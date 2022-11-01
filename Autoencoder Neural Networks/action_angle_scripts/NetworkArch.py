"""Creates the network architecture for linearizing PDEs."""
import numpy as np
import tensorflow as tf
from tensorflow import keras

from action_angle_scripts.DenseResBlock import DenseResBlock
from action_angle_scripts.RelMSE import RelMSE

@tf.function
def print_fn(x):
  # Note that `print_trace` is a TensorFlow Op. See the next section for details
  print_trace = tf.print(
      "`input` has type", type(x), "and shape", tf.shape(x))


class NetworkArch(keras.Model):
    """Subclass the Keras Model class."""

    def __init__(self,
                 n_inputs=128,
                 n_latent=21,
                 len_time=51,
                 num_shifts=50,
                 num_shifts_middle=50,
                 outer_encoder=DenseResBlock(),
                 outer_decoder=DenseResBlock(),
                 inner_config=dict(),
                 inner_loss_weights=[1, 1],
                 L_diag=False,
                 L_param=False,
                 auxiliary_network=None,
                 train_autoencoder_only=False,
                 delta_t=None,
                 **kwargs):
        """
        Create network architecture for linearizing PDEs.

        Arguments:
            n_inputs -- the number of inputs to the network
                (spatial discretization of the PDE)
            n_latent -- the dimensionality of the latent space
                (i.e. number of Koopman eigenfunctions in expansion)
            len_time -- number of time steps for each trajectory in data
            num_shifts -- the number of time steps in the future each network
                input will predict when calculating prediction loss
            num_shifts_middle -- the number of time steps in the future each
                network input will predict when calculating the linearity loss
            outer_encoder -- a Keras Layer or Model with the architecture for
                the outer encoder
            outer_decoder -- a Keras Layer or Model with the architecture for
                the outer decoder (typically the same as the outer encoder)
            inner_config -- Python dictionary with keyword arguments for
                the inner encoder and decoder
            inner_loss_weights -- Python list of length 2 with weights for the
                inner autoencoder loss and linearity loss
            L_diag -- Boolean that determines whether the dynamics matrix L is
                constrained to be diagonal
            L_param -- Boolean that determines whether the dynamics matrix L is 
                parameterized by the encoded space
            auxiliary_network -- a Keras Layer or Model with the architecture for
                the auxiliary network (needed if L_param == True)
            train_autoencoder_only -- Boolean that determines whether only the
                autoencoder losses are used for training. It is recommendeded
                that you do several epochs of pretraining with only the
                autoencoder losses and then set this option to False to include
                the prediction and linearity losses.
            **kwargs -- additional keyword arguments. Can be used to name the
                Model.
        """
        super().__init__(**kwargs)

        self.n_inputs = n_inputs
        self.n_latent = n_latent
        self.len_time = len_time
        self.num_shifts = num_shifts
        self.num_shifts_middle = num_shifts_middle
        self.outer_encoder = outer_encoder
        self.outer_decoder = outer_decoder

        # Create the inner encoder layer
        self.inner_encoder = keras.layers.Dense(
            n_latent,
            name='inner_encoder',
            activation=None,
            use_bias=False,
            kernel_initializer=identity_init,
            **inner_config)

        # The dynamics matrix, initialized as identity
        self.L = tf.Variable(tf.eye(n_latent), trainable=True)

        self.L_param = L_param
        # if the dynamics matrix L is parameterized, we need an auxiliary network
        # that maps from the encoded state to L
        if self.L_param:
            self.auxiliary_network = auxiliary_network
            self.delta_t = delta_t

        # Create the inner decoder layer
        self.inner_decoder = keras.layers.Dense(
            n_inputs,
            name='inner_decoder',
            activation=None,
            use_bias=False,
            kernel_initializer=identity_init,
            **inner_config)

        self.RelMSE = RelMSE(name='RelMSE')
        self.inner_loss_weights = inner_loss_weights
        self.L_diag = L_diag
        
        self.train_autoencoder_only = train_autoencoder_only

    def call(self, inputs):
        """
        Run given inputs through the newtork.

        Arguments:
            inputs:
            inputs -- A Numpy array or Tensorflow tensor with shape
                (number_of_trajectories, len_time, n_inputs)
            outputs:
            autoencoder_output -- The output of running each time
                step of each trajectory through the autoencoder
            outer_auto_output -- The output of running each time step
                of each trajectory through the outer autoencoder
            predictions -- The predictions for num_steps steps into the
                future for all trajectories and all time steps for which
                num_steps steps into the future are in the data

        """
        # For the shape of the outputs for prediction and linearity loss
        len_pred = self.len_time - self.num_shifts
        len_lin = self.len_time - self.num_shifts_middle

        # Create arrays for prediction and linearity
        # inputs shape was (number_of_trajectories, len_time, n_inputs)
        # pred_inputs shape is (number_of_trajectories, len_pred, n_inputs)
        pred_inputs = inputs[:, :len_pred, :]
        # lin_inputs shape is (number_of_trajectories, len_lin, n_inputs)
        lin_inputs = inputs[:, :len_lin, :]

        # Inputs for "exact" solutions for prediction and linearity loss
        pred_exact = stack_predictions(inputs, self.num_shifts)
        lin_advanced = stack_predictions(inputs, self.num_shifts_middle)

        #  Reshape inputs as 2D arrays
        auto_inputs, pred_inputs, lin_inputs, lin_advanced = reshape_inputs(
            (inputs, pred_inputs, lin_inputs, lin_advanced))
        # now pred_inputs is reshaped to (number_of_trajectories x len_pred, n_inputs)
        # and lin_inputs is reshaped to (number_of_trajectories x len_lin, n_inputs)

        # Autoencoder
        partially_encoded = self.outer_encoder(auto_inputs)
        fully_encoded = self.inner_encoder(partially_encoded)
        partially_decoded = self.inner_decoder(fully_encoded)
        autoencoder_output = self.outer_decoder(partially_decoded)

        autoencoder_output = tf.reshape(autoencoder_output,
                                        [-1, self.len_time, self.n_inputs])

        # Outer Autoencoder
        outer_auto_output = self.outer_decoder(partially_encoded)

        outer_auto_output = tf.reshape(outer_auto_output,
                                       [-1, self.len_time, self.n_inputs])

        # Inner Autoencoder Loss
        self.add_loss(self.inner_loss_weights[0]
                      * self.RelMSE(partially_encoded, partially_decoded))

        # If training autoencoder only, output results
        if self.train_autoencoder_only:
            predictions = 0 * pred_exact
            return autoencoder_output, outer_auto_output, predictions

        # Set dynamics matrix L
        if self.L_diag:
            Lmat = tf.linalg.diag(tf.linalg.diag_part(self.L))
        else:
            Lmat = self.L

        # Prediction
        predictions_list = []
        # pred_inputs has shape (number_of_trajectories x len_pred, n_inputs) so part_encoded_pred also has that shape
        part_encoded_pred = self.outer_encoder(pred_inputs)
        
        # current_encoded has shape [number_of_trajectories x len_pred, n_latent]
        current_encoded = self.inner_encoder(part_encoded_pred)

        # if dynamics matrix L is parameterized by the encoding, set that here
        # post-AA TODO: generalize the auxiliary network: currently assuming fixed L for each trajectory & no mu
        # TODO: generalize structure when L is parameterized: currently only allowing 2x2 complex conjugate block

        #print_fn(current_encoded)
        if self.L_param:
            # current_encoded has shape [number_of_trajectories x len_pred, n_latent]
            # output of auxiliary_network has shape [number_of_trajectories x len_pred, n_latent/2] (no mus)
            # Lmat has shape [number_of_trajectories x len_pred, n_latent, n_latent]
            print("using parameterized L")
            # make auxiliary network function of y_1^2 + y_2^2, etc. instead of y directly
            # TODO: allow more than one oscillator 
            radius = tf.reduce_sum(tf.square(current_encoded), axis=1, keepdims=True)
            # radius should have shape [number_of_trajectories x len_pred, n_latent/2] (for now: [number_of_trajectories x len_pred, 1])
            Lmat = form_complex_conjugate_block(self.auxiliary_network(radius), self.delta_t)
        #print_fn(Lmat)
        for shift in range(self.num_shifts):
            # in case of L_param, matmul on [number_of_trajectories x len_pred, n_latent] times [number_of_trajectories x len_pred, n_latent, n_latent] 
            # should mean advanced_encoded has shape [number_of_trajectories x len_pred, n_latent]
            advanced_encoded = tf.linalg.matvec(Lmat, current_encoded)
            #print_fn(advanced_encoded)
            adv_part_decoded = self.inner_decoder(advanced_encoded)
            advanced_decoded = self.outer_decoder(adv_part_decoded)
            predictions_list.append(tf.reshape(advanced_decoded,
                                               [-1, len_pred, self.n_inputs]))
            current_encoded = tf.identity(advanced_encoded)
        predictions = tf.concat(predictions_list, axis=1)

        # Linearity predictions
        linearity_list = []
        # lin_inputs has shape (number_of_trajectories x len_lin, n_inputs) so part_encoded_lin also has that shape
        part_encoded_lin = self.outer_encoder(lin_inputs)

        # current_encoded has shape [number_of_trajectories x len_lin, n_latent]
        current_encoded = self.inner_encoder(part_encoded_lin)

        # if dynamics matrix L is parameterized by the encoding, set that here
        # Note: Lmat is potentially different than in above case because of num_shifts vs. num_shifts middle
        if self.L_param:
            # slightly different than above: Lmat has shape [number_of_trajectories x len_lin, n_latent, n_latent]
            # make auxiliary network function of y_1^2 + y_2^2, etc. instead of y directly
            # TODO: allow more than one oscillator 
            radius = tf.reduce_sum(tf.square(current_encoded), axis=1, keepdims=True)
            Lmat = form_complex_conjugate_block(self.auxiliary_network(radius), self.delta_t)

        for shift in range(self.num_shifts_middle):
            advanced_encoded = tf.linalg.matvec(Lmat, current_encoded)
            current_encoded = tf.identity(advanced_encoded)
            linearity_list.append(tf.reshape(current_encoded,
                                             [-1, len_lin, self.n_latent]))
        lin_pred = tf.concat(linearity_list, axis=1)

        lin_part_encoded = self.outer_encoder(lin_advanced)
        lin_exact = self.inner_encoder(lin_part_encoded)
        lin_exact = tf.reshape(
            lin_exact,
            [-1, self.num_shifts_middle * len_lin, self.n_latent])

        # Add Linearity loss
        self.add_loss(self.inner_loss_weights[1]
                      * self.RelMSE(lin_exact, lin_pred))

        return autoencoder_output, outer_auto_output, predictions


def identity_init(shape, dtype=tf.float32):
    """Initialize weight matrices as identity-like matrices."""
    n_rows = shape[0]
    n_cols = shape[1]
    if n_rows >= n_cols:
        A = np.zeros((n_rows, n_cols), dtype=np.float32)
        for col in range(n_cols):
            for row in range(col, col + n_rows - n_cols + 1):
                A[row, col] = 1.0 / (n_rows - n_cols + 1)
    else:
        A = np.zeros((n_rows, n_cols), dtype=np.float32)
        for row in range(n_rows):
            for col in range(row, row + n_cols - n_rows + 1):
                A[row, col] = 1.0 / (n_cols - n_rows + 1)
    return A


def reshape_inputs(inputs):
    """Reshape inputs to be 2D arrays."""
    input_list = []
    for data in inputs:
        input_list.append(tf.reshape(data, [-1, data.shape[-1]]))
    return tuple(input_list)


def stack_predictions(data, num_shifts):
    """Stack inputs for linearity or prediction loss."""
    len_pred = data.shape[1] - num_shifts
    prediction_list = []
    for j in range(num_shifts):
        prediction_list.append(data[:, j + 1:j + 1 + len_pred, :])
    prediction_tensor = tf.concat(prediction_list, axis=1)

    return prediction_tensor

def form_complex_conjugate_block(omegas, delta_t):
    """Form a 2x2 block for a complex conj. pair of eigenvalues, but for each example, so dimension [None, 2, 2]
    2x2 Block is
    [cos(omega * delta_t), -sin(omega * delta_t)
    sin(omega * delta_t), cos(omega * delta_t)]
    Arguments:
        omegas -- array of parameters for blocks: freq. (omega), size [None, 1] (no mu for action angles)
        delta_t -- time step in trajectories from input data
    Returns:
        stack of 2x2 blocks, size [None, 2, 2], where first dimension matches first dimension of omegas
    Side effects:
        None
    """
    entry11 = tf.cos(omegas[:, 0] * delta_t)
    entry12 = tf.sin(omegas[:, 0] * delta_t)
    row1 = tf.stack([entry11, -entry12], axis=1)  # [None, 2]
    row2 = tf.stack([entry12, entry11], axis=1)  # [None, 2]
    return tf.stack([row1, row2], axis=2)  # [None, 2, 2] put one row below other
