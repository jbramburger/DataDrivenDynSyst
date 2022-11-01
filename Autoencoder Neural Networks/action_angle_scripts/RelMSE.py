"""Calculate relative mean squared error."""
from tensorflow import keras
import tensorflow as tf


class RelMSE(keras.losses.Loss):
    """Subclass the keras Loss class."""

    def __init__(self, denom_nonzero=1e-5, **kwargs):
        """Compute relative mean squared error between labels and preds.

        Arguments:
            denom_nonzero -- a small nonzero term to add to the denominator
                to avoid dividing by zero
            **kwargs -- additional keyword arguments
        """
        super().__init__(**kwargs)

        self.denom_nonzero = denom_nonzero

    def call(self, y_true, y_pred):
        """Calculate relative MSE given true and predicted values."""
        # Compute the MSE of each example
        mse = tf.reduce_mean(tf.square(y_pred - y_true), axis=-1)

        # Compute the mean of squares of the true values
        true_norm = tf.reduce_mean(tf.square(y_true), axis=-1)

        # Ensure there are no 'zero' values in the denominator before division
        true_norm += self.denom_nonzero

        # Compute relative MSE of each example
        err = tf.truediv(mse, true_norm)

        # Compute mean over batch
        err = tf.reduce_mean(err, axis=-1)

        # Return the error
        return err
