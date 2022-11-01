"""Experiment Helper Functions."""
import random as r
import json
import sys
import os

import numpy as np
import tensorflow as tf
from tensorflow import keras

# Add the architecture path for the networkarch and RelMSE
from action_angle_scripts.NetworkArch import NetworkArch


def construct_network(**architecture_config):
    """Construct neural network architecture."""
    return NetworkArch(**architecture_config)


def getdatasize(data_file_prefix):
    """Get size of data files."""
    return np.load("{}_train1_x.npy".format(data_file_prefix)).shape


def get_data(data_file_prefix, data_train_len, num_shifts):
    """Load and organize data."""
    # Load in the data
    data_train = np.load("{}_train1_x.npy".format(data_file_prefix))
    for k in range(2, data_train_len + 1):
        new_data = np.load("{}_train{}_x.npy".format(data_file_prefix, k))
        data_train = np.vstack([data_train, new_data])
    data_val = np.load("{}_val_x.npy".format(data_file_prefix))

    # Stack data for prediction loss
    train_pred = stack_predictions(data_train, num_shifts)
    val_pred = stack_predictions(data_val, num_shifts)

    # Tensors of zeros for training autoencoder only
    train_zeros = np.zeros(train_pred.shape)
    val_zeros = np.zeros(val_pred.shape)

    return (data_train, data_val, train_zeros,
            val_zeros, train_pred, val_pred)


def evaluate_initial_models(save_prefix, all_data,
                            train_opts, network_config):
    """Train 20 models and choose best one."""
    # Extract the training/validation data
    (data_train, data_val, train_zeros,
     val_zeros, train_pred, val_pred) = all_data

    # Gather the relevant training options
    aec_only_epochs = train_opts['aec_only_epochs']
    init_full_epochs = train_opts['init_full_epochs']
    num_init_models = train_opts['num_init_models']
    loss_fn = train_opts['loss_fn']
    opt = train_opts['optimizer']
    optimizer_opts = train_opts['optimizer_opts']
    batch_size = train_opts['batch_size']
    loss_weights = train_opts['loss_weights']
    inner_loss_weights = [loss_weights[i] for i in [2, 4]]
    outer_loss_weights = [loss_weights[i] for i in [0, 1, 3]]

    # Set up results dictionary
    results = {'full_hist': [],
               'aec_hist': [],
               'lr': [],
               'best_loss': [],
               'model_path': []}

    # For loop for generating, training, and evaluating the initial models
    for i in range(num_init_models):
        # Randomly selected learning rate
        lr = 10**(-r.uniform(3, 6))

        # Create a model, initially only train autoencoders
        model = construct_network(train_autoencoder_only=True,
                                  inner_loss_weights=inner_loss_weights,
                                  **network_config)
        # Compile the model
        model.compile(loss=3 * [loss_fn],
                      optimizer=opt(lr=lr, **optimizer_opts),
                      loss_weights=outer_loss_weights)

        # Use checkpointing
        checkpoint_path_aec = save_prefix + 'checkpoint_aec_{}'.format(i)
        cbs_aec = [keras.callbacks.ModelCheckpoint(checkpoint_path_aec,
                                                   save_weights_only=True,
                                                   monitor='val_loss',
                                                   save_best_only=True)]
        #  Fit autoencoder-only model
        aec_hist = model.fit(x=data_train,
                             y=[data_train, data_train, train_zeros],
                             validation_data=(data_val,
                                              [data_val, data_val, val_zeros]),
                             callbacks=cbs_aec, batch_size=batch_size,
                             epochs=aec_only_epochs, verbose=True)

        # Re-load weights with best validation loss
        model.load_weights(checkpoint_path_aec)

        # Now set the model to train with prediction losses
        model.train_autoencoder_only = False

        # Re-compile the model
        model.compile(loss=3 * [loss_fn],
                      optimizer=opt(lr=lr, **optimizer_opts),
                      loss_weights=outer_loss_weights)

        # Train full model
        checkpoint_path_full = save_prefix + 'checkpoint_{}'.format(i)
        cbs = [keras.callbacks.ModelCheckpoint(checkpoint_path_full,
                                               save_weights_only=True,
                                               monitor='val_loss',
                                               save_best_only=True)]
        # Fit the full model
        full_hist = model.fit(x=data_train,
                              y=[data_train, data_train, train_pred],
                              validation_data=(data_val,
                                               [data_val, data_val, val_pred]),
                              callbacks=cbs, batch_size=batch_size,
                              epochs=init_full_epochs)

        # Re-load weights with best validation loss
        model.load_weights(checkpoint_path_full)

        # Evaluate the model to get final validation loss
        best_loss = model.evaluate(x=data_val,
                                   y=[data_val, data_val, val_pred],
                                   verbose=False)

        # Save the model
        model_path = save_prefix + "model_{}".format(i)
        model.save(model_path)

        # Append the results to the results list
        results['full_hist'].append(full_hist.history.copy())
        results['aec_hist'].append(aec_hist.history.copy())
        results['lr'].append(lr)
        results['best_loss'].append(best_loss[0])
        results['model_path'].append(model_path)

        # Delete the model variable and clear_session to remove any graph
        del model
        tf.keras.backend.clear_session()

    # Select the best model from the loop
    best_model_idc = np.argmin(results['best_loss'])
    best_model_path = results['model_path'][best_model_idc]

    # Return the best model's path
    return results, best_model_path


def train_final_model(model_path,
                      save_prefix,
                      all_data,
                      train_opts,
                      custom_objects):
    """Load best initial model and train until convergence."""
    # Gather the relevant training options
    best_model_epochs = train_opts['best_model_epochs']
    batch_size = train_opts['batch_size']

    # Extract the training/validation data
    (data_train, data_val, train_zeros,
     val_zeros, train_pred, val_pred) = all_data

    # Load the model
    model = tf.keras.models.load_model(model_path,
                                       custom_objects=custom_objects)

    # Set the place to save the checkpoint model weights
    checkpoint_model_path = save_prefix + 'checkpoint_final'

    # Use checkpointing
    cbs = [keras.callbacks.ModelCheckpoint(checkpoint_model_path,
                                           save_weights_only=True,
                                           monitor='val_loss',
                                           save_best_only=True)]

    # Train the model
    hist = model.fit(x=data_train,
                     y=[data_train, data_train, train_pred],
                     validation_data=(data_val,
                                      [data_val, data_val, val_pred]),
                     callbacks=cbs, batch_size=batch_size,
                     epochs=best_model_epochs)

    # Re-load the best model weights
    model.load_weights(checkpoint_model_path)

    # Save the model
    model_path = save_prefix + 'final_model'
    model.save(model_path)

    return hist.history, model_path


def save_results(results_path, random_seed,
                 model_path, custom_objects,
                 final_hist, init_hist):
    """Save the results."""
    # Load and save the model
    model = tf.keras.models.load_model(model_path,
                                       custom_objects=custom_objects)
    model.save(results_path + 'final_model')
    print("Best model saved to:", model_path)

    # Export the initial model training dictionary
    hist_filepath = results_path + "initial_pool_results.json"
    json.dump(init_hist, open(hist_filepath, 'w'))

    # Export the final model training dictionary
    final_hist['random_seed'] = random_seed
    hist_filepath = results_path + "final_model_history.json"
    json.dump(final_hist, open(hist_filepath, 'w'))

    print("Exported training dictionaries to: ", results_path)


def check_for_directories(expt_name):
    """Create necessary directories if they do not exist."""
    pardir = os.path.abspath(os.pardir)
    dirs = ['logs',
            'model_weights',
            'model_weights' + os.sep + expt_name,
            'results',
            'results' + os.sep + expt_name,
            ]
    for dirname in dirs:
        os.makedirs(pardir + os.sep + dirname, exist_ok=True)


def stack_predictions(data, num_shifts):
    """Create tensors to be used as inputs for prediction/linearity losses."""
    len_pred = data.shape[1] - num_shifts
    prediction_list = []
    for j in range(num_shifts):
        prediction_list.append(data[:, j + 1:j + 1 + len_pred, :])
    prediction_tensor = np.concatenate(prediction_list, axis=1)

    return prediction_tensor


def run_experiment(random_seed, expt_name, data_file_prefix,
                   training_options, network_config, custom_objects):
    """Run experiment for Koopman autoencoder."""
    # Assign a random number generator seed for learning rates
    r.seed(random_seed)

    # Create necessary directories
    check_for_directories(expt_name)

    # Get the training data
    all_data = get_data(data_file_prefix,
                        training_options['data_train_len'],
                        network_config['num_shifts'])

    # Set the prefix for where to save the results/checkpointed models
    save_prefix = '../model_weights/{}/'.format(expt_name)

    # Step 1 -- Train a collection of initial models
    # Autoencoders-only, then full model
    # This method returns the file path to the best model:
    init_hist, model_path = evaluate_initial_models(save_prefix,
                                                    all_data,
                                                    training_options,
                                                    network_config)

    # Step 2 -- Load the best model, and train for the full time
    # Load the best model
    final_hist, model_path = train_final_model(model_path, save_prefix,
                                               all_data,
                                               training_options,
                                               custom_objects)

    # Step 3 -- Save the results
    results_path = '../results/{}/'.format(expt_name)
    save_results(results_path, random_seed,
                 model_path, custom_objects,
                 final_hist, init_hist)
