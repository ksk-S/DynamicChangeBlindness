import ball_reader
import csv
#from easyesn import PredictionESN
from FeedbackPredictionESN import PredictionESN
#import FeedbackPredictionESN
import matplotlib.pyplot as plt
import numpy as np
import os
import random
import argparse
import optuna


''' Point process function '''
def point_window(n_bins=4):
    y = np.zeros((n_bins)).astype("float")
    y[int(n_bins/2)] = 1.
    return y


''' Gaussian window function '''
def gaussian_window(n_bins=4):
    sigma = 1.
    x = np.arange(-1, 1, 2./float(n_bins))
    y = np.exp((-1)*np.power(x, 2)/sigma**2)
    return y


def create_data(dim_input=3, num_train_epis=1000, num_test_epis=100, min_epis_interval=5,
                max_epis_interval=5, len_response_interval=3, stimulus_bins=7, window_func=gaussian_window):
    # assert min_epis_interval <= max_epis_interval, "The maximum of episode-interval must be larger tha the minimum."
    # assert min_epis_interval > stimulus_bins, "Stimulus bins must be larger than the minimum of episode-interval."
    # assert num_train_epis > num_test_epis, "The number of training episodes must be larger than the one of testing."

    train_data = None
    test_data = None

    # Training data
    for i in range(num_train_epis):
        response_id = np.random.randint(0, dim_input)
        # response_id = 0
        len_initial_interval = random.randint(min_epis_interval, max_epis_interval)
        total_len_epis = len_initial_interval + stimulus_bins*2 + len_response_interval
        epis_arr = np.zeros((total_len_epis, dim_input)).astype("float")
        for _idx in range(dim_input):
            if _idx == response_id:
                epis_arr[total_len_epis-stimulus_bins:, _idx] = window_func(n_bins=stimulus_bins) # Response
            else:
                epis_arr[len_initial_interval:len_initial_interval+stimulus_bins, _idx] = window_func(n_bins=stimulus_bins) # Stimulus
        if train_data is None:
            train_data = epis_arr
        else:
            train_data = np.concatenate([train_data, epis_arr], axis=0)
    
    # Test data
    for i in range(num_test_epis):
        stimulus_id = np.random.randint(0, dim_input)
        # stimulus_id = 1
        len_initial_interval = random.randint(min_epis_interval, max_epis_interval)
        total_len_epis = len_initial_interval + stimulus_bins*2 + len_response_interval
        epis_arr = np.zeros((total_len_epis, dim_input)).astype("float")
        for _idx in range(dim_input):
            if _idx == stimulus_id:
                epis_arr[len_initial_interval:len_initial_interval+stimulus_bins, _idx] = window_func(n_bins=stimulus_bins) # Stimulus
            else:
                epis_arr[total_len_epis-stimulus_bins:, _idx] = window_func(n_bins=stimulus_bins) # Response
        if test_data is None:
            test_data = epis_arr
        else:
            test_data = np.concatenate([test_data, epis_arr], axis=0)
    
    train_X = train_data[:-1, :]
    train_Y = train_data[1:, :]
    test_X = test_data[:-1, :]
    test_Y = test_data[1:, :]

    return train_X, train_Y, test_X, test_Y


# def create_data(dim_input=3, num_train_epis=1000, num_test_epis=100, min_epis_interval=5,
#                 max_epis_interval=5, len_response_interval=3, stimulus_bins=7, window_func=gaussian_window):
#     # assert min_epis_interval <= max_epis_interval, "The maximum of episode-interval must be larger tha the minimum."
#     # assert min_epis_interval > stimulus_bins, "Stimulus bins must be larger than the minimum of episode-interval."
#     # assert num_train_epis > num_test_epis, "The number of training episodes must be larger than the one of testing."

#     train_X = None
#     train_Y = None
#     test_X = None
#     test_Y = None

#     # Training data
#     for i in range(num_train_epis):
#         response_id = np.random.randint(0, dim_input)
#         # response_id = 0
#         len_initial_interval = random.randint(min_epis_interval, max_epis_interval)
#         total_len_epis = len_initial_interval + stimulus_bins*2 + len_response_interval
#         epis_arr_x = np.zeros((total_len_epis, dim_input)).astype("float")
#         epis_arr_y = np.zeros((total_len_epis, dim_input)).astype("float")
#         for _idx in range(dim_input):
#             if _idx == response_id:
#                 epis_arr_y[total_len_epis-stimulus_bins:, _idx] = window_func(n_bins=stimulus_bins) # Response
#             else:
#                 epis_arr_x[len_initial_interval:len_initial_interval+stimulus_bins, _idx] = window_func(n_bins=stimulus_bins) # Stimulus
#         if train_X is None:
#             train_X = epis_arr_x
#             train_Y = epis_arr_y
#         else:
#             train_X = np.concatenate([train_X, epis_arr_x], axis=0)
#             train_Y = np.concatenate([train_Y, epis_arr_y], axis=0)
    
#     # Test data
#     for i in range(num_test_epis):
#         stimulus_id = np.random.randint(0, dim_input)
#         # stimulus_id = 1
#         len_initial_interval = random.randint(min_epis_interval, max_epis_interval)
#         total_len_epis = len_initial_interval + stimulus_bins*2 + len_response_interval
#         epis_arr_x = np.zeros((total_len_epis, dim_input)).astype("float")
#         epis_arr_y = np.zeros((total_len_epis, dim_input)).astype("float")
#         for _idx in range(dim_input):
#             if _idx == stimulus_id:
#                 epis_arr_x[len_initial_interval:len_initial_interval+stimulus_bins, _idx] = window_func(n_bins=stimulus_bins) # Stimulus
#             else:
#                 epis_arr_y[total_len_epis-stimulus_bins:, _idx] = window_func(n_bins=stimulus_bins) # Response
#         if test_X is None:
#             test_X = epis_arr_x
#             test_Y = epis_arr_y
#         else:
#             test_X = np.concatenate([test_X, epis_arr_x], axis=0)
#             test_Y = np.concatenate([test_Y, epis_arr_y], axis=0)
    
#     return train_X, train_Y, test_X, test_Y


def train_tracker(args, plot_path, dim_input=3):
    plot_name = ''
    if not os.path.exists(plot_path+"networks"):
        os.makedirs(plot_path+"networks")
    
    # Choosing a window function for stimulus and response signals
    if args.use_gauss_window:
        window_func = gaussian_window
    else:
        window_func = point_window
    
    window_func = gaussian_window
    
    # Creating data
    train_X, train_Y, test_X, test_Y = create_data(
        dim_input=dim_input,
        num_train_epis=args.num_train_epis,
        num_test_epis=args.num_test_epis,
        min_epis_interval=args.min_epis_interval,
        max_epis_interval=args.max_epis_interval,
        len_response_interval=args.len_response_interval,
        stimulus_bins=args.stimulus_bins,
        window_func=window_func
    )

    # Creating data for training validation
    valid_X, valid_Y, __x, __y = create_data(
        dim_input=dim_input,
        num_train_epis=5,
        num_test_epis=1,
        min_epis_interval=args.min_epis_interval,
        max_epis_interval=args.max_epis_interval,
        len_response_interval=args.len_response_interval,
        stimulus_bins=args.stimulus_bins,
        window_func=window_func
    )

    #save original data
    header = ""
    header = header + "presence"
    np.savetxt(plot_path + "track_original.csv", test_Y, delimiter=",", header=header, fmt='%3f', comments='')

    for i in range(args.networks_count):
        esn = PredictionESN(
            n_input=dim_input, n_output=dim_input,
            n_reservoir=args.num_hiddens,
            leakingRate=args.leaking_rate,
            spectralRadius=args.spectral_radius,
            regressionParameters=[1e-2],
            feedback=args.feedback
        )
        # Fitting training data
        esn.fit(train_X, train_Y, transientTime=0, verbose=1)
        # Confirming whether the training is succeeded or not
        y_valid_pred = esn.predict(valid_X, verbose=1)
        # Predicting with test data
        y_test_pred = esn.predict(test_X, verbose=1)

        # Saving prediction data
        np.savetxt(plot_path + "track_" + str(i) + ".csv", y_test_pred, delimiter=",", header=header, fmt='%3f', comments='')
        #save network
        esn.save(path=plot_path + "networks/" + str(i))
        
        # Prediction Plot
        plt.figure('Training validation', figsize=(14,21)).clear()
        #plt.yscale('log')
        plt.subplot(3, 1, 1)
        plt.plot(valid_Y[:,0], color='red', linewidth=5, label='Target Value')
        plt.plot(y_valid_pred[:,0], color='blue', linestyle="--", linewidth=1, label='ESN Prediction')
        plt.ylim([-2.0, 2.0])
        plt.legend()
        plt.subplot(3, 1, 2)
        plt.plot(valid_Y[:,1], color='red', linewidth=5, label='Target Value')
        plt.plot(y_valid_pred[:,1], color='blue', linestyle="--", linewidth=1, label='ESN Prediction')
        plt.ylim([-2.0, 2.0])
        plt.legend()
        plt.subplot(3, 1, 3)
        plt.plot(valid_Y[:,2], color='red', linewidth=5, label='Target Value')
        plt.plot(y_valid_pred[:,2], color='blue', linestyle="--", linewidth=1, label='ESN Prediction')
        plt.ylim([-2.0, 2.0])
        plt.legend()
        figname = plot_path + plot_name + "_validation_" + str(i)  + ".png"
        plt.savefig(figname)
        print('\t[+]Plot saved in', figname)
        
        plt.figure('Prediction', figsize=(14,21)).clear()
        #plt.yscale('log')
        plt.subplot(3, 1, 1)
        plt.plot(test_Y[:,0], color='red', linewidth=5, label='Target Value')
        plt.plot(y_test_pred[:,0], color='blue', linestyle="--", linewidth=1, label='ESN Prediction')
        plt.ylim([-2.0, 2.0])
        plt.legend()
        plt.subplot(3, 1, 2)
        plt.plot(test_Y[:,1], color='red', linewidth=5, label='Target Value')
        plt.plot(y_test_pred[:,1], color='blue', linestyle="--", linewidth=1, label='ESN Prediction')
        plt.ylim([-2.0, 2.0])
        plt.legend()
        plt.subplot(3, 1, 3)
        plt.plot(test_Y[:,2], color='red', linewidth=5, label='Target Value')
        plt.plot(y_test_pred[:,2], color='blue', linestyle="--", linewidth=1, label='ESN Prediction')
        plt.ylim([-2.0, 2.0])
        plt.legend()
        figname = plot_path + plot_name + "_prediction_" + str(i)  + ".png"
        plt.savefig(figname)
        print('\t[+]Plot saved in', figname)


def search_func(trial):
    plot_path='cb_results/minimal_model/'
    if not os.path.exists(plot_path+"networks"):
        os.makedirs(plot_path+"networks")
    
    # window_func = gaussian_window
    window_func = point_window
    dim_input = 3
    num_train_epis = 1000
    num_test_epis = 100
    min_epis_interval = 10
    stimulus_bins = 1

    num_hiddens = 200
    feedback = False
    networks_count = 1

    leaking_rate = trial.suggest_uniform('leaking_rate', 0.0, 1.0)
    spectral_radius = trial.suggest_uniform('spectral_radius', 0.0, 1.0)
    max_epis_interval = trial.suggest_int("max_epis_interval", 10, 20)
    len_response_interval = trial.suggest_int("len_response_interval", 3, 5)

    train_X, train_Y, test_X, test_Y = create_data(
        dim_input=dim_input,
        num_train_epis=num_train_epis,
        num_test_epis=num_test_epis,
        min_epis_interval=min_epis_interval,
        max_epis_interval=max_epis_interval,
        len_response_interval=len_response_interval,
        stimulus_bins=stimulus_bins,
        window_func=window_func
    )

    esn = PredictionESN(
        n_input=dim_input,
        n_output=dim_input,
        n_reservoir=num_hiddens,
        leakingRate=leaking_rate,
        spectralRadius=spectral_radius,
        regressionParameters=[1e-2],
        feedback=feedback
    )

    # Fitting training data
    esn.fit(train_X, train_Y, transientTime=0, verbose=1)
    # Predicting with test data
    y_test_pred = esn.predict(test_X, verbose=1)
    
    # Saving prediction data
    np.savetxt(plot_path + "track_" + str(trial) + ".csv", y_test_pred, delimiter=",", header="presence", fmt='%3f', comments='')
    #save network
    esn.save(path=plot_path + "networks/" + str(trial))
    
    cur_score = test_Y - y_test_pred
    cur_score = np.power(cur_score, 2)
    cur_score = np.mean(cur_score)

    return cur_score


if __name__=='__main__':
    '''
    {'leaking_rate': 0.9958270662907154, 'spectral_radius': 0.9307680421663365, 'max_epis_interval': 20, 'len_response_interval': 4}
    '''
    parser = argparse.ArgumentParser()
    # About data
    parser.add_argument('--num_train_epis', type=int, default=1000)
    parser.add_argument('--num_test_epis', type=int, default=2)
    parser.add_argument('--min_epis_interval', type=int, default=4)
    parser.add_argument('--max_epis_interval', type=int, default=4)
    parser.add_argument('--len_response_interval', type=int, default=2)
    parser.add_argument('--stimulus_bins', type=int, default=7)
    parser.add_argument('--use_gauss_window', action='store_true', default=False)
    # About echo state networks
    parser.add_argument('--num_hiddens', type=int, default=200)
    parser.add_argument('--leaking_rate', type=float, default=0.9958270662907154)
    parser.add_argument('--spectral_radius', type=float, default=0.9307680421663365)
    parser.add_argument('--feedback', action='store_true', default=False)
    parser.add_argument('--networks_count', type=int, default=1)
    args = parser.parse_args()
    args.use_gauss_window = True
    
    train_tracker(args, plot_path='cb_results/minimal_model/', dim_input=3)

    # study = optuna.create_study()
    # study.optimize(search_func, n_trials=1000)