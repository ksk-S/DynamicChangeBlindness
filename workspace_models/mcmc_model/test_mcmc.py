import numpy as np
import pandas as pd
import pystan
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

import warnings
warnings.filterwarnings('ignore')

import csv
import glob

def create_data(ts, f1, f2, theta1, theta2, A1, A2, sigma=0.02):
    outs = []
    for t in ts:
        val = A1*math.sin(2*math.pi*f1*t+theta1) + A2*math.sin(2*math.pi*f2*t+theta2) + np.random.normal(0.0, sigma)
        outs.append(val)
    return np.array(outs)


def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))


def relu(x):
    if x > 0:
        return x
    else:
        return 0.


sm = pystan.StanModel(file='./stan/wave_minimal_mixture_phase_v2.stan')

def main(file_id=0):
    # w_dict = {'w1': 2.129664768615264, 'w2': 8.852169933587668, 'w3': 7.246868900669156, 'w4': 6.140351315929692, 'w5': -6.893550507176817, 'w6': -3.536844111212398, 'w7': 0.4368130392565739, 'w8': -9.984198843831692, 'w9': -2.555181532914446, 'w10': 8.102549737405944, 'w11': 7.006321079909512, 'w12': 9.564507619066738, 'w13': 7.629611432175047, 'w14': 6.872460616740369, 'w15': -7.6622867363527325}
    # w_dict = {'w1': -1.321608147937682, 'w2': -1.5152116491564795, 'w3': -9.966715996596518, 'w4': 6.59422839639349, 'w5': 9.815788344261179, 'w6': -0.011891300816384232, 'w7': 9.648864825269271, 'w8': -4.918440864817081, 'w9': -1.2781061267131406, 'w10': 2.222622486704873, 'w11': 5.170732571748707, 'w12': -9.472978755216273, 'w13': 7.5546316778520435, 'w14': 7.078069839919416, 'w15': -7.685538254667767}
    w_dict = {'w1': 2.9594431735720113, 'w2': 2.9834750286183813, 'w3': -2.0373151373032354, 'w4': -2.9562080390382732, 'w5': 1.3535656270456662, 'w6': 2.9842366488694174, 'w7': 2.727711056806083, 'w8': 2.9747890067311733, 'w9': 0.8716217311959202, 'w10': 2.4465328655326286, 'w11': -2.938165083407799, 'w12': 2.376187463010098, 'w13': 1.6636081984702027, 'w14': 2.969931340836359, 'w15': 2.9834357626984325, 'w_sigma': 7.544885306559859}
    w1 = w_dict['w1']
    w2 = w_dict['w2']
    w3 = w_dict['w3']
    w4 = w_dict['w4']
    w5 = w_dict['w5']
    w6 = w_dict['w6']
    w7 = w_dict['w7']
    w8 = w_dict['w8']
    w9 = w_dict['w9']
    w10 = w_dict['w10']
    w11 = w_dict['w11']
    w12 = w_dict['w12']
    w13 = w_dict['w13']
    w14 = w_dict['w14']
    w15 = w_dict['w15']
    w_sigma = w_dict['w_sigma']
    W_1 = np.array([w1, w2, w3, w4])
    W_2 = np.array([w5, w6, w7, w8])
    W_3 = np.array([w9, w10, w11, w12])
    W_4 = np.array([w13, w14, w15])
    
    N = 50
    dt = 0.001
    sigma = 0.005
    theta1 = 0.
    theta2 = 0.25
    f1 = 10
    f2 = 100
    A1 = 1.
    A2 = 0.3
    ts = np.arange(0, N)*dt*2.
    
    nb_steps = 26
    nb_t = int(N/nb_steps)
    memory_size = 15
    theta1_prior_mean = 0.
    theta2_prior_mean = 0.
    theta1_prior_sigma = 100.
    theta2_prior_sigma = 100.
    f2_prior_mean = 0.5
    f2_prior_sigma = 100.
    y_pred = np.zeros(nb_t)
    
    total_error = 0.0
    prev_error = 0.0

    params_list = np.array([
        [0., 0.25, 100, math.pi, -0.87, 100]
    ])
    k = 0

    file_id = len(glob.glob('./data/csv/*.csv'))
    f = open('./data/csv/test_'+str(file_id)+'.csv', 'w')
    writer = csv.writer(f, lineterminator='\n')

    theta1 = params_list[k, 0]
    theta2 = params_list[k, 1]
    f2 = params_list[k, 2]
    outs = create_data(ts, f1, f2, theta1, theta2, A1, A2, sigma=sigma)

    ### Before chaning parameters
    # Initial prediction
    _ys = create_data(ts, f1, f2, theta1, theta2, A1, A2, sigma=sigma)
    _ts = ts + np.random.randn(N)*0.001
    y_pred = create_data(_ts, f1, f2_prior_mean, theta1_prior_mean, theta2_prior_mean, A1, A2, sigma=sigma)
    error = np.mean(np.power(_ys - np.array(y_pred), 2))
    total_error = total_error + error
    stan_data = {
        "N": len(_ts),
        "t": _ts,
        "y": _ys,
        "sigma": sigma,
        "theta1_prior_mean": theta1_prior_mean,
        "theta2_prior_mean": theta2_prior_mean,
        "theta1_prior_sigma": theta1_prior_sigma,
        "theta2_prior_sigma": theta2_prior_sigma,
        "f2_prior_mean": f2_prior_mean,
        "f2_prior_sigma": f2_prior_sigma,
        "f1": f1,
        "A1": A1,
        "A2": A2
    }
    fit = sm.sampling(data=stan_data, iter=4000, chains=1)
    alpha_inputs = np.array([error, float(np.std(theta1_prior_sigma)/100.), float(np.std(theta2_prior_sigma)/100.), float(np.std(f2_prior_sigma)/100.)])
    h1 = np.tanh(np.dot(alpha_inputs, W_1))
    h2 = np.tanh(np.dot(alpha_inputs, W_2))
    h3 = np.tanh(np.dot(alpha_inputs, W_3))
    alpha = sigmoid(np.dot(np.array([h1, h2, h3]), W_4))
    theta1_prior_mean = (1-alpha)*theta1_prior_mean + alpha*np.mean(fit.extract(permuted=True)['theta1'])
    theta2_prior_mean = alpha*theta2_prior_mean + (1-alpha)*np.mean(fit.extract(permuted=True)['theta2'])
    f2_prior_mean = alpha*f2_prior_mean + (1-alpha)*np.mean(fit.extract(permuted=True)['f2'])
    prev_error = error
    writer.writerow([1, alpha, error, theta1_prior_sigma, theta2_prior_sigma, f2_prior_sigma, theta1_prior_mean, theta2_prior_mean, f2_prior_mean, theta1, theta2, f2])

    # Second prediction
    _ys = create_data(ts, f1, f2, theta1, theta2, A1, A2, sigma=sigma)
    _ts = ts + np.random.randn(N)*0.001
    y_pred = create_data(_ts, f1, f2_prior_mean, theta1_prior_mean, theta2_prior_mean, A1, A2, sigma=sigma)
    error = np.mean(np.power(_ys - np.array(y_pred), 2))
    total_error = total_error + error
    stan_data = {
        "N": len(_ts),
        "t": _ts,
        "y": _ys,
        "sigma": sigma,
        "theta1_prior_mean": theta1_prior_mean,
        "theta2_prior_mean": theta2_prior_mean,
        "theta1_prior_sigma": theta1_prior_sigma,
        "theta2_prior_sigma": theta2_prior_sigma,
        "f2_prior_mean": f2_prior_mean,
        "f2_prior_sigma": f2_prior_sigma,
        "f1": f1,
        "A1": A1,
        "A2": A2
    }
    fit = sm.sampling(data=stan_data, iter=4000, chains=1)
    alpha_inputs = np.array([error, float(np.std(theta1_prior_sigma)/100.), float(np.std(theta2_prior_sigma)/100.), float(np.std(f2_prior_sigma)/100.)])
    h1 = np.tanh(np.dot(alpha_inputs, W_1))
    h2 = np.tanh(np.dot(alpha_inputs, W_2))
    h3 = np.tanh(np.dot(alpha_inputs, W_3))
    alpha = sigmoid(np.dot(np.array([h1, h2, h3]), W_4))
    theta1_prior_mean = (1-alpha)*theta1_prior_mean + alpha*np.mean(fit.extract(permuted=True)['theta1'])
    theta2_prior_mean = alpha*theta2_prior_mean + (1-alpha)*np.mean(fit.extract(permuted=True)['theta2'])
    f2_prior_mean = alpha*f2_prior_mean + (1-alpha)*np.mean(fit.extract(permuted=True)['f2'])
    diff_error = error - prev_error
    theta1_prior_sigma = (1-alpha)*theta1_prior_sigma + alpha*theta1_prior_sigma*np.exp(w_sigma*diff_error)
    theta2_prior_sigma = alpha*theta2_prior_sigma + (1-alpha)*theta2_prior_sigma*np.exp(w_sigma*diff_error)
    f2_prior_sigma = alpha*f2_prior_sigma + (1-alpha)*f2_prior_sigma*np.exp(w_sigma*diff_error)
    prev_error = error
    writer.writerow([2, alpha, error, theta1_prior_sigma, theta2_prior_sigma, f2_prior_sigma, theta1_prior_mean, theta2_prior_mean, f2_prior_mean, theta1, theta2, f2])

    # 3rd prediction
    _ys = create_data(ts, f1, f2, theta1, theta2, A1, A2, sigma=sigma)
    _ts = ts + np.random.randn(N)*0.001
    y_pred = create_data(_ts, f1, f2_prior_mean, theta1_prior_mean, theta2_prior_mean, A1, A2, sigma=sigma)
    error = np.mean(np.power(_ys - np.array(y_pred), 2))
    total_error = total_error + error
    stan_data = {
        "N": len(_ts),
        "t": _ts,
        "y": _ys,
        "sigma": sigma,
        "theta1_prior_mean": theta1_prior_mean,
        "theta2_prior_mean": theta2_prior_mean,
        "theta1_prior_sigma": theta1_prior_sigma,
        "theta2_prior_sigma": theta2_prior_sigma,
        "f2_prior_mean": f2_prior_mean,
        "f2_prior_sigma": f2_prior_sigma,
        "f1": f1,
        "A1": A1,
        "A2": A2
    }
    fit = sm.sampling(data=stan_data, iter=4000, chains=1)
    alpha_inputs = np.array([error, float(np.std(theta1_prior_sigma)/100.), float(np.std(theta2_prior_sigma)/100.), float(np.std(f2_prior_sigma)/100.)])
    h1 = np.tanh(np.dot(alpha_inputs, W_1))
    h2 = np.tanh(np.dot(alpha_inputs, W_2))
    h3 = np.tanh(np.dot(alpha_inputs, W_3))
    alpha = sigmoid(np.dot(np.array([h1, h2, h3]), W_4))
    theta1_prior_mean = (1-alpha)*theta1_prior_mean + alpha*np.mean(fit.extract(permuted=True)['theta1'])
    theta2_prior_mean = alpha*theta2_prior_mean + (1-alpha)*np.mean(fit.extract(permuted=True)['theta2'])
    f2_prior_mean = alpha*f2_prior_mean + (1-alpha)*np.mean(fit.extract(permuted=True)['f2'])
    diff_error = error - prev_error
    theta1_prior_sigma = (1-alpha)*theta1_prior_sigma + alpha*theta1_prior_sigma*np.exp(w_sigma*diff_error)
    theta2_prior_sigma = alpha*theta2_prior_sigma + (1-alpha)*theta2_prior_sigma*np.exp(w_sigma*diff_error)
    f2_prior_sigma = alpha*f2_prior_sigma + (1-alpha)*f2_prior_sigma*np.exp(w_sigma*diff_error)
    prev_error = error
    writer.writerow([3, alpha, error, theta1_prior_sigma, theta2_prior_sigma, f2_prior_sigma, theta1_prior_mean, theta2_prior_mean, f2_prior_mean, theta1, theta2, f2])

    ### After chaning parameters
    theta1 = params_list[k, 3]
    theta2 = params_list[k, 4]
    f2 = params_list[k, 5]
    outs = create_data(ts, f1, f2, theta1, theta2, A1, A2, sigma=sigma)
    for i in range(10):
        _ys = create_data(ts, f1, f2, theta1, theta2, A1, A2, sigma=sigma)
        _ts = ts + np.random.randn(N)*0.001
        y_pred = create_data(_ts, f1, f2_prior_mean, theta1_prior_mean, theta2_prior_mean, A1, A2, sigma=sigma)
        error = np.mean(np.power(_ys - np.array(y_pred), 2))
        total_error = total_error + error
        stan_data = {
            "N": len(_ts),
            "t": _ts,
            "y": _ys,
            "sigma": sigma,
            "theta1_prior_mean": theta1_prior_mean,
            "theta2_prior_mean": theta2_prior_mean,
            "theta1_prior_sigma": theta1_prior_sigma,
            "theta2_prior_sigma": theta2_prior_sigma,
            "f2_prior_mean": f2_prior_mean,
            "f2_prior_sigma": f2_prior_sigma,
            "f1": f1,
            "A1": A1,
            "A2": A2
        }
        fit = sm.sampling(data=stan_data, iter=4000, chains=1)
        alpha_inputs = np.array([error, float(np.std(theta1_prior_sigma)/100.), float(np.std(theta2_prior_sigma)/100.), float(np.std(f2_prior_sigma)/100.)])
        h1 = np.tanh(np.dot(alpha_inputs, W_1))
        h2 = np.tanh(np.dot(alpha_inputs, W_2))
        h3 = np.tanh(np.dot(alpha_inputs, W_3))
        alpha = sigmoid(np.dot(np.array([h1, h2, h3]), W_4))
        theta1_prior_mean = (1-alpha)*theta1_prior_mean + alpha*np.mean(fit.extract(permuted=True)['theta1'])
        theta2_prior_mean = alpha*theta2_prior_mean + (1-alpha)*np.mean(fit.extract(permuted=True)['theta2'])
        f2_prior_mean = alpha*f2_prior_mean + (1-alpha)*np.mean(fit.extract(permuted=True)['f2'])
        diff_error = error - prev_error
        theta1_prior_sigma = (1-alpha)*theta1_prior_sigma + alpha*theta1_prior_sigma*np.exp(w_sigma*diff_error)
        theta2_prior_sigma = alpha*theta2_prior_sigma + (1-alpha)*theta2_prior_sigma*np.exp(w_sigma*diff_error)
        f2_prior_sigma = alpha*f2_prior_sigma + (1-alpha)*f2_prior_sigma*np.exp(w_sigma*diff_error)
        prev_error = error
        writer.writerow([i+4, alpha, error, theta1_prior_sigma, theta2_prior_sigma, f2_prior_sigma, theta1_prior_mean, theta2_prior_mean, f2_prior_mean, theta1, theta2, f2])
    
    f.close()
    return total_error


main()