import numpy as np
import pandas as pd
import pystan
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

import warnings
warnings.filterwarnings('ignore')

import optuna


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

def objective(trial):
    w1 = trial.suggest_uniform('w1', -3., 3.)
    w2 = trial.suggest_uniform('w2', -3., 3.)
    w3 = trial.suggest_uniform('w3', -3., 3.)
    w4 = trial.suggest_uniform('w4', -3., 3.)
    w5 = trial.suggest_uniform('w5', -3., 3.)
    w6 = trial.suggest_uniform('w6', -3., 3.)
    w7 = trial.suggest_uniform('w7', -3., 3.)
    w8 = trial.suggest_uniform('w8', -3., 3.)
    w9 = trial.suggest_uniform('w9', -3., 3.)
    w10 = trial.suggest_uniform('w10', -3., 3.)
    w11 = trial.suggest_uniform('w11', -3., 3.)
    w12 = trial.suggest_uniform('w12', -3., 3.)
    w13 = trial.suggest_uniform('w13', -3., 3.)
    w14 = trial.suggest_uniform('w14', -3., 3.)
    w15 = trial.suggest_uniform('w15', -3., 3.)
    w_sigma = trial.suggest_uniform('w_sigma', -10., 10.)
    W_1 = np.array([w1, w2, w3, w4])
    W_2 = np.array([w5, w6, w7, w8])
    W_3 = np.array([w9, w10, w11, w12])
    W_4 = np.array([w13, w14, w15])
    
    N = 50
    dt = 0.001
    sigma = 0.01
    theta1 = 0.
    theta2 = 0.25
    f1 = 10
    f2 = 100
    A1 = 1.
    A2 = 0.3
    ts = np.arange(0, N)*dt*2.
    # outs = create_data(ts, f1, f2, theta1, theta2, A1, A2, sigma=sigma)
    
    # nb_steps = 100
    nb_steps = 26
    nb_t = int(N/nb_steps)
    memory_size = 15
    theta1_prior_mean = 0.
    theta2_prior_mean = 0.
    theta1_prior_sigma = 20.
    theta2_prior_sigma = 20.
    f2_prior_mean = 0.5
    f2_prior_sigma = 20.
    y_pred = np.zeros(nb_t)
    
    total_error = 0.0

    params_list = np.array([
        [0.44, 0.25, 100, math.pi, -0.87, 88],
        [-0.54, 0.0, 80, -0.54, math.pi, 80],
        [math.pi*0.25, 0.5, 100, -math.pi*0.25, 0.5, 100],
        [0.6, math.pi*0.5, 60, 0.6, math.pi*0.5, 100],
        [0.3, 0.75, 60, -0.7, 0.0, 100],
    ])

    for k in range(len(params_list)):
        prev_error = 0.0
        theta1 = params_list[k, 0]
        theta2 = params_list[k, 1]
        f2 = params_list[k, 2]
        # outs = create_data(ts, f1, f2, theta1, theta2, A1, A2, sigma=sigma)

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
        fit = sm.sampling(data=stan_data, iter=5000, chains=1)
        alpha_inputs = np.array([error, float(theta1_prior_sigma/100.), float(theta2_prior_sigma/100.), float(f2_prior_sigma/100.)])
        h1 = np.tanh(np.dot(alpha_inputs, W_1))
        h2 = np.tanh(np.dot(alpha_inputs, W_2))
        h3 = np.tanh(np.dot(alpha_inputs, W_3))
        alpha = sigmoid(np.dot(np.array([h1, h2, h3]), W_4))
        theta1_prior_mean = np.mean(fit.extract(permuted=True)['theta1'])
        theta2_prior_mean = np.mean(fit.extract(permuted=True)['theta2'])
        f2_prior_mean = np.mean(fit.extract(permuted=True)['f2'])
        prev_error = error
        print("1-1 error:", error)

        # Second prediction
        for i in range(3):
            _ys = create_data(ts, f1, f2, theta1, theta2, A1, A2, sigma=sigma)
            _ts = ts + np.random.randn(N)*0.0005
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
            fit = sm.sampling(data=stan_data, iter=5000, chains=1)
            alpha_inputs = np.array([error, float(theta1_prior_sigma/20.), float(theta2_prior_sigma/20.), float(f2_prior_sigma/20.)])
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
            print("1-"+str(i+2)+" error:", error)

        ### After chaning parameters
        theta1 = params_list[k, 3]
        theta2 = params_list[k, 4]
        f2 = params_list[k, 5]
        outs = create_data(ts, f1, f2, theta1, theta2, A1, A2, sigma=sigma)
        for i in range(5):
            _ys = create_data(ts, f1, f2, theta1, theta2, A1, A2, sigma=sigma)
            _ts = ts + np.random.randn(N)*0.0005
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
            fit = sm.sampling(data=stan_data, iter=5000, chains=1)
            alpha_inputs = np.array([error, float(theta1_prior_sigma/20.), float(theta2_prior_sigma/20.), float(f2_prior_sigma/20.)])
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
            print("2-"+str(i+1)+" error:", error)
    return total_error


study = optuna.create_study()
study.optimize(objective, n_trials=100)

print("params_{}".format(study.best_params))
print("value_{}".format(study.best_value))


## First result (8/23) 
'''
{'w1': 5.018517220060774, 'w2': -3.6255715455977686, 'w3': -3.724743569108436, 'w4': 4.687505381875587, 'w5': 6.3858834752867075, 'w6': 3.7940853161519374, 'w7': -5.769327360288206, 'w8': 5.598481284654969, 'w9': 9.166909160777632, 'w10': -3.9337044952088127, 'w11': 5.112566876509584, 'w12': 2.2696441148754314, 'w13': 3.971480589411998, 'w14': 8.547504462553892, 'w15': -4.873755732191103, 'w_sigma': 6.361999199168867}
'''
