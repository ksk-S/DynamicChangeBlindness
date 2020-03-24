import numpy as np
import pandas as pd
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

import warnings
warnings.filterwarnings('ignore')

import optuna


def create_data(f1, f2, A1, A2, sigma=0.02):
    outs = []
    ts = 4500
    theta1 = 0.7
    theta2 = 0.4
    for t in range(ts):
        if t == 500:
            theta1 = 1.0
            theta2 = -1.0
        elif t == 1000:
            theta1 = 1.0
            theta2 = 0.0
        elif t == 1500:
            theta1 = 0.3
            theta2 = 0.0
        elif t == 2000:
            theta1 = -0.5
            theta2 = 0.9
        elif t == 2500:
            theta1 = 0.3
            theta2 = 0.0
        elif t == 3000:
            theta1 = 0.9
            theta2 = 0.0
        elif t == 3500:
            theta1 = 0.9
            theta2 = -1.3
        elif t == 4000:
            theta1 = -1.2
            theta2 = 1.3
        n_f1 = np.random.normal(0.0, 0.05)
        n_f2 = np.random.normal(0.0, 0.05)
        val = A1*math.sin(f1*t+theta1+n_f1) + A2*math.sin(f2*t+theta2+n_f2) + np.random.normal(0.0, sigma)
        outs.append(val)
    return np.array(outs)


def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))


def relu(x):
    if x > 0:
        return x
    else:
        return 0.


### EKF
def predict_phase(x_vec, P_mat, J_s=np.eye(2), dw=np.array([0.01, 0.1]), Q_t=np.ones((2,2))):
    # J_s: Jacobian
    x_hat = x_vec + dw
    P_hat = np.matmul(np.matmul(J_s,P_mat),J_s.T) + Q_t
    return x_hat, P_hat

def update_phase(obs, x_hat, P_hat, x_vec, P_mat, w_vec, R_t=np.eye(2)):
    y_error = obs - (np.sin(x_hat[0])+0.3*np.sin(x_hat[1]))
    w_err = np.array([np.tanh(y_error*w_vec[0]), np.tanh(y_error*w_vec[1]), np.tanh(y_error*w_vec[2]), np.tanh(y_error*w_vec[3])])
    alpha = sigmoid(np.dot(w_err, w_vec[4:]))
    J_o = np.array([np.cos(x_hat[0]), 0.3*np.cos(x_hat[1])]) # Jacobian
    S_t = np.matmul(np.matmul(J_o, P_hat), J_o.T) + R_t
    K_t = np.matmul(np.matmul(P_hat, J_o.T), np.linalg.inv(S_t)) # Kalman Gain
    K_t = K_t*np.array([alpha, 1.-alpha])
    new_x_vec = x_vec + K_t*y_error
    new_P_mat = np.matmul((np.eye(2) - np.matmul(K_t, J_o)), P_hat)
    return new_x_vec, new_P_mat, y_error


ys = create_data(f1=0.01, f2=0.1, A1=1.0, A2=0.3, sigma=0.05)


def objective(trial):
    w1 = trial.suggest_uniform('w1', -3., 3.)
    w2 = trial.suggest_uniform('w2', -3., 3.)
    w3 = trial.suggest_uniform('w3', -3., 3.)
    w4 = trial.suggest_uniform('w4', -3., 3.)
    w5 = trial.suggest_uniform('w5', -3., 3.)
    w6 = trial.suggest_uniform('w6', -3., 3.)
    w7 = trial.suggest_uniform('w7', -3., 3.)
    w8 = trial.suggest_uniform('w8', -3., 3.)
    W_1 = np.array([w1, w2, w3, w4, w5, w6, w7, w8])

    x_vec = np.array([0.0, 0.0])
    P_mat = np.eye(2)

    total_err = 0.0
    for _y in ys[1:]:
        x_hat, P_hat = predict_phase(x_vec, P_mat)
        new_x_vec, new_P_mat, y_error = update_phase(_y, x_hat, P_hat, x_vec, P_mat, W_1)
        x_vec = new_x_vec
        P_mat = new_P_mat
        total_err = total_err + np.sqrt(y_error*y_error)
    total_err = total_err/float(len(ys[1:]))
    return total_err


study = optuna.create_study()
study.optimize(objective, n_trials=1000)

print("params_{}".format(study.best_params))
print("value_{}".format(study.best_value))


# Current best value is 0.09289854978248953 with parameters: {'w1': 0.8868339845276003, 'w2': -2.4239527390853723, 'w3': 2.5663446991064536, 'w4': -1.835679959314501, 'w5': 2.668697875044799, 'w6': -0.578802425496894, 'w7': -2.3135794565999737, 'w8': -0.9460572459969298}