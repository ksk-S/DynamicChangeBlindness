import numpy as np
import pandas as pd
import math
import json
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

import warnings
warnings.filterwarnings('ignore')

import optuna


def create_data(f1, f2, A1, A2, sigma=0.02):
    outs = []
    ts = 1000
    theta1 = 1.4
    theta2 = 1.0
    for t in range(ts):
        # if t == 500:
        #     theta1 = 1.4
        #     theta2 = -0.5
        # elif t == 1500:
        #     theta1 = 0.7
        #     theta2 = 0.0
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
    return new_x_vec, new_P_mat, y_error, alpha, K_t

ys = create_data(f1=0.01, f2=0.1, A1=1.0, A2=0.3, sigma=0.05)

# w_dict = {'w1': 0.8654948627671226, 'w2': -1.7444762795695032, 'w3': -1.256158244213108, 'w4': 2.9877172040880846, 'w5': 0.7674940690302532, 'w6': -0.5751565428986629, 'w7': -2.1525316155059886, 'w8': -1.593668210140296}
w_dict = {'w1': 0.8868339845276003, 'w2': -2.4239527390853723, 'w3': 2.5663446991064536, 'w4': -1.835679959314501, 'w5': 2.668697875044799, 'w6': -0.578802425496894, 'w7': -2.3135794565999737, 'w8': -0.9460572459969298}
w1 = w_dict['w1']
w2 = w_dict['w2']
w3 = w_dict['w3']
w4 = w_dict['w4']
w5 = w_dict['w5']
w6 = w_dict['w6']
w7 = w_dict['w7']
w8 = w_dict['w8']
W_1 = np.array([w1, w2, w3, w4, w5, w6, w7, w8])

x_vec = np.array([0.0, 0.0])
P_mat = np.eye(2)

total_err = 0.0
alphas = []
y_errors = []
preds = []
k_gains = []
ttt = 1
for _y in ys[1:]:
    x_hat, P_hat = predict_phase(x_vec, P_mat)
    new_x_vec, new_P_mat, y_error, _alpha, K_t = update_phase(_y, x_hat, P_hat, x_vec, P_mat, W_1)
    x_vec = new_x_vec
    P_mat = new_P_mat
    total_err = total_err + np.sqrt(y_error*y_error)
    alphas.append(_alpha)
    y_errors.append(np.abs(y_error))
    preds.append(np.sin(x_vec[0])+0.3*np.sin(x_vec[1]))
    k_gains.append(K_t.tolist())
    # print(ttt, y_error, _alpha)
    ttt = ttt + 1
total_err = total_err/float(len(ys[1:]))

with open("./data/json/test_ekf_no_alpha5.json", "w") as f:
    out_dict = {
        "k_gain": k_gains,
        # "alphas": alphas,
        "y_errors": y_errors,
        "ys": ys[1:].tolist(),
        "preds": preds
    }
    json.dump(out_dict, f)

# print(alphas)
# print(y_errors)
# print(ys[1:].tolist())
# print(preds)