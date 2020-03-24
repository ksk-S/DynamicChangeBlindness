import ball_reader
import csv
#from easyesn import PredictionESN
from FeedbackPredictionESN import PredictionESN
#import FeedbackPredictionESN
import matplotlib.pyplot as plt
import numpy as np
import os


def train_tracker(ball_data_folder, plot_path, feedback, networks_count, 
						outputSize=3, train_cycles=300, test_cycles=100, 
						alpha=0.8, resSize=100, plot_show=True):

	#where is the prediction? add speed estimation, or make it predict the next frame, or values for the next frame (delay output by 1)
	#48
	inputSize=4*4*2
	plot_name='tracker'
	if not os.path.exists(plot_path+"networks"):
		os.makedirs(plot_path+"networks")

	#how many time steps each video frame is shown the the ESN
	total_n = (train_cycles+test_cycles)
	#output duplication: how many perceptrons to train for each output
	#n_perceptrons = 10

	# Get Data
	#return nx2x4x4 pixel color data, nx2 [main ball vertical direction, second ball presence]
	multi_data, y = ball_reader.getDataIgnoreBalls(ball_data_folder, train_cycles+test_cycles)
	
	# reshape the data
	data = multi_data.reshape((total_n,inputSize))

	data_train, y_train = data[:train_cycles], y[:train_cycles]
	data_test, y_test = data[train_cycles:], y[train_cycles:]
	# time_delay = 5
	# #shift input t+1 to get predictions
	# y_train = np.roll(y_train,-time_delay, axis=0)
	# y_train[train_cycles-time_delay:] = np.zeros((time_delay,outputSize))
	# y_test = np.roll(y_test,-time_delay, axis=0)
	# y_test[test_cycles-time_delay:] = np.zeros((time_delay,outputSize))


	#save original data
	header = ""
	n_balls = 2
	for ball in range(n_balls):
		header = header + "x" + str(ball) + ",y" + str(ball) + ","
	header = header + "presence"
	np.savetxt(plot_path + "track_original.csv", y_test, delimiter=",", header=header, fmt='%3f', comments='')

	for i in range(networks_count):
		esn = PredictionESN(n_input=inputSize, n_output=outputSize, n_reservoir=10, leakingRate=0.2, spectralRadius=0.9, regressionParameters=[1e-2], feedback = feedback)
		esn.fit(data_train, y_train, transientTime=0, verbose=1)
		y_test_pred = esn.predict(data_test, verbose=1)

		
		np.savetxt(plot_path + "track_" + str(i) + ".csv", y_test_pred, delimiter=",", header=header, fmt='%3f', comments='')
		#save network
		esn.save(path=plot_path + "networks/" + str(i))

		# Prediction Plot
		plt.figure('Direction Prediction', figsize=(14,7)).clear()
		#plt.yscale('log')
		plt.plot(y_test[:,0], color='red', linewidth=5, label='Target Value')
		plt.plot(y_test_pred[:,0], color='blue', linestyle="--", linewidth=1, label='ESN Tracking Prediction')
		plt.legend()
		figname = plot_path + plot_name + "_tracking_" + str(i)  + ".png"
		plt.savefig(figname)
		print('\t[+]Plot saved in', figname)
		#if plot_show:
		#	plt.show()

		plt.figure('Presence Prediction', figsize=(14,7)).clear()
		#plt.yscale('log')
		plt.plot(y_test[:,outputSize-1], color='red', linewidth=5, label='Target Value')
		plt.plot(y_test_pred[:,outputSize-1], color='blue', linestyle="--", linewidth=1, label='ESN Presence Prediction')
		plt.legend()
		figname = plot_path + plot_name + "_presence_" + str(i) + ".png"
		plt.savefig(figname)
		print('\t[+]Plot saved in', figname)
		#if plot_show:
		#	plt.show()

def evaluate_tracker(ball_data_folder, esn_path, steps, offset, output_size):
	esn = PredictionESN.load(esn_path)

	data_test, y_test = ball_reader.getDataSwitchBalls(ball_data_folder, steps, offset=offset)
	y_test_pred = esn.predict(data_test, verbose=1)
	plt.plot(y_test[:,0], color='red', linewidth=5, label='Target Value')
	plt.plot(y_test_pred[:,0], color='blue', linestyle="--", linewidth=1, label='ESN Direction Prediction')
	plt.legend()
	plt.show()

	plt.plot(y_test[:,output_size-1], color='red', linewidth=5, label='Target Value')
	plt.plot(y_test_pred[:,output_size-1], color='blue', linestyle="--", linewidth=1, label='ESN Presence Prediction')
	plt.legend()
	plt.show()


if __name__=='__main__':
	ball_data_folder = "../data/ignore_balls/" #"../data/changing_ball/"
	data_path = 'cb_results/ignore_balls/'
	main_balls_count = 2
	#test()
	train_tracker(ball_data_folder=ball_data_folder, plot_path=data_path, feedback=False, networks_count=1, outputSize = main_balls_count*2+1)
	#evaluate_tracker(ball_data_folder=ball_data_folder, esn_path=data_path+"networks/0", steps=1000, offset=0, output_size = main_balls_count*2+1)

