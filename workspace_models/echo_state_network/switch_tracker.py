import ball_reader
import csv
from StackedESN import StackedESN
import matplotlib.pyplot as plt
import numpy as np
import os


def train_tracker(ball_data_folder, plot_path, networks_count, n_reservoir,
						train_cycles=20, test_cycles=10, 
						alpha=0.8, resSize=100, plot_show=True):

	plot_name='tracker'
	if not os.path.exists(plot_path+"networks"):
		os.makedirs(plot_path+"networks")

	#number of time steps
	total_n = train_cycles+test_cycles
	
	# Get Data
	#return nx2x4x4 pixel color data
	# x, y
	# switched or not
	input_size, output_size1, output_size2, pix_data, output1, output2 = ball_reader.getDataSwitchBalls(ball_data_folder, train_cycles+test_cycles)

	# reshape the data to remove dimensions
	input_data = pix_data.reshape((total_n,input_size))

	data_train, y1_train, y2_train = input_data[:train_cycles], output1[:train_cycles], output2[:train_cycles]
	data_test, y1_test, y2_test = input_data[train_cycles:], output1[train_cycles:], output2[train_cycles:]

	#save input data
	header = "x0,y0"
	np.savetxt(plot_path + "input.csv", y1_test, delimiter=",", header=header, fmt='%3f', comments='')

	for i in range(networks_count):
		esn = StackedESN(n_input=input_size, n_output1=output_size1, n_output2=output_size2, n_reservoir=n_reservoir, leakingRate=0.2, spectralRadius=0.9, regressionParameters=[1e-2])

		esn.fit(inputData=data_train, outputData1=y1_train, outputData2=y2_train, transientTime=0, verbose=1)
		y1_test_pred, y2_test_pred = esn.predict(inputData=data_test, outputData1=y1_test, outputData2=y2_test, verbose=1) #continuation can be set to false
		
		np.savetxt(plot_path + "track_" + str(i) + ".csv", y1_test_pred, delimiter=",", header=header, fmt='%3f', comments='')
		#save network
		esn.save(path=plot_path + "networks/" + str(i))

		# Prediction Plot
		plt.figure('Coordinates Estimation', figsize=(14,7)).clear()
		plt.plot(y1_test[:,0], color='red', linewidth=5, label='Target Value')
		plt.plot(y1_test_pred[:,0], color='blue', linestyle="--", linewidth=1, label='Estimation')
		plt.legend()
		figname = plot_path + plot_name + "_tracking_" + str(i)  + ".png"
		plt.savefig(figname)
		print('\t[+]Plot saved in', figname)
		#if plot_show:
		#	plt.show()

		plt.figure('Switch Estimation', figsize=(14,7)).clear()
		#plt.yscale('log')
		plt.plot(y2_test[:,0], color='red', linewidth=5, label='Target Value')
		plt.plot(y2_test_pred[:,0], color='blue', linestyle="--", linewidth=1, label='Estimation')
		plt.legend()
		figname = plot_path + plot_name + "_switch_" + str(i) + ".png"
		plt.savefig(figname)
		print('\t[+]Plot saved in', figname)
		#if plot_show:
		#	plt.show()

def evaluate_tracker(plot_path, ball_data_folder, esn_path, steps, offset):
	plot_name='tracker'
	esn = StackedESN.load(esn_path)

	input_size, output_size1, output_size2, pix_data, output1, output2 = ball_reader.getDataChangeBall(ball_data_folder, steps, offset=offset)
	# reshape the data to remove dimensions
	input_data = pix_data.reshape((steps,input_size))
	dummy_switch_info = np.zeros((output2.shape))

	y_test_pred, _ = esn.predict(inputData=input_data, outputData1=output1, outputData2=dummy_switch_info, verbose=1, continuation=False)

	plt.figure('Change Estimation', figsize=(14,7)).clear()
	plt.plot(output1[:,0], color='red', linewidth=5, label='Target Value')
	plt.plot(y_test_pred[:,0], color='blue', linestyle="--", linewidth=1, label="Estimation")
	plt.plot(output2[:,0]*100, color='black', linewidth=2, label="Ball color")

	plt.legend()
	figname = plot_path + plot_name + "_change.png"
	plt.savefig(figname)

	#save in csv file
	header = "x0,y0"
	np.savetxt(plot_path + "evaluation_data.csv", output1, delimiter=",", header=header, fmt='%3f', comments='')

	header = "x0,y0"
	np.savetxt(plot_path + "evaluated_data.csv", y_test_pred, delimiter=",", header=header, fmt='%3f', comments='')
	plt.show()


if __name__=='__main__':
	ball_data_folder = "../data/switch_balls/" #"../data/changing_ball/"
	data_path = 'cb_results/switch_balls/'
	main_balls_count = 2
	train_cycles = 800
	test_cycles = 200
	n_reservoir = 50

	#train_tracker(ball_data_folder=ball_data_folder, plot_path=data_path, networks_count=1, train_cycles = train_cycles, test_cycles=test_cycles, n_reservoir=n_reservoir)
	evaluate_tracker(plot_path=data_path, ball_data_folder="../data/changing_ball/", esn_path=data_path+"networks/0", steps=500, offset=0)

