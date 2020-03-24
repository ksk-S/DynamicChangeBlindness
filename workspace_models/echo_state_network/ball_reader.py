import csv
import numpy as np
#import numpy.random as rnd
import random
from PIL import Image


def getDataYContinuous(folder, n, offset=0):
	#how to coarse grain the input image (eg 4x4 grid)
	cells_by_axis = 4
	#nx2x4x4 table: n images, green and blue channel, 4x4 cells
	pixel_data = np.zeros( (n,2,cells_by_axis,cells_by_axis) )
	#number of main balls
	main_balls_count = 3
	#nx2
	balls_data = np.zeros( (n,main_balls_count*2+1) )

	for i in range(n):
		image_name = str(i+offset).zfill(5) + ".png"
		#print "image ", image_name
		image = Image.open(folder + "images/"+image_name)
		pixels = image.convert('RGBA')

		#h,w
		image_size = image.size
		#print "image_size ", image_size 
		cell_size = tuple(ti/cells_by_axis for ti in image_size)
		cell_pixels = cell_size[0]*cell_size[1]
		#print "cell_pixels ", cell_pixels

		for row in range(image_size[0]):
			cell_row = int(row/cell_size[0])
			for col in range(image_size[1]):
				cell_col = int(col/cell_size[1])
				#print "cell_row ", cell_row, " cell_col ", cell_col
				r, g, b, a = pixels.getpixel((col, row))
				#green
				pixel_data[i, 0, cell_row, cell_col] += g
				#blue
				pixel_data[i, 1, cell_row, cell_col] += b
				
		#divide by cell size
		for x in range(cells_by_axis):
			for y in range(cells_by_axis):
				pixel_data[i, 0, x, y] = pixel_data[i, 0, x, y]/cell_pixels
				pixel_data[i, 1, x, y] = pixel_data[i, 1, x, y]/cell_pixels

	

	#direction of blue ball
	with open(folder + 'main_ball.csv') as csvfile:
		reader = csv.DictReader(csvfile)
		i = 0
		for row in reader:
			for ball in range(main_balls_count):
				index = ball*2
				balls_data[i,index] = float(row["x" + str(ball)])
				balls_data[i,index+1] = float(row["y" + str(ball)])
			i = i+1
			if i>=n:
				break

	with open(folder + 'secondary_ball.csv') as csvfile:
		reader = csv.DictReader(csvfile)
		i = 0
		for row in reader:
			#presence of green ball
			balls_data[i,main_balls_count*2] = int(row["visible"])
			i = i+1
			if i>=n:
				break
	csvfile.close()

	return pixel_data, balls_data

def getDataIgnoreBallsChange(folder, n, offset=0):
	#how to coarse grain the input image (eg 4x4 grid)
	cells_by_axis = 4
	#nxcolorsx4x4 table: n images, color channels, 4x4 cells
	pixel_data = np.zeros( (n, 3, cells_by_axis, cells_by_axis) )
	#number of main balls not to ignore
	main_balls_count = 2
	#n * (x,y)
	balls_data = np.zeros( (n,main_balls_count*2+1) )

	for i in range(n):
		image_name = str(i+offset).zfill(5) + ".png"
		#print "image ", image_name
		image = Image.open(folder + "images/"+image_name)
		pixels = image.convert('RGBA')

		#h,w
		image_size = image.size
		#print "image_size ", image_size 
		cell_size = tuple(ti/cells_by_axis for ti in image_size)
		cell_pixels = cell_size[0]*cell_size[1]
		#print "cell_pixels ", cell_pixels

		for row in range(image_size[0]):
			cell_row = int(row/cell_size[0])
			for col in range(image_size[1]):
				cell_col = int(col/cell_size[1])
				#print "cell_row ", cell_row, " cell_col ", cell_col
				r, g, b, a = pixels.getpixel((col, row))
				#red
				pixel_data[i, 0, cell_row, cell_col] += r
				#green
				pixel_data[i, 1, cell_row, cell_col] += g
				#blue
				pixel_data[i, 2, cell_row, cell_col] += b

				
		#divide by cell size
		for x in range(cells_by_axis):
			for y in range(cells_by_axis):
				pixel_data[i, 0, x, y] = pixel_data[i, 0, x, y]/cell_pixels
				pixel_data[i, 1, x, y] = pixel_data[i, 1, x, y]/cell_pixels
				pixel_data[i, 2, x, y] = pixel_data[i, 2, x, y]/cell_pixels

	#Y
	with open(folder + 'main_ball.csv') as csvfile:
		reader = csv.DictReader(csvfile)
		i = 0
		for row in reader:
			for ball in range(main_balls_count):
				index = ball*2
				balls_data[i,index] = float(row["x" + str(ball)])
				balls_data[i,index+1] = float(row["y" + str(ball)])
			i = i+1
			if i>=n:
				break

	with open(folder + 'secondary_ball.csv') as csvfile:
		reader = csv.DictReader(csvfile)
		i = 0
		for row in reader:
			#presence of secondary object
			balls_data[i,main_balls_count*2] = int(row["visible"])
			i = i+1
			if i>=n:
				break
	csvfile.close()

	return pixel_data, balls_data

def getDataChangeBall(folder, n, offset=0):
	#how to coarse grain the input image (eg 4x4 grid)
	cells_by_axis = 4
	input_size = 2*4*4
	#nx2x4x4 table: n images, red and blue channels, 4x4 cells
	pixel_data = np.zeros( (n,2,cells_by_axis,cells_by_axis) )
	switch_report = np.zeros((n))
	#number of main balls not to ignore
	main_balls_count = 1

	#expected output
	#n * (x,y,x2,y2.,...)
	output_size1 = main_balls_count*2
	balls_data = np.zeros( (n, output_size1) )
	#switch will store the color of the ball
	output_size2 = 1
	switch_data = np.zeros( (n, output_size2) )

	for i in range(n):
		image_name = str(i+offset).zfill(5) + ".png"
		image = Image.open(folder + "images/"+image_name)
		pixels = image.convert('RGBA')

		#h,w
		image_size = image.size
		cell_size = tuple(ti/cells_by_axis for ti in image_size)
		cell_pixels = cell_size[0]*cell_size[1]

		for row in range(image_size[0]):
			cell_row = int(row/cell_size[0])
			for col in range(image_size[1]):
				cell_col = int(col/cell_size[1])
				#print "cell_row ", cell_row, " cell_col ", cell_col
				r, g, b, a = pixels.getpixel((col, row))
				#red
				pixel_data[i, 0, cell_row, cell_col] += r
				#blue
				pixel_data[i, 1, cell_row, cell_col] += b
				
		#divide by cell size
		for x in range(cells_by_axis):
			for y in range(cells_by_axis):
				pixel_data[i, 0, x, y] = pixel_data[i, 0, x, y]/cell_pixels
				pixel_data[i, 1, x, y] = pixel_data[i, 1, x, y]/cell_pixels

	ballIndex = 0

	#x and y
	with open(folder + 'balls.csv') as csvfile:
		reader = csv.DictReader(csvfile)
		i = 0
		for row in reader:		
			#print(row)
			index = ballIndex*2
			balls_data[i,0] = float(row["x" + str(ballIndex)])
			balls_data[i,1] = float(row["y" + str(ballIndex)])
			switch_data[i,0] = int(row["color"])

			i = i+1
			if i>=n:
				break

	csvfile.close()

	return input_size, output_size1, output_size2, pixel_data, balls_data, switch_data

def getDataSwitchBalls(folder, n, offset=0):
	#how to coarse grain the input image (eg 4x4 grid)
	cells_by_axis = 4
	input_size = 2*4*4
	#nx2x4x4 table: n images, red and blue channels, 4x4 cells
	pixel_data = np.zeros( (n,2,cells_by_axis,cells_by_axis) )
	switch_report = np.zeros((n))
	#number of main balls not to ignore
	main_balls_count = 1

	#expected output
	#n * (x,y,x2,y2.,...)
	output_size1 = main_balls_count*2
	balls_data = np.zeros( (n, output_size1) )
	#switch report
	output_size2 = 1
	switch_data = np.zeros( (n, output_size2) )

	for i in range(n):
		image_name = str(i+offset).zfill(5) + ".png"
		image = Image.open(folder + "images/"+image_name)
		pixels = image.convert('RGBA')

		#h,w
		image_size = image.size
		cell_size = tuple(ti/cells_by_axis for ti in image_size)
		cell_pixels = cell_size[0]*cell_size[1]

		for row in range(image_size[0]):
			cell_row = int(row/cell_size[0])
			for col in range(image_size[1]):
				cell_col = int(col/cell_size[1])
				#print "cell_row ", cell_row, " cell_col ", cell_col
				r, g, b, a = pixels.getpixel((col, row))
				#red
				pixel_data[i, 0, cell_row, cell_col] += r
				#blue
				pixel_data[i, 1, cell_row, cell_col] += b
				
		#divide by cell size
		for x in range(cells_by_axis):
			for y in range(cells_by_axis):
				pixel_data[i, 0, x, y] = pixel_data[i, 0, x, y]/cell_pixels
				pixel_data[i, 1, x, y] = pixel_data[i, 1, x, y]/cell_pixels

	#number of switches
	#meanSwitchTime = 100
	#switches = int(n/meanSwitchTime)
	ballIndex = 0
	#switchTiming = rnd.rand(switches)*meanSwitchTime
	switchTiming = random.randint(100,150)
	#x and y
	with open(folder + 'balls.csv') as csvfile:
		reader = csv.DictReader(csvfile)
		i = 0
		for row in reader:
			
			if(switchTiming<0):
				switchTiming = random.randint(100,150)
				#print("switch at ", i, "ballIndex ", ballIndex)

				if (ballIndex == 0):
					ballIndex = 1
				else :
					ballIndex = 0

			#report change 1 time step later
			# if(switchTiming<0 and i+1<n):
			# 	switch_data[i+1,0] = 1

			#report which ball should be tracked
			if(i+1<n):
				switch_data[i+1,0] = ballIndex*2 - 1 #-1 or 1 
			
			index = ballIndex*2
			balls_data[i,0] = float(row["x" + str(ballIndex)])
			balls_data[i,1] = float(row["y" + str(ballIndex)])
			
			i = i+1
			switchTiming = switchTiming -1;
			if i>=n:
				break

	csvfile.close()

	return input_size, output_size1, output_size2, pixel_data, balls_data, switch_data

def getDataIgnoreBalls(folder, n, offset=0):
	#how to coarse grain the input image (eg 4x4 grid)
	cells_by_axis = 4
	#nx2x4x4 table: n images, green and blue channels, 4x4 cells
	pixel_data = np.zeros( (n,2,cells_by_axis,cells_by_axis) )
	#number of main balls not to ignore
	main_balls_count = 2
	#n * (x,y)
	balls_data = np.zeros( (n,main_balls_count*2+1) )

	for i in range(n):
		image_name = str(i+offset).zfill(5) + ".png"
		#print "image ", image_name
		image = Image.open(folder + "images/"+image_name)
		pixels = image.convert('RGBA')

		#h,w
		image_size = image.size
		#print "image_size ", image_size 
		cell_size = tuple(ti/cells_by_axis for ti in image_size)
		cell_pixels = cell_size[0]*cell_size[1]
		#print "cell_pixels ", cell_pixels

		for row in range(image_size[0]):
			cell_row = int(row/cell_size[0])
			for col in range(image_size[1]):
				cell_col = int(col/cell_size[1])
				#print "cell_row ", cell_row, " cell_col ", cell_col
				r, g, b, a = pixels.getpixel((col, row))
				#green
				pixel_data[i, 0, cell_row, cell_col] += g
				#blue
				pixel_data[i, 1, cell_row, cell_col] += b

				
		#divide by cell size
		for x in range(cells_by_axis):
			for y in range(cells_by_axis):
				pixel_data[i, 0, x, y] = pixel_data[i, 0, x, y]/cell_pixels
				pixel_data[i, 1, x, y] = pixel_data[i, 1, x, y]/cell_pixels

	#Y
	with open(folder + 'main_ball.csv') as csvfile:
		reader = csv.DictReader(csvfile)
		i = 0
		for row in reader:
			for ball in range(main_balls_count):
				index = ball*2
				balls_data[i,index] = float(row["x" + str(ball)])
				balls_data[i,index+1] = float(row["y" + str(ball)])
			i = i+1
			if i>=n:
				break
	csvfile.close()

	with open(folder + 'secondary_ball.csv') as csvfile:
		reader = csv.DictReader(csvfile)
		i = 0
		for row in reader:
			#presence of secondary object
			balls_data[i,main_balls_count*2] = int(row["visible"])
			i = i+1
			if i>=n:
				break
	csvfile.close()

	return pixel_data, balls_data


def getDataColorChange(folder, n, offset=0):
	#how to coarse grain the input image (eg 4x4 grid)
	cells_by_axis = 4
	#nx3x4x4 table: n images, green, blue and red channels, 4x4 cells
	pixel_data = np.zeros( (n,3,cells_by_axis,cells_by_axis) )
	#number of main balls
	main_balls_count = 1
	#n * (x,y)
	balls_data = np.zeros( (n,main_balls_count*2+1) )

	for i in range(n):
		image_name = str(i+offset).zfill(5) + ".png"
		#print "image ", image_name
		image = Image.open(folder + "images/"+image_name)
		pixels = image.convert('RGBA')

		#h,w
		image_size = image.size
		#print "image_size ", image_size 
		cell_size = tuple(ti/cells_by_axis for ti in image_size)
		cell_pixels = cell_size[0]*cell_size[1]
		#print "cell_pixels ", cell_pixels

		for row in range(image_size[0]):
			cell_row = int(row/cell_size[0])
			for col in range(image_size[1]):
				cell_col = int(col/cell_size[1])
				#print "cell_row ", cell_row, " cell_col ", cell_col
				r, g, b, a = pixels.getpixel((col, row))
				#red
				pixel_data[i, 0, cell_row, cell_col] += r
				#green
				pixel_data[i, 1, cell_row, cell_col] += g
				#blue
				pixel_data[i, 2, cell_row, cell_col] += b

				
		#divide by cell size
		for x in range(cells_by_axis):
			for y in range(cells_by_axis):
				pixel_data[i, 0, x, y] = pixel_data[i, 0, x, y]/cell_pixels
				pixel_data[i, 1, x, y] = pixel_data[i, 1, x, y]/cell_pixels
				pixel_data[i, 2, x, y] = pixel_data[i, 2, x, y]/cell_pixels

	
	with open(folder + 'main_ball.csv') as csvfile:
		reader = csv.DictReader(csvfile)
		i = 0
		for row in reader:
			for ball in range(main_balls_count):
				index = ball*2
				balls_data[i,index] = float(row["x" + str(ball)])
				balls_data[i,index+1] = float(row["y" + str(ball)])
			i = i+1
			if i>=n:
				break

	with open(folder + 'secondary_ball.csv') as csvfile:
		reader = csv.DictReader(csvfile)
		i = 0
		for row in reader:
			#presence of green ball
			balls_data[i,main_balls_count*2] = int(row["visible"])
			i = i+1
			if i>=n:
				break
	csvfile.close()

	return pixel_data, balls_data

#n<12,000
#return nx2x4x4 pixel color data, nx2 [main ball vertical direction, second ball presence]
#duplicate: how many times to repeat each frame in the returned dataset (default 1, no repetition)
def getDataYBinary(folder, n, duplicate=1, offset=0):
	#how to coarse grain the input image (eg 4x4 grid)
	cells_by_axis = 4
	#nx2x4x4 table: n images, green and blue channel, 4x4 cells
	pixel_data = np.zeros( (n*duplicate,2,cells_by_axis,cells_by_axis) )
	#nx2
	balls_data = np.zeros( (n*duplicate,2) )

	for i in range(n):
		image_name = str(i+offset).zfill(5) + ".png"
		#print "image ", image_name
		image = Image.open(folder + "images/"+image_name)
		pixels = image.convert('RGBA')

		#h,w
		image_size = image.size
		#print "image_size ", image_size 
		cell_size = tuple(ti/cells_by_axis for ti in image_size)
		cell_pixels = cell_size[0]*cell_size[1]
		#print "cell_pixels ", cell_pixels

		for row in range(image_size[0]):
			cell_row = int(row/cell_size[0])
			for col in range(image_size[1]):
				cell_col = int(col/cell_size[1])
				#print "cell_row ", cell_row, " cell_col ", cell_col
				r, g, b, a = pixels.getpixel((col, row))
				for rep in range(duplicate):
					index = i*duplicate + rep
					#green
					pixel_data[index, 0, cell_row, cell_col] += g
					#blue
					pixel_data[index, 1, cell_row, cell_col] += b
				
		#divide by cell size
		for x in range(cells_by_axis):
			for y in range(cells_by_axis):
				pixel_data[i, 0, x, y] = pixel_data[i, 0, x, y]/cell_pixels
				pixel_data[i, 1, x, y] = pixel_data[i, 1, x, y]/cell_pixels

	

	#direction of blue ball
	with open(folder + 'main_ball.csv') as csvfile:
		reader = csv.DictReader(csvfile)
		i = 0
		for row in reader:
			for rep in range(duplicate):
				index = i*duplicate + rep
				balls_data[index,0] = int(row["direction"])
			i = i+1
			if i>=n:
				break

	with open(folder + 'secondary_ball.csv') as csvfile:
		reader = csv.DictReader(csvfile)
		i = 0
		for row in reader:
			for rep in range(duplicate):
				index = i*duplicate + rep
				#presence of green ball
				balls_data[index,1] = int(row["visible"])
			i = i+1
			if i>=n:
				break
	csvfile.close()

	return pixel_data, balls_data

