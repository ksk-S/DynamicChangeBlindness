import math
import numpy as np
import random
import psychopy.visual
import psychopy.event
import psychopy.core
import gabor_ball

from enum import IntEnum

save_video = False

class ChangeType(IntEnum):
    Rotation = 0
    Shift = 1
    AllRotation = 2

# experiment parameter
# 0 or 1
is_control = 0
control_time = 1.0

change_type = ChangeType.Rotation

# space params
# TODO conssitent case
ScreenSize =[800, 800]
speed = 4
change_angle = 15

# init position
central_pos = [ScreenSize[0]/2,  ScreenSize[1]/2]

win = psychopy.visual.Window(
    size=ScreenSize,
    units="pix",
    color = [.5, .5, .5],
    fullscr=False,
)

stim = gabor_ball.init(central_pos, ScreenSize, win)
gratings = stim["gratings"]
fixation_dot = stim["fixation_dot"]
x_pos = stim["x_pos"]
y_pos = stim["y_pos"]

# store the target positions
target_pos = []
for i in range(gabor_ball.n_patches):
    target_pos.append(gratings[i].pos)

print(target_pos[0][0]) # I don't know why but the array is one level deeper

    
# Define initial positions for each gabor balls
for i in range(gabor_ball.n_patches):
    expand_pos = random.randrange(0, ScreenSize[0]*2 +ScreenSize[1]*2)
    if expand_pos < ScreenSize[0]:
        gratings[i].pos = [expand_pos, 0]
    elif expand_pos < ScreenSize[0]+ScreenSize[1]:
        gratings[i].pos = [ScreenSize[0], expand_pos - ScreenSize[0]]
    elif expand_pos < ScreenSize[0]*2+ScreenSize[1]:
        gratings[i].pos = [expand_pos - ScreenSize[0]-ScreenSize[1], ScreenSize[1]]    
    else:
        gratings[i].pos = [0, expand_pos - ScreenSize[0]*2 - ScreenSize[1]]
        
    gratings[i].pos -= central_pos
    gratings[i].draw()
    
# Calcurate the speeds from the intial positions and the target positions, 
speeds = []
for i in range(gabor_ball.n_patches):
    dist_x = gratings[i].pos[0] - target_pos[i][0][0]
    dist_y = gratings[i].pos[1] - target_pos[i][0][1]
    dist = math.sqrt(dist_x * dist_x + dist_y*dist_y)

    speed =[]
    speed.append(dist_x / 200) # need to adjust the speed
    speed.append(dist_y / 200) # need to adjust the speed

    speeds.append(speed)

def Change():
    index=random.randrange(0,gabor_ball.n_patches)

    if change_type == ChangeType.Rotation:
        gratings[index].ori = gratings[index].ori + change_angle

    elif change_type == ChangeType.AllRotation:
        for index in range(gabor_ball.n_patches):
            gratings[index].ori = gratings[index].ori + change_angle
   
    elif change_type == ChangeType.Shift:
        x_pos[index] = x_pos[index] + 2


def IsNearBy(posA, posB):
    dist_x = posA[0] - posB[0]
    dist_y = posA[1] - posB[1]
    
    dist = math.sqrt(dist_x * dist_x + dist_y*dist_y)
    if dist < 0.1:
        return True
    else:
        return False
    
# loop

clock = psychopy.core.Clock()
keep_going = True

status = 0
target_time = 0

while keep_going:
    
    if ( status == 0 and IsNearBy(gratings[0].pos, target_pos[0][0]) ):
        status = 1

        target_time = clock.getTime()
        if is_control == 0:
            Change()

    elif ( status == 1 and clock.getTime() - target_time > control_time ):
        status = 2

        if is_control == 1:
            Change()
        
    elif ( status == 2 and (gratings[0].pos[0]< 0 or gratings[0].pos[0] > ScreenSize[1] )):
        keep_going = False        

    
    # move coordinates    
    fixation_dot.pos = [central_pos[0] - ScreenSize[0]/2, central_pos[1] - ScreenSize[1]/2]
    fixation_dot.draw()

    for i in range(gabor_ball.n_patches):
        gratings[i].pos = [gratings[i].pos[0] - speeds[i][0], gratings[i].pos[1] - speeds[i][1]] 
        gratings[i].draw()        

    win.flip()

    if save_video:
        win.getMovieFrame()

    #escape
    keys = psychopy.event.getKeys()    
    if len(keys) > 0:
        keep_going = False

if save_video:
    win.saveMovieFrames('original' + ['rotation','shift','allRotation'][int(change_type)] + '_' + ['','control'][is_control] + '.mp4', fps=40)

win.close()
