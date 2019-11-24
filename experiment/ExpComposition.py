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

change_type = ChangeType.Rotation

# experiment parameter
# 0 or 1
is_control = 0

control_time = 1.0
end_time = 2.5

# space params
n_patches = 6
speed = 4
change_angle = 30

# init position
clock = psychopy.core.Clock()
status = 0
keep_going = False

def Change():
    index=random.randrange(0, stim.n_patches)

    if change_type == ChangeType.Rotation:
        stim.gratings[index].ori = stim.gratings[index].ori + change_angle

    elif change_type == ChangeType.AllRotation:
        for index in range(stim.n_patches):
            stim.gratings[index].ori = stim.gratings[index].ori + change_angle
   
    elif change_type == ChangeType.Shift:
        stim.x_pos[index] = stim.x_pos[index] + 2


def Init(w, s):

    global win, ScreenSize
    
    win = w
    ScreenSize = s

    ResetTrial()
    
def InitialisePos():
    
    global target_pos, speeds

    target_pos = []
    for i in range(stim.n_patches):
        target_pos.append(stim.gratings[i].pos)

    
    # Define initial positions for each gabor balls
    for i in range(stim.n_patches):
        expand_pos = random.randrange(0, ScreenSize[0]*2 +ScreenSize[1]*2)
        if expand_pos < ScreenSize[0]:
            stim.gratings[i].pos = [expand_pos, 0]
        elif expand_pos < ScreenSize[0]+ScreenSize[1]:
            stim.gratings[i].pos = [ScreenSize[0], expand_pos - ScreenSize[0]]
        elif expand_pos < ScreenSize[0]*2+ScreenSize[1]:
            stim.gratings[i].pos = [expand_pos - ScreenSize[0]-ScreenSize[1], ScreenSize[1]]    
        else:
            stim.gratings[i].pos = [0, expand_pos - ScreenSize[0]*2 - ScreenSize[1]]
        
        stim.gratings[i].pos -= central_pos
        stim.gratings[i].draw()

    speeds = []
    for i in range(stim.n_patches):
        dist_x = stim.gratings[i].pos[0] - target_pos[i][0][0]
        dist_y = stim.gratings[i].pos[1] - target_pos[i][0][1]
        dist = math.sqrt(dist_x * dist_x + dist_y*dist_y)
    
        speed =[]
        speed.append(dist_x / 200) # need to adjust the speed
        speed.append(dist_y / 200) # need to adjust the speed
    
        speeds.append(speed)
    
    
def ResetTrial():
    
    global central_pos, stim, clock, status, keep_going
    

    central_pos = [ScreenSize[0]/2,ScreenSize[1]/2]
    stim = gabor_ball.init(central_pos, ScreenSize, win, n_patches)
    
    InitialisePos()
    
    clock = psychopy.core.Clock()
    status = 0
    keep_going = True
    

def StartTrial(condition):

    global is_control

    is_control = condition

    ResetTrial()

    while keep_going:
        Update()
        
    if save_video:
        win.saveMovieFrames('original' + ['rotation','shift','allRotation'][int(change_type)] + '_' + ['','control'][is_control] + '.mp4', fps=40)
 

def StimTurnAround():
    for i in range(stim.n_patches):
        x = speeds[i][0]
        y = speeds[i][1]
        
        if(random.randrange(0, 1) == 0):
            speeds[i][0] = y
            speeds[i][1] = -x
        else:
            speeds[i][0] = -y
            speeds[i][1] = x
    

def Update():
    
    global status, keep_going, central_pos, stim, target_pos, target_time
    
    if ( status == 0 and IsNearBy(stim.gratings[0].pos, target_pos[0][0]) ):
        status = 1

        StimTurnAround()

        target_time = clock.getTime()
        if is_control == 0:
            Change()

    elif ( status == 1 and clock.getTime() - target_time > control_time ):
        status = 2

        if is_control == 1:
            Change()
        
    elif ( status == 2 and clock.getTime() - target_time > end_time):
        keep_going = False        
    
    #update positions
    stim.fixation_dot.pos = [central_pos[0] - ScreenSize[0]/2, central_pos[1] - ScreenSize[1]/2]
    stim.fixation_dot.draw()
    
    for i in range(stim.n_patches):
        stim.gratings[i].pos = [stim.gratings[i].pos[0] - speeds[i][0], stim.gratings[i].pos[1] - speeds[i][1]] 
        stim.gratings[i].draw()
    
    win.flip()
    
    if save_video:
        win.getMovieFrame()

    #escape
    keys = psychopy.event.getKeys()

    if len(keys) > 0:
        keep_going = False


def IsNearBy(posA, posB):
    dist_x = posA[0] - posB[0]
    dist_y = posA[1] - posB[1]
    
    dist = math.sqrt(dist_x * dist_x + dist_y*dist_y)
    if dist < 0.1:
        return True
    else:
        return False