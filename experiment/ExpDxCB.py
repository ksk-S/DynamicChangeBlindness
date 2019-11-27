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

change_type = ChangeType.Shift

# experiment parameter
# 0 or 1
is_control = 0

# space params
n_patches = 6

rotation_speed = 10  ### Rotation speed
shift_rate = 0.1

# init position
clock = psychopy.core.Clock()
status = 0
tmp_t = 0
keep_going = False

def Change():

    index=random.randrange(0, stim.n_patches)

    if change_type == ChangeType.Rotation:
        stim.gratings[index].ori = stim.gratings[index].ori + change_angle

    elif change_type == ChangeType.AllRotation:
        for index in range(stim.n_patches):
            stim.gratings[index].ori = stim.gratings[index].ori + change_angle
   
    elif change_type == ChangeType.Shift:
        index=random.randrange(0, stim.n_patches)
        if random.random() > 0.5:
            stim.x_pos[index] = stim.x_pos[index]*(1.+shift_rate)
            stim.y_pos[index] = stim.y_pos[index]*(1.+shift_rate)
        else:
            stim.x_pos[index] = stim.x_pos[index]*(1.-shift_rate)
            stim.y_pos[index] = stim.y_pos[index]*(1.-shift_rate)


def Init(w, s):

    global win, ScreenSize
    
    win = w
    ScreenSize = s

    ResetTrial()
    

def ResetTrial():
    
    global central_pos, stim, clock, status, keep_going, ScreenSize, tmp_t
    
    central_pos = [ScreenSize[0]/2, ScreenSize[1]/2]
    stim = gabor_ball.init(central_pos, ScreenSize, win, n_patches)
    clock = psychopy.core.Clock()
    status = 0
    tmp_t = 0
    keep_going = True
    

def StartTrial(condition):

    global is_control

    is_control = condition

    ResetTrial()

    while keep_going:
        Update()
        
    if save_video:
        win.saveMovieFrames('dx_cb' + ['rotation','shift','allRotation'][int(change_type)] + '_' + ['','control'][is_control] + '.mp4', fps=40)
 

def Update():
    
    global status, keep_going, central_pos, stim, tmp_t, ScreenSize, clock

    # status changes
    if ( status == 0 and tmp_t > 60):
            
        status = 1
        for index in range(stim.n_patches):
            stim.gratings[index].ori = stim.gratings[index].ori - rotation_speed
        if is_control == 0:
            Change()
            
    elif ( status == 1 and tmp_t > 120):
        
        status = 2
        
        if is_control == 1:
            Change()
        
    elif ( status == 2 and tmp_t > 180):
        
        keep_going = False
    
    tmp_t = tmp_t + 1

    # move coordinates
    if status == 0:
        for index in range(stim.n_patches):
            stim.gratings[index].ori = stim.gratings[index].ori + rotation_speed
        # gratings[1].ori = gratings[1].ori + rotation_speed
    elif status >= 1:
        for index in range(stim.n_patches):
            stim.gratings[index].ori = stim.gratings[index].ori - rotation_speed
        # gratings[1].ori = gratings[1].ori - rotation_speed
    #    for index in range(gabor_ball.n_patches):
        #    gratings[index].ori = gratings[index].ori + rotation_speed
    
    #update positions
    stim.fixation_dot.pos = [central_pos[0] - ScreenSize[0]/2, central_pos[1] - ScreenSize[1]/2]
    stim.fixation_dot.draw()
    
    for i in range(stim.n_patches):
        stim.gratings[i].pos = [stim.x_pos[i] + central_pos[0] - ScreenSize[0]/2, stim.y_pos[i] + central_pos[1] - ScreenSize[1]/2]
        stim.gratings[i].draw()
    
    win.flip()
    if save_video:
        win.getMovieFrame()

    #escape
    keys = psychopy.event.getKeys()

    if len(keys) > 0:
        keep_going = False