import math
import numpy as np
import random
import psychopy.visual
import psychopy.event
import psychopy.core
import gabor_ball

from enum import IntEnum

save_video = True

class ChangeType(IntEnum):
    
    Rotation = 0
    Shift = 1
    AllRotation = 2

change_type = ChangeType.Rotation


num_waiting_patch = 1

# experiment parameter
# 0 or 1
is_control = 0

# space params
n_patches = 6
speed = 4
change_angle = 15

# init position
clock = psychopy.core.Clock()
status = 0
keep_going = False

# waiting index
waiting_indice = []
moving_indice = []
    

def Change():
    indice = []
    if is_control == 1:
        indice = waiting_indice;
    else:
        indice = moving_indice;
        
    for i in indice:
    
        if change_type == ChangeType.Rotation:
            stim.gratings[i].ori = stim.gratings[i].ori + change_angle

        elif change_type == ChangeType.AllRotation:
            for index in range(stim.n_patches):
                stim.gratings[i].ori = stim.gratings[i].ori + change_angle
   
        elif change_type == ChangeType.Shift:
            stim.x_pos[i] = stim.x_pos[i] + 2

def Init(w, s):

    global win, ScreenSize
    
    win = w
    ScreenSize = s

    ResetTrial()
    

def InitialisePos():
    
    global stim, waiting_indice, moving_indice

    waiting_indice = []
    moving_indice =  []
    
    for i in range(stim.n_patches):
        moving_indice.append(i)
    
    for i in range(num_waiting_patch):
        
        index = random.randrange(0, stim.n_patches)
        while index in waiting_indice:
            index = random.randrange(0, stim.n_patches)
        waiting_indice.append(index)
        moving_indice.remove(index)

    #print(waiting_indice)
    #print(moving_indice)
    
    # set the corner position
    central_pos = [gabor_ball.total_diameter/2, gabor_ball.total_diameter/2]
    stim = gabor_ball.init(central_pos, ScreenSize, win, n_patches)
    
    # store the patch positions at the corner
    target_pos = []
    for i in range(len(waiting_indice)):
        target_pos.append(stim.gratings[waiting_indice[i]].pos)

    # rest to the initial pos
    central_pos = [gabor_ball.total_diameter/2, ScreenSize[1] - (gabor_ball.total_diameter/2)]
    stim = gabor_ball.init(central_pos, ScreenSize, win, n_patches)

    # restore the patches position for 'waiting' patches
    for i in range(len(waiting_indice)):
        stim.gratings[waiting_indice[i]].pos = target_pos[i]

def ResetTrial():
    
    global central_pos, stim, clock, status, keep_going
    
    
    central_pos = [gabor_ball.total_diameter/2, ScreenSize[1] - (gabor_ball.total_diameter/2)]
    stim = gabor_ball.init(central_pos, ScreenSize, win, n_patches)
    
    InitialisePos();

    
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
        win.saveMovieFrames('separation_' + ['rotation','shift','allRotation'][int(change_type)] + '_' + ['','control'][is_control] + '.mp4', fps=40)
 

def Update():
    
    global status, keep_going, central_pos, stim
    
    # status changes
    if ( status == 0 and central_pos[1] < (gabor_ball.total_diameter/2) ):
            
        status = 1
        
        Change()
        
        #if is_control == 0:
            
    elif ( status == 1 and central_pos[0] >  ScreenSize[0]/ 2 ):
        
        status = 2
        
        #if is_control == 1:
        #    Change()
        
    elif ( status == 2 and central_pos[0] > - (gabor_ball.total_diameter/2)  + ScreenSize[0] ):
        
        keep_going = False


    #control movement direction
    if (status == 0):
        central_pos = [central_pos[0], central_pos[1] - speed]
    elif (status >= 1):
        central_pos = [central_pos[0] + speed, central_pos[1]]
    
    #update positions
    stim.fixation_dot.pos = [central_pos[0] - ScreenSize[0]/2, central_pos[1] - ScreenSize[1]/2]
    stim.fixation_dot.draw()
    
    for i in range(stim.n_patches):
        if not (status == 0 and i in waiting_indice):
            stim.gratings[i].pos = [stim.x_pos[i]+ central_pos[0] - ScreenSize[0]/2, stim.y_pos[i] + central_pos[1] - ScreenSize[1]/2]
        stim.gratings[i].draw()
    
    win.flip()
    
    if save_video:
        win.getMovieFrame()

    #escape
    keys = psychopy.event.getKeys()

    if len(keys) > 0:
        keep_going = False

