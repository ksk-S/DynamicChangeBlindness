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



class MovementType(IntEnum):
    
    RandomToRandom = 0
    RandomToCircle = 1
    CircleToRandom = 2

movement_type = MovementType.RandomToRandom


# experiment parameter
# 0 or 1
is_control = 0

control_time = 0.5
end_time = 2.0

# space params
n_patches = 6

grating_speed = 4

init_patches_distance = 400

change_angle = 15

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
    
    global stim, target_pos, fixation_speed, speeds

    target_pos = []
    for i in range(stim.n_patches):
        target_pos.append(stim.gratings[i].pos)
  
    for i in range(stim.n_patches):
        angle = random.randrange(0.0, 360.0)
        stim.gratings[i].pos = [init_patches_distance * math.cos(angle), init_patches_distance * math.sin(angle)]

    fixation_speed =  (ScreenSize[1] - gabor_ball.total_diameter) * ( grating_speed / ScreenSize[0] )

    speeds = []
    for i in range(stim.n_patches):
        dist_x = stim.gratings[i].pos[0] - target_pos[i][0][0]
        dist_y = stim.gratings[i].pos[1] - target_pos[i][0][1]
        dist = math.sqrt(dist_x * dist_x + dist_y*dist_y)
    
        speed =[]
        speed.append(dist_x * (grating_speed / ScreenSize[0]) ) # need to adjust the speed
        speed.append(dist_y * (grating_speed / ScreenSize[1]) ) # need to adjust the speed
    
        speeds.append(speed)
    
    
    
def ResetTrial():
    
    global central_pos, stim, clock, status, keep_going
    
    #target pos
    central_pos = [gabor_ball.total_diameter, gabor_ball.total_diameter]

    
    stim = gabor_ball.init(central_pos, ScreenSize, win, n_patches)
    
    
    InitialisePos()
    
    #initial pos
    central_pos = [gabor_ball.total_diameter, ScreenSize[1]]
    
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
    
    
    if (status == 0):
        central_pos = [central_pos[0], central_pos[1] - fixation_speed]
    else:
        central_pos = [central_pos[0]+ fixation_speed, central_pos[1] ]
        
    stim.fixation_dot.pos = [central_pos[0] - ScreenSize[0]/2, central_pos[1] - ScreenSize[1]/2]
    stim.fixation_dot.draw()
    
    if movement_type == MovementType.RandomToRandom:
    
        for i in range(stim.n_patches):
            stim.gratings[i].pos = [stim.gratings[i].pos[0] - speeds[i][0], stim.gratings[i].pos[1] - speeds[i][1]] 
            stim.gratings[i].draw()
    
    elif movement_type == MovementType.RandomToCircle:
        if status == 0:
            for i in range(stim.n_patches):
                stim.gratings[i].pos = [stim.gratings[i].pos[0] - speeds[i][0], stim.gratings[i].pos[1] - speeds[i][1]] 
                stim.gratings[i].draw()
        else:
            for i in range(stim.n_patches):
                stim.gratings[i].pos = [stim.x_pos[i]+ central_pos[0] - ScreenSize[0]/2, stim.y_pos[i] + central_pos[1] - ScreenSize[1]/2]
                stim.gratings[i].draw()
    
    elif movement_type == MovementType.CircleToRandom:
        if status == 0:
            for i in range(stim.n_patches):
                stim.gratings[i].pos = [stim.x_pos[i]+ central_pos[0] - ScreenSize[0]/2, stim.y_pos[i] + central_pos[1] - ScreenSize[1]/2]
                stim.gratings[i].draw()
        else:
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