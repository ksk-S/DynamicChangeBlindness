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

# space params
n_patches = 6
speed = 4
change_angle = 15

# init position
clock = psychopy.core.Clock()
status = 0
keep_going = False

# for creating 2 subgroups of gabors
beta_idx = 0

def Change():
    global beta_idx
    index = beta_idx
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
    

def ResetTrial():
    
    global central_pos, stim, clock, status, keep_going, beta_idx
    
    central_pos = [gabor_ball.total_diameter/2, ScreenSize[1] - (gabor_ball.total_diameter/2)]
    stim = gabor_ball.init(central_pos, ScreenSize, win, n_patches)
    clock = psychopy.core.Clock()
    status = 0
    keep_going = True

    for index in range(stim.n_patches):
        stim.gratings[index].ori = 0
    beta_idx = random.randrange(0, stim.n_patches-1)
    stim.gratings[beta_idx].ori = random.randrange(30, 150)
    stim.gratings[beta_idx+1].ori = random.randrange(30, 150)
    

def StartTrial(condition):

    global is_control

    is_control = condition

    ResetTrial()

    while keep_going:
        Update()
        
    if save_video:
        win.saveMovieFrames('original' + ['rotation','shift','allRotation'][int(change_type)] + '_' + ['','control'][is_control] + '.mp4', fps=40)
 

def Update():
    
    global status, keep_going, central_pos, stim
    
    # status changes
    if ( status == 0 and central_pos[1] < (gabor_ball.total_diameter/2) ):
            
        status = 1
        
        if is_control == 0:
            Change()
            
    elif ( status == 1 and central_pos[0] >  ScreenSize[0]/ 2 ):
        
        status = 2
        
        if is_control == 1:
            Change()
        
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
        stim.gratings[i].pos = [stim.x_pos[i]+ central_pos[0] - ScreenSize[0]/2, stim.y_pos[i] + central_pos[1] - ScreenSize[1]/2]
        stim.gratings[i].draw()
    
    win.flip()
    
    if save_video:
        win.getMovieFrame()

    #escape
    keys = psychopy.event.getKeys()

    if len(keys) > 0:
        keep_going = False