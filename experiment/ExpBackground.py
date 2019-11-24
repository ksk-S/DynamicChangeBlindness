import math
import numpy as np
import random
import psychopy.visual
import psychopy.event
import psychopy.core
import gabor_ball

from enum import Enum

save_video = False

class ChangeType(Enum):
    
    GratingMotion = 0
    GratingColor = 1
    BackgroundColor = 2

change_type = ChangeType.GratingMotion

# experiment parameter
# 0 or 1
is_control = 0

n_patches = 6
speed = 10
change_angle = 30

# init position
clock = psychopy.core.Clock()
status = 0
keep_going = False


    
def Init(w, s):

    global win, ScreenSize
    
    win = w
    ScreenSize = s

    ResetTrial()


def CreateDots():
    
    global dot_stim
    
    n_dots = 1000
    

    dot_stim = psychopy.visual.DotStim(
        win=win,
        units="pix",
        nDots=n_dots,
        coherence = 1,
        fieldPos=(0, 0),
        fieldSize=(800, 800),
        dotSize=5.0,
        dotLife=2,
        dir=0.0,
        speed=speed,
        color=(0.0, 0.0, 0.0),
        opacity=1.0
    )
    
def ResetTrial():
    
    global central_pos, stim, clock, status, keep_going, background, start_time
    
    central_pos = [ScreenSize[0]/2,ScreenSize[1]/2]
    stim = gabor_ball.init(central_pos, ScreenSize, win, n_patches)
    clock = psychopy.core.Clock()
    status = 0
    keep_going = True
    
    start_time = clock.getTime()  

    CreateDots()
    
    
    
    if change_type == ChangeType.BackgroundColor:
        window.color = [1,0,0]

    elif change_type == ChangeType.GratingColor:
        background.color = [1,0,0]


def StartTrial(condition):

    global is_control

    is_control = condition

    ResetTrial()

    while keep_going:
        Update()
        
    if save_video:
        win.saveMovieFrames('original' + ['rotation','shift','allRotation'][int(change_type)] + '_' + ['','control'][is_control] + '.mp4', fps=40)
 

def changeBackground():

    if change_type == ChangeType.GratingMotion:
        dot_stim.dir = 90

    elif change_type == ChangeType.GratingColor:
        dot_stim.color = [0,0,1]

    elif change_type == ChangeType.BackgroundColor:
        window.color = [0,0,1]
        
        
def Update():
    
    global status, keep_going, central_pos, stim
    
    elapsed_time = clock.getTime() - start_time 
    if ( status == 0 and elapsed_time > 3 ):
            
        status = 1
        
        changeBackground()

        if is_control == 0:
            index=random.randrange(0,6)
            stim.gratings[index].ori = stim.gratings[index].ori + change_angle
            
            
    elif ( status == 1 and  elapsed_time > 5  ):
        
        status = 2
        
        if is_control == 1:
            index=random.randrange(0,6)
            stim.gratings[index].ori = stim.gratings[index].ori + change_angle
        
    elif ( status == 2 and elapsed_time > 8  ):
        
        keep_going = False
    
    
    dot_stim.draw()

    
    #update positions
    #stim.fixation_dot.pos = [central_pos[0] - ScreenSize[0]/2, central_pos[1] - ScreenSize[1]/2]
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


