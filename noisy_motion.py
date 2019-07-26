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

change_type = ChangeType.Rotation

# space params
# TODO conssitent case
ScreenSize =[800, 800]
speed = 4
change_angle = 15

# initial position
central_pos = [gabor_ball.total_diameter/2, ScreenSize[1] - (gabor_ball.total_diameter/2)]
# initial value for sin wave
wave_angle = 0
wave_amplitude = gabor_ball.stimulus_diameter/2

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


def Change():
    index=0 #random.randrange(0,gabor_ball.n_patches)
    print("rotated ", index)
    print("old ", gratings[index].ori)

    if change_type == ChangeType.Rotation:
        gratings[index].ori = gratings[index].ori + change_angle


    elif change_type == ChangeType.AllRotation:
        for index in range(gabor_ball.n_patches):
            gratings[index].ori = gratings[index].ori + change_angle
   
    elif change_type == ChangeType.Shift:
        x_pos[index] = x_pos[index] + 2

    print("new ", gratings[index].ori)


clock = psychopy.core.Clock()
keep_going = True

status = 0

while keep_going:
    wave_angle = wave_angle + math.pi/200
#    grating.phase = np.mod(clock.getTime() / 0.5, 1)

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


    # move coordinates
    # use original trajectory + sine wave
    if (status == 0):
        central_pos = [central_pos[0], central_pos[1] - speed]
        # add sine wave on gabor 1 and cos on gabor 2
    elif (status >= 1):
        central_pos = [central_pos[0] + speed, central_pos[1]]
    
    #update positions
    fixation_dot.pos = [central_pos[0] - ScreenSize[0]/2, central_pos[1] - ScreenSize[1]/2]
    fixation_dot.draw()

    
    for i in range(gabor_ball.n_patches):
        angle = math.pi*2*i/gabor_ball.n_patches
        if(i == 0):
            x_pos[i] = wave_amplitude*math.cos(angle + wave_angle)
            y_pos[i] = wave_amplitude*math.sin(angle + wave_angle)
        else:
            x_pos[i] = wave_amplitude*math.cos(angle - wave_angle)
            y_pos[i] = wave_amplitude*math.sin(angle - wave_angle)

        gratings[i].pos = [x_pos[i]+ central_pos[0] - ScreenSize[0]/2, y_pos[i] + central_pos[1] - ScreenSize[1]/2]
    
        # if (i == 0):
        #     gratings[i].pos[0] = gratings[i].pos[0] + wave_amplitude * math.sin(wave_angle)
        # else:
        #     gratings[i].pos[0] = gratings[i].pos[0] + wave_amplitude * math.cos(wave_angle)
        # line = psychopy.visual.Line(win, start=fixation_dot.pos, end=gratings[i].pos, lineColor = [0,0,0])
        # line.draw()
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