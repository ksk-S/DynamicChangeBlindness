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

# experiment parameter
# 0 or 1
is_control = 0

change_type = ChangeType.Shift

# space params
# TODO conssitent case
ScreenSize =[800, 800]
rotation_speed = 4  ### Rotation speed
shift_size = 2.5
change_angle = 15

# init position
#central_pos = [gabor_ball.total_diameter/2, ScreenSize[1] - (gabor_ball.total_diameter/2)]

central_pos = [ScreenSize[0]/2, ScreenSize[1]/2]

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
    index=random.randrange(0,int(gabor_ball.n_patches/2))
    gratings[index].ori = gratings[index].ori - 10
    # x_pos[index] = x_pos[index] + shift_size
    # y_pos[index] = y_pos[index] + shift_size


clock = psychopy.core.Clock()
keep_going = True

status = 0
tmp_t = 0

while keep_going:
#    grating.phase = np.mod(clock.getTime() / 0.5, 1)

    # status changes
    if ( status == 0 and tmp_t > 60):
            
        status = 1
        # for index in range(gabor_ball.n_patches):
        #     gratings[index].ori = gratings[index].ori - rotation_speed
        if is_control == 0:
            Change()
            
    elif ( status == 1 and tmp_t > 100):
        
        status = 2
        
        if is_control == 1:
            Change()
        
    elif ( status == 2 and tmp_t > 150):
        
        keep_going = False
    
    tmp_t = tmp_t + 1

    # move coordinates
    if status == 0:
        # gratings[0].ori = gratings[0].ori + rotation_speed*2.0
        for index in range(int(gabor_ball.n_patches/2.)):
            gratings[index+3].ori = gratings[index+3].ori + rotation_speed*2.0
    elif status >= 1:
        # gratings[0].ori = gratings[0].ori - rotation_speed*2.0
        for index in range(int(gabor_ball.n_patches/2.)):
            gratings[index+3].ori = gratings[index+3].ori - rotation_speed*2.0
    
    #update positions
    fixation_dot.pos = [central_pos[0] - ScreenSize[0]/2, central_pos[1] - ScreenSize[1]/2]
    fixation_dot.draw()
    
    for i in range(gabor_ball.n_patches):
        gratings[i].pos = [x_pos[i]+ central_pos[0] - ScreenSize[0]/2, y_pos[i] + central_pos[1] - ScreenSize[1]/2]
        gratings[i].draw()
    
    win.flip()
    if save_video:
        win.getMovieFrame()

    #escape
    keys = psychopy.event.getKeys()

    if len(keys) > 0:
        keep_going = False

if save_video:
    win.saveMovieFrames('thetas_only_' + ['rotation','shift','allRotation'][int(change_type)] + '_' + ['','control'][is_control] + '.mp4', fps=40)

win.close()