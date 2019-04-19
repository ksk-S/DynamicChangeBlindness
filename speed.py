import math
import numpy as np
import random
import psychopy.visual
import psychopy.event
import psychopy.core

from enum import Enum

class ChangeType(Enum):
    
    Rotation = 0
    Shift = 1

# experiment parameter
# 0 or 1
is_control = 0

change_type = ChangeType.Rotation

#space params
ScreenSize =[1000, 500]

speed = 2
change_angle = 20

r_stimulus = 60
r_grating  = 15
r_total = r_stimulus + r_grating


# define shape
cos60 = math.sqrt(3) / 2
sin60 = 0.5

x_pos = [cos60*r_stimulus,          0, -cos60*r_stimulus, -cos60*r_stimulus,           0,  cos60*r_stimulus]
y_pos = [sin60*r_stimulus, r_stimulus,  sin60*r_stimulus, -sin60*r_stimulus, -r_stimulus, -sin60*r_stimulus]


# init position
central_pos = [r_total, (ScreenSize[1] - r_total)/2]

# timing of gabor change
max_steps = (int) (ScreenSize[0]/speed)
change_step = random.randrange((int)(max_steps/3),(int)(2*max_steps/3))

win = psychopy.visual.Window(
    size=ScreenSize,
    units="pix",
    fullscr=False
)

gratings = []
for i in range(6):
    #print (x_pos[i] , " " , y_pos[i])
    gratings.append(psychopy.visual.GratingStim(
        win=win,
        size=[r_grating*2, r_grating*2],
        pos = [x_pos[i]+ central_pos[0] - ScreenSize[0]/2, y_pos[i] + central_pos[1] - ScreenSize[1]/2],
        mask="circle",
        units="pix",
        ori=random.randrange(0,360),
        sf=1.0 / (r_grating *2)
    )
    )
    
fixation_dot = psychopy.visual.Circle(
    win=win,
    pos = [central_pos[0] - ScreenSize[0]/2, central_pos[1] - ScreenSize[1]/2],
    units="pix",
    radius=3,
    fillColor=[-1] * 3,
    lineColor=[-1] * 3,
    edges=128
)


def change():
    index=random.randrange(0,6)

    if change_type == ChangeType.Rotation:
        
        gratings[index].ori = gratings[index].ori + change_angle

    elif change_type == ChangeType.Shift:
    
        x_pos[index] = x_pos[index] + 2


clock = psychopy.core.Clock()
keep_going = True

step = 0

while keep_going:

    if (step==change_step):
        change()

        if is_control == 0:
            speed = speed*2
    
    if (central_pos[0] > ScreenSize[0] - r_total):
        keep_going = False

    # move coordinates
    central_pos = [central_pos[0] + speed, central_pos[1]]
    
    #update positions
    fixation_dot.pos = [central_pos[0] - ScreenSize[0]/2, central_pos[1] - ScreenSize[1]/2]
    fixation_dot.draw()
    
    for i in range(6):
        gratings[i].pos = [x_pos[i]+ central_pos[0] - ScreenSize[0]/2, y_pos[i] + central_pos[1] - ScreenSize[1]/2]
        gratings[i].draw()
    
    win.flip()

    #escape
    keys = psychopy.event.getKeys()

    if len(keys) > 0:
        keep_going = False

    step = step + 1

win.close()