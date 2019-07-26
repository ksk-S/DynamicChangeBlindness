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
ScreenSize =[1200, 800]

speed = 2
change_angle = 20

turn_height = 100

r_stimulus = 60
r_grating  = 15
r_total = r_stimulus + r_grating


# define shape
cos60 = math.sqrt(3) / 2
sin60 = 0.5

x_pos = [cos60*r_stimulus,          0, -cos60*r_stimulus, -cos60*r_stimulus,           0,  cos60*r_stimulus]
y_pos = [sin60*r_stimulus, r_stimulus,  sin60*r_stimulus, -sin60*r_stimulus, -r_stimulus, -sin60*r_stimulus]


# init position
central_pos = [r_total, ScreenSize[1] /2]


win = psychopy.visual.Window(
    size=ScreenSize,
    units="pix",
    fullscr=False
)

# init gratings
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


def Change():
    print("change")
    index=random.randrange(0,6)

    if change_type == ChangeType.Rotation:
        
        gratings[index].ori = gratings[index].ori + change_angle

    elif change_type == ChangeType.Shift:
    
        x_pos[index] = x_pos[index] + 2


clock = psychopy.core.Clock()
keep_going = True

status = 0

turn_count = 0

while keep_going:
#    grating.phase = np.mod(clock.getTime() / 0.5, 1)

    # status changes
    if ( status == 0 and central_pos[1] < ScreenSize[1]/ 2 - turn_height  ):
        
        if(turn_count ==3):
            Change()
            status = 2
            continue
        
        status = 1
        turn_count += 1
        
        Change()
            
    elif ( status == 1 and central_pos[1] > ScreenSize[1]/ 2 + turn_height ):

        if(turn_count ==3):
            Change()
            status = 2
            continue
        
        status = 0
        turn_count += 1

        Change()
        
        
    elif ( central_pos[0] > - r_total  + ScreenSize[0] ):
        
        keep_going = False


    # move coordinates
    if (status == 0):
        central_pos = [central_pos[0] + speed, central_pos[1] - speed]
    elif (status >= 1):
        central_pos = [central_pos[0] + speed, central_pos[1] + speed]
    
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

win.close()