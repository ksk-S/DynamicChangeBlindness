import numpy as np
import random
import psychopy.visual
import psychopy.event
import psychopy.core

from enum import Enum

class ChangeType(Enum):
    
    GratingMotion = 0
    GratingColor = 1
    BackgroundColor = 2

#params

# 0 or 1
is_control = 0

change_type = ChangeType.GratingMotion

speed = 1.5
change_angle = 30

ScreenSize =[500, 500]



x_pos = [43,  0, -43, -43,   0,  43]
y_pos = [25, 50, 25, -25, -50, -25]

r_stimulus = 50
r_grating  = 15
r_total = r_stimulus + r_grating


# init position
central_pos = [ScreenSize[0]/2,ScreenSize[1]/2]


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

background = psychopy.visual.GratingStim(
        win=win,
        size=ScreenSize,
        pos = [0,0],
        units="pix",
        ori=-90,
        sf=5.0 / ScreenSize[0],
)


def changeBackground():

    if change_type == ChangeType.GratingMotion:
        background.ori = 0

    elif change_type == ChangeType.GratingColor:
        background.color = [0,0,1]

    elif change_type == ChangeType.BackgroundColor:
        win.color = [0,0,1]
    

clock = psychopy.core.Clock()

keep_going = True


status = 0
start_time = clock.getTime()  

if change_type == ChangeType.BackgroundColor:
    win.color = [1,0,0]

elif change_type == ChangeType.GratingColor:
    background.color = [1,0,0]

while keep_going:
    elapsed_time = clock.getTime() - start_time 
    #print elapsed_time
    # status changes
    if ( status == 0 and elapsed_time > 3 ):
            
        status = 1
        
        changeBackground()

        if is_control == 0:
            index=random.randrange(0,6)
            gratings[index].ori = gratings[index].ori + change_angle
            
            
    elif ( status == 1 and  elapsed_time > 6  ):
        
        status = 2
        
        if is_control == 1:
            index=random.randrange(0,6)
            gratings[index].ori = gratings[index].ori + change_angle
        
    elif ( status == 2 and elapsed_time > 9  ):
        
        keep_going = False


    # move stimulus
    
    #if (status == 0):
    #    central_pos = [central_pos[0], central_pos[1] - speed]
    #elif (status >= 1):
    #    central_pos = [central_pos[0] + speed, central_pos[1]]
    
    
    
    #update pos
    
    background.phase = np.mod(clock.getTime() * speed, 1)
    if change_type == ChangeType.GratingColor or change_type == ChangeType.GratingMotion:
        background.draw()
    
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