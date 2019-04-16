import numpy as np
import math
import random
import psychopy.visual
import psychopy.event
import psychopy.core

#params

# speed in rad/step
rotation_speed = math.pi/100

# circle radius
motion_radius = 100

# gabor change (degrees)
change_angle = 20

ScreenSize =[500, 500]

x_pos = [43,  0, -43, -43,   0,  43]
y_pos = [25, 50, 25, -25, -50, -25]

r_stimulus = 50
r_grating  = 15
r_total = r_stimulus + r_grating

# duration of experiment
max_steps = 500

# timing of gabor change
change_step = random.randrange(max_steps/4,3*max_steps/4)

# init position
central_pos = [(ScreenSize[0]/2) - motion_radius, (ScreenSize[1]/2)]
sitmulus_angle = 0

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
    print(index)    
    gratings[index].ori = gratings[index].ori + change_angle


clock = psychopy.core.Clock()
keep_going = True
step = 0

#??
# 1 
# class ChangeType(Enum):
    
#     GratingMotion = 0
#     GratingColor = 1
#     BackgroundColor = 2


while keep_going:

    if (step >= max_steps):
        keep_going = False

    # update positions
    sitmulus_angle = sitmulus_angle + rotation_speed;
    central_pos = [ScreenSize[0]/2 + motion_radius * math.cos(sitmulus_angle), ScreenSize[1]/2 + motion_radius * math.sin(sitmulus_angle)]
    
    fixation_dot.pos = [central_pos[0] - ScreenSize[0]/2, central_pos[1] - ScreenSize[1]/2]
    fixation_dot.draw()
    
    for i in range(6):
        gratings[i].pos = [x_pos[i] + central_pos[0] - ScreenSize[0]/2, y_pos[i] + central_pos[1] - ScreenSize[1]/2]
        gratings[i].draw()

    # change gabor
    if (step == change_step ):
        change()
    
    win.flip()

    #escape
    keys = psychopy.event.getKeys()

    if len(keys) > 0:
        keep_going = False

    step = step + 1

win.close()