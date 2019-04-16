import numpy as np
import random
import psychopy.visual
import psychopy.event
import psychopy.core

#params

# speed = 2
rotation_speed = 
# degrees
change_angle = 20

ScreenSize =[500, 500]

x_pos = [43,  0, -43, -43,   0,  43]
y_pos = [25, 50, 25, -25, -50, -25]

r_stimulus = 50
r_grating  = 15
r_total = r_stimulus + r_grating

# duration of experiment
max_steps = 100

# timing of gabor change
change_step = random.randrange(max_steps/4,3*max_steps/4)

# init position
central_pos = [r_total, ScreenSize[1] - r_total]
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

    if (step >= max_steps)
        keep_going = False

    # move coordinates
    sitmulus_angle = sitmulus_angle + rotationSpeed);

            double d = mainVariables[3];
            float ballX = 200 + (float) ( d*Math.cos(newAngle));
            float ballY = 200 + (float) ( d*Math.sin(newAngle));

            float[] newPosition = {ballX,ballY,newAngle,mainVariables[3]};
            return newPosition;

    if (status == 0):
        central_pos = [central_pos[0], central_pos[1] - speed]
    elif (status >= 1):
        central_pos = [central_pos[0] + speed, central_pos[1]]
    
    #update positions
    fixation_dot.pos = [central_pos[0] - ScreenSize[0]/2, central_pos[1] - ScreenSize[1]/2]
    fixation_dot.draw()
    
    for i in range(6):
        gratings[i].pos = [x_pos[i] + central_pos[0] - ScreenSize[0]/2, y_pos[i] + central_pos[1] - ScreenSize[1]/2]
        gratings[i].draw()
    
    win.flip()

    #escape
    keys = psychopy.event.getKeys()

    if len(keys) > 0:
        keep_going = False

    step = step + 1

win.close()