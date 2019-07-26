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
    
    Rotation = 0
    Shift = 1

class CausalityType(Enum):
    Causal = 0
    SpatialGap = 1
    TemporalGap = 2
    PassThrough = 3


# EXP setting
change_type = ChangeType.Rotation

#causality_type = CausalityType.Causal
causality_type = CausalityType.SpatialGap
#causality_type = CausalityType.TemporalGap
#causality_type = CausalityType.PassThrough

left_visible = False

# Causality params
space_A_and_B = 0
if causality_type == CausalityType.SpatialGap:
    space_A_and_B = 80

    
frames_for_gap = 60


#space params
ScreenSize =[1200, 600]

speed = 2
change_angle = 20

turn_height = 100

r_stimulus = 40
r_grating  = 15
r_total = r_stimulus + r_grating

if causality_type == CausalityType.PassThrough:
    space_A_and_B = -r_total * 2

# define shape
cos60 = math.sqrt(3) / 2
sin60 = 0.5

x_pos = [cos60*r_stimulus,          0, -cos60*r_stimulus, -cos60*r_stimulus,           0,  cos60*r_stimulus]
y_pos = [sin60*r_stimulus, r_stimulus,  sin60*r_stimulus, -sin60*r_stimulus, -r_stimulus, -sin60*r_stimulus]


# init position
central_pos_A = [r_total, ScreenSize[1] /2]
central_pos_B = [ScreenSize[0]/2, ScreenSize[1] /2]

win = psychopy.visual.Window(
    size=ScreenSize,
    units="pix",
    fullscr=False
)

# init gratings
gratings_A = []
for i in range(6):
    gratings_A.append(psychopy.visual.GratingStim(
        win=win,
        size=[r_grating*2, r_grating*2],
        pos = [x_pos[i]+ central_pos_A[0] - ScreenSize[0]/2, y_pos[i] + central_pos_A[1] - ScreenSize[1]/2],
        mask="circle",
        units="pix",
        ori=random.randrange(0,360),
        sf=1.0 / (r_grating *2)
    )
    )

fixation_dot_A = psychopy.visual.Circle(
    win=win,
    pos = [central_pos_A[0] - ScreenSize[0]/2, central_pos_A[1] - ScreenSize[1]/2],
    units="pix",
    radius=3,
    fillColor=[-1] * 3,
    lineColor=[-1] * 3,
    edges=128
)


# gratings_B = []
# for i in range(6):
#     gratings_B.append(psychopy.visual.GratingStim(
#         win=win,
#         size=[r_grating*2, r_grating*2],
#         pos = [x_pos[i]+ central_pos_B[0] - ScreenSize[0]/2, y_pos[i] + central_pos_B[1] - ScreenSize[1]/2],
#         mask="circle",
#         units="pix",
#         ori=random.randrange(0,360),
#         sf=1.0 / (r_grating *2)
#     )
#     )

# fixation_dot_B = psychopy.visual.Circle(
#     win=win,
#     pos = [central_pos_B[0] - ScreenSize[0]/2, central_pos_B[1] - ScreenSize[1]/2],
#     units="pix",
#     radius=3,
#     fillColor=[-1] * 3,
#     lineColor=[-1] * 3,
#     edges=128
# )


stim = gabor_ball.init(central_pos_B, ScreenSize, win)
gratings_B = stim["gratings"]
fixation_dot_B = stim["fixation_dot"]
x_pos = stim["x_pos"]
y_pos = stim["y_pos"]
r_total = gabor_ball.total_diameter/2


def Change():
    index=random.randrange(0,gabor_ball.n_patches)

    gratings_B[index].ori = gratings_B[index].ori + change_angle



clock = psychopy.core.Clock()
keep_going = True

status = 0

frame_count = 0


while keep_going:
#    grating.phase = np.mod(clock.getTime() / 0.5, 1)

    # status changes
    if ( status == 0 and central_pos_A[0] > ScreenSize[0]/ 2 - 2* r_total - space_A_and_B ):

        if causality_type == CausalityType.TemporalGap:
            status = 1
            frame_count = 0
            
        elif causality_type == CausalityType.PassThrough:
            status = 3
            frames_for_gap = 200000
            Change()
            
        else:
            status = 2
            Change()
            
    elif ( status == 1 and frame_count > frames_for_gap ):

        status = 2
        Change()
        
    elif ( central_pos_B[0] > - r_total/2 + ScreenSize[0]  or central_pos_A[0] > - r_total  + ScreenSize[0] ):
        
        keep_going = False


    # move coordinates
    if (status == 0 or status == 3):
        central_pos_A = [central_pos_A[0] + speed, central_pos_A[1] ]
    elif (status == 2):
        central_pos_B = [central_pos_B[0] + speed, central_pos_B[1] ]
    
    #update positions
    if(left_visible):
        fixation_dot_A.pos = [central_pos_A[0] - ScreenSize[0]/2, central_pos_A[1] - ScreenSize[1]/2]
        fixation_dot_A.draw()

    fixation_dot_B.pos = [central_pos_B[0] - ScreenSize[0]/2, central_pos_B[1] - ScreenSize[1]/2]
    fixation_dot_B.draw()



    for i in range(gabor_ball.n_patches):
        if(left_visible):
            gratings_A[i].pos = [x_pos[i]+ central_pos_A[0] - ScreenSize[0]/2, y_pos[i] + central_pos_A[1] - ScreenSize[1]/2]
            gratings_A[i].draw()
        
        gratings_B[i].pos = [x_pos[i]+ central_pos_B[0] - ScreenSize[0]/2, y_pos[i] + central_pos_B[1] - ScreenSize[1]/2]
        gratings_B[i].draw()
    
    win.flip()
    if save_video:
            win.getMovieFrame()

    frame_count += 1
    #escape
    keys = psychopy.event.getKeys()
    
    if len(keys) > 0:
        keep_going = False

if save_video:
    win.saveMovieFrames('causality.mp4', fps=40)

win.close()