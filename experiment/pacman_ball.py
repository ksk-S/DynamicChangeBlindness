import math
import psychopy.visual
import random

# matlab parameters
# gaborbgcolor = [.5 .5 .5 0];
# bgcolor = gaborbgcolor.*256-1;
# gaborSize = 150;
# freq = gaborSize/5000;
# sc = gaborSize/10; %Correlated with patch size

# fixation
# crossSize = 10;
# crossThickness = 2;

# change
# movementIncrement = 10;
# rotationSize = 15;


# the main stimulus
# use this to control the size of all stimuli
size_factor = 3

contrast = 1

gabor_diameter  = 35 * size_factor
stimulus_diameter = gabor_diameter * 2
total_diameter = gabor_diameter + stimulus_diameter
gabor_freq = gabor_diameter/5000 # r_grating/800,
# gauss spatial constant
sc = gabor_diameter/10


def init(central_pos, ScreenSize, win, n_patches = 6, size_factor=1.5):

    global gabor_diameter, stimulus_diameter, total_diameter, gabor_freq, sc

    # variables for stimulus
    # radii
    gabor_diameter  = 35 * size_factor
    stimulus_diameter = gabor_diameter * 3.0
    total_diameter = gabor_diameter + stimulus_diameter
    gabor_freq = gabor_diameter/5000 # r_grating/800,
    # gauss spatial constant
    sc = gabor_diameter/10

    # define shape
    sd = stimulus_diameter
    x_pos = []
    y_pos = []
    gratings = []

    ini_phase = math.pi*2*random.random()
    for i in range(n_patches):
        angle = math.pi*2*i/n_patches #+ 0.25*math.pi
        ttt = math.pi - angle
        x_pos.append(math.cos(angle+ini_phase)*sd/2)
        y_pos.append(math.sin(angle+ini_phase)*sd/2)

        # grating_pos = [x_pos[i] - ScreenSize[0]/2, y_pos[i] - ScreenSize[1]/2]
        grating_pos = [x_pos[i]+ central_pos[0] - ScreenSize[0]/2, y_pos[i] + central_pos[1] - ScreenSize[1]/2]
        gratings.append(psychopy.visual.ImageStim(
            win = win,
            image="imgs/dynamic_cb.png",
            size = [gabor_diameter, gabor_diameter],
            pos = grating_pos,
            # blurr edge only
            maskParams = {"sd": 5}, #{"fringeWidth":0.2},
            units = "pix",
            ori = ttt*180.0/math.pi,
            # ori = random.randrange(0,360),
            # sf = gabor_freq,
            color = [1, 1, 1],
            contrast = contrast,
        )
        )
    
    fixation_dot = psychopy.visual.Circle(
        win=win,
        pos = [central_pos[0], central_pos[1] - ScreenSize[1]/2],
        units="pix",
        radius=3,
        fillColor=[-1] * 3,
        lineColor=[-1] * 3,
        edges=128
    )
    ball = type('', (), {})
    ball.gratings = gratings
    ball.fixation_dot = fixation_dot
    ball.x_pos = x_pos
    ball.y_pos = y_pos
    ball.n_patches = n_patches

    return ball