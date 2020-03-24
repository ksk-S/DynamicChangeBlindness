## Setup
Open a terminal.

Install your favorite python env. eg, anaconda, virtualenv...

### Run test code

Create environment
`conda create -n cblindness37 python=3.7
source activate cblindness37
cd path/to/echo_state_network`

Install missing packages
`conda install -n cblindness37 Matplotlib Numpy Scipy`
`pip install easyesn Pillow`


If necessary, create `images` and `csv_files` in `path/to/echo_state_network`

Run
`python switch_tracker.py`


## References

Model of prediction-related change blindness: partial access here
https://books.google.co.uk/books?id=clVkDQAAQBAJ&pg=PA387&lpg=PA387&dq=NEURAL+MODELS+OF+PREDICTION+AND+SUSTAINED+INATTENTIONAL+BLINDNESS&source=bl&ots=9x_pkVTxli&sig=ACfU3U2uQbIYlzfovTGK-jdmACCKcbrdqg&hl=en&sa=X&ved=2ahUKEwi3g9-ojcXgAhXBXhUIHQpJAGMQ6AEwAnoECAcQAQ#v=onepage&q=NEURAL%20MODELS%20OF%20PREDICTION%20AND%20SUSTAINED%20INATTENTIONAL%20BLINDNESS&f=false

## Summary of model

### Input
Video of 1 blue ball bouncing in 2D: stimulus that should be tracked.
Unexpected stimulus is a green ball that randomly crosses the screen from left to right. 
Visual area divided in 4x4 grid for the video input (average pixels in each cell)
Green or blue channels only -> 4*4*2
Ball size seems to be about 1/3 of one cell.
Unexpected stimulus appears on the left, y = random.

### ESN 
Echo State Netowrk (ESN) with 100 reservoir neurons (https://cloudfront.escholarship.org/dist/prd/content/qt7z89g7p6/qt7z89g7p6.pdf) takes 32 inputs from the grid
//no reference for this claim: Video is updated every 10 steps of the ESN (1 frame every 10 steps)
2 types of perceptrons readout from the ESN: detection perceptrons DP and tracking perceptrons TP (all fully connected to the ESN)
Average activity of the TP is fed back to the ESN through one input neuron.

### Task
Average of DP must be >=0.5 (<0.5) to indicate presence (absence) of the unexpected stimulus, average of TP must be >=0.5 (<0.5) to indicate upwards (downwards) motion of the tracked stimulus. Right and left motions are ignored.

100 DP and 100 TP are trained for each ESN, but during evaluation only the best DP and the best TP are kept.
100 ESN are trained.
Training time: 10,000 time steps, an average of 40 changes of direction and 10 appearances of the unexpected object.
Evaluation: 1,000 time steps, unexpected object present 50% of the time. Performance = percentage of time steps where the perceptron gave the correct answer. Keep best perceptron.

Reapeat with feedback from tracking.

## Code

Deprecated
Code based on https://github.com/Devrim-Celik/simple_echo_state_network

New
https://github.com/kalekiu/easyesn

### Input

See reasonning below.
Training:
ball size = 30px
speed = 15 pix/frame
initial position = (40,40)
unexpected input speed = (15,0)
unexpected apparition time = random(0,75) reset every 100 frames.
unexpected position = (0,rand(400))


"Ball size seems to be about 1/3 of one cell"
400x400 image, ball size = (400/4)/3 = 30 px

//I cannot find reference for the repetition of frames and it may not be critical
"1 frame every 10 steps"
"10,000 time steps, an average of 40 changes of direction"
video frames nf = 10,000/10 = 1,000
change direction after x frames = nf/40 = 1,000/40 = 25
speed s pix/frame ~= scene size / x = 400/25 = 16, let's say 15 p/f

 
"10 appearances of the unexpected object"
frequency = 1,000/10 =  1 per 100 frames
takes 25 frames to cross
apparition time = random(0,75) reset every 100 frames.

"unexpected object present 50% of the time"
100/2 = 50 frames


