import psychopy.event
import random

import DataCtl
import ExpOriginal
import ExpBackground
import ExpComposition
import ExpDxCB

num_trials = 20

ScreenSize =[800, 800]

subId = 0

ExperimentTypes = {
  "original"                : False,
  "background"              : False,
  "compositionRandToRand"   : True,
  "compositionRandToCircle" : True,
  "compositionCicleToRand"  : True,
  "dx_cb"           : True,
}

def StartExp(info):
    
    global subId
    
    subId = info["SubId"]
    
    DataCtl.Init(info)
    DataCtl.SaveHeader()
    
    
    for exp in ExperimentTypes:
        if ExperimentTypes[exp]:
            print(exp)
            StartCondition(exp)

def CreateTrialList():
    
    list1 = [0 for i in range(int(num_trials/2))]
    list2 = [1 for i in range(int(num_trials/2))]
    
    list1.extend(list2)
    
    random.shuffle(list1)

    return list1


def ExpInitHandler(win, ScreenSize, callback):
    return func(*args)
    

def ExpRunHandler(func, *arg):
    return func(*args)




def StartCondition(exp):
    
    global win
    
    trial_list = CreateTrialList()
    print(trial_list)
    
    win = psychopy.visual.Window(
        size=ScreenSize,
        units="pix",
        color = [.5, .5, .5],
        fullscr=False,
    )
    
    if exp =='original':
    
        ExpOriginal.Init(win, ScreenSize)
        
    elif exp == 'background':
        
        ExpBackground.Init(win, ScreenSize)
       
    elif exp == 'compositionRandToRand':
        
        ExpComposition.Init(win, ScreenSize)
        
    elif exp == 'compositionRandToCircle':
        
        ExpComposition.Init(win, ScreenSize)
        
    elif exp == 'compositionCicleToRand':
        
        ExpComposition.Init(win, ScreenSize)
<<<<<<< HEAD
=======
    
    elif exp == 'dx_cb':

        ExpDxCB.Init(win, ScreenSize)

>>>>>>> 1a437d9f967978681753b99325b3239ce36fe100
        
    ShowInstruction()
    
    
    for trialId in range(num_trials):
        condition = trial_list[trialId]
        
        if exp =='original':
    
            ExpOriginal.StartTrial(condition)
        
        elif exp == 'background':
            
            ExpBackground.StartTrial(condition)
        
        elif exp == 'compositionRandToRand':
            ExpComposition.movement_type = ExpComposition.MovementType.RandomToRandom
            ExpComposition.StartTrial(condition)
        
        elif exp == 'compositionRandToCircle':
            
            ExpComposition.movement_type = ExpComposition.MovementType.RandomToCircle
            ExpComposition.StartTrial(condition)
        
        elif exp == 'compositionCicleToRand':
            
            ExpComposition.movement_type = ExpComposition.MovementType.CircleToRandom
            ExpComposition.StartTrial(condition)
        
        elif exp == 'dx_cb':

            ExpDxCB.StartTrial(condition)
        
        resKey = GetResponse()
        
        print( 'target' if condition == 0 else 'control' )
        
        response = 1 if resKey[0] == "right" else 0
        
        DataCtl.SaveData([subId, exp, trialId, condition, response])
    
    win.close()
    

def ShowInstruction():
   
    text = psychopy.visual.TextStim(
        win=win,
        pos=(0.0, 0.0),
        text="In this session, you will see a moving object consisting of multiple grating patches. You will need to fixate the central dot at all times. Your task is to detect any changes in one of the grating patches.\n\nPress any key to start",
        color=[-1, -1, -1]
    )
    
    text.draw()
    win.flip()
    
    keys = psychopy.event.waitKeys()



def GetResponse():
    
    win.flip()

    text = psychopy.visual.TextStim(
        win=win,
        text="Did you detect the change?\nNo : Left Arrow : Yes : Right Arrow",
        color=[-1, -1, -1]
    )

    text.draw()
    win.flip()

    keys = psychopy.event.waitKeys()
    return keys


#info = {'SubId':'1', 'gender':'male', 'age':'10', 'ExpVersion': 1.0}
#StartExp(info)