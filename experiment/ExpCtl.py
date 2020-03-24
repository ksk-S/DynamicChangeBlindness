import psychopy.event
import random

import DataCtl
import ExpOriginal
import ExpBackground
import ExpComposition
import ExpSeparation

import ExpA
import ExpA_prime
import ExpB
import ExpB_prime
import ExpC
import ExpC_prime
import ExpD
import ExpD_prime


num_trials = 20

ScreenSize =[800, 800]

subId = 0

# ExperimentTypes = {
#   "original"                : False,
#   "background"              : False,
#   "compositionRandToRand"   : False,
#   "compositionRandToCircle" : False,
#   "compositionCicleToRand"  : False,
#   "separation"              : True,
# }

ExperimentTypes = {
    # "a"                       : True,
    # "a_prime"                 : True,
    # "b"                       : True,
    "b_prime"                 : True,
    "c"                       : True,
    "c_prime"                 : True,
    "d"                       : True,
    "d_prime"                 : True,
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

    elif exp == 'separation':
        
        ExpSeparation.Init(win, ScreenSize)

    elif exp == "a":

        ExpA.Init(win, ScreenSize)
    
    elif exp == "a_prime":

        ExpA_prime.Init(win, ScreenSize)
    
    elif exp == "b":

        ExpB.Init(win, ScreenSize)
    
    elif exp == "b_prime":

        ExpB_prime.Init(win, ScreenSize)
    
    elif exp == "c":

        ExpC.Init(win, ScreenSize)
    
    elif exp == "c_prime":

        ExpC_prime.Init(win, ScreenSize)
    
    elif exp == "d":

        ExpD.Init(win, ScreenSize)
    
    elif exp == "d_prime":

        ExpD_prime.Init(win, ScreenSize)


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
        
        elif exp == 'separation':
            
            ExpSeparation.StartTrial(condition)
        
        elif exp == "a":
            
            ExpA.StartTrial(condition)
    
        elif exp == "a_prime":

            ExpA_prime.StartTrial(condition)
    
        elif exp == "b":
            
            ExpB.StartTrial(condition)
    
        elif exp == "b_prime":
            
            ExpB_prime.StartTrial(condition)
            
        elif exp == "c":
            
            ExpC.StartTrial(condition)
        
        elif exp == "c_prime":
            
            ExpC_prime.StartTrial(condition)
        
        elif exp == "d":
            
            ExpD.StartTrial(condition)
            
        elif exp == "d_prime":
            
            ExpD_prime.StartTrial(condition)
        
        resKey = GetResponse()
        
        print( 'target' if condition == 0 else 'control' )
        
        response = 1 if resKey[0] == "right" else 0
        
        DataCtl.SaveData([subId, exp, trialId, condition, response])
    
    win.close()
    

def ShowInstruction():
   
    text = psychopy.visual.TextStim(
        win=win,
        pos=(0.0, 0.0),
        text="In this session, you will see a moving object cosists of multipul grating patches. You will need fixate the central dot all the time. Your task is to detect any changes in one of the grating patches.\n\nPress any key to start",
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