import psychopy.event
import random

import DataCtl
import ExpOriginal
import ExpBackground
import ExpComposition

num_trials = 10

ScreenSize =[800, 800]

subId = 0

ExperimentTypes = {
  "original"        : True,
  "background"      : True,
  "composition"     : True,
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
       
    elif exp == 'composition':
        
        ExpComposition.Init(win, ScreenSize)

        
    ShowInstruction()
    
    
    for trialId in range(num_trials):
        condition = trial_list[trialId]
        
        if exp =='original':
    
            ExpOriginal.StartTrial(condition)
        
        elif exp == 'background':
            
            ExpBackground.StartTrial(condition)
        
        elif exp == 'composition':
            
            ExpComposition.StartTrial(condition)
        
        resKey = GetResponse()
        
        response = 1 if resKey[0] == "right" else 0
        
        DataCtl.SaveData([subId, exp, trialId, condition, response])
    
    win.close()
    

def ShowInstruction():
   
    text = psychopy.visual.TextStim(
        win=win,
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