import csv
from datetime import datetime
import numpy as np
import os
import random



data = []

filename = ""
output_dir = "Results"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

def Init(info):
    
    global filename
    
    filename = 'Sub' + str(info["SubId"]) + "-" + info["gender"] + "-age" + str(info["age"]) + "-ver" + str(info["ExpVersion"])
    filename +=  '-' + datetime.now().strftime("%Y-%m-%d-%H-%M-%S")+".csv"

def StoreData(trial):
    data.append(trial)


def SaveHeader():
    
    with open(output_dir + "/" + filename, 'a', newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["SubId", "ExpType", "trialId", "condition", "response"])

def SaveData(trial):
    
    with open(output_dir + "/" + filename, 'a', newline="") as f:
        writer = csv.writer(f)
        writer.writerow(trial)