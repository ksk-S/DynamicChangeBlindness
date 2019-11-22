#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Dynamic ChangeBlindness Experiment.
"""

import ExpCtl

#from __future__ import absolute_import, division, print_function

from psychopy import gui

# create a DlgFromDict
info = {'SubId':'1', 'gender':['male', 'female'], 'age':'', 'ExpVersion': 1.0}

infoDlg = gui.DlgFromDict(dictionary=info, title='TestExperiment', order=['SubId', 'gender', 'age'], fixed=['ExpVersion']) 

if infoDlg.OK:  # this will be True (user hit OK) or False (cancelled)
    ExpCtl.StartExp(info)
    
else:
    print('User Cancelled')

