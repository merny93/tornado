# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 17:31:41 2021

@author: David
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

std=[]
for i in range(9):
    data = pd.read_csv('../csv/Na-22_Calibration_00{}.csv'.format(i), skiprows = 6)
    std.append([])
    for j in range(0,2037):
        std[i].append(np.std(data['Counts'][j:j+10]))
    for j in range(10):
        std[i].append(std[i][2036])
plt.figure()
plt.plot(np.linspace(0,2047,2047), np.mean(np.array(std), axis =0))
plt.show()