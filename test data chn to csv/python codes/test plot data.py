# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 17:31:41 2021

@author: David
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data_df = pd.read_csv('Ba-133 Calibration 009.csv', skiprows = 6)

x = np.array(data_df['Channel'])
y = np.array(data_df['Counts'])

plt.plot(x, y)



