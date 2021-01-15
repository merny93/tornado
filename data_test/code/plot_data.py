# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 17:31:41 2021
@author: David
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

data_df = pd.read_csv('Na-22_Calibration_009.csv', skiprows = 6)

x = np.array(data_df['Channel'])
y = np.array(data_df['Counts'])

def gaussian_fit(x, a, mean, sigma):
    return a*np.exp(-(x - mean)**2/(2*sigma**2))

popt,pcov = curve_fit(gaussian_fit, x, y, p0 = [1, 1350, 2])
print(pcov)

plt.plot(x, y)
plt.plot(x, gaussian_fit(x, *popt))
plt.xlim(1150,1550)
plt.ylim(0,30)
plt.xlabel('Channel')
plt.ylabel('Counts')
plt.show()
