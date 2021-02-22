import numpy as np
import file_tools as ft
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


# Load the absorption efficiency and interpolate it
data = ft.csv_generic('../calibration_data/Efficiency.csv')
energy = data['Energy']
absorption = data['Absorption']
detector = interp1d(energy, absorption, kind='cubic')

a,b = [],[]
for angle in [55,75,95,105,220,135,230,310]:
    data = np.load('../data/tungsten/Angles/{}/line_coefs.npz'.format(angle))
    a.append(data['coefs'][0])
    b.append(data['coefs'][1])

channels = np.linspace(0,2047,2048)
absorption = []
for i in range(len(channels)):
    if np.mean(a)*channels[i]+np.mean(b) < 227:
        absorption.append(float(detector(227)/100))
    else:
        absorption.append(float(detector(np.mean(a)*channels[i]+np.mean(b))/100))
np.savez('../data/absorption.npz', absorption = absorption)

plt.figure()
plt.plot(channels,absorption)
plt.show()