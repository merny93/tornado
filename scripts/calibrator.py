import sys 
import gaussian_fitter as g_f 
import matplotlib.pyplot as plt 
import numpy as np
import file_tools as ft
from os import path

def calibrator(sup_path, plotting = False):
    '''
    Reads in the path of a data set. This path folder must contain folders called
    00_Other_Sources
    and
    04_Other_Sources
    calibrates with both folder and then returns the params and uncertainties.
    Can make a plot of the lines.
    '''
    elements1, elements2 = {},{}
    # for element in ft.SOURCE_NAMES:
    print(path.join(sup_path,"00_Other_Sources", 'Cs-137'))
    elements1['Cs-137'] = ft.get_data(path.join(sup_path,"00_Other_Sources", 'Cs-137'), source_name = 'Cs-137')
    elements2['Cs-137'] = ft.get_data(path.join(sup_path, "04_Other_Sources", 'Cs-137'), source_name = 'Cs-137')
    # params1, unc1 = g_f.calibrator_fit(elements1)
    # params2, unc2 = g_f.calibrator_fit(elements2)

    if plotting: 
        plt.clf() 
        plt.figure()
        x = np.linspace(0,2047,2048)
        y1 = g_f.collapse_data(elements1['Cs-137'])
        y2 = g_f.collapse_data(elements2['Cs-137'])
        plt.scatter(x, y1, marker = '.', color = 'blue', label='Before', alpha=0.1)
        plt.scatter(x, y2, marker = '.', color = 'red', label='After', alpha=0.1)
        plt.savefig(path.join(sup_path,'calibration.png'))

    # np.savez(path.join(sup_path,'calibration.npz',element), params1 = params1, params2 = params2, unc1 = unc1, unc2 = unc2)
    # return params1,unc1,params2,unc2

if __name__ == '__main__':
    for angle in [55,75,95,105,135,220,230,310]:
        calibrator('../data/tungsten/Angles/{}'.format(angle), plotting=True)

