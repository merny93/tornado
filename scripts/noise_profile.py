import file_tools as ft
import gaussian_fitter as gf
from os import path
import numpy as np 
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit 

def polynome(x,a0,a1,a2,a3,a4,a5,a6,a7,a8):
    y = a0+a1*x+ a2*x**2+a3*x**3+a4*x**4+a5*x**5+a6*x**6+a7*x**7+a8*x**8
    return y

def fitter(sup_path):
    '''
    Performs a calibration on a data set of size (4x2x2048)
    and returns the parameters a,b of the line fit and their uncertainties
    '''
    data = gf.collapse_data(ft.get_data(path.join(sup_path, "01_No_Scatterer")))[220:]
    x = np.linspace(220,2047,2048-220)
    uncert = np.sqrt(data)
    print(data)
    res = curve_fit(polynome, x, data, p0 = [1,1,1,1,1,1,1,1,1]) 
    popt = res[0]
    uncertainty = np.sqrt(np.diag(res[1]))
    print(popt, uncertainty)
    plt.figure()
    plt.scatter(x,data, marker = '.', alpha = 0.1)
    plt.plot(x, polynome(x, *popt))
    plt.savefig(path.join(sup_path + 'background_noise.png'))
    np.savez(path.join(sup_path + 'background_noise.npz'), coefs = popt)


fitter('../data/tungsten/Angles/55/')
