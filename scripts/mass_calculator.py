import file_tools as ft
import noise_profile as nsp
import gaussian_fitter as gf
from os import path
import sys
import numpy as np 
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit 

def mass_calculator(sup_path, plot = True):
    '''
    Loads rod data from specified path, substracts the noise profile and plot.
    The formula used to calculate the mass is the following:
    m = (1/c^2) * 1/( 1/(energy[channel]*10^3*1.602*10^{-19}) - 1/E_i*10^3*1.602*10^{-19} ) * (1-cos(theta)),
    where c = 3*10^8 m/s, E_i = 661.6 keV.
    '''
    data = gf.collapse_data(ft.get_data(path.join(sup_path, "02_Tungsten")))
    x = np.linspace(0,2047,2048)
    uncert_data = np.sqrt(data)
    noise_coefs = nsp.fitter(path.join('../data/', sys.argv[1]))
    noise = nsp.polynome(x, *noise_coefs)
    # NOISE UNCERTAINTY
    data_corrected = data[220:] - noise[220:]
    
    if plot: 
        plt.figure()
        plt.scatter(x[220:],data_corrected, marker = '.', alpha = 1)
        plt.show()
    
    


if __name__ == '__main__':
    mass_calculator(path.join('../data', sys.argv[1]))