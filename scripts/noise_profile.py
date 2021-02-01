import file_tools as ft
import gaussian_fitter as gf
from os import path
import sys
import numpy as np 
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit 

def polynome(*argv):
    y = sum([argv[0]**i * argv[i+1] for i in range(len(argv)-1)])
    return y

def fitter(sup_path, plot = True):
    '''
    Performs a calibration on a data set of size (4x2x2048)
    and returns the parameters a,b of the line fit and their uncertainties
    '''
    data = gf.collapse_data(ft.get_data(path.join(sup_path, "01_No_Scatterer")))[220:]
    
    padded_noise = np.pad(calib_noise,(window//2, window//2), mode="edge")
    array_view = np.lib.stride_tricks.as_strided(padded_noise, shape=(calib_noise.size,window), strides=[padded_noise.strides[0] for i in range(2)])
    noise_model =np.mean(array_view, axis=-1)
    noise_model = np.maximum(np.ones(noise_model.size), noise_model)
    noise_model = np.sqrt(noise_model)
    print(np.min(noise_model))
    x = np.linspace(220,2047,2048-220)
    uncert = np.sqrt(data)
    print(data)
    res = curve_fit(polynome, x, data, p0 = [-3.20194279e+00,-3.38704699e-02,1.48655540e-02,-1.15786819e-04,
                                            3.68163557e-07, -5.20736574e-10,  1.64862236e-14,  1.13040755e-15,
                                            -1.95524233e-18,  1.75095304e-21, -9.49046487e-25,  3.13656100e-28,
                                            -5.83843311e-32,  4.70517366e-36]) 
    popt = res[0]
    uncertainty = np.sqrt(np.diag(res[1]))
    print(popt, uncertainty)
    if plot:
        plt.figure()
        plt.scatter(x,data, marker = '.', alpha = 0.1)
        plt.plot(x, polynome(x, *popt))
        plt.savefig(path.join(sup_path + 'background_noise.png'))
    np.savez(path.join(sup_path + 'background_noise.npz'), coefs = popt)
    return popt

if __name__ == '__main__':
    fitter(path.join('../data/', sys.argv[1]))