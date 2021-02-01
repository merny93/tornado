import file_tools as ft
import noise_profile as nsp
import gaussian_fitter as gf
from os import path
import numpy as np 
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit 

def peak_finder(sup_path, bins, plot = True):
    '''
    Loads rod data from specified path, substracts the noise profile and plot.
    '''
    data = gf.collapse_data(ft.get_data(path.join(sup_path, "02_Tungsten")))[bins[0]:bins[1]]
    x = np.linspace(0,2047,2048)[bins[0]:bins[1]]
    noise_coefs = nsp.fitter(path.join(sup_path), plot = False)
    noise = nsp.polynome(x, *noise_coefs)
    # NOISE UNCERTAINTY
    data_corrected = data - noise
    new_x = x[(data_corrected > 0)]
    data_corrected = data_corrected[(data_corrected > 0)]
    data_uncertainty = np.sqrt(data_corrected)
    res = curve_fit(gf.total_fit, new_x, data_corrected, 
            p0 = [np.mean(bins), 10, 100, 10, 1, 1], sigma = data_uncertainty, bounds = (0,1e5))
    popt = res[0]
    unc = np.sqrt(np.diag(res[1]))
    
    if plot: 
        plt.figure()
        plt.scatter(new_x,data_corrected, marker = '.', alpha = 1)
        plt.plot(x, gf.total_fit(x, *popt))
        plt.xlim(bins[0],bins[1])
        # plt.show()
        plt.savefig(path.join(sup_path + 'gaussian.png'))
    return popt, unc
    
    


if __name__ == '__main__':
    peaks=[]
    uncertainty = []
    bins = [[500,700],[600,800],[700,950],[750,1050],[1200,1550]]
    angles = [55,75,95,105,220]
    for i in range(len(angles)):
        res = peak_finder(path.join('../data', 'tungsten/Angles/{}/'.format(angles[i])), [bins[i][0], bins[i][1]])
        peaks.append(res[0][0])
        uncertainty.append(res[1][0])
    np.savez(path.join('../data', 'tungsten/Angles/peaks.npz'), angles = angles, peaks = peaks, uncertainty = uncertainty)
    
    
    
    