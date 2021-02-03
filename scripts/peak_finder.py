import file_tools as ft
import noise_profile as nsp
import gaussian_fitter as gf
from os import path
import numpy as np 
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit 
from scipy.interpolate import interp1d


def peak_finder(sup_path, bins, plot = True):
    '''
    Loads rod data from specified path, substracts the noise profile and plot.
    '''
    data = gf.collapse_data(ft.get_data(path.join(sup_path, "02_Tungsten")))
    x = np.linspace(0,2047,2048)
    background = nsp.fitter(path.join(sup_path), plot = True)
    corr_data = data-(background*30/12)

    # noise_function = interp1d(x, noise)

    # def denoised_gaus(x_,mean,sigma,a0,a2):
    #     return gf.total_fit(x_,mean,sigma,a0,a2)# + a3*noise_function(x_)
    
    # NOISE UNCERTAINTY
    res = curve_fit(gf.total_fit, x[bins[0]:bins[1]], corr_data[bins[0]:bins[1]], 
            p0 = [np.mean(bins), 10, 10, 1], sigma = np.sqrt(data[bins[0]:bins[1]] + 30/12*background[bins[0]:bins[1]]), bounds = (0,1e5))
    popt = res[0]
    unc = np.sqrt(np.diag(res[1]))
    
    if plot: 
        plt.figure()
        plt.scatter(x, corr_data, marker = '.', alpha = 1)
        plt.plot(x, gf.total_fit(x, *popt))
        plt.ylim(0,180)
        # plt.ylim(min(data[bins[0]:bins[1]]), max(data[bins[0]:bins[1]]))
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
    
    
    
    