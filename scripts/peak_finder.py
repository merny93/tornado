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
    background = nsp.fitter(path.join(sup_path), plot = False)
    corr_data = data-(background*30/12)

    # noise_function = interp1d(x, noise)

    # def denoised_gaus(x_,mean,sigma,a0,a2):
    #     return gf.total_fit(x_,mean,sigma,a0,a2)# + a3*noise_function(x_)
    
    # NOISE UNCERTAINTY
    res = curve_fit(gf.total_fit, x[bins[0]:bins[1]], corr_data[bins[0]:bins[1]], 
            p0 = [np.mean(bins), 10, 10, 1], sigma = np.sqrt(data[bins[0]:bins[1]] +\
                30/12*background[bins[0]:bins[1]]), bounds = (0,1e5))
    popt = res[0]
    unc = np.sqrt(np.diag(res[1]))
    chi_sqd = np.sum((corr_data[bins[0]:bins[1]]-\
        gf.total_fit(x[bins[0]:bins[1]], *popt))**2/(data[bins[0]:bins[1]] +\
            30/12*background[bins[0]:bins[1]])) / (len(corr_data[bins[0]:bins[1]]) - 4)
    
    if plot: 
        plt.clf()
        plt.rcParams.update({'font.size': 18})
        fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw = {'height_ratios': [3,1],'wspace':0, 'hspace':0})

        plt.ylim(0,180)
        plt.xlim(bins[0],bins[1])
        axs[0].scatter(x, corr_data, marker = '.', label = 'data')
        axs[0].plot(x, gf.total_fit(x, *popt), color = 'red', label = 'Gaussian fit')
        max = popt+unc
        min = np.clip(popt-unc,1e-15,None)
        max[0], max[2], min[0], min[2] = popt[0], popt[2], popt[0], popt[2]
        # axs[0].plot(x, gf.total_fit(x, *max), color = 'red',linestyle = '--')
        # axs[0].plot(x, gf.total_fit(x, *min), color = 'red',linestyle = '--')
        axs[0].errorbar(x, corr_data, yerr = np.sqrt(data + 30/12*background), 
                        linestyle = "None",capsize=0)
        axs[0].set_ylabel('Counts')
        axs[0].set_ylim(-20,150)
        axs[0].set_yticks([0,100])
        axs[0].legend(fontsize=14)

        #residual plot
        axs[1].scatter(x, corr_data-gf.total_fit(x, *popt), marker = '.')
        axs[1].errorbar(x, corr_data-gf.total_fit(x, *popt), yerr = np.sqrt(data + 30/12*background), linestyle = "None",capsize=0)
        axs[1].set_ylabel('Residuals') #, position = (0,0))
        plt.xlabel('Channel number')
        axs[1].plot(x,np.zeros(len(x)), color='grey', linestyle = '--')
        # axs[1].fill_between(x, gf.total_fit(x, *popt)-gf.total_fit(x, *max),
        #                      gf.total_fit(x, *popt)-gf.total_fit(x, *min), color = 'red', alpha = 0.5)
        axs[1].set_yticks([-20,0,20])
        axs[1].set_ylim(-50,50)
        plt.subplots_adjust(wspace=0, hspace=0)
        plt.tight_layout()
        # plt.xlim(bins[0],bins[1])

        print(chi_sqd)
        plt.savefig(path.join(sup_path + 'gaussian.png'))
        # plt.show()
        plt.close()
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
    
    
    
    