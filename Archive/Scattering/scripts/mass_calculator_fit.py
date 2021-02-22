import numpy as np
from uncertainties import ufloat
from uncertainties.umath import cos, sqrt
from scipy import constants as con
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def mass_electron(theta, m, g):
    '''
    Fit Compton's formula for the mass of the electron. If we solve for Ef in Compton's formula,
    we get Ef = 1/( (1/(mc^2))*(1-cos(theta)) + 1/Ei ). We fit y = Ef and x = theta, and 
    we know Ei and work in units such that c^2 disappears. 
    '''
    Ei = 661.657 # keV
    Ef = 1/( (1/(m*1e3))*(1 - np.cos(theta*np.pi/180 - g))  + 1/Ei ) 
    # c^2 cancels with the c^2 from [m]=MeV/c^2.
    # factor of 1e3 is to convert MeV to keV.
    
    return Ef


def fitter():
    plt.rcParams.update({'font.size': 18})
    peaks_data = np.load('../data/tungsten/Angles/peaks.npz')
    plt.clf()
    plt.figure()
    mass = []
    data = []
    data_uncertainties = []
    for i in range(len(peaks_data['angles'])):
        line_data = np.load("../data/tungsten/Angles/{}/line_coefs.npz".format(peaks_data['angles'][i]))
        coef, unc = line_data['coefs'], line_data['unc']
        a = ufloat(coef[0], unc[0])
        b = ufloat(coef[1], unc[1])
        theta = peaks_data['angles'][i]

        value = a*ufloat(peaks_data['peaks'][i], peaks_data['uncertainty'][i]) + b
        # print(a,b)
        data.append(value.n)
        # print(value.n, value.std_dev)
        data_uncertainties.append(value.std_dev)
        
    res = curve_fit(mass_electron, peaks_data['angles'], data, p0 = [0.511, np.pi], sigma = data_uncertainties, absolute_sigma=True)
    popt = res[0]
    unc = np.sqrt(np.diag(res[1]))
    chi_sqd = np.sum((data-mass_electron(peaks_data['angles'], *popt))**2/(data_uncertainties)) / (len(data) - 2)
    linspace = np.linspace(0,360,1000)

    
    # print(popt, unc)
    mass = ufloat(popt[0],unc[0])
    print(mass, popt, chi_sqd)

    plt.clf()
    plt.rcParams.update({'font.size': 15})
    fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw = {'height_ratios': [3,1],'wspace':0, 'hspace':0})
    axs[0].scatter(peaks_data['angles'],data, marker = '.', label = 'data', color = 'black')
    axs[0].plot(linspace, mass_electron(linspace, *popt), label = 'Compton\'s fit', color='blue')

    axs[0].errorbar(peaks_data['angles'],data, yerr = data_uncertainties, linestyle = "None",capsize=0, color='black')
    axs[0].legend(fontsize=14)

    #residual plot
    axs[1].scatter(peaks_data['angles'], data-mass_electron(peaks_data['angles'], *popt), marker = '.', color='blue')
    axs[1].errorbar(peaks_data['angles'], data-mass_electron(peaks_data['angles'], *popt), yerr = data_uncertainties, linestyle = "None",capsize=0, color='blue')
    axs[1].set_ylabel('Residuals', fontsize=15)
    axs[1].plot(linspace,np.zeros(len(linspace)), color='grey', linestyle = '--')
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.tight_layout()

    plt.xlabel(r'$\theta \ (^\circ)$', fontsize = 14)
    axs[0].set_ylabel(r'$E_f$ (keV)', fontsize = 14)

    plt.savefig('../figures/mass_fit.png')
    plt.show()
    plt.close()

if __name__ == '__main__':
    fitter()
