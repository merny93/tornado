import numpy as np
from uncertainties import ufloat
from uncertainties.umath import cos, sqrt
from scipy import constants as con
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def mass_electron(theta, m, g):
    '''
    Calculate the predicted mass of the electron from a single peak in a photon spectrum.
    Convert a peak location (in channel number) into energy via a linear relation y=ax+b, 
    where a and b must be provided. The following formula is then used to calculate the 
    mass of the electron:
    m = (1/c^2) * 1/( 1/Ef - 1/E_i ) * (1 - cos(theta)),
    where c = 3*10^8 m/s, E_i = 661.6 keV (energy of the source) and theta is in rad.
    '''

    Ei = 661.657e3*con.e
    Ef = 1/(((1-np.cos(theta*np.pi/180 - g))/(m/5.609588357e29)/con.c**2) + 1/Ei)/1e3/con.e

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

        value = (a*ufloat(peaks_data['peaks'][i], peaks_data['uncertainty'][i]) + b)
        data.append(value.n)
        data_uncertainties.append(value.std_dev)

    res = curve_fit(mass_electron, peaks_data['angles'], data, p0=[0.511,np.pi], sigma = data_uncertainties)
    popt = res[0]
    unc = np.sqrt(np.diag(res[1]))
    linspace = np.linspace(0,360,1000)

    
    # print(popt, unc)
    mass = ufloat(popt[0],unc[0])
    print(mass)

    plt.clf()
    plt.rcParams.update({'font.size': 18})
    fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw = {'height_ratios': [3,1],'wspace':0, 'hspace':0})
    axs[0].scatter(peaks_data['angles'],data, marker = '.', label = 'data')
    axs[0].plot(linspace, mass_electron(linspace, *popt), color = 'red', label = 'Gaussian fit')

    axs[0].errorbar(peaks_data['angles'],data, yerr = data_uncertainties, linestyle = "None",capsize=0)
    axs[0].legend(fontsize=14)

    #residual plot
    axs[1].scatter(peaks_data['angles'], data-mass_electron(peaks_data['angles'], *popt), marker = '.')
    axs[1].errorbar(peaks_data['angles'], data-mass_electron(peaks_data['angles'], *popt), yerr = data_uncertainties, linestyle = "None",capsize=0)
    axs[1].set_ylabel('Residuals')
    axs[1].plot(linspace,np.zeros(len(linspace)), color='grey', linestyle = '--')
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.tight_layout()

    plt.xlabel(r'$\theta \ (^\circ)$', fontsize = 14)
    plt.ylabel(r'$E_f$', fontsize = 14)

    plt.savefig('../figures/mass_fit.png')
    plt.show()
    plt.close()

if __name__ == '__main__':
    fitter()
