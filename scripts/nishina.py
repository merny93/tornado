import numpy as np 
import matplotlib.pyplot as plt 
import gaussian_fitter as gf 
import file_tools as ft
import peak_finder as pf
from os import path
import uncertainties as un
import scipy.constants as con
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d 
from scipy import integrate

def thomson(theta, a,b,c):
    angle = theta*np.pi/180 - b
    r0 = con.physical_constants['classical electron radius'][0]
    delta_omega = 1 # we will fit for it
    return (1/2)*r0**2*(1 + np.cos(angle)**2)*a + c

def nishina(theta, a,b,c):
    thom = thomson(theta, a,b,c)
    angle = theta*np.pi/180 - b
    Ei = 661.657e3*con.e
    alpha = Ei/(con.m_e*con.c**2)
    k_n = (1+(alpha**2*(1-np.cos(angle))**2)/((1+np.cos(angle)**2)*(1+alpha*(1-np.cos(angle)))))/(1+alpha*(1-np.cos(angle)))**2
    value = a*k_n*thom
    return value

    

if __name__ == '__main__':
    peaks=[]
    uncertainty = []
    bins = [[500,700],[600,800],[700,950],[750,1050],[1200,1550],[1100,1450],[1100,1350],[540,650], [950,1300],[900,1350],[800,1200]]
    angles = [55,75,95,105,220,135,230,310,125,240,250]
    rates = []
    rates_uncertainties = []


    for i in range(len(angles)):
        res = pf.peak_finder(path.join('../data', 'tungsten/Angles/{}/'.format(angles[i])), [bins[i][0], bins[i][1]], plot=False)
        a_0 = un.ufloat(res[0][2], res[1][2])
        sigma = un.ufloat(res[0][1],res[1][1])
        total_time = un.ufloat(30,2)*30
        rate = (a_0*sigma*np.sqrt(2*np.pi)/total_time)
        rates.append(rate.nominal_value)
        rates_uncertainties.append(rate.std_dev)
    
    x = np.linspace(0,360,1000)
    res1 = curve_fit(nishina, angles, rates, p0 = [1e15, np.pi,0], sigma = rates_uncertainties)
    angles_fit = []
    rates_fit = []
    unc_fit = []
    for i in range(len(angles)):
        if np.abs(180-angles[i])<120:
            angles_fit.append(angles[i])
            rates_fit.append(rates[i])
            unc_fit.append(rates_uncertainties[i])
    res2 = curve_fit(thomson, angles_fit, rates_fit, p0 = [1e30, np.pi,0], sigma = unc_fit)
    print(res1[0],res2[0], len(unc_fit))
    

    # plt.clf()
    plt.rcParams.update({'font.size': 15})
    fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw = {'height_ratios': [3,1],'wspace':0, 'hspace':0})
    axs[0].plot(x, nishina(x, *res1[0]), label = 'Klein-Nishina', color='blue')  
    axs[0].plot(x, thomson(x, *res2[0]), label = 'Thomson', color='red')  
    axs[0].scatter(angles,  rates, label = 'data', color='black')
    axs[0].errorbar(angles,  rates, yerr =  rates_uncertainties, xerr = 1.5*np.ones(len(angles)), ls = 'none', color = 'black')
    axs[0].legend(fontsize=14)
    axs[0].set_ylabel('Detection Rate (photons/s)')
    axs[0].legend(['Klein-Nishina','Thomson','Data'])

    axs[1].scatter(angles, rates-nishina(np.array(angles),*res1[0]), marker = 'o', color='blue')
    axs[1].errorbar(angles, rates-nishina(np.array(angles),*res1[0]), yerr = rates_uncertainties, linestyle = "None",capsize=2, color='blue')
    axs[1].scatter(angles, rates-thomson(np.array(angles),*res2[0]), marker = 'o', color='red')
    axs[1].errorbar(angles, rates-thomson(np.array(angles),*res2[0]), yerr = rates_uncertainties, linestyle = "None",capsize=2, color='red')
    axs[1].set_ylabel('Residuals', fontsize=15)
    axs[1].plot(x,np.zeros(len(x)), color='grey', linestyle = '--')
    axs[1].set_yticks([-10,0,10])
    plt.xlim(0, 360)
    plt.xlabel(r'Angle $(^\circ)$')
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.tight_layout()

    plt.savefig('../figures/nishina.png')
    plt.show()
    plt.close()