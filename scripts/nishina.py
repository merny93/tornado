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

def thomson(theta, a):
    theta = theta*np.pi/180 - np.pi
    r0 = con.physical_constants['classical electron radius'][0]
    delta_omega = 1 # we will fit for it
    return (1/2)*r0**2*(1 + np.cos(theta)**2)*a

def nishina(theta, a):
    thom = thomson(theta, a)
    theta = theta*np.pi/180 - np.pi
    Ei = 661.657e3*con.e
    alpha = Ei/(con.m_e*con.c**2)
    k_n = (1+(alpha**2*(1-np.cos(theta))**2)/((1+np.cos(theta)**2)*(1+alpha*(1-np.cos(theta)))))/(1+alpha*(1-np.cos(theta)))**2
    value = a*k_n*thom
    return value

    

if __name__ == '__main__':
    peaks=[]
    uncertainty = []
    bins = [[500,700],[600,800],[700,950],[750,1050],[1200,1550],[1100,1450],[1100,1350],[540,650]]
    angles = [55,75,95,105,220,135,230,310]
    rates = []
    rates_uncertainties = []

    # Load the absorption effeciency and interpolate it
    data = ft.csv_generic('../calibration_data/Efficiency.csv')
    energy = data['Energy']
    absorption = data['Absorption']
    print(energy)
    detector = interp1d(energy, absorption)

    for i in range(len(angles)):
        res = pf.peak_finder(path.join('../data', 'tungsten/Angles/{}/'.format(angles[i])), [bins[i][0], bins[i][1]], plot=False)
        a_0 = un.ufloat(res[0][2], res[1][2])
        sigma = un.ufloat(res[0][1],res[1][1])
        total_time = un.ufloat(30,2)*30
        rate = (a_0*sigma*np.sqrt(2*np.pi)/total_time) / detector(res[0][0])
        rates.append(rate.nominal_value)
        rates_uncertainties.append(rate.std_dev)
    
    x = np.linspace(0,360,1000)
    res1 = curve_fit(nishina, angles, rates, p0 = [1e30], sigma = rates_uncertainties)
    res2 = curve_fit(thomson, angles, rates, p0 = [1e30], sigma = rates_uncertainties)
    
    plt.figure()
    plt.plot(x, nishina(x, res1[0][0]), label = 'Klein-Nishina')  
    plt.plot(x, thomson(x, res2[0][0]), label = 'Thomson')  
    plt.scatter(angles,  rates, label = 'data')
    plt.errorbar(angles,  rates, yerr =  rates_uncertainties, xerr = 1.5*np.ones(len(angles)), ls = 'none', color = 'b')
    plt.xlim(0, 360)
    plt.xlabel(r'Angle $(^\circ)$')
    plt.ylabel('Detection Rate (photons/s)')
    plt.legend()
    plt.savefig('../figures/nishina.png')
    plt.show()