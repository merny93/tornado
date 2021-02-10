import numpy as np
from uncertainties import ufloat
from uncertainties.umath import cos, sqrt
from scipy import constants as con
import matplotlib.pyplot as plt

def mass_electron(peak_location, a, b, theta):
    '''
    Calculate the predicted mass of the electron from a single peak in a photon spectrum.
    Convert a peak location (in channel number) into energy via a linear relation y=ax+b, 
    where a and b must be provided. The following formula is then used to calculate the 
    mass of the electron:
    m = (1/c^2) * 1/( 1/Ef - 1/E_i ) * (1 - cos(theta)),
    where c = 3*10^8 m/s, E_i = 661.6 keV (energy of the source) and theta is in rad.
    '''
    
    
    Ei = ufloat(661.657e3*con.e, 0.003*con.e) # 661.6 keV converted to Joules
    
    Ef = (a*peak_location + b)*10**3*con.e #if these numbers are wrong thats linda fault
    
    mass = ((1/con.c**2) * 1/( 1/Ef - 1/Ei ) * ( 1 - cos( theta*np.pi/180-np.pi) ))*con.c**2/con.e/1e6
    
    print(1/ (1/Ef - 1/Ei ))
    return mass
    

if __name__ == '__main__':
    plt.rcParams.update({'font.size': 18})
    peaks_data = np.load('../data/tungsten/Angles/peaks.npz')
    plt.clf()
    plt.figure()
    mass = []
    for i in range(len(peaks_data['angles'])):
        line_data = np.load("../data/tungsten/Angles/{}/line_coefs.npz".format(peaks_data['angles'][i]))
        coef, unc = line_data['coefs'], line_data['unc']
        a = ufloat(coef[0], unc[0])
        b = ufloat(coef[1], unc[1])

        temp_mass = mass_electron(ufloat(peaks_data['peaks'][i], peaks_data['uncertainty'][i]), a, b,
                 ufloat(peaks_data['angles'][i], 1.5) )
        temp_mass = dict(zip(["nominal_value", "std_dev"], [temp_mass.nominal_value, temp_mass.std_dev]))
        # print(peaks_data['uncertainty'][i])
        mass.append(temp_mass) 
        print('theta = {} deg -> mass = {} kg'.format( peaks_data['angles'][i], temp_mass))

        plt.scatter(i+1, temp_mass["nominal_value"], marker = 's')
        plt.errorbar(i+1, temp_mass["nominal_value"], yerr = temp_mass["std_dev"], linewidth = 3, capsize=10)

    avg = np.average([x["nominal_value"] for x in mass], weights= [x["std_dev"]**-2 for x in mass])
    # std = np.std([x["nominal_value"] for x in mass])
    std = 1/np.sqrt(np.sum([x["std_dev"]**-2 for x in mass]))
    plt.plot(np.linspace(0,359,360), np.ones(360)*avg, linestyle = '--')
    plt.plot(np.linspace(0,359,360), 
            np.ones(360)*con.physical_constants['electron mass energy equivalent in MeV'][0], 
            linestyle = '--')
    plt.fill_between(np.linspace(0,359,360), np.ones(360)*avg-std,np.ones(360)*avg+std, color = 'red', alpha = 0.3)
    plt.ylim(0.47, 0.56)
    plt.xlim(0,9)
    plt.legend(['Avg', 'Lit',r'$55^\circ$',r'$75^\circ$',r'$95^\circ$',r'$125^\circ$',
                r'$220^\circ$',r'$135^\circ$',r'$230^\circ$',r'$310^\circ$'], fontsize = 12, loc='lower left')
    plt.ylabel(r'Mass (MeV/$c^2$)')
    plt.xticks([])
    plt.tight_layout()
    plt.savefig('../figures/mass.png')
    plt.show()
        
    print('Average is given by: ', avg, "pm", std)
