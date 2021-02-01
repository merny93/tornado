import numpy as np
from uncertainties import ufloat
from uncertainties.umath import cos, sqrt
from scipy import constants as con

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
    
    mass = (1/con.c**2) * 1/( 1/Ef - 1/Ei ) * ( 1 - cos( theta*np.pi/180) )
    
    return mass
    

if __name__ == '__main__':
    peaks_data = np.load('../data/tungsten/Angles/peaks.npz')
    for i in range(len(peaks_data['angles'])):
        line_data = np.load("../data/tungsten/Angles/{}/line_coefs.npz".format(peaks_data['angles'][i]))
        coef, unc = line_data['coefs'], line_data['unc']
        a = ufloat(coef[0], unc[0])
        b = ufloat(coef[1], unc[1])
        mass = []
        
        mass.append(mass_electron(ufloat(peaks_data['peaks'][i], peaks_data['uncertainty'][i]), a, b, ufloat(peaks_data['angles'][i], 0.1) ))
        print(i)
        print(len(mass))
        print('theta = {} deg -> mass = {} kg'.format( peaks_data['angles'][i], mass[i] ))
        
    print('Average is given by: ', np.average([x.nominal_value for x in mass], weights = [x.std_dev for x in mass]), "pm", 1/np.sqrt(np.sum([x.std_dev for x in mass])))