import numpy as np


def mass_electron(peak_location, a, b, theta):
    '''
    Calculate the predicted mass of the electron from a single peak in a photon spectrum.
    Convert a peak location (in channel number) into energy via a linear relation y=ax+b, 
    where a and b must be provided. The following formula is then used to calculate the 
    mass of the electron:
    m = (1/c^2) * 1/( 1/Ef - 1/E_i ) * (1 - cos(theta)),
    where c = 3*10^8 m/s, E_i = 661.6 keV (energy of the source) and theta is in rad.
    '''
    
    c = 3e8
    Ei = 661.6e3*1.602e-19 # 661.6 keV converted to Joules
    
    Ef = ( a*peak_location + b )*10**3*1.602e-19
    
    mass = (1/c**2) * 1/( 1/Ef - 1/Ei ) * ( 1 - np.cos( theta*np.pi/180) )
    
    return mass
    
    
peaks_data = np.load('../data/tungsten/Angles/peaks.npz')
a = 0.3856 
b = -14.9

mass = []
for i in range(len(peaks_data['angles'])):
    mass.append( mass_electron(peaks_data['peaks'][i], a, b, peaks_data['angles'][i]) )
    print('theta = {} deg -> mass = {} kg'.format( peaks_data['angles'][i], round(mass[i], 35) ))
    
print('Standard deviation of masses = ', round(np.std(mass), 35) )

    
