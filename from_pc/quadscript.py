# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 20:59:20 2021

@author: remotelab
"""

q1  = d[2]+1j*d[3]             # Complex phasor for Channel 1
q2  = d[4]+1j*d[5]             # Complex phasor for Channel 2
q1r = q1*abs(q2)/q2            # q1 but with phase relative to that of q2
x   = ( d[1] )                 # Frequency for the x-data
y   = ( abs(q1r), angle(q1r) ) # Scaled magnitude and phase for Channel 1

xlabels = 'Frequency (Hz)'
ylabels = ( 'Magnitude', 'Relative Phase' )