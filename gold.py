import numpy as np 
import x_sim
import x 
import spify 
import matplotlib.pyplot as plt
import os

if __name__ == '__main__':
    for i,item in enumerate(os.listdir('./X-Ray/data/Cu3Au')):
        gold_bins = [ [41.75, 48.61, 71.11, 85.99, 90.81],[23.80, 33.89, 41.82, 48.6, 54.84, 60.62, 71.22, 76.23, 81.2, 86.07, 90.91, 95.71]]
        print(item)
        if i==0:
            x_sim.full_fitter('./X-Ray/data/Cu3Au/{}'.format(item),
                    gold_bins[i],
                    item[:-10], [4.5e-10, 0], mil_in= 'fcc' )
        else:
            x_sim.full_fitter('./X-Ray/data/Cu3Au/{}'.format(item),
                    gold_bins[i],
                    item[:-10], [4.5e-10, 0], mil_in= 'simple' )