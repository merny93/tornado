import numpy as np 
import x_sim
import x 
import spify 
import matplotlib.pyplot as plt
import os

if __name__ == '__main__':
    for i,item in enumerate(os.listdir('./X-Ray/data/mystery')):
        mystery_bins = [[40.55, 58.69, 73.74, 87.76],
                        [ 38.52, 55.60, 69.70, 82.60, 94.95 ],
                        [ 28.55, 47.40, 56.17, 69.17, 76.49, 88.18, 95.08 ],
                        [ 27.35, 45.40, 66.04, 72.96, 83.79 ]]
        print(item)
        if item[0] == 'M':
            x_sim.full_fitter('./X-Ray/data/mystery/{}'.format(item),
                    mystery_bins[i],
                    item[:-10], [4.5e-10, 0], mil_in= 'bcc' )
        elif item[0] == 'S':
            x_sim.full_fitter('./X-Ray/data/mystery/{}'.format(item),
                    mystery_bins[i],
                    item[:-10], [4.5e-10, 0], mil_in= 'diamond' )
        else:
            pass