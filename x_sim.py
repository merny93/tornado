import spinmob as s
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import voigt_profile
import spify
import os


class SmartFit:
    def __init__(self, data):
        #data is a dict of numpy arrays:
        #"angle" give the angles
        #"count" gives the counts
        self.d = data

    def set_peak_pos(self, peak_pos, peak_width):
        """
        We want to create a set of windows where our peaks overlap so we can fit simulatniously there
        """
        
        #we want in degrees for fitting later!
        self.deg_pos = peak_pos
        
        #convert into the array index instead of degrees
        peak_pos = [int(np.argwhere(self.d["angle"] > i)[0]) for i in peak_pos]
        bin_size = self.d["angle"][1] - self.d["angle"][0]
        peak_width = [int(i/bin_size) for i in peak_width]
        
        ##we want to create windows that we will be fitting on!
        window_list = [[pos-width//2, pos + width//2] for pos, width in zip(peak_pos,peak_width)]

        #now the windows might overlap so we want to combine
        window_list.sort(key=lambda x: x[0]) #should be totally not needed
        from functools import reduce
        window_merge = reduce(lambda x,y: x + [y] if y[0] > x[-1][1] else x[:-1] + [[x[-1][0] , max(x[-1][1], y[1])]], window_list, [window_list[0]])
        #i think this works :)
        print(window_merge)
        self.merged= window_merge
        self.peaks = peak_pos
        self.extents = window_list

    def full_fit(self):
        '''
            loops through the windows and simulationsly fits returning all the peaks and unceratinties
        '''

        #create the fitting voigth function we will use
        def fitting_func(x, *argvs):
            return argvs[2]* voigt_profile(x - argvs[0], argvs[1], argvs[3], out=None) + np.polyval(argvs[4:], x) 
       
        #init the resutlts
        peaks_res = [[],[]]

        #use the merged windows
        for window in self.merged:
            ##this is one fitting section
            #lets get the extents and the centers

            #this grabs the indecies of peaks involved
            inds_peaks = [i for i,x in enumerate(self.peaks) if (x>window[0] and x<window[1]) ]
            N = len(inds_peaks)
            poly_deg = 2
            fit_deg = 4
            total_deg = poly_deg + fit_deg

            ##get the relavant data
            x = self.d["angle"][window[0]:window[1]] 
            y = self.d["count"][window[0]:window[1]] 
            n = np.sqrt(y)
            n[n<3] = 3
            #this is how far we gotta scale
            offset = window[0]

            #generate the full fitting function
            def multi_fit(x, *argvs):
                total_counts = np.zeros_like(y)
                for i, ind in enumerate(inds_peaks):
                    #convert the full array index to this window index
                    extents = [self.extents[ind][0]-offset, self.extents[ind][1]-offset]
                    #grab the fitting function
                    vals = fitting_func(x[extents[0]:extents[1]], *argvs[i*total_deg: (i+1)*total_deg])
                    total_counts[extents[0]: extents[1]] += vals #superimpose 
                return total_counts
            

            #create the good first guess
            p0 = [[self.deg_pos[i], 0.1, 30, 0.1, 0,0] for i in inds_peaks]
            p0 = [x for i in p0 for x in i]
            
            #ploting debug
            '''
            y_try = multi_fit(x,*p0)
            plt.scatter(x,y)
            plt.plot(x,y_try)
            plt.show()
            '''

            #do the fit and get the reuslts
            popt,pcov = curve_fit(multi_fit,x,y,p0=p0,sigma=n, absolute_sigma=True, maxfev=100000)
            peaks_res[0] += [popt[i] for i in range(len(popt)) if i%total_deg == 0]
            peaks_res[1] += [np.sqrt(np.diag(pcov))[i] for i in range(len(popt)) if i%total_deg == 0]

            #more plotting
            '''
            goddem = multi_fit(x,*popt)
            plt.scatter(x,y)
            plt.plot(x,goddem)
            plt.show()
            '''
        return peaks_res



if __name__=="__main__":
    d = s.data.load("X-Ray/data/lead_tin_series/Pb25Sn75_09-09-20.UXD")
    data = {"angle": d[0], "count": d[1]}
    plt.plot(data["angle"], data["count"])
    plt.show()
    fit = SmartFit(data)
    peak_pos = [30.7,31.3,32.1, 36.5, 36.3,44,45,52.5,55.4]
    peak_width = [2 for i in range(len(peak_pos))]
    fit.set_peak_pos(peak_pos, peak_width)
    # plt.plot(data["count"])
    # plt.show()
    fit.full_fit()