import spinmob as s
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import voigt_profile
import spify
import os
import x as x_py


class SmartFit:
    """
        Smart overlapping fitting routine

        Init it with the data as a dictionary with angle and count

        Call set_peak_pos with the positions and widths in 2theta degrees

        then Call full fit for fun and profit. Will return the peak pos and uncertainties
    """
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
        # print(window_merge)
        self.merged= window_merge
        self.peaks = peak_pos
        self.extents = window_list

    def full_fit(self, plot=False, fname=None):
        '''
            loops through the windows and simulationsly fits returning all the peaks and unceratinties
        '''

        #create the fitting voigth function we will use
        def fitting_func(x, *argvs):
            #NEWWWW THE first argument is now the strong peak and no longer the average!
            #will do double peak
            wavelengths = [1.540562e-10,1.544390e-10]
            power_ratio = [1, 0.5]
            #set theta_i = arcsin(lambda_i/lambda_1 sin(theta_1))
            #with theta and not 2theta
            thetas = [np.rad2deg(np.arcsin((wavelength/wavelengths[0]) * np.sin(np.deg2rad(argvs[0]/2))))*2 for wavelength in wavelengths]
            #now eval all the peaks
            voigts = [power*argvs[2]* voigt_profile(x - theta, argvs[1], argvs[3], out=None) for theta,power in zip(thetas, power_ratio)]
            #sum them and return!
            return sum(voigts) + np.polyval(argvs[4:], x) 
       
        #init the resutlts
        peaks_res = [[],[]]

        #use the merged windows
        for count, window in enumerate(self.merged):
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
                total_counts = np.zeros_like(y) #fix this cunt
                for i, ind in enumerate(inds_peaks):
                    #convert the full array index to this window index
                    extents = [self.extents[ind][0]-offset, self.extents[ind][1]-offset]
                    #grab the fitting function
                    vals = fitting_func(x[extents[0]:extents[1]], *argvs[i*total_deg: (i+1)*total_deg])
                    total_counts[extents[0]: extents[1]] += vals #superimpose 
                return total_counts
            

            #create the good first guess
            p0 = [[self.deg_pos[i], 0.05, 20, 0.05, 0,0] for i in inds_peaks]
            p0 = [x for i in p0 for x in i]
            
            """
            #ploting debug
            plt.clf()
            y_try = multi_fit(x,*p0)
            plt.scatter(x,y)
            plt.plot(x,y_try)
            plt.show()
            """

            #do the fit and get the reuslts
            popt,pcov = curve_fit(multi_fit,x,y,p0=p0,sigma=n, absolute_sigma=True, maxfev=100000)
            peaks_res[0] += [popt[i] for i in range(len(popt)) if i%total_deg == 0]
            peaks_res[1] += [np.sqrt(np.diag(pcov))[i] for i in range(len(popt)) if i%total_deg == 0]

            #more plotting
            if plot:
                spify.residual_plot(x, y, n,multi_fit, popt, 
                                    r"$2\theta$", "Counts",
                                    "Residuals", 
                                    fname+"_"+str(count)+"_"+"voight")
        return peaks_res

def fcc_anal(peak_data, p0, header):
    # mil_in could be either 'tin' or 'fcc' or others
    miller = np.load('./miller/fcc_indices.npz')['indices']
    # if header == '4_Pb75Sn25_09_09_20':
    #     miller.pop(0)
    #     miller.pop(5)
    # For face centered, 1dim x array suffices in which we combine the hkl with lambda
    x_fit = (miller[:,0][:len(peak_data[0])]**2+\
            miller[:,1][:len(peak_data[0])]**2+\
            miller[:,2][:len(peak_data[0])]**2)*0.1540562e-9**2
    # Fit
    popt,pcov = curve_fit(x_py.shape_func, x_fit, peak_data[0], p0= p0, sigma=peak_data[1])
    unc = np.sqrt(np.diag(pcov))
    chi_sqd = np.sum((peak_data[0] - x_py.shape_func(x_fit, *popt))**2/peak_data[1]**2)/len(peak_data[0]-2)
    # print(chi_sqd)
    # plot
    spify.residual_plot(x_fit, peak_data[0], peak_data[1], x_py.shape_func, [popt[0], popt[1]], 
                    r"$\lambda^2(h^2+k^2+l^2) \ (10^{-19} \ \mathrm{m}^2)$",
                        r"$2\theta \ (^{\circ} )$",
                        "Residuals", 
                        header + "_fcc_lattice", renorm=False)
    # Return the juice
    return popt, unc, chi_sqd


def bcc_anal(peak_data, p0, header):
    # mil_in could be either 'tin' or 'fcc' or others
    miller = np.load('./miller/bcc_indices.npz')['indices']
    # if header == '4_Pb75Sn25_09_09_20':
    #     miller.pop(0)
    #     miller.pop(5)
    # For face centered, 1dim x array suffices in which we combine the hkl with lambda
    x_fit = (miller[:,0][:len(peak_data[0])]**2+\
            miller[:,1][:len(peak_data[0])]**2+\
            miller[:,2][:len(peak_data[0])]**2)*0.1540562e-9**2
    # Fit
    popt,pcov = curve_fit(x_py.shape_func, x_fit, peak_data[0], p0= p0, sigma=peak_data[1])
    unc = np.sqrt(np.diag(pcov))
    chi_sqd = np.sum((peak_data[0] - x_py.shape_func(x_fit, *popt))**2/peak_data[1]**2)/len(peak_data[0]-2)
    # print(chi_sqd)
    # plot
    spify.residual_plot(x_fit, peak_data[0], peak_data[1], x_py.shape_func, [popt[0], popt[1]], 
                    r"$\lambda^2(h^2+k^2+l^2) \ (10^{-19} \ \mathrm{m}^2)$",
                        r"$2\theta \ (^{\circ} )$",
                        "Residuals", 
                        header + "_bcc_lattice", renorm=False)
    # Return the juice
    return popt, unc, chi_sqd


def diamond_anal(peak_data, p0, header):
    # mil_in could be either 'tin' or 'fcc' or others
    miller = np.load('./miller/diamond_indices.npz')['indices']
    if header == 'SC32_02apr':
        miller = np.delete(miller, [2], axis=0)
    # For face centered, 1dim x array suffices in which we combine the hkl with lambda
    x_fit = (miller[:,0][:len(peak_data[0])]**2+\
            miller[:,1][:len(peak_data[0])]**2+\
            miller[:,2][:len(peak_data[0])]**2)*0.1540562e-9**2
    # Fit
    popt,pcov = curve_fit(x_py.shape_func, x_fit, peak_data[0], p0= p0, sigma=peak_data[1])
    unc = np.sqrt(np.diag(pcov))
    chi_sqd = np.sum((peak_data[0] - x_py.shape_func(x_fit, *popt))**2/peak_data[1]**2)/len(peak_data[0]-2)
    # print(chi_sqd)
    # plot
    spify.residual_plot(x_fit, peak_data[0], peak_data[1], x_py.shape_func, [popt[0], popt[1]], 
                    r"$\lambda^2(h^2+k^2+l^2) \ (10^{-19} \ \mathrm{m}^2)$",
                        r"$2\theta \ (^{\circ} )$",
                        "Residuals", 
                        header + "_diamond_lattice", renorm=False)
    # Return the juice
    return popt, unc, chi_sqd

def tin_anal(peak_data, p0, header):
    # mil_in could be either 'tin' or 'fcc' or others
    miller = np.load('./miller/tin_indices.npz')['indices']
    if header == '4_Pb75Sn25_09_09_20':
        miller = np.delete(miller, [0,5], axis=0)
    # For the tetragonal structure, we make the hkl and lambda independant for the fitting function
    # We make the x array a 4-dim made of [h,k,l,lambda]
    x_fit = [miller[:,0][:len(peak_data[0])],
            miller[:,1][:len(peak_data[0])],
            miller[:,2][:len(peak_data[0])],
            np.ones(len(peak_data[0]))*0.1540562e-9]

    # fit the thing
    popt,pcov = curve_fit(x_py.tin_func, x_fit, peak_data[0], p0= p0, sigma=peak_data[1])
    unc = np.sqrt(np.diag(pcov))
    chi_sqd = np.sum((peak_data[0] - x_py.tin_func(x_fit, *popt))**2/peak_data[1]**2)/len(peak_data[0]-2)
    # print(chi_sqd)

    # Plot the thing
    x_plot = 1/4*((x_fit[0]**2+x_fit[1]**2)/popt[0]**2 + x_fit[2]**2/popt[1]**2)*x_fit[3]**2
    spify.residual_plot(x_plot, peak_data[0], peak_data[1], x_py.tin_func_2, [popt[2]], 
                    r"$\lambda^2(h^2+k^2+l^2) \ (10^{-19} \ \mathrm{m}^2)$",
                        r"$2\theta \ (^{\circ} )$",
                        "Residuals", 
                        header + "_tin_lattice", renorm=False)
    # Return the juice
    return popt, unc, chi_sqd


def simple_anal(peak_data, p0, header):
    # mil_in could be either 'tin' or 'fcc' or others
    miller = np.load('./miller/simple_indices.npz')['indices']
    # if header == '4_Pb75Sn25_09_09_20':
    #     miller.pop(0)
    #     miller.pop(5)
    # For face centered, 1dim x array suffices in which we combine the hkl with lambda
    x_fit = (miller[:,0][:len(peak_data[0])]**2+\
            miller[:,1][:len(peak_data[0])]**2+\
            miller[:,2][:len(peak_data[0])]**2)*0.1540562e-9**2
    # Fit
    popt,pcov = curve_fit(x_py.shape_func, x_fit, peak_data[0], p0= p0, sigma=peak_data[1])
    unc = np.sqrt(np.diag(pcov))
    chi_sqd = np.sum((peak_data[0] - x_py.shape_func(x_fit, *popt))**2/peak_data[1]**2)/len(peak_data[0]-2)
    # print(chi_sqd)
    # plot
    spify.residual_plot(x_fit, peak_data[0], peak_data[1], x_py.shape_func, [popt[0], popt[1]], 
                    r"$\lambda^2(h^2+k^2+l^2) \ (10^{-19} \ \mathrm{m}^2)$",
                        r"$2\theta \ (^{\circ} )$",
                        "Residuals", 
                        header + "_simple_lattice", renorm=False)
    # Return the juice
    return popt, unc, chi_sqd



def full_fitter(path_, peak_pos, header, p0, mil_in = 'tin', p02=None, tin_list=None, lead_list = None):
    '''
    This method takes in a data set, and returns lattice parameters

    path_: the path to the UXD file
    peak_pos: the peak positions to help the fitter (can be looked at from plot)
    header: header for saving plots
    p0: guess for the fit
    mil_in: can be tin, fcc or other and decides which fitting function to use among other
    '''
    # Load up the data
    d = s.data.load(path_)
    data = {"angle": d[0], "count": d[1]}
    #Plot it
    plt.figure()
    plt.plot(data["angle"], data["count"])
    # plt.show()
    # return [data["angle"], data["count"]]
    
    # plt.savefig('./figures/{}_full.png'.format(header))

    # Fit all the peaks
    fit = SmartFit(data)
    peak_width = [2.5 for i in range(len(peak_pos))]
    if header == 'Cu3Au_B_14':
        peak_width[7] = 4
        print(peak_width)
    fit.set_peak_pos(peak_pos, peak_width)
    # Recover peak position and uncertainty
    peak_data = np.array(fit.full_fit(plot = True, fname=header))
    # print(peak_data)
    # Function that orders the choices properly using a and c from litterature
    # print(peak_data)
    if mil_in == 'tin':
        popt, unc, chi_sqd = tin_anal(peak_data, p0, header)
        print(popt,unc, chi_sqd)

    elif mil_in=='fcc':
        popt, unc, chi_sqd = fcc_anal(peak_data, p0, header)
        print(popt,unc, chi_sqd)

    elif mil_in=='bcc':
        popt, unc, chi_sqd = bcc_anal(peak_data, p0, header)
        print(popt,unc, chi_sqd)

    elif mil_in=='diamond':
        popt, unc, chi_sqd = diamond_anal(peak_data, p0, header)
        print(popt,unc, chi_sqd)
    
    elif mil_in=='simple':
        popt, unc, chi_sqd = simple_anal(peak_data, p0, header)
        print(popt,unc, chi_sqd)
    
    elif mil_in == 'mix':
        peak_data_lead = [[],[]]
        peak_data_tin = [[],[]]
        for i in range(len(peak_data[0])):
            if i in tin_list:
                peak_data_tin[0].append(peak_data[0][i])
                peak_data_tin[1].append(peak_data[1][i])
            if i in lead_list:
                peak_data_lead[0].append(peak_data[0][i])
                peak_data_lead[1].append(peak_data[1][i])
        # print(peak_data_lead)
        popt_tin, unc_tin, chi_sqd = tin_anal(np.array(peak_data_tin), p02, header)
        print(popt_tin,unc_tin, chi_sqd)
        popt_lead, unc_lead, chi_sqd_fcc = fcc_anal(np.array(peak_data_lead), p0, header)
        print(popt_lead,unc_lead, chi_sqd_fcc)
    
    

        

if __name__=="__main__":
    # Doing Pure tin
    # full_fitter("X-Ray/data/lead_tin_series/1_Sn08-09-20.UXD",
    #             [30.76,32.08,43.86,45.05,55.58,62.61,63.91,64.53],
    #             'tin',[5.83e-10, 3.18e-10, 0])

    # # # Doing Pure Lead
    # full_fitter("X-Ray/data/lead_tin_series/5_Pb_08_09_20D.UXD",
    #             [31.6,36.3,52.5,61.9,65.2,77.1,85.8,88.2,99.4,107.8],
    #             'lead', [4.5e-10, 0], mil_in= 'fcc' )

    # Doing the copper nickel alloy
    copper_bins = [[43.3,50.6,74.5,89.5,95.3],
                [43.5,50.8,74.8,90.5,95.8],
                [43.8,51.3,75.7,91.5,97.0],
                [44,51.1,75.7,92.7,97.7],
                [44.5,52.0,76.2,92.9,98.7]]
    # for i,item in enumerate(os.listdir('./X-Ray/data/copper_nickel_series')):
    #     print(item)
    #     full_fitter('./X-Ray/data/copper_nickel_series/{}'.format(item),
    #                 copper_bins[i],
    #                 item[:-4], [3.5e-10, 0], mil_in= 'fcc' )

    # tin_lead_bins = [[31.4, 32.2,36.44,44.1, 45.0, 52.4, 55.5, 62.3, 64.8, 65.4],
    #                 [30.70,31.4, 32.2,36.44,44.1, 45.0, 52.4, 55.5, 62.3, 64.8, 65.4],
    #                 [30.70,31.4, 32.2,36.44,44.1, 45.0, 52.4, 55.5, 62.3, 64.8, 65.4]]
    # tin_lists = [[1,3,4,6,7,8],[0,2,4,5,7,8,9,10],[0,2,4,5,7,8,9,10]]
    # lead_lists = [[0,2,5,7,8],[1,3,6,8,10],[1,3,6,8,10]]
    # full_fitter('./X-Ray/data/lead_tin_series/4_Pb75Sn25_09_09_20.UXD', tin_lead_bins[0],
    #             '4_Pb75Sn25_09_09_20', [4.94e-10, 0], mil_in= 'mix', 
    #             p02 = [5.83e-10, 3.18e-10, 0], tin_list=tin_lists[0], lead_list=lead_lists[0])
    # full_fitter('./X-Ray/data/lead_tin_series/3_Pb50Sn50_09-09-20.UXD', tin_lead_bins[1],
    #             '3_Pb50Sn50_09-09-20', [4.5e-10, 0], mil_in= 'mix' , p02 = [5.83e-10, 3.18e-10, 0], tin_list=tin_lists[1], lead_list=lead_lists[1])
    # full_fitter('./X-Ray/data/lead_tin_series/2_Pb25Sn75_09-09-20.UXD', tin_lead_bins[2],
    #             '2_Pb25Sn75_09-09-20', [4.5e-10, 0], mil_in= 'mix' , p02 = [5.83e-10, 3.18e-10, 0], tin_list=tin_lists[1], lead_list=lead_lists[2])

    data = []
    legend = []
    for i,item in enumerate(os.listdir('./X-Ray/data/Cu3Au')):
        print(item)
        data.append(full_fitter('./X-Ray/data/Cu3Au/{}'.format(item),
                    copper_bins[i],
                    item[:-4], [3.5e-10, 0], mil_in= 'fcc' ))
        legend.append(item[:4])
    
    plt.figure()
    plt.clf()
    plt.rcParams.update({'font.size': 30})
    for i,item in enumerate(data):
        if i == 0:
            plt.plot(item[0],item[1] + 900*(4-i) + 500, color = 'magenta')
        else:
            plt.plot(item[0],item[1] + 900*(4-i) + 500)
    # miller = np.load('./miller/fcc_indices.npz')['indices']
    # plt.legend(['100% Cu', '75% Cu 25% Ni', '50% Cu 50% Ni', '25% Cu 75% Ni', '100% Ni'], fontsize=24, loc = 'upper left')
    # positions = [44.5-1,52.0,76.2,92.9-1,98.7+0.5]
    # plt.text(25, 200, 'Miller Indices:', fontsize = 28)
    # for i in range(5):
    #     plt.text(positions[i] - 3, 200, '({} {} {})'.format(miller[i][0],miller[i][1],miller[i][2]),
    #      fontsize = 28)
    plt.ylim(0,5000)
    plt.xlim(20,105)
    plt.yticks([])
    plt.xlabel(r'$2\theta \ (^\circ)$')
    plt.ylabel('Counts')
    plt.tight_layout()
    plt.show()


    
