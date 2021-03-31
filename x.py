import spinmob as s
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import voigt_profile
import spify
import os
#load pure copper data!

#the new fitting fit needs to work with different h,k,l
# Make a wrapper that will generate fitting functions with different h,k,l 
def shape_func(lmda_sqd_times_combined_hkl, a, theta_0):
    return np.rad2deg(np.arcsin(np.sqrt(lmda_sqd_times_combined_hkl/ (4* a**2)) )) * 2 + 2*theta_0

# Same but for tin
# x is [h,k,l,lambda]
def tin_func(x, a, c, theta_0):
    return np.rad2deg(np.arcsin(np.sqrt(1/4*((x[0]**2+x[1]**2)/a**2 + x[2]**2/c**2)*x[3]**2))) * 2 + 2*theta_0

def tin_func_2(x, theta_0):
    return np.rad2deg(np.arcsin(np.sqrt(x))) * 2 + 2*theta_0

#GAUSSIAN FIRST TRY
def fit_func(x, *argvs):
    #start with gaussian
    #args is x then [mean, sigma,  a, c]
    return argvs[2] * np.exp(-((x-argvs[0])/argvs[1])**2) + argvs[3]

##DOUGBLE PEAK MIT FIT
def double_fit(x, *argvs):
    #now what if we tried something spicy ;)
    #double 
    # args is x then [mean1, sigma1, a1, gamma1, mean2, sigma2, a2, gamma2, *poly]
    return sum([argvs[2+i]* voigt_profile(x - argvs[0+i], argvs[1+i], argvs[3+i], out=None) for i in [0]]) + np.polyval(argvs[4:], x) 
    #that line sums over the two voigt profiles and then adds a polynomial as discribed by mit

def line(x, a,b):
    return a*x+b



def fitter(filename, windows, path_, hkls, p0_fitter):
    d = s.data.load(os.path.join(path_,"{}.UXD".format(filename)))
    data = {"angle": d[0], "count": d[1]}


    #plot the data for some good vibes and easy checks
    plt.plot(data["angle"], data["count"])
    plt.show()


    ##DEBUGING OF FITTING FUNCTION
    ''' 
    x = np.linspace(-2,2)
    y = double_fit(x, *[-1, 0.5, 250, 0.7,1, 0.5, 0, 0.7, 0,0,0])
    plt.plot(x,y)
    plt.show()
    '''


    #these windows are where the peaks are in angle units (NOT ARRAY INDEX)
    # factor shifts the window to the right because of alloys




    ##loop through the windows and try to fit:

    theta_fit= []
    peaks_data = {}
    for i,window in enumerate(windows):
        #window is in angle so lets get it in index
        pos = np.argwhere(np.logical_and(data["angle"]>window[0], data["angle"]<window[1]))
        x = np.array(data["angle"][pos]).flatten()
        y = np.array(data["count"][pos]).flatten()
        #noise as sqrt of count but dont let it get too small
        noise = np.sqrt(y)
        noise[noise < 2] = 2

        #optimize***************
        #first get a easy run-through using the gaussian
        popt,pcov = curve_fit(fit_func, x,y, p0=[sum(window)/2,1, 250, 0.7], sigma=noise, absolute_sigma=True)

        ##PLOT THE FIRST FIT WITH RESIDUALS
        '''
        fig = plt.figure()
        gs = fig.add_gridspec(2, hspace = 0, height_ratios = [3,1])
        axs = gs.subplots(sharex= True)
        axs[0].scatter(x,y)
        axs[0].plot(x,fit_func(x,*popt))
        axs[1].scatter(x, y- fit_func(x, *popt))
        plt.show()
        '''

        #this gives a really good first guess for the next step
        p0 = [popt[0], 0.5 *popt[1], 2*popt[2], 0.1, 
            0,0]
        
        ##IF U WANT TO PLOT THE FIRST GEUSS
        # plt.scatter(x,y)
        # plt.plot(x,double_fit(x,*p0))
        # plt.show()


        ##WAS BEING USED BEFORE TO IMPROVE FIT
        '''
        bounds = ([0, 0,0,0,
                   popt[0],0,0,0, 
                   0,0,0,0,
                   -np.inf, -np.inf],
                   [popt[0], np.inf,np.inf,np.inf,
                   np.inf,np.inf,np.inf,np.inf, 
                   popt[0], np.inf,np.inf,np.inf,
                   np.inf,np.inf])
        '''

        popt, pcov = curve_fit(double_fit, x,y, p0=p0, sigma=noise, absolute_sigma=True, maxfev=1000000) #, bounds=bounds)
        # print(popt)
        #now get a prediction for residual calculations!
        y_pred = double_fit(x,*popt)
        
        #compute chisqd
        chi_sqd = np.sum(((y-y_pred)/noise)**2)/x.size
        
        #output
        # print("chisqd is about,",  chi_sqd)
        # print("position is,", popt[0], "+/-", np.sqrt(np.diag(pcov))[0] )
        
        #save to the results
        #order is critical here!
        theta_fit.append((popt[0], np.sqrt(np.diag(pcov))[0]))

        #get spify
        spify.residual_plot(x,y,noise, double_fit, popt,
                            r"$2\theta \ (^\circ)$","Counts","Residuals", 
                            filename+"_"+str(i+1)+"_voigt")
        
        peaks_data['peak_{}'.format(i)] = {}
        peaks_data['peak_{}'.format(i)]['popt'] = popt[0]
        peaks_data['peak_{}'.format(i)]['unc'] = np.sqrt(np.diag(pcov))[0]
        peaks_data['peak_{}'.format(i)]['chi'] = chi_sqd

    ##we fit for theta_fit
    ## its an list with 2 lists in it. Each of the inner lists constains tupples giving the peak position and 

    lmdas = sorted([0.1540562e-9])
    ##Create the x_s by combining all the hkls and lambdas
    from itertools import product
    x_s = np.array(list(map(lambda x: x[0] * x[1]**2, product(hkls, lmdas))))
    #collaps the data
    y_s = np.array([points[0] for points in theta_fit])
    n_s = np.array([points[1] for points in theta_fit])

    #fit for the parameter a
    popt,pcov = curve_fit(shape_func,x_s, y_s, p0= p0_fitter, sigma=n_s)

    chi_sqd = np.sum((y_s - shape_func(x_s, *popt))**2/n_s**2)/len(y_s-2)

    x_high = np.linspace(x_s[0], x_s[-1], num=250)
    y_high = shape_func(x_high, *popt)
    y_pred = shape_func(x_s, *popt)


    spify.residual_plot(x_s, y_s, n_s,shape_func, popt, 
                        r"$\lambda^2(h^2+k^2+l^2) \ (10^{-19} \ \mathrm{m}^2)$",
                         r"$2\theta \ (^{\circ} )$",
                         "Residuals", 
                         filename + "_lattice", renorm=True)
    peaks_data['line_popt'] = popt
    peaks_data['line_unc'] = np.sqrt(np.diag(pcov))
    peaks_data['line_chi'] = chi_sqd
    return peaks_data

def full_anal(path_, bins, hkls, print_=False, line_plot=False, p0_fitter = [3.5e-10, 0]):
    filenames = os.listdir(path_)
    filenames = list(map(lambda x: x.split(".")[0], filenames))
    # filenames = ["Cu_03_09_20", "Cu75Ni25","Cu50Ni50","Cu25Ni75", "Ni_03_09_20"]
    total_data = {}
    for i in range(len(filenames)):
        total_data['{}'.format(filenames[i])] = (fitter(filenames[i], bins[i], path_,hkls, p0_fitter))
    # print(total_data)
    if print_:
        try:
            import pprint 
            pp = pprint.PrettyPrinter(indent=4)
            pp.pprint(total_data)
        except:
            print(total_data)
    
    if line_plot:
        y_data = [[],[]]
        for name in filenames:
            y_data[0].append(total_data[name]["line_popt"][0])
            y_data[1].append(total_data[name]["line_unc"][0])
        
        y_data[0] = np.array(y_data[0], dtype = float)
        y_data[1] = np.array(y_data[1], dtype = float)
        # plot_litterature = [[0,100],
        #                     [361.49e-12, 352.4e-12],
        #                     ['Litterature']]
        # spify.lattice_alloy_plot(plot_data, plot_litterature, 'Cu-Ni')
        x = np.array([0,75,50,25, 100] ,dtype = float)
        popt,pcov = curve_fit(line,x,y_data[0], p0= [1,1], sigma=y_data[1])

        chi_sqd = np.sum((y_data[0] -line(np.array(x), *popt))**2/y_data[1]**2)/len(y_data[0]-2)

        spify.residual_plot(x, y_data[0], y_data[1],
                            line, popt,"Nickel Concentration", "Lattice Parameter","Residuals", 
                            path_[-15:] + "_vegard", legend_loc = "upper right")

if __name__ == '__main__':
    copper_bins = [[[41.5,45.5],[48,53],[72,77],[88,93], [93,98]],
            [[41.5,49],[48,56],[72,80],[88,95], [95,101]],
            [[41.5,48],[48,56],[72,80],[88,95], [94,100]],
            [[41.5,47],[48,56],[72,80],[87,94], [95,99]],
            [[41.5,50],[48,56],[72,80],[90,96], [97,101]]]
    copper_nickel_hkls = [3,4,8,11,12]
    # full_anal('X-Ray/data/copper_nickel_series', copper_bins, copper_nickel_hkls, print_=True)
    lead_bins = [[[30,35],[36,38],[50,56],[60,65],[63,67],[74,82],[82,87],[87,90],[96,103],[104,112]]]
    lead_hkls = [3,4,8,11,12,16,19,20,24,27]
    full_anal('X-Ray/data/pb', lead_bins, lead_hkls , print_=True, p0_fitter = [1e-9, 0])

    


