import spinmob as s
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import voigt_profile
#load pure copper data!
d = s.data.load("X-Ray/data/copper_nickel_series/Cu_03_09_20.UXD")
data = {"angle": d[0], "count": d[1]}


#GAUSSIAN FIRST TRY
def fit_func(x, *argvs):
    #start with gaussian
    #args is x then [mean, sigma,  a, c]
    return np.array(argvs[2] * np.exp(-((x-argvs[0])/argvs[1])**2) + argvs[3], dtype=float)

##DOUGBLE PEAK MIT FIT
def double_fit(x, *argvs):
    #now what if we tried something spicy ;)
    #double 
    # args is x then [mean1, sigma1, a1, gamma1, mean2, sigma2, a2, gamma2, *poly]
    return sum([argvs[2+i]* voigt_profile(x - argvs[0+i], argvs[1+i], argvs[3+i], out=None) for i in [0,4]]) + np.polyval(argvs[8:], x) 
    #that line sums over the two voigt profiles and then adds a polynomial as discribed by mit


#plot the data for some good vibes and easy checks
# plt.plot(data["angle"], data["count"])
# plt.show()

##DEBUGING OF FITTING FUNCTION 
# x = np.linspace(-2,2)
# y = double_fit(x, *[-1, 0.5, 250, 0.7,1, 0.5, 0, 0.7, 0,0,0])
# plt.plot(x,y)
# plt.show()


#these windows are where the peaks are in angle units (NOT ARRAY INDEX)
windows = [[41,46],[48,53],[72,77],[88,93],[93,98]]



##loop through the windows and try to fit:
for window in windows:
    #window is in angle so lets get it in index
    pos = np.argwhere(np.logical_and(data["angle"]>window[0], data["angle"]<window[1]))
    x = np.array(data["angle"][pos]).flatten()
    y = np.array(data["count"][pos]).flatten()
    #noise as sqrt of count but dont let it get too small
    noise = np.sqrt(y)
    noise[noise < 5] = 5

    #optimize***************
    #first get a easy run-through using the gaussian
    popt,pcov = curve_fit(fit_func, x,y, p0=[sum(window)/2,1, 250, 0.7], sigma=noise, absolute_sigma=True)

    #this gives a really good first guess for the next step
    popt,pcov = curve_fit(double_fit, x,y, p0=[popt[0] + 0.1,popt[1], popt[2]/2, 0,popt[0] - 0.1,popt[1], popt[2]/2, 0, 0,0,0], sigma=noise, absolute_sigma=True, maxfev=100000)
    #now get a prediction for residual calculations!
    y_pred = double_fit(x,*popt)
    
    #compute chisqd
    chi_sqd = np.sum(((y-y_pred)/noise)**2)/x.size
    
    #output
    print("chisqd is about,",  chi_sqd)
    print("position is,", popt[0], "+/-", np.sqrt(np.diag(pcov))[0] )
    
    #for plot get high res
    x_high = np.linspace(x[0], x[-1], num=250)
    y_high = double_fit(x_high, *popt)

    #get spify
    plt.rcParams.update({'font.size': 16})


    fig = plt.figure()
    gs = fig.add_gridspec(2, hspace = 0, height_ratios = [3,1])
    axs = gs.subplots(sharex= True)
    #fig.suptitle("Sample fits")

    axs[0].plot(x_high,y_high, label="fit", color="magenta")
    axs[0].scatter(x, y, label="data", marker="x", s=25, c="black" , linewidth= 2)
    
    axs[0].legend()
    axs[1].scatter(x, y_pred-y, marker="x", s=25, c="black" , linewidth= 2)
    axs[1].axhline(y=0,c="magenta", linestyle="--")
    
    
    #get rid of overlap
    to_kills = axs[0].yaxis.get_ticklabels()
    to_kills[0].set_visible(False)

    #reduce density of residual plot
    for n, label in enumerate(axs[1].yaxis.get_ticklabels()):
        if n % 2 != 0:
            label.set_visible(False)

    #add the y axis label
    axs[0].set(ylabel="Counts")
    axs[1].set(ylabel="Residuals")
    #do something im not sure what
    for ax in axs:
        ax.label_outer()
        ax.set(xlabel=r"Diffraction Scattering Angle ($2\theta$)")

    plt.show()

