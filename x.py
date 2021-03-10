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
    return sum([argvs[2+i]* voigt_profile(x - argvs[0+i], argvs[1+i], argvs[3+i], out=None) for i in [0]]) + np.polyval(argvs[4:], x) 
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
windows = [[41.5,45.5],[48,53],[72,77],[88,93], [93,98]]



##loop through the windows and try to fit:

theta_fit= []

for window in windows:
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
    fig = plt.figure()
    gs = fig.add_gridspec(2, hspace = 0, height_ratios = [3,1])
    axs = gs.subplots(sharex= True)

    axs[0].scatter(x,y)
    axs[0].plot(x,fit_func(x,*popt))
    axs[1].scatter(x, y- fit_func(x, *popt))
    plt.show()
    #this gives a really good first guess for the next step
    p0 = [popt[0], 0.5 *popt[1], 2*popt[2], 0.1, 
          0,0]
    
    ##IF U WANT TO PLOT THE FIRST GEUSS
    # plt.scatter(x,y)
    # plt.plot(x,double_fit(x,*p0))
    # plt.show()
    # bounds = ([0, 0,0,0,
    #            popt[0],0,0,0, 
    #            0,0,0,0,
    #            -np.inf, -np.inf],
    #            [popt[0], np.inf,np.inf,np.inf,
    #            np.inf,np.inf,np.inf,np.inf, 
    #            popt[0], np.inf,np.inf,np.inf,
    #            np.inf,np.inf])
    popt,pcov = curve_fit(double_fit, x,y, p0=p0, sigma=noise, absolute_sigma=True, maxfev=1000000) #, bounds=bounds)
    print(popt)
    #now get a prediction for residual calculations!
    y_pred = double_fit(x,*popt)
    
    #compute chisqd
    chi_sqd = np.sum(((y-y_pred)/noise)**2)/x.size
    
    #output
    print("chisqd is about,",  chi_sqd)
    print("position is,", popt[0], "+/-", np.sqrt(np.diag(pcov))[0] )
    
    #save to the results
    #order is critical here!
    theta_fit.append((popt[0], np.sqrt(np.diag(pcov))[0]))
   

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
    axs[0].scatter(x, y, label="data", marker="x", s=25, c="black", linewidth=2)
    axs[0].errorbar(x, y, yerr = noise, linestyle="")
    
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
        ax.set(xlabel=r"$2\theta$ (in degrees)")

    plt.show()

##we fit for theta_fit
## its an list with 2 lists in it. Each of the inner lists constains tupples giving the peak position and 

#the new fitting fit needs to work with different h,k,l
# Make a wrapper that will generate fitting functions with different h,k,l 
def shape_func(lmda_sqd_times_combined_hkl, a):
    return np.rad2deg(np.arcsin(np.sqrt(lmda_sqd_times_combined_hkl/ (4* a**2)))) * 2

hkls = [3,4,8,11,12]
lmdas = sorted([0.1540562e-9])

from itertools import product
x_s = np.array(list(map(lambda x: x[0] * x[1]**2, product(hkls, lmdas))))
y_s = np.array([points[0] for points in theta_fit])
n_s = np.array([points[1] for points in theta_fit])
popt,pcov = curve_fit(shape_func,x_s, y_s, p0= [3.5e-10], sigma=n_s)
print(n_s)
print(popt[0], "+/-", np.sqrt(np.diag(pcov))[0])
print(np.sum((y_s - shape_func(x_s, *popt)/n_s)**2))

x_high = np.linspace(x_s[0], x_s[-1], num=250)
y_high = shape_func(x_high, *popt)
y_pred = shape_func(x_s, *popt)


#get spify
plt.rcParams.update({'font.size': 16})


fig = plt.figure()
gs = fig.add_gridspec(2, hspace = 0, height_ratios = [3,1])
axs = gs.subplots(sharex= True)
#fig.suptitle("Sample fits")

axs[0].plot(x_high,y_high, label="fit", color="magenta")
axs[0].scatter(x_s, y_s, label="data", marker="x", s=25, c="black" , linewidth= 2)

axs[0].legend()
axs[1].scatter(x_s, y_pred-y_s, marker="x", s=25, c="black" , linewidth= 2)
axs[1].axhline(y=0,c="magenta", linestyle="--")
axs[1].errorbar(x_s,  y_pred-y_s, yerr = n_s, linestyle="")
axs[1].set_ylim(-10,10)

#get rid of overlap
to_kills = axs[0].yaxis.get_ticklabels()
to_kills[0].set_visible(False)

#reduce density of residual plot
for n, label in enumerate(axs[1].yaxis.get_ticklabels()):
    if n % 2 != 0:
        label.set_visible(False)

#add the y axis label
axs[0].set(ylabel=r"$\theta$", fontsize = 15)
axs[1].set(ylabel="Residuals")
#do something im not sure what
for ax in axs:
    ax.label_outer()
    ax.set(xlabel=r"$\lambda^2(h^2+k^2+l^2)$", fontsize = 15)

plt.show()