import numpy as np
import file_tools as ft
import matplotlib.pyplot as plt


##high level plan
##get a basic idea of the map from channels to energy assuming that channel*alpha = energy
##do this with a single fit of one of the calibration runs of sodium as the peak is well defined

##next we will transition into log space where gaussians are simply inverted parabolas

##Now lets establish a noise model in the log space which is done by taking the std of the white noise on small sections
##this is important as we see that not all data is treated equal and thus we must weigh them as such

##Now we will try to find a model that fits more accuratly the channel -> energy.
##to do this we cut out peaks of interest in the data (select them based on the strong peaks we see in literature)
##cut out some data around them and then add this to our model fitting problem
##we will be effectivly constructing a huge model to fit all at once and then only extract the parameters relevant to us
##this will ensure that errors are treated fairly and there are no compounding mistakes .

## critically we need to make a model that takes in a bunch of parameters (of interest is the map from energy to count)
## and outputs the predicted peaks (for all the sets of data) and then we compute chi^2 to minimize with neutons method


#start by fitting a sample of sodium to get a good idea of what to expect:
#code copied from DABID

calib_na = ft.read_csv("../data/calibration/Na-22_Calibration_009.csv")

from scipy.optimize import curve_fit
def gaussian_fit(x, a, mean, sigma):
    return a*np.exp(-(x - mean)**2/(2*sigma**2))

res = curve_fit(gaussian_fit, x := calib_na["Channel"], y := calib_na["Counts"], p0 = [1, 1350, 2]) #look at this walrus
popt = res[0]
plt.clf()
plt.plot(x, y)
plt.plot(x, gaussian_fit(x, *popt))
plt.show()

#popt gives a great starting point since we know that we expected that peak to be at 511Kev
#for the sake of argument we can assume that a*energy = channel and we can calculate this a to be:

inital_scaling = popt[1]/511 #i hate the hardcoded number as much as u do dont hate

## we can also get an estimate of what counts we can expect!
## while this will break all units we can do as follows
##we can convert dose into count by doing something i have yet to figure out

inital_count = 1 ##FIXXXXXX

#this allows us to get a good estimate of the higest and lowest peaks we are willing to fit. 
#and strongest and weakest 

##lets side step and read in the literature:
lit_vals = ft.read_lit("../calibration_data/litterature")


##next we want a model that takes in a bunch of parameters and gives the "prediction"


def model(params): #params is gonna be a list of all the parameters we care about



