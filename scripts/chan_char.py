import numpy as np
import file_tools as ft
import matplotlib.pyplot as plt


##high level plan
##get a basic idea of the map from channels to energy assuming that channel*alpha = energy
##do this with a single fit of one of the calibration runs of sodium as the peak is well defined


##Now lets establish a noise model which is done by taking the std of the white noise on small sections
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
x = calib_na["Channel"]
y = calib_na["Counts"]
res = curve_fit(gaussian_fit, x, y, p0 = [1, 1350, 2]) #look at this walrus
popt = res[0]
# plt.clf()
# plt.plot(x, y)
# plt.plot(x, gaussian_fit(x, *popt))
# plt.show()

#popt gives a great starting point since we know that we expected that peak to be at 511Kev
#for the sake of argument we can assume that a*energy = channel and we can calculate this a to be:

inital_scaling = popt[1]/511 #i hate the hardcoded number as much as u do dont hate

## we can also get an estimate of what counts we can expect!
## while this will break all units we can do as follows
##we can convert dose into count by doing something i have yet to figure out

inital_count = 1 ##FIXXXXXX

#this allows us to get a good estimate of the higest and lowest peaks we are willing to fit. 
#and strongest and weakest 

##now we pick out the noise floor as follows:

window=30
#switch to noise later
calib_data = ft.read_csv("../data/calibration/Na-22_Calibration_009.csv")
calib_noise = calib_data["Counts"].copy()
niter = 150 #number of iterations (this is a guess but works ok)
for _ in range(niter):
    #I started with a python loop implementation and it took forever so sorry for the hard code 
    # I want to make a rolling window array. So 3 wide such that the first row is 1st ,2nd ,3rd element 
    #second row is 2nd, 3rd, 4th element and so on so forth 
    #numpy does not have a funciton to do this so lets tell it to read the array differently with strides!

    stride = (calib_noise.strides[0], calib_noise.strides[0])
    #this line tells python that a move to the right or a move to down is the same equivalent 
    #since the elemnt to the right and the element bellow are both indexed 1 away in the original array

    calib_rolling = np.lib.stride_tricks.as_strided(calib_noise, shape = (calib_noise.size - 2, 3), strides=stride, writeable=False)
    ##that line generated the array representation (not actually writable)

    calib_argmin = np.argmin(calib_rolling, axis=-1)#check which is the minimum 

    calib_reset = np.where(calib_argmin == 1, (calib_rolling[:,0] + calib_rolling[:,2])/2, calib_rolling[:,1])
    #another hard line. This one checks if the min happened in the middle and if that is the case fills an array
    #of size ft_noise.size -2 with the average of the points around or just copies the value if it wasnt the min

    calib_noise[1:-1] = calib_reset #finally adding it to the array for the next itteration


padded_noise = np.pad(calib_noise,(window//2, window//2), mode="edge")
array_view = np.lib.stride_tricks.as_strided(padded_noise, shape=(calib_noise.size,window), strides=[padded_noise.strides[0] for i in range(2)])
noise_model =np.mean(array_view, axis=-1)
noise_model = np.maximum(np.ones(noise_model.size), noise_model)
noise_model = np.sqrt(noise_model)
print(np.min(noise_model))

##we also want to do a running average
# plt.clf()
# plt.plot(calib_data["Counts"])
# plt.plot(noise_model)

# plt.show()
##lets side step and read in the literature:


##next we want a model that takes in a bunch of parameters and gives the "prediction"


def gaussian_model(x, mean, a, sigma):
    return a*np.exp(-(x - mean)**2/(2*sigma**2))

def model(params): #params is gonna be a list of all the parameters we care about
    #params is [a,b, alpha_1, sigma_1, alpha_2, etc]
    params_alphas = list(params[2::2])
    params_beta = list(params[3::2])
    lit_vals = ft.read_lit("../calibration_data/litterature")
    prediction = {}
    for key in lit_vals: #key=='Ba-133'
        peaks = lit_vals[key] #list of dictionaries
        prediction[key] = {"mask":np.zeros(2048, dtype=bool), "Counts": np.zeros(2048)}
        for peak in peaks:
            #peak is dictionary with important value being "energy"
            mean_energy = peak["energy"]
            #using initial fit to decide where to plot
            central_channel = int(np.round(inital_scaling*mean_energy))
            #width to fit in channels
            num_chan = 200
            if (num_chan//2 + central_channel) > 2000:
                continue
            if (central_channel -num_chan//2)< 0:
                continue
            #generate the channels
            channels = np.arange(central_channel-(num_chan//2), central_channel + (num_chan//2))
            #model which goes from params to prediction
            counts = gaussian_model(params[0]*channels + params[1], mean_energy, params_alphas.pop(0), params_beta.pop(0)) 
            # where_add = np.where(prediction[key]["mask"][central_channel-(num_chan//2): central_channel + (num_chan//2)] == False)[0]
            # where_add_full = where_add + central_channel-(num_chan//2)
            # prediction[key]["Counts"][where_add_full] = counts[where_add]
            prediction[key]["Counts"][central_channel-(num_chan//2): central_channel + (num_chan//2)] += counts
            # plt.clf()
            # plt.plot(counts)
            # plt.show()
            prediction[key]["mask"][central_channel-(num_chan//2): central_channel + (num_chan//2)] = True

    return prediction

## prediction = {'Ba-133': {"mask": nparraybools, "Counts": nparrayfloats}}


def model_dif(modelA, modelB):
    md = []
    for key in modelA:
        mask = np.where(modelA[key]["mask"]==True)[0]
        md.append(modelA[key]["Counts"][mask] - modelB[key]["Counts"][mask])
    return np.concatenate(md)

def get_resid(prediction):
    ##take prediction from the model and data and compute chisqud
    res  = []
    for key in prediction:
        #key is name of sample
        data_sets = ft.get_data("../data/calibration", source_name=key)
        full_counts = []
        for data_set in data_sets:
            #we compute chi^2
            full_counts.append(data_set["Counts"])
        data_counts = np.mean(full_counts, axis=0)
        mask = np.where(prediction[key]["mask"]==True)[0]
        #noise_model from above
        diff = (data_counts[mask] - prediction[key]["Counts"][mask])
        res.append(diff)
    return np.concatenate(res)

def get_noise(for_model):
    # lit_vals = ft.read_lit("../calibration_data/litterature")
    ##incorporate the uncertainties later?

    md = []
    for key in for_model:
        mask = np.where(for_model[key]["mask"]==True)[0]
        md.append(noise_model[mask])
    return np.concatenate(md)

    

def chi_sqd(prediction, do_plot = False):
    ##take prediction from the model and data and compute chisqud
    total_chi = 0
    for key in prediction:
        #key is name of sample
        data_sets = ft.get_data("../data/calibration", source_name=key)
        # print(data_sets)
        for data_set in data_sets:
            #we compute chi^2
            data_counts = data_set["Counts"]
            mask = np.where(prediction[key]["mask"]==True)[0]
            #noise_model from above
            diff = (prediction[key]["Counts"] - data_counts)**2
            chi_vec = (diff/noise_model)[mask]
            chi_cur = np.sum(chi_vec)
            total_chi += chi_cur
        if do_plot:
            plt.clf()
            plt.plot(data_counts)
            plt.plot(prediction[key]["Counts"])
            plt.savefig("../data_test/calibration_sources/fit_example"+key+".png")
    return total_chi


lit = ft.read_lit("../calibration_data/litterature")
from functools import reduce
num_par = sum([len(lit[key]) for key in lit])
print(num_par)
init_guess = [1/inital_scaling, 0] + [popt[0] if x%2==0 else popt[2] for x in range(2*num_par)]
import time as t
# print(init_guess)

t1 = t.time()
res_model = model(init_guess)
my_chi = chi_sqd(res_model)
diff = model_dif(res_model, res_model)
print(diff.shape)
resid = get_resid(res_model)
print(resid.shape)
t2 = t.time()

print(t2-t1)

import newton_solver as ns
print(init_guess)

params, cov = ns.newton_solve(model, get_resid, get_noise, chi_sqd, model_dif, init_guess)
print(params[:2], init_guess[:2])

chi_sqd(model(params), do_plot=True)