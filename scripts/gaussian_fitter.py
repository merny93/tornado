import numpy as np
import file_tools as ft
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def collapse_data(data_list):
    # load the data and collapses the runs in one big array
    from functools import reduce
    return reduce(lambda x,y: x+y["Counts"], data_list, np.zeros_like(data_list[0]["Counts"]))

def gaussian_fit(x, a, mean, sigma, c,d):
    return a*np.exp(-(x - mean)**2/(2*sigma**2)) + c +d*x

def plotter(num, key, guess_mean, guess_width = 2, guess_height=100, plot = False):
    '''
    Takes in
    num: limits of the x axis to fit
    key: element to fit
    guesses: guess for the gaussian fit
    plot: output a plot file
    '''
    x = np.linspace(0,2047,2048)[num[0]:num[1]] # channel numbers
    y = collapse_data(elements[key])[num[0]:num[1]] # counts
    uncert = np.sqrt(y) # uncertainty of counting
    #fit the gaussian
    res = curve_fit(gaussian_fit, x, y, p0 = [guess_height, guess_mean, guess_width, 1,1], sigma = uncert) 
    # recover the parameter values
    popt = res[0]
    # The uncertainty of the fit
    uncertainty = np.sqrt(np.diag(res[1]))
    # Chi Squared
    chi_sqd = np.sum((y-gaussian_fit(x, *popt))**2/np.sqrt(gaussian_fit(x, *popt)))
    # Plotting the fit
    if plot:
        plt.clf() #clear the figure
        # Make 2 plots, one for the fit and one for the residuals
        fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw = {'hspace': 0})
        fig.suptitle('{}'.format(key))

        # Scatter of the data and gaussian line
        axs[0].scatter(x, y, marker = '.', label = 'data')
        axs[0].plot(x, gaussian_fit(x, *popt), color = 'orange', label = 'Gaussian fit')
        axs[0].errorbar(x, y, yerr = uncert, linestyle = "None")
        axs[0].set_ylabel('Counts', position = (0,0))
        axs[0].legend()

        #residual plot
        axs[1].scatter(x, y-gaussian_fit(x, *popt), marker = '.')
        axs[1].errorbar(x, y-gaussian_fit(x, *popt), yerr = uncert, linestyle = "None")
        plt.xlabel('Channel number')
        plt.savefig('../figures/{}_gaussian.png'.format(key))
    return [popt[1], uncertainty[1]]

def line(x,a,b):
    return a*x+b

def line_fit(points_y, litterature,plot = True):
    '''
    Takes in
    points_y: calibration data and uncertainty
    litterature: litterature data
    plot: output a plot file
    ''' 
    res = curve_fit(line, points_y[:,0], litterature[0],  p0 = [0.29, -25], sigma = litterature[1])
    popt = res[0]
    uncertainty = np.sqrt(np.diag(res[1]))
    chi_sqd = np.sum((litterature[0]-line(points_y[:,0], *popt))**2/litterature[1])
    if plot:
        # plot the line
        plt.clf()
        fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw = {'hspace': 0})
        fig.suptitle('Line Fit')
        axs[0].scatter(points_y[:,0], litterature[0], marker = '.', label = 'data')
        axs[0].plot(points_y[:,0], line(points_y[:,0], *popt), color = 'orange', label = 'Linear Fit')
        axs[0].errorbar(points_y[:,0], line(points_y[:,0], *popt), yerr = litterature[1], linestyle = "None")
        axs[0].set_ylabel('Energy', position = (0,0))
        axs[0].legend()
        #plot the residuals
        axs[1].scatter(points_y[:,0], litterature[0]-line(points_y[:,0], *popt), marker = '.')
        axs[1].errorbar(points_y[:,0], litterature[0]-line(points_y[:,0], *popt), yerr = litterature[1], linestyle = "None")
        plt.xlabel('Channel Number')
        plt.savefig('../figures/line.png')
    print(popt, uncertainty, chi_sqd)

if __name__ == "__main__":
    elements = {}
    for element in ft.SOURCE_NAMES:
        elements[element] = ft.get_data("../data/calibration", source_name = element)
    line_points = []
    line_points.append(plotter(num = [1200,1550], key = 'Na-22', guess_mean = 1360))
    line_points.append(plotter(num = [900,1200], key = 'Ba-133', guess_mean = 970))
    line_points.append(plotter([1550,1910], 'Cs-137', 1742, guess_width = 58, guess_height=77))
    line_points.append(plotter([315,400], 'Co-57', 350,  guess_width = 2, guess_height=1600))
    line_fit(np.array(line_points), [[511.0, 356.0129, 661.657, 122.06065],[5, 7, 3, 12]])