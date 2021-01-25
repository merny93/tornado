import numpy as np
import file_tools as ft
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import erfc
from matplotlib import rc


def collapse_data(data_list):
    # load the data and collapses the runs in one big array
    from functools import reduce
    return reduce(lambda x,y: x+y["Counts"], data_list, np.zeros_like(data_list[0]["Counts"]))


def total_fit(x,mean,sigma,beta,a0,a1,a2):
    gaus = a0*np.exp(-(x - mean)**2/(2*sigma**2))
    expo = a1*np.exp((x-mean)/beta)*erfc((x-mean)/(np.sqrt(2)*sigma)+sigma/(np.sqrt(2)*beta))
    step = a2*erfc((x-mean)/(np.sqrt(2)*sigma))
    return gaus + expo + step

def plotter(num, key, guess_mean, guess_width = 2, guess_height=1, plot = True):
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
    res = curve_fit(total_fit, x, y, p0 = [guess_mean, guess_width, 100,guess_height,1,1], sigma = uncert, bounds = (0,1e5)) 
    # recover the parameter values
    popt = res[0]
    # The uncertainty of the fit
    uncertainty = np.sqrt(np.diag(res[1]))
    # Chi Squared
    chi_sqd = np.sum((y-total_fit(x, *popt))**2/np.sqrt(y))
    # Plotting the fit
    if plot:
        # Make 2 plots, one for the fit and one for the residuals
        plt.rcParams.update({'font.size': 18})
        fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw = {'height_ratios': [3,1],'wspace':0, 'hspace':0})
        # fig.suptitle('{}'.format(key))

        # Scatter of the data and gaussian line
        axs[0].scatter(x, y, marker = '.', label = 'data')
        axs[0].plot(x, total_fit(x, *popt), color = 'orange', label = 'Gaussian fit')
        max = popt+uncertainty
        min = np.clip(popt-uncertainty,1e-15,None)
        max[0], max[2], min[0], min[2] = popt[0], popt[2], popt[0], popt[2]
        axs[0].plot(x, total_fit(x, *max), color = 'orange',linestyle = '--')
        axs[0].plot(x, total_fit(x, *min), color = 'orange',linestyle = '--')
        axs[0].errorbar(x, y, yerr = uncert, linestyle = "None")
        axs[0].set_ylabel('Counts')
        axs[0].legend()

        #residual plot
        axs[1].scatter(x, y-total_fit(x, *popt), marker = '.')
        axs[1].errorbar(x, y-total_fit(x, *popt), yerr = uncert, linestyle = "None")
        axs[1].set_ylabel('Residuals') #, position = (0,0))
        axs[1].plot(x,np.zeros(len(x)), color='grey', linestyle = '--')
        axs[1].fill_between(x, total_fit(x, *popt)-total_fit(x, *max), total_fit(x, *popt)-total_fit(x, *min), color = 'orange', alpha = 0.5)
        plt.xlabel('Channel number')
        plt.subplots_adjust(wspace=0, hspace=0)
        plt.tight_layout()
        plt.savefig('../figures/{}_gaussian.pdf'.format(key))
        # plt.show()
        plt.close()
    # print(chi_sqd)
    return [popt[0], uncertainty[1]]

def reverse_line(x,a,b):
    return (x-b)/a

def line(x,a,b):
    return a*x+b

def line_fit(points_y, litterature, plot = True):
    '''
    Takes in
    points_y: calibration data and uncertainty
    litterature: litterature data
    plot: output a plot file
    ''' 
    res = curve_fit(reverse_line, litterature[0], points_y[:,0],  p0 = [0.29, -25], sigma = points_y[:,1])
    popt = res[0]
    uncertainty = np.sqrt(np.diag(res[1]))
    chi_sqd = np.sum((points_y[:,0]-reverse_line(litterature[0], *popt))**2/points_y[:,1])
    if plot:
        # plot the line
        fig, axs = plt.subplots(2,1, sharex=True, sharey=False, gridspec_kw = {'height_ratios': [3,1],'wspace':0, 'hspace':0})

        max = uncertainty + popt
        min = -uncertainty + popt
        axs[0].scatter(litterature[0], points_y[:,0],marker = '.', label = 'data')
        axs[0].errorbar(litterature[0],points_y[:,0],  yerr = points_y[:,1], linestyle = "None")
        axs[0].plot(line(points_y[:,0], *popt), points_y[:,0], color = 'orange', label = 'Linear Fit')
        n = ['Ba-133','Co-57','Na-22', 'Cs-137']
        positions = [(litterature[0][1]+10, points_y[:,0][1]-50),(litterature[0][3]+20,points_y[:,0][3]-30),
                     (litterature[0][2]-20, points_y[:,0][2]-100),(litterature[0][0], points_y[:,0][0]-50)]
        for i, txt in enumerate(n):
            axs[0].annotate(txt, positions[i], fontsize=10)
        axs[0].set_ylabel('Channel Number')
        axs[0].legend()
        #plot the residuals
        linspace = np.linspace(np.min(litterature[0]),np.max(litterature[0]),1000)
        axs[1].scatter(litterature[0],points_y[:,0]-reverse_line(litterature[0], *popt),  marker = '.')
        axs[1].errorbar(litterature[0],points_y[:,0]-reverse_line(litterature[0], *popt),  yerr = points_y[:,1], linestyle = "None")
        axs[1].fill_between(linspace, reverse_line(linspace, *popt)-reverse_line(linspace, *max),
                             reverse_line(linspace, *popt)-reverse_line(linspace, *min), color = 'orange', alpha = 0.5)
        axs[1].plot(linspace,np.zeros(len(linspace)), color='grey', linestyle = '--')
        axs[1].set_ylabel('Residuals')
        plt.xlabel('Energy (keV)')

        plt.tight_layout()
        plt.savefig('../figures/line.pdf')
        plt.show()
    print(popt, uncertainty, chi_sqd)

if __name__ == "__main__":
    elements = {}
    for element in ft.SOURCE_NAMES:
        elements[element] = ft.get_data("../data/calibration", source_name = element)
    line_points = []
    # line_points.append(plotter([1200,1450], 'Na-22', guess_mean = 1360))
    line_points.append(plotter([900,1050], 'Ba-133', guess_mean = 970))
    # line_points.append(plotter([1550,1910], 'Cs-137', 1742, guess_width = 58, guess_height=77))
    # line_points.append(plotter([315,380], 'Co-57', 350,  guess_width = 2, guess_height=1600))
    # line_fit(np.array(line_points), [[511.0, 356.0129, 661.657, 122.06065],[5, 7, 3, 12]])