import matplotlib.pyplot as plt
import numpy as np


def residual_plot(x,data,noise,func,params, xlabel, ylabel1,ylabel2, filename):
    #get spify
    plt.clf()
    plt.rcParams.update({'font.size': 16})


    fig = plt.figure()
    gs = fig.add_gridspec(2, hspace = 0, height_ratios = [3,1])
    axs = gs.subplots(sharex= True)
    #fig.suptitle("Sample fits")

    #generate the high
    x_high = np.linspace(x[0],x[-1], num=250)
    y_high = func(x_high, *params)
    y_pred = func(x,*params)

    axs[0].plot(x_high,y_high, label="fit", color="magenta")
    axs[0].scatter(x, data, label="data", marker="x", s=25, c="black", linewidth=2)
    axs[0].errorbar(x, data, yerr = noise, linestyle="", c="black")
    
    axs[0].legend()
    axs[1].scatter(x, y_pred-data, marker="x", s=25, c="black" , linewidth= 2)
    axs[1].axhline(y=0,c="magenta", linestyle="--")
    
    
    #get rid of overlap
    to_kills = axs[0].yaxis.get_ticklabels()
    to_kills[0].set_visible(False)

    #reduce density of residual plot
    for n, label in enumerate(axs[1].yaxis.get_ticklabels()):
        if n % 2 != 0:
            label.set_visible(False)

    #add the y axis label
    axs[0].set(ylabel=ylabel1)
    axs[1].set(ylabel=ylabel2)
    #do something im not sure what
    for ax in axs:
        ax.label_outer()
        ax.set(xlabel=r"{}".format(xlabel))

    plt.savefig('./figures/{}_voigt.png'.format(filename))


def lattice_alloy_plot(data, litterature, header):
    '''
    Makes a plot comparing the lattice numbers of alloys to litterature numbers
    
    data[0] is an array of the indexing in x of the points
    data[1] is an array of the lattice numbers
    data[2] is an array of the uncertainties on the numbers
    data[3] is an array of the names that should be displayed in the legend
    
    litterature[0] is an array of the indexing in x of the points
    litterature[1] is an array of the litterature lattice numbers
    litterature[2] is an array of the names that should be displayed in the legend

    header is the name to use when saving the file
    '''
    plt.clf()
    plt.figure()
    plt.scatter(data[0], data[1], color = 'black')
    plt.scatter(litterature[0], litterature[1], color = 'blue')
    plt.legend(data[3]+litterature[2])
    plt.errorbar(data[0], data[1], yerr=data[2], color = 'black', linestyle = "")
    plt.xlabel('% Nickel')
    plt.ylabel('Lattice Number')
    plt.savefig('./figures/{}_lattice_alloy.png'.format(header))