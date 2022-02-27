## This script reads in 2D histogram data in csv.format,
## i.e. number of samples per bin, where number of bins is given by number of rows
## times number of cols.The rows give bins in range (-xrange, xrange)+bias_x, 
## the cols bins in range (-yrange, yrange)+bias_y. xrange, yrange and the biasses, 
## as well as temperature needs to be specified below.
## Read in data, normalize histogram, print histogram and true density, and print
## distributional error vs h (if activated).

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import LogLocator
import csv
from cycler import cycler
# from integrators_tests import dist_diff
from GM_functions import posterior_plot
from GM_functions import gaussmix
from matplotlib import cm

# plt.close('all')
# plt.rcParams.update({'font.size': 33})
# plt.rc('legend', fontsize=33)    # legend fontsize
# plt.rcParams['axes.grid'] = True
# plt.rc('lines', linewidth=3)
# plt.rcParams['axes.titlepad'] = 25

#%%
 
x_range = 1  # configuration interval (-x_range, x_range) where densities are notably > 0
bias_x = 5
y_range = 1
bias_y = 0

T=1
h="0.025"
L=1
# rho = config_density_HO

plotrange = (-0.07,0.07)
name = "/home/rene/PhD/Research/code/Integrators/GM/histo_OBABO_"   # file name head to read
title = r"  Harmonic Oscillator Sampling with OBABO"

files = ["h0.050"]
labels = [r"$h$=0.05"]


# colors = ["b", "orange",  "g", "c"]
colors = ["c"]

dist_mse = []
fig, ax = plt.subplots()

for (file, label, c) in zip(files,labels, colors):
    print(name+file)
    with open(name+file) as csv_file:

        csv_reader = csv.reader(csv_file, delimiter=' ')
        # rowct = 0
        # for row in csv_reader:  # count rows
        #     rowct += 1
        # bin_ctrs = np.zeros((rowct,rowct))
        bin_ctrs = []
        for row in csv_reader:
            bin_ctrs += [[int(i) for i in row[0:-1]]]
    
    nr_bins = len(bin_ctrs) - 2
    bin_ctrs = np.array(bin_ctrs)
    # dist_mse += [dist_diff(bin_ctrs, nr_bins, (-1*x_range, x_range), rho, (T))]    

    # plot histograms
    delta_x = 2*x_range / nr_bins  # width of bins
    delta_y = 2*y_range / nr_bins
    histo = bin_ctrs / (np.sum(bin_ctrs)*delta_x*delta_y)  # normalize
    # midx = np.arange(-x_range+bias_x + 0.5*delta_x, x_range+bias_x, delta_x)  # center of bins in x
    # midy = np.arange(-y_range+bias_y + 0.5*delta_y, y_range+bias_y, delta_y)  # center of bins in y 
    # ax.plot(midx, histo[1:-1], c, label=label )
    plt.imshow(histo[1:-1,1:-1], cmap=cm.jet, extent=(bias_x-x_range,bias_x+x_range,
                                                      bias_y+y_range,bias_y-y_range))
    plt.colorbar()


# plot density and finalize plot    
# rho = rho(midx, T)
# # rho = rho(midx, T)
# ax.plot(midx, rho,  linestyle= "dashdot", c="k", alpha=0.7, label=r"$\rho_{{\beta}}$")
# ax.set_xlabel(r"$\theta$")
# ax.set_ylabel(r"Occurence Frequency")
# ax.set_xlim(plotrange)
# fig.suptitle(title, y=0.95)
# ax.legend()
#%%
# plot distribution error
# fig2, ax2 = plt.subplots()
# # h_axis = [float(h) for h in hs]
# N_axis = [float(N[2:]) for N in Ns]
# c="r"
# ax2.scatter(np.log10(N_axis), np.log10(dist_mse), s=80, label="OBABO", c=c)
# ax2.plot(np.log10(N_axis), np.log10(dist_mse), c=c)
# # ax2.scatter(N_axis, np.log(dist_mse), s=80, label="distributional mse")
# ax2.set_ylabel(r"log$_{{10}}$(mse)")
# ax2.set_xlabel(r"log$_{{10}}$(N)")
# fig2.suptitle(r"Method Convergence on HO, Additive Noise $\sigma$=0.1, $h$=0.05") 
# fig2.legend()
   

plt.show()