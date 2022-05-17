# This script generates a data set of a 1D Gaussian Mixture (GM) model, consisting
# of the sum of two independent Gaussians to specified parameters. 
# Furthermore, it computes the log-posterior of the means of the two Gaussians given the data and a prior.
# It can plot the true density, the data histogram, and a heatmap of the log-posterior.
# It relies on the functions in the script GM_functions.py.


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pylab import figure, cm
from matplotlib.ticker import LogLocator
import csv
from cycler import cycler
import time
import GM_functions as GM

# plt.close('all')
plt.rcParams.update({'font.size': 44})
plt.rc('legend', fontsize=40)    # legend fontsize
plt.rcParams['axes.grid'] = True
plt.rc('lines', linewidth=5)
plt.rcParams['axes.titlepad'] = 26


        
seed = 2
np.random.seed(seed)
seeds = [1]

## create data to true density
mus_true = [-2.5,2.5]
sigs = [3,0.5]
a_s = [0.8, 0.2]
N = 5000
X = GM.create_data(mus_true, sigs, a_s, N)
GM.write_to_file(X, "GM_data_5000.csv")

## plot density and data histogram
# xrange = (-12.5,12.5)
# plt.hist(X,100,xrange, density=True, label="Data Distribution") 
# x = np.linspace(xrange[0], xrange[1],1000)
# y = GM.gaussmix(*mus_true, x, sigs,a_s)
# plt.plot(x,y, label="Underlying Density")

# prior specifics
prior_func = GM.prior_gauss
(m1,m2,sig0) = (0,0,5)
prior_params = (m1,m2,sig0)

## evaluate the posterior on grid
mu1, mu2 = (np.linspace(-5,-1,1000), np.linspace(1,5,1000))  # grid
POST = GM.posterior_plot(mu1[:,None], mu2[None,:], X, GM.gaussmix, (sigs, a_s), prior_func, prior_params)  # compute log-posterior

plt.imshow(POST, cmap=cm.jet, extent=((mu2[0], mu2[-1], mu1[-1], mu1[0])))  # plot result
plt.clim(POST.max()-20,POST.max())
plt.colorbar()
plt.xlabel(r"$\mu_{2}$")
plt.ylabel(r"$\mu_{1}$")

mu1_star = mu1[ np.where( abs(POST-POST.max())<0.0001 )[0][0] ]  # find modes of log-posterior
mu2_star = mu2[ np.where( abs(POST-POST.max())<0.0001 )[1][0] ]
