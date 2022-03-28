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



# %% MAIN
        
seed = 2
np.random.seed(seed)
seeds = [1]

## create data to true density
mus_true = [-2.5,2.5]
sigs = [3,0.5]
a_s = [0.8, 0.2]
N = 500
X = GM.create_data(mus_true, sigs, a_s, N)
# GM.write_to_file(X, "GM_data_5000.csv")

# xrange = (-12.5,12.5)
# plt.hist(X,100,xrange, density=True, label="Data Distribution") 
# x = np.linspace(xrange[0], xrange[1],1000)
# y = GM.gaussmix(*mus_true, x, sigs,a_s)
# plt.plot(x,y, label="Underlying Density")

#prior specifics
prior_func = GM.prior_gauss
(m1,m2,sig0) = (0,0,5)
prior_params = (m1,m2,sig0)

# # evaluate the posterior on grid
# mu1, mu2 = (np.linspace(-7.5,7.5,1000), np.linspace(-7.5,7.5,1000))
mu1, mu2 = (np.linspace(-5,-1,1000), np.linspace(1,5,1000))
POST = GM.posterior_plot(mu1[:,None], mu2[None,:], X, GM.gaussmix, (sigs, a_s), prior_func, prior_params)
plt.imshow(POST, cmap=cm.jet, extent=((mu2[0], mu2[-1], mu1[-1], mu1[0])))
plt.clim(POST.max()-20,POST.max())
plt.colorbar()
plt.xlabel(r"$\mu_{2}$")
plt.ylabel(r"$\mu_{1}$")
mu1_star = mu1[ np.where( abs(POST-POST.max())<0.0001 )[0][0] ]
mu2_star = mu2[ np.where( abs(POST-POST.max())<0.0001 )[1][0] ]


#%% simulation params
theta0 = np.array([2.,3.])
p0 = np.zeros(2, float)
N_0 = 1e3
# Ns = [100, 500, 100000]
h = 0.05
# hs = [0.1]
# hs = [ 1, 0.5, 0.1]
# T = [1e-3, 1e-2, 1e-1]
T = 1
gamma = 1

L=20
SF=1
# label_MH_SF = r"MH-only, $L$={} $h=${}".format(L,h)
label_MH_SF = r"SF, $L$={}".format(L)

xrange = (-5, 10)
title = r"GM Posterior Sampling"


# %%
fig, ax = plt.subplots(1,2)

force_params = (T, GM.gaussmix, sigs, a_s, X)
# for N_0 in Ns:
# for h in hs:
for seed in seeds:
    np.random.seed(seed)
    TAVG_OBABO = GM.OBABO_simu(theta0, p0, int(N_0), h, T, gamma, GM.noisy_force, force_params)
    # TAVG_MH = OBABO_simu_MH(theta0, p0, int(N_0/(L+1)), h, T, gamma, L, SF, noisy_force, force_params, posterior,
                        # X, gaussmix, (sigs, a_s))
    
    # histo_OBABO, bins = np.histogram(theta_OBABO[:,0], bins=200, range=xrange, density=True)
    # midx = (bins[0:-1]+bins[1:])/2
    # ax[0].plot(midx, histo_OBABO, label=r"$N$={}".format(N_0))
    # histo_OBABO, bins = np.histogram(theta_OBABO[:,1], bins=200, range=xrange, density=True)
    # midx = (bins[0:-1]+bins[1:])/2
    # ax[1].plot(midx, histo_OBABO, label=r" $N$={}".format(N_0))

    ax[0].plot(np.arange(0,N_0+1), np.abs(TAVG_OBABO[:,0]-mu1_star), 
                        label=r"OBABO, $h=${}".format(h))
    ax[1].plot(np.arange(0,N_0+1), np.abs(TAVG_OBABO[:,1]-mu2_star), 
                        label=r"OBABO, $h=${}".format(h))
    
    # ax[0].plot(np.arange(0,N_0+1, (L+1)), np.abs(TAVG_MH[:,0]-mu1_star), 
    #                     label=label_MH_SF+r", $h$={}".format(h))
    # ax[1].plot(np.arange(0,N_0+1, (L+1)), np.abs(TAVG_MH[:,1]-mu2_star), 
    #                     label=label_MH_SF+r", $h$={}".format(h))


    
#%%
# fig, ax = plt.subplots()
# histo_OBABO, bins = np.histogram(theta_OBABO[:,0], bins=200, range=xrange, density=True)
# midx = (bins[0:-1]+bins[1:])/2
# ax.plot(midx, histo_OBABO, label=r"mu1")
# histo_OBABO, bins = np.histogram(theta_OBABO[:,1], bins=200, range=xrange, density=True)
# midx = (bins[0:-1]+bins[1:])/2
# ax.plot(midx, histo_OBABO, label=r"mu2")

# plt.title(title)
# plt.ylabel("|$\hat{{\mu}}_{{1}}-\mu_{{1}}^{{*}}$|")
# plt.xlabel(r"$N$")
# plt.legend()



fig.suptitle(title)
ax[0].set_ylabel("|$\hat{{\mu}}_{{1}}-\mu_{{1}}^{{*}}$|")
ax[1].set_ylabel("|$\hat{{\mu}}_{{2}}-\mu_{{2}}^{{*}}$|")
ax[0].set_xlabel(r"$N$")
ax[1].set_xlabel(r"$N$")
ax[0].set_yscale("log")
ax[0].set_xscale("log")
ax[1].set_yscale("log")
ax[1].set_xscale("log")
plt.legend()
