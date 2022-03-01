# This script reads and plots configurational temperature data from csv file
# (single column of floats).

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import LogLocator
import csv
from cycler import cycler
from matplotlib import cm

plt.close('all')
plt.rcParams.update({'font.size': 33})
plt.rc('legend', fontsize=33)    # legend fontsize
plt.rcParams['axes.grid'] = True
plt.rc('lines', linewidth=3)
plt.rcParams['axes.titlepad'] = 0

#%%


# name = "/home/rene/PhD/Research/code/Integrators/GaussianMixtureSampling/Tconf_OBABO_"   # file name head to read
# title = r"GM Sampling, Absolute Error in $T_{{conf}}$, Gradient Noise"
name = "/home/rene/PhD/Research/code/Integrators/GaussianMixtureSampling/Traj"   # file name head to read
title = r"GM Sampling, Absolute Error in $T_{{conf}}$"

# files = ["h0.050", "MH_SF0_L1_h0.050", "MH_SF0_L10_h0.050", "MH_SF0_L50_h0.050", 
#          "MH_SF1_L1_h0.050", "MH_SF1_L10_h0.050", "MH_SF1_L50_h0.050", "h0.010"]
# labels = [r"OBABO, $h$=0.05",  r"MH-only, $h$=0.05, $L$=1", r"MH-only, $h$=0.05, $L$=10", r"MH-only, $h$=0.05, $L$=50", 
#           r"SF, $h$=0.05, $L$=1", r"SF, $h$=0.05, $L$=10", r"SF, $h$=0.05, $L$=50", r"OBABO, $h$=0.01",]
# files = ["h0.010","h0.010_gradnoiseB495", "h0.010_gradnoiseB490", "h0.010_gradnoiseB450"]
# labels = [r"OBABO, $h$=0.01", r"OBABO, $h$=0.01, $B$=99%",r"OBABO, $h$=0.01, $B$=98%", r"OBABO, $h$=0.01, $B$=90%"]
# files = ["h0.010", "h0.010_gradnoiseB490", "h0.050", "h0.050_gradnoiseB490"]
# labels = [r"OBABO, $h$=0.01", r"OBABO, $h$=0.01, $B$=98%", r"OBABO, $h$=0.05", r"OBABO, $h$=0.05, $B$=98%"]
files = ["1 Tconf_OBABO_h0.050", "2 Tconf_OBABO_h0.050", "3 Tconf_OBABO_h0.050"]
labels = [r"OBABO, $h$=0.05, seed 1", r"OBABO, $h$=0.05, seed 2", r"OBABO, $h$=0.05, seed 3"]

# colors = ["b", "orange",  "g", "c", "m","pink" , "k","r"]
colors = [ "c", "m", "b"]

dist_mse = []
fig, ax = plt.subplots()

for (file, label, c) in zip(files,labels, colors):
    print(name+file)
    with open(name+file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=' ')
        Tconf = []
        n_axis = []
        for row in csv_reader:
            n_axis += [int(row[0])]
            Tconf += [float(row[1])]
        
        # ax.plot(np.arange(0,len(Tconf)), Tconf)
        Tconf = np.array(Tconf)
        # Tconf_avg = [np.mean(Tconf[0:i]) for i in np.arange(1, len(Tconf), n_avg)]
        # ax.plot(np.arange(1, len(Tconf), n_avg), Tconf_avg, c=c, label=label)
        ax.plot(n_axis[1::], np.abs(Tconf[1::]-2), c=c, label = label)



#%%

ax.set_xlabel("No. of Samples")
ax.set_ylabel(r"Error")
fig.suptitle(title)
plt.legend()
plt.show()