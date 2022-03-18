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


name = "/home/rene/PhD/Research/Integrators/GM/few_data/data/Tconf_"   # file name head to read
title = r"GM Sampling, Absolute Error in $T_{{conf}}$, Correction of Disc. Bias"
# name = "/home/rene/PhD/Research/code/Integrators/GaussianMixtureSampling/Traj"   # file name head to read
# title = r"GM Sampling, Absolute Error in $T_{{conf}}$"
h = "0.030"
L="1"
files = ["OBABO_h"+h, "MOBABO_SF0_L"+L+"_h"+h, "MOBABO_SFA_L"+L+"_h"+h, 
          "MOBABO_SFR_L"+L+"_h"+h, "OMBABO_SF0_L"+L+"_h"+h, "OMBABO_SFA_L"+L+"_h"+h, "OMBABO_SFR_L"+L+"_h"+h]
labels = [r"OBABO, $h$="+h,  r"MOBABO SF0, $L$="+L, r"MOBABO SFA, $L$="+L, "MOBABO SFR, $L$="+L, 
          r"OMBABO SF0, $L$="+L, r"OMBABO SFA, $L$="+L, r"OMBABO SFR, $L$="+L]


# files = ["OMBABO_L10_h0.050"]
# labels = ["test"]

# colors = ["b", "orange",  "g", "c", "m","pink" , "k","r"]
colors = ["r", "g", "c", "m", "orange", "b", "yellow"]

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
        # ax.plot(n_axis[1::], np.abs(Tconf[1::]-2), c=c, label = label)
        ax.plot(n_axis, np.abs(Tconf-2), c=c, label = label)



#%%

ax.set_xlabel(r"$N_{{Samples}}$")
ax.set_ylabel(r"Error")
fig.suptitle(title)
plt.legend()
plt.show()