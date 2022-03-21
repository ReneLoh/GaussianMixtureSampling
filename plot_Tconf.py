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


name = "/home/rene/PhD/Research/Integrators/GM/few_data/data/Tconf_SRC_"   # file name head to read
title = r"GM Sampling, Absolute Error in $T_{{conf}}$, Correction of Disc. Bias"
title2 = r"GM Sampling, Acceptance Probability, Correction of Disc. Bias"
# name = "/home/rene/PhD/Research/code/Integrators/GaussianMixtureSampling/Traj"   # file name head to read
# title = r"GM Sampling, Absolute Error in $T_{{conf}}$"
h = "0.080"
L="1"
files = ["OBABO_h"+h, "MOBABO_SF0_L"+L+"_h"+h, "MOBABO_SFA_L"+L+"_h"+h, 
          "MOBABO_SFR_L"+L+"_h"+h, "OMBABO_SF0_L"+L+"_h"+h, "OMBABO_SFA_L"+L+"_h"+h, "OMBABO_SFR_L"+L+"_h"+h]
labels = [r"OBABO, $h$="+h,  r"MOBABO SF0, $L$="+L, r"MOBABO SFA, $L$="+L, "MOBABO SFR, $L$="+L, 
          r"OMBABO SF0, $L$="+L, r"OMBABO SFA, $L$="+L, r"OMBABO SFR, $L$="+L]

# files = ["MOBABO_SF0_L"+L+"_h"+h+"_newINIT",]
# labels = ["MOBABO SF0, $L$="+L]

# colors = ["b", "orange",  "g", "c", "m","pink" , "k","r"]
colors = ["r", "g", "c", "m", "orange", "b", "yellow"]

dist_mse = []
fig, ax = plt.subplots()
fig2, ax2 = plt.subplots()

for (file, label, c) in zip(files,labels, colors):
    print(name+file)
    with open(name+file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=' ')
        Tconf = []
        acceptP = []
        n_axis = []
        for row in csv_reader:
            n_axis += [int(row[0])]
            Tconf += [float(row[1])]
            acceptP += [float(row[2])]
        
        Tconf = np.array(Tconf)
        acceptP = np.array(acceptP)
  
        ax.plot(n_axis, np.abs(Tconf-2), c=c,  label = label)
        ax2.plot(n_axis, acceptP, c=c, linestyle="--", label = label)


#%%

ax.set_xlabel(r"$N_{{Samples}}$")
ax.set_ylabel(r"Error")
fig.suptitle(title)
ax.legend()
# fig.legend()

ax2.set_xlabel(r"$N_{{Samples}}$")
ax2.set_ylabel(r"$P_{{Acceptance}}$")
fig2.suptitle(title2)
ax2.legend()
# fig2.legend()
