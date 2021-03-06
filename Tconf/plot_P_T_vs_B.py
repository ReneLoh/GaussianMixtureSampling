# This script reads and plots sampling data from csv file.
# 3 columns. 1) Iteration, 2) Configurational Temperature, 3) Metropolis Acceptance Probability.
# Specifically, it plots the Metropolis acceptance probability and the error in the configurational
# temperature against subsampling gradient noise (i.e. against batch size B), both averaged over the
# last n rows in the files.

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import LogLocator
import csv
from cycler import cycler
from matplotlib import cm
import matplotlib as mp

plt.close('all')
plt.rcParams.update({'font.size': 35})
plt.rc('legend', fontsize=35)    # legend fontsize
plt.rcParams['axes.grid'] = True
plt.rc('lines', linewidth=3)
plt.rcParams['axes.titlepad'] = 0
mp.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

#%%

name = "/home/rene/PhD/Research/Integrators/GM/much_data/data/Tconf_"    # head of file name to be read
title = r"GM Sampling, 5000 Data Points, Grad. Noise"    # plot titles

# parameters specifying file names to be read and how to process data
h = "0.005"
Xsize = 5000
Tconf_star = 1
L="10"
Bs = ["1", "50", "500", "5000"]  # batch sizes for which files are present to be plotted
grad = "_gradnoiseB"
avg = "99"
n = 200

files = ["OBABO_h"+h, "MOBABO_SF0_L"+L+"_h"+h, "MOBABO_SFR_L"+L+"_h"+h, "OMBABO_SF0_L"+L+"_h"+h, "OMBABO_SFR_L"+L+"_h"+h]
labels = [r"OBABO", r"MOBABO SF0, $L$="+L,"MOBABO SFR, $L$="+L, r"OMBABO SF0, $L$="+L, r"OMBABO SFR, $L$="+L]

colors = ["k", "g", "m", "orange", "yellow"]  

# arrays storing error in Tconf and acceptance probabilities for each files and batch size
T_err = np.empty( (len(files), len(Bs)) )
Probs = np.empty( (len(files), len(Bs)) )

for i in range(0, len(files)):

    for j in range(0, len(Bs)):
        
        # if batch size equals size of data set, the file name structure is different
        if int(Bs[j])==Xsize:   
            file = name + files[i] + "_avg" + avg  # full batch file name
        else: 
            file = name + files[i] + grad + Bs[j] + "_avg" + avg   # grad. noise file name
        
        with open(file) as csv_file:
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
            
            # average over the last n rows
            T_err[i,j] = np.mean(np.abs(Tconf-Tconf_star)[-n::])
            Probs[i,j] = np.mean(acceptP[-n::])




#%% plot results
fig, ax = plt.subplots() 
fig2, ax2 = plt.subplots() 

for i in range(0, len(files)):
    if i==0:
        continue
    ax.plot(Bs, T_err[i,:], "-o", c=colors[i] , label=labels[i])
    ax2.plot(Bs, Probs[i,:], "-o", c=colors[i], label=labels[i])

# fig.suptitle(title)
plt.tight_layout()
ax.set_xlabel(r"$B$")
ax.set_ylabel(r"Error")
# ax.legend(loc="upper right", bbox_to_anchor=(1.3,1))
ax.legend()
ax.set_title(r" Error in Config. Temp. vs. Batch Size")
ax.set_yscale("log")

ax2.set_xlabel(r"$B$")
ax2.set_ylabel(r"$P_{{acc}}$")
# ax2.legend(loc="upper right", bbox_to_anchor=(1.3,1))
ax2.legend()
ax2.set_title(r"Metropolis Acceptance Probability vs. Batch Size")
ax2.set_yscale("log")

