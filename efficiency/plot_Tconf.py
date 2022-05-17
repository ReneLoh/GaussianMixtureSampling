# This script reads and plots sampling data from csv file.
# 3 columns. 1) Iteration, 2) Configurational Temperature, 3) Metropolis Acceptance Probability

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


name = "/home/rene/PhD/Research/Integrators/GM/few_data/efficiency/Tconf_"   # head of file name to be read
# name = "/home/rene/PhD/Research/Code/Integrators/GM/Tconf/Tconf_"   # head of file name to be read

title = r"GM Sampling, 500 Data Points, Grad. Noise"   # plot titles
title2 = r"GM Sampling, 500 Data Points, Grad. Noise"


# simulation parameters specifying file names to be read
h = "0.010"
Tconf_star = 5
L="1"
B="25"
grad = "_gradnoiseB"+B
Bperc = "50%"
avg = "66"

# ## manual entered file names and plot labels
files = ["OBABO_h0.010", "OBABO_h0.010_gradnoiseB250", "OBABO_h0.010_gradnoiseB25","OBABO_h0.010_gradnoiseB5", 
         "OMBABO_SFR_L1_h0.010_gradnoiseB250","OMBABO_SFR_L1_h0.010_gradnoiseB25", "OMBABO_SFR_L1_h0.010_gradnoiseB5"]
labels = [r"OBABO Full Gradient", "OBABO B=50%", "OBABO B=5%", "OBABO B=1%", "OMBABO SFR, L=1, B=50%", "OMBABO SFR, L=1, B=5%", "OMBABO SFR, L=1, B=1%"]

# files = ["OBABO_h0.075_avg99", "OBABO_h0.010_avg99", "OBABO_h0.0050_avg99",  "OBABO_h0.0010_avg99"]
# labels = [r"OBABO, h=0.075", "OBABO, h=0.01", "OBABO, h=0.005", "OBABO, h=0.001"]


## file names based on simulation parameters above
#2 OBABO files, one with full gradient, one with partial gradient
# files = ["OBABO_h"+h, "OBABO_h"+h, "OMBABO_SFR_L"+L+"_h"+h]  
# files = [i +grad+"_avg"+avg for i in files]
# files[0]="OBABO_h"+h+"_avg66" # manually enter name of OBABO full gradient
# labels = [r"OBABO, $h$="+h + ", Full Gradient", "OBABO, $h$="+h+", $B=$"+Bperc, r"OMBABO SFR, $L$="+L+", $B=$"+Bperc]
# colors = ["r", "k", "g", "m", "orange", "yellow"] 

# # 1 OBABO file (only with partial gradient)
# files = ["OBABO_h"+h, "MOBABO_SF0_L"+L+"_h"+h, "MOBABO_SFR_L"+L+"_h"+h, "OMBABO_SF0_L"+L+"_h"+h, "OMBABO_SFR_L"+L+"_h"+h] 
# files = [i +grad+"_avg"+avg for i in files]
# labels = [r"OBABO, $h$="+h, r"MOBABO SF0, $L$="+L, r"MOBABO SFA, $L$="+L, "MOBABO SFR, $L$="+L, 
#           r"OMBABO SF0, $L$="+L, r"OMBABO SFA, $L$="+L, r"OMBABO SFR, $L$="+L]


colors = ["r", "g", "m", "b", "orange", "cyan", "k"]

#%% plot files

fig, ax = plt.subplots()  # for Tconf
fig2, ax2 = plt.subplots()  # for P_accept

for (file, label, c) in zip(files,labels, colors):
    # if c=="k":
    #     continue
    print(name+file)
    with open(name+file+"_avg"+avg) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=' ')
        Tconf = []
        acceptP = []
        n_axis = []
        for row in csv_reader:
            n_axis += [float(row[0])]
            Tconf += [float(row[1])]
            acceptP += [float(row[2])]

        Tconf = np.array(Tconf)
        acceptP = np.array(acceptP)            
        ax.plot(n_axis, np.abs(Tconf-Tconf_star), c=c,  label = label) 
        if not np.all(acceptP == 1):
            ax2.plot(n_axis, acceptP, c=c, linestyle="--", label = label)


#%% customize plots

ax.set_xlabel(r"Time [$\Delta t_{{\mathrm{Full \: Gradient}}}$]")
ax.set_ylabel(r"Error")
# ax.set_ylim(-0.05, 0.9)
fig.suptitle(title)
plt.tight_layout()
ax.legend(loc="center right")
ax.set_title(r" $ |-0.5\langle \theta \cdot \nabla U(\theta) \rangle - T_{{conf}}|$")
# fig.legend()

ax2.set_xlabel(r"Time [$\Delta t_{{\mathrm{Full \: Gradient}}}$]")
ax2.set_ylabel(r"$P_{{\mathrm{Acceptance}}}$")
fig2.suptitle(title2)
plt.tight_layout()
ax2.legend()
# fig2.legend()
