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


name = "/home/rene/PhD/Research/Integrators/GM/few_data/Tconf/data/Tconf_"   # head of file name to be read
title = r"GM Sampling, 500 Data Points, Grad. Noise"   # plot titles
title2 = r"GM Sampling, 500 Data Points, Grad. Noise"


# simulation parameters specifying file names to be read
h = "0.010"
Tconf_star = 0.1
L="10"
B="250"
grad = "_gradnoiseB"+B
Bperc = "50%"
avg = "99"

# ## manual entered file names and plot labels
# files = ["OBABO_h0.075_avg99", "OBABO_h0.010_avg99", "OBABO_h0.001_avg99", "OBABO_h0.0005_avg99"]
# labels = [r"OBABO, h=0.075", "OBABO, h=0.01", "OBABO, h=0.001", "OBABO, h=0.0005"]

# files = ["OBABO_h0.075_avg99", "OBABO_h0.010_avg99", "OBABO_h0.0050_avg99",  "OBABO_h0.0010_avg99"]
# labels = [r"OBABO, h=0.075", "OBABO, h=0.01", "OBABO, h=0.005", "OBABO, h=0.001"]


## file names based on simulation parameters above
# 2 OBABO files, one with full gradient, one with partial gradient
files = ["OBABO_h"+h, "OBABO_h"+h, "MOBABO_SF0_L"+L+"_h"+h, "MOBABO_SFR_L"+L+"_h"+h, "OMBABO_SF0_L"+L+"_h"+h, "OMBABO_SFR_L"+L+"_h"+h]  
files = [i +grad+"_avg"+avg for i in files]
files[0]="OBABO_h"+h+"_avg33" # manually enter name of OBABO full gradient
labels = [r"OBABO, $h$="+h + ", Full Gradient", "OBABO, $h$="+h+", $B=$"+Bperc, r"MOBABO SF0, $L$="+L+", $B=$"+Bperc, 
          "MOBABO SFR, $L$="+L+", $B=$"+Bperc, r"OMBABO SF0, $L$="+L+", $B=$"+Bperc, r"OMBABO SFR, $L$="+L+", $B=$"+Bperc]
colors = ["r", "k", "g", "m", "orange", "yellow"] 

# # 1 OBABO file (only with partial gradient)
# files = ["OBABO_h"+h, "MOBABO_SF0_L"+L+"_h"+h, "MOBABO_SFR_L"+L+"_h"+h, "OMBABO_SF0_L"+L+"_h"+h, "OMBABO_SFR_L"+L+"_h"+h] 
# files = [i +grad+"_avg"+avg for i in files]
# labels = [r"OBABO, $h$="+h, r"MOBABO SF0, $L$="+L, r"MOBABO SFA, $L$="+L, "MOBABO SFR, $L$="+L, 
#           r"OMBABO SF0, $L$="+L, r"OMBABO SFA, $L$="+L, r"OMBABO SFR, $L$="+L]
# colors = ["r", "g", "m", "orange", "yellow"]

#%% plot files

fig, ax = plt.subplots()  # for Tconf
fig2, ax2 = plt.subplots()  # for P_accept

for (file, label, c) in zip(files,labels, colors):
    # if c=="k":
    #     continue
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
        ax.plot(n_axis, np.abs(Tconf-Tconf_star), c=c,  label = label) 
        ax2.plot(n_axis, acceptP, c=c, linestyle="--", label = label)


#%% customize plots

ax.set_xlabel(r"$N_{{Samples}}$")
ax.set_ylabel(r"Error")
# ax.set_ylim(-0.05, 0.9)
fig.suptitle(title)
plt.tight_layout()
ax.legend(loc="center right")
ax.set_title(r" $ |-0.5\langle \theta \cdot \nabla U(\theta) \rangle - T_{{conf}}|$,  $T_{{conf}}=$" + str(Tconf_star))
# fig.legend()

ax2.set_xlabel(r"$N_{{Samples}}$")
ax2.set_ylabel(r"$P_{{Acceptance}}$")
fig2.suptitle(title2)
plt.tight_layout()
ax2.legend()
# fig2.legend()
