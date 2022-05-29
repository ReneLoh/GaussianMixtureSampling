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


name = "/home/rene/PhD/Research/Integrators/GM/much_data/efficiency/Tconf_"   # head of file name to be read
# name = "/home/rene/PhD/Research/Code/Integrators/GM/Tconf/Tconf_"   # head of file name to be read

title = r"GM Sampling, 5000 Data Points, Grad. Noise"   # plot titles
title2 = r"GM Sampling, 5000 Data Points, Grad. Noise"


# simulation parameters specifying file names to be read
h = "0.010"
Tconf_star = 5
L="1"
B="25"
grad = "_gradnoiseB"+B
Bperc = "50%"
avg = "66"

# ## manual entered file names and plot labels
# files = ["OBABO_h0.010",  "OBABO_h0.010_gradnoiseB250",  "OBABO_h0.010_gradnoiseB25",  "OBABO_h0.010_gradnoiseB5", 
#          "OMBABO_SFR_L1_h0.010_gradnoiseB250",  "OMBABO_SFR_L1_h0.010_gradnoiseB25",  "OMBABO_SFR_L1_h0.010_gradnoiseB5"]

# labels = [r"OBABO Full Gradient",  "OBABO B=50%",  "OBABO B=5%",  "OBABO B=1%", 
#           "OMBABO SFR, L=1, B=50%",  "OMBABO SFR, L=1, B=5%",  "OMBABO SFR, L=1, B=1%"]


files = ["OBABO_h0.010",  "OBABO_h0.010_gradnoiseB250",  "OBABO_h0.010_gradnoiseB50",  "OBABO_h0.010_gradnoiseB1", 
         "OMBABO_SFR_L1_h0.010_gradnoiseB250",  "OMBABO_SFR_L1_h0.010_gradnoiseB50",  "OMBABO_SFR_L1_h0.010_gradnoiseB1"]

labels = [r"OBABO Full Gradient",  "OBABO B=5%",  "OBABO B=1%",  "OBABO B=1pt", 
          "OMBABO SFR, L=1, B=5%",  "OMBABO SFR, L=1, B=1%",  "OMBABO SFR, L=1, B=1pt"]

# files = ["OBABO_h0.010",  "OBABO_h0.010_gradnoiseB250",  "OBABO_h0.010_gradnoiseB50",  "OBABO_h0.010_gradnoiseB1", 
#          "OMBABO_SFR_L10_h0.010_gradnoiseB250",  "OMBABO_SFR_L10_h0.010_gradnoiseB50",  "OMBABO_SFR_L10_h0.010_gradnoiseB1"]

# labels = [r"OBABO Full Gradient",  "OBABO B=5%",  "OBABO B=1%",  "OBABO B=1pt", 
#           "OMBABO SFR, L=10, B=5%",  "OMBABO SFR, L=10, B=1%",  "OMBABO SFR, L=10, B=1pt"]

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
