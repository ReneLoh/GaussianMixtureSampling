# This script reads and plots configurational temperature data from csv file
# (single column of floats).

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


# name = "/home/rene/PhD/Research/Integrators/GM/few_data/data/Tconf_"   # file name head to read
name = "/home/rene/PhD/Research/Integrators/GM/few_data/data/mus_"   # file name head to read
# title = r"GM Sampling, Absolute Error in $T_{{conf}}$, Correction of Disc. Bias"
# title2 = r"GM Sampling, Acceptance Probability, Correction of Disc. Bias"
# title = r"GM Sampling, Absolute Error in $T_{{conf}}$, Correction of Grad. Noise"
# title = r"GM Sampling, $|\theta \cdot \nabla U(\theta) -N_{{d}}T|$"
# title = r"GM Sampling, t-aveaged $|\mu_{{i}}-\mu_{{i}}^{{*}}|$, Disc. Bias."
title = r"GM Sampling, 500 Data Points, Grad. Noise"
title2 = r"GM Sampling, 500 Data Points, Grad. Noise"


mus_true = np.array([-2.5, 2.5])
post_argmax = np.array([-2.7817817817817816, 2.6536536536536537 ])  # argmax of posterior
h = "0.010"
avg = "99"
Tconf_star = 1
L="10"

# files = ["OBABO_h"+h, "OBABO_h"+h, "MOBABO_SF0_L"+L+"_h"+h, "MOBABO_SFA_L"+L+"_h"+h, 
#           "MOBABO_SFR_L"+L+"_h"+h, "OMBABO_SF0_L"+L+"_h"+h, "OMBABO_SFA_L"+L+"_h"+h, "OMBABO_SFR_L"+L+"_h"+h]
# grad = "_gradnoiseB375"
# files = [i +grad+"_avg"+avg for i in files]
# files[0]="OBABO_h"+h+"_avg99"
# # labels = [r"OBABO, $h$="+h, r"MOBABO SF0, $L$="+L, r"MOBABO SFA, $L$="+L, "MOBABO SFR, $L$="+L, 
# #           r"OMBABO SF0, $L$="+L, r"OMBABO SFA, $L$="+L, r"OMBABO SFR, $L$="+L]
# labels = [r"OBABO, $h$="+h + ", Full Gradient", "OBABO, $h$="+h+", $B=$75%", r"MOBABO SF0, $L$="+L+", $B=$75%", 
#           r"MOBABO SFA, $L$="+L+", $B=$75%", "MOBABO SFR, $L$="+L+", $B=$75%", 
#           r"OMBABO SF0, $L$="+L+", $B=$75%", r"OMBABO SFA, $L$="+L+", $B=$75%", r"OMBABO SFR, $L$="+L+", $B=$75%"]

### WITHOUT SFA
files = ["OBABO_h"+h, "OBABO_h"+h, "MOBABO_SF0_L"+L+"_h"+h, "MOBABO_SFR_L"+L+"_h"+h, "OMBABO_SF0_L"+L+"_h"+h, "OMBABO_SFR_L"+L+"_h"+h]
grad = "_gradnoiseB250"
files = [i +grad+"_avg"+avg for i in files]
files[0]="OBABO_h"+h+"_avg99"
# labels = [r"OBABO, $h$="+h, r"MOBABO SF0, $L$="+L, r"MOBABO SFA, $L$="+L, "MOBABO SFR, $L$="+L, 
#           r"OMBABO SF0, $L$="+L, r"OMBABO SFA, $L$="+L, r"OMBABO SFR, $L$="+L]
labels = [r"OBABO, $h$="+h + ", Full Gradient", "OBABO, $h$="+h+", $B=$50%", r"MOBABO SF0, $L$="+L+", $B=$50%", 
          "MOBABO SFR, $L$="+L+", $B=$50%", r"OMBABO SF0, $L$="+L+", $B=$50%", r"OMBABO SFR, $L$="+L+", $B=$50%"]
###

# files = ["small_OBABO_h"+h, "small_OBABO_h"+h+"_gradnoiseB490", "large_OBABO_h"+h, "large_OBABO_h"+h+"_gradnoiseB4900"]
# labels = [r"500 Data Points, Full Gradient", "500 Data Points, B=98%", 
#           "5k Data Points, Full Gradient", "5k Data Points, B=98%"  ]
# files = ["OBABO_h"+h, "OBABO_h"+h+"_gradnoiseB490", "OBABO_h"+h+"_gradnoiseB450",
#          "OBABO_h"+h+"_gradnoiseB250","OBABO_h"+h+"_gradnoiseB100"]
# labels = [r"B=100%", "B=98%", "B=90%", "B=50%", "B=20%"]

# hs = ["0.001", "0.010", "0.100"]
# files = ["OBABO_h"+h + "_avg99" for h in hs]
# labels = [r"OBABO, $h=0.001$", r"OBABO, $h=0.01$", r"OBABO, $h=0.1$"]

# files = ["MOBABO_SF0_L"+L+"_h"+h+"_newINIT",]
# labels = ["MOBABO SF0, $L$="+L]

# colors = ["b", "orange",  "g", "c", "m","pink" , "k","r"]
# colors = ["r", "g", "c", "m", "orange", "b", "yellow"]
colors = ["r", "k", "g", "c", "m", "orange", "b", "yellow"]  # in case of two OBABO (second color)
colors = ["r", "k", "g", "m", "orange", "yellow"]  # in case of two OBABO (second color), w/o SFA

# fig, ax = plt.subplots()  # for Tconf
fig, ax = plt.subplots(1,2)    # for mu_i
fig2, ax2 = plt.subplots()

# for (file, label, c) in zip(files,labels, colors):
for (file, label, c) in zip(files,labels, colors):
    # if c=="r":
    #     continue
    print(name+file)
    with open(name+file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=' ')
        Tconf = []
        mu1 = []
        mu2 = []
        acceptP = []
        n_axis = []
        for row in csv_reader:
            n_axis += [int(row[0])]
            # Tconf += [float(row[1])]
            mu1 += [float(row[1])]
            mu2 += [float(row[2])]
            acceptP += [float(row[3])]
        
        # Tconf = np.array(Tconf)
        # TconfAVG = np.array( [np.mean(Tconf[i-30000:i]) for i in range(300000, len(Tconf))]  )
        mu1 = np.array(mu1)
        # mu1AVG = np.array( [np.mean(mu1[i-n:i]) for i in range(n, len(mu1), ndist)]  )
        mu2 = np.array(mu2)
        # mu2AVG = np.array( [np.mean(mu2[i-n:i]) for i in range(n, len(mu2), ndist)]  )
        acceptP = np.array(acceptP)
   
        # ax.plot(n_axis, Tconf, c=c, label = label)                 
        # ax[0].plot(n_axis[n::ndist], np.abs(mu1AVG-mus_true[0]), c=c,  label = label) 
        # ax[1].plot(n_axis[n::ndist], np.abs(mu2AVG-mus_true[1]), c=c,  label = label) 
        # ax[0].plot(n_axis[n::ndist], np.abs(mu1AVG-post_argmax[0]), c=c,  label = label) 
        # ax[1].plot(n_axis[n::ndist], np.abs(mu2AVG-post_argmax[1]), c=c,  label = label) 
        ax[0].plot(n_axis, np.abs(mu1-post_argmax[0]), c=c,  label = label) 
        ax[1].plot(n_axis, np.abs(mu2-post_argmax[1]), c=c,  label = label) 
        # ax.plot(n_axis, np.abs(Tconf-Tconf_star), c=c,  label = label)
        # ax.plot(n_axis[300000::], np.abs(TconfAVG-2), c=c,  label = label)
        ax2.plot(n_axis, acceptP, c=c, linestyle="--", label = label)


#%%

# ax.set_xlabel(r"$N_{{Samples}}$")
ax[0].set_xlabel(r"$N_{{Samples}}$")
ax[0].set_ylabel(r"Error")
ax[1].set_xlabel(r"$N_{{Samples}}$")
fig.suptitle(title)
ax[1].legend(loc="upper right", bbox_to_anchor=(1.3,1))
ax[0].set_ylim(0, 0.01)
ax[1].set_ylim(0, 0.0035)
ax[0].set_title(r" $\langle |\mu_{{1}}-\mathrm{argmax}_{\mu_{1}}\mathrm{P}(\theta|X)|\rangle$")
ax[1].set_title(r" $\langle |\mu_{{2}}-\mathrm{argmax}_{\mu_{2}}\mathrm{P}(\theta|X)|\rangle$")
# fig.legend()

ax2.set_xlabel(r"$N_{{Samples}}$")
ax2.set_ylabel(r"$P_{{Acceptance}}$")
fig2.suptitle(title2)
ax2.legend()
# fig2.legend()
