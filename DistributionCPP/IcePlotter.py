from mpl_toolkits.mplot3d import Axes3D
from pylab import *
import sys,os
import numpy
import re
import pylab
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter, ScalarFormatter, MaxNLocator
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import figure, axes, plot, xlabel, ylabel, title, grid, savefig, show

inches_per_pt = 1.0/72.27
width = 700
height = 300
fig_size = [width*inches_per_pt,height*inches_per_pt]

params = {   #'axes.labelsize': 30,
             #'text.fontsize': 40,
             #'legend.fontsize': 20,
             'xtick.labelsize': 18,
             'ytick.labelsize': 18,             
             'figure.figsize':fig_size,
             #'figure.markersize': 50}
}
pylab.rcParams.update(params)

# ******************************************************
# RUN (COMPUTATIONAL CODE)
# ******************************************************

RUN = 421; AP = "";
NUM = 6;
chord = 0.5334;
runsdir = "/home/adegenna/LagrangianIcingCode/DistributionCPP/RUNS";

for i in range(1,NUM+1):
    XY = genfromtxt(runsdir + "/" + str(RUN) + AP +"/T" + str(i) + "/XY_NEW.out", delimiter = "\t");
    figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);

# ******************************************************
# RUN (LEWICE COMPARISON)
# ******************************************************

lewdir = "/home/adegenna/LagrangianIcingCode/Validations/LewiceIceshapes";
lewrun = genfromtxt(lewdir + "/" + "Run" + str(RUN) + ".csv", delimiter = ",");
plt.scatter(lewrun[:,0]/21.0*chord,lewrun[:,1]/21.0*chord,c='b',s=50);

axis('equal'); plt.xlim([-0.05,0.2])
plt.grid(b=True)

# ******************************************************
# CLEAN AIRFOIL
# ******************************************************

basedir = "/home/adegenna/LagrangianIcingCode/Validations/Ice/Run405Rime/";
XY = genfromtxt("./Grid/NACA0012/NACA0012-SP", delimiter = "\t");
figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'k',linewidth=3);


# UPPER = genfromtxt("./RUNS/409/T0/THERMO_SOLN_UPPER.out", delimiter = "\t");
# LOWER = genfromtxt("./RUNS/409/T0/THERMO_SOLN_LOWER.out", delimiter = "\t");
# BETA  = genfromtxt("./RUNS/409/T0/BETA.out", delimiter = "\t");
# LWC = 0.55e-3; Uinf = 102.8; 
# figure(2);
# subplot(311); plot(UPPER[:,0],UPPER[:,1],'b.-'); plot(LOWER[:,0],LOWER[:,1],'r.-'); #plt.xlim([-0.06,0.03])
# subplot(312); plot(UPPER[:,0],UPPER[:,2],'b.-'); plot(LOWER[:,0],LOWER[:,2],'r.-'); #plt.xlim([-0.06,0.03])
# subplot(313); plot(UPPER[:,0],UPPER[:,3],'b.-'); plot(LOWER[:,0],LOWER[:,3],'r.-'); plot(BETA[:,0],LWC*Uinf*BETA[:,1],'--',c='k'); #plt.xlim([-0.06,0.03])

show()
