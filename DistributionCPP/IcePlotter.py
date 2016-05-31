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
ratio = 1.62
height = 300
width  = ratio*height
fig_size = [width*inches_per_pt,height*inches_per_pt]

params = {   #'axes.labelsize': 30,
             #'text.fontsize': 40,
             #'legend.fontsize': 20,
             'xtick.labelsize': 14,
             'ytick.labelsize': 14,             
             'figure.figsize':fig_size,
             #'figure.markersize': 50}
}
pylab.rcParams.update(params)

# ******************************************************
# RUN (COMPUTATIONAL CODE)
# ******************************************************

RUN = ["206","207","212","213","308","314","316","402","404","405","409","421","422","423","424","425","426","427","428","429"];
AP = "";
NUM = [12,    12,   9,    8,     4,    6,    4,    7,    7,    7,    6,    6,    6,    6,    6,    6,    6,    6,    6,    5  ];
chord = 0.5334;
runsdir = "/home/adegenna/LagrangianIcingCode/DistributionCPP/VALIDATIONS_2";

for k in range(0,np.size(RUN)):
    figure(1);
    K = k+1;
    for i in range(1,NUM[k]+1):
        try:
            XY = genfromtxt(runsdir + "/" + str(RUN[k]) + AP +"/T" + str(i) + "/XY_NEW.out", delimiter = "\t");
            subplot(4,5,K); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3,alpha=0.75);
            plt.xticks([-0.02, 0.10])
            plt.yticks([-0.04, 0.04])
            plt.title(str(RUN[k]),fontweight='bold');
        except:
            print "Could not open" + str(RUN[k]) + AP + "/T" + str(i);
    figure(2); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3); axis('equal'); plt.xlim([-0.02,0.1]);plt.tight_layout()

    # ******************************************************
    # RUN (LEWICE COMPARISON)
    # ******************************************************
    
    figure(1);
    lewdir = "/home/adegenna/LagrangianIcingCode/Validations/LewiceIceshapes";
    lewrun = genfromtxt(lewdir + "/" + "Run" + str(RUN[k]) + ".csv", delimiter = ",");
    plt.scatter(lewrun[:,0]/21.0*chord,lewrun[:,1]/21.0*chord,c='r',s=50);
    axis('equal'); plt.xlim([-0.02,0.1])

    # ******************************************************
    # CLEAN AIRFOIL
    # ******************************************************
    
    basedir = "/home/adegenna/LagrangianIcingCode/Validations/Ice/Run405Rime/";
    XY = genfromtxt("./Grid/NACA0012/NACA0012-SP", delimiter = "\t");
    plot(XY[:,0]*chord,XY[:,1]*chord,'k',linewidth=3);


# UPPER = genfromtxt("./RUNS/409/T0/THERMO_SOLN_UPPER.out", delimiter = "\t");
# LOWER = genfromtxt("./RUNS/409/T0/THERMO_SOLN_LOWER.out", delimiter = "\t");
# BETA  = genfromtxt("./RUNS/409/T0/BETA.out", delimiter = "\t");
# LWC = 0.55e-3; Uinf = 102.8; 
# figure(2);
# subplot(311); plot(UPPER[:,0],UPPER[:,1],'b.-'); plot(LOWER[:,0],LOWER[:,1],'r.-'); #plt.xlim([-0.06,0.03])
# subplot(312); plot(UPPER[:,0],UPPER[:,2],'b.-'); plot(LOWER[:,0],LOWER[:,2],'r.-'); #plt.xlim([-0.06,0.03])
# subplot(313); plot(UPPER[:,0],UPPER[:,3],'b.-'); plot(LOWER[:,0],LOWER[:,3],'r.-'); plot(BETA[:,0],LWC*Uinf*BETA[:,1],'--',c='k'); #plt.xlim([-0.06,0.03])

show()
