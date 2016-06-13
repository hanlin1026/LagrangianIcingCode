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
ratio = 1.0
height = 700
width  = ratio*height
fig_size = [width*inches_per_pt,height*inches_per_pt]

params = {   #'axes.labelsize': 30,
             #'text.fontsize': 40,
             #'legend.fontsize': 20,
             'xtick.labelsize': 24,
             'ytick.labelsize': 24,             
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
runsdir  = "/home/adegenna/LagrangianIcingCode/DistributionCPP/VALIDATIONS_2";
runsdir2 = "/home/adegenna/LagrangianIcingCode/DistributionCPP/TMP";

for k in range(0,np.size(RUN)):
    figure(1);
    K = k+1;
    for i in range(1,NUM[k]+1):
        try:
            XY = genfromtxt(runsdir + "/" + str(RUN[k]) + AP +"/T" + str(i) + "/XY_NEW.out", delimiter = "\t");
            subplot(5,4,K); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3,alpha=0.75);
            plt.xticks([-0.02, 0.10])
            plt.yticks([-0.04, 0.04])
            plt.title(str(RUN[k]),fontweight='bold',fontsize=30);
        except:
            print "Could not open" + str(RUN[k]) + AP + "/T" + str(i);
            # XY = genfromtxt(runsdir2 + "/" + str(RUN[k]) + AP +"/T" + str(i) + "/XY_NEW.out", delimiter = "\t");
            # subplot(5,4,K); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3,alpha=0.75);
            # plt.xticks([-0.02, 0.10])
            # plt.yticks([-0.04, 0.04])
            # plt.title(str(RUN[k]),fontweight='bold',fontsize=30);
    figure(2); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3); axis('equal'); plt.xlim([-0.02,0.1]); plt.tight_layout()

    # ******************************************************
    # RUN (LEWICE COMPARISON)
    # ******************************************************
    
    figure(1);
    lewdir = "/home/adegenna/LagrangianIcingCode/Validations/LewiceIceshapes";
    if (RUN[k] == "308_2"):
        lewrun = genfromtxt(lewdir + "/" + "Run308.csv", delimiter = ",");
    else:
        lewrun = genfromtxt(lewdir + "/" + "Run" + str(RUN[k]) + ".csv", delimiter = ",");
    plt.scatter(lewrun[:,0]/21.0*chord,lewrun[:,1]/21.0*chord,c='r',s=50);
    axis('equal'); plt.xlim([-0.02,0.1])

    # ******************************************************
    # CLEAN AIRFOIL
    # ******************************************************
    
    basedir = "/home/adegenna/LagrangianIcingCode/Validations/Ice/Run405Rime/";
    XY = genfromtxt("./Grid/NACA0012/NACA0012-SP", delimiter = "\t");
    plot(XY[:,0]*chord,XY[:,1]*chord,'k',linewidth=3);

# UPPER = genfromtxt("TMP/308_2/T3/THERMO_SOLN_UPPER.out", delimiter = "\t");
# LOWER = genfromtxt("TMP/308_2/T3/THERMO_SOLN_LOWER.out", delimiter = "\t");
# BETA  = genfromtxt("TMP/308_2/T3/BETA.out", delimiter = "\t");
# LWC = 1.0e-3; Uinf = 102.8; 
# figure(2);
# idx = np.linspace(0,1,50);
# for j in range(0,50):
#     ind1 = j*1000;
#     ind2 = (j+1)*1000-1;
#     for jj in range(0,8):
#         LOWER[ind1:ind2,jj] = np.flipud(LOWER[ind1:ind2,jj]);
#     LOWER[ind1:ind2,0] = -LOWER[ind1:ind2,0];
#     subplot(311); plot(UPPER[ind1:ind2,0],UPPER[ind1:ind2,1],color=plt.cm.coolwarm(idx[j]),lw=2); 
#     plot(LOWER[ind1:ind2,0],LOWER[ind1:ind2,1],color=plt.cm.coolwarm(idx[j]),lw=2); 
#     plt.xlim([-0.06,0.03]); plt.gca().set_xticks([-0.04,-0.02,0,0.02]); 
#     plt.gca().set_yticks([0,4e-6,8e-6,12e-6]); plt.gca().set_yticklabels([0,4,8,12]);
#     plt.gca().set_ylabel(r'$h_f$ ($\mu$m)',fontsize=24);
#     subplot(312); plot(UPPER[ind1:ind2,0],UPPER[ind1:ind2,2],color=plt.cm.coolwarm(idx[j]),lw=2); 
#     plot(LOWER[ind1:ind2,0],LOWER[ind1:ind2,2],color=plt.cm.coolwarm(idx[j]),lw=2); 
#     plt.xlim([-0.06,0.03]); plt.gca().set_xticks([-0.04,-0.02,0,0.02]);
#     plt.ylim([-20,1.0]);    plt.gca().set_yticks([-20,0]);
#     plt.gca().set_ylabel(r'$T$ ($C$)', fontsize=24);
#     subplot(313); plot(UPPER[ind1:ind2,0],UPPER[ind1:ind2,3],color=plt.cm.coolwarm(idx[j]),lw=2); 
#     plot(LOWER[ind1:ind2,0],LOWER[ind1:ind2,3],color=plt.cm.coolwarm(idx[j]),lw=2);
#     plt.xlim([-0.06,0.03]); plt.gca().set_xticks([-0.04,-0.02,0,0.02]);
#     plt.gca().set_yticks([0,0.04,0.08,0.12]);
#     plt.gca().set_ylabel(r'$\dot{m}_{ice}$ ($kg/(m^2 s)$)',fontsize=24);
#     plt.gca().set_xlabel(r'$s$ ($m$)',fontsize=24);
# plot(BETA[:,0],LWC*Uinf*BETA[:,1],'--',c='k',lw=2);
# plt.tight_layout();

show()
#savefig('VALIDATIONS/IceShapeValidations.png');
