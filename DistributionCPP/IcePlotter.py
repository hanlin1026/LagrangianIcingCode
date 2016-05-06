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

# Length scaling and basedirectory
chord = 0.5334;
#chord = 21;
basedir = "/home/adegenna/LagrangianIcingCode/Validations/Ice/Run405Rime/";
# Import time history of ice accretion
# for i in range(0,7):
#     filename = basedir + "XY_NEW" + str(i+1) + ".out";
#     XY = genfromtxt(filename, delimiter = "\t")
#     figure(1);
#     plot(XY[:,0]*chord,XY[:,1]*chord,lw=3,c='b')

XY = genfromtxt("./Grid/NACA0012/NACA0012-SP", delimiter = "\t");
figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'k',linewidth=3);

# ******************************************************
# RUN 404
# ******************************************************

# XY = genfromtxt("./Grid/RUN404/T1/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);
# XY = genfromtxt("./Grid/RUN404/T2/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);
# XY = genfromtxt("./Grid/RUN404/T3/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);
# XY = genfromtxt("./Grid/RUN404/T4/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);
# XY = genfromtxt("./Grid/RUN404/T5/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);
# XY = genfromtxt("./Grid/RUN404/T6/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);
# XY = genfromtxt("./Grid/RUN404/T7/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);

# XY = genfromtxt("./Grid/RUN404/T_SIMUL_ROE1/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'m',linewidth=3);
# XY = genfromtxt("./Grid/RUN404/T_SIMUL_ROE2/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'m',linewidth=3);
# XY = genfromtxt("./Grid/RUN404/T_SIMUL_ROE3/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'m',linewidth=3);
# XY = genfromtxt("./Grid/RUN404/T_SIMUL_ROE4/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'m',linewidth=3);
# XY = genfromtxt("./Grid/RUN404/T_SIMUL_ROE5/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'m',linewidth=3);
# XY = genfromtxt("./Grid/RUN404/T_SIMUL_ROE6/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'m',linewidth=3);
# XY = genfromtxt("./Grid/RUN404/T_SIMUL_ROE7/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'m',linewidth=3);

# ******************************************************
# RUN 308
# ******************************************************

# XY = genfromtxt("./Grid/RUN308/T1/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'r',linewidth=3);
# XY = genfromtxt("./Grid/RUN308/T2/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'r',linewidth=3);
# XY = genfromtxt("./Grid/RUN308/T3/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'r',linewidth=3);
# XY = genfromtxt("./Grid/RUN308/T4/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'r',linewidth=3);

# XY = genfromtxt("./Grid/RUN308/T_ROE1/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);
# XY = genfromtxt("./Grid/RUN308/T_ROE2/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);
# XY = genfromtxt("./Grid/RUN308/T_ROE3/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);
# XY = genfromtxt("./Grid/RUN308/T_ROE4/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);

# XY = genfromtxt("./Grid/RUN308/T_CFPOS_ROE1/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'m',linewidth=3);


# XY = genfromtxt("./Grid/RUN308/T_SIMUL_ROE1/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);
# XY = genfromtxt("./Grid/RUN308/T_SIMUL_ROE2/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);
# XY = genfromtxt("./Grid/RUN308/T_SIMUL_ROE3/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);
# XY = genfromtxt("./Grid/RUN308/T_SIMUL_ROE4/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);

# ******************************************************
# RUN 405
# ******************************************************

# XY = genfromtxt("./Grid/RUN405/T1/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'g',linewidth=3);
# XY = genfromtxt("./Grid/RUN405/T2/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'g',linewidth=3);
# XY = genfromtxt("./Grid/RUN405/T3/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'g',linewidth=3);
# XY = genfromtxt("./Grid/RUN405/T4/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'g',linewidth=3);
# XY = genfromtxt("./Grid/RUN405/T5/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'g',linewidth=3);
# XY = genfromtxt("./Grid/RUN405/T6/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'g',linewidth=3);
# XY = genfromtxt("./Grid/RUN405/T7/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'g',linewidth=3);

# XY = genfromtxt("./Grid/RUN405/T_SIMUL_ROE1/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'m',linewidth=3);
# XY = genfromtxt("./Grid/RUN405/T_SIMUL_ROE2/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'m',linewidth=3);
# XY = genfromtxt("./Grid/RUN405/T_SIMUL_ROE3/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'m',linewidth=3);
# XY = genfromtxt("./Grid/RUN405/T_SIMUL_ROE4/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'m',linewidth=3);
# XY = genfromtxt("./Grid/RUN405/T_SIMUL_ROE5/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'m',linewidth=3);
# XY = genfromtxt("./Grid/RUN405/T_SIMUL_ROE6/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'m',linewidth=3);
# XY = genfromtxt("./Grid/RUN405/T_SIMUL_ROE7/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'m',linewidth=3);

# ******************************************************
# RUN 402
# ******************************************************

# XY = genfromtxt("./Grid/RUN402/T_SIMUL_ROE1/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);
# XY = genfromtxt("./Grid/RUN402/T_SIMUL_ROE2/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);
# XY = genfromtxt("./Grid/RUN402/T_SIMUL_ROE3/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);
# XY = genfromtxt("./Grid/RUN402/T_SIMUL_ROE4/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);
# XY = genfromtxt("./Grid/RUN402/T_SIMUL_ROE5/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);
# XY = genfromtxt("./Grid/RUN402/T_SIMUL_ROE6/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);
# XY = genfromtxt("./Grid/RUN402/T_SIMUL_ROE7/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);

# ******************************************************
# RUN 409
# ******************************************************

# XY = genfromtxt("./Grid/RUN409/T_SIMUL_ROE1/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);
# XY = genfromtxt("./Grid/RUN409/T_SIMUL_ROE2/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);
# XY = genfromtxt("./Grid/RUN409/T_SIMUL_ROE3/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);
# XY = genfromtxt("./Grid/RUN409/T_SIMUL_ROE4/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);
# XY = genfromtxt("./Grid/RUN409/T_SIMUL_ROE5/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);
# XY = genfromtxt("./Grid/RUN409/T_SIMUL_ROE6/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);

# XY = genfromtxt("./Grid/RUN409/T_singleShot1/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'m',linewidth=3);
# XY = genfromtxt("./Grid/RUN409/T_singleShot2/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'m',linewidth=3);
# XY = genfromtxt("./Grid/RUN409/T_singleShot3/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'m',linewidth=3);
# XY = genfromtxt("./Grid/RUN409/T_singleShot4/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'m',linewidth=3);
# XY = genfromtxt("./Grid/RUN409/T_singleShot5/XY_NEW.out", delimiter = "\t");
# figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'m',linewidth=3);

# ******************************************************
# TEMPORARY RUNS
# ******************************************************

for k in range(1,6):
    for i in range(1,6):
        XY = genfromtxt("DAKOTA/workdir." + str(i) + "/T" + str(k) + "/XY_NEW.out", delimiter = "\t");
        figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,color=cm.coolwarm(0.2*i),linewidth=3);


plt.xlim([-0.05,0.2])
axis('equal')
plt.grid(b=True)

# Compare to LEWICE results
# RUN404 = genfromtxt("/home/adegenna/LagrangianIcingCode/Validations/LewiceIceshapes/Run404.csv", delimiter = ",");
# plt.scatter(RUN404[:,0]/21.0*chord,RUN404[:,1]/21.0*chord,c='b',s=50);
# RUN308 = genfromtxt("/home/adegenna/LagrangianIcingCode/Validations/LewiceIceshapes/Run308.csv", delimiter = ",");
# plt.scatter(RUN308[:,0]/21.0*chord,RUN308[:,1]/21.0*chord,c='r',s=50);
# RUN3082 = genfromtxt("/home/adegenna/LagrangianIcingCode/Validations/LewiceIceshapes/Run308Lewice.csv", delimiter = ", ");
# plt.plot(RUN3082[:,0]/21.0*chord,RUN3082[:,1]/21.0*chord,c='g',lw=3);
# RUN405 = genfromtxt("/home/adegenna/LagrangianIcingCode/Validations/LewiceIceshapes/Run405.csv", delimiter = ",");
# plt.scatter(RUN405[:,0]/21.0*chord,RUN405[:,1]/21.0*chord,c='g',s=50);
# RUN402 = genfromtxt("/home/adegenna/LagrangianIcingCode/Validations/LewiceIceshapes/Run402.csv", delimiter = ",");
# plt.scatter(RUN402[:,0]/21.0*chord,RUN402[:,1]/21.0*chord,c='b',s=50);
# RUN409 = genfromtxt("/home/adegenna/LagrangianIcingCode/Validations/LewiceIceshapes/Run409.csv", delimiter = ",");
# plt.scatter(RUN409[:,0]/21.0*chord,RUN409[:,1]/21.0*chord,c='b',s=50);

# UPPER = genfromtxt("./Grid/RUN409/T_singleShot0/THERMO_SOLN_UPPER.out", delimiter = "\t");
# LOWER = genfromtxt("./Grid/RUN409/T_singleShot0/THERMO_SOLN_LOWER.out", delimiter = "\t");
# BETA = genfromtxt("./Grid/RUN409/T_singleShot0/BETA.out", delimiter = "\t");
# LWC = 1.3e-3; Uinf = 67.1; 
# figure(2);
# subplot(311); plot(UPPER[:,0],UPPER[:,1],'b.-'); plot(LOWER[:,0],LOWER[:,1],'r.-'); #plt.xlim([-0.06,0.03])
# subplot(312); plot(UPPER[:,0],UPPER[:,2],'b.-'); plot(LOWER[:,0],LOWER[:,2],'r.-'); #plt.xlim([-0.06,0.03])
# subplot(313); plot(UPPER[:,0],UPPER[:,3],'b.-'); plot(LOWER[:,0],LOWER[:,3],'r.-'); plot(BETA[:,0],LWC*Uinf*BETA[:,1],'--',c='k'); #plt.xlim([-0.06,0.03])

show()
