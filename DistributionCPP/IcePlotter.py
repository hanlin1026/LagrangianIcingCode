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

XY = genfromtxt("./Grid/RUN404/T1area2/XY_NEW.out", delimiter = "\t");
figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);
XY = genfromtxt("./Grid/RUN404/T2area2/XY_NEW.out", delimiter = "\t");
figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);
XY = genfromtxt("./Grid/RUN404/T3area2/XY_NEW.out", delimiter = "\t");
figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);
XY = genfromtxt("./Grid/RUN404/T4area2/XY_NEW.out", delimiter = "\t");
figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);
XY = genfromtxt("./Grid/RUN404/T5area2/XY_NEW.out", delimiter = "\t");
figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);
XY = genfromtxt("./Grid/RUN404/T6area2/XY_NEW.out", delimiter = "\t");
figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);
XY = genfromtxt("./Grid/RUN404/T7area2/XY_NEW.out", delimiter = "\t");
figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'b',linewidth=3);

XY = genfromtxt("./Grid/RUN308/T1area2/XY_NEW.out", delimiter = "\t");
figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'r',linewidth=3);
XY = genfromtxt("./Grid/RUN308/T2area2/XY_NEW.out", delimiter = "\t");
figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'r',linewidth=3);
XY = genfromtxt("./Grid/RUN308/T3area2/XY_NEW.out", delimiter = "\t");
figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'r',linewidth=3);
XY = genfromtxt("./Grid/RUN308/T4area2/XY_NEW.out", delimiter = "\t");
figure(1); plot(XY[:,0]*chord,XY[:,1]*chord,'r',linewidth=3);



# Current iteration
# XY = genfromtxt("XY_NEW.out", delimiter="\t");
# plot(XY[:,0]*chord,XY[:,1]*chord,lw=3,c='r');
plt.xlim([-0.05,0.2])
axis('equal')
plt.grid(b=True)
#legend(['230 K','240 K','250 K','260 K','270 K','NACA0012'])
# Compare to LEWICE results
RUN404 = genfromtxt("/home/adegenna/LagrangianIcingCode/Validations/LewiceIceshapes/Run404.csv", delimiter = ",");
plt.scatter(RUN404[:,0]/21.0*chord,RUN404[:,1]/21.0*chord,c='b',s=50);
RUN308 = genfromtxt("/home/adegenna/LagrangianIcingCode/Validations/LewiceIceshapes/Run308.csv", delimiter = ",");
plt.scatter(RUN308[:,0]/21.0*chord,RUN308[:,1]/21.0*chord,c='r',s=50);

#habashi = genfromtxt("/home/adegenna/LagrangianIcingCode/Validations/LewiceIceshapes/Habashi405.csv",delimiter=",");
#plt.scatter(habashi[:,0]/21.0*chord,habashi[:,1]/21.0*chord,c='r');

UPPER = genfromtxt("THERMO_SOLN_UPPER.out", delimiter = "\t");
LOWER = genfromtxt("THERMO_SOLN_LOWER.out", delimiter = "\t");
BETA = genfromtxt("Grid/RUN404/T4/BETA.out", delimiter = "\t");
LWC = 0.55e-3; Uinf = 102.8; 
figure(2);
subplot(311); plot(UPPER[:,0],UPPER[:,1],'b.-'); plot(LOWER[:,0],LOWER[:,1],'r.-'); plt.xlim([-0.1,0.1])
subplot(312); plot(UPPER[:,0],UPPER[:,2],'b.-'); plot(LOWER[:,0],LOWER[:,2],'r.-'); plt.xlim([-0.1,0.1])
subplot(313); plot(UPPER[:,0],UPPER[:,3],'b.-'); plot(LOWER[:,0],LOWER[:,3],'r.-'); plot(BETA[:,0],LWC*Uinf*BETA[:,1],'--',c='k'); plt.xlim([-0.1,0.1])

show()
