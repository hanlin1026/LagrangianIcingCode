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
height = 500
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

XY = genfromtxt("XY_OLD_NEW.out", delimiter = "\t")
figure(1);
plot(XY[:,0],XY[:,1],c='k')
plot(XY[:,2],XY[:,3],c='r')
axis('equal')

UPPER = genfromtxt("THERMO_SOLN_UPPER.out", delimiter = "\t");
LOWER = genfromtxt("THERMO_SOLN_LOWER.out", delimiter = "\t");
BETA = genfromtxt("ThermoEqns/BetaXY.dat", delimiter = ",");
LWC = 0.55e-3; Uinf = 100; 
figure(2);
subplot(311); plot(UPPER[:,0],UPPER[:,1],c='b'); plot(LOWER[:,0],LOWER[:,1],c='r'); 
subplot(312); plot(UPPER[:,0],UPPER[:,2],c='b'); plot(LOWER[:,0],LOWER[:,2],c='r');
subplot(313); plot(UPPER[:,0],UPPER[:,3],c='b'); plot(LOWER[:,0],LOWER[:,3],c='r'); plot(BETA[:,0],LWC*Uinf*BETA[:,1],'--',c='k')

show()
