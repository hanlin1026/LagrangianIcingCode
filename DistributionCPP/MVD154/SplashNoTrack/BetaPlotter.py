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
width = 1000
height = 800
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
fig = plt.figure()
ax = fig.gca()
sizeX = 35;
bins = 27;
X = numpy.zeros((sizeX,bins));
Y = numpy.zeros((sizeX,bins));
weight = genfromtxt("Weights.dat", delimiter = '\n');
xyExp = genfromtxt("BetaMVD154Experiment.dat", delimiter = ',');
xInterp = numpy.linspace(-0.3,0.2,1000);
yInterp = numpy.zeros(1000);
# Interpolate beta
for i in range(0,27):
    xName = "workdir." + str(i+1) + "/BetaBins.out";
    yName = "workdir." + str(i+1) + "/Beta.out";
    x = genfromtxt(xName, delimiter = '\n');
    y = genfromtxt(yName, delimiter = '\n');
    yInterp = yInterp + weight[i]*numpy.interp(xInterp,x,y,left=0,right=0)/100.0;
plt.plot(xInterp,yInterp,lw=5);
plt.plot(xyExp[:,0],xyExp[:,1],c='r',lw=5);
legend(["Computational","Experimental"])
# Write to file
XY = np.vstack((xInterp,yInterp));
XY = XY.transpose();
np.savetxt('BetaXY.dat', XY, delimiter=',')

#plt.scatter(X,Y,c="r",edgecolor='',lw=0,s=15)
#plt.scatter(Xc,Yc,edgecolor='',lw=0,s=15,c='g')
#axis('equal')

plt.show()
