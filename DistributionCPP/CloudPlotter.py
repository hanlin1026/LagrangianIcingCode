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
X = genfromtxt("CloudX.out", delimiter = '\n')
Y = genfromtxt("CloudY.out", delimiter = '\n')
XYa = genfromtxt("AirfoilXY.out", delimiter = "\t")
Xc = genfromtxt("CloudCELLX.out", delimiter = "\n")
Yc = genfromtxt("CloudCELLY.out", delimiter = "\n")
plt.scatter(X,Y,c="r",edgecolor='',lw=0,s=2)
#plt.scatter(Xc,Yc,edgecolor='',lw=0)
plt.plot(XYa[:,0],XYa[:,1])

plt.show()
