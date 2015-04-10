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
XY = genfromtxt("QuadTreeXY.dat", delimiter = '\t')
#plt.scatter(XY[:,0],XY[:,1])
for i in xrange(0,size(XY,0)/4):
    ind = 4*i
    X = np.append(XY[ind:ind+4,0], XY[ind,0]);
    Y = np.append(XY[ind:ind+4,1], XY[ind,1]);
    plt.plot(X,Y,color='k');
XYdata = genfromtxt("DataSet.dat", delimiter = "\t")
plt.scatter(XYdata[:,0],XYdata[:,1]);

#pylab.savefig('temp2.png',bbox_inches=0)

plt.show()

#import Image
#import numpy as np

#image=Image.open('temp2.png')
# image=Image.open('temp2.eps')
#image.load()

#print image.size

#crop_image = image.crop([370,220,1150,950])

#crop_image.save('temp2_cropped.png')
# crop_image.save('temp2_cropped.eps')
