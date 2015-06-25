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
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.mlab as mlab

inches_per_pt = 1.0/72.27
width = 1150
height = 800
fig_size = [width*inches_per_pt,height*inches_per_pt]

params = {   #'axes.labelsize': 30,
             #'text.fontsize': 40,
             #'legend.fontsize': 20,
             'xtick.labelsize': 40,
             'ytick.labelsize': 40,             
             'figure.figsize':fig_size,
             #'figure.markersize': 50}
}
pylab.rcParams.update(params)

fig = plt.figure()
ax = plt.gca()
#ax = fig.add_subplot(111)

X = genfromtxt("Xgrid.out", delimiter = '\n');
Y = genfromtxt("Ygrid.out", delimiter = '\n');
Z = genfromtxt("Zgrid.out", delimiter = '\n');

ax.set_xlim([0,1]);
ax.set_ylim([0,1]);
ax.set_xlabel(r'$X$',fontsize=50);
ax.set_ylabel(r'$Y$',fontsize=50);
ax.set_aspect('equal');

Zgrid = mlab.griddata(X,Y,Z,X,Y);
plt.contourf(X,Y,Zgrid);


#plt.contour(POD1x,POD2y,PODZ)
# im = plt.imshow(PODZ, interpolation='bilinear', origin='lower',
#                 cmap=cm.coolwarm, extent=(-2,2,-2,2))
# levels = np.arange(0.95, 1.02, 0.01)
# CS = plt.contour(PODZ,
#                  origin='lower',
#                  linewidths=5,
#                  extent=(-2,2,-2,2),colors='k')
# ax.set_xticks([-2,-1,0,1,2]);
# ax.set_yticks([-2,-1,0,1,2]);
# ax.set_xlabel('POD 1',fontsize=50)
# ax.set_ylabel('POD 2',fontsize=50)
# # We can still add a colorbar for the image, too.
# CBI = plt.colorbar(im, orientation='vertical',ticks=[PODZ.min(),PODZ.max()])


#plt.plot([0,100],[0,0],'k--',linewidth=3)
#plt.plot(Span,coeffs[0,:],'b-o',linewidth=5,markersize=12,label="POD Mode 1")
#plt.plot(Span,coeffs[1,:],'r-o',linewidth=5,markersize=12,label="POD Mode 2")
#plt.legend()

# for i in xrange(0, 10):
#    plt.plot(i+1,XY[i],'bo',markeredgecolor='black',markersize=15)
# lns1 = ax.plot(xGLOB,MEAN,color='black',linewidth=5,label="Mean")
# ax2 = ax.twinx()
# lns2 = ax2.plot(xGLOB,POD[:,0],color='blue',linewidth=5,label="POD Mode 1")
# lns3 = ax2.plot(xGLOB,POD[:,1],'--',color='blue',linewidth=5,label="POD Mode 2")
# for tl in ax2.get_yticklabels():
#    tl.set_color('b')
# ax2.set_ylim([-0.15, 0.15])
# ax.set_xlabel('$s$',fontsize=50)
# ax.set_ylabel('$N(s)$',fontsize=50)
# ax2.set_ylabel('$N(s)$',fontsize=50,color='blue')
# lns = lns1+lns2+lns3
# labs = [l.get_label() for l in lns]
# ax.legend(lns, labs, loc=0, prop={'size':30})
# ax.set_ylim([0, 2])
# ax.set_xlim([-3,6])
# ax.set_xticks([-3,0,3,6])
# ax.set_xticklabels([-3,0,3,6])
# ax.set_yticks([0, 0.5, 1, 1.5, 2])

# ax.set_xlim([0.5, 10.5])
# ax.set_ylim([2, 50])
# ax.set_yticklabels(['2e0', '1e1', '5e1'])
# ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

#Colorbar code ************
# min, max = (0,100)
# step = 1

# #Using contourf to provide my colorbar info, then clearing the figure
# Z = [[0,0],[0,0]]
# levels = range(min,max+step,step)
# CS3 = plt.contourf(Z, levels, cmap=cm.coolwarm)
# plt.clf()

# #Plotting what I actually want
# plt.xlabel(r'$Y/c$',fontsize=50)
# plt.ylabel(r'$X/c$',fontsize=50)
# for i in xrange(0, 10):
#     plt.plot(IceX[:,i], IceY[:,i], color=cmap(N[i]), linewidth=2)
# #plt.plot(IceXScaleH[:,8], IceYScaleH[:,8], color=cmap(1.), linewidth=2)
# plt.plot(mean[:,0], mean[:,1], 'k--', linewidth = 2)
# plt.plot(clean[:,0], clean[:,1], 'k-', linewidth = 2)
# #MEAN = genfromtxt("MeanShiftedIce.dat", delimiter = ',')
# #plt.plot(xGLOB, MEAN, '--', color='black', linewidth=2)
# ax2 = plt.gca()
# ax2.set_xlim([-0.015, 0.040])
# ax2.set_ylim([-0.020, 0.020])
# ax2.set_aspect('equal')
# ax2.xaxis.set_ticks([0, 0.02, 0.04])
# ax2.yaxis.set_ticks([-0.02, 0, 0.02])
# # Inset subplot
# # eps1 = [-2,-2,-2,0,0,0,2,2,2]
# # eps2 = [-2,0,2,-2,0,2,-2,0,2]
# # a = axes([.55, .352, .4, .4], axisbg='0.75')
# # plot(PODcoords[0,:], PODcoords[1,:],'o',color='k')
# # for i in xrange(0,8):
# #    plot(eps1[i],eps2[i],'o',color=cmap(N[i]),markersize=15)
# # plot(eps1[8],eps2[8],'o',color=cmap(N[8]),markersize=15)
# # title('POD Coordinates',fontsize=25)
# # a.set_aspect('equal')
# # setp(a, xticks=[-10,-5, 0,5, 10], yticks=[-10,-5, 0,5, 10])
# # for tick in a.xaxis.get_major_ticks():
# #    tick.label.set_fontsize(25)
# # for tick in a.yaxis.get_major_ticks():
# #    tick.label.set_fontsize(25) 

# cbar = plt.colorbar(CS3, ticks=[0, 100],shrink=0.77)
# cbar.ax.set_yticklabels(['0.75', '1.25'])
#*********************

# plt.xlabel(r'Mode 1',fontsize=25)
# plt.ylabel(r'Mode 2',fontsize=25)
# plt.grid()


# loc parameter: upper right = 1; upper left = 2
#plt.legend(["gPC", "Monte Carlo"], loc=2, prop={'size':30})
ax.set_axisbelow(True)
#ax.set_aspect('equal')
gcf().subplots_adjust(left=0.16)
gcf().subplots_adjust(bottom=0.15)
ax.xaxis.labelpad = 20
#ax.yaxis.labelpad = 20
gcf().tight_layout()

#plt.grid(None)

pylab.savefig('2Dplot.eps',bbox_inches='tight')

#plt.gca().tight_layout()
plt.draw()
plt.show()

#import Image
#import numpy as np

#image=Image.open('temp2.png')
#image.load()

#print image.size

#crop_image = image.crop([370,220,1150,950])

#crop_image.save('temp2_cropped.png')
