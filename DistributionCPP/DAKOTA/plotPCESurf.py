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
width = 162*5
height = 100*5
fig_size = [width*inches_per_pt,height*inches_per_pt]

params = {   'axes.labelsize': 30,
             'text.fontsize': 50,
             #'legend.fontsize': 20,
             'xtick.labelsize': 30,
             'ytick.labelsize': 30,             
             'figure.figsize':fig_size,
             #'figure.markersize': 50}
}
pylab.rcParams.update(params)

# **********************************
# SUBROUTINE TO PARSE DAKOTA FILE
# **********************************

def parseDakotaOutfile(filename):
    
    # Open Dakota.out file and find where the PCE coefficients for output_1 begin
    datafile = open(filename,'r');
    found = False;
    iter = 0;
    while True:
        line = datafile.readline();
        if 'Coefficients of Polynomial Chaos Expansion' in line:
            found = True;
            break;
    for i in range(1,4):
        line = datafile.readline();

    # Start reading in output_1 PCE coefficients
    coeffs = [];
    order1 = [];
    order2 = [];
    while True:
        if 'Coefficients of Polynomial Chaos Expansion' in line:
            break;
        else:
            linedata = line.split();
            coeffs.append( float(linedata[0]) );
            order1.append( int(re.findall('\d+',linedata[1])[0]) );
            order2.append( int(re.findall('\d+',linedata[2])[0]) );
            line = datafile.readline();
    COEFF1 = [coeffs, order1, order2];

    # Start reading in output_2 PCE coefficients
    for i in range(1,4):
        line = datafile.readline();
    coeffs = [];
    order1 = [];
    order2 = [];
    while True:
        if '-----------------------------------------------------------------------------' in line:
            break;
        else:
            linedata = line.split();
            coeffs.append( float(linedata[0]) );
            order1.append( int(re.findall('\d+',linedata[1])[0]) );
            order2.append( int(re.findall('\d+',linedata[2])[0]) );
            line = datafile.readline();
    COEFF2 = [coeffs, order1, order2];

    return np.squeeze(np.array([[COEFF1, COEFF2]]));

# **********************************
# SUBROUTINE FOR LEGENDRE POLYNOMIALS
# **********************************

def legendre(x,N):
    M = np.size(x);
    Y = np.zeros(M);
    L = np.zeros((N+1,M));
    
    if (N==0):
        for i in range(0,M):
            Y[i] = 1.0;
    elif (N==1):
        for i in range(0,M):
            Y[i] = x[i];
    else:
        L[0,:] = 1.0;
        L[1,:] = x;
        for i in range(2,N+1):
                L[i,:] = (1.0/i)*((2*i-1)*x[:]*L[i-1,:] - (i-1)*L[i-2,:]);
        for i in range(0,M):
            Y[i] = L[N,i];
    
    return Y;

# **********************************
# SUBROUTINE FOR SAMPLING PCE
# **********************************

def samplePCE(coeffs,x):
    # coeffs = [C_ij, i, j]
    d      = np.shape(x)[0];
    N      = np.shape(x)[1];
    Y      = np.zeros(N);
    order  = np.shape(coeffs)[1];
    for i in range(0,order):
        Pi = np.ones(N);
        for j in range(0,d):
            Pi = Pi*legendre(x[j,:],int(coeffs[j+1,i]));
        Y = Y + coeffs[0,i]*Pi;

    return Y;

# **********************************
# MAIN SCRIPT
# **********************************

# Load PCE coefficients
basedir = "/home/adegenna/LagrangianIcingCode/DistributionCPP/DAKOTA";
PCEfile = basedir + "/ICE.out";
coeffs  = parseDakotaOutfile(PCEfile);

# Load clean airfoil coordinates
xy0 = genfromtxt(basedir + "/../Grid/NACA0012/NACA0012-SP", delimiter = "\t");

# Load quadrature point
data = genfromtxt(basedir + "/ICE.dat");
xq   = data[1:,1];
yq   = data[1:,2];
zq   = data[1:,3];

# Test legendre function
x       = np.linspace(-1,1,100);
[XX,YY] = np.meshgrid(x,x);
X       = np.squeeze(np.array([[np.hstack(XX),np.hstack(YY)]]));
Cij     = np.squeeze(coeffs[0,:,:]);
y       = samplePCE(Cij,X);

# Contour plot of 2-D PCE map
# fig  = plt.figure();
# ax   = fig.gca(projection='3d');
# surf = ax.plot_surface((XX+1)/2*(270-250)+250, (YY+1)/2*(1-0.3)+0.3, y.reshape(100,100), rstride=1, cstride=1,
#                        cmap=cm.coolwarm,vmin=0.34,vmax=0.44,
#                        linewidth=0, antialiased=False);
# ax.set_xlim([250,270]);
# ax.set_ylim([0.3,1.0]);
# ax.set_zlim([0.35,0.44]);
# ax.view_init(elev=30, azim=120)
figure();
plt.contourf((XX+1)/2*(270-250)+250, (YY+1)/2*(1-0.3)+0.3, y.reshape(100,100),levels=np.linspace(0.32,0.50,100),
             cmap=cm.coolwarm,vmin=0.34,vmax=0.44,
             linewidth=0, antialiased=False);
plt.colorbar(ticks=[0.32,0.50])
plt.xticks([250,260,270])
plt.yticks([0.3,0.65,1.0])
plt.xlabel(r'$T_{\infty}$ (K)');
plt.ylabel(r'LWC ($g/m^3$)');
plt.tight_layout()

# Plot quadrature points
# fig2 = plt.figure();
# ax2  = fig2.gca(projection='3d');
# ax2.scatter(xq,yq*(1e3),zq,s=50,c=zq,depthshade=False,cmap=cm.coolwarm,vmin=0.34,vmax=0.44);
# ax2.set_xlim([250,270]);
# ax2.set_ylim([0.3,1.0]);
# ax2.set_zlim([0.35,0.44]);
# ax2.view_init(elev=30, azim=120)
figure();
plt.scatter(xq,yq*(1e3),s=50,c=zq,cmap=cm.coolwarm,vmin=0.34,vmax=0.44);
#plt.colorbar(ticks=[0.32,0.50])
gca().set_xlim([250,270]);
gca().set_ylim([0.3,1.0]);
plt.xticks([250,260,270])
plt.yticks([0.3,0.65,1.0])
plt.xlabel(r'$T_{\infty}$ (K)');
plt.ylabel(r'LWC ($g/m^3$)');
plt.tight_layout()

# Approximate statistics
xsamp = np.random.uniform(-1,1,1e5);
ysamp = np.random.uniform(-1,1,1e5);
XY    = np.squeeze(np.array([[xsamp],[ysamp]]));
Z     = samplePCE(Cij,XY);

# Plot statistics histogram
fig3 = plt.figure();
plt.hist(Z,100,normed=1,facecolor='b',alpha=1);
gca().set_xlim([0.32,0.48]);
plt.xticks([0.32,0.40,0.48])
plt.xlabel(r'$C_L$');
plt.ylabel(r'PDF($C_L$)');
plt.tight_layout()

# Divide shapes into groups
zhigh = []; zlow = [];
limhigh = np.percentile(Z,75);
limlow  = np.percentile(Z,25);
print(limhigh); print(limlow);
figure(12); plt.scatter(xq,yq*(1e3),s=50,facecolors='none',edgecolors='k'); 
gca().set_xlim([250,270]); gca().set_ylim([0.3,1.0]);
plt.xticks([250,260,270]); plt.yticks([0.3,0.65,1.0])
plt.xlabel(r'$T_{\infty}$ (K)');
plt.ylabel(r'LWC ($g/m^3$)');
plt.tight_layout()
for i in range(0,np.size(zq)):
    if (zq[i] >= limhigh):
        xy = genfromtxt(basedir + '/RESULTS/workdir.' + str(i+1) + '/T5/XY_NEW.out');
        figure(11); plt.plot(xy[:,0],xy[:,1],'r',lw=1); 
        figure(12); plt.scatter(xq[i],yq[i]*(1e3),s=50,c='r');
    elif (zq[i] <= limlow):
        xy = genfromtxt(basedir + '/RESULTS/workdir.' + str(i+1) + '/T5/XY_NEW.out');
        figure(11); plt.plot(xy[:,0],xy[:,1],'b',lw=1);
        figure(12); plt.scatter(xq[i],yq[i]*(1e3),s=50,c='b');

figure(11); plot(xy0[:,0],xy0[:,1],'k',lw=2);
figure(11); axis('equal'); plt.xlim([-0.05,0.05]);
plt.xlabel(r'$X/c$'); plt.ylabel(r'$Y/c$'); plt.tight_layout()

show();
