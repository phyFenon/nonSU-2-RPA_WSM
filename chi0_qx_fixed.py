#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 16:14:24 2019
@author: feng
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.colors as colors
import matplotlib
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
#from matplotlib.mlab import griddata

#
cm='jet'
cm2='rainbow'
cmap = plt.get_cmap(cm2)
gamma=0

N=32 #beta=100 eta=0.01
NX=16   #1 2 3 4 (-pi/2 0 pi/2 pi) 
print(NX)
data=np.genfromtxt("gama%s_full_real.txt"%(gamma))


x=data[:,0]/np.pi
y=data[:,1]/np.pi
z=data[:,2]/np.pi

z1=np.real(data[:,3])
z2=np.real(data[:,4])
z3=np.real(data[:,5])

z4=np.real(data[:,6])
z5=np.real(data[:,7])
z6=np.real(data[:,8])

X=np.reshape(x,(N,N,N))
Y=np.reshape(y,(N,N,N))
Z=np.reshape(z,(N,N,N))
Z1=np.reshape(z1,(N,N,N))
Z2=np.reshape(z2,(N,N,N))
Z3=np.reshape(z3,(N,N,N))

Z4=np.reshape(z4,(N,N,N))
Z5=np.reshape(z5,(N,N,N))
Z6=np.reshape(z6,(N,N,N))

fig=plt.figure(figsize=(6,6),dpi=400)

xy_labelsize=18 

my_x_ticks=[-1,-0.5,0,0.5,1]
my_y_ticks=[-1,-0.5,0,0.5,1]
ticks1=('','-0.5','0','0.5','1')
ticks2=('-1','-0.5','0','0.5','')


index=-1.0+2*NX/N
Yaxis=Y[NX-1,:,:]
Zaxis=Z[NX-1,:,:]
Z1axis=Z1[NX-1,:,:]
Z2axis=Z2[NX-1,:,:]
Z3axis=Z3[NX-1,:,:]

Z4axis=Z4[NX-1,:,:]
Z5axis=Z5[NX-1,:,:]
Z6axis=Z6[NX-1,:,:]
print(Z1.max(),Z1.min())
print(Z2.max(),Z2.min())
print(Z3.max(),Z3.min())

#this is for qx=constant,qy&qz
#ax1
levels = MaxNLocator(nbins=200).tick_values(Z1axis.min(), Z1axis.max())
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
ax1=fig.add_subplot(321)

c1 = ax1.contourf(Yaxis, Zaxis, Z1axis,cmap=cmap,levels=levels, norm=norm)

##ax1.set_title(r'$q_x=0$',FONTSIZE=18)

ax1.set_ylabel(r'$q_{z}/\pi$',FONTSIZE=xy_labelsize)
plt.xlim(-1,1)
plt.ylim(-1,1)
plt.yticks(my_y_ticks)
ax1.set_xticklabels([])


ax1.xaxis.set_major_locator(MultipleLocator(0.5))
ax1.yaxis.set_major_locator(MultipleLocator(0.5))
ax1.tick_params(which='major',direction='in',top='True',right='True',length=8,width=2,labelsize=14)


ax1.spines['top'].set_linewidth(2)
ax1.spines['bottom'].set_linewidth(2)
ax1.spines['left'].set_linewidth(2)
ax1.spines['right'].set_linewidth(2)

ticks=np.linspace(Z1axis.min(),Z1axis.max(),2,endpoint=True)
cbar=plt.colorbar(c1,ax=ax1,ticks=ticks,shrink=0.75)
cbar.ax.set_yticklabels(["{:.2f}".format(i) for i in ticks])
cbar.ax.tick_params(labelsize=14)

##ax1.text(-0.9, 0.8,'$(a1)$',color='white', fontsize=15)
ax1.text(-0.1, 0.7,r'$\chi^{0}_{dd}$',color='black', fontsize=16)




#ax2
levels = MaxNLocator(nbins=200).tick_values(Z2axis.min(), Z2axis.max())
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

ax2=fig.add_subplot(323)
c2= ax2.contourf(Yaxis, Zaxis, Z2axis,cmap=cmap,levels=levels, norm=norm)
plt.xlim(-1,1)
plt.ylim(-1,1)
ax2.set_ylabel(r'$q_{z}/\pi$',FONTSIZE=xy_labelsize)

ax2.set_xticklabels([])
ticks=np.arange(-1,1,0.5)
plt.yticks(ticks)

ax2.xaxis.set_major_locator(MultipleLocator(0.5))
ax2.yaxis.set_major_locator(MultipleLocator(0.5))
ax2.tick_params(which='major',direction='in',top='True',right='True',length=8,width=2,labelsize=14)
ax2.spines['top'].set_linewidth(2)
ax2.spines['bottom'].set_linewidth(2)
ax2.spines['left'].set_linewidth(2)
ax2.spines['right'].set_linewidth(2)

ticks=np.linspace(Z2axis.min(),Z2axis.max(),2,endpoint=True)
cbar=plt.colorbar(c2,ax=ax2,ticks=ticks,shrink=0.75)
cbar.ax.set_yticklabels(["{:.2f}".format(i) for i in ticks])
cbar.ax.tick_params(labelsize=14)

ax2.text(-0.1, 0.7,r'$\chi^{0}_{zz}$',color='black', fontsize=16)

#ax3
ax3=fig.add_subplot(325)
levels = MaxNLocator(nbins=200).tick_values(Z3axis.min(), Z3axis.max())
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
c3 = ax3.contourf(Yaxis, Zaxis, Z3axis,cmap=cmap,levels=levels, norm=norm)
plt.xlim(-1,1)
plt.ylim(-1,1)
ax3.set_xlabel(r'$q_{y}/\pi$',FONTSIZE=xy_labelsize)
ax3.set_ylabel(r'$q_{z}/\pi$',FONTSIZE=xy_labelsize)

ax3.xaxis.set_major_locator(MultipleLocator(0.5))
ax3.yaxis.set_major_locator(MultipleLocator(0.5))
ax3.tick_params(which='major',direction='in',top='True',right='True',length=8,width=2,labelsize=14)
ax3.spines['top'].set_linewidth(2)
ax3.spines['bottom'].set_linewidth(2)
ax3.spines['left'].set_linewidth(2)
ax3.spines['right'].set_linewidth(2)
ticks=np.arange(-1,1,0.5)
plt.xticks(my_x_ticks)
plt.yticks(ticks)

ticks=np.linspace(Z3axis.min(),Z3axis.max(),2,endpoint=True)
cbar=plt.colorbar(c3,ax=ax3,ticks=ticks,shrink=0.75)
cbar.ax.set_yticklabels(["{:.2f}".format(i) for i in ticks])
cbar.ax.tick_params(labelsize=14)
ax3.text(-0.1, 0.7,r'$\chi^{0}_{yy}$',color='black', fontsize=16)


###second column
#ax1
levels = MaxNLocator(nbins=200).tick_values(Z4axis.min(), Z4axis.max())
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
ax1=fig.add_subplot(322)

c1 = ax1.contourf(Yaxis, Zaxis, Z4axis,cmap=cmap,levels=levels, norm=norm)

##ax1.set_title(r'$q_x=0$',FONTSIZE=18)

##ax1.set_ylabel(r'$q_{z}/\pi$',FONTSIZE=xy_labelsize)
plt.xlim(-1,1)
plt.ylim(-1,1)
ax1.set_xticklabels([])
ax1.set_yticklabels([])

ax1.xaxis.set_major_locator(MultipleLocator(0.5))
ax1.yaxis.set_major_locator(MultipleLocator(0.5))
ax1.tick_params(which='major',direction='in',top='True',right='True',length=8,width=2,labelsize=14)


ax1.spines['top'].set_linewidth(2)
ax1.spines['bottom'].set_linewidth(2)
ax1.spines['left'].set_linewidth(2)
ax1.spines['right'].set_linewidth(2)

ticks=np.linspace(Z4axis.min(),Z4axis.max(),2,endpoint=True)
cbar=plt.colorbar(c1,ax=ax1,ticks=ticks,shrink=0.75)
cbar.ax.set_yticklabels(["{:.2f}".format(i) for i in ticks])
cbar.ax.tick_params(labelsize=14)

##ax1.text(-0.9, 0.8,'$(a1)$',color='white', fontsize=15)
ax1.text(-0.1, 0.7,r'$\chi^{0}_{xx}$',color='black', fontsize=16)




#ax2
levels = MaxNLocator(nbins=200).tick_values(Z5axis.min(), Z5axis.max())
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

ax2=fig.add_subplot(324)
c2= ax2.contourf(Yaxis, Zaxis, Z5axis,cmap=cmap,levels=levels, norm=norm)
plt.xlim(-1,1)
plt.ylim(-1,1)
##ax2.set_ylabel(r'$q_{z}/\pi$',FONTSIZE=xy_labelsize)
ax2.set_xticklabels([])
ax2.set_yticklabels([])
##plt.yticks(my_y_ticks,ticks2)


ax2.xaxis.set_major_locator(MultipleLocator(0.5))
ax2.yaxis.set_major_locator(MultipleLocator(0.5))
ax2.tick_params(which='major',direction='in',top='True',right='True',length=8,width=2,labelsize=14)
ax2.spines['top'].set_linewidth(2)
ax2.spines['bottom'].set_linewidth(2)
ax2.spines['left'].set_linewidth(2)
ax2.spines['right'].set_linewidth(2)

ticks=np.linspace(Z5axis.min(),Z5axis.max(),2,endpoint=True)
cbar=plt.colorbar(c2,ax=ax2,ticks=ticks,shrink=0.75)
cbar.ax.set_yticklabels(["{:.2f}".format(i) for i in ticks])
cbar.ax.tick_params(labelsize=14)

ax2.text(-0.1, 0.7,r'$\chi^{0}_{+-}$',color='black', fontsize=16)

#ax3
ax3=fig.add_subplot(326)
levels = MaxNLocator(nbins=200).tick_values(Z6axis.min(), Z6axis.max())
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
c3 = ax3.contourf(Yaxis, Zaxis, Z6axis,cmap=cmap,levels=levels, norm=norm)
plt.xlim(-1,1)
plt.ylim(-1,1)
ax3.set_xlabel(r'$q_{y}/\pi$',FONTSIZE=xy_labelsize)
##ax3.set_ylabel(r'$q_{z}/\pi$',FONTSIZE=xy_labelsize)

ax3.xaxis.set_major_locator(MultipleLocator(0.5))
ax3.yaxis.set_major_locator(MultipleLocator(0.5))
ax3.tick_params(which='major',direction='in',top='True',right='True',length=8,width=2,labelsize=14)
ax3.spines['top'].set_linewidth(2)
ax3.spines['bottom'].set_linewidth(2)
ax3.spines['left'].set_linewidth(2)
ax3.spines['right'].set_linewidth(2)

ax3.set_yticklabels([])
plt.xticks(my_x_ticks)


ticks=np.linspace(Z6axis.min(),Z6axis.max(),2,endpoint=True)
cbar=plt.colorbar(c3,ax=ax3,ticks=ticks,shrink=0.75)
cbar.ax.set_yticklabels(["{:.2f}".format(i) for i in ticks])
cbar.ax.tick_params(labelsize=14)
ax3.text(-0.1, 0.7,r'$\chi^{0}_{-+}$',color='black', fontsize=16)



fig.tight_layout()#调整整体空白
plt.subplots_adjust(wspace =0.2, hspace =0.)#调整子图间距
plt.savefig('gamma%sqx%s_chi0.png'%(gamma,NX),bbox_inches='tight')


plt.show()
plt.close()






