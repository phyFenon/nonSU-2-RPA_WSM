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


cmap = plt.get_cmap('bwr')
N=64 #64 beta=10

gama=1 #3
NX=16*4 #1 2 3 4(-0.5pi 0 0.5pi pi)
print(NX)

data=np.genfromtxt("gama%s_full_real.txt"%(gama))
x=data[:,0]/np.pi
y=data[:,1]/np.pi
z=data[:,2]/np.pi

z1=np.real(data[:,3])
z2=np.real(data[:,4])
z3=np.real(data[:,5])


X=np.reshape(x,(N,N,N))
Y=np.reshape(y,(N,N,N))
Z=np.reshape(z,(N,N,N))
Z1=np.reshape(z1,(N,N,N))
Z2=np.reshape(z2,(N,N,N))
Z3=np.reshape(z3,(N,N,N))
alpha=0.8


Xaxis=X[:,:,NX-1]
Yaxis=Y[:,:,NX-1]
Z2axis=Z2[:,:,NX-1]
Z3axis=Z3[:,:,NX-1]

if Z2axis.max()>Z3axis.max():
    Zmax=Z2axis.max()
else:
    Zmax=Z3axis.max()
    
print(Zmax)
norm =colors.Normalize(vmin=0,vmax=Zmax, clip=False)
xy_labelsize=15
my_x_ticks=np.arange(-1,1.01,1)
my_y_ticks=np.arange(-1,1.01,1)
my_ticks=[-1,0,1]

fig=plt.figure(figsize=(5.5,3),dpi=150)

ax=fig.add_subplot(121,projection='3d')

ax.plot_surface(Xaxis,Yaxis,Z2axis,cmap=cmap,cstride=1,rstride=1,alpha=alpha,linewidth=0,norm=norm,shade=False)
ax.set_title(r'$\chi^{RPA}_{zz}$',FONTSIZE=16)

##ax.contourf(Xaxis,Yaxis,Z2axis,zdir='y',offset=-1.2,cmap=cmap)
##ax.contourf(Xaxis,Yaxis,Z2axis,zdir='x',offset=-1.2,cmap=cmap)
ax.axis([Xaxis.min(), Xaxis.max(), Yaxis.min(), Yaxis.max()])
ax.set_xlabel(r'$q_{x}/\pi$',FONTSIZE=xy_labelsize)
ax.set_ylabel(r'$q_{y}/\pi$',FONTSIZE=xy_labelsize)
plt.xticks(my_x_ticks,my_ticks)
plt.yticks(my_y_ticks,my_ticks)
plt.tick_params(labelsize=12)
ax.set_xlim(-1.,1.)
ax.set_ylim(-1.,1.)
ax.view_init(10,45)
ax.set_zlim(0,Zmax)

ax=fig.add_subplot(122,projection='3d')

ax.plot_surface(Xaxis,Yaxis,Z3axis,cmap=cmap,cstride=1,rstride=1,alpha=alpha,linewidth=0, norm=norm,shade=False)
ax.set_title(r'$\chi^{RPA}_{+-}$',FONTSIZE=16)
ax.set_xlabel(r'$q_{x}/\pi$',FONTSIZE=xy_labelsize)
ax.set_ylabel(r'$q_{y}/\pi$',FONTSIZE=xy_labelsize)
ax.contourf(Xaxis,Yaxis,Z3axis,zdir='y',offset=-1.2,cmap=cmap)
ax.contourf(Xaxis,Yaxis,Z3axis,zdir='x',offset=-1.2,cmap=cmap)
plt.xticks(my_x_ticks)
plt.yticks(my_y_ticks)
plt.xticks(my_x_ticks,my_ticks)
plt.yticks(my_y_ticks,my_ticks)
plt.tick_params(labelsize=12)
ax.set_zlim(0,Zmax)
ax.set_xlim(-1.2,1.)
ax.set_ylim(-1.2,1.)
ax.view_init(10,45)

fig.tight_layout()#调整整体空白
plt.subplots_adjust(wspace =0.05, hspace =0.)#调整子图间距



plt.savefig('gama%s_qz=%s_3d.png'%(gama,NX))
plt.show()






