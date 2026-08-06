#! /usr/bin/env python3
###################################################################################################
#
# Calcul de la vitesse de groupe analytique de Lax-Wendroff en fonction de k
#
# checked with Kiwi
#
###################################################################################################

from math import *
import numpy as np
import scipy as sp
from pylab import *
from matplotlib import *
import matplotlib.pyplot as plt

def racines(cfl,omega):
    epsilon=1e-10
    N=omega.size
    tmp1=np.zeros((N),'complex')
    tmp2=np.zeros((N),'complex')
    
    for i in range(N):
        mu=2.*tan(omega[i]/2.*cfl)
        tmp1[i]=-1j*mu+np.sqrt(1-mu*mu)
        tmp2[i]=-1j*mu-np.sqrt(1-mu*mu)
        
    return tmp1,tmp2


def phases(r1,r2):
    epsilon=1e-10
    N=r1.size
    tmp1=np.zeros((N),'d')
    tmp2=np.zeros((N),'d')

    for i in range(N):
        tmp1[i]=atan2(r1[i].imag,(r1[i].real+epsilon))
        tmp2[i]=atan2(r2[i].imag,(r2[i].real+epsilon))

    return tmp1,tmp2


def vitesses_phase(omega,p1,p2):
    epsilon=1e-10
    N=omega.size
    tmp1=np.ones((N),'d')
    tmp2=np.ones((N),'d')

    for i in range(1,N):
        tmp1[i]=pow(-p1[i]/(epsilon+omega[i]),-1)
        tmp2[i]=pow(-p2[i]/(epsilon+omega[i]),-1)

    return tmp1,tmp2   

#####################

toto=0
N = 101
CFL = 0.999
omega=np.zeros((N),'d') # omega est \omega \Delta / c ou c est le CFL
for i in range(N):
    tmp = i/100. *pi
    omega[i] = tmp

un=np.ones((N),'d')
ep1=np.zeros((N),'complex')
ep2=np.zeros((N),'complex')
phi1=np.zeros((N),'d')
phi2=np.zeros((N),'d')
v1=np.zeros((N),'d')
v2=np.zeros((N),'d')

ep1,ep2=racines(CFL,omega)
phi1,phi2=phases(ep1,ep2)
v1,v2=vitesses_phase(omega,phi1,phi2)

plot(omega,v1,'bx')
plot(omega,v2,'gx')
xlim((0,pi))
##ylim((0,1))
show()

if toto==1:
    fig1=plt.figure(figsize=(8,8))

# Affichage
    sp1=fig1.add_subplot(111)
    sp1.set_xticks((0,pi/4.,pi/2.,3.*pi/4.,pi))
    sp1.set_yticks((0,pi/4.,1,pi/2.,2,3.*pi/4.,pi))
    sp1.set_xticklabels((r'$0$',r'$\pi/4$',r'$\pi/2$',r'$3\pi/4$',r'$\pi$'),size=20)
    sp1.set_yticklabels((r'$0$',r'$\pi/4$',r'$1$',r'$\pi/2$',r'$2$',r'$3\pi/4$',r'$\pi$'),size=20)
    sp1.plot(omega,omega,'k-',lw=2,label=r"Exact")
    sp1.legend(shadow=True,fancybox=True,loc='upper left')
    sp1.set_xlabel(r"$\frac{\omega \Delta t}{c}$",fontsize=20)
    sp1.set_ylabel(r"$v_{phi}^{\star} / u",fontsize=20)
    sp1.set_xlim((0,pi))
    sp1.set_ylim((0,pi))
    sp1.grid(True)

    plt.show()

