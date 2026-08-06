#! /usr/bin/env python3
###################################################################################################
#
# Code de calcul permettant de resoudre l'equation de convection lineaire 1D
# sur un domaine periodique q_t + q_x = 0 , x \in (0,1) 
#
# checked with Kiwi
#
###################################################################################################

from math import *
from numpy import *
from scipy import *
from pylab import *
from matplotlib import *
import matplotlib.pyplot as plt

def SolInit(L,np,choix):
    x=zeros((np),'d')
    init=zeros((np),'d')

    for i in range(np):
        x[i]=L*float(i)/(np-1)

    if choix=='zero':
        print(choix)
        
    if choix=='gauss':
        print(choix)
        sigma=0.05
        centre=0.3
        for i in range(np):
            init[i]=exp(-(x[i]-centre)*(x[i]-centre)/sigma/sigma)

    if choix=='wiggle':
        print(choix)
        for i in range(np):
            init[i]=pow(-1.,i)
            print(init[i])

    if choix=='packet':
        print(choix)
        sigma=0.05
        periode=0.05
        centre=0.3
        for i in range(np):
            init[i]=sin(2*pi/periode*x[i])*exp(-(x[i]-centre)*(x[i]-centre)/sigma/sigma)    

    if choix=='packet_wig':
        print(choix)
        sigma=0.1
        for i in range(np):
            init[i]=pow(-1.,i)*exp(-(x[i]-0.5)*(x[i]-0.5)/sigma/sigma)    

    if choix=='sine':
        print(choix)
        periode=0.2
        for i in range(np):
            init[i]=sin(2*pi/periode*x[i])
            
    return x,init


def LaxWendroff(cfl,q):
    np=len(q)
    tmp=zeros((np),'d')
    for i in range(1,np-1):
        tmp[i]=-0.5*cfl*(q[i+1]-q[i-1]) + 0.5*cfl*cfl*(q[i+1]-2.*q[i]+q[i-1])

    tmp[0]=tmp[np-1]=-0.5*cfl*(q[1]-q[np-2]) + 0.5*cfl*cfl*(q[1]-2.*q[0]+q[np-2])

    return tmp

def SourceTerm(q,dt,iter,periode,localisation):
    np=len(q)
    time=iter*dt

    for i in localisation:
        q[i]=sin(2*pi/periode*time+pi/2.)

    return q

#####################################################################################

#####################
# Programme principal
#####################

# Parametres - declarations - calculs
L=4.
np=1001
dt=0.0004
tfinal=1.
choix='zero'
legende='oui'
image_base=choix
perio_source=2*pi*dt/(pi*0.1)
i_source=[500]

exact=zeros((np),'d')
x=zeros((np),'d')
q=zeros((np),'d')
q0=zeros((np),'d')
dq=zeros((np),'d')

dx=L/float(np-1)
niter=int(tfinal/dt)
#niter=100
print("niter=",niter)
cfl = dt/dx
print('CFL=',cfl)
    

# Solution initiale et maillage
x,q0=SolInit(L,np,choix)
q=array(q0,copy=True)
for i in range(len(i_source)):
    print("Terme source a x=",x[i_source[i]])


# Boucle iterative
for i in range(1,niter+1):
    if i%100==0:
        print("n=",i)
        
    dq=LaxWendroff(cfl,q)
    q+=dq

    if len(i_source)>0:
        q=SourceTerm(q,dt,i,perio_source,i_source)

str='$t=%2.1f$'%(tfinal)
print(str)
fig=plt.figure()
plt.plot(x,q0,'k-',label=r'$Initial$')
plt.plot(x,q,'r-o',label=str)
plt.axis([0, 4, -1.5e-2, 1.5e-2])
if legende=="oui":
    plt.legend(shadow=True,fancybox=True,loc='upper right')
plt.xlabel(r'$x$',fontsize=20)
plt.ylabel(r'$q(x)$',fontsize=20)
plt.grid(True)

plt.show()

