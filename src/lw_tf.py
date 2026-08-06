#! /usr/bin/env python3
###################################################################################################
#
# Calcul de la vitesse de groupe analytique de Lax-Wendroff en fonction de k
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

# Improve LaTeX rendering in plots
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern']
plt.rcParams['font.size'] = 11
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['legend.fontsize'] = 11
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 150

def racines(cfl,omega):
    epsilon=1e-10
    np=size(omega)
    tmp1=zeros((np),'complex')
    tmp2=zeros((np),'complex')
    
    for i in range(np):
        b=2.*(1-cfl*cfl-exp(1j*omega[i]*cfl))/(cfl*(cfl-1))
        delta = b*b - 4.*(cfl+1)/(cfl-1)

        tmp1[i]=0.5*(-b+sqrt(delta))
        tmp2[i]=0.5*(-b-sqrt(delta))
        
    return tmp1,tmp2


def phases(r1,r2):
    epsilon=1e-10
    np=size(r1)
    tmp1=zeros((np),'d')
    tmp2=zeros((np),'d')

    for i in range(np):
        tmp1[i]=atan2(r1[i].imag,(r1[i].real+epsilon))
        tmp2[i]=atan2(r2[i].imag,(r2[i].real+epsilon))

    return tmp1,tmp2


def vitesses_phase(omega,p1,p2):
    epsilon=1e-10
    np=size(omega)
    tmp1=ones((np),'d')
    tmp2=zeros((np),'d')

    for i in range(1,np):
        tmp1[i]=pow(-p1[i]/(epsilon+omega[i]),-1)
        tmp2[i]=pow(-p2[i]/(epsilon+omega[i]),-1)

    return tmp1,tmp2


def vitesses_groupe(omega,p1,p2):
    epsilon=1e-10
    np=size(omega)
    tmp1=ones((np),'d')
    tmp2=zeros((np),'d')

    tmp1[0]=-1./((p1[1]-p1[0])/(omega[1]-omega[0]))
    tmp2[0]=-1./((p2[1]-p2[0])/(omega[1]-omega[0]))
    for i in range(1,np-1):
        tmp1[i]=-1./((p1[i+1]-p1[i-1])/(omega[i+1]-omega[i-1]))
        tmp2[i]=-1./((p2[i+1]-p2[i-1])/(omega[i+1]-omega[i-1]))

    tmp1[np-1]=-1./((p1[np-1]-p1[np-2])/(omega[np-1]-omega[np-2]))
    tmp2[np-1]=-1./((p2[np-1]-p2[np-2])/(omega[np-1]-omega[np-2]))

    return tmp1,tmp2    


def longueurs_onde(omega,p1,p2):
    epsilon=1e-10
    np=size(omega)
    tmp1=ones((np),'d')
    tmp2=zeros((np),'d')
    diff=zeros((np),'d')

    tmp1[0]=1000.
    tmp2[0]=2.
    diff[0]=tmp1[0]-tmp2[0]
    
    for i in range(1,np):
        tmp1[i]= 2*pi/(-p1[i]+epsilon)
        tmp2[i]= 2*pi/(-p2[i]+epsilon)
        diff[i]=abs(tmp1[i]-tmp2[i])

    mini=argmin(diff)
        
    return tmp1,tmp2,mini

def error(omega,p1):
    epsilon=1e-10
    np=size(omega)
    p=0

    for i in range(1,np):
        error=abs((-p1[i] - omega[i])/(omega[i]+epsilon))
        print((i,error))
        if error>0.01:
            p=i-1
            break

    return omega[p],2*pi/omega[p]

#####################

def plot_lw_tf(toto='lambda', n_pts=501, CFL=0.7):
    omega = linspace(0, pi, n_pts)

    lambada = zeros(n_pts)
    lambada[0] = 1000.
    lambada[1:] = 2*pi / omega[1:]

    un = ones(n_pts)

    ep1, ep2 = racines(CFL, omega)
    phi1, phi2 = phases(ep1, ep2)
    v1, v2 = vitesses_phase(omega, phi1, phi2)
    vg1, vg2 = vitesses_groupe(omega, phi1, phi2)
    l1, l2, mini = longueurs_onde(omega, phi1, phi2)

    if toto == 'none':
        print(omega[mini], 2*pi/omega[mini], -phi1[mini], -phi2[mini])
        omegaa, Ta = error(omega, phi1)
        print(omegaa, Ta)
        return None, None

    xticks_vals   = (0, pi/4., pi/2., 3.*pi/4., pi)
    xticks_labels = (r'$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$')
    xlabel_str    = r"$\frac{\omega k}{c}$"

    if toto == 'vphi':
        fig, sp = plt.subplots(figsize=(8, 10))
        sp.set_xticks(xticks_vals)
        sp.set_yticks((0, 0.25, 0.5, 0.75, 1, 1.25))
        sp.set_xticklabels(xticks_labels, size=15)
        sp.set_yticklabels((r'$0$', r'$0.25$', r'$0.5$', r'$0.75$', r'$1$', r'$1.25$'), size=15)
        sp.plot(omega, un, 'k-', lw=2, label=r"Exact")
        sp.plot(omega[:mini], v1[:mini], 'g-', lw=1, label=r"p-wave")
        sp.plot(omega[:mini], v2[:mini], 'r-', lw=1, label=r"q-wave")
        sp.plot(omega[mini+1:], v1[mini+1:], 'g--', lw=0.5)
        sp.plot(omega[mini+1:], v2[mini+1:], 'r--', lw=0.5)
        sp.legend(shadow=True, fancybox=True, loc='lower right')
        sp.set_xlabel(xlabel_str, fontsize=24)
        sp.set_ylabel(r"$v_{\varphi}^{\star} / u$", fontsize=20)
        sp.set_xlim((0, pi))
        sp.set_ylim((0, 1.25))
        sp.grid(True)

    elif toto == 'lambda':
        fig, sp = plt.subplots(figsize=(10, 10))
        sp.set_xticks(xticks_vals)
        sp.set_yticks((0, 5, 10, 15, 20, 25, 30))
        sp.set_xticklabels(xticks_labels, size=15)
        sp.set_yticklabels((r'$0$', r'$5$', r'$10$', r'$15$', r'$20$', r'$25$', r'$30$'), size=15)
        sp.plot(omega, lambada, 'k-', lw=2, label=r"Exact")
        sp.plot(omega[:mini], l1[:mini], 'g-', lw=1, label=r"p-wave")
        sp.plot(omega[:mini], l2[:mini], 'r-', lw=1, label=r"q-wave")
        sp.plot(omega[mini+1:], l1[mini+1:], 'g--', lw=0.5)
        sp.plot(omega[mini+1:], l2[mini+1:], 'r--', lw=0.5)
        sp.legend(shadow=True, fancybox=True, loc='upper right')
        sp.set_xlabel(xlabel_str, fontsize=24)
        sp.set_ylabel(r"$\lambda^{\star}/h$", fontsize=20)
        sp.set_xlim((0, pi))
        sp.set_ylim((0, 30))
        sp.grid(True)

    elif toto == 'vg':
        fig, sp = plt.subplots(figsize=(10, 10))
        sp.set_xticks(xticks_vals)
        sp.set_xticklabels(xticks_labels, size=15)
        sp.plot(omega, un, 'k-', lw=2, label=r"Exact")
        sp.plot(omega[:mini], vg1[:mini], 'g-', lw=1, label=r"p-wave")
        sp.plot(omega[:mini], vg2[:mini], 'r-', lw=1, label=r"q-wave")
        sp.plot(omega[mini+1:], vg1[mini+1:], 'g--', lw=0.5)
        sp.plot(omega[mini+1:], vg2[mini+1:], 'r--', lw=0.5)
        sp.legend(shadow=True, fancybox=True, loc='upper right')
        sp.set_xlabel(xlabel_str, fontsize=24)
        sp.set_ylabel(r"$\mathcal{V}_g^{\star}/ u$", fontsize=20)
        sp.set_xlim((0, pi))
        sp.set_ylim((-10, 10))
        sp.grid(True)

    elif toto == 'abs':
        fig, sp = plt.subplots(figsize=(10, 10))
        sp.set_xticks(xticks_vals)
        sp.set_xticklabels(xticks_labels, size=15)
        sp.plot(omega, un, 'k-', lw=2, label=r"Exact")
        sp.plot(omega[:mini], abs(ep1[:mini]), 'g-', lw=1, label=r"p-wave")
        sp.plot(omega[:mini], abs(ep2[:mini]), 'r-', lw=1, label=r"q-wave")
        sp.plot(omega[mini+1:], abs(ep1[mini+1:]), 'g--', lw=0.5)
        sp.plot(omega[mini+1:], abs(ep2[mini+1:]), 'r--', lw=0.5)
        sp.legend(shadow=True, fancybox=True, loc='upper right')
        sp.set_xlabel(xlabel_str, fontsize=24)
        sp.set_ylabel(r"$\left| E \right|$", fontsize=24)
        sp.set_xlim((0, pi))
        sp.set_ylim((0, 5))
        sp.grid(True)

    elif toto == 'omega':
        fig, sp = plt.subplots(figsize=(10, 10))
        sp.set_xticks(xticks_vals)
        sp.set_xticklabels(xticks_labels, size=15)
        sp.plot(omega, omega, 'k-', lw=2, label=r"Exact")
        sp.plot(omega[:mini], -phi1[:mini], 'g-', lw=1, label=r"p-wave")
        sp.plot(omega[:mini], -phi2[:mini], 'r-', lw=1, label=r"q-wave")
        sp.plot(omega[mini+1:], -phi1[mini+1:], 'g--', lw=0.5)
        sp.plot(omega[mini+1:], -phi2[mini+1:], 'r--', lw=0.5)
        sp.legend(shadow=True, fancybox=True, loc='lower right')
        sp.set_xlabel(xlabel_str, fontsize=24)
        sp.set_ylabel(r"$\frac{\omega^{\star} k}{c}$", fontsize=24)
        sp.set_xlim((0, pi))
        sp.set_ylim((0, pi))
        sp.grid(True)

    return fig, sp


if __name__ == "__main__":
    # choose plot: 'vphi', 'lambda', 'vg', 'abs', 'omega', or 'none' (print only)
    toto  = 'omega'
    n_pts = 1001    # number of pulsation samples
    CFL   = 0.1    # CFL number

    fig, ax = plot_lw_tf(toto=toto, n_pts=n_pts, CFL=CFL)
    if fig is not None:
        fig.savefig('doc/Images/lw_tf_{}.png'.format(toto), bbox_inches='tight')
        plt.show()
