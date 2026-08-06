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
from matplotlib.patches import Rectangle

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

def kdx_lw(cfl,k):
    epsilon=1e-10
    np=size(k)
    kstar=zeros((np),'d')

    for i in range(np):
        num = cfl*sin(k[i])
        denom = 1 - cfl*cfl*(1-cos(k[i]))
        kstar[i]=(1/(cfl+epsilon))*atan2(num,denom)

    return kstar


def vg_lw(cfl,k):
    epsilon=1e-10
    np=size(k)
    vg=zeros((np),'d')
    ##vg2=zeros((np),'d')

    for i in range(np):
        num1 = cos(k[i])*(1-cfl*cfl*(1.-cos(k[i]))) + sin(k[i])*(cfl*cfl*sin(k[i]))

        denom1 = (1-cfl*cfl*(1.-cos(k[i]))) * (1-cfl*cfl*(1.-cos(k[i])))

        num2 = num1 / denom1

        denom2 = 1. + cfl*cfl*sin(k[i])*sin(k[i])/denom1

        vg[i] = num2 / denom2

        ##vg2[i]=((1.-cfl*cfl)*cos(k[i])+cfl*cfl)/(cfl*cfl*sin(k[i])*sin(k[i])+(cfl*cfl*cos(k[i])-cfl*cfl+1)*(cfl*cfl*cos(k[i])-cfl*cfl+1))

        ##print vg[i]-vg2[i]

    return vg

def modg_lw(cfl,k):
    epsilon=1e-10
    np=size(k)
    modg=zeros((np),'d')

    for i in range(np):
        g=1-1j*cfl*sin(k[i])-cfl*cfl*(1-cos(k[i]))
        modg[i]=abs(g)
        
    return modg    

#####################

def plot_lw(toto="ampli", n_pts=101, cfl0=0.1, dcfl=0.2, ncfl=6):
    k = linspace(0, pi, n_pts)
    un = ones(n_pts)

    cfls = cfl0 + arange(ncfl) * dcfl
    liste = []
    for c in cfls:
        print("cfl=", c)
        if toto == "kdx":
            liste.append(kdx_lw(c, k))
        elif toto == "vg":
            liste.append(vg_lw(c, k))
        elif toto == "ampli":
            liste.append(modg_lw(c, k))

    if toto == "kdx":
        fig, sp = plt.subplots(figsize=(8, 8))
        sp.set_xticks((0, pi/4., pi/2., 3.*pi/4., pi))
        sp.set_yticks((0, pi/4., 1, pi/2., 2, 3.*pi/4., pi))
        sp.set_xticklabels((r'$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$'), size=20)
        sp.set_yticklabels((r'$0$', r'$\pi/4$', r'$1$', r'$\pi/2$', r'$2$', r'$3\pi/4$', r'$\pi$'), size=20)
        sp.plot(k, k, 'k-', lw=2, label=r"Exact")
        for c, data in zip(cfls, liste):
            sp.plot(k, data, lw=1, label="%3.1f" % c)
        sp.legend(shadow=True, fancybox=True, loc='upper left')
        sp.set_xlabel(r"$\xi h$", fontsize=20)
        sp.set_ylabel(r"$\xi^{\star} h$", fontsize=20)
        sp.set_xlim((0, pi))
        sp.set_ylim((0, pi))
        sp.grid(True)

    elif toto == "vg":
        fig, sp = plt.subplots(figsize=(12, 7))
        sp.set_xticks((0, pi/4., pi/2., 3.*pi/4., pi))
        sp.set_yticks((-2, -1, 0, 1))
        sp.set_xticklabels((r'$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$'), size=20)
        sp.set_yticklabels((r'$-2$', r'$-1$', r'$0$', r'$1$'), size=20)
        sp.plot(k, un, 'k-', lw=2, label=r"Exact")
        for c, data in zip(cfls, liste):
            sp.plot(k, data, lw=1, label="%3.1f" % c)
        sp.legend(shadow=True, fancybox=True, loc='lower left')
        sp.set_xlabel(r"$\xi h$", fontsize=20)
        sp.set_ylabel(r"$\mathcal{V}^{\star}_{g} / u$", fontsize=20)
        sp.set_xlim((0, pi))
        sp.set_ylim((-2, 1.5))
        sp.add_patch(Rectangle((0, -2), pi, 2, facecolor="#aaaaaa"))
        sp.grid(True)

    elif toto == "ampli":
        fig, sp = plt.subplots(figsize=(10, 7))
        sp.set_xticks((0, pi/4., pi/2., 3.*pi/4., pi))
        sp.set_yticks((0, 0.5, 1))
        sp.set_xticklabels((r'$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$'), size=20)
        sp.set_yticklabels((r'$0$', r'$0.5$', r'$1$'), size=20)
        sp.plot(k, un, 'k-', lw=2, label=r"Exact")
        for c, data in zip(cfls, liste):
            sp.plot(k, data, lw=1, label="%3.1f" % c)
        sp.legend(shadow=True, fancybox=True, loc='lower left')
        sp.set_xlabel(r"$\xi h$", fontsize=20)
        sp.set_ylabel(r"$\left|g(\xi) \right|$", fontsize=20)
        sp.set_xlim((0, pi))
        sp.set_ylim((0, 1.2))
        sp.add_patch(Rectangle((0, 1), pi, 0.2, facecolor="#aaaaaa"))
        sp.grid(True)

    return fig, sp


if __name__ == "__main__":
    # choose plot: "kdx" (modified wavenumber), "vg" (group velocity), "ampli" (amplification factor)
    toto  = "ampli"
    n_pts = 501    # number of wavenumber samples
    cfl0  = 0.1    # first CFL number
    dcfl  = 0.2    # CFL step
    ncfl  = 6      # number of CFL values

    fig, ax = plot_lw(toto=toto, n_pts=n_pts, cfl0=cfl0, dcfl=dcfl, ncfl=ncfl)
    
    # Display or save the figure
    fig.savefig('../doc/Images/{}_lw.png'.format(toto), dpi=150, bbox_inches='tight')
    # plt.show()



