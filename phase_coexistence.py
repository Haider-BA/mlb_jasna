#
#  Maxwell construction for phase coexistence
#
#  Copyright (c) 2015 Ulf D. Schiller <mail@ulfschiller.de>
#  All rights reserved.


import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plot, ticker
from scipy import integrate, optimize, misc


titlefontsize=24
labelfontsize=20
tickfontsize=12
textfontsize=16
linewidth=2

mpl.rc('font', size=titlefontsize)
mpl.rc('axes', labelsize=labelfontsize)
mpl.rc(('xtick','ytick'), labelsize=tickfontsize)


R1 = 0.5
R2 = 1.0
R3 = 2.0*R2-R1
S1 = 0.6
S2 = 0.75
W  = 0.25*(R2-R1)

Tmin = 0.5
Tmax = 1.0

npoints   = 100
delta_rho = (R3-R1)/npoints
epsilon   = 1.e-2


def cs2(x, S2, S1=S1, W=W):
    cs2 = (S2-S1)*np.exp(-(x-R2)*(x-R2)/(2.0*W*W))+S1
    return cs2


def pressure(x, S2, S1=S1, W=W):
    p = x*cs2(x, S2, S1, W)
    return p


def psi(x, S2, S1=S1, W=W):    
    func = lambda x: np.log(cs2(x, S2, S1, W))
    psi = misc.derivative(func, x, delta_rho)
    return psi


def derP(x, S2, S1=S1, W=W):
    func = lambda x: pressure(x, S2, S1, W)
    dP   = misc.derivative(func, x, delta_rho)
    return dP


def der2P(x, S2, S1=S1, W=W):
    func = lambda x: derP(x, S2, S1, W)
    d2P  = misc.derivative(func, x, delta_rho)
    return d2P


def free_energy(rhomin, rhomax, S2, S1=S1, W=W):
    func = lambda x: pressure(1./x, S2, S1, W)
    res = integrate.quad(func, 1./rhomin, 1./rhomax)
    energy = res[0]
    return energy


def critical_density(S2, S1=S1, W=W):
    func = lambda x: psi(x, S2, S1, W)
    res = optimize.minimize_scalar(func) # may not be exact
    rhocrit = res.x
    func = lambda x: 2.*derP(x, S2, S1, W) + x*der2P(x, S2, S1, W)
    res = optimize.fsolve(func, rhocrit)
    rhocrit = res
    return rhocrit


def stability_limits(S2, S1=S1, W=W):
    func = lambda x: psi(x, S2, S1, W) + 1./x
    rhocrit = critical_density(S2, S1, W)
    if (psi(rhocrit, S2, S1, W) < -1./rhocrit):
        rhomin = optimize.bisect(func, R2, rhocrit)
        rhomax = optimize.bisect(func, rhocrit, 10.)
    else:
        rhomin = rhocrit
        rhomax = rhocrit
    return rhomin, rhomax


def densities(p, S2, S1=S1, W=W):
    func = lambda x: pressure(x, S2, S1, W) - p
    rhomin, rhomax = stability_limits(S2, S1, W)
    rho1 = optimize.fsolve(func, rhomin-delta_rho)[0]
    rho2 = optimize.fsolve(func, 0.5*(rhomin+rhomax))[0]
    rho3 = optimize.fsolve(func, rhomax+delta_rho)[0]
    return rho1, rho2, rho3


def maxwell_construction(p, S2, S1=S1, W=W):
    rho1, rho2, rho3 = np.sort(densities(p, S2, S1, W))
    area = free_energy(rho1,rho3,S2,S1,W) - p * (1./rho3 - 1./rho1)
    return area


if __name__ == '__main__':

    fig = plot.figure(figsize=(8.,8.))
    ax1 = plot.subplot(224)
    ax2 = plot.subplot(223)
    ax3 = plot.subplot(221)
    ax4 = plot.subplot(222)

    rs = np.arange(R1,R3,delta_rho)
    rg = []
    rl = []
    ts = []

    for T in np.arange(Tmin,Tmax,(Tmax-Tmin)/npoints):

        rhomin, rhomax = stability_limits(S2=T)

        pmin  = pressure(rhomax, S2=T)
        pmax  = pressure(rhomin, S2=T)
    
        a1 = maxwell_construction(pmin, S2=T)
        a2 = maxwell_construction(pmax, S2=T)

        if (a1*a2 < 0.):

            pv = optimize.bisect(maxwell_construction, pmin, pmax, (T,S1,W))

            rhos = np.array(densities(pv,S2=T))
            vps  = [ pressure(x, S2=T) for x in rhos ]
            
            rhogas    = rhos[0]
            rholiquid = rhos[2]

            rg.append(rhogas)
            rl.append(rholiquid)
            ts.append(T)

    T = S2
    rhocrit = critical_density(S2=T)
    rhomin, rhomax = stability_limits(S2=T)
    pcrit = pressure(rhocrit, S2=T)
    pmin  = pressure(rhomin, S2=T)
    pmax  = pressure(rhomax, S2=T)
    a1 = maxwell_construction(pmin, S2=T)
    a2 = maxwell_construction(pmax, S2=T)
    if (a1*a2 < 0.):
        pv = optimize.bisect(maxwell_construction, pmin, pmax, (T,S1,W))
    else:
        pv = pcrit
    rhos  = np.array(densities(pv, S2=T))
    vps  = [ pressure(x, S2=T) for x in rhos ]
    cs    = cs2(rs, S2=T)
    ps    = pressure(rs, S2=T)
    psis  = psi(rs, S2=T)
    #es = np.array( [ free_energy(R1, rho, S2=T) for rho in rs ] )

    ax1.plot(rg,ts,".b")
    ax1.plot(rl,ts,".b")

    ax1.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax1.set_xlabel(r'$\rho$')
    ax1.set_ylabel(r'$A$')

    ax2.plot(1/rs, ps, lw=linewidth)
    ax2.plot(1/rhos, vps,"o")
    ax2.plot(1/rhos, np.zeros_like(rhos)+pv)
    #ax2.plot(1/rhomin, pmin, "o")
    #ax2.plot(1/rhomax, pmax, "o")
    #ax2.plot(1/rhocrit, pcrit, "o")
    
    ax2.set_xlabel(r'$1/\rho$')
    ax2.set_ylabel(r'$p(\rho)$')

    ax3.plot(rs, cs, lw=linewidth)
    
    ax3.set_xlabel(r'$\rho$')
    ax3.set_ylabel(r'$c_s^2(\rho)$')

    ax4.plot(rs, psis, lw=linewidth)
    ax4.plot(rs, 1/rs, "k:", lw=linewidth)
    ax4.plot(rs, -1/rs, "k:", lw=linewidth)
    ax4.plot(rhocrit, psi(rhocrit, S2=T), "o")
    #ax4.plot(rhomin, psi(rhomin, S2=T), "o")
    #ax4.plot(rhomax, psi(rhomax, S2=T), "o")
    ax4.set_xlabel(r'$\rho$')
    ax4.set_ylabel(r'$\psi(\rho)$')

    fig.tight_layout()

    fig.savefig('phase_coexistence.png')

    data = zip(ts,rg,rl)

    np.savetxt('phasdiag_data.txt', data)

    plot.show()
