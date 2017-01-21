from scipy import *
from geometry import *
import solver
import approximations
import scikits.bvp1lg as bvp
import scipy.interpolate as interp
from scipy.optimize import fmin
from scipy.optimize import brentq
import copy

def optimV(sys, Vmin, Vmax):
    V = linspace(Vmin, Vmax, 200)
    E = linspace(0.1,2,50)
    x = linspace(0, sys.L, 5)
    best = -1.
    bestV = -1
    current = []
    for v in V:
        sys.V = v
        solv = solver(sys, E, x)
        current.append(solv.heatCurrent())
        if current[-1] > best:
            best = current[-1]
            bestV = v
        print v, current[-1]

    print (bestV, max(current))
    return (bestV, max(current))


def optimV2(sys, approx_theta=False, approx_hc=True):
    E = linspace(0.1,2,50)
    x = linspace(0,sys.L,5)
    solv=solver(sys, E, x, approx_theta=approx_theta, approx_hc=approx_hc)
    def minfun(v):
        sys2 = sys
        sys2.V = v
        return -solver(sys2, E, x, approx_theta=approx_theta, approx_hc=approx_hc).heatCurrent()
    sol = fmin(minfun,0.85,xtol=1e-3)
    sys.V = sol[0]
    fsol = solver(sys,E,x, approx_theta=approx_theta, approx_hc=approx_hc).heatCurrent()
    print sys.RL, sys.L, sol, fsol
    return array((sol, fsol), dtype=float)

def optimL2(sys, approx_theta=False, approx_hc = False):
    E=linspace(0.01,3,150)
    x=linspace(0,sys.L,5)
    Ndef=350
    solv = solver.solver(sys, approx_theta=approx_theta, approx_hc=approx_hc,N=Ndef)
    def minfun(l):
        sys2 = sys
        sys2.L = l
        return -solver.solver(sys2, approx_theta=approx_theta, approx_hc=approx_hc,N=Ndef).heatCurrent()
    sol = fmin(minfun,4.5,xtol=1e-2, disp=False)
    sys.L = sol[0]
    fsol = solver.solver(sys, approx_theta=approx_theta, approx_hc=approx_hc,N=Ndef).heatCurrent()
    #print sys.RL, sys.V, sol, fsol
    return array((sol, fsol), dtype=float)

def optimLV(sys):
    E=linspace(0.1,2,50)
    def minfun(l):
        x=linspace(0,l,5)
        sys.L = l
        sys.V = optimV2(sys)[0]
        return -solver(sys, E, x).heatCurrent()

    sol = fmin(minfun,4,xtol=1e-2)
    sys.L = sol[0]
    x=linspace(0,sys.L,5)
    fsol = solver(sys,E,x).heatCurrent()
    print sys.RL, sys.V, sol, fsol
    return array((sys.RL,sys.V, sys.L, fsol), dtype=float)

def minT(sys, approx_theta=False, approx_hc = False,optV=False, short_approx=False):
    def fun(tn):
        V = sys.V
        if optV: V = 1-0.66*tn
        syss = geometry(tn, sys.TR, sys.RL, sys.RR, sys.L, V, sys.sys)
        hc = solver.solver(syss, N=300, approx_theta=approx_theta, approx_hc=approx_hc, short_approx=short_approx).heatCurrent()
        print syss.V,syss.TL,syss.TR, hc
        return hc
    
    return brentq(fun, 0.001, 0.5)


def minTshort(sys, optV=False, approx=False):
    def fun(tn):
        sys2 = sys
        sys2.TL = tn
        if optV: sys2.V = 1-0.66*tn
        return approximations.QNIS_accurate(sys2)
    
    return brentq(fun, 0.001, 0.5)

def optimRR(sys, optV=False):
    def minfun(rr):
        sys.RR = abs(rr)
        if abs(rr) > sys.RL: sys.RR = sys.RL
        if abs(rr) < 1e-10: sys. RR = 1e-10
        try:
            mint = minT(sys, optV)
        except ValueError:
            mint = 0.5
        print sys.RR, sys.RL, mint
        return mint
    sol = fmin(minfun,1.e-2*sys.RL,xtol=1., disp=False)
    sys.RR=sol[0]
    fsol = minT(sys,optV)
    return array((sol, fsol), dtype=float)


def optimL(sys, approx_hc=False, approx_theta=False, short_approx=False, optV=False):
    def minfun(l):
        sys.L = abs(l)
        try: mint = minT(sys, approx_theta=approx_theta, approx_hc=approx_hc, short_approx=short_approx, optV=optV)
        except ValueError: mint = .5
        
    sol = fmin(minfun, 1., xtol=.05, disp=False)
    sys.L=sol[0]
    fsol = minT(sys, approx_theta=approx_theta, approx_hc=approx_hc, short_approx=short_approx, optV=optV)
    return array((sol, fsol),dtype=float)

def optimRR(sys, approx_hc=False, approx_theta=False, short_approx=False, optV=False):
    def minfun(r):
        sys.RR = abs(r)
        try: mint = minT(sys, approx_theta=approx_theta, approx_hc=approx_hc, short_approx=short_approx, optV=optV)
        except ValueError: mint = .5
        
    sol = fmin(minfun, 100., xtol=10., disp=False)
    sys.RR=sol[0]
    fsol = minT(sys, approx_theta=approx_theta, approx_hc=approx_hc, short_approx=short_approx, optV=optV)
    return array((sol, fsol),dtype=float)
