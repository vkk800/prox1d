# This code fits delta and R_L (the cooler resistance) to the data from sample w43. 
# We use the short wire approximation and the approximations for the distribution functions in the fit because with exact numerics it takes a long time and the for this short wire these approximations should be accurate enough.

from scipy import *
from matplotlib.pyplot import *
from scipy.io import loadmat
from scipy.optimize import curve_fit
import geometry
import solver
import coarse

data=loadmat('sinis_AlMn.mat')

Vref = linspace(0, 0.0004, 20)

# We load the raw data from sinis_AlMn.mat to variables ending 'raw'. We then coarsen the voltage grid (coarseArr function) to get values that are as close as possible to ones given in Vref.
w43Tind=20
TRw43raw = data['w43']['T'][0,w43Tind][:,0]
TLw43raw = data['w43']['thermometer'][0,w43Tind][0,:]
Vw43raw = data['w43']['Vc'][0,w43Tind][:,0]
Iw43raw = data['w43']['Ic'][0,w43Tind][:,0]
arr, ind = coarse.coarseArr(Vw43raw, Vref)
TRw43 = TRw43raw[ind] 
TLw43 = TLw43raw[ind]
Vw43 = Vw43raw[ind]
Iw43 = Iw43raw[ind]

# Unit charge and Boltzmann constant
ce = 1.602e-19
kb = 1.38e-23

# System parameters; gap, resistance of a one coherence length long patch of aluminium in normal state, length of the superconductor, area of the contact and normal state resistance (calculated from 4K curves). 
delta0 = .55e-22 # in SI units
Rxi0 = 9.7e-6 # in SI units
L=1.22 # in the units of the coherence length
A=15808 # not really used for anything
RLw43 = 0.345 # not really used for anything

# Takes in voltage grid and some parameters in SI units and returns the calculated charge current in SI units
def IVnum(V, delta=delta0, RL = RLw43, RR=RLw43/20., Rxi=Rxi0):
    RR = abs(RR); RL = abs(RL); delta = abs(delta); Rxi = abs(Rxi)
    # Transform to dimensionless units
    v = ce * V / delta
    tr = TRw43 * kb / delta
    tl = TLw43 * kb / delta
    RR = 0.00646/RL # from using RR value found from low-bias fit
    # Calculate the I-V curve (note: we use the short wire approximation and the approximation for the distribution functions
    current = zeros(len(v))
    for vv in enumerate(v):
        sys = geometry.geometry(tl[vv[0]],tr[vv[0]], RL/Rxi0, RR/Rxi0, L, vv[1],'NSN')
        current[vv[0]] = solver.solver(sys, N=300, shortappr=True, approx_hc=True).chargeCurrent()
    print delta, RL/Rxi
    # And transform back to dimensionful units for the output
    return current*delta/(ce*Rxi)


# Fitting algorithm from scipy.optimize. We fit delta and RL and use for RR value RR = RL/100
paramsw43 = curve_fit(IVnum, Vw43, Iw43, p0=[.5e-22, RLw43])
currentw43 = IVnum(Vw43,delta=paramsw43[0][0],RL=paramsw43[0][1])


# New fit from approximations: RL*RR = 0.00646, RL = 245400, delta = 5.69e-23

# Plot the actual I-V curve and the fit
plot(Vw43, currentw43, 'b', Vw43, Iw43, 'r')
show()
