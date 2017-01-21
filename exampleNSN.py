# Here we calculate the heat current from the left normal metal as a function of the bias voltage in NSN system.
# Import all the relevant libraries:
from scipy import *
from matplotlib.pyplot import *
import geometry # Class that describes the geometry (+parameters) of the system
import solver # Class that solves the Usadel equation
reload(solver)

# A function that constructs the system with bias voltage and S wire length given by the arguments
# and following other parameters:
# Temperatures of the reservoirs T_L = T_R = 0.1,
# Resistances of the barriers R_L = 1000, R_R = 100
def generate_geometry(l, v): return geometry.geometry(TL=0.1, TR=0.1, RL = 10000., RR=100., L=l, V=v, sys='NSN')

# Calculate the heat current with different values of bias voltage between 0 and 1.2 and length 0.5
# This takes a while since we are solving all three Usadel equations numerically Nxlen(V) = 200 x 40 = 12000 times. N here comes from the energy discretization.
V = linspace(0., 1.2, 10)
L=.5
#heatCurrent = array([solver.solver(generate_geometry(L, v), N=200,approx_hc=False).heatCurrent() for v in V])

# Lets also do the same but use the short wire approximation and analytic expressions for the integrand in the heat current (this is MUCH faster):
heatCurrentApprox = array([solver.solver(generate_geometry(L, v), bulkappr=True, approx_hc=True).heatCurrent() for v in V])

# Plot the results:
#plot(V,heatCurrent[:],'b',V,heatCurrentApprox[:],'b--')
plot(V,heatCurrentApprox[:],'b--')
xlabel("$e V / \Delta$")
ylabel("$ \dot{Q} e^2 R_\\xi^2 / \Delta^2$")
show()

# Note: there are three approximation schemes built into the solver that can be turned on when constructing solver instance):
# shortappr - assumes the short wire approximation for \theta
# longappr - assumes the long wire approximation with RR=0 for \theta
# approx_hc - uses analytic expressions for distribution function that are calculated assuming the spectral quantities don't vary in space.
