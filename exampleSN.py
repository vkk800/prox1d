# Here we calculate the heat current from the left normal metal as a function of the bias voltage in NSN system.
# Import all the relevant libraries:
from scipy import *
from matplotlib.pyplot import *
# Class that describes the geometry (+parameters) of the system
import geometry
import solver  # Class that solves the Usadel equation
reload(solver)

# Define the geometry. It's really an NSN system, but if we set the left resistance to ~infinity and keep in equilibrium, it should be the same as just SN system.

# The variables are: TL, TR = temperatures of the left and right leads (assumed to be in equilibrium), RL, RR = resistances of the left and right interfaces, L = length of the superconductor, V=bias voltage, sys = describes the system components, doesn't work very reliably for anything else than NSN, but in principle could be also SNS, SSS or whatever.

geom1 = geometry.geometry(TL=.1, TR=.1, RL=10000.,
                          RR=0., L=5., V=0., sys='NSN')
geom2 = geometry.geometry(TL=.1, TR=.1, RL=10000.,
                          RR=1., L=5., V=0., sys='NSN')

# Create solver object and solve for theta
solv = solver.solver(geom1)
solv.solveTheta()
solv2 = solver.solver(geom2)
solv2.solveTheta()


# Plot the density of states at different points. The interface is at x=5.
# For zero interface resistance
subplot(2, 1, 1)
plot(solv.E, real(cosh(solv.theta[:, :, 0])))
title("$R = 0$")
xlabel("$E / \Delta$")
ylabel("$DOS$")
labels = []
for xx in solv.x:
    labels.append("$x/\\xi = "+str(xx)+"$")
legend(labels)

# For resistance = 1
subplot(2, 1, 2)
plot(solv2.E, real(cosh(solv2.theta[:, :, 0])))
title("$R = 1$")
xlabel("$E / \Delta$")
ylabel("$DOS$")

show()
