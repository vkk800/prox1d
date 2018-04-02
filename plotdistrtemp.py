from scipy import *
from matplotlib.pyplot import *
# Class that describes the geometry (+parameters) of the system
import geometry
import solver  # Class that solves the Usadel equation
reload(solver)

geom = geometry.geometry(TL=0.03, TR=0.03, RL=1e-6,
                         RR=1e-6, L=10., V=2., sys='NSN')

solution = solver.solver(geom, E=linspace(-2, 2, 200))
solution.solveAll()

plot(solution.E, 1-solution.fl[:, :, 0]-solution.ft[:, :, 0])
show()
