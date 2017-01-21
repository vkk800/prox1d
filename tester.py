from scipy import *
from matplotlib.pyplot import *
import solver
import geometry

def gensys(v): return geometry.geometry(0.1,0.1,1000.,1.e-10,4.,v,'NSN')

#V=linspace(0.1,2,6)
V=1.

s1 = solver.solver(gensys(V))
s1.solveAll()

#s2 = solver.solver(gensys(V))
#s2.solveAll()


plot(s1.E,s1.fl[:,0,0])
