from scipy import *
from matplotlib.pyplot import *
import solver
import geometry
from scipy.integrate import *
reload(solver)


def gengeom(v): return geometry.geometry(0.15, 0.15, 100, 1e-10, 3.0, v, 'NSN')


x = linspace(0.0, 3., 200)

V = linspace(0., 1.2, 20)
jq = zeros((len(V),))
jq2 = zeros((len(V),))

solution1 = solver.solver(gengeom(0.7), x=x, N=400)
solution1.solveAll()
solution1.solveDelta()
solution1.solveAll()
solution2 = solver.solver(gengeom(0.7), x=x, N=400)
solution2.solveAll()
solution1.heatCurrent()
solution2.heatCurrent()

for v in enumerate(V):
    solution1 = solver.solver(gengeom(v[1]), x=x, N=400)
    solution1.solveAll()
    solution1.solveDelta()
    solution1.solveAll()
    jq[v[0]] = solution1.heatCurrent()
    solution2 = solver.solver(gengeom(v[1]), x=x, N=400)
    solution2.solveAll()
    jq2[v[0]] = solution2.heatCurrent()

plot(V, jq, V, jq2)
show()
