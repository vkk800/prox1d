from scipy import *
from matplotlib.pyplot import *
import optims
import geometry
import pickle
RL = 100000


def gensys(rl, tr): return geometry.geometry(
    0.1, tr, rl, 1e-10, 13, 0.9, 'NSN')


RL = array([1000, 5000, 10000, 50000])
RR = array([1e-10])

TR = linspace(0.01, 0.4, 10)
COMP = True

if COMP:
    minT = zeros((len(RL), len(RR), len(TR)))
    minTa = zeros((len(RL), len(RR), len(TR)))
    for rl in enumerate(RL):
        for rr in enumerate(RR):
            for tr in enumerate(TR):
                sys = geometry.geometry(
                    tr[1], tr[1], rl[1], rr[1], 13., 0.4, 'NSN')
                minT[rl[0], rr[0], tr[0]] = optims.minT(
                    sys, approx_hc=False, approx_theta=False, optV=True)
                minTa[rl[0], rr[0], tr[0]] = optims.minT(
                    sys, approx_hc=True, approx_theta=False, optV=True)
                print minT
    out = open('minTplot.pkl', 'wb')
    pickle.dump(minT, out)
    out.close()

inp = open('minTplot.pkl', 'rb')
minT = pickle.load(inp)
inp.close()


plot(TR, minT[0, 0, :], 'g-', TR, minT[1, 0, :], 'b-.', TR, minT[2,
                                                                 0, :], 'r--', TR, minT[3, 0, :], 'c-o', TR, TR, 'k--', markersize=3)
xlabel("$k_B T_R / \Delta_0$", fontsize=18)
ylabel("$k_B T_L / \Delta_0$", fontsize=18)
legend(["$R_L / R_\\xi = 1000$", "$R_L / R_\\xi = 5000 $",
        "$R_L / R_\\xi = 10000$", "$R_L / R_\\xi = 50000$"], loc=2)
axis([0.01, 0.31, 0, 0.31])
gcf().subplots_adjust(left=0.13, bottom=0.13)
savefig("minTL2.png")
show()
