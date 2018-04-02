from scipy import *
import geometry
from scipy.integrate import quad
from scipy.integrate import trapz

# Polylogarithm


def polylog(n, x, N=40):
    sum = 0.
    for k in range(1, N):
        sum = sum + x**k / k**n
    return sum
# Fermi integral:


def g(n, x, N=40):
    return polylog(n, -exp(-x), N)

# Ideal NIS heat current:


def QNIS(sys, numerical=True):
    if not numerical:
        return 0.25*sqrt(sys.TL*pi/2.)*(4*(sys.V-1.)*g(0.5, (1.-sys.V)/sys.TL)+sys.TL/2 * (3*sys.V-7)*g(1.5, (1-sys.V)/sys.TL))
    return


def QNIS2(sys, numerical=True):
    if not numerical:
        return -sqrt(pi*sys.TL**3)/sqrt(2) * (0.5*g(1.5, (1-sys.V)/sys.TL)+(1-sys.V)/sys.TL * g(0.5, (1-sys.V)/sys.TL))+sqrt(2*pi*sys.TR)*exp(-1/sys.TR)
    return


def INIS(sys, numerical=True):
    if not numerical:
        return sqrt(pi/2)*(sqrt(sys.TL)*exp(-(1.-sys.V)/sys.TL)+sqrt(sys.TR)*exp(-1./sys.TR))
    else:
        def integrand(e): return abs(real(e/sqrt(e**2-1.))) * \
            (1/(exp((e-sys.V)/sys.TL)+1)-1/(exp(e/sys.TR)+1))
    return quad(integrand, -sys.V-2., sys.V+2.)[0]


# Based on the half-infinite approimation:
# Nonequilibrium correction
def dQneq(sys):
    cutoff = (1.*sys.L/(1.*sys.RL))**2
    return -0.5*(sys.L/sys.RL)*exp(-(1-sys.V)/sys.TL)*log(sys.TL/(cutoff*1.78))

# Proxmity effect correction


def dQpe(sys):
    const2 = (2.-sqrt(2.))*exp(-sys.L*sqrt(2.))*(1.+sqrt(2.)+2.*sys.L)/2.
    # -const2*(pi**2 * sys.TL**2 * sys.V**2 / 6 - sys.V**4 / 12)
    return -(tan(pi/8.)*exp(-sys.L*sqrt(2.)))*(2*sys.V**2-2*pi**2 / 3 * (sys.TL**2 - sys.TR**2))


def dIneq(sys):
    return sys.TL*exp(-(1.-sys.V)/sys.TL)*sys.L/sys.RL


def dIpe(sys):
    return 4*tan(pi/8)*exp(-sys.L*sqrt(2))*sys.V


def INIS_accurate(sys):
    Rred = sys.RL*sys.RR/(sys.RL+sys.RR)
    gamma = 1./(2*sys.L*Rred)
    maxE = sys.V+2.
    N = 1000
    de = 3./N
    E = concatenate((linspace(0.000001, 0.95, round(
        N/3)), linspace(0.95+de, 1.05-de, round(N/3)), linspace(1.05, maxE, round(N/3))))
    Ns = real((E+1j*gamma)/sqrt((E+1j*gamma)**2-1.))
    jt = zeros(len(E))
    Dt = 0.5*(1+(E**2+gamma**2+1.)/((E**2+gamma**2)*abs(1-1./(E+1j*gamma)**2)))
    R = -2*1j*(E-1j*gamma-(E+1j*gamma)*sign(1-1./(E+1j*gamma)**2)) / \
        ((E**2+gamma**2)*sqrt(1-1./(E+1j*gamma)**2))
    K = sqrt(R/Dt)
    # jt=Ns*(sys.ftR(E)*Ns*K-sys.ftL(E)*(Ns*K*cosh(sys.L*K)+R*sys.RR*sinh(sys.L*K)))/(Dt*Ns*K*(sys.RL+sys.RR)*cosh(sys.L*K)+(Ns**2+Dt*R*sys.RL*sys.RR)*sinh(sys.L*K))
    jt = Dt*(K*sys.ftL(E)*Ns*(cosh(K*sys.L)*Ns+Dt*K*sys.RR*sinh(K*sys.L))-sys.ftR(E)*Ns**2 * K) / \
        (Dt*K*cosh(K*sys.L)*(sys.RR*Ns+sys.RL*Ns) +
         (Dt**2 * K**2 * sys.RR*sys.RL+Ns**2)*sinh(K*sys.L))

    return trapz(E, jt)


def QNIS_accurate(sys):
    Rred = sys.RL*sys.RR/(sys.RL+sys.RR)
    gamma = 1./(2*sys.L*Rred)
    maxE = sys.V+2.
    N = 400
    de = 3./N
    E = concatenate((linspace(0.000001, 0.95, round(
        N/3)), linspace(0.95+de, 1.05-de, round(N/3)), linspace(1.05, maxE, round(N/3))))
    Ns = real((E+1j*gamma)/sqrt((E+1j*gamma)**2-1.))
    jt = zeros(len(E))
    Dt = 0.5*(1+(E**2+gamma**2+1.)/((E**2+gamma**2)*abs(1-1./(E+1j*gamma)**2)))
    Dl = 0.5*(1+(E**2+gamma**2-1.)/((E**2+gamma**2)*abs(1-1./(E+1j*gamma)**2)))
    R = -1j*(E-1j*gamma-(E+1j*gamma)*sign(1-1./(E+1j*gamma)**2)) / \
        ((E**2+gamma**2)*sqrt(1-1./(E+1j*gamma)**2))
    K = sqrt(R/Dt)
    jl = (sys.flR(E)-sys.flL(E))*Ns/((sys.RR+sys.RL)+Ns*sys.L / Dl)
    jt = Dt*(K*sys.ftL(E)*Ns*(cosh(K*sys.L)*Ns+Dt*K*sys.RR*sinh(K*sys.L))-sys.ftR(E)*Ns**2 * K) / \
        (Dt*K*cosh(K*sys.L)*(sys.RR*Ns+sys.RL*Ns) +
         (Dt**2 * K**2 * sys.RR*sys.RL+Ns**2)*sinh(K*sys.L))
    return -trapz(E, E*jl+sys.VL*jt)


def INISdynes(sys, gamma):
    maxE = sys.V+2.
    N = 400
    de = 3./N
    E = concatenate((linspace(0.000001, 0.95, round(
        N/3)), linspace(0.95+de, 1.05-de, round(N/3)), linspace(1.05, maxE, round(N/3))))
    R = sys.RL+sys.RR
    Ns = real((E+1j*gamma)/sqrt((E+1j*gamma)**2-1.))
    jspectral = Ns*(sys.fL(E)-sys.fR(E))/R
    return trapz(E, jspectral)
