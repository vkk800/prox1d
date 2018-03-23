from scipy import *
from geometry import *
import scikits.bvp1lg as bvp
import scipy.interpolate as interp
import scipy.integrate as integrate


class solver:
    def __init__(self, geometry, E=None, x=None, N=500, shortappr=False,
                 longappr=False, approx_hc=False, approx_theta=False,
                 short_approx=False, bulkappr=False):
        self.geometry = geometry
        self.E = E
        self.x = x

        if self.E is None:
            maxE = self.geometry.V + 2.
            de = 3./N
            self.E = concatenate((linspace(0.000001, 0.95, round(N/3)),
                                  linspace(0.95+de, 1.05-de, round(N/3)),
                                  linspace(1.05, maxE, round(N/3))))
        if self.x is None:
            self.x = linspace(0, self.geometry.L, 11)

        self.delta = geometry.Delta0
        self.shortappr = shortappr
        self.shortappr = short_approx
        self.approx_hc = approx_hc
        self.longappr = longappr
        self.longappr = approx_theta
        self.bulkappr = bulkappr

        self.theta = zeros((self.E.size, self.x.size, 2), dtype=complex)
        self.theta_func = []
        self.fl = zeros((self.E.size, self.x.size, 2))
        self.fl_func = []
        self.ft = zeros((self.E.size, self.x.size, 2))
        self.ft_func = []

        self.solved = {'theta': False, 'fl': False, 'ft': False}

    # Let's first solve the pairing angle

    # The solver for theta:
    def solveTheta(self):
        # Short wire approximation
        if self.shortappr:
            for e in enumerate(self.E):
                def thetafun(x):
                    x = asarray(x)
                    Rpar = ones(len(x))*self.geometry.RL * \
                        self.geometry.RR/(self.geometry.RL+self.geometry.RR)
                    return array([arctanh(self.delta(x)/(e[1]+1.*1j/(2*self.geometry.L*Rpar))),
                                  0.*ones(len(x))]).transpose()
                self.theta_func.append(thetafun)
                self.theta[e[0], :, :] = thetafun(self.x)
            self.solved['theta'] = True
            return

        if self.bulkappr:
            for e in enumerate(self.E):
                def thetafun(x):
                    x = asarray(x)
                    return array([ones(len(x))*arctanh(self.delta(x)/(e[1])),
                                  0.*ones(len(x))]).transpose()
                self.theta_func.append(thetafun)
                self.theta[e[0], :, :] = thetafun(self.x)
            self.solved['theta'] = True
            # Long wire approximation (with RR = 0)
        if self.longappr:
            for e in enumerate(self.E):
                def thetafun(x):
                    x = asarray(x)
                    ts = arctanh(self.delta(x)/e[1])
                    alpha = sqrt(self.delta(x)**2-e[1]**2)
                    return array([ts-4*arctanh(exp((x-self.geometry.L)*sqrt(2*alpha))*tanh(ts/4)),
                                  -2.*sqrt(2*alpha)*sinh(ts/2)/(cosh(sqrt(2*alpha)*(x-self.geometry.L))-cosh(ts/2)*sinh(sqrt(2*alpha)*(x-self.geometry.L)))]).transpose()
                self.theta_func.append(thetafun)
                self.theta[e[0], :, :] = thetafun(self.x)
            self.solved['theta'] = True
            return
        # # In some cases we can use the analytical half-infinite approximation for theta. Valid when R_L >> 1, R_R = 0 and L >~ 4.
        # if self.approx_theta:
        #     for e in enumerate(self.E):
        #         def thetafun(x): x = asarray(x); return array([arctanh(self.delta(x)/e[1])-4*arctanh(exp((x-self.geometry.L)*sqrt(2*sqrt(self.delta(x)**2-e[1]**2)))*tanh(arctanh(self.delta(x)/e[1])/4)), -2*sqrt(2)*sqrt(sqrt(self.delta(x)**2-e[1]**2))*sinh(arctanh(self.delta(x)/e[1])/2)/(cosh(sqrt(2)*(x-self.geometry.L)*sqrt(sqrt(self.delta(x)**2-e[1]**2)))-cosh(arctanh(self.delta(x)/e[1])/2)*sinh(sqrt(2)*(x-self.geometry.L)*sqrt(sqrt(self.delta(x)**2-e[1]**2))))]).transpose()

        #         self.theta_func.append(thetafun)
        #         self.theta[e[0],:,:] = thetafun(self.x)
        #     self.solved['theta'] = True
        #     return

        # Here's the exact numerical solution:
        boundaries = [0, self.geometry.L]
        tol = [1e-6, 1e-6]
        self.E = self.E[::-1]
        for e in enumerate(self.E):
            def fsub(x, z):
                t, dt = z
                return array([-2*1j*e[1]*sinh(t)+2*1j*self.delta(x)*cosh(t)])

            def dfsub(x, z):
                t, dt = z
                return array([[-2*1j*e[1]*cosh(t)+2*1j*self.delta(x)*sinh(t)]])

            def gsub(z):
                t, dt = z
                return array([self.geometry.RL*dt[0]+sinh(self.geometry.thetaL(e[1])-t[0]),
                              self.geometry.RR*dt[1]+sinh(t[1]-self.geometry.thetaR(e[1]))])

            def dgsub(z):
                t, dt = z
                return array([[-cosh(self.geometry.thetaL(e[1])-t[0]), self.geometry.RL],
                              [cosh(t[1]-self.geometry.thetaR(e[1])), self.geometry.RR]])

            def guess_fun(x):
                #z=array([[real(arctanh(1/e[1])*(self.geometry.L-xx)/self.geometry.L+1),real(arctanh(1/e[1])/self.geometry.L),imag(arctanh(1/e[1])*(self.geometry.L-xx)/self.geometry.L), imag(arctanh(1/e[1])/self.geometry.L)] for xx in x]).transpose()
                #dm = zeros((2,) + asarray(x).shape)
                z = zeros((4,) + asarray(x).shape)
                dm = zeros((2,) + asarray(x).shape)
                z[0] = 1

                return z, dm

            if e[0] > 0:
                def guess(x):
                    dm = zeros((2,) + asarray(x).shape)
                    return array([real(self.theta_func[e[0]-1](x)[:, 0]), real(self.theta_func[e[0]-1](x)[:, 1]),
                                  imag(self.theta_func[e[0]-1](x)[:, 0]), imag(self.theta_func[e[0]-1](x)[:, 1])]), dm
            else:
                guess = guess_fun

            # ncomp = 1, mstar = 2

            solution = bvp.colnew.solve(boundaries, [2], fsub, gsub, tolerances=tol,
                                        is_linear=False, maximum_mesh_size=2000, is_complex=True, vectorized=True, coarsen_initial_guess_mesh=True, initial_guess=guess)

            self.theta[e[0], :, :] = solution(self.x)
            self.theta_func.append(solution)
        self.solved['theta'] = True
        self.E = self.E[::-1]
        self.theta[:, :, :] = self.theta[::-1, :, :]
        self.theta_func = self.theta_func[::-1]


# The solver for fL:

    def solvefl(self):
        boundaries = [0, self.geometry.L]
        tol = [1e-8, 1e-8]
        fl_func = []

        for e in enumerate(self.E):
            theta = self.theta_func[e[0]]

            def Dl(x): return cos(imag(theta(array(x, ndmin=1))[:, 0]))**2

            def dDl(x): return -2*imag(theta(array(x, ndmin=1))[:, 1])*sin(
                imag(theta(array(x, ndmin=1)))[:, 0])*cos(imag(theta(array(x, ndmin=1)))[:, 0])

            def Ns(x): return real(cosh(theta(array(x, ndmin=1))[:, 0]))

            def fsub(x, z):
                fl, dfl = z
                return array([-dDl(x)*dfl])

            def gsub(z):
                fl, dfl = z
                return array([dfl[0]*Dl(0)*self.geometry.RL-Ns(0)*(fl[0]-self.geometry.flL(e[1])),
                              Dl(self.geometry.L)*dfl[1]*self.geometry.RR-Ns(self.geometry.L)*(self.geometry.flR(e[1])-fl[1])])

            def guess(x):
                z = array([[self.geometry.flR(e[1])*xx/self.geometry.L+self.geometry.flL(
                    e[1])*(self.geometry.L-xx)/self.geometry.L, 0.1] for xx in x]).transpose()
                #dm = zeros((2,) + asarray(x).shape)
#                z = zeros((4,) + asarray(x).shape)
                dm = zeros((2,) + asarray(x).shape)
#                z[0] = 1
                return z, dm

            solution = bvp.colnew.solve(boundaries, [
                                        2], fsub, gsub, tolerances=tol, is_linear=True, maximum_mesh_size=10000, is_complex=False, initial_guess=guess)

            fl_func.append(solution)
            self.fl[e[0], :, 0] = solution(self.x)[:, 0]
            self.fl[e[0], :, 1] = solution(self.x)[:, 1]*Dl(self.x)

        self.fl_func = fl_func
        self.solved['fl'] = True

#######
    def solveft(self):
        boundaries = [0, self.geometry.L]
        tol = [1e-8, 1e-8]

        ft_func = []

        for e in enumerate(self.E):
            theta = self.theta_func[e[0]]

            def Dt(x): return cosh(real(theta(array(x, ndmin=1))[:, 0]))**2

            def dDt(x): return 2*real(theta(array(x, ndmin=1))[:, 1])*cosh(
                real(theta(array(x, ndmin=1))[:, 0]))*sinh(real(theta(array(x, ndmin=1))[:, 0]))

            def R(x): return -imag(sinh(theta(array(x, ndmin=1))[:, 0]))

            def Ns(x): return real(cosh(theta(array(x, ndmin=1))[:, 0]))

            def guess(x):
                z = array([[self.geometry.ftR(e[1])*xx/self.geometry.L+self.geometry.ftL(
                    e[1])*(self.geometry.L-xx)/self.geometry.L, 0.1] for xx in x]).transpose()
                #dm = zeros((2,) + asarray(x).shape)
#                z = zeros((4,) + asarray(x).shape)
                dm = zeros((2,) + asarray(x).shape)
#                z[0] = 1
                return z, dm

            def fsub(x, z):
                ft, dft = z
                return array([ft*2.*self.delta(x)*R(x)/Dt(x) - dDt(x)*dft/Dt(x)])

            def gsub(z):
                ft, dft = z
                return array([dft[0]*Dt(0)-Ns(0)*(ft[0]-self.geometry.ftL(e[1]))/self.geometry.RL,
                              Dt(self.geometry.L)*dft[1]-Ns(self.geometry.L)*(self.geometry.ftR(e[1])-ft[1])/self.geometry.RR])

            solution = bvp.colnew.solve(boundaries, [
                                        2], fsub, gsub, tolerances=tol, is_linear=True, maximum_mesh_size=10000, is_complex=False)

            self.ft[e[0], :, 0] = solution(self.x)[:, 0]
            self.ft[e[0], :, 1] = solution(self.x)[:, 1]*Dt(self.x)
            ft_func.append(solution)

        self.ft_func = ft_func
        self.solved['ft'] = True

#####

    def solveAll(self):
        self.solveTheta()
        self.solvefl()
        self.solveft()
        return

    def solveDelta(self):

        delta = zeros((len(self.x),))
        Estar = max(self.E)
        print shape(delta)
        print shape(self.x)
        for x in enumerate(self.x):
            delta[x[0]] = trapz(self.fl[:, x[0], 0] *
                                real(sinh(self.theta[:, x[0], 0])), self.E)
            delta[x[0]] = delta[x[0]]/log(Estar+sqrt(Estar**2-1.))

        self.delta = interp.interp1d(self.x, delta)
        # def f(d,e): return d/(sqrt(e**2-d**2))

        # def deltaIntegral(d, x):
        #     fl_h = array([fll(x)[0] for fll in self.fl_func])
        #     ft_h = array([ftt(x)[0] for ftt in self.ft_func])

        #     integrand = real(f(d,self.E))*(1-fl_h)-1j*imag(f(d,self.E))*ft_h

        #     return integrate.trapz(integrand, self.E)

        # delta = []
        # for x in enumerate(self.x):
        #     delta.append(exp(-deltaIntegral(self.delta(x[1]), x[1])))

        # print array(delta).shape
        # print x.shape
        # return interp.interp1d(self.x, delta)

    def heatCurrent(self):

        if not self.solved['theta']:
            self.solveTheta()
        if self.approx_hc:
            # Spectral coefficients at the interface (assumed to stay constant)
            Dt = cosh(real(self.theta[:, 0, 0]))**2
            Dl = cos(imag(self.theta[:, 0, 0]))**2
            R = - imag(sinh(self.theta[:, 0, 0]))
            # DOS at the interface and at the end of the wire
            Ns0 = real(cosh(self.theta[:, 0, 0]))
            NsL = real(cosh(self.theta[:, -1, 0]))
            K = sqrt(R/Dt)
            # Approximate heat currents from these
            jl = (self.geometry.flR(self.E) - self.geometry.flL(self.E)) * Ns0*NsL / \
                (self.geometry.RR*Ns0 + self.geometry.RL *
                 NsL + Ns0*NsL*self.geometry.L/Dl)
            jt = Dt*(self.geometry.ftR(self.E)*NsL*Ns0*K - K*self.geometry.ftL(self.E)*Ns0 * (cosh(K*self.geometry.L)*NsL + Dt*K*self.geometry.RR*sinh(K*self.geometry.L))) / \
                (Dt*K*cosh(K*self.geometry.L)*(self.geometry.RR*Ns0+self.geometry.RL*NsL)+(Dt **
                                                                                           2 * K**2 * self.geometry.RL*self.geometry.RR+Ns0*NsL)*sinh(K*self.geometry.L))
            # heat current = integral over E*jl + VL*jt
            integrand = self.E*jl-self.geometry.VL*jt
            return trapz(integrand, self.E)

        if not self.solved['fl']:
            self.solvefl()
        if not self.solved['ft']:
            self.solveft()

        integrand = self.fl[:, 0, 1]*self.E-self.geometry.VL*self.ft[:, 0, 1]
        return trapz(integrand, self.E)

    def chargeCurrent(self):
        if not self.solved['theta']:
            self.solveTheta()
        if self.approx_hc:
            # Spectral coefficients at the interface (assumed to stay constant)
            Dt = cosh(real(self.theta[:, 0, 0]))**2
            Dl = cos(imag(self.theta[:, 0, 0]))**2
            R = - imag(sinh(self.theta[:, 0, 0]))
            # DOS at the interface and at the end of the wire
            Ns0 = real(cosh(self.theta[:, 0, 0]))
            NsL = real(cosh(self.theta[:, -1, 0]))
            K = sqrt(R/Dt)
            jt = Dt*(self.geometry.ftR(self.E)*NsL*Ns0*K - K*self.geometry.ftL(self.E)*Ns0 * (cosh(K*self.geometry.L)*NsL + Dt*K*self.geometry.RR*sinh(K*self.geometry.L))) / \
                (Dt*K*cosh(K*self.geometry.L)*(self.geometry.RR*Ns0+self.geometry.RL*NsL)+(Dt **
                                                                                           2 * K**2 * self.geometry.RL*self.geometry.RR+Ns0*NsL)*sinh(K*self.geometry.L))
            return trapz(jt, self.E)

        if not self.solved['ft']:
            self.solveft()

        integrand = self.ft[:, 0, 1]
        return trapz(integrand, self.E)
