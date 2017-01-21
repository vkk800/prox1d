from scipy import *

class geometry:
    def __init__(self, TL=0.1, TR=0.1, RL=1000, RR=1000, L=4, V=0, sys='NSN'):
        self.TL=abs(TL)
        self.TR=abs(TR)
        self.RL=abs(RL)
        self.RR=abs(RR)
        self.L=abs(L)
        self.V=V
        self.sys=sys
        self.VL = self.RL*self.V/(self.RL+self.RR)
        self.VR = -self.RR*self.V/(self.RL+self.RR)

    def Delta0(self,X):
        if self.sys[1] == 'N': return 0.
        elif self.sys[1] == 'S': return 1.
        else:
            print 'Error: unknown label for wire'
            return None

    def thetaL(self, E):
        if self.sys[0] == 'N': return 0.
        elif self.sys[0] == 'S': return arctanh(1/E)
        else:
            print 'Error: unknown label for left reservoir'
            return None

    def thetaR(self, E):
        if self.sys[2] == 'N': return 0.
        elif self.sys[2] == 'S': return arctanh(1/E)
        else:
            print 'Error: unknown label for right reservoir'
            return None

    def flL(self, E): return 0.5*(tanh((E-self.VL)/(2*self.TL))+tanh((E+self.VL)/(2*self.TL)))
    def flR(self, E): return 0.5*(tanh((E-self.VR)/(2*self.TR))+tanh((E+self.VR)/(2*self.TR)))
#    def flR(self, E): return tanh(E/(2*self.TR))
    def ftL(self, E): return 0.5*(tanh((E-self.VL)/(2*self.TL))-tanh((E+self.VL)/(2*self.TL)))
    def ftR(self, E): return 0.5*(tanh((E-self.VR)/(2*self.TR))-tanh((E+self.VR)/(2*self.TR)))
#    def ftR(self, E): return 0.
    def fL(self, E): return 1/(exp((E-self.VL)/self.TL)+1)
    def fR(self, E): return 1/(exp((E-self.VR)/self.TR)+1)
