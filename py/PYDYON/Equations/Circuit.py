# Implementation of circuit equation


from .. import Conductivity


class CircuitEquation:
    

    def __init__(self, quantities, ions, Vloop, Lp=5.4e-6, LMK2=9.1e-6, M=2.49e-6, RMK2=7.5e-4):
        """
        Constructor.
        """
        self.quantities = quantities
        self.ions = ions

        self.L_p = Lp
        self.L_MK2 = LMK2
        self.M = M
        self.R_MK2 = RMK2

        if callable(Vloop):
            self.Vloop = Vloop
        else:
            self.Vloop = lambda t : Vloop


    def __call__(self, x):
        return self.eval(x)


    def evalIp(self, t, x):
        """
        Evaluate dIp/dt.
        """
        pf = self.M - self.L_MK2*self.L_p / self.M
        I_MK2 = self.quantities['IMK2']
        I_p = self.quantities['Ip']
        R_p = Conductivity.evalResistance(self.quantities)

        return (self.Vloop(t) - self.R_MK2*I_MK2 - self.L_MK2/self.M * (Vloop(t) - R_p*I_p)) / pf
    

    def evalIMK2(self, x):
        """
        Evaluate dIMK2/dt.
        """
        pf = self.M - self.L_MK2*self.L_p / self.M
        I_MK2 = self.quantities['IMK2']
        I_p = self.quantities['Ip']
        R_p = Conductivity.evalResistance(self.quantities)
        V = self.Vloop(t)

        return (V - R_p*I_p - self.L_p/self.M * (V - self.R_MK2*I_MK2)) / pf


