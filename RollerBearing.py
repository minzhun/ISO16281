import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy import special

# bearing = RollerBearing()
# bearing.geometry()
# bearing.material()
# bearing.stiffness()
# bearing.internal_clearance()
# bearing.info() : print information
# bearing.disp() : calculate relative displacement
# bearing.rating() : print basic ref rating life


class RollerBearing:

    def __init__(self, i, Z, Dwe, Dpw, Lwe, alpha, phi0, bearing_type, ns=30):

        # i : number of rows
        self.i = i
        # Z : number of rolling elements
        self.Z = Z
        # Dwe : roller diameter , mm
        self.Dwe = Dwe
        # Dpw : Pitch Circle Diameter , mm
        self.Dpw = Dpw
        # Lwe : effective roller length , mm
        self.Lwe = Lwe
        # alpha : nominal contact angle , degree
        self.alpha = alpha
        self.alpha_rad = math.radians(alpha)
        # phi0 : first element angle , degree
        self.phi0 = phi0
        # type : 0 - radial ball bearing ; 1- thrust ball bearing
        self.type = bearing_type
        # ns : number of laminae
        self.ns = ns

        # calculated in geometry()
        self.gamma = 0.0

        # calculated in stiffness()
        self.cl = 0.0
        self.cs = 0.0

        # defined in internal_clearance()
        self.s = 0.0

        # defined and calculated in disp()
        self.Fr = 0.0
        self.Mz = 0.0
        self.Delta_r = 0.0
        self.Delta_phi = 0.0

        # defined and calculated in capacity()
        self.C = 0.0
        self.Qci = 0.0
        self.Qce = 0.0
        self.qci = 0.0
        self.qce = 0.0

        # calculated in element()
        self.Delta_j_k = 0.0
        self.Q_j_k = 0.0

        # calculated in load()
        self.qkei = 0.0
        self.qkee = 0.0

        # calculated in basic_ref_life()
        self.L10r = 0.0
        self.Pref_r = 0.0

    def geometry(self):
        pass

    def stiffness(self):
        pass

    def roller_profile(self):
        pass

    def internal_clearance(self):
        pass

    def info(self):
        pass

    def disp(self):
        pass

    def rating(self, C):
        pass

    def capacity(self, C):
        pass

    def concentration_edge(self):
        pass

    def element(self):
        pass

    def load(self):
        pass

    def basic_ref_life(self):
        pass

    def equilibrium(self):
        pass
