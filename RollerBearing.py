import math
import numpy as np
from scipy import optimize

# bearing = RollerBearing()
# bearing.geometry()
# bearing.stiffness()
# bearing.internal_clearance()
# bearing.roller_profile()  // no profile in BallBearing()
# bearing.info() : print information
# bearing.disp() : calculate relative displacement
# bearing.rating() : print basic ref rating life

# bearing type
# 1 - cylindrical roller bearing, needle roller bearing
# 2 - tapered roller bearing
# 3 - thrust cylindrical roller bearing, thrust needle roller bearing
# 4 - thrust tapered roller bearing


class RollerBearing:

    # Input for disp() : Fx, Fy, M
    # Output of disp() : Delta_r, Delta_psi
    # Input for capacity() : Capacity
    # Output of capacity() : Qci, Qce, qci, qce
    # Output for element() : Delta_j, Psi_j, Delta_j_k, q_j_k
    # Output for load() : qkei, qkee
    # Output for basic_ref_life() : L10r, Pref

    def __init__(self, i, Z, Dwe, Dpw, Lwe, alpha_deg, phi0, bearing_type, ns=30):

        # i : number of rows
        self.i = i
        # Z : number of rolling elements
        self.Z = Z
        # Dwe : element diameter , mm
        self.Dwe = Dwe
        # Dpw : Pitch Circle Diameter , mm
        self.Dpw = Dpw
        # Lwe : effective roller length , mm
        self.Lwe = Lwe
        # alpha : nominal contact angle , degree
        self.alpha_deg = alpha_deg
        self.alpha = math.radians(alpha_deg)
        # phi0 : first element angle , degree
        self.phi0 = phi0
        # type : 0 - radial roller bearing ; 1- thrust roller bearing
        self.type = bearing_type
        # ns : number of laminae
        self.ns = ns

        # calculated in geometry()
        self.gamma = 0.0
        self.phi_deg = np.zeros(self.Z)
        self.phi = np.zeros(self.Z)

        # calculated in stiffness()
        self.cl = 0.0
        self.cs = 0.0

        # defined in internal_clearance()
        self.s = 0.0

        # defined in roller_profile()
        self.xk = np.zeros(self.ns)
        self.profile = np.zeros(self.ns)

        # defined and calculated in disp()
        self.Fr = 0.0
        self.M = 0.0
        self.Delta_r = 0.0
        self.Delta_psi = 0.0
        self.angle_fr = 0.0
        self.phi_cal = np.zeros(self.Z)
        self.phi_cal_deg = np.zeros(self.Z)

        # defined and calculated in capacity()
        self.C = 0.0
        self.Qci = 0.0
        self.Qce = 0.0
        self.qci = 0.0
        self.qce = 0.0

        # defined in concentration of edge stress
        self.f_i = np.zeros(self.ns)
        self.f_e = np.zeros(self.ns)

        # calculated in element()
        self.Delta_j = np.zeros(self.Z)
        self.Psi_j = np.zeros(self.Z)
        self.Delta_j_k = np.zeros((self.Z, self.ns))
        self.q_j_k = np.zeros((self.Z, self.ns))

        # calculated in load()
        self.qkei = np.zeros(self.ns)
        self.qkee = np.zeros(self.ns)

        # calculated in basic_ref_life()
        self.L10r = 0.0
        self.Pref = 0.0

    # bearing geometry
    def geometry(self):
        # phi : angular position of elements , deg
        self.phi_deg = self.phi0 + np.linspace(0.0, 360.0, self.Z, endpoint=False)
        self.phi = np.radians(self.phi_deg)
        # gamma : auxiliary parameter , 1
        self.gamma = self.Dwe * math.cos(self.alpha) / self.Dpw

    # spring constant fro steel material
    def stiffness(self):
        # N / mm**10/9
        self.cl = 35948 * math.pow(self.Lwe, 8/9)
        # N / mm**8/9
        self.cs = 35948 * math.pow(self.Lwe, 8/9) / self.ns

    # define internal clearance
    def internal_clearance(self, s):
        self.s = s

    # roller profile
    # Equation 42 - 44
    def roller_profile(self, use_or_not=0):
        self.xk = (np.arange(self.ns) - self.ns / 2 + 0.5) * self.Lwe / self.ns  # ns 为偶数
        if use_or_not == 1:
            if self.type == 1 or self.type == 3:
                if self.Lwe <= 2.5 * self.Dwe:
                    self.profile = 0.00035 * self.Dwe * np.log(1.0 / (1.0 - np.power(2 * self.xk / self.Lwe, 2)))
                else:
                    for i in list(range(self.ns)):
                        if abs(self.xk[i]) <= (self.Lwe - 2.5 * self.Dwe) / 2.0:
                            self.profile[i] = 0.0
                        else:
                            self.profile[i] = 0.0005 * self.Dwe * math.log(1.0 / 1.0 - math.pow((2.0 * abs(self.xk[i]) - self.Lwe - 2.5 * self.Dwe) / 2.5 / self.Dwe, 2))
            elif self.type == 2 or self.type == 4:
                self.profile = 0.00045 * self.Dwe * np.log(1.0 / (1.0 - np.power(2 * self.xk / self.Lwe, 2)))
        elif use_or_not == 0:
            self.profile = np.zeros(self.ns)

    def info(self):
        print("phi = ", self.phi_deg)
        print("gamma = %.4f" % self.gamma)
        print("cl = %.4f" % self.cl)
        print("cs = %.4f" % self.cs)
        print("s = %.4f" % self.s)
        print("bearing type :", self.type)

    # calculate displacement from equilibrium condition
    def disp(self, Fx, Fy, M):
        # Fr : radial load , N
        self.Fr = math.sqrt(Fx * Fx + Fy * Fy)
        # radial load angle , rad
        if Fx != 0.0:
            self.angle_fr = math.atan(Fy / Fx)
        else:
            self.angle_fr = math.pi / 2.0
        #
        self.phi_cal = self.phi - self.angle_fr
        self.phi_cal_deg = np.degrees(self.phi_cal)
        # M : moment , N*mm
        self.M = M
        # delta0 : initial value
        delta0 = np.array([0.01, 0.01])    # 待考查初值的影响
        # solve the equilibrium function
        sol = optimize.root(self.equilibrium, delta0, method='hybr', tol=1e-50)
        # Delta_r : relative radial displacement , mm
        self.Delta_r = sol.x[0]
        # Delta_psi : total misalignment , rad
        self.Delta_psi = sol.x[1]
        print("Delta_r = %.4f" % self.Delta_r)
        print("Delta_psi = %.4f" % self.Delta_psi)

    # bearing rating
    def rating(self, C):
        self.capacity(C)
        self.concentration_edge_stress()
        self.element()
        self.load()
        self.basic_ref_life()

    # Equation 47-57 : calculate basic dynamic load rating for inner and outer ring
    def capacity(self, C):
        # C : basic dynamic load rating , N
        self.C = C
        temp_1 = math.pow(1 + math.pow(1.038 * pow((1-self.gamma)/(1+self.gamma), 143/108), 9/2), 2/9)
        temp_2 = math.pow(1 + math.pow(1.038 * pow((1-self.gamma)/(1+self.gamma), 143/108), -9/2), 2/9)
        temp_3 = math.pow(1 + math.pow(pow((1-self.gamma)/(1+self.gamma), 143/108), 9/2), 2/9)
        temp_4 = math.pow(1 + math.pow(pow((1-self.gamma)/(1+self.gamma), 143/108), -9/2), 2/9)
        temp_5 = pow(2, 2/9)
        temp_6 = pow(2, 2/9)
        # Qci : rolling element load for the basic dynamic load rating of inner ring or shaft washer , N
        # Qce : rolling element load for the basic dynamic load rating of outer ring or housing washer , N
        if self.type == 1 or self.type == 2:
            const_rambda_v = 0.83
            self.Qci = self.C / const_rambda_v / 0.378 / self.Z / math.cos(self.alpha) / math.pow(self.i, 7/9) * temp_1
            self.Qce = self.C / const_rambda_v / 0.364 / self.Z / math.cos(self.alpha) / math.pow(self.i, 7/9) * temp_2
        if (self.type == 3 or self.type == 4) and self.alpha != 90.0:
            const_rambda_v = 0.73
            self.Qci = self.C / const_rambda_v / self.Z / math.sin(self.alpha) * temp_3
            self.Qce = self.C / const_rambda_v / self.Z / math.sin(self.alpha) * temp_4
        if (self.type == 3 or self.type == 4) and self.alpha == 90.0:
            const_rambda_v = 0.73
            self.Qci = self.C / const_rambda_v / self.Z * temp_5
            self.Qce = self.C / const_rambda_v / self.Z * temp_6
        # qci : basic dynamic load rating of a lamina of the inner ring , N
        # qce : basic dynamic load rating of a lamina of the outer ring , N
        self.qci = self.Qci * pow(1 / self.ns, 7/9)
        self.qce = self.Qce * pow(1 / self.ns, 7/9)
        print("Qci = %.4f" % self.Qci)
        print("Qce = %.4f" % self.Qce)
        print("qci = %.4f" % self.qci)
        print("qce = %.4f" % self.qce)

    # concentration of edge stress // not in BallBearing()
    # Equation 60
    def concentration_edge_stress(self, use_or_not=0):
        if use_or_not == 1:
            self.f_i = np.arange(self.ns) + 1
            self.f_e = np.arange(self.ns) + 1
            self.f_i = 1.0 - 0.01 / np.log(1.985 * np.abs((2 * self.f_i - self.ns - 1) / (2 * self.ns - 2)))
            self.f_e = 1.0 - 0.01 / np.log(1.985 * np.abs((2 * self.f_i - self.ns - 1) / (2 * self.ns - 2)))
        else:
            self.f_i = np.zeros(self.ns) + 1.0
            self.f_e = np.zeros(self.ns) + 1.0

    # Equation 38-41 : calculate deflection and load of rolling element
    # Delta_j : array , elastic deflection of the rolling element j , mm
    # Psi_j : array , total misalignment angle between the raceways of rolling element j , rad
    # Delta_j_k : matrix , elastic deflection of the lamina k of rolling element j , mm
    # q_j_k : matrix , load on lamina k of rolling element j , N
    def element(self):
        self.Delta_j = self.Delta_r * np.cos(self.phi_cal) - 0.5 * self.s
        self.Psi_j = np.arctan(np.tan(self.Delta_psi) * np.cos(self.phi_cal))
        for j in list(range(self.Z)):
            self.Delta_j_k[j, :] = self.Delta_j[j] - self.xk * np.tan(self.Psi_j[j]) - 2 * self.profile
        self.Delta_j_k = np.maximum(self.Delta_j_k, 0.0)
        self.q_j_k = self.cs * np.power(self.Delta_j_k, 10 / 9)
        print("Delta_j_k = ", self.Delta_j_k)
        print("q_j_k", self.q_j_k)

    # Equation 61-64 : calculate dynamic equivalent load on a lamina
    # state = 1 , inner is rotating
    # state = 0 , outer is rotating
    # qkei : dynamic equivalent load on a lamina k of the inner ring , N
    # qkee : dynamic equivalent load on a lamina k of the outer ring , N
    def load(self, state=1):
        temp_1 = np.zeros(self.ns)
        temp_2 = np.zeros(self.ns)
        temp_3 = np.zeros(self.ns)
        temp_4 = np.zeros(self.ns)
        for k in list(range(self.ns)):
            sum_1 = np.power(np.power(self.f_i[k] * self.q_j_k[:, k], 4).sum() / self.Z, 0.25)
            sum_2 = np.power(np.power(self.f_i[k] * self.q_j_k[:, k], 4.5).sum() / self.Z, 1 / 4.5)
            sum_3 = np.power(np.power(self.f_e[k] * self.q_j_k[:, k], 4.5).sum() / self.Z, 1 / 4.5)
            sum_4 = np.power(np.power(self.f_e[k] * self.q_j_k[:, k], 4).sum() / self.Z, 0.25)
            temp_1[k] = sum_1
            temp_2[k] = sum_2
            temp_3[k] = sum_3
            temp_4[k] = sum_4
        if state == 1:
            self.qkei = temp_1
            self.qkee = temp_3
        if state == 0:
            self.qkei = temp_2
            self.qkee = temp_4
        print("qkei = ", self.qkei)
        print("qkee = ", self.qkee)

    # Equation 65 : calculate basic reference rating life
    # L10r : basic reference rating life , 1e6 cycle
    # Pref_r : dynamic equivalent refrence load , N
    def basic_ref_life(self):
        temp = np.power(self.qkei / self.qci, 4.5) + np.power(self.qkee / self.qce, 4.5)
        temp = temp.sum()
        self.L10r = np.power(temp, -8/9)
        self.Pref = self.C / np.power(self.L10r, 3/10)
        print("L10r = %.4f" % self.L10r)
        print("Pref_r = %.4f" % self.Pref)

    # Equation 45-46
    def equilibrium(self, delta):
        temp_sum_1 = 0.0
        temp_sum_2 = 0.0
        #
        Delta_j = delta[0] * np.cos(self.phi_cal) - 0.5 * self.s
        Psi_j = np.arctan(math.tan(delta[1]) * np.cos(self.phi_cal))
        #
        Delta_j_k = np.zeros((self.Z, self.ns))
        for j in list(range(self.Z)):
            Delta_j_k[j, :] = Delta_j[j] - self.xk[:] * np.tan(Psi_j[j]) - 2.0 * self.profile
        Delta_j_k = np.maximum(Delta_j_k, 0.0)
        #
        for j in list(range(self.Z)):
            temp_sum_1 = temp_sum_1 + np.sum(np.cos(self.phi_cal[j]) * np.power(Delta_j_k[j, :], 10/9))
            temp_sum_2 = temp_sum_2 + np.sum(np.cos(self.phi_cal[j]) * np.power(Delta_j_k[j, :], 10/9) * self.xk)
        #
        return [self.Fr - self.cl / self.ns * temp_sum_1,
                self.M - self.cl / self.ns * temp_sum_2]

    def residual(self):
        print("Equilibrium Residual = ", self.equilibrium([self.Delta_r, self.Delta_psi]))


# Example 1 :
# RollerBearing ( i, Z, Dwe, Dpw, Lwe, alpha_deg, phi0, bearing_type, ns=30 )
print("----------------------------------------------------------------")
print("Example 1")
bearing_1 = RollerBearing(1, 6, 9.0, 43.5, 13.6, 0.0, 30.0, 1, 30)
bearing_1.geometry()
bearing_1.stiffness()
bearing_1.internal_clearance(0.0)
bearing_1.roller_profile()
bearing_1.info()
bearing_1.disp(1000.0, 0.0, 0.0)
bearing_1.rating(31500)
bearing_1.residual()
print(bearing_1.phi_cal_deg)
print(bearing_1.Delta_j)
