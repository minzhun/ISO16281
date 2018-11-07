import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy import special

# bearing = BallBearing()
# bearing.geometry()
# bearing.material()
# bearing.stiffness()
# bearing.internal_clearance()
# bearing.info() : print information
# bearing.disp() : calculate relative displacement
# bearing.rating() : print basic ref rating life


class BallBearing:

    # Input for disp() : Fx, Fy, Fz, Mz
    # Output of disp() : Delta_r, Delta_a, Delta_psi
    # Input for capacity() : Capacity, Bearing type
    # Output of capacity() : Qci, Qce
    # Output for element() : delta_j, alpha_j, Q_j
    # Output for load() : Qei, Qee
    # Output for basic_ref_life() : L10r, Pref_r

    def __init__(self, i, Z, Dw, Dpw, ri, re, alpha, phi0, bearing_type):

        # i : number of rows
        self.i = i
        # Z : number of rolling elements
        self.Z = Z
        # Dw : Element Diameter , mm
        self.Dw = Dw
        #  Dpw : Pitch Circle Diameter , mm
        self.Dpw = Dpw
        # ri : inner groove diameter , mm
        self.ri = ri
        # re: outer groove diameter , mm
        self.re = re
        # alpha : nominal contact angle , degree
        self.alpha = alpha
        self.alpha_rad = math.radians(alpha)
        # phi0 : first element angle , degree
        self.phi0 = phi0
        # type : 0 - radial ball bearing ; 1- thrust ball bearing
        self.type = bearing_type

        # calculated in geometry()
        self.phi = list(range(self.Z))
        self.phi_cal = list(range(self.Z))
        self.A = 0.0
        self.gamma = 0.0
        self.sigma_rho_i = 0.0
        self.sigma_rho_e = 0.0
        self.Fi_rho = 0.0
        self.Fe_rho = 0.0
        self.chi_i = 0.0
        self.chi_e = 0.0

        # defined in material()
        self.mat_E = 0.0
        self.mat_nu = 0.0

        # calculated in stiffness()
        self.cp = 0.0

        # defined in internal clearance()
        self.s_r = 0.0
        self.s_a = 0.0
        self.alpha0 = 0.0
        self.alpha0_deg = 0.0
        self.Ri = 0.0

        # defined and calculated in disp()
        self.Fr = 0.0
        self.Fa = 0.0
        self.Mz = 0.0
        self.angle_fr = 0.0
        self.Delta_r = 0.0
        self.Delta_a = 0.0
        self.Delta_a_MASTA = 0.0
        self.Delta_psi = 0.0

        # defined and calculated in capacity()
        self.C = 0.0
        self.Qci = 0.0
        self.Qce = 0.0

        # calculated in element()
        self.Delta_Element = list(range(Z))
        self.Alpha_Element = list(range(Z))
        self.Q_Element = list(range(Z))

        # calculated in load()
        self.Qei = 0.0
        self.Qee = 0.0

        # calculated in basic_ref_life()
        self.L10r = 0.0
        self.Pref_r = 0.0

    # bearing geometry
    def geometry(self):
        # phi : angular position of elements , deg
        for i in self.phi:
            self.phi[i] = self.phi0 + self.phi[i] * 360.0 / self.Z
        # A : distance between raceway groove curvature center , mm
        self.A = self.ri + self.re - self.Dw
        # gamma : auxiliary parameter , 1
        self.gamma = self.Dw * math.cos(self.alpha) / self.Dpw
        # sigma_rho_i : curvature sum at inner ring contact , 1/mm
        self.sigma_rho_i = 2.0 / self.Dw * (2.0 + self.gamma / (1 - self.gamma) - self.Dw / 2.0 / self.ri)
        # sigma_rho_e : curvature sum at outer ring contact , 1/mm
        self.sigma_rho_e = 2.0 / self.Dw * (2.0 - self.gamma / (1 + self.gamma) - self.Dw / 2.0 / self.re)
        # Fi_rho : relative curvature difference at the inner ring contact , 1
        self.Fi_rho = (self.gamma / (1 - self.gamma) + self.Dw / 2.0 / self.ri) / (2.0 + self.gamma / (1 - self.gamma) - self.Dw / 2.0 / self.ri)
        # Fe_rho : relative curvature difference at the outer ring contact , 1
        self.Fe_rho = (-1.0 * self.gamma / (1 + self.gamma) + self.Dw / 2.0 / self.re) / (2.0 - self.gamma / (1 + self.gamma) - self.Dw / 2.0 / self.re)
        # chi_i : ratio of semi-major to semi-minor of the contact ellipse at inner ring , 1
        self.chi_i = optimize.newton(self.eq_2_1, 0.1)
        # chi_i : ratio of semi-major to semi-minor of the contact ellipse at outer ring , 1
        self.chi_e = optimize.newton(self.eq_2_2, 0.1)

    # bearing material
    def material(self, mat_E=210000.0, mat_nu=0.3):
        # mat_E : modulus of elasticity , MPa
        self.mat_E = mat_E
        # mat_nu : poisson ratio
        self.mat_nu = mat_nu

    # bearing stiffness
    def stiffness(self):
        # cp : spring contact , N / mm**1.5
        self.cp = 1.48 * self.mat_E / (1 - self.mat_nu ** 2) * np.power(
            self.k(self.chi_i) * np.power(self.sigma_rho_i / self.chi_i ** 2 / self.e(self.chi_i), 1/3) +
            self.k(self.chi_e) * np.power(self.sigma_rho_e / self.chi_e ** 2 / self.e(self.chi_e), 1/3), -1.5)

    # calculate internal clearance
    def internal_clearance(self, s_r):
        # s_r : radial operating clearance , mm
        self.s_r = s_r
        # alpha0 : initial contact angle
        self.alpha0 = math.acos(1.0 - self.s_r / 2.0 / self.A)
        self.alpha0_deg = math.degrees(math.acos(1.0 - self.s_r / 2.0 / self.A))
        # s : axial operating clearance , mm
        self.s_a = 2.0 * self.A * math.sin(self.alpha0)
        # Ri : distance between the center of curvature of the inner race groove and the axis of rotation , mm
        self.Ri = 0.5 * self.Dpw + (self.ri - 0.5 * self.Dw) * math.cos(self.alpha0)

    # print bearing information
    def info(self):
        print("phi = ", self.phi)
        print("A = %.4f" % self.A)
        print("alpha0_deg = %.4f" % self.alpha0_deg)
        print("s_a = %.4f" % self.s_a)
        print("gamma = %.4f" % self.gamma)
        print("sigma_rho_i = %.4f" % self.sigma_rho_i)
        print("sigma_rho_e = %.4f" % self.sigma_rho_e)
        print("Fi_rho = %.4f" % self.Fi_rho)
        print("Fe_rho = %.4f" % self.Fe_rho)
        print("chi_i = %.4f" % self.chi_i)
        print("chi_e = %.4f" % self.chi_e)
        print("Ri = %.4f" % self.Ri)
        print("cp = %.4f" % self.cp)
        print("bearing type :", self.type)

    # calculate displacement from equilibrium condition
    def disp(self, Fx, Fy, Fz, Mz):
        # Fr : radial load , N
        self.Fr = math.sqrt(Fx * Fx + Fy * Fy)
        # radial load angle , deg
        self.angle_fr = math.degrees(math.atan(Fy / Fx))
        #
        for i in self.phi_cal:
            self.phi_cal[i] = self.phi[i] - self.angle_fr
        # Fa : axial load , N
        self.Fa = Fz
        # Mz : moment , N*mm
        self.Mz = Mz
        # delta0 : initial value
        delta0 = np.array([0.01, 0.01, 0.01])    # 初值不取零；若取零，径向间隙为零、受轴向力情况将得不出结果
        # solve the equilibrium function
        sol = optimize.root(self.equilibrium, delta0, method='hybr', tol=1e-50)
        # Delta_r : relative radial displacement , mm
        self.Delta_r = sol.x[0]
        # Delta_a : relative axial displacement , mm
        self.Delta_a = sol.x[1]
        self.Delta_a_MASTA = self.Delta_a + self.s_a * 0.5
        # Delta_psi : total misalignment , rad
        self.Delta_psi = sol.x[2]
        print("Delta_r = %.4f" % self.Delta_r)
        print("Delta_a = %.4f" % self.Delta_a_MASTA)
        print("Delta_psi = %.4f" % self.Delta_psi)

    # bearing rating
    def rating(self, C):
        self.capacity(C)
        self.element()
        self.load()
        self.basic_ref_life()

    # Equation 19-24 : calculate basic dynamic load rating for inner and outer ring
    def capacity(self, C):
        # C : basic dynamic load rating , N
        self.C = C
        temp_1 = math.pow(1+ math.pow(1.044 * pow((1-self.gamma)/(1+self.gamma), 1.72) * pow(self.ri/self.re*(2*self.re-self.Dw)/(2*self.ri-self.Dw), 0.41), 10/3), 0.3)
        temp_2 = math.pow(1+ math.pow(1.044 * pow((1-self.gamma)/(1+self.gamma), 1.72) * pow(self.ri/self.re*(2*self.re-self.Dw)/(2*self.ri-self.Dw), 0.41), -10/3), 0.3)
        temp_3 = math.pow(1+ math.pow(pow((1-self.gamma)/(1+self.gamma), 1.72) * pow(self.ri/self.re*(2*self.re-self.Dw)/(2*self.ri-self.Dw), 0.41), 10/3), 0.3)
        temp_4 = math.pow(1+ math.pow(pow((1-self.gamma)/(1+self.gamma), 1.72) * pow(self.ri/self.re*(2*self.re-self.Dw)/(2*self.ri-self.Dw), 0.41), -10/3), 0.3)
        temp_5 = math.pow(1+ math.pow(pow(self.ri/self.re*(2*self.re-self.Dw)/(2*self.ri-self.Dw), 0.41), 10/3), 0.3)
        temp_6 = math.pow(1+ math.pow(pow(self.ri/self.re*(2*self.re-self.Dw)/(2*self.ri-self.Dw), 0.41), -10/3), 0.3)
        # Qci : rolling element load for the basic dynamic load rating of inner ring or shaft washer , N
        # Qce : rolling element load for the basic dynamic load rating of outer ring or housing washer , N
        if self.type == 0:
            self.Qci = self.C / 0.407 / self.Z / math.cos(self.alpha) / math.pow(self.i, 0.7) * temp_1
            self.Qce = self.C / 0.389 / self.Z / math.cos(self.alpha) / math.pow(self.i, 0.7) * temp_2
        if self.type == 1 and self.alpha != 90.0:
            self.Qci = self.C / self.Z / math.sin(self.alpha_rad) * temp_3
            self.Qce = self.C / self.Z / math.sin(self.alpha_rad) * temp_4
        if self.type == 1 and self.alpha == 90.0:
            self.Qci = self.C / self.Z * temp_5
            self.Qce = self.C / self.Z * temp_6
        print("Qci = %.4f" % self.Qci)
        print("Qce = %.4f" % self.Qce)

    # Equation 12 and Equation 10 : calculate deflection and load of rolling element
    # Delta_Element : elastic deflection of rolling element , mm
    # Alpha_Element : operating contact angle of rolling element , deg
    # Q_Element : load of rolling element , N
    def element(self):
        for i in list(range(self.Z)):
            self.Delta_Element[i] = 0.0
            self.Alpha_Element[i] = 0.0
            self.Q_Element[i] = 0.0
        for j in list(range(self.Z)):
            delta_j = math.sqrt((self.A * math.cos(self.alpha0) + self.Delta_r * math.cos(math.radians(self.phi_cal[j]))) ** 2 +
                                (self.A * math.sin(self.alpha0) + self.Delta_a + self.Ri * math.sin(self.Delta_psi) * math.cos(
                                    math.radians(self.phi_cal[j]))) ** 2) - self.A
            if delta_j < 0.0:
                delta_j = 0.0
            alpha_j = math.atan(
                (self.A * math.sin(self.alpha0) + self.Delta_a + self.Ri * math.sin(self.Delta_psi) * math.cos(math.radians(self.phi_cal[j]))) /
                (self.A * math.cos(self.alpha0) + self.Delta_r * math.cos(math.radians(self.phi_cal[j]))))
            self.Delta_Element[j] = delta_j
            self.Alpha_Element[j] = math.degrees(alpha_j)
            self.Q_Element[j] = self.cp * math.pow(delta_j, 1.5)
        print("Delta_Element = ", self.Delta_Element)
        print("Alpha_Element = ", self.Alpha_Element)
        print("Q_Element = ", self.Q_Element)

    # Equation 25-28 : calculate dynamic equivalent load
    # state = 1 , inner is rotating
    # state = 0 , outer is rotating
    # Qei : dynamic equivalent rolling element load on inner ring or shaft washer , N
    # Qee : dynamic equivalent rolling element load on outer ring or housing washer , N
    def load(self, state=1):
        sum_1 = 0.0
        sum_2 = 0.0
        for Qj in self.Q_Element:
            sum_1 = sum_1 + math.pow(Qj, 3)
            sum_2 = sum_2 + math.pow(Qj, 10/3)
        t_1 = math.pow(sum_1 / self.Z, 1/3)
        t_2 = math.pow(sum_2 / self.Z, 3/10)
        if state == 1:
            self.Qei = t_1
            self.Qee = t_2
        if state == 0:
            self.Qei = t_2
            self.Qee = t_1
        print("Qei = %.4f" % self.Qei)
        print("Qee = %.4f" % self.Qee)

    # Equation 29 : calculate basic reference rating life
    # L10r : basic reference rating life , 1e6 cycle
    # Pref_r : dynamic equivalent refrence load , N
    def basic_ref_life(self):
        self.L10r = math.pow(math.pow(self.Qci/self.Qei, -10/3) + math.pow(self.Qce/self.Qee, -10/3), -9/10)
        self.Pref_r = self.C / math.pow(self.L10r, 1/3)
        print("L10r = %.4f" % self.L10r)
        print("Pref_r = %.4f" % self.Pref_r)

    # Equation 3 : complete elliptic integral of the first kind
    def k(self, chi):
        return special.ellipk(1 - 1.0 / chi ** 2)

    # Equation 4 : complete elliptic integral of the second kind
    def e(self, chi):
        return special.ellipe(1 - 1.0 / chi ** 2)

    # Equation 2_1 : calculate chi_i
    def eq_2_1(self, chi):
        return 1 - 2.0 / (chi ** 2 - 1.0) * (self.k(chi) / self.e(chi) - 1.0) - self.Fi_rho

    # Equation 2_1 : calculate chi_e
    def eq_2_2(self, chi):
        return 1 - 2.0 / (chi ** 2 - 1.0) * (self.k(chi) / self.e(chi) - 1.0) - self.Fe_rho

    # Equation 16-18 : equilibrium condition
    def equilibrium(self, delta):
        temp_sum_1 = 0.0
        temp_sum_2 = 0.0
        temp_sum_3 = 0.0
        for phi_j in self.phi_cal:
            delta_j = math.sqrt((self.A * math.cos(self.alpha0) + delta[0] * math.cos(math.radians(phi_j))) ** 2 +
                                (self.A * math.sin(self.alpha0) + delta[1] + self.Ri * math.sin(delta[2]) * math.cos(
                                    math.radians(phi_j))) ** 2) - self.A
            if delta_j < 0.0:
                delta_j = 0.0
            alpha_j = math.atan(
                (self.A * math.sin(self.alpha0) + delta[1] + self.Ri * math.sin(delta[2]) * math.cos(math.radians(phi_j))) /
                (self.A * math.cos(self.alpha0) + delta[0] * math.cos(math.radians(phi_j))))
            temp_sum_1 = temp_sum_1 + pow(delta_j, 1.5) * math.cos(alpha_j) * math.cos(math.radians(phi_j))
            temp_sum_2 = temp_sum_2 + pow(delta_j, 1.5) * math.sin(alpha_j)
            temp_sum_3 = temp_sum_3 + pow(delta_j, 1.5) * math.sin(alpha_j) * math.cos(math.radians(phi_j))
        return [self.Fr - self.cp * temp_sum_1,
                self.Fa - self.cp * temp_sum_2,
                self.Mz - self.Dpw / 2.0 * self.cp * temp_sum_3]

    def residual(self):
        print("Equilibrium Residual = ", self.equilibrium([self.Delta_r, self.Delta_a, self.Delta_psi]))


# input_data = np.arange(0.01, 20.0, 0.1)
# output_data_1 = bearing.eq_2_1(input_data)
# output_data_2 = bearing.eq_2_2(input_data)
# plt.plot(input_data, output_data_1, color='red')
# plt.plot(input_data, output_data_2, color='blue')
# plt.grid()
# plt.show()

# Example 1 : 深沟球轴承，径向间隙为零，只受径向力
# Delta : 16.9905 um
# L10r : 12008.9256
# BallBearing ( i, Z, Dw, Dpw, ri, re, alpha, phi0, type )
print("----------------------------------------------------------------")
print("Example 1")
bearing_1 = BallBearing(1, 6, 11.1, 43.5, 5.772, 5.883, 0.0, 0.0, 0)
bearing_1.geometry()
bearing_1.material()
bearing_1.stiffness()
bearing_1.internal_clearance(0.0)
bearing_1.info()
bearing_1.disp(500.0, 866.0, 0.0, 0.0)
bearing_1.rating(23400)
bearing_1.residual()
Err1 = (bearing_1.Delta_r * 1000 - 16.9905) / 16.9905 * 100
Err2 = (bearing_1.L10r - 12008.9256) / 12008.9256 * 100
print("Err 1 = %.4f" % Err1, "%")
print("Err 2 = %.4f" % Err2, "%")
print(bearing_1.phi_cal)
