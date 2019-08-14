import math
import numpy as np
from scipy import optimize
from scipy import special
from Head import Material, Load, BasicCapacity

# material = Material()
# bearing = BallBearing()
# bearing.geometry(s_r)
# bearing.stiffness(material)
# bearing.info() : print information
# load = Load()
# bearing.disp(load) : calculate relative displacement
# basic_capacity = BasicCapacity()
# bearing.rating(basic_capacity) : print basic ref rating life
# bearing.residual()

# bearing_type
# 1 - deep groove ball bearing, angular contact ball bearing, separable ball bearing
# 2 - spherical roller bearing
# 3 - self aligning ball bearing
# 4 - thrust ball bearing, thrust angular contact ball bearing
# 5 - thrust spherical roller bearing

class BallBearing:

    # Input for disp() : Load
    # Output of disp() : Delta_r, Delta_a, Delta_psi
    # Input for capacity() : BasicCapacity
    # Output of capacity() : Qci, Qce
    # Output for element() : Delta_j, Alpha_j, Q_j
    # Output for load() : Qei, Qee
    # Output for basic_ref_life() : L10r, Pref

    def __init__(self, i, Z, Dw, Dpw, alpha_deg, phi0, bearing_type):

        # i : number of rows
        self.i = i
        # Z : number of rolling elements
        self.Z = Z
        # Dw : Element Diameter , mm
        self.Dw = Dw
        #  Dpw : Pitch Circle Diameter , mm
        self.Dpw = Dpw
        # alpha : nominal contact angle , degree
        self.alpha_deg = alpha_deg
        self.alpha = math.radians(alpha_deg)
        # phi0 : first element angle , degree
        self.phi0 = phi0
        # type
        self.type = bearing_type

        # geometry()
        # input for geometry()
        self.s_r = 0.0
        # calculated in geometry()
        self.ri = 0.0
        self.re = 0.0
        self.phi_deg = np.zeros(self.Z)
        self.phi = np.zeros(self.Z)
        self.A = 0.0
        # calculated in internal clearance()
        self.s_a = 0.0
        self.alpha0 = 0.0
        self.alpha0_deg = 0.0
        self.Ri = 0.0
        # continue calculated in geometry()
        self.gamma = 0.0
        self.sigma_rho_i = 0.0
        self.sigma_rho_e = 0.0
        self.Fi_rho = 0.0
        self.Fe_rho = 0.0
        self.chi_i = 0.0
        self.chi_e = 0.0

        # stiffness()
        # defined in stiffness()
        self.mat_E = 0.0
        self.mat_nu = 0.0
        # calculated in stiffness()
        self.cp = 0.0

        # defined and calculated in disp()
        self.phi_cal = np.zeros(self.Z)
        self.phi_cal_deg = np.zeros(self.Z)
        self.Fr = 0.0
        self.Fa = 0.0
        self.M = 0.0
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
        self.Delta_Element = np.zeros(self.Z)
        self.Alpha_Element = np.zeros(self.Z)
        self.Q_Element = np.zeros(self.Z)

        # calculated in load()
        self.Qei = 0.0
        self.Qee = 0.0

        # calculated in basic_ref_life()
        self.L10r = 0.0
        self.Pref = 0.0

    # bearing geometry
    def geometry(self, s_r):
        # ri : inner groove diameter , mm
        # re: outer groove diameter , mm
        if self.type == 1:
            self.re = 0.53 * self.Dw
            self.ri = 0.52 * self.Dw
        elif self.type == 2:
            self.re = self.Dpw / 2.0 / math.cos(self.alpha) + 0.5 * self.Dw
            self.ri = self.re
            self.Dw = 2 * 0.97 * self.re
        elif self.type == 3:
            self.re = 0.5 * self.Dw * (1.0 + 1.0 / self.gamma)
            self.ri = 0.53 * self.Dw
        elif self.type == 4:
            self.re = 0.54 * self.Dw
            self.ri = 0.54 * self.Dw
        elif self.type == 5:
            self.re = (self.Dpw + self.Dw * math.cos(math.radians(45.0))) / 2 / math.cos(self.alpha)
            self.ri = self.re
            self.Dw = 2 * 0.97 * self.re
        # phi : angular position of elements , deg
        self.phi_deg = self.phi0 + np.linspace(0.0, 360.0, self.Z, endpoint=False)
        self.phi = np.radians(self.phi_deg)
        # A : distance between raceway groove curvature center , mm
        self.A = self.ri + self.re - self.Dw
        # s_r : radial operating clearance , mm
        self.s_r = s_r
        self.__internal_clearance()
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
        self.chi_i = optimize.newton(self.__eq_2_1, 0.1)
        # chi_i : ratio of semi-major to semi-minor of the contact ellipse at outer ring , 1
        self.chi_e = optimize.newton(self.__eq_2_2, 0.1)

    # calculate internal clearance
    def __internal_clearance(self):
        # alpha0 : initial contact angle
        self.alpha0 = math.acos(1.0 - self.s_r / 2.0 / self.A)
        self.alpha0_deg = math.degrees(math.acos(1.0 - self.s_r / 2.0 / self.A))
        # s : axial operating clearance , mm
        self.s_a = 2.0 * self.A * math.sin(self.alpha0)
        # Ri : distance between the center of curvature of the inner race groove and the axis of rotation , mm
        self.Ri = 0.5 * self.Dpw + (self.ri - 0.5 * self.Dw) * math.cos(self.alpha0)

    # bearing stiffness
    def stiffness(self, material):
        # mat_E : modulus of elasticity , MPa
        self.mat_E = material.mat_E
        # mat_nu : poisson ratio
        self.mat_nu = material.mat_nu
        # cp : spring contact , N / mm**1.5
        self.cp = 1.48 * self.mat_E / (1 - self.mat_nu ** 2) * np.power(
            self.__k(self.chi_i) * np.power(self.sigma_rho_i / self.chi_i ** 2 / self.__e(self.chi_i), 1 / 3) +
            self.__k(self.chi_e) * np.power(self.sigma_rho_e / self.chi_e ** 2 / self.__e(self.chi_e), 1 / 3), -1.5)

    # print bearing information
    def info(self):
        print("i  =  %2d" % self.i)
        print("Z  =  %2d" % self.Z)
        print("Dw  =  %.2f " % self.Dw)
        print("Dpw  =  %.2f" % self.Dpw)
        print("alpha_deg  =  %.2f" % self.alpha_deg)
        print("phi0  =  %.2f" % self.phi0)
        print("bearing type  =  %1d" % self.type)
        print("s_r  =  %.4f" % self.s_r)
        print("ri  =  %.4f" % self.ri)
        print("re  =  %.4f" % self.re)
        print("phi  =  ", self.phi_deg)
        print("A  =  %.4f" % self.A)
        print("alpha0_deg  =  %.4f" % self.alpha0_deg)
        print("s_a  =  %.4f" % self.s_a)
        print("Ri  =  %.4f" % self.Ri)
        print("gamma  =  %.4f" % self.gamma)
        print("sigma_rho_i  =  %.4f" % self.sigma_rho_i)
        print("sigma_rho_e  =  %.4f" % self.sigma_rho_e)
        print("Fi_rho  =  %.4f" % self.Fi_rho)
        print("Fe_rho  =  %.4f" % self.Fe_rho)
        print("chi_i  =  %.4f" % self.chi_i)
        print("chi_e  =  %.4f" % self.chi_e)
        print("mat_E  =  %.2f" % self.mat_E)
        print("mat_nu  =  %.2f" % self.mat_nu)
        print("cp  =  %.4f" % self.cp)

    # calculate displacement from equilibrium condition
    def disp(self, load):
        # Fr : radial load , N
        self.Fr = load.Fr
        # radial load angle , rad
        if load.Fx != 0.0:
            self.angle_fr = math.atan(load.Fy / load.Fx)
        else:
            self.angle_fr = math.pi / 2.0
        #
        self.phi_cal = self.phi - self.angle_fr
        self.phi_cal_deg = np.degrees(self.phi_cal)
        # Fa : axial load , N
        self.Fa = load.Fa
        # M : moment , N*mm
        self.M = load.Mb
        # delta0 : initial value
        delta0 = np.array([0.01, 0.01, 0.01])    # 初值不取零；若取零，径向间隙为零、受轴向力情况将得不出结果
        # solve the equilibrium function
        sol = optimize.root(self.__equilibrium, delta0, method='hybr', tol=1e-50)
        # Delta_r : relative radial displacement , mm
        self.Delta_r = sol.x[0]
        # Delta_a : relative axial displacement , mm
        self.Delta_a = sol.x[1]
        self.Delta_a_MASTA = self.Delta_a + self.s_a * 0.5
        # Delta_psi : total misalignment , rad
        self.Delta_psi = sol.x[2]
        print("Delta_r = %.8f" % self.Delta_r)
        print("Delta_a = %.8f" % self.Delta_a_MASTA)
        print("Delta_psi = %.8f" % self.Delta_psi)

    # bearing rating
    def rating(self, basic_capacity):
        self.__capacity(basic_capacity)
        self.__element()
        self.__load()
        self.__basic_ref_life()

    # Equation 19-24 : calculate basic dynamic load rating for inner and outer ring
    def __capacity(self, basic_capacity):
        # C : basic dynamic load rating , N
        self.C = basic_capacity.C
        temp_1 = math.pow(1 + math.pow(1.044 * pow((1-self.gamma)/(1+self.gamma), 1.72) * pow(self.ri/self.re*(2*self.re-self.Dw)/(2*self.ri-self.Dw), 0.41), 10/3), 0.3)
        temp_2 = math.pow(1 + math.pow(1.044 * pow((1-self.gamma)/(1+self.gamma), 1.72) * pow(self.ri/self.re*(2*self.re-self.Dw)/(2*self.ri-self.Dw), 0.41), -10/3), 0.3)
        temp_3 = math.pow(1 + math.pow(pow((1-self.gamma)/(1+self.gamma), 1.72) * pow(self.ri/self.re*(2*self.re-self.Dw)/(2*self.ri-self.Dw), 0.41), 10/3), 0.3)
        temp_4 = math.pow(1 + math.pow(pow((1-self.gamma)/(1+self.gamma), 1.72) * pow(self.ri/self.re*(2*self.re-self.Dw)/(2*self.ri-self.Dw), 0.41), -10/3), 0.3)
        temp_5 = math.pow(1 + math.pow(pow(self.ri/self.re*(2*self.re-self.Dw)/(2*self.ri-self.Dw), 0.41), 10/3), 0.3)
        temp_6 = math.pow(1 + math.pow(pow(self.ri/self.re*(2*self.re-self.Dw)/(2*self.ri-self.Dw), 0.41), -10/3), 0.3)
        # Qci : rolling element load for the basic dynamic load rating of inner ring or shaft washer , N
        # Qce : rolling element load for the basic dynamic load rating of outer ring or housing washer , N
        if self.type == 1 or self.type == 2 or self.type == 3:
            self.Qci = self.C / 0.407 / self.Z / math.pow(math.cos(self.alpha), 0.7) / math.pow(self.i, 0.7) * temp_1
            self.Qce = self.C / 0.389 / self.Z / math.pow(math.cos(self.alpha), 0.7) / math.pow(self.i, 0.7) * temp_2
        if (self.type == 4 or self.type == 5) and self.alpha_deg != 90.0:
            self.Qci = self.C / self.Z / math.sin(self.alpha) * temp_3
            self.Qce = self.C / self.Z / math.sin(self.alpha) * temp_4
        if (self.type == 4 or self.type == 5) and self.alpha_deg == 90.0:
            self.Qci = self.C / self.Z * temp_5
            self.Qce = self.C / self.Z * temp_6
        print("Qci = %.4f" % self.Qci)
        print("Qce = %.4f" % self.Qce)

    # Equation 12 and Equation 10 : calculate deflection and load of rolling element
    # Delta_Element : elastic deflection of rolling element , mm
    # Alpha_Element : operating contact angle of rolling element , deg
    # Q_Element : load of rolling element , N
    def __element(self):
        self.Delta_Element = np.sqrt((self.A * math.cos(self.alpha0) + self.Delta_r * np.cos(self.phi_cal)) ** 2 +
                                     (self.A * math.sin(self.alpha0) + self.Delta_a + self.Ri * math.sin(self.Delta_psi) *
                                      np.cos(self.phi_cal)) ** 2) - self.A
        self.Delta_Element = np.maximum(self.Delta_Element, 0.0)
        self.Alpha_Element = np.arctan(
                (self.A * math.sin(self.alpha0) + self.Delta_a + self.Ri * math.sin(self.Delta_psi) * np.cos(self.phi_cal)) /
                (self.A * math.cos(self.alpha0) + self.Delta_r * np.cos(self.phi_cal)))
        self.Alpha_Element = np.degrees(self.Alpha_Element)
        self.Q_Element = self.cp * np.power(self.Delta_Element, 1.5)
        print("Delta_Element = ", self.Delta_Element)
        print("Alpha_Element = ", self.Alpha_Element)
        print("Q_Element = ", self.Q_Element)

    # Equation 25-28 : calculate dynamic equivalent load
    # state = 1 , inner is rotating
    # state = 0 , outer is rotating
    # Qei : dynamic equivalent rolling element load on inner ring or shaft washer , N
    # Qee : dynamic equivalent rolling element load on outer ring or housing washer , N
    def __load(self, state=1):
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
    def __basic_ref_life(self):
        self.L10r = math.pow(math.pow(self.Qci/self.Qei, -10/3) + math.pow(self.Qce/self.Qee, -10/3), -9/10)
        self.Pref = self.C / math.pow(self.L10r, 1/3)
        print("L10r = %.4f" % self.L10r)
        print("Pref_r = %.4f" % self.Pref)

    # Equation 3 : complete elliptic integral of the first kind
    def __k(self, chi):
        return special.ellipk(1 - 1.0 / chi ** 2)

    # Equation 4 : complete elliptic integral of the second kind
    def __e(self, chi):
        return special.ellipe(1 - 1.0 / chi ** 2)

    # Equation 2_1 : calculate chi_i
    def __eq_2_1(self, chi):
        return 1 - 2.0 / (chi ** 2 - 1.0) * (self.__k(chi) / self.__e(chi) - 1.0) - self.Fi_rho

    # Equation 2_1 : calculate chi_e
    def __eq_2_2(self, chi):
        return 1 - 2.0 / (chi ** 2 - 1.0) * (self.__k(chi) / self.__e(chi) - 1.0) - self.Fe_rho

    # Equation 16-18 : equilibrium condition
    def __equilibrium(self, delta):
        delta_j = np.sqrt((self.A * math.cos(self.alpha0) + delta[0] * np.cos(self.phi_cal)) ** 2 +
                          (self.A * math.sin(self.alpha0) + delta[1] + self.Ri * math.sin(delta[2]) * np.cos(
                           self.phi_cal)) ** 2) - self.A
        delta_j = np.maximum(delta_j, 0)
        alpha_j = np.arctan(
            (self.A * math.sin(self.alpha0) + delta[1] + self.Ri * math.sin(delta[2]) * np.cos(self.phi_cal)) /
            (self.A * math.cos(self.alpha0) + delta[0] * np.cos(self.phi_cal)))
        temp_sum_1 = np.power(delta_j, 1.5) * np.cos(alpha_j) * np.cos(self.phi_cal)
        temp_sum_1 = temp_sum_1.sum()
        temp_sum_2 = np.power(delta_j, 1.5) * np.sin(alpha_j)
        temp_sum_2 = temp_sum_2.sum()
        temp_sum_3 = np.power(delta_j, 1.5) * np.sin(alpha_j) * np.cos(self.phi_cal)
        temp_sum_3 = temp_sum_3.sum()
        return [self.Fr - self.cp * temp_sum_1,
                self.Fa - self.cp * temp_sum_2,
                self.M - self.Dpw / 2.0 * self.cp * temp_sum_3]

    def residual(self):
        print("Equilibrium Residual = ", self.__equilibrium([self.Delta_r, self.Delta_a, self.Delta_psi]))


# mode = 1 test local
mode = 1
if mode == 1:
    # Example 1 : 深沟球轴承，径向间隙为零，只受径向力
    # BallBearing ( i, Z, Dw, Dpw, alpha_deg, phi0, type )
    # Test with 6305
    print("----------------------------------------------------------------")
    print("Example")
    material_1 = Material("Steel", 210000.0, 0.3)
    bearing_1 = BallBearing(1, 7, 11.5, 43.5, 0.0, 0.0, 1)
    bearing_1.geometry(0.01)
    bearing_1.stiffness(material_1)
    bearing_1.info()
    print("###")
    load_1 = Load(1000.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    load_2 = Load(0.0, 0.0, 1000.0, 0.0, 0.0, 0.0)
    load_3 = Load(1000.0, 0.0, 1000.0, 0.0, 0.0, 0.0)
    load_4 = Load(0.0, 0.0, 0.0, 1000.0, 0.0, 0.0)
    load_5 = Load(1000.0, 0.0, 1000.0, 1000.0, 0.0, 0.0)
    bearing_1.disp(load_5)
    print("###")
    basic_capacity_1 = BasicCapacity(23400)
    bearing_1.rating(basic_capacity_1)
    print("###")
    bearing_1.residual()

