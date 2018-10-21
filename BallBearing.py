import math
import numpy as np
from scipy import optimize
from scipy import special


class BallBearing:

    mat_E = 210000.0
    mat_nu = 0.3

    def __init__(self, Z, Dw, Dpw, ri, re, alpha, phi0, s):
        self.Z = Z
        self.Dw = Dw
        self.Dpw = Dpw
        self.ri = ri
        self.re = re
        self.alpha = alpha
        self.alpha_rad = math.radians(alpha)
        self.phi0 = phi0
        self.s = s
        self.Fr = 0.0
        self.Fa = 0.0
        self.Mz = 0.0
        self.Delta_r = 0.0
        self.Delta_a = 0.0
        self.Delta_psi = 0.0
        self.Cr = 0.0
        self.Ca = 0.0
        self.Qci = 0.0
        self.Qce = 0.0
        self.Qei = 0.0
        self.Qee = 0.0
        self.L10r = 0.0

        self.phi = list(range(Z))
        for i in self.phi:
            self.phi[i] = self.phi0 + self.phi[i] * 360.0 / self.Z
        print("phi = ", self.phi)

        self.Delta_Element = list(range(Z))
        self.Q_Element = list(range(Z))
        for i in list(range(Z)):
            self.Delta_Element[i] = 0.0
            self.Q_Element[i] = 0.0

        self.A = self.ri + self.re - self.Dw
        print("A = %.4f" % self.A)

        self.alpha0 = math.acos(1.0 - self.s / 2.0 / self.A)
        self.alpha0_deg = math.degrees(math.acos(1.0 - self.s / 2.0 / self.A))
        print("alpha0 = %.4f" % self.alpha0_deg)

        self.gamma = self.Dw * math.cos(self.alpha) / self.Dpw
        print("gamma = %.4f" % self.gamma)

        self.sigma_rho_i = 2.0 / self.Dw * (2.0 + self.gamma / (1 - self.gamma) - self.Dw / 2.0 / self.ri)
        print("sigma_rho_i = %.4f" % self.sigma_rho_i)

        self.sigma_rho_e = 2.0 / self.Dw * (2.0 - self.gamma / (1 + self.gamma) - self.Dw / 2.0 / self.re)
        print("sigma_rho_e = %.4f" % self.sigma_rho_e)

        self.Fi_rho = (self.gamma / (1 - self.gamma) + self.Dw / 2.0 / self.ri) / (2.0 + self.gamma / (1 - self.gamma) - self.Dw / 2.0 / self.ri)
        print("Fi_rho = %.4f" % self.Fi_rho)

        self.Fe_rho = (-1.0 * self.gamma / (1 + self.gamma) + self.Dw / 2.0 / self.re) / (2.0 - self.gamma / (1 + self.gamma) - self.Dw / 2.0 / self.re)
        print("Fe_rho = %.4f" % self.Fe_rho)

        self.chi_i = optimize.newton(self.eq_2_1, 0.1)
        print("chi_i = %.4f" % self.chi_i)

        self.chi_e = optimize.newton(self.eq_2_2, 0.1)
        print("chi_e = %.4f" % self.chi_e)

        self.Ri = 0.5 * self.Dpw + (self.ri - 0.5 * self.Dw) * math.cos(self.alpha0)
        print("Ri = %.4f" % self.Ri)

        self.cp = 1.48 * BallBearing.mat_E / (1 - BallBearing.mat_nu ** 2) * np.power(
            self.k(self.chi_i) * np.power(self.sigma_rho_i / self.chi_i ** 2 / self.e(self.chi_i), 1/3) +
            self.k(self.chi_e) * np.power(self.sigma_rho_e / self.chi_e ** 2 / self.e(self.chi_e), 1/3), -1.5)
        print("cp = %.4f" % self.cp)

    def k(self, chi):
        return special.ellipk(1 - 1.0 / chi ** 2)

    def e(self, chi):
        return special.ellipe(1 - 1.0 / chi ** 2)

    def eq_2_1(self, chi):
        return 1 - 2.0 / (chi ** 2 - 1.0) * (self.k(chi) / self.e(chi) - 1.0) - self.Fi_rho

    def eq_2_2(self, chi):
        return 1 - 2.0 / (chi ** 2 - 1.0) * (self.k(chi) / self.e(chi) - 1.0) - self.Fe_rho

    def delta(self, Fr, Fa, Mz):
        self.Fr = Fr
        self.Fa = Fa
        self.Mz = Mz
        delta0 = np.array([0.0, 0.0, 0.0])
        sol = optimize.root(self.delta_func, delta0)
        self.Delta_r = sol.x[0] * 1000
        self.Delta_a = sol.x[1] * 1000
        self.Delta_psi = sol.x[2] * 1000
        print("Fr = %.4f" % self.Fr)
        print("Fa = %.4f" % self.Fa)
        print("Mz = %.4f" % self.Mz)
        print("Delta_r = %.4f" % self.Delta_r)
        print("Delta_a = %.4f" % self.Delta_a)
        print("Delta_psi = %.4f" % self.Delta_psi)

    def element(self):
        for j in list(range(self.Z)):
            delta_j = math.sqrt((self.A * math.cos(self.alpha0) + self.Delta_r * math.cos(math.radians(self.phi[j]))) ** 2 +
                                (self.A * math.sin(self.alpha0) + self.Delta_a + self.Ri * math.sin(self.Delta_psi) * math.cos(
                                    math.radians(self.phi[j]))) ** 2) - self.A
            if delta_j < 0.0:
                delta_j = 0.0
            self.Delta_Element[j] = delta_j
            self.Q_Element[j] = self.cp * math.pow(delta_j/1000, 1.5)
        print("Delta_Element = ", self.Delta_Element)
        print("Q_Element = ", self.Q_Element)

    def delta_func(self, delta):
        temp_sum_1 = 0.0
        temp_sum_2 = 0.0
        temp_sum_3 = 0.0
        for phi_j in self.phi:
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

    def capacity(self, Cr, Ca):
        self.Cr = Cr
        self.Ca = Ca
        temp_1 = math.pow(1+ math.pow(1.044 * pow((1-self.gamma)/(1+self.gamma), 1.72) * pow(self.ri/self.re*(2*self.re-self.Dw)/(2*self.ri-self.Dw), 0.41), 10/3), 0.3)
        temp_2 = math.pow(1+ math.pow(1.044 * pow((1-self.gamma)/(1+self.gamma), 1.72) * pow(self.ri/self.re*(2*self.re-self.Dw)/(2*self.ri-self.Dw), 0.41), -10/3), 0.3)
        temp_3 = math.pow(1+ math.pow(pow((1-self.gamma)/(1+self.gamma), 1.72) * pow(self.ri/self.re*(2*self.re-self.Dw)/(2*self.ri-self.Dw), 0.41), 10/3), 0.3)
        temp_4 = math.pow(1+ math.pow(pow((1-self.gamma)/(1+self.gamma), 1.72) * pow(self.ri/self.re*(2*self.re-self.Dw)/(2*self.ri-self.Dw), 0.41), -10/3), 0.3)
        temp_5 = math.pow(1+ math.pow(pow(self.ri/self.re*(2*self.re-self.Dw)/(2*self.ri-self.Dw), 0.41), 10/3), 0.3)
        temp_6 = math.pow(1+ math.pow(pow(self.ri/self.re*(2*self.re-self.Dw)/(2*self.ri-self.Dw), 0.41), -10/3), 0.3)
        if self.Ca == 0.0:
            self.Qci = self.Cr / 0.407 / self.Z / math.pow(math.cos(self.alpha_rad), 0.7) * temp_1
            self.Qce = self.Cr / 0.389 / self.Z / math.pow(math.cos(self.alpha_rad), 0.7) * temp_2
        if self.Cr == 0.0 and self.alpha != 90.0:
            self.Qci = self.Ca / self.Z / math.sin(self.alpha_rad) * temp_3
            self.Qce = self.Ca / self.Z / math.sin(self.alpha_rad) * temp_4
        if self.Cr == 0.0 and self.alpha == 90.0:
            self.Qci = self.Ca / self.Z * temp_5
            self.Qce = self.Ca / self.Z * temp_6
        print("Qci = %.4f" % self.Qci)
        print("Qce = %.4f" % self.Qce)

    def load(self, state):
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

    def basic_ref_life(self):
        self.L10r = math.pow(math.pow(self.Qci/self.Qei, -10/3) + math.pow(self.Qce/self.Qee, -10/3), -9/10)
        print("L10r = %.4f" % self.L10r)


bearing = BallBearing(7, 11.1, 43.5, 5.772, 5.883, 0.0, 0.0, 0.0)
bearing.delta(1000.0, 0.0, 0.0)
bearing.element()
bearing.capacity(23400, 0.0)
bearing.load(0)
bearing.basic_ref_life()
