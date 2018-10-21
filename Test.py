import math
import numpy as np
from scipy import optimize
from scipy import special
import matplotlib.pyplot as plt

# Z : number of rolling elements
Z = 7

# mat_E : modulus of elasticity , MPa
mat_E = 210000.0

# mat_nu : poisson ratio
mat_nu = 0.3

# Dw : Element Diameter , mm
Dw = 11.1

#  Dpw : Pitch Circle Diameter , mm
Dpw = 43.5

# ri : inner groove diameter , mm
ri = 5.772

# re: outer groove diameter , mm
re = 5.883

# alpha : nominal contact angle , degree
alpha = 0.0

# s : radial operating clearance , mm
s = 0.0

# phi : angular position , degree
phi = [0.0, 360.0/7, 360.0/7*2, 360.0/7*3, 360.0/7*4, 360.0/7*5, 360.0/7*6]

# Fr : radial load , N
Fr = 1000.0

# Fa : axial load , N
Fa = 0.0

# Mz : moment , N*mm
# Mz = math.sqrt(2.0675**2+0.4658**2) * 1000
Mz = 0.0

# A : distance between the curvature centres of the raceway groove radii , mm
A = ri + re - Dw
print("A = %12.4f" % A)

# alpha_0 : initial contact angle , degree
alpha0_deg = math.degrees(math.acos(1.0 - s / 2.0 / A))
alpha0 = math.acos(1.0 - s / 2.0 / A)
print("alpha_0 = %12.4f" % alpha0_deg)

# Ri : distance between the centre of curvature of the inner raceway groove and the axis of rotation , mm
Ri = 0.5 * Dpw + (ri - 0.5 * Dw) * math.cos(alpha0)
print("Ri = %12.4f" % Ri)

# gamma : auxiliary parameter
gamma = Dw * math.cos(math.radians(alpha)) / Dpw
print("gamma = %12.4f" % gamma)

# sigma_rho_i : curvature sum at inner ring contact
sigma_rho_i = 2.0 / Dw * (2.0 + gamma / (1 - gamma) - Dw / 2.0 / ri)
print("sigma_rho_i = %12.4f" % sigma_rho_i)

# sigma_rho_e : curvature sum at outer ring contact
sigma_rho_e = 2.0 / Dw * (2.0 - gamma / (1 + gamma) - Dw / 2.0 / re)
print("sigma_rho_e = %12.4f" % sigma_rho_e)

# Fi_rho : relative curvature difference at the inner ring contact
Fi_rho = (gamma / (1 - gamma) + Dw / 2.0 / ri) / (2.0 + gamma / (1 - gamma) - Dw / 2.0 / ri)
print("Fi_rho = %12.4f" % Fi_rho)

# Fe_rho : relative curvature difference at the outer ring contact
Fe_rho = (-1.0 * gamma / (1 + gamma) + Dw / 2.0 / re) / (2.0 - gamma / (1 + gamma) - Dw / 2.0 / re)
print("Fe_rho = %12.4f" % Fe_rho)


# K(chi) : Equation 3
def k(chi):
    return special.ellipk(1 - 1.0 / chi**2)


# E(chi) : Equation 4
def e(chi):
    return special.ellipe(1 - 1.0 / chi**2)


# Equation 2
def eq_2_1(chi):
    return 1 - 2.0 / (chi**2 - 1.0) * (k(chi) / e(chi) - 1.0) - Fi_rho


def eq_2_2(chi):
    return 1 - 2.0 / (chi**2 - 1.0) * (k(chi) / e(chi) - 1.0) - Fe_rho


# Solve Equation 2
chi_i = optimize.newton(eq_2_1, 0.1)
chi_e = optimize.newton(eq_2_2, 0.1)
print("chi_i = %12.4f" % chi_i)
print("chi_e = %12.4f" % chi_e)

# cp : spring constant
cp = 1.48 * mat_E / (1 - mat_nu ** 2) * np.power(k(chi_i) * np.power(sigma_rho_i / chi_i ** 2 / e(chi_i), 1 / 3) +
                                                 k(chi_e) * np.power(sigma_rho_e / chi_e ** 2 / e(chi_e), 1 / 3), -1.5)
print("cp = %12.4f" % cp)

# input_data = np.arange(0.01, 20.0, 1.0)
# output_data = Eq_2_2(input_data)
# plt.plot(input_data, output_data)
# plt.show()


# Equation 16
def eq_16(delta_r, delta_a, psi):
    temp_sum = 0.0
    for phi_j in phi:
        delta_j = math.sqrt((A * math.cos(alpha0) + delta_r * math.cos(math.radians(phi_j))) ** 2 +
                            (A * math.sin(alpha0) + delta_a + Ri * math.sin(psi) * math.cos(math.radians(phi_j))) ** 2) - A
        if delta_j < 0.0:
            delta_j = 0.0
        alpha_j = math.atan((A * math.sin(alpha0) + delta_a + Ri * math.sin(psi) * math.cos(math.radians(phi_j))) /
                            (A * math.cos(alpha0) + delta_r * math.cos(math.radians(phi_j))))
        temp_sum = temp_sum + pow(delta_j, 1.5) * math.cos(alpha_j) * math.cos(math.radians(phi_j))
    return Fr - cp * temp_sum


# Equation 17
def eq_17(delta_r, delta_a, psi):
    temp_sum = 0.0
    for phi_j in phi:
        delta_j = math.sqrt((A * math.cos(alpha0) + delta_r * math.cos(math.radians(phi_j))) ** 2 +
                            (A * math.sin(alpha0) + delta_a + Ri * math.sin(psi) * math.cos(math.radians(phi_j))) ** 2) - A
        if delta_j < 0.0:
            delta_j = 0.0
        alpha_j = math.atan((A * math.sin(alpha0) + delta_a + Ri * math.sin(psi) * math.cos(math.radians(phi_j))) /
                            (A * math.cos(alpha0) + delta_r * math.cos(math.radians(phi_j))))
        temp_sum = temp_sum + pow(delta_j, 1.5) * math.sin(alpha_j)
    return Fa - cp * temp_sum


# Equation 18
def eq_18(delta_r, delta_a, psi):
    temp_sum = 0.0
    for phi_j in phi:
        delta_j = math.sqrt((A * math.cos(alpha0) + delta_r * math.cos(math.radians(phi_j))) ** 2 +
                            (A * math.sin(alpha0) + delta_a + Ri * math.sin(psi) * math.cos(math.radians(phi_j))) ** 2) - A
        if delta_j < 0.0:
            delta_j = 0.0
        alpha_j = math.atan((A * math.sin(alpha0) + delta_a + Ri * math.sin(psi) * math.cos(math.radians(phi_j))) /
                            (A * math.cos(alpha0) + delta_r * math.cos(math.radians(phi_j))))
        temp_sum = temp_sum + pow(delta_j, 1.5) * math.sin(alpha_j) * math.cos(math.radians(phi_j))
    return Mz - Dpw / 2.0 * cp * temp_sum


# Solve Equation 16 - 18
def displacement(delta):
    temp_sum_1 = 0.0
    temp_sum_2 = 0.0
    temp_sum_3 = 0.0
    for phi_j in phi:
        delta_j = math.sqrt((A * math.cos(alpha0) + delta[0] * math.cos(math.radians(phi_j))) ** 2 +
                            (A * math.sin(alpha0) + delta[1] + Ri * math.sin(delta[2]) * math.cos(math.radians(phi_j))) ** 2) - A
        if delta_j < 0.0:
            delta_j = 0.0
        alpha_j = math.atan((A * math.sin(alpha0) + delta[1] + Ri * math.sin(delta[2]) * math.cos(math.radians(phi_j))) /
                            (A * math.cos(alpha0) + delta[0] * math.cos(math.radians(phi_j))))
        temp_sum_1 = temp_sum_1 + pow(delta_j, 1.5) * math.cos(alpha_j) * math.cos(math.radians(phi_j))
        temp_sum_2 = temp_sum_2 + pow(delta_j, 1.5) * math.sin(alpha_j)
        temp_sum_3 = temp_sum_3 + pow(delta_j, 1.5) * math.sin(alpha_j) * math.cos(math.radians(phi_j))
    return [Fr - cp * temp_sum_1,
            Fa - cp * temp_sum_2,
            Mz - Dpw / 2.0 * cp * temp_sum_3]


delta0 = np.array([0.0, 0.0, 0.0])
sol = optimize.root(displacement, delta0)

# Results , um
Delta_r = sol.x[0] * 1000
Delta_a = sol.x[1] * 1000
Delta_psi = sol.x[2] * 1000
print("Delta_r = %12.4f" % Delta_r)
print("Delta_a = %12.4f" % Delta_a)
print("Delta_psi = %12.4f" % Delta_psi)
