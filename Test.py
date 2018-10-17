import math
import numpy as np
from scipy import optimize
from scipy import special
import matplotlib.pyplot as plt

# Z : number of rolling elements
Z = 6

# mat_E : modulus of elasticity , MPa
mat_E = 210000.0

# mat_nu : poisson ratio
mat_nu = 0.3

# Dw : Element Diameter , mm
Dw = 7.5

#  Dpw : Pitch Circle Diameter , mm
Dpw = 24.5

# ri : inner groove diameter , mm
ri = 3.9

# re: outer groove diameter , mm
re = 3.975

# alpha : nominal contact angle , degree
alpha = 0.0

# s : radial operating clearance , mm
s = 0.0105

# phi : angular position , degree
phi = np.array([0.0, 60.0, 120.0, 180.0, 240.0, 300.0])

# A : distance between the curvature centres of the raceway groove radii , mm
A = ri + re - Dw
print("A = ", A)

# alpha_0 : initial contact angle , degree
alpha_0 = math.degrees(math.acos(1.0 - s / 2.0 / A))
print("alpha_0 = ", alpha_0)

# Ri : distance between the centre of curvature of the inner raceway groove and the axis of rotation , mm
Ri = 0.5 * Dpw + (ri - 0.5 * Dw) * math.cos(math.radians(alpha_0))
print("Ri = ", Ri)

# gamma : auxiliary parameter
gamma = Dw * math.cos(math.radians(alpha)) / Dpw
print("gamma = ", gamma)

# sigma_rho_i : curvature sum at inner ring contact
sigma_rho_i = 2.0 / Dw * (2.0 + gamma / (1 - gamma) - Dw / 2.0 / ri)
print("sigma_rho_i = ", sigma_rho_i)

# sigma_rho_e : curvature sum at outer ring contact
sigma_rho_e = 2.0 / Dw * (2.0 - gamma / (1 + gamma) - Dw / 2.0 / re)
print("sigma_rho_e = ", sigma_rho_e)

# Fi_rho : relative curvature difference at the inner ring contact
Fi_rho = (gamma / (1 - gamma) + Dw / 2.0 / ri) / (2.0 + gamma / (1 - gamma) - Dw / 2.0 / ri)
print("Fi_rho = ", Fi_rho)

# Fe_rho : relative curvature difference at the outer ring contact
Fe_rho = (-1.0 * gamma / (1 + gamma) + Dw / 2.0 / re) / (2.0 - gamma / (1 + gamma) - Dw / 2.0 / re)
print("Fe_rho = ", Fe_rho)

# K(chi) : Equation 3
def K(chi):
    return special.ellipk(1 - 1.0 / chi**2)

# E(chi) : Equation 4
def E(chi):
    return special.ellipe(1 - 1.0 / chi**2)

# Equation 2
def Eq_2_1(chi):
    return (1 - 2.0 / (chi**2 - 1.0) * (K(chi) / E(chi) - 1.0) - Fi_rho)
def Eq_2_2(chi):
    return (1 - 2.0 / (chi**2 - 1.0) * (K(chi) / E(chi) - 1.0) - Fe_rho)

# Solve Equation 2
chi_i = optimize.newton(Eq_2_1, 0.1)
chi_e = optimize.newton(Eq_2_2, 0.1)
print("chi_i = ", chi_i)
print("chi_e = ", chi_e)

# cp : spring constant
cp = 1.48 * mat_E / (1 - mat_nu ** 2) * (K(chi_i) * pow((sigma_rho_i / chi_i ** 2 / E(chi_i)), 1 / 3) + K(chi_e) * pow((sigma_rho_e / chi_e ** 2 / E(chi_e)), 1 / 3)) ** (-1.5)
print("cp = ", cp)


# input_data = np.arange(0.01, 20.0, 1.0)
# output_data = Eq_2_2(input_data)
# plt.plot(input_data, output_data)
# plt.show()