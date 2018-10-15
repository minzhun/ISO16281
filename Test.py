import math
from scipy import optimize
from scipy import special
import numpy as np
import matplotlib.pyplot as plt

# E: modulus of elasticity, MPa
mat_E = 210000.0

# nu_E: poisson ratio
nu_E = 0.3

# Dw: Element Diameter, mm
Dw = 7.5

#  Dpw: Pitch Circle Diameter, mm
Dpw = 24.5

# ri: Inner Groove Diameter, mm
ri = 3.9

# re: Outer Groove Diameter, mm
re = 3.975

# alpha: Nominal Contact Angle
alpha = 0.0

# gamma: auxiliary parameter
gamma = Dw * math.cos(alpha) / Dpw

# sigma_rho_i: curvature sum at inner ring contact
sigma_rho_i = 2.0 / Dw * (2.0 + gamma / (1 - gamma) - Dw / 2.0 / ri)

# sigma_rho_e: curvature sum at outer ring contact
sigma_rho_e = 2.0 / Dw * (2.0 - gamma / (1 + gamma) - Dw / 2.0 / re)

# Fi_rho: relative curvature difference at the inner ring contact
Fi_rho = (gamma / (1 - gamma) + Dw / 2.0 / ri) / (2.0 + gamma / (1 - gamma) - Dw / 2.0 / ri)

# Fe_rho: relative curvature difference at the outer ring contact
Fe_rho = (-1.0 * gamma / (1 + gamma) + Dw / 2.0 / re) / (2.0 - gamma / (1 + gamma) - Dw / 2.0 / re)

# K(chi): Equation 3
def K(chi):
    return special.ellipk(1 - 1.0 / chi**2)

# E(chi): Equation 4
def E(chi):
    return special.ellipe(1 - 1.0 / chi**2)

# Equation 2
def Eq_2_1(chi):
    return (1 - 2.0 / (chi**2 - 1.0) * (K(chi) / E(chi) - 1.0) - Fi_rho)
def Eq_2_2(chi):
    return (1 - 2.0 / (chi**2 - 1.0) * (K(chi) / E(chi) - 1.0) - Fe_rho)

# Solve Equation 2
chi_i = optimize.newton(Eq_2_1, 5.0)
chi_e = optimize.newton(Eq_2_2, 5.0)

print("chi = ", chi_i, chi_e)

# cp: spring constant
cp = 1.48 * mat_E / (1 - nu_E**2) * (K(chi_i) * pow((sigma_rho_i / chi_i**2 / E(chi_i)), 1/3) + K(chi_e) * pow((sigma_rho_e / chi_e**2 / E(chi_e)),1/3))**(-1.5)

# input_data = np.arange(0.01, 20.0, 1.0)
# output_data = Eq_2_2(input_data)
# plt.plot(input_data, output_data)
# plt.show()