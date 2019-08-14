from BallBearing import *
# from RollerBearing import *

material_1 = Material("Steel", 210000.0, 0.3)
#
load_1 = Load(1000.0, 0.0, 0.0, 0.0, 0.0, 0.0)
load_2 = Load(0.0, 0.0, 1000.0, 0.0, 0.0, 0.0)
load_3 = Load(1000.0, 0.0, 1000.0, 0.0, 0.0, 0.0)
load_4 = Load(0.0, 0.0, 0.0, 1000.0, 0.0, 0.0)
load_5 = Load(1000.0, 0.0, 1000.0, 1000.0, 0.0, 0.0)
#
basic_capacity_1 = BasicCapacity(23400)

# Example 1 : 深沟球轴承，径向间隙为零，只受径向力
# Delta : 15.1101 um
# L10r : 13163.019
# BallBearing ( i, Z, Dw, Dpw, alpha_deg, phi0, type )
print("----------------------------------------------------------------")
print("Example 1")
bearing_1 = BallBearing(1, 7, 11.5, 43.5, 0.0, 0.0, 1)
bearing_1.geometry(0.0)
bearing_1.stiffness(material_1)
bearing_1.info()
bearing_1.disp(load_1)
bearing_1.rating(basic_capacity_1)
bearing_1.residual()
Err1 = (bearing_1.Delta_r * 1000 - 15.1101) / 15.1101 * 100
Err2 = (bearing_1.L10r - 13163.019) / 13163.019 * 100
print("Err 1 = %.4f" % Err1, "%")
print("Err 2 = %.4f" % Err2, "%")

# Example 2 : 深沟球轴承，径向间隙为零，只受轴向力
# Delta : 133.5381 um
# L10r : 2222.7318
# BallBearing ( i, Z, Dw, Dpw, alpha_deg, phi0, type )
print("----------------------------------------------------------------")
print("Example 2")
bearing_2 = BallBearing(1, 7, 11.5, 43.5, 0.0, 0.0, 1)
bearing_2.geometry(0.0)
bearing_2.stiffness(material_1)
bearing_2.info()
bearing_2.disp(load_2)
bearing_2.rating(basic_capacity_1)
bearing_2.residual()
Err1 = (bearing_2.Delta_a_MASTA * 1000 - 133.5381) / 133.5381 * 100
Err2 = (bearing_2.L10r - 2222.7318) / 2222.7318 * 100
print("Err 1 = %.4f" % Err1, "%")
print("Err 2 = %.4f" % Err2, "%")

# Example 3 : 深沟球轴承，径向间隙不为零，只受径向力
# Delta : 21.223 um
# L10r : 10795.7948
# BallBearing ( i, Z, Dw, Dpw, alpha_deg, phi0, type )
print("----------------------------------------------------------------")
print("Example 3")
bearing_3 = BallBearing(1, 7, 11.5, 43.5, 0.0, 0.0, 1)
bearing_3.geometry(0.01)
bearing_3.stiffness(material_1)
bearing_3.info()
bearing_3.disp(load_1)
bearing_3.rating(basic_capacity_1)
bearing_3.residual()
Err1 = (bearing_3.Delta_r * 1000 - 21.223) / 21.223 * 100
Err2 = (bearing_3.L10r - 10795.7948) / 10795.7948 * 100
print("Err 1 = %.4f" % Err1, "%")
print("Err 2 = %.4f" % Err2, "%")

# Example 4 : 深沟球轴承，径向间隙不为零，只受轴向力
# Delta : 149.1717 um
# L10r : 3115.8781
# BallBearing ( i, Z, Dw, Dpw, alpha_deg, phi0, type )
print("----------------------------------------------------------------")
print("Example 4")
bearing_4 = BallBearing(1, 7, 11.5, 43.5, 0.0, 0.0, 1)
bearing_4.geometry(0.01)
bearing_4.stiffness(material_1)
bearing_4.info()
bearing_4.disp(load_2)
bearing_4.rating(basic_capacity_1)
bearing_4.residual()
Err1 = (bearing_4.Delta_a_MASTA * 1000 - 149.1717) / 149.1717 * 100
Err2 = (bearing_4.L10r - 3115.8781) / 3115.8781 * 100
print("Err 1 = %.4f" % Err1, "%")
print("Err 2 = %.4f" % Err2, "%")

# Example 5 : 深沟球轴承，径向间隙不为零，受径向力与轴向力
# Delta : 19.5923 um , 146.4603 um , 2.6872 mrad
# L10r : 1630.5367
# BallBearing ( i, Z, Dw, Dpw, alpha_deg, phi0, type )
print("----------------------------------------------------------------")
print("Example 5")
bearing_5 = BallBearing(1, 7, 11.5, 43.5, 0.0, 0.0, 1)
bearing_5.geometry(0.01)
bearing_5.stiffness(material_1)
bearing_5.info()
bearing_5.disp(load_3)
bearing_5.rating(basic_capacity_1)
bearing_5.residual()
Err1 = (bearing_5.Delta_r * 1000 - 19.5923) / 19.5923 * 100
Err2 = (bearing_5.Delta_a_MASTA * 1000 - 146.4603) / 146.4603 * 100
Err3 = (bearing_5.L10r - 1630.5367) / 1630.5367 * 100
print("Err 1 = %.4f" % Err1, "%")
print("Err 2 = %.4f" % Err2, "%")
print("Err 3 = %.4f" % Err3, "%")

# Example 6 : 深沟球轴承，径向间隙为零，受弯矩
# Delta : 3.6957 mrad
# L10r : 785144.7521
print("----------------------------------------------------------------")
print("Example 6")
bearing_6 = BallBearing(1, 7, 11.5, 43.5, 0.0, 0.0, 1)
bearing_6.geometry(0.0)
bearing_6.stiffness(material_1)
bearing_6.info()
bearing_6.disp(load_4)
bearing_6.rating(basic_capacity_1)
bearing_6.residual()
Err1 = (bearing_6.Delta_psi * 1000 - 3.6957) / 3.6957 * 100
Err2 = (bearing_6.L10r - 785144.7521) / 785144.7521 * 100
print("Err 1 = %.4f" % Err1, "%")
print("Err 2 = %.4f" % Err2, "%")

# Example 7 : 深沟球轴承，径向间隙不为零，受径向力，受轴向力，受弯矩
# Delta : 19.7144 um, 146.3084 um, 2.7239 mrad
# L10r : 1627.7479
print("----------------------------------------------------------------")
print("Example 7")
bearing_7 = BallBearing(1, 7, 11.5, 43.5, 0.0, 0.0, 1)
bearing_7.geometry(0.01)
bearing_7.stiffness(material_1)
bearing_7.info()
bearing_7.disp(load_5)
bearing_7.rating(basic_capacity_1)
bearing_7.residual()
Err1 = (bearing_7.Delta_r * 1000 - 19.7144) / 19.7144 * 100
Err2 = (bearing_7.Delta_a_MASTA * 1000 - 146.3084) / 146.3084 * 100
Err3 = (abs(bearing_7.Delta_psi) * 1000 - 2.7239) / 2.7239 * 100
Err4 = (bearing_7.L10r - 1627.7479) / 1627.7479 * 100
print("Err 1 = %.4f" % Err1, "%")
print("Err 2 = %.4f" % Err2, "%")
print("Err 3 = %.4f" % Err3, "%")
print("Err 4 = %.4f" % Err4, "%")