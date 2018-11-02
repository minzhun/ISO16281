from BallBearing import *

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
bearing_1.disp(1000.0, 0.0, 0.0)
bearing_1.rating(23400)
bearing_1.residual()
Err1 = (bearing_1.Delta_r * 1000 - 16.9905) / 16.9905 * 100
Err2 = (bearing_1.L10r - 12008.9256) / 12008.9256 * 100
print("Err 1 = %.4f" % Err1, "%")
print("Err 2 = %.4f" % Err2, "%")

# Example 2 : 深沟球轴承，径向间隙为零，只受轴向力
# Delta : 138.9128 um
# L10r : 2485.2612
# BallBearing ( i, Z, Dw, Dpw, ri, re, alpha, phi0, type )
print("----------------------------------------------------------------")
print("Example 2")
bearing_2 = BallBearing(1, 6, 11.1, 43.5, 5.772, 5.883, 0.0, 0.0, 0)
bearing_2.geometry()
bearing_2.material()
bearing_2.stiffness()
bearing_2.internal_clearance(0.0)
bearing_2.info()
bearing_2.disp(0.0, 1000.0, 0.0)
bearing_2.rating(23400)
bearing_2.residual()
Err1 = (bearing_2.Delta_a_MASTA * 1000 - 138.9128) / 138.9128 * 100
Err2 = (bearing_2.L10r - 2485.2612) / 2485.2612 * 100
print("Err 1 = %.4f" % Err1, "%")
print("Err 2 = %.4f" % Err2, "%")

# Example 3 : 深沟球轴承，径向间隙不为零，只受径向力
# Delta : 23.224 um
# L10r : 9323.8946
# BallBearing ( i, Z, Dw, Dpw, ri, re, alpha, phi0, type )
print("----------------------------------------------------------------")
print("Example 3")
bearing_3 = BallBearing(1, 6, 11.1, 43.5, 5.772, 5.883, 0.0, 0.0, 0)
bearing_3.geometry()
bearing_3.material()
bearing_3.stiffness()
bearing_3.internal_clearance(0.01)
bearing_3.info()
bearing_3.disp(1000.0, 0.0, 0.0)
bearing_3.rating(23400)
bearing_3.residual()
Err1 = (bearing_3.Delta_r * 1000 - 23.224) / 23.224 * 100
Err2 = (bearing_3.L10r - 9323.8946) / 9323.8946 * 100
print("Err 1 = %.4f" % Err1, "%")
print("Err 2 = %.4f" % Err2, "%")

# Example 4 : 深沟球轴承，径向间隙不为零，只受轴向力
# Delta : 153.9587 um
# L10r : 3402.641
# BallBearing ( i, Z, Dw, Dpw, ri, re, alpha, phi0, type )
print("----------------------------------------------------------------")
print("Example 4")
bearing_4 = BallBearing(1, 6, 11.1, 43.5, 5.772, 5.883, 0.0, 0.0, 0)
bearing_4.geometry()
bearing_4.material()
bearing_4.stiffness()
bearing_4.internal_clearance(0.01)
bearing_4.info()
bearing_4.disp(0.0, 1000.0, 0.0)
bearing_4.rating(23400)
bearing_4.residual()
Err1 = (bearing_4.Delta_a_MASTA * 1000 - 153.9587) / 153.9587 * 100
Err2 = (bearing_4.L10r - 3402.641) / 3402.641 * 100
print("Err 1 = %.4f" % Err1, "%")
print("Err 2 = %.4f" % Err2, "%")

# Example 5 : 深沟球轴承，径向间隙不为零，受径向力与轴向力
# Delta : 21.4149 um , 151.1059 um , 2.8266 mrad
# L10r : 1728.1527
# BallBearing ( i, Z, Dw, Dpw, ri, re, alpha, phi0, type )
print("----------------------------------------------------------------")
print("Example 5")
bearing_5 = BallBearing(1, 6, 11.1, 43.5, 5.772, 5.883, 0.0, 0.0, 0)
bearing_5.geometry()
bearing_5.material()
bearing_5.stiffness()
bearing_5.internal_clearance(0.01)
bearing_5.info()
bearing_5.disp(1000.0, 1000.0, 0.0)
bearing_5.rating(23400)
bearing_5.residual()
Err1 = (bearing_5.Delta_r * 1000 - 21.4149) / 21.4149 * 100
Err2 = (bearing_5.Delta_a_MASTA * 1000 - 151.1059) / 151.1059 * 100
Err3 = (bearing_5.L10r - 1728.1527) / 1728.1527 * 100
print("Err 1 = %.4f" % Err1, "%")
print("Err 2 = %.4f" % Err2, "%")
print("Err 3 = %.4f" % Err3, "%")