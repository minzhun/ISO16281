from ISO_TS_16281 import Ui_MainWindow
from PyQt5 import QtWidgets
from BallBearing import *
from Head import *


def geometry_clicked():

    temp1 = int(ui.lineEdit.text())
    temp2 = int(ui.lineEdit_2.text())
    temp3 = float(ui.lineEdit_3.text())
    temp4 = float(ui.lineEdit_4.text())
    temp5 = float(ui.lineEdit_5.text())
    temp6 = float(ui.lineEdit_6.text())
    temp7 = float(ui.lineEdit_7.text())
    temp8 = ui.comboBox.currentText()
    if temp8 == "Deep Groove":
        temp81 = 1
    else:
        temp81 = 2

    temp_ball_bearing = BallBearing(temp1, temp2, temp3, temp4, temp5, temp6, temp81)
    temp_ball_bearing.check()
    temp_ball_bearing.geometry(temp7)
    res1 = temp_ball_bearing.s_a
    res2 = temp_ball_bearing.alpha0_deg
    str1 = "{:.4f}".format(res1)
    str2 = "{:.4f}".format(res2)

    ui.label_19.setText(str1)
    ui.label_20.setText(str2)

    return temp_ball_bearing


def stiffness_clicked():

    temp9 = str(ui.lineEdit_9.text())
    temp10 = float(ui.lineEdit_10.text())
    temp11 = float(ui.lineEdit_11.text())

    temp_material = Material(temp9, temp10, temp11)
    temp_ball_bearing = geometry_clicked()
    temp_ball_bearing.stiffness(temp_material)

    res3 = temp_ball_bearing.cp
    str3 = "{:.2f}".format(res3)

    ui.label_21.setText(str3)

    return temp_ball_bearing


def disp_clicked():

    temp12 = float(ui.lineEdit_12.text())
    temp13 = float(ui.lineEdit_13.text())
    temp14 = float(ui.lineEdit_14.text())
    temp15 = float(ui.lineEdit_15.text())
    temp16 = float(ui.lineEdit_16.text())
    temp17 = float(ui.lineEdit_17.text())

    temp_load = Load(temp12, temp13, temp14, temp15, temp16, temp17)
    temp_ball_bearing = stiffness_clicked()
    temp_ball_bearing.disp(temp_load)

    res4 = temp_ball_bearing.Delta_r
    res5 = temp_ball_bearing.Delta_a_MASTA
    res6 = temp_ball_bearing.Delta_psi
    str4 = "{:.4f}".format(res4)
    str5 = "{:.4f}".format(res5)
    str6 = "{:.4f}".format(res6)

    ui.label_22.setText(str4)
    ui.label_25.setText(str5)
    ui.label_36.setText(str6)

    return temp_ball_bearing


def rating_clicked():

    temp18 = float(ui.lineEdit_18.text())

    temp_basic_capacity = BasicCapacity(temp18)
    temp_ball_bearing = disp_clicked()
    temp_ball_bearing.rating(temp_basic_capacity)

    res7 = temp_ball_bearing.Qci
    res8 = temp_ball_bearing.Qce
    res9 = temp_ball_bearing.Qei
    res10 = temp_ball_bearing.Qee
    res11 = temp_ball_bearing.L10r
    res12 = temp_ball_bearing.Pref
    str7 = "{:.2f}".format(res7)
    str8 = "{:.2f}".format(res8)
    str9 = "{:.2f}".format(res9)
    str10 = "{:.2f}".format(res10)
    str11 = "{:.2f}".format(res11)
    str12 = "{:.2f}".format(res12)

    ui.label_37.setText(str7)
    ui.label_38.setText(str8)
    ui.label_39.setText(str9)
    ui.label_40.setText(str10)
    ui.label_41.setText(str11)
    ui.label_42.setText(str12)

    return temp_ball_bearing


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()

    # # 注意位置
    # MainWindow.slot1 = slot1

    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)

    ui.pushButton_2.clicked.connect(geometry_clicked)
    ui.pushButton_3.clicked.connect(stiffness_clicked)
    ui.pushButton_4.clicked.connect(disp_clicked)
    ui.pushButton_5.clicked.connect(rating_clicked)

    MainWindow.show()
    sys.exit(app.exec_())

