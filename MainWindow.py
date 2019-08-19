from ISO_TS_16281 import Ui_MainWindow
from PyQt5 import QtWidgets
from BallBearing import *


def geometry_clicked():
    temp1 = int(ui.lineEdit.text())
    temp2 = int(ui.lineEdit_2.text())
    temp3 = float(ui.lineEdit_3.text())
    temp4 = float(ui.lineEdit_4.text())
    temp5 = float(ui.lineEdit_5.text())
    temp6 = float(ui.lineEdit_6.text())
    temp7 = float(ui.lineEdit_7.text())
    temp8 = int(ui.lineEdit_8.text())

    temp_ball_bearing = BallBearing(temp1, temp2, temp3, temp4, temp5, temp6, temp8)
    temp_ball_bearing.geometry(temp7)
    res1 = temp_ball_bearing.s_a
    res2 = temp_ball_bearing.alpha0_deg
    print(res1)
    print(res2)
    print(str(res1))
    print(str(res2))
    str1 = str(res1)[0:6]
    str2 = str(res2)[0:6]

    ui.label_19.setText(str1)
    ui.label_20.setText(str2)


def stiffness_clicked():
    pass


def disp_clicked():
    pass


def rating_clicked():
    pass


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

