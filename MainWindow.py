from ISO_TS_16281 import Ui_MainWindow
from PyQt5 import QtWidgets
from BallBearing import *


def button_clicked():
    temp = float(ui.lineEdit.text())
    ui.label_30.setText(str(temp))


def slot1():
    temp = float(ui.lineEdit.text())
    ui.label_30.setText(str(temp))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()

    MainWindow.slot1 = slot1 ## 注意位置

    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)


    #ui.pushButton_2.clicked.connect(button_clicked)

    MainWindow.show()
    sys.exit(app.exec_())

