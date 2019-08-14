import math

class Material:

    def __init__(self, mat_name, mat_E, mat_nu):

        self.mat_name = mat_name
        self.mat_E = mat_E
        self.mat_nu = mat_nu


class Load:

    def __init__(self, Fx, Fy, Fz, Mx, My, Mz):

        self.Fx = Fx
        self.Fy = Fy
        self.Fz = Fz
        self.Mx = Mx
        self.My = My
        self.Mz = Mz
        self.Fr = math.sqrt(Fx*Fx+Fy*Fy)
        self.Fa = Fz
        self.Mb = math.sqrt(Mx*Mx+My*My)


class BasicCapacity:

    def __init__(self, C):

        self.C = C