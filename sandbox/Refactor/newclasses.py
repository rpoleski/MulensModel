class PointLensMagnificationCurve():
    def __init__(self, trajectory=None, parameters=None):
        self.trajectory = trajectory
        if parameters is None:
            self.parameters = {'t_0': None, 'u_0': None, 't_E': None}

    def _check_parameters(self):
        pass

    def get_magnification(self):
        return point_lens_magnification(trajectory)

    def get_d_A_d_params(self):
        return d_A_d_params

class FiniteSourceGould94MagnificationCurve(PointLensMagnificationCurve):

    def __init__(self, **kwargs):
        PointLensMagnificationCurve.__init__(**kwargs)
        self.pspl_magnification = None
        self.B_0 = None

    def _check_parameters(self):
        PointLensMagnificationCurve._check_parameters()
        # Check for rho

    def get_pspl_magnification(self):
        self.pspl_magnification = PointLensMagnificationCurve.get_magnification()
        return self.pspl_magnification

    def get_B_0(self):
        return self.B_0

    def get_magnification(self):
        if self.pspl_magnification is None:
            self.get_pspl_magnification()

        if self.B_0 is None:
            self.get_B_0()

        self.fspl_magnification = self.pspl_magnification() * self.B_0

    def get_d_A_d_params(self):
        PointLensMagnificationCurve.get_d_A_d_params()
        # Modifications for FSPL
        return d_A_d_params

class FiniteSourceYoo04MagnificationCurve(
    FiniteSourceGould94MagnificationCurve):

    def __init__(self, **kwargs):
        FiniteSourceGould94MagnificationCurve.__init__(**kwargs)
        self.B_1 = None

    def _check_parameters(self):
        FiniteSourceGould94MagnificationCurve._check_parameters()
        # Check for gamma

    def get_B_1(self):
        return self.B_1

    def get_magnification(self):
        FiniteSourceGould94MagnificationCurve.get_magnification()
        if self.B_1 is None
            self.get_B_1()

        self.fspl_magnification += self.parameters.gamma * self.B_1

    def get_d_A_d_params(self):
        FiniteSourceGould94MagnificationCurve.get_d_A_d_params()
        # Modifications for B_1 term
        return d_A_d_params

