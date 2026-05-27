import numpy as np

import VBMicrolensing

from MulensModel.pointlens import _AbstractMagnification
from MulensModel.binarylens import _LimbDarkeningForMagnification, _FiniteSource


class _MultipleLensPointSourceMagnification(_AbstractMagnification):
    """
    Equations for calculating point-source--multiple-lens magnification.
    This is a placeholder class to establish the basic methods and attributes
    and over-write methods from
    :py:class:`~MulensModel.pointlens.PointSourcePointLensMagnification`
    that do not apply to binary lenses.

    Arguments :
        trajectory: :py:class:`~MulensModel.trajectory.Trajectory`
            Including trajectory.parameters =
            :py:class:`~MulensModel.modelparameters.ModelParameters`
    """

    def __init__(self, **kwargs):
        super().__init__(trajectory=kwargs['trajectory'])
        # This speeds-up code for np.float input.
        # Can be manually changed to 'numpy'.
        self._solver = 'Skowron_and_Gould_12'
        self._vbm = VBMicrolensing.VBMicrolensing()

        self._source_x = self.trajectory.x
        self._source_y = self.trajectory.y
        self._geometry = self.trajectory.parameters.get_lens_geometry(
            self.trajectory.times)
        if len(self._geometry) == 1:
            self._geometry = self._geometry * len(self._source_x)
        self._zip_kwargs = None

    def get_magnification(self):
        """
        Calculate the magnification

        Parameters : None

        Returns :
            magnification: *np.ndarray*
                The magnification for each point in :py:attr:`~trajectory`.
        """
        zip_args = [self._source_x, self._source_y, self._geometry]

        out = []
        if self._zip_kwargs is None:
            for (x, y, geometry) in zip(*zip_args):
                out.append(self._get_1_magnification(x, y, geometry))
        else:
            zip_args += [self._zip_kwargs]
            for (x, y, geometry, kwargs_) in zip(*zip_args):
                out.append(self._get_1_magnification(
                    x, y, geometry, **kwargs_))
        self._magnification = np.array(out)
        return self._magnification


class MultipleLensPointSourceVBMMagnification(_MultipleLensPointSourceMagnification):
    """
    Equations for calculating point-source--multiple-lens magnification using VBM for point sources.
    Arguments :
        trajectory: :py:class:`~MulensModel.trajectory.Trajectory`
            Including trajectory.parameters =
            :py:class:`~MulensModel.modelparameters.ModelParameters`
    """

    def _get_1_magnification(self, x, y, geometry):
        """
        Calculate 1 magnification using VBM.
        """
        return self._get_1_magnification_point_source(float(x), float(y), float(geometry))

    def _get_1_magnification_point_source(self, x, y, geometry):
        """
        Call VBM to get 1 magnification for point source.
        This function is also called by child classes.
        """
        self._vbm.SetLensGeometry(geometry)
        return self._vbm.MultiMag0(x, y)


class MultipleLensVBMMagnification(_MultipleLensPointSourceMagnification, _LimbDarkeningForMagnification,
                                   _FiniteSource):
    """
    Multiple lens finite source magnification calculated using VBMicrolensing library that implements advanced contour
    integration algorithm presented by
    `Bozza 2010 MNRAS, 408, 2188 <https://ui.adsabs.harvard.edu/abs/2010MNRAS.408.2188B/abstract>`_.
    See also `VBBL website by Valerio Bozza <http://www.fisica.unisa.it/GravitationAstrophysics/VBBinaryLensing.htm>`_.

    For coordinate system convention see
    :py:class:`BinaryLensQuadrupoleMagnification`

    Arguments :
        trajectory: :py:class:`~MulensModel.trajectory.Trajectory`
            Including trajectory.parameters =
            :py:class:`~MulensModel.modelparameters.ModelParameters`

        gamma: *float*
            Linear limb-darkening coefficient in gamma convention.

        u_limb_darkening: *float*
            Linear limb-darkening coefficient in u convention.
            Note that either *gamma* or *u_limb_darkening* can be
            set.  If neither of them is provided then limb
            darkening is ignored.

        accuracy: *float*
            Requested accuracy of the result.

        relative_accuracy: *float*
            Requested relative accuracy of the result.
    """

    def __init__(self, gamma=None, u_limb_darkening=None, accuracy=0.001, relative_accuracy=None, algorithm=None,
                 **kwargs):
        super().__init__(**kwargs)
        self._set_LD_coeffs(u_limb_darkening=u_limb_darkening, gamma=gamma)
        self._set_and_check_rho()
        self._accuracy = self._parse_accuracy(accuracy)
        self._relative_accuracy = self._parse_accuracy(relative_accuracy)

        self._vbm = VBMicrolensing.VBMicrolensing()
        self._set_algorithm(algorithm)
        if self._u_limb_darkening is None:
            self._vbm_function = self._vbm.MultiMag
        else:
            self._vbm_function = self._vbm.MultiMag2

    def _parse_accuracy(self, accuracy):
        """check if None and run float()"""
        if accuracy is None:
            return None

        if accuracy <= 0.:
            raise ValueError("VBM/VBBL requires accuracy >= 0")

        return float(accuracy)

    def _set_algorithm(self, algorithm):
        """Set algorithm for VBM, should be called before SetLensGeometry"""
        if algorithm is None:
            return
        if algorithm == 'Nopoly':
            self._vbm.SetMethod(self._vbm.Method.Nopoly)
        if algorithm == 'Multipoly':
            self._vbm.SetMethod(self._vbm.Method.Multipoly)
        if algorithm == 'Singlepoly':
            self._vbm.SetMethod(self._vbm.Method.Singlepoly)

    def _get_1_magnification(self, x, y, geometry):
        """MultipleLensPointSourceMagnification
        Calculate 1 magnification using VBM.
        """
        args = [float(x), float(y), self._rho]
        if self._u_limb_darkening is None:
            args.append(self._accuracy)
        else:
            self._vbm.a1 = self._u_limb_darkening
            if self._accuracy is not None:
                self._vbm.Tol = self._accuracy
            if self._relative_accuracy is not None:
                self._vbm.RelTol = self._relative_accuracy
        self._vbm.SetLensGeometry(geometry)
        return self._vbm_function(*args)
