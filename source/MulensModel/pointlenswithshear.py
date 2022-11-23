import numpy as np

from MulensModel.pointlens import PointLens
from MulensModel.trajectory import Trajectory


class PointLensWithShear(PointLens):
    """
    Lens of point mass plus shear and convergence, i.e., Chang-Refsdal lens.

    This lens was first considered by
    `Chang and Refsdal (1979; Nature, 282, 561)
    <https://ui.adsabs.harvard.edu/abs/1979Natur.282..561C/abstract>`_.

    Keywords :
        parameters: :py:class:`~MulensModel.modelparameters.ModelParameters`
            Parameters of the model.
    """
    def get_point_source_magnification(self, trajectory):
        """
        Calculate point source magnification for the lens composed of
        a single mass plus a mass sheat.

        Arguments :
            trajectory: :py:class:`~MulensModel.trajectory.Trajectory` object
                Trajectory of the source.

        Returns :
            pspl_magnification: *np.ndarray*
                The point-source--point-lens magnification for each point
                specified by `trajectory`.
        """
        if not isinstance(trajectory, Trajectory):
            raise TypeError("trajectory must be a Trajectory object.")

        shear_G = self.parameters.shear_G
        if shear_G is None:
            shear_G = 0.
        convergence_K = self.parameters.convergence_K

        shear_G_conj = shear_G.conjugate()
        zeta = trajectory.x + trajectory.y * 1j
        zeta_conj = zeta.conjugate()

        temp = [shear_G * shear_G_conj**2 - convergence_K**2 * shear_G_conj +
                2 * convergence_K * shear_G_conj - shear_G_conj]
        coeffs_array = np.stack(
            [[shear_G] * len(zeta),
             2 * shear_G * zeta_conj - zeta * convergence_K + zeta,
             2 * shear_G * shear_G_conj + shear_G * zeta_conj**2 -
             zeta * zeta_conj * convergence_K + zeta * zeta_conj,
             2 * shear_G * shear_G_conj * zeta_conj -
             zeta * convergence_K * shear_G_conj + zeta * shear_G_conj -
             convergence_K**2 * zeta_conj + 2 * convergence_K * zeta_conj -
             zeta_conj,
             temp * len(zeta)], axis=1)

        magnification = []
        const = (1. - convergence_K)**2
        for coeffs in coeffs_array:
            mag = 0.
            roots = np.polynomial.polynomial.polyroots(coeffs)
            for root in roots:
                if root == 0:
                    continue
                root_conj = np.conjugate(root)
                mag += 1. / abs(
                    const -
                    (root_conj**-2 - shear_G) * (root**-2 - shear_G_conj))
            magnification.append(mag)

        magnification = np.array(magnification)
        self._test_maginification_values(magnification, trajectory)

        return magnification

    def _test_maginification_values(self, magnification, trajectory):
        """
        Test if all magnifications are > 1 and raise ValueError if not.
        """
        if np.all(magnification >= 1.):
            return

        error = ('PointLensWithShear.get_point_source_magnification() '
                 'failed for input:\nconvergence_K : {:}\nshear_G : {:}\n\n'
                 'x y mag\n')
        error = error.format(
            self.parameters.convergence_K, self.parameters.shear_G)
        for (mag, x, y) in zip(magnification, trajectory.x, trajectory.y):
            if mag < 1.:
                error += "{:} {:} {:}\n".format(x, y, mag)

        raise ValueError(error)

    def get_point_lens_finite_source_magnification(self, *args, **kwargs):
        """Not implemented for Chang-Refsdal"""
        raise NotImplementedError("not implemented for Chang-Refsdal")

    def get_point_lens_limb_darkening_magnification(self, *args, **kwargs):
        """Not implemented for Chang-Refsdal"""
        raise NotImplementedError("not implemented for Chang-Refsdal")

    def get_point_lens_uniform_integrated_magnification(self, *args, **kwargs):
        """Not implemented for Chang-Refsdal"""
        raise NotImplementedError("not implemented for Chang-Refsdal")

    def get_point_lens_LD_integrated_magnification(self, *args, **kwargs):
        """Not implemented for Chang-Refsdal"""
        raise NotImplementedError("not implemented for Chang-Refsdal")

    def get_point_lens_large_finite_source_magnification(self, *args,
                                                         **kwargs):
        """Not implemented for Chang-Refsdal"""
        raise NotImplementedError("not implemented for Chang-Refsdal")

    def get_point_lens_large_LD_integrated_magnification(self, *args,
                                                         **kwargs):
        """Not implemented for Chang-Refsdal"""
        raise NotImplementedError("not implemented for Chang-Refsdal")
