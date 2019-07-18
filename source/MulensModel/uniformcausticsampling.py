import numpy as np
import math
import matplotlib.pyplot as plt  # XXX probably remove that

import MulensModel as MM


class UniformCausticSampling(object):
    """
    Uniform sampling of a binary lens caustic.
    Note that calculations take some time for given (s, q).
    Keep that in mind, when optimizing your fitting routine.

    Arguments :
        s: *float*
            Separation of the two lens components relative to
            the Einstein ring size.

        q: *float*
            Mass ratio of the two lens components.

        n_points: *int*
            Number of points used for internal integration.
            Default value should work find.

    Instead of standard parameters (*t_0*, *u_0*, *t_E*, *alpha*), here
    we use four other parameters: two epochs of caustic crossing
    (*t_caustic_in*, *t_caustic_out*) and two curvelinear coordinates of
    caustic crossing (*x_caustic_in*, *x_caustic_out*). The curvelinear
    coordinates are defined so that going from 0 to 1 draws all caustics
    for given separation and mass ratio (this is different convention than
    in papers cited below). For a wide topology (i.e., 2 caustics), there
    is a value between 0 and 1 which separates the caustics (:math:`{x_sep}`)
    and a trajectory
    exists only if *x_caustic_in* and *x_caustic_out* correspond to
    the same caustic, i.e., both are smaller than :math:`x_{sep}` or
    both are larger than :math:`x_{sep}`. For a close topology
    (i.e., 3 caustics), there are two such separating values.

    For description of the curvelinear coordinates, see:

    `Cassan A. 2008 A&A 491, 587 "An alternative parameterisation for
    binary-lens caustic-crossing events"
    <https://ui.adsabs.harvard.edu/abs/2008A%26A...491..587C/abstract>`_

    `Cassan A. et al. 2010 A&A 515, 52
    "Bayesian analysis of caustic-crossing microlensing events"
    <https://ui.adsabs.harvard.edu/abs/2010A%26A...515A..52C/abstract>`_

    In order to visualize the curvelinear coordinates,
    you can run a code like that:

    .. code-block:: python

      import matplotlib.pyplot as plt
      sampling = UniformCausticSampling(s=1.1, q=0.3)
      color = np.linspace(0., 1., 200)
      points = [sampling.caustic_point(c) for c in color]
      x = [p.real for p in points]
      y = [p.imag for p in points]
      plt.scatter(x, y, c=color)
      plt.axis('equal')
      plt.colorbar()
      plt.show()

    This will show an intermediate topology. Change *s=1.1* to *s=2.*
    to plot a wide topology, or to *s=0.7* to plot a close topology.

    To be specific, the central caustics are plotted counter-clockwise
    and *x_caustic=0.* corresponds to right-hand point where the caustic
    crosses the X-axis. For a wide topology, the planetary caustic is
    plotted in a similar way. For a close topology, the lower planetary
    caustic is plotted counter-clockwise and the upper planetary caustic
    is symmetric, thus plotted clockwise. For planetary caustics in
    a close topology, the zero-point of *x_caustic* values is defined
    in a very complicated way, however it is a smooth function of
    *s* and *q*.
    """
    def __init__(self, s, q, n_points=10000):
        self._s = s
        self._q = q
        self._n_points = n_points

        self._get_n_caustics()
        self._get_phi()
        self._integrate()
        self._find_inflections_and_correct()

    def _get_n_caustics(self):
        """
        Get number of caustics: 1, 2, or 3.
        """
        limit = (1. + self.q) / (1. + self.q**(1./3.))**3
        if self.s > 1. / math.sqrt(limit):
            self._n_caustics = 2
        elif self.s < math.pow(limit, 0.25):
            self._n_caustics = 3
        else:
            self._n_caustics = 1

    def _get_phi(self):
        """
        Prepare internal variables:
            self._phi - gives all phi values used for integration
            self._d_phi - step between adjacent phi values
        """
        phi_begin = 0.
        phi_end = 2. * np.pi - 1e-14
        if self._n_caustics == 1:
            phi_end = 4. * np.pi - 1e-14
        self._d_phi = (phi_end - phi_begin) / self._n_points
        self._phi = np.linspace(phi_begin, phi_end, self._n_points)

    def _get_indexes_of_inflection_points(self, values):
        """
        Find inflection points in give tabulated function.
        """
        diff_ = values[1:] - values[:-1]
        diff = np.concatenate(([diff_[-2], diff_[-1]], diff_))
        out = []
        for i in range(1, len(diff)-1):
            if diff[i-1] > diff[i] and diff[i+1] > diff[i]:
                parabola = np.polyfit([-1., 0., 1.], diff[i-1:i+2], 2)
                shift = -0.5 * parabola[1] / parabola[0]  # XXX 1) use it
                # XXX 2) if shift > 0.5 or < -0.5 than use other
                # triple to calculate it
                out.append(i)
        return out

    def _zeta(self, z):
        """
        Apply lens equation in complex coordinates and
        shift to center of mass coordinates.
        """
        z_bar = z.conjugate()
        zeta = -z_bar + (1./z + self.q/(z+self.s)) / (1. + self.q)
        zeta -= self.s * self.q / (1. + self.q)
        return zeta

    def _find_nearest_index(self, array, value):
        """
        returns element of array that is closest to value
        """
        idx = (np.abs(array - value)).argmin()
        return array[idx]

    def _critical_curve(self, phi):
        """
        Calculate a point on critical curve - see eq. 6 in Cassan 2008.
        """
        coeffs = [0., 0., 0., 2.*self.s, 1.]
        exp_i_phi = np.exp(1j * phi)
        coeffs[0] = -self.s * self.s * exp_i_phi / (1. + self.q)
        coeffs[1] = -2. * self.s * exp_i_phi / (1. + self.q)
        coeffs[2] = self.s * self.s - exp_i_phi

        roots = np.polynomial.polynomial.polyroots(np.array(coeffs))

        if self._n_caustics == 3:
            if self._critical_curve_previous is not None:
                # This makes sure we're following right branch.
                new_roots = np.array([0.+0.j, 0.+0.j, 0.+0.j, 0.+0.j])
                for (i, v) in enumerate(self._critical_curve_previous):
                    new_roots[i] = self._find_nearest_index(roots, v)
                roots = new_roots
            self._critical_curve_previous = roots

        return roots

    def _dz_dphi(self, z):
        """
        Eq. 11 from Cassan (2008)
        """
        z_plus_d = z + self.s
        z_plus_d_2 = z_plus_d**2
        z_plus_d_3 = z_plus_d_2 * z_plus_d
        q_z_2 = self.q * z * z
        value = (z_plus_d_2 + q_z_2) * z_plus_d * z / (z_plus_d_3 + z * q_z_2)
        return 0.5j * value

    def _dz_bar_dphi(self, z_bar):
        """
        almost the same as eq. 11 from Cassan (2008)
        """
        return -self._dz_dphi(z_bar)

    def _dzeta_dphi(self, z, phi):
        """
        Eq. (9) and (11) from Cassan (2008)
        """
        dz_dphi_ = self._dz_dphi(z)
        dz_bar_dphi_ = self._dz_bar_dphi(np.conjugate(z))
        return dz_dphi_ + np.exp(1j * phi) * dz_bar_dphi_

    def _caustic_and_trajectory(self, zetas, u_0, alpha,
                                sum_use, flip, caustic):
        """
        check if caustic crosses the line defined by 2 points
        """
        cos_a = math.cos(alpha * np.pi / 180.)
        sin_a = math.sin(alpha * np.pi / 180.)

        x_1 = zetas[:-1].real
        y_1 = zetas[:-1].imag
        x_2 = zetas[1:].real
        y_2 = zetas[1:].imag
        dx = x_2 - x_1
        dy = y_2 - y_1
        tau = (x_1 * y_2 - x_2 * y_1 + u_0 * (dy * sin_a + dx * cos_a))
        tau /= dy * cos_a - dx * sin_a
        x_cross = tau * cos_a - u_0 * sin_a

        # Check if crossing point is between x_1 and x_2:
        index = np.where((x_cross - x_1) * (x_cross - x_2) < 0.)[0]

        fraction = (x_cross[index] - x_1[index]) / (x_2[index] - x_1[index])
        sum_ = sum_use[index] * (1. - fraction) + sum_use[index+1] * fraction
        in_caustic = sum_ / sum_use[-1]
        if self._n_caustics == 3 and caustic > 1:
            begin = self._which_caustic[caustic-1]
            end = self._which_caustic[caustic]
        else:
            end = 0.5 * (self._which_caustic[caustic] +
                         self._which_caustic[caustic-1])
            if flip:
                begin = self._which_caustic[caustic]
            else:
                begin = self._which_caustic[caustic-1]
        x_caustic = begin + in_caustic * (end - begin)
        return x_caustic.tolist()

    def _integrate(self):
        """
        Main integration for Cassan (2008) parameterization.
        It sets internal variables:
        - self._z_all
        - self._sum_1
        - self._sum_2
        - self._z_sum_1
        - self._z_sum_2
        - self._z_index_sum_1
        - self._z_index_sum_2
        """
        size = (self._n_points, 4)
        self._z_all = np.zeros(size, dtype=np.complex128)
        self._sum_1 = np.zeros(self._n_points)
        self._z_sum_1 = np.zeros(self._n_points, dtype=np.complex128)
        self._z_index_sum_1 = np.zeros(self._n_points, dtype=int)
        if self._n_caustics > 1:
            self._sum_2 = np.zeros(self._n_points)
            self._z_sum_2 = np.zeros(self._n_points, dtype=np.complex128)
            self._z_index_sum_2 = np.zeros(self._n_points, dtype=int)
        self._critical_curve_previous = None

        for (i, phi) in enumerate(self._phi):
            self._z_all[i] = self._critical_curve(phi)

            if self._n_caustics == 1:
                self._z_index_sum_1[i] = int(phi / np.pi)
                self._z_sum_1[i] = self._z_all[i, self._z_index_sum_1[i]]
                abs_1 = abs(self._dzeta_dphi(self._z_sum_1[i], phi))
                self._sum_1[i] = self._sum_1[i-1] + abs_1 * self._d_phi
            if self._n_caustics > 1:
                if self._n_caustics == 2:
                    self._z_index_sum_2[i] = int(phi / np.pi)
                    self._z_index_sum_1[i] = self._z_index_sum_2[i] + 2
                if self._n_caustics == 3:
                    signs = 1. / np.conjugate(self._z_all[i])
                    signs = (signs - self._z_all[i]).imag
                    if signs[1] * signs[2] >= 0.:
                        args = [
                            self.s, self.q, self._n_points, i, self._z_all[i]]
                        raise ValueError("Critical error: {:}".format(args))
                    if signs[1] < 0.:
                        self._z_index_sum_2[i] = 2
                    else:
                        self._z_index_sum_2[i] = 1
                self._z_sum_1[i] = self._z_all[i, self._z_index_sum_1[i]]
                self._z_sum_2[i] = self._z_all[i, self._z_index_sum_2[i]]
                abs_1 = abs(self._dzeta_dphi(self._z_sum_1[i], phi))
                abs_2 = abs(self._dzeta_dphi(self._z_sum_2[i], phi))
                self._sum_1[i] = self._sum_1[i-1] + abs_1 * self._d_phi
                self._sum_2[i] = self._sum_2[i-1] + abs_2 * self._d_phi

    def _find_inflections_and_correct(self):
        """
        Find inflection points of s(phi). In the case of close configuration
        also correct the phase of planetary caustic - Do we do it in current
        version XXX ??? - yes because of self._which_caustic
        """
        indexes = self._get_indexes_of_inflection_points(self._sum_1)
        value_1 = [float(i)/self._n_points for i in indexes]
        self._inflections_fractions = {1: value_1}
        if self._n_caustics == 1:
            self._which_caustic = np.array([0., 1.])
            return

        cusps_z_1 = [self._z_sum_1[i] for i in indexes]
        if self._n_caustics == 2:
            add = self._z_all[indexes[1], 2]
        else:
            add = self._z_all[indexes[1], 3]
        cusps_z_1 = [add] + cusps_z_1
        cusps_z_1 += [cusps_z_1[0].conjugate()]

        cusps_zeta_1 = [self._zeta(z) for z in cusps_z_1]

        indexes = self._get_indexes_of_inflection_points(self._sum_2)
        value_2 = [float(i)/self._n_points for i in indexes]
        self._inflections_fractions[2] = value_2
        cusps_z_2 = [self._z_sum_2[i] for i in indexes]
        if self._n_caustics == 2:
            cusps_z_2 = [self._z_all[indexes[1], 0]] + cusps_z_2
            cusps_z_2 += [cusps_z_2[0].conjugate()]

        cusps_zeta_2 = [self._zeta(z) for z in cusps_z_2]

        #if self._n_caustics == 3:
            #arg = np.argmin([np.abs(z) for z in cusps_zeta_2])
            #self._shift = self._sum_2[indexes[arg] - 1]
# HERE
# XXX - currently we skip this part because it produces an artifact
# in self.plot_caustic()
        if False:
# XXX - THE PART BELOW WAS NOT CHECKED FOR Y-AXIS CHANGE
        # Only for close configuration there is a need to modify self._sum_2
        # in order to make inflection point closer to 0.
        #  if self._n_caustics == 3:
            arg = np.argmin([np.abs(z) for z in cusps_zeta_2])
            shift = indexes[arg] - 1
            # shift = 0
            # shift = indexes[arg]
            # print("SHIFT", shift)
            d_1 = self._sum_2[-1] - self._sum_2[-2]
            d_2 = self._sum_2[1] - self._sum_2[0]
            value_shift = self._sum_2[-1] + (d_1 + d_2) / 2.
            self._sum_2 = np.roll(self._sum_2, -shift)
            self._z_all_2 = np.roll(self._z_all, -shift, axis=0)
            # XXX self._z_sum_2
            self._z_index_sum_2 = np.roll(self._z_index_sum_2, shift)
            self._sum_2[:-shift] -= value_shift
            self._sum_2 -= self._sum_2[0]

            indexes = self._get_indexes_of_inflection_points(self._sum_2)
            value_2 = [float(i)/self._n_points for i in indexes]
            self._inflections_fractions[2] = value_2
            value_3 = [self._inflections_fractions[2][i] for i in [0, 2, 1]]
            self._inflections_fractions[3] = value_3
            cusps_z_2 = [self._z_sum_2[i] for i in indexes]
            cusps_zeta_2 = [self._zeta(z) for z in cusps_z_2]

        # XXX This seems unused
        if self._n_caustics == 3:
            indexes = [indexes[i] for i in [0, 2, 1]]
            cusps_z_3 = [self._z_sum_2[i].conjugate() for i in indexes]
            cusps_zeta_3 = [self._zeta(z) for z in cusps_z_3]

        length_1 = 2. * self._sum_1[-1]
        lengths_sum = length_1
        lengths = [length_1]
        if self._n_caustics > 1:
            length_2 = self._sum_2[-1]
            if self._n_caustics == 2:
                length_2 *= 2.
            lengths_sum += length_2
            lengths += [length_1 + length_2]
            if self._n_caustics == 3:
                length_3 = length_2
                lengths_sum += length_2
                lengths += [lengths[-1] + length_3]
        self._which_caustic = np.array([0.] + lengths) / lengths_sum
# XXX  what we want to remember from above:
# - length_
# - indexes
# - cusps_z_
# - cusps_zeta_

    def get_standard_parameters(self, x_caustic_in, x_caustic_out,
                                t_caustic_in, t_caustic_out):
        """
        Get standard binary lens parameters (i.e., ``t_0``, ``u_0``, ``t_E``,
        ``alpha``; see
        :py:class:`~MulensModel.modelparameters.ModelParameters`)
        based on provided curvelinear parameters.

        Note that this function quite frequently raises ``ValueError``
        exception. This is because not all
        (``s``, ``q``, ``x_caustic_in``, and ``x_caustic_out``)
        correspond to real trajectories. The returned values are in
        conventions used by :py:class:`~MulensModel.model.Model`.

        Keywords :
            x_caustic_in: *float*
                Curvelinear coordinate of caustic entrance.
                Must be in (0, 1) range.

            x_caustic_out: *float*
                Curvelinear coordinate of caustic exit.
                Must be in (0, 1) range.

            t_caustic_in: *float*
                Epoch of caustic entrance.

            t_caustic_out: *float*
                Epoch of caustic exit.

        Returns :
            parameters: *dict*
                Dictionary with standard binary parameters, i.e, keys are
                ``t_0``, ``u_0``, ``t_E``, and ``alpha``.
        """
        if (x_caustic_in < 0. or x_caustic_in > 1. or
                x_caustic_out < 0. or x_caustic_out > 1. or
                t_caustic_in >= t_caustic_out):
            msg = 'Wrong input in get_standard_parameters(): {:} {:} {:} {:}'
            raise ValueError(msg.format(x_caustic_in, x_caustic_out,
                                        t_caustic_in, t_caustic_out))

        caustic_in = self.which_caustic(x_caustic_in)
        caustic_out = self.which_caustic(x_caustic_out)
        if caustic_in != caustic_out:
            message = (
                "Function get_standard_parameters() got curvelinear caustic " +
                "coordinates on different caustics.\n" +
                "x_caustic_in = {:} is on caustic ".format(x_caustic_in) +
                "{:} and\n".format(caustic_in) +
                "x_caustic_out = {:} is on caustic ".format(x_caustic_out) +
                "{:}".format(caustic_out))
            raise ValueError(message)

        zeta_in = self.caustic_point(x_caustic_in)
        zeta_out = self.caustic_point(x_caustic_out)

        u_0 = zeta_out.real*zeta_in.imag - zeta_out.imag*zeta_in.real
        u_0 /= abs(zeta_out - zeta_in)
        if zeta_out.real == zeta_in.real:
            alpha = 0.5 * np.pi * np.sign(zeta_out.imag - zeta_in.imag)
        else:
            diff = zeta_out.imag - zeta_in.imag
            alpha = np.arctan(diff / (zeta_out.real - zeta_in.real))
            alpha += np.pi * np.heaviside(zeta_in.real-zeta_out.real, .5)
        alpha *= 180. / np.pi
        if alpha < 0.:
            alpha += 360.
        t_E = (t_caustic_out - t_caustic_in) / abs(zeta_out - zeta_in)
        t_0 = ((zeta_out + zeta_in) / (zeta_out - zeta_in)).real
        t_0 *= 0.5 * (t_caustic_in - t_caustic_out)
        t_0 += 0.5 * (t_caustic_out + t_caustic_in)

        return {'t_0': t_0, 'u_0': u_0, 't_E': t_E, 'alpha': alpha}

    def get_x_in_x_out(self, u_0, alpha):  # XXX name of the function?
        """
        Calculate where given trajectory crosses the caustic.

        Parameters :
            u_0: *float*
                The parameter u_0 of source trajectory, i.e., impact parameter.

            alpha: *float*
                Angle defining the source trajectory.

        Returns :
            x_caustic_points: *list* of *float*
                Caustic coordinates of points where given trajectory crosses
                the caustic. The length is between 0 and 6.
        """
        zetas_1 = self._zeta(self._z_sum_1)
        if self._n_caustics > 1:
            zetas_2 = self._zeta(self._z_sum_2)

        points = self._caustic_and_trajectory(
            zetas_1, u_0, alpha, self._sum_1, flip=False, caustic=1)
        points += self._caustic_and_trajectory(
            zetas_1.conjugate(), u_0, alpha, self._sum_1, flip=True, caustic=1)
        if self._n_caustics > 1:
            points += self._caustic_and_trajectory(
                zetas_2, u_0, alpha, self._sum_2, flip=False, caustic=2)
            points += self._caustic_and_trajectory(
                zetas_2.conjugate(), u_0, alpha, self._sum_2, flip=True,
                caustic=self._n_caustics)
        return points

    def get_uniform_sampling(self, n_points, n_min_for_caustic=10,
                             caustic=None):
        """
        Sample uniformly (x_caustic_in, x_caustic_out) space according to
        Jacobian and requirement that x_caustic_in corresponds to caustic
        entrance, and x_caustic_out corresponds to caustic exit. The Jacobian
        is defined by Eq. 23 of:

        `Cassan A. et al. 2010 A&A 515, 52
        "Bayesian analysis of caustic-crossing microlensing events"
        <https://ui.adsabs.harvard.edu/abs/2010A%26A...515A..52C/abstract>`_

        and the above requirement is defined under Eq. 27 of that paper.

        Relative number of points per caustic is not yet specified.
        Points do not repeat. This function is useful for sampling starting
        distribution for model fitting. For example sampling see
        `Cassan et al. (2010)
        <https://ui.adsabs.harvard.edu/abs/2010A%26A...515A..52C/abstract>`_
        bottom panel of Fig. 1.

        Parameters :
            n_points: *int*
                number of points to be returned

            n_min_for_caustic: *int*
                minimum number of points in each caustic

            caustic: *int* or *None*
                Select which caustic will be sampled. *None* means all
                caustics. Can be *1*, *2*, or *3* but has to be
                <= :py:attr:`n_caustics`.

        Returns :
            x_caustic_in: *np.ndarray*
                Randomly drawn entrance points.

            x_caustic_out: *np.ndarray*
                Corresponding randomly drawn exit points.
        """
        if caustic is not None:
            instance = isinstance(caustic, int)
            if not instance or caustic < 1 or caustic > self._n_caustics:
                raise ValueError(
                    'Wrong caustic in get_uniform_sampling(): ' + str(caustic))
        if n_min_for_caustic * self._n_caustics > n_points:
            msg = 'wrong input for {:} caustics: {:} {:}'
            raise ValueError(msg.format(
                self._n_caustics, n_points, n_min_for_caustic))

        increase_factor = 10  # How many more points we will have internally.

        if caustic is not None:
            n_points_for_caustic = [n_points]
            caustics = [caustic]
        else:
            n_points_for_caustic = self._select_n_points(n_points,
                                                         n_min_for_caustic)
            caustics = [i+1 for i in range(self._n_caustics)]

        x_1_out = []
        x_2_out = []
        for (caustic_, n_points_) in zip(caustics, n_points_for_caustic):
            out = [False]
            factor = 1.
            while not out[0]:
                out = self._get_uniform_sampling_one_caustic(
                    caustic_, n_points_, increase_factor*factor)
                if not out[0]:
                    factor *= 1.1 * out[1]
            x_1_out += out[1].tolist()
            x_2_out += out[2].tolist()
        return (np.array(x_1_out), np.array(x_2_out))

    def _get_uniform_sampling_one_caustic(self, caustic, n_points,
                                          increase_factor=10):
        """
        Get uniform sampling for a single caustic.
        """
        min_factor = 3

        n_all = int(n_points * increase_factor) + 1
        begin = self._which_caustic[caustic-1]
        scale = self._which_caustic[caustic] - begin
        x_1 = np.random.rand(n_all) * scale + begin
        x_2 = np.random.rand(n_all) * scale + begin
        jacobian = np.zeros(n_all)
        for (i, (x_1_, x_2_)) in enumerate(zip(x_1, x_2)):
            jacobian[i] = self._jacobian(x_1_, x_2_)
        index = np.where(jacobian > 0.)[0]
        if len(index) < n_points * min_factor:
            return (False, n_points * min_factor / len(index))
        jacobian_masked = jacobian[index]
        probabilities = jacobian_masked / np.sum(jacobian_masked)
        out = np.random.choice(index, size=n_points, replace=False,
                               p=probabilities)
        return (True, x_1[out], x_2[out])

# XXX the function below should be public
    def _jacobian(self, x_caustic_in, x_caustic_out):
        """
        Evaluates Eq. 23 from Cassan et al. (2010) with condition under Eq. 27.

        Parameters :
            x_caustic_in: *float*
                Point of caustic entrance.
            x_caustic_out: *float*
                Point of caustic exit.

        Returns :
            jacobian: *float*
                Value of Jacobian. Returns *0.* if trajectory does not exist.
        """
        check = self._orientation_check(x_caustic_in, x_caustic_out)
        if not check[0]:
            return 0.
        (zeta_in, zeta_out, dzeta_dphi_in, dzeta_dphi_out) = check[1:]
        dzeta = zeta_out - zeta_in
        wedge_in = (dzeta.real * dzeta_dphi_in.imag -
                    dzeta.imag * dzeta_dphi_in.real)
        wedge_out = (dzeta.real * dzeta_dphi_out.imag -
                     dzeta.imag * dzeta_dphi_out.real)
        jacobian = np.abs(wedge_in) * np.abs(wedge_out) / np.abs(dzeta)**4
        return jacobian

    def _select_n_points(self, n_points, n_min_for_caustic):
        """
        divide n_points into caustics
        """
        if self._n_caustics == 1:
            out = [n_points]
        elif self._n_caustics == 2:
            area_1 = (self._which_caustic[1] - self._which_caustic[0])**2
            area_2 = (self._which_caustic[2] - self._which_caustic[1])**2
            fraction = area_1 / (area_1 + area_2)
            n_1 = int(fraction * n_points + 0.5)
            if n_1 < n_min_for_caustic:
                n_1 = n_min_for_caustic
            n_2 = n_points - n_1
            if n_2 < n_min_for_caustic:
                n_2 = n_min_for_caustic
                n_1 = n_points - n_2
            out = [n_1, n_2]
        elif self._n_caustics == 3:
            area_1 = (self._which_caustic[1] - self._which_caustic[0])**2
            area_2 = (self._which_caustic[2] - self._which_caustic[1])**2
            area_3 = area_2
            fraction_2 = area_2 / (area_1 + area_2 + area_3)
            n_2 = int(fraction_2 * n_points + 0.5)
            if n_2 < n_min_for_caustic:
                n_2 = n_min_for_caustic
            n_3 = n_2
            n_1 = n_points - n_2 - n_3
            if n_1 < n_min_for_caustic:
                n_1 = n_min_for_caustic
                n_2 = (n_points - n_1) // 2
                n_3 = n_points - n_1 - n_2
            out = [n_1, n_2, n_3]
        else:
            raise ValueError('strange error: {:}'.format(self._n_caustics))
        return out

    def orientation_check(self, x_caustic_in, x_caustic_out):  # XXX
        """
        **TO DO: this function should have different name because we check
        the other condition i.e., if the same caustics are hit;
        maybe just name it check() or check_x_in_x_out()**

        Check if given (x_caustic_in, x_caustic_out) define an existing
        trajectory. An obvious case, when they don't is when both caustic
        points are on the same fold, but other cases exists.

        Parameters :
            x_caustic_in: *float*
                Coordinate of putative caustic entrance.

            x_caustic_out: *float*
                Coordinate of putative caustic exit.

        Returns :
            check: *bool*
                *True* if input defines a trajectory, *False* if it does not.
        """
        return self._orientation_check(x_caustic_in, x_caustic_out)[0]

    def _orientation_check(self, x_caustic_in, x_caustic_out):
        """
        Check if given parameters define real trajectory.

        Returns a list with first element being True/False
        """
        if self._n_caustics > 1:
            caustic_in = self.which_caustic(x_caustic_in)
            caustic_out = self.which_caustic(x_caustic_out)
            if caustic_in != caustic_out:
                return [False]

        zeta_in = self.caustic_point(x_caustic_in)
        dzeta_dphi_in = self._last_dzeta_dphi
        zeta_out = self.caustic_point(x_caustic_out)
        dzeta_dphi_out = self._last_dzeta_dphi
        # Eq. 27 from Cassan+2010:
        n_t = (zeta_out - zeta_in) / np.abs(zeta_out - zeta_in)
        # Eq. 26 from Cassan+2010:
        n_c_in = 1j * dzeta_dphi_in / np.abs(dzeta_dphi_in)
        n_c_out = 1j * dzeta_dphi_out / np.abs(dzeta_dphi_out)

        condition_in = n_c_in.real * n_t.real + n_c_in.imag * n_t.imag
        condition_out = n_c_out.real * n_t.real + n_c_out.imag * n_t.imag

        if condition_in < 0. and condition_out > 0.:
            return [True, zeta_in, zeta_out, dzeta_dphi_in, dzeta_dphi_out]
        else:
            return [False]

    def _mirror_normalize_to_0_1(self, x, x_min=0., x_max=1.):
        """
        Normalizes input to 0-1 range but in special way,
        which considers the middle point.
        """
        if x < x_min or x > x_max:
            msg = "problem in _mirror_normalize: {:} {:} {:}"
            raise ValueError(msg.format(x, x_min, x_max))
        middle = (x_min + x_max) / 2.
        if x < middle:
            xx = (x - x_min) / (middle - x_min)
            return (xx, False)
        else:
            xx = (x_max - x) / (x_max - middle)
            return (xx, True)

    def caustic_point(self, x_caustic):
        """
        Calculate caustic position corresponding to given x_caustic.

        Keywords :
            x_caustic: *float*
                Curvelinear coordinate of the point considered.
                Has to be in 0-1 range.

        Returns :
            point: *numpy.complex128*
                Caustic point in complex coordinates.
        """
        caustic = self.which_caustic(x_caustic)

        if self._n_caustics < 3 or caustic == 1:
            (fraction_in_caustic, flip) = self._mirror_normalize_to_0_1(
                x_caustic, self._which_caustic[caustic-1],
                self._which_caustic[caustic])
        else:
            in_caustic = x_caustic - self._which_caustic[caustic-1]
            diff = (
                self._which_caustic[caustic] - self._which_caustic[caustic-1])
            fraction_in_caustic = in_caustic / diff
            flip = False

        if caustic == 1:
            sum_use = self._sum_1
            z_use = self._z_sum_1
        else:
            sum_use = self._sum_2
            z_use = self._z_sum_2
        sum_ = fraction_in_caustic * sum_use[-1]
        phi_interp = np.interp([sum_], sum_use, self._phi)[0]
        zeta = self._zeta(np.interp([phi_interp], self._phi, z_use)[0])
        # XXX the calculation of dzeta_dphi should not use index
        # XXX and should only be done if really needed
        index = np.argsort(np.abs(phi_interp-self._phi))[0]
        zeta_1 = self._zeta(z_use[index])
        index_ = index + 1
        if index_ == self._n_points:
            index_ = 0
        zeta_2 = self._zeta(z_use[index_])
        if flip or caustic == 3:
            zeta = zeta.conjugate()
            dzeta = zeta_2.conjugate() - zeta_1.conjugate()
        else:
            dzeta = zeta_1 - zeta_2
        dzeta /= np.abs(dzeta)
        self._last_dzeta_dphi = dzeta
        return zeta

    def which_caustic(self, x_caustic):
        """
        Indicates on which caustic given point is.

        Keywords :
            x_caustic: *float*
                Curvelinear coordinate to be checked

        Returns :
            i_caustic: *int*
                Number indicating the caustic:

                ``1`` - central caustic,

                ``2`` - planetary caustic; for close configuration it is
                the lower of the two planetary caustics,

                ``3`` - upper planetary caustic.
        """
        if x_caustic > 1.:
            raise ValueError('Got x_caustic > 1 : {:}'.format(x_caustic))
        if x_caustic < 0.:
            raise ValueError('Got x_caustic < 0 : {:}'.format(x_caustic))
        if self._n_caustics == 1:
            return 1
        caustic = np.searchsorted(self._which_caustic, x_caustic)
        if caustic == 0:
            if x_caustic == 0.:
                return 1
            else:
                fmt = 'which_caustic() got {:} and internally had {:}'
                raise ValueError(fmt.format(x_caustic, self._which_caustic))
        return caustic

    @property
    def n_caustics(self):
        """
        *int*

        Number of caustics: *1* for resonant topology, *2* for wide topology,
        or *3* for close topology.
        """
        return self._n_caustics

    @property
    def s(self):
        """
        *float*

        separation of the two lens components relative to Einstein ring size
        """
        return self._s

    @property
    def q(self):
        """
        *float*

        Mass ratio.
        """
        return self._q

    def _plot_caustic(self, n_points=200):
        """
        Plot caustic using uniform sampling and color scale
        """
        color = np.linspace(0, 1, n_points)
        points = [self.caustic_point(c) for c in color]
        x = [p.real for p in points]
        y = [p.imag for p in points]
        plt.scatter(x, y, c=color)
        plt.axis('equal')
        plt.colorbar()

    def _plot_full(self, x_caustic_in, x_caustic_out, t_caustic_in,
                   t_caustic_out, n_points=200):
        """
        Plot caustic and trajectory - useful for checks of model parameter
        calculations.
        """
        self._plot_caustic(n_points=n_points)
        params = self.get_standard_parameters(x_caustic_in, x_caustic_out,
                                              t_caustic_in, t_caustic_out)
        params['s'] = self.s
        params['q'] = self.q
        model = MM.Model(params)  # remove import, if you remove this line
        model.plot_trajectory(t_start=t_caustic_in, t_stop=t_caustic_out)
        model.plot_trajectory(np.array([params['t_0']]),
                              marker='X', arrow=False)

        c_in = self.caustic_point(x_caustic_in)
        c_out = self.caustic_point(x_caustic_out)
        plt.scatter(c_in.real, c_in.imag, marker='X', c='pink')
        plt.scatter(c_out.real, c_out.imag, marker='X', c='red')

        txt = "s = {:} q = {:} x_caustic_in = {:} x_caustic_out = {:}"
        plt.title(txt.format(self.s, self.q, x_caustic_in, x_caustic_out))
        plt.tight_layout()


if __name__ == "__main__":
    # caustic = UniformCausticSampling(s=1.1, q=0.1)
    if False:  # Check basic calculation
        params = caustic.get_standard_parameters(0.13, 0.04, 0., 1.)
        caustic.get_x_in_x_out(u_0=params['u_0'], alpha=params['alpha'])

    if True:  # Test get_uniform_sampling()
        for s in [1.1, 2., 0.7]:
            caustic = UniformCausticSampling(s=s, q=0.3)
            caustic._plot_caustic()
            # 1) Counter-clockwise; starting at largest X
            # 2) Counter-clockwise Counter-clockwise; starting at largest X - BOTH
            # 3) Counter-clockwise Counter-clockwise CLOCKWISE; starting at largest X - AND - two not specified
            plt.show()

    if False:  # Get x_caustic for many (u_0, alpha)
        n_points = 10000
        u_0_ = np.random.rand(n_points) * 3.4 - 1.7
        alpha_ = np.random.rand(n_points) * 360.
        for (u_0, alpha) in zip(u_0_, alpha_):
            x_caustic = caustic.get_x_in_x_out(u_0=u_0, alpha=alpha)
            for i in range(len(x_caustic)):
                for j in range(i+1, len(x_caustic)):
                    check = caustic.orientation_check(
                        x_caustic[i], x_caustic[j])
                    if not check:
                        continue
                    print(x_caustic[i], x_caustic[j], len(x_caustic))
                    print(x_caustic[j], x_caustic[i], len(x_caustic))

    if False:  # Test direction_check()
        n_points = 100000
        x_in = np.random.rand(n_points)
        x_out = np.random.rand(n_points)
        for (x_in_, x_out_) in zip(x_in, x_out):
            if caustic.which_caustic(x_in_) != caustic.which_caustic(x_out_):
                continue
            directions = caustic.direction_check(x_in_, x_out_)
            print(x_in_, x_out_, directions[0], directions[1])

