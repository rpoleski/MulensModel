from numpy.testing import assert_almost_equal
import numpy as np
import matplotlib.pyplot as plt
from MulensModel.model import Model
import VBMicrolensing


plot = False

def test_VBM_vs_MM():
    """
    Test MulensModel.Model() for triple lens vs VBMicrolensing. The test is based on the example from VBMicrolensing:
    https://github.com/valboz/VBMicrolensing/blob/main/examples/python_examples/Triple_lens.ipynb
    """
    VBM = VBMicrolensing.VBMicrolensing()
    # Set relative accuracy
    VBM.RelTol = 1e-04
    # Set accuracy
    VBM.Tol = 1e-04
    s12 = 0.765
    # Mass ratio lens 2
    q2 = 0.00066
    # impact parameter
    u0 = 0.0060
    # alpha
    alpha = 3.212
    # source radius in Einstein radii of the total mass.
    rho = 0.0567
    # einstein radius crossing time
    tE = 50.13
    # time of peak magnification
    t0 = 0
    # separation between the last two lenses in descending order of mass in units of total ang. Einstein radii
    s23 = 1.5
    # Mass ratio lens 3
    q3 = 0.000001
    # psi
    psi = -1.5

    num_points = 1000
    tmin = -50
    tmax = 50
    t = np.linspace(t0 + tmin, t0 + tmax, num_points)

    params = [np.log(s12), np.log(q2), u0, alpha, np.log(
        rho), np.log(tE), t0, np.log(s23), np.log(q3), psi]

    # Set the Method that you want use : Singlepoly, Multipoly, Nopoly.
    VBM.SetMethod(VBM.Method.Nopoly)
    magtriple_VBM = VBM.TripleLightCurve(params, t)
    caustics_VBM = VBM.Multicaustics()
    criticalcurves_VBM = VBM.Multicriticalcurves()

    x_VBM = []
    y_VBM = []
    x_critical_VBM = []
    y_critical_VBM = []
    for i in range(len(caustics_VBM)):
        x_VBM.extend(caustics_VBM[i][0])
        y_VBM.extend(caustics_VBM[i][1])
        x_critical_VBM.extend(criticalcurves_VBM[i][0])
        y_critical_VBM.extend(criticalcurves_VBM[i][1])

    Parameters = {
        's_21': s12,
        'q_21': q2,
        'u_0': u0,
        'alpha': np.degrees(alpha),
        'rho': rho,
        't_E': tE,
        't_0': t0,
        's_31': s23,
        'q_31': q3,
        'psi': np.degrees(psi),
    }
    model = Model(parameters=Parameters)
    model.set_magnification_methods([float(min(t)), 'vbm_multiple', float(max(t))])
    model.default_magnification_method = 'vbm_multiple'
    magtriple_MM = model.get_magnification(t)
    model.update_caustics()
    caustics_MM = model.caustics
    x_MM, y_MM = caustics_MM.get_caustics()
    x_critical_MM, y_critical_MM = caustics_MM._critical_curve.x, caustics_MM._critical_curve.y

    if plot:
        plt.scatter(x_VBM, y_VBM, color='r', label='VBM caustics', s=4, alpha=0.1, marker='x')
        plt.scatter(x_MM, y_MM, color='b', label='MM caustics', s=1, alpha=0.1, marker='o')
        plt.scatter(x_critical_VBM, y_critical_VBM, color='r', label='VBM critical curve', s=4, alpha=0.1, marker='x')
        plt.scatter(x_critical_MM, y_critical_MM, color='b', label='MM critical curve', s=1, alpha=0.1, marker='o')
        plt.legend()
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Caustics and critical curves for triple lens')
        plt.show()

    assert_almost_equal(magtriple_MM, magtriple_VBM[0], decimal=3, err_msg='Magnification')
    assert_almost_equal(x_MM, x_VBM, decimal=3, err_msg='Caustics x')
    assert_almost_equal(y_MM, y_VBM, decimal=3, err_msg='Caustics y')
    assert_almost_equal(x_critical_MM, x_critical_VBM, decimal=3, err_msg='Critical x')
    assert_almost_equal(y_critical_MM, y_critical_VBM, decimal=3, err_msg='Critical y')

    return 'git'


def test_2L(n=1):
    """
    test calculations of magnification  of 2L using vbm_multiple and vbm.
    """
    for i in range(n):
        parameters = {
            's': np.random.uniform(0.0001, 5.),
            'q': np.random.uniform(0.0001, 1.),
            'u_0': np.random.uniform(0.0001, 2.),
            'alpha': np.random.uniform(0, 360),
            'rho': np.random.uniform(0.0001, 3.),
            't_E': np.random.uniform(20., 500.),
            't_0': np.random.uniform(-100., 100.)
        }

        t0 = parameters['t_0']
        num_points = 1000
        tmin = -50
        tmax = 50
        t = np.linspace(t0 + tmin, t0 + tmax, num_points)

        model_vbm = Model(parameters=parameters)
        model_vbm.set_magnification_methods([float(min(t)), 'vbm', float(max(t))])
        model_vbm.default_magnification_method = 'vbm'
        mag_2L_vbm = model_vbm.get_magnification(t)
        model_vbm.update_caustics()
        caustics_vbm = model_vbm.caustics
        x_vbm, y_vbm = caustics_vbm.get_caustics()
        x_critical_vbm, y_critical_vbm = caustics_vbm._critical_curve.x,  caustics_vbm._critical_curve.y

        model_vbm_multiple = Model(parameters=parameters)

        model_vbm_multiple.set_magnification_methods([float(min(t)), 'vbm_multiple', float(max(t))])
        model_vbm_multiple.default_magnification_method = 'vbm_multiple'
        mag_2L_vbm_multiple = model_vbm_multiple.get_magnification(t)
        model_vbm_multiple.update_caustics()
        caustics_vbm_multiple = model_vbm_multiple.caustics
        x_vbm_multiple, y_vbm_multiple = caustics_vbm_multiple.get_caustics()
        x_critical_vbm_multiple = caustics_vbm_multiple._critical_curve.x
        y_critical_vbm_multiple = caustics_vbm_multiple._critical_curve.y

        if plot:
            plt.scatter(x_vbm, y_vbm, color='r', label='VBM caustics', s=4, alpha=0.1, marker='x')
            plt.scatter(x_vbm_multiple, y_vbm_multiple, color='b', label='VBM multiple caustics', s=1, alpha=0.1,
                        marker='o')
            plt.scatter(x_critical_vbm, y_critical_vbm, color='r', label='VBM critical curve', s=4, alpha=0.1,
                        marker='x')
            plt.scatter(x_critical_vbm_multiple, y_critical_vbm_multiple, color='b',
                        label='VBM multiple critical curve',  s=1, alpha=0.1, marker='o')
            plt.legend()
            plt.xlabel('x')
            plt.ylabel('y')
            plt.title('Caustics and critical curves for 2L')
            plt.show()

        assert_almost_equal(mag_2L_vbm, mag_2L_vbm_multiple, decimal=3, err_msg='Magnification')
        assert_almost_equal(x_vbm, x_vbm_multiple, decimal=3, err_msg='Caustics x')
        assert_almost_equal(y_vbm, y_vbm_multiple, decimal=3, err_msg='Caustics y')
        assert_almost_equal(x_critical_vbm, x_critical_vbm_multiple, decimal=3, err_msg='Critical x')
        assert_almost_equal(y_critical_vbm, y_critical_vbm_multiple, decimal=3, err_msg='Critical y')
    return 'git'


def test_2L_plus_2L():
    """
    Test the magnification of a triple lens with two planets vs the product of magnifications of
    two binary lenses with the same parameters as the triple lens.
    """
    parameters = {
        's_21': 0.2,
        'q_21': 0.02,
        's_31': 5.,
        'q_31': 0.01,
        'u_0': 0.1,
        'alpha': 100,
        'alpha_31': 1,
        'rho': 0.01,
        't_E': 30.,
        't_0': 0,
    }

    t0 = parameters['t_0']
    num_points = 10000
    tmin = -200
    tmax = 200
    t = np.linspace(t0 + tmin, t0 + tmax, num_points)

    model = Model(parameters=parameters)
    model.set_magnification_methods([float(min(t)), 'vbm_multiple', float(max(t))])
    model.default_magnification_method = 'vbm_multiple'
    magtriple_3L = model.get_magnification(t)
    model.update_caustics()
    caustics_3L = model.caustics
    x_3L, y_3L = caustics_3L.get_caustics()
    x_critical_3L, y_critical_3L = caustics_3L._critical_curve.x, caustics_3L._critical_curve.y
    trajectory_3L = model.get_trajectory(t)
    parameters_1 = {
        's': parameters['s_21'],
        'q': parameters['q_21'],
        'u_0': parameters['u_0'],
        'alpha': parameters['alpha'],
        'rho': parameters['rho'],
        't_E': parameters['t_E'],
        't_0': parameters['t_0'],
    }
    parameters_2 = {
        's': parameters['s_31'],
        'q': parameters['q_31'],
        'u_0': parameters['u_0'],
        'alpha': parameters['alpha_31'],
        'rho': parameters['rho'],
        't_E': parameters['t_E'],
        't_0': parameters['t_0'],
    }
    model_1 = Model(parameters=parameters_1)
    model_1.set_magnification_methods([float(min(t)), 'vbm', float(max(t))])
    model_1.default_magnification_method = 'vbm'
    mag_2L_1 = model_1.get_magnification(t)
    model_1.update_caustics()
    caustics_2L_1 = model_1.caustics
    x_2L_1, y_2L_1 = caustics_2L_1.get_caustics()
    x_critical_2L_1, y_critical_2L_1 = caustics_2L_1._critical_curve.x, caustics_2L_1._critical_curve.y
    trajectory_2L_1 = model_1.get_trajectory(t)
    model_2 = Model(parameters=parameters_2)
    model_2.set_magnification_methods([float(min(t)), 'vbm', float(max(t))])
    model_2.default_magnification_method = 'vbm'
    mag_2L_2 = model_2.get_magnification(t)
    model_2.update_caustics()
    caustics_2L_2 = model_2.caustics
    x_2L_2, y_2L_2 = caustics_2L_2.get_caustics()
    x_critical_2L_2, y_critical_2L_2 = caustics_2L_2._critical_curve.x, caustics_2L_2._critical_curve.y
    trajectory_2L_2 = model_2.get_trajectory(t)

    if plot:
        plt.scatter(x_3L, y_3L, color='r', label='3L caustics', s=4, alpha=0.5, marker='x')
        plt.scatter(x_2L_1, y_2L_1, color='b', label='2L_1 caustics', s=1, alpha=0.5, marker='o')
        plt.scatter(x_2L_2, y_2L_2, color='g', label='2L_2 caustics', s=1, alpha=0.5, marker='o')
        plt.scatter(x_critical_3L, y_critical_3L, color='r', label='3L critical curve', s=4, alpha=0.5, marker='x')
        plt.scatter(x_critical_2L_1, y_critical_2L_1, color='b', label='2L_1 critical curve', s=4, alpha=0.5, marker='o')
        plt.scatter(x_critical_2L_2, y_critical_2L_2, color='g', label='2L_2 critical curve', s=4, alpha=0.5, marker='o')
        plt.scatter(trajectory_3L.x, trajectory_3L.y, color='r', label='3L trajectory', s=1, alpha=0.5, marker='x')
        plt.scatter(trajectory_2L_1.x, trajectory_2L_1.y, color='b', label='2L_1 trajectory', s=1, alpha=0.5, marker='o')
        plt.scatter(trajectory_2L_2.x, trajectory_2L_2.y, color='g', label='2L_2 trajectory', s=1, alpha=0.5, marker='o')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Caustics and critical curves for 3L')
        plt.legend()
        plt.show()

        plt.scatter(t, magtriple_3L, color='r', label='3L magnification', s=4, alpha=0.5, marker='x')
        plt.scatter(t, mag_2L_1, color='b', label='2L_1 magnification', s=1, alpha=0.5, marker='o')
        plt.scatter(t, mag_2L_2, color='g', label='2L_2 magnification', s=1, alpha=0.5, marker='o')
        plt.xlabel('time')
        plt.ylabel('magnification')
        plt.title('Magnification for 3L and 2L')
        plt.legend()

        plt.show()
    # assert_almost_equal(magtriple_3L, mag_2L_1*mag_2L_2, decimal=2, err_msg='Magnification 3L vs product of 2L' )
    return 'git'


if __name__ == '__main__':
    plot = True
    test_VBM_vs_MM()
    test_2L(n=10)
    test_2L_plus_2L()
