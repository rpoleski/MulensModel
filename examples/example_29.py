"""
Code that translates VBMicrolensing/RTModel xallarap parameters (called BO there) to MulensModel conventions.

Written with Claude Code
"""
import numpy as np
import matplotlib.pyplot as plt

import MulensModel as mm
import VBMicrolensing


def rtmodel_to_mulensmodel_xallarap(tE, FR, u01, u02, t01, t02, gamma1, gamma2, gammaz,
                                    mass_luminosity_exponent=4.0):
    """
    Transform xallarap orbital parameters from VBM/RTM to MM.
    """
    q = FR ** (1.0 / mass_luminosity_exponent)

    dt = t02 - t01
    t02_eff = t02 + gamma1 * dt * tE
    u02_eff = u02 - gamma2 * dt

    s = np.array([(t01 - t02_eff) / tE, u02_eff - u01, 0.0])
    if abs(gammaz) > 1e-15:
        s[2] = -(s[0]*gamma1 + s[1]*gamma2) / gammaz

    g = np.array([gamma1, gamma2, gammaz])
    A = np.linalg.norm(s)
    g_norm = np.linalg.norm(g)

    if A < 1e-15:
        raise ValueError("Sources coincident (|s| ~ 0).")

    if g_norm < 1e-15:
        raise ValueError("All gamma components zero; no orbital motion.")

    period = 2.0 * np.pi * A / g_norm  # [days]
    L = np.cross(s, g)  # Calculate angular momentum.
    L_norm = np.linalg.norm(L)
    L_sky_norm = np.hypot(L[0], L[1])
    i_rad = np.arccos(np.clip(L[2] / L_norm, -1., 1.)) if L_norm > 1e-15 else 0.0
    Omega = np.array([-L[1], L[0]]) / L_sky_norm if L_sky_norm > 1e-12 else np.array([1., 0.])
    Omega_rad = np.arctan2(Omega[1], Omega[0]) % (2.0 * np.pi)

    phi0_2 = np.arccos(np.clip(np.dot(s[:2], Omega) / A, -1., 1.))
    if s[2] < 0.:
        phi0_2 = -phi0_2
    phi0_1 = (phi0_2 + np.pi) % (2.0 * np.pi)
    a1 = A * q / (1.0 + q)

    return {
        'xi_period': period,
        'xi_semimajor_axis': a1,
        'xi_inclination': np.degrees(i_rad),
        'xi_Omega_node': np.degrees(Omega_rad),
        'xi_argument_of_latitude_reference': np.degrees(phi0_1),
        't_0_xi': t01,
        'q_source': q}


def rtmodel_to_mulensmodel_for_get_trajectory(tE, FR, u01, u02, t01, t02, piN, piE, rho, gamma1, gamma2, gammaz,
                                              mass_luminosity_exponent=4.0,):
    """
    Transform VBM/RTM parameters of BO model to get MM parameters that result in corresponding trajectory.
    Returns parameters *dict* and angle needed for trajectory rotation.
    """
    out = rtmodel_to_mulensmodel_xallarap(tE, FR, u01, u02, t01, t02, gamma1, gamma2, gammaz,
                                          mass_luminosity_exponent=mass_luminosity_exponent)
    q = FR**(1.0 / mass_luminosity_exponent)
    vt = gamma1 * q / (1.0 + q) + 1.0 / tE
    vu = gamma2 * q / (1.0 + q)
    tE_eff = 1.0 / np.sqrt(vt**2 + vu**2)
    alpha_eff_deg = np.degrees(np.arctan2(vu, vt))
    out['xi_Omega_node'] -= alpha_eff_deg
    mm_params = {'t_0': t01, 'u_0': u01, 't_E': tE_eff, 'pi_E_N': piN, 'pi_E_E': piE, 'rho_1': rho, **out}
    return mm_params, alpha_eff_deg


def source1_mulensmodel_get_trajectory(model, alpha_eff_deg, u01, times):
    """
    Get trajectory of source 1 and rotate it by alpha_eff_deg around u01,
    finally, mirror the coordinates to be in VBM/RTModel convention
    """
    trajs = model.get_trajectory(times)
    tau_mm = np.array(trajs[0].x)
    u_mm = np.array(trajs[0].y)

    ca = np.cos(np.radians(alpha_eff_deg))
    sa = np.sin(np.radians(alpha_eff_deg))
    tau_rt = ca * tau_mm - sa * (u_mm - u01)
    u_rt = sa * tau_mm + ca * (u_mm - u01) + u01
    return (-tau_rt, -u_rt)


def get_dict_from_VBM_BO_list(list_VBM):
    """Change VBM/RTM BO list to a dict"""
    out = {'tE': np.exp(list_VBM[0]), 'FR': np.exp(list_VBM[1]), 'u01': list_VBM[2], 'u02': list_VBM[3],
           't01': list_VBM[4] + 2450000., 't02': list_VBM[5] + 2450000., 'rho': np.exp(list_VBM[6]),
           'piN': list_VBM[7], 'piE': list_VBM[8],
           'gamma1': list_VBM[9], 'gamma2': list_VBM[10], 'gammaz': list_VBM[11]}
    return out


if __name__ == '__main__':
    # Values for tests:
    coords = "18:00:00.0 -30:00:00.0"
    p_VBM = [4.605, -0.693, 0.1, 0.3, 10000.0, 9990.0, -6.908, 0.4, -0.3, 0.01, -0.005, 0.008]

    print("VBM/RTModel list:")
    print(p_VBM)
    print("\ndict format:")
    p = get_dict_from_VBM_BO_list(p_VBM)
    print(p)
    t_grid = np.linspace(p['t01'] - 100, p['t01'] + 200, 600)

    (mm_params, alpha_eff_deg) = rtmodel_to_mulensmodel_for_get_trajectory(**p)

    VBM = VBMicrolensing.VBMicrolensing()
    VBM.SetObjectCoordinates(coords)
    (_, tau_rt, u_rt, _, _) = VBM.BinSourceExtLightCurveXallarap(p_VBM, t_grid-2450000.)

    model = mm.Model(mm_params, coords=coords)
    print("\nMM dict:")
    print(model.parameters.parameters)
    (tau_aligned, u_aligned) = source1_mulensmodel_get_trajectory(model, alpha_eff_deg, p['u01'], t_grid)

    plt.plot(tau_rt, u_rt, lw=5, label='RTModel')
    plt.plot(tau_aligned, u_aligned, lw=2, ls='--', label='MulensModel rotated')
    plt.legend()
    plt.xlabel("x [theta_E]")
    plt.ylabel("y [theta_E]")
    plt.show()
