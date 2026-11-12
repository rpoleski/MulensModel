"""
Code that translates VBMicrolensing/RTModel xallarap parameters (called BO there) to MulensModel conventions.

TO DO:
- going directly from VBM to MM, not via current p dict
- parallax; also VBM +-2450000
- combine rtmodel_to_mulensmodel_xallarap() with rtmodel_to_mulensmodel_for_get_trajectory()
- docstrings

Written with Claude Code
"""
import numpy as np
import matplotlib.pyplot as plt

import MulensModel as mm
import VBMicrolensing


def rtmodel_to_mulensmodel_xallarap(tE, FR, u01, u02, t01, t02, gamma1, gamma2, gammaz,
                                    mass_luminosity_exponent=4.0):
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


def rtmodel_to_mulensmodel_for_get_trajectory(tE, FR, u01, u02, t01, t02, gamma1, gamma2, gammaz,
                                              mass_luminosity_exponent=4.0,):

    out = rtmodel_to_mulensmodel_xallarap(tE, FR, u01, u02, t01, t02, gamma1, gamma2, gammaz,
                                          mass_luminosity_exponent=mass_luminosity_exponent)
    q = FR**(1.0 / mass_luminosity_exponent)
    vt = gamma1 * q / (1.0 + q) + 1.0 / tE
    vu = gamma2 * q / (1.0 + q)
    tE_eff = 1.0 / np.sqrt(vt**2 + vu**2)
    alpha_eff_deg = np.degrees(np.arctan2(vu, vt))
    out['xi_Omega_node'] -= alpha_eff_deg
    mm_params = {'t_0': t01, 'u_0': u01, 't_E': tE_eff, **out}
    return mm_params, alpha_eff_deg


def source1_mulensmodel_get_trajectory(model, alpha_eff_deg, u01, times):
    trajs = model.get_trajectory(times)
    tau_mm = np.array(trajs[0].x)
    u_mm = np.array(trajs[0].y)

    ca = np.cos(np.radians(alpha_eff_deg))
    sa = np.sin(np.radians(alpha_eff_deg))
    tau_rt = ca * tau_mm - sa * (u_mm - u01)
    u_rt = sa * tau_mm + ca * (u_mm - u01) + u01
    return (-tau_rt, -u_rt)


if __name__ == '__main__':
    p = dict(
        tE=100.0, u01=0.1,  u02=0.3, t01=2460000.0, t02=2459990.0, FR=0.5,
        gamma1=0.01, gamma2=-0.005, gammaz=0.008,
    )
    rho = 0.001

    t_grid = np.linspace(p['t01'] - 100, p['t01'] + 200, 600)

    (mm_params, alpha_eff_deg) = rtmodel_to_mulensmodel_for_get_trajectory(**p)

    p_VBM = [np.log(p['tE']), np.log(p['FR']), p['u01'], p['u02'], p['t01']-2450000., p['t02']-2450000.,
             np.log(rho), 0.0, 0.0, p['gamma1'], p['gamma2'], p['gammaz']]
    VBM = VBMicrolensing.VBMicrolensing()
    VBM.SetObjectCoordinates("18:00:00.0 -30:00:00.0")
    (_, tau_rt, u_rt, _, _) = VBM.BinSourceExtLightCurveXallarap(p_VBM, t_grid-2450000.)

    model = mm.Model(mm_params)
    (tau_aligned, u_aligned) = source1_mulensmodel_get_trajectory(model, alpha_eff_deg, p['u01'], t_grid)

    plt.plot(tau_rt, u_rt, lw=5, label='RTModel')
    plt.plot(tau_aligned, u_aligned, lw=2, ls='--', label='MulensModel rotated')
    plt.legend()
    plt.xlabel("x [theta_E]")
    plt.ylabel("y [theta_E]")
    plt.show()
