"""Screening-level kinetics for stage 1b (no tblite needed; unit-tested).

partial_hessian_freqs : harmonic frequencies of a block (partial) Hessian over the moving protons, for any isotope masses
kie_semiclassical     : ZPE-only (Bigeleisen-Mayer-type) H/D KIE from R and TS partial-Hessian frequencies
wigner, bell, eckart  : tunnelling factors (Wigner 1932; Bell parabolic barrier, see Skodje & Truhlar 1981;
                        unsymmetrical Eckart barrier, Eckart 1930 / Johnston & Heicklen 1962)
Partial Hessians over the transferring protons are an approximation (screening only); the DFT stage uses full Hessians.
"""
import numpy as np

KB = 1.380649e-23; H = 6.62607015e-34; C_CM = 2.99792458e10; NA = 6.02214076e23
EV = 1.602176634e-19; AMU = 1.66053906660e-27; KCAL = 4184.0 / NA
M_H, M_D = 1.00782503207, 2.0141017778


def partial_hessian_freqs(hess_ev_A2, masses_amu):
    """hess: (3m x 3m) eV/Angstrom^2 over m atoms; returns wavenumbers in cm^-1 (imaginary as negative), ascending."""
    m = np.repeat(np.asarray(masses_amu, float), 3)
    Hm = hess_ev_A2 / np.sqrt(np.outer(m, m))
    lam = np.linalg.eigvalsh(0.5 * (Hm + Hm.T)) * EV / (1e-20 * AMU)     # s^-2
    return np.sign(lam) * np.sqrt(np.abs(lam)) / (2 * np.pi * C_CM)


def zpe_kcal(wn):
    wn = np.asarray([w for w in wn if w > 0.0])
    return float(0.5 * H * C_CM * wn.sum() / KCAL)


def kie_semiclassical(wR_H, wTS_H, wR_D, wTS_D, T=298.15):
    """k_H/k_D from ZPE differences only (real modes; the TS imaginary mode is excluded)."""
    dH = zpe_kcal(wTS_H) - zpe_kcal(wR_H); dD = zpe_kcal(wTS_D) - zpe_kcal(wR_D)
    return float(np.exp(-(dH - dD) * KCAL / (KB * T)))


def wigner(nu_imag_cm, T=298.15):
    u = H * C_CM * abs(nu_imag_cm) / (KB * T)
    return float(1 + u * u / 24.0)


def bell(nu_imag_cm, T=298.15):
    """Bell's infinite-parabola factor (u/2)/sin(u/2); undefined for u >= 2*pi (returns nan)."""
    u = H * C_CM * abs(nu_imag_cm) / (KB * T)
    return float((u / 2) / np.sin(u / 2)) if u < 2 * np.pi else float("nan")


def eckart(nu_imag_cm, V1_kcal, V2_kcal, T=298.15, npts=4000):
    """Unsymmetrical Eckart transmission factor. V1/V2: forward/reverse barrier heights (kcal/mol, > 0)."""
    if V1_kcal <= 0 or V2_kcal <= 0:
        return float("nan")
    hnu = H * C_CM * abs(nu_imag_cm); kT = KB * T
    V1, V2 = V1_kcal * KCAL, V2_kcal * KCAL
    a1, a2 = 2 * np.pi * V1 / hnu, 2 * np.pi * V2 / hnu
    E0 = max(0.0, V1 - V2)                            # lowest energy (from the reactant) at which products exist
    E = np.linspace(E0, V1 + 40 * kT, npts)
    xi = E / V1
    den = 1 / np.sqrt(a1) + 1 / np.sqrt(a2)
    tpa = 2 * np.sqrt(a1 * xi) / den
    tpb = 2 * np.sqrt(np.abs((xi - 1) * a1 + a2)) / den
    dd = a1 * a2 - 4 * np.pi ** 2 / 16
    # P = 1 - (cosh(a-b) + cosh(d)) / (cosh(a+b) + cosh(d)), evaluated with a common scale to avoid overflow
    am, ap = tpa - tpb, tpa + tpb
    if dd >= 0:
        dv = 2 * np.sqrt(dd)
        M = np.maximum(np.maximum(np.abs(am), ap), dv)
        cdn = np.exp(dv - M) + np.exp(-dv - M)
    else:
        M = np.maximum(np.abs(am), ap)
        cdn = 2 * np.cos(2 * np.sqrt(-dd)) * np.exp(-M)
    num = np.exp(am - M) + np.exp(-am - M) + cdn
    den = np.exp(ap - M) + np.exp(-ap - M) + cdn
    P = 1 - num / den
    integrand = P * np.exp(-(E - V1) / kT) / kT
    return float(np.trapezoid(integrand, E) if hasattr(np, "trapezoid") else np.trapz(integrand, E))
