"""
The helium dimer is weakly bound and has a single vibrational level.
Calculate its binding energy and compare it to the literature value.
"""

import math
import os.path
import sys
sys.path.append("..")
import cooley


def he2_potential_energy(r):
    """Potential energy of two He atoms (in units of Kelvin)
    as a function of internuclear separation (in units of nm).

    Computed using the parameters in Table I of Ref. [1]
    (see References below), together with Eq.'s 2-4 of Ref. [2].

    This is the so-called "HFD-B3-FCI1" potential, which is compared against
    many other potentials in Ref. [3] and reputed to be the "best",
    by for example, Ref. [4].

    References:
    [1] Aziz, Janzen and Moldover, Phys. Rev. Lett, 74, 1586 (1995).
        http://dx.doi.org/10.1103/PhysRevLett.74.1586
    [2] Aziz and Slaman, Metrologia 27. 211 219 (1990).
        http://dx.doi.org/10.1088/0026-1394/27/4/005
    [3] Janzen and Aziz, J. Chem. Phys. v. 103, pg 9626 (1995).
        http://dx.doi.org/10.1063/1.469978
    [4] Grisenti et al., Phys. Rev. Lett, v. 85, 2284 (2000).
        http://dx.doi.org/10.1103/PhysRevLett.85.2284
    """
    r_m = 0.2968300
    d = 1.4380000
    a_star = 1.86924404e5
    alpha_star = 10.5717543
    beta_star = -2.07758779
    c6 = 1.35186623
    c8 = 0.41495143
    c10 = 0.17151143
    epsilon = 10.956

    x = r / r_m
    if x < d:
        f = math.exp(-(d / x - 1.0)**2)
    else:
        f = 1.0
    v_star = (a_star * math.exp(-alpha_star * x + beta_star * x**2) -
              f * (c6 / x**6 + c8 / x**8 + c10 / x**10))
    v = epsilon * v_star

    return v


def calculate_energy():
    reduced_mass_amu = 0.5 * 4.002602
    he2_Potential = cooley.Potential(
        he2_potential_energy,
        "nm",
        "K",
        reduced_mass_amu,
        "AMU",
        0)
    guess = -0.001
    npoints = 10000
    rmax = 100.0
    rmin = 0.01
    tolerance = 1.0e-9
    results = cooley.find_single_eigen(he2_Potential,
                                       guess,
                                       rmin,
                                       rmax,
                                       npoints,
                                       tolerance,
                                       diagnostics=False
                                       )
    return results["energy"]


if __name__ == "__main__":
    # literature value from Table I of Janzen and Aziz (see references above).
    literature_value = -1.594e-3
    energy = calculate_energy()
    print("He2 energy (K): ", energy)
    print("compared to literature value (K): ", literature_value)
