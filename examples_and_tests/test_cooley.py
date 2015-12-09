"""Tests for the Cooley package.

To run all tests:
  python test_cooley.py
To run a specific test from the command-line:
  python test_cooley.py test_he2.test
"""

import unittest

import os.path, sys
sys.path.append("..")
from cooley import *

class morse_potential():
    """A "Morse" potential."""
    def __init__(self,d,a,re):
        """Potential parameter names are the same as given in Cooley's paper.
        """
        self.d=d
        self.re=re
        self.a=a
    def v(self,r):
        return -self.d+self.d*(1-math.exp(-self.a*(r-self.re)))**2

class Test_cashion_compare(unittest.TestCase):
    """Compare with some results of:
    J. K. Cashion, J. Chem. Phys. 39, 1872 (1963),
    http://dx.doi.org/10.1063/1.1734545
    (It would be better to compare against more results of this paper,
    especially those involving rotation).
    """
    def test(self):
        d=605.559
        a=0.988879
        re=2.40873
        class test_potential():
            def __init__(self,d,a,re):
                self.v=morse_potential(d,a,re).v
                self.hbar2_div_2m=1.0
        Potential=test_potential(d,a,re)

        guess=-580.0
        npoints=1000
        rmax=10.0
        rmin=re-2
        tolerance=1.0e-9

        results=find_single_eigen(Potential,
                                  guess,
                                  rmin,
                                  rmax,
                                  npoints,
                                  tolerance)
        # see Table I of Cashion:
        literature_value=-581.46902
        print
        print "Calculated v=0 energy: ", results["energy"]
        print "Cashion's energy (Table I): ", literature_value
        self.assertAlmostEqual(results["energy"], literature_value, places=4)

class Test_cooley_compare(unittest.TestCase):
    """Compare with some results of:
    Cooley,  Mathematics of Computation, v. 15, pg 363 (1961).
    http://dx.doi.org/10.2307/2003025
    (It would be better to compare against more results of this paper.)
    """
    def test(self):
        class test_potential():
            def __init__(self,d,a,re):
                self.v=morse_potential(d,a,re).v
                self.hbar2_div_2m=1.0
        # Cooley's Morse parameters:
        d=188.4355
        a=0.711248
        re=1.9975
        Potential=test_potential(d,a,re)

        guess=-180.0
        npoints=200
        rmax=10.0
        rmin=0.000
        tolerance=1.0e-9

        results=find_single_eigen(Potential,
                                  guess,
                                  rmin,
                                  rmax,
                                  npoints,
                                  tolerance
                              )
        # see Table 2 of Cooley:
        literature_value=-178.79857
        print
        print "Calculated v=0 energy: ", results["energy"]
        print "Cooley's energy (Table I): ", literature_value
        self.assertAlmostEqual(results["energy"], literature_value, places=4)

class test_he2(unittest.TestCase):
    """Check Helium dimer binding energy."""
    def test(self):
        import helium_dimer_example # see this module for documentation.
        energy=helium_dimer_example.calculate_energy()
        # compare against value in Table I of Janzen and Aziz, 
        # J. Chem. Phys. v. 103, pg 9626 (1995) 
        # http://dx.doi.org/10.1063/1.469978
        literature_value=-1.594e-3
        print
        print "He2 energy (K): ", energy
        print "compared to literature value (K): ", literature_value
        self.assertAlmostEqual(energy, literature_value, places=6)

class test_show(unittest.TestCase):
    """Check SHO ground state energy"""
    def test(self):
        class SHOPotential():
            hbar2_div_2m = 0.5
            def v(self, r):
                return 0.5*r*r

        guess = 0.5
        rmin = -10.0
        rmax = 10.0
        npoints = 1000
        tolerance = 1.0e-6
        results=find_single_eigen(SHOPotential(),
                                  guess,
                                  rmin,
                                  rmax,
                                  npoints,
                                  tolerance,
                                  diagnostics = True
                              )
        print results["energy"]

if __name__ == "__main__":
    unittest.main()
