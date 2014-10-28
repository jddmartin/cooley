"""Tests for the Cooley package.
XX I need to actually add checks that returned results are correct.

To run all tests:
  python test_cooley.py
To run a specific test from the command-line:
  python test_cooley.py test_he2.test

XX add references for cooley and cashion papers and any other tests.
"""

import unittest

import os.path, sys
sys.path.append("..")
from cooley import *

class morse_potential():
    def __init__(self,d,a,re):
        self.d=d
        self.re=re
        self.a=a
    def v(self,r):
        return -self.d+self.d*(1-math.exp(-self.a*(r-self.re)))**2

class Test_cashion_compare(unittest.TestCase):
    def test(self):
        d=605.559
        a=0.988879
        re=2.40873
        class test_potential():
            def __init__(self,d,a,re):
                self.v=morse_potential(d,a,re).v
                self.hbar2_div_2m=1.0
        Potential=test_potential(d,a,re)

        guess=-530.0
        npoints=1000
        rmax=10.0
        rmin=re-2

        results=driver(Potential,
                       guess,
                       rmin,
                       rmax,
                       npoints,
                       1.0)
        # see: Table I of Cashion (XX I should check a few different levels).
        print "Cashion energy: ", results["energy"]


class Test_cooley_compare(unittest.TestCase):
    """Tests from Cooley's paper (Morse potential)."""
    def test(self):
        kg_per_amu=1.660538921e-27
        kg_per_atomic_unit=9.10938291e-31
        reduced_mass_amu=0.5
        # goofy scaling that Cooley uses:
        class test_potential():
            def __init__(self,d,a,re):
                self.v=morse_potential(d,a,re).v
                self.hbar2_div_2m=1.0

        # Cooley's Morse parameters:
        d=188.4355
        a=0.711248
        re=1.9975
        Potential=test_potential(d,a,re)

        guess=-110.8
#        guess=-160.4
        npoints=200
        rmax=10.0
        rmin=0.000

        results=driver(Potential,
                       guess,
                       rmin,
                       rmax,
                       npoints
                   )
        # see: Table 2 of Cooley (XX I should a few different levels, etc...)
        print results["energy"]

class test_he2(unittest.TestCase):
    def test(self):
        root=os.path.expanduser("~")
        # should define he2 potential here.
        sys.path.append(root+"/2014/fall2014/calculations/"
                        "he2_potential_energy/20140927/")
        import he2_potential_energy
        reduced_mass_amu=0.5*4.002602
        test_scaled_potential=ScaledPotential(
            he2_potential_energy.he2_potential_energy,
            "nm",
            "K",
            reduced_mass_amu,
            "AMU",
            0)                                   
        guess=-0.001
        npoints=10000
        rmax=100.0
        rmin=0.01

        results=driver(test_scaled_potential,
                       guess,
                       rmin,
                       rmax,
                       npoints
                   )
        print "He2 binding energy: ", results["energy"]
    
if __name__ == "__main__":
    unittest.main()
