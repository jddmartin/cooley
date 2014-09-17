"""Implementation of Cooley-Cashion technique for solving for vibrational
energy levels of a diatomic molecule.

See Cashion: http://dx.doi.org/10.1063/1.1734545
"""

import math

class ScaledPotential():
    # helper class to take a potential function and generate scaled
    # version of it suitable for Numerov integration.
    # possibly allow addition of centrigual term
    def __init__(self, potential_function, 
                 potential_distance_units,
                 potential_energy_units,
                 reduced_mass, 
                 reduced_mass_units,
                 j):
        energy_units_in_joules={"eV":1.602176565e-19,
                                "cm^{-1}":1.986445685e-23,    
                                "K":1.3806488e-23,            
                                "AU": 4.35974434e-18,}        
        distance_units_in_meters={"A":1.0e-10,      
                                  "nm":1.0e-9,                
                                  "AU":0.52917721092e-10,}    
        mass_units_in_kilograms={"AMU":1.660538921e-27,
                                 "AU":9.10938291e-31,}   
        hbar=1.054571726e-34 # J s
        
        self.hbar2_div_2m_in_user_units=(
            hbar**2/(mass_units_in_kilograms[reduced_mass_units]
                     *reduced_mass)/2.0
            /energy_units_in_joules[potential_energy_units]
            /(distance_units_in_meters[potential_distance_units])**2)
                                      
        self.potential_function=potential_function
        self.j=j

    def scaled_v(self,r):
        val=self.potential_function(r) / self.hbar2_div_2m_in_user_units
        if self.j == 0:
            return val
        else:
            return val+self.j*(self.j+1.0)/r**2

def integrate(energy, potential, rstart, rstep, max_points):
    """
    rstep < 0.0 specified inwards integration with special stopping criteria
    """

    h2=rstep*rstep

    
    rcurrent=rstart
    rnext=rstart+rstep
    if rstep < 0: 
        # start inwards integration based on WKB (A3 of Cashion):
        psi=[1.0e-6,]
        psi.append(psi[0]
                   /math.exp(rnext*math.sqrt(potential(rnext)-energy)-
                             rcurrent*math.sqrt(potential(rcurrent)-energy)))
    else:
        # start outwards integration as recommended by Cashion:
        psi=[0,1.0e-6]

    for i in range(2,max_points):
        rprev=rcurrent
        rcurrent=rnext
        rnext=rstart+rstep*i

        yprev=(1.0-h2/12.0*(potential(rprev)-energy))*psi[-2]
        ycurrent=(1.0-h2/12.0*(potential(rcurrent)-energy))*psi[-1]
        ynew=2.0*ycurrent-yprev+h2*(potential(rcurrent)-energy)*psi[-1]
        psi.append(ynew/(1.0-h2/12.0*(potential(rnext)-energy)))
        if rstep < 0.0: # inwards integration
            if ynew < ycurrent: # inwards integration stopping criteria
                break
                # XX this may never be satisfied if energy is too low
        if abs(psi[-1]) > 1.0e6: # rescale to prevent overflow:
            rescale(psi, abs(psi[-1]))
    return psi

def rescale(psi, factor):
    for i, psival in enumerate(psi):
        psi[i]=psival/factor

def update_energy(potential, energy, rmin, rmax, npoints):
    h=(rmax-rmin)/(npoints-1)

    # integrate in:
    psi_inwards=integrate(energy, potential, rmax, -h, npoints)
    rescale(psi_inwards,psi_inwards[-1])
    m_index=npoints-len(psi_inwards)

    # integrate out:
    psi_outwards=integrate(energy, potential, rmin, h, m_index+1)
    rescale(psi_outwards,psi_outwards[-1])

    # splice two solutions together:
    psi=psi_outwards[:-1]+psi_inwards[::-1]

    f=open("test_wavefunction.dat","w")
    for i,apsi in enumerate(psi):
        f.write("%f %f\n" % (rmin+i*h,apsi))
    f.close()

    # calculate correction to energy, based on discontinuity at meeting point 
    # (A4 of Cashion):
    sum_psi2=0.0
    for apsi in psi:
        sum_psi2 += apsi**2

    rm=rmin+h*m_index
    rinner=rm-h
    router=rm+h

    h2=h**2
    yinner=psi[m_index-1]*(1-h2/12.0*(potential(rinner)-energy))
    ym=psi[m_index]*(1-h2/12.0*(potential(rm)-energy))
    youter=psi[m_index+1]*(1-h2/12.0*(potential(router)-energy))
    
    energy_correction=(((-yinner+2.0*ym-youter)/h2
                        +(potential(rm)-energy)*psi[m_index])
                       /sum_psi2)

    new_energy=energy+energy_correction
    return new_energy, psi

def driver(potential, energy_guess, rmin, rmax, npoints, corr_fact):
    max_iterations=10
    tolerance=1.0e-9
    count=0
    energies=[energy_guess,]
    while count < max_iterations:
        count += 1
        print "Trying energy: ", energies[-1]*corr_fact
        new_energy, psi = update_energy(potential, energies[-1], 
                                        rmin, rmax, npoints)
        energies.append(new_energy)
        if abs(energies[-1]-energies[-2]) < tolerance:
            success_code=1
            break
    else: # failure, the "while" loop condition was no longer satisfied:
        success_code=0
        print "Failure."
        return success_code
    return success_code
    
class morse_potential():
    def __init__(self,d,a,re):
        self.d=d
        self.re=re
        self.a=a
    def v(self,r):
        return -self.d+self.d*(1-math.exp(-self.a*(r-self.re)))**2

import unittest

class Test_cashion_compare(unittest.TestCase):
    def test(self):
        d=605.559
        a=0.988879
        re=2.40873
        test_potential=morse_potential(d,a,re)
        guess=-530.0
        npoints=1000
        rmax=10.0
        rmin=re-2

        driver(test_potential.v,
               guess,
               rmin,
               rmax,
               npoints,
               1.0)

class Test_cooley_compare(unittest.TestCase):
    def test(self):
        kg_per_amu=1.660538921e-27
        kg_per_atomic_unit=9.10938291e-31
        reduced_mass_amu=0.5
        # goofy scaling that Cooley uses:
        corr_fact=2.0*reduced_mass_amu*kg_per_amu/kg_per_atomic_unit
        # Cooley's Morse parameters:
        d=188.4355/corr_fact
        a=0.711248
        re=1.9975
        test_potential=morse_potential(d,a,re)
        test_scaled_potential=ScaledPotential(test_potential.v,
                                              "AU",
                                              "AU",
                                              reduced_mass_amu,
                                              "AMU",
                                              0)
        guess=-110.8/corr_fact/test_scaled_potential.hbar2_div_2m_in_user_units
        npoints=200
        rmax=10.0
        rmin=0.000

        driver(test_scaled_potential.scaled_v,
               guess,
               rmin,
               rmax,
               npoints,
               (corr_fact*test_scaled_potential.hbar2_div_2m_in_user_units)
           )

class test_he2(unittest.TestCase):
    def test(self):
        import sys
        sys.path.append("../../helium_dimer")
        import he2_potential_energy
        reduced_mass_amu=0.5*4.002602
        test_scaled_potential=ScaledPotential(
            he2_potential_energy.he2_potential_energy,
            "nm",
            "K",
            reduced_mass_amu,
            "AMU",
            0)
        guess=-0.001/test_scaled_potential.hbar2_div_2m_in_user_units
        print guess
        npoints=10000
        rmax=100.0
        rmin=0.01

        driver(test_scaled_potential.scaled_v,
               guess,
               rmin,
               rmax,
               npoints,
               test_scaled_potential.hbar2_div_2m_in_user_units
           )
    

if __name__ == "__main__":
    # to run a specific test from the command-line:
    #   python cooley.py test_he2.test
    unittest.main()
