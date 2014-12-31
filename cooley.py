"""A Python implementation of the Cooley technique for solving for the 
vibrational energy levels and wavefunctions of a diatomic molecule.

For a description of the techniques see:
J. K. Cashion, J. Chem. Phys. 39, 1872 (1963),
http://dx.doi.org/10.1063/1.1734545 (and references within).

For examples and tests see the "examples_and_tests" subdirectory.

This is written for readability and simplicity, not efficiency.

Written by J. Martin, Fall 2014, GPLv3.
"""

import math

# base class of exceptions for module:
class Error(Exception):
    pass

class Potential():
    """Helper class for generating objects describing potential and units
    for potential in a form that the integration functions of this module
    can use.

    The integration functions of this module take an object with the
    member function ".v(r)" and attribute ".hbar2_div_2m".  Given
    a potential function accepting distance and returning energy in user 
    specified units, together with the reduced mass of the system, 
    this class allows easy generation of the required object.  A centrifugal
    contribution to the potential may also be added.
    """

    def __init__(self, potential_function, 
                 potential_distance_units,
                 potential_energy_units,
                 reduced_mass, 
                 reduced_mass_units,
                 j):
        """Arguments:
             potential_function:  a function of a single variable representing
                                  internuclear separation that
                                  returns the value of the potential at that
                                  separation.
             potential_distance_units: string representing distance units
                                       used in "potential_function".
                                       e.g. "nm", "m"
             potential_energy_units: string representing energy units
                                     returned by "potential_function".
             reduced_mass: reduced mass in units to be specified by:
             reduced_mass_units: units used to specify reduced mass.
                                 e.g. "AU" (atomic units)
                                      "AMU" (atomic mass units).
        """
        energy_units_in_joules={"eV":1.602176565e-19,
                                "meV":1.602176565e-22,
                                "cm^{-1}":1.986445685e-23,    
                                "K":1.3806488e-23,            
                                "AU": 4.35974434e-18,
                                "J":1.0}        
        distance_units_in_meters={"A":1.0e-10,      
                                  "nm":1.0e-9,       
                                  "m":1.0,         
                                  "AU":0.52917721092e-10,}    
        mass_units_in_kilograms={"AMU":1.660538921e-27,
                                 "AU":9.10938291e-31,
                                 "kg":1.0,}   
        hbar=1.054571726e-34 # J s
        
        # define hbar^2/(2m) in the units of user_energy*user_length^2, where
        # user_energy is the energy unit for the returned value of the 
        # potential function, and user_length is the unit for the argument of the
        # potential function:
        self.hbar2_div_2m=(
            hbar**2/(mass_units_in_kilograms[reduced_mass_units]
                     *reduced_mass)/2.0
            /energy_units_in_joules[potential_energy_units]
            /(distance_units_in_meters[potential_distance_units])**2)
                                      
        self.potential_function=potential_function
        self.j=j

    def v(self,r):
        val=self.potential_function(r)
        if self.j == 0:
            return val
        else:
            return (val
                    +self.hbar2_div_2m*self.j*(self.j+1.0)/r**2)

class InwardsError(Error):
    pass

class EnergyError(Error):
    pass

def integrate(Potential, energy, rstart, rstep, max_points, 
              start_linear_from_origin=False):
    """Perform either an ingoing or outgoing integration (depending
    on sign of "rstep") of the 1d Schrodinger equation using Numerov's
    method (as described by Cashion).  The wavefunction returned
    is *not normalized*, and should be normalized by the user as required.

    Arguments:
      Potential: an object with a member function:   ".v(r)" 
                 and the attribute:                ".hbar2_div_2m" 
                 The member function ".v(r)" should return the value of
                 the potential energy as a function of internuclear distance.
                 The attribute ".hbar2_div_2m" should return hbar^2/(2m) 
                 using the same energy and distance units that ".v(r)" uses
                 for its return value and argument respectively.
                 (The class "Potential" in this  module is
                 a convenient way to make this object from a potential function).
      rstart: starting point for integration 
              (in same distance units as argument of "Potential.v(r)")
      rstep: increment for integration; positive for outgoing integration, 
             negative for ingoing integration.
      max_points: maximum number of points to integrate.
    Returns: a list with the wavefunction.  Note that the first element
             corresponds to rstart, and subsequent elements correspond
             to rstart+rstep, rstart+2.0*rstep, etc...
             i.e. for inwards integration (rstep < 0.0), the wavefunction will
             be listed in order of descending r.

    Notes:
    This function can be used to: 
    1) compute continuum wavefunctions (in which case it will be called once)
     or 
    2) used as part of an iterative procedure to compute a bound state
       energy eigenvalue.  In this case it will be called twice for
       each trial energy: once for an inwards integration and once 
       for an outwards integration.
    In case 1) this function may be called directly by the user, but
    in case 2) it will normally be called by a driver function which
    takes care of matching inwards and outwards integrations and determining
    a new energy eigenvalue.

    See module docstring for reference.
    """

    hbar2_div_2m=Potential.hbar2_div_2m
    h2=rstep*rstep    
    rcurrent=rstart
    rnext=rstart+rstep
    if rstep < 0: 
        # start inwards integration based on WKB (A3 of Cashion):
        psi=[1.0e-6,]
        try:
            psi.append(psi[0]
                       /math.exp(rnext*math.sqrt(Potential.v(rnext)-energy)
                                 /hbar2_div_2m-
                                 rcurrent*math.sqrt(Potential.v(rcurrent)-energy)
                                 /hbar2_div_2m))
        except ValueError as e:
            message=("Energy exceeds Potential.v(rnext) "
                     "or Potential.v(rcurrent): %17g" % energy)
            raise EnergyError(message)
    else:
        # start outwards integration as recommended by Cashion:
        psi=[0,1.0e-6]

        if start_linear_from_origin is True:
            # this option is possibly useful when looking at s-wave
            # states in Coulomb potential.  This was introduced to
            # help solve e^- - H scattering length problem.
            psi=[1.0e-6, 1.0e-6*(rstart+rstep)/rstart]

    for i in range(2,max_points):
        rprev=rcurrent
        rcurrent=rnext
        rnext=rstart+rstep*i

        yprev=(1.0-h2/12.0*(Potential.v(rprev)-energy)/hbar2_div_2m)*psi[-2]
        ycurrent=(1.0-h2/12.0*(Potential.v(rcurrent)-energy)/hbar2_div_2m)*psi[-1]
        ynew=(2.0*ycurrent-yprev+h2*(Potential.v(rcurrent)-energy)
              /hbar2_div_2m*psi[-1])
        psi.append(ynew/(1.0-h2/12.0*(Potential.v(rnext)-energy)/hbar2_div_2m))
        if rstep < 0.0: # inwards integration
            if ynew < ycurrent: # inwards integration stopping criteria
                break
        if abs(psi[-1]) > 1.0e6: # rescale to prevent overflow:
            rescale(psi, abs(psi[-1]))
    else: # we have iterated to end:
        if rstep < 0.0: # inwards integration stopping criteria *not met*:
            # energy is possibly too low.
            raise InwardsError("Inwards integration stopping criteria not met.")
    return psi

def rescale(psi, factor):
    for i, psival in enumerate(psi):
        psi[i]=psival/factor

def update_energy(Potential, energy, rmin, rmax, npoints):
    """Perform a single iteration of energy updating technique described
    by Cashion (originally due to Cooley).
    This function will normally not be called directly by the user.

    For a description of arguments see "integrate".  Returns the updated
    energy (in energy units of "Potential.v(r)") and wavefunction.
    
    See module docstring for reference.
    """

    hbar2_div_2m=Potential.hbar2_div_2m
    h=(rmax-rmin)/(npoints-1)
    h2=h**2

    # integrate in:
    psi_inwards=integrate(Potential, energy, rmax, -h, npoints)
    rescale(psi_inwards,psi_inwards[-1])
    m_index=npoints-len(psi_inwards)

    # integrate out:
    psi_outwards=integrate(Potential, energy, rmin, h, m_index+1)
    rescale(psi_outwards,psi_outwards[-1])

    # splice two solutions together:
    psi=psi_outwards[:-1]+psi_inwards[::-1]

    # calculate correction to energy, based on discontinuity at meeting point 
    # (A4 of Cashion):
    sum_psi2=0.0
    for apsi in psi:
        sum_psi2 += apsi**2

    rm=rmin+h*m_index
    rinner=rm-h
    router=rm+h

    yinner=psi[m_index-1]*(1-h2/12.0*(Potential.v(rinner)-energy)/hbar2_div_2m)
    ym=psi[m_index]*(1-h2/12.0*(Potential.v(rm)-energy)/hbar2_div_2m)
    youter=psi[m_index+1]*(1-h2/12.0*(Potential.v(router)-energy)/hbar2_div_2m)
    
    energy_correction=hbar2_div_2m*(((-yinner+2.0*ym-youter)/h2
                                     +(Potential.v(rm)-energy)
                                     /hbar2_div_2m*psi[m_index])
                                    /sum_psi2)

    new_energy=energy+energy_correction
    return new_energy, psi

        
class MaxIterationsError(Error):
    pass

def find_single_eigen(Potential, 
                      energy_guess, 
                      rmin, rmax, npoints, 
                      tolerance,
                      max_iterations=100,
                      diagnostics=False):
    """Iterative procedure to find a single energy eigenvalue and wavefunction,
    using Cooley's procedure.  For many cases this will be the only function 
    directly called by the user.

    Arguments:
      Potential: an object with a member function:   ".v(r)" 
                 and the attribute:                ".hbar2_div_2m" 
                 The member function ".v(r)" should return the value of
                 the potential energy as a function of internuclear distance.
                 The attribute ".hbar2_div_2m" should return hbar^2/(2m) 
                 using the same energy and distance units that ".v(r)" uses
                 for its return value and argument respectively.
                 (The class "Potential" in this  module is
                 a convenient way to make this object from a potential function).
      energy_guess: initial guess value for energy 
                    (same energy units as the return value of "Potential.v(r)")
      rmin, rmax: integration limits
                  (same distance units as argument of "Potential.v(r)").
      npoints: number of integration points 
               i.e. grid spacing will be: dr=(rmax-rmin)/(npoints-1.0)
      tolerance: the energy eigenvalue is considered found when two
                 consecutive updates give an absolute change in energy
                 less than tolerance.
                 (same energy units as the return value of "Potential.v(r)")
      max_iterations: maximum number of times to obtain an update of
                      the energy.
      diagnostics: whether or not diagnostics information should
                   be obtained.
    Returns:
      A dictionary with various results of the calculation indexed by keys.
      The most important corresponds is the "success_code" value, indicating
      whether or not the calculation was a success (=1) or not.  You should
      *always* check this value.  
    
    See module docstring for reference.
    """

    count=0
    energies=[energy_guess,]
    h=(rmax-rmin)/(npoints-1)
    success_code=0
    while count < max_iterations:
        count += 1
        if diagnostics:
            print "Trying energy: ", energies[-1]
        new_energy, psi = update_energy(Potential, energies[-1], 
                                        rmin, rmax, npoints)
        if diagnostics:
            f=open("diagnostic_wavefunction.dat","w")
            for i,apsi in enumerate(psi):
                f.write("%f %f\n" % (rmin+i*h,apsi))
            f.close()

        energies.append(new_energy)
        if ((abs(energies[-1]-energies[-2]) <= tolerance) and
            ((abs(energies[-2]-energies[-3]) <= tolerance))):
            success_code=1
            break
    else: # "max_iterations" reached without success:
        message=("Maximum iterations exceeded: %d" % max_iterations)
        raise MaxIterationsError(message)
    # normalize wavefunction:
    factor=1.0/math.sqrt(h*sum([apsi*apsi for apsi in psi]))
    normalized_psi=[apsi*factor for apsi in psi]
    return {"energy":energies[-2], # this is energy corresponding to psi
            # not the most recent
            "psi":normalized_psi,
            "rmin":rmin,
            "h":h}
