from cooley import *

def main():
    import sys
    sys.path.append("../../helium_dimer")
    import he2_potential_energy
    reduced_mass_amu=0.5*4.002602
#    reduced_mass_amu=(3.0*4.0/(3.0+4.0))
    test_scaled_potential=ScaledPotential(
        he2_potential_energy.he2_potential_energy,
        "nm",
        "K",
        reduced_mass_amu,
        "AMU",
        0)
    npoints=10000
    rmax=100.0
    rmin=0.01
    h=(rmax-rmin)/(npoints-1)
    psi=integrate(0.0, test_scaled_potential.scaled_v, rmin, h, 
                  npoints)
    f=open("cont.dat","w")
    for i, apsi in enumerate(psi):
        f.write("%f %f %f\n" % (i*h+rmin,apsi,
                                test_scaled_potential.scaled_v(i*h+rmin)))

if __name__ == "__main__":
    # to run a specific test from the command-line:
    #   python cooley.py test_he2.test
    main()
