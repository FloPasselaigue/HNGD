#ifndef Precipitation_hpp
#define Precipitation_hpp

/**
    This is a superclass for hydride nucleation and growth.
    It inherits from Mechanism and contains functions used by both nucleation
    and growth: hydride formation energy, volume fraction, lever rule.
 */

#include <stdio.h>
#include "Mechanism.hpp"
#include "PhysicsConstants.h"

using namespace std ;

class Precipitation : public Mechanism
{
    public:
        Precipitation(Sample* sample) ;
    
    // Define the coefficients used to compute the formation energy
        static void defineEnergyPolynomial(double Eth0, double Eth1, double Eth2, double Eth3) ;
    
    // Compute the coefficent used to compute nucleation and growth kinetids (volume fraction, lever rule, chi)
        static void computeCoeffs(Sample* sample) ;
        
    protected:
    
        static vector<double> _lever_rule ;
        static vector<double> _f_alpha ;
        
        static void compute_leverRule(Sample* sample) ;
        static void compute_f_alpha(Sample* sample) ;
        static void compute_chi(Sample* sample) ;

    // Compute the atom fraction of hydrogen at the alpha/alpha+delta frontier of the phase diagram
        static double compute_xdelta(double T) ;
    
    // Convert a concentration c from wt.ppm to atom fraction
        static double convertToAtomFrac(double c) ;
        
    // Compute the hydride formation energy (temperature polnomial)
        static double formation_energy(double T) ;
    
        static double _Eth0, _Eth1, _Eth2, _Eth3 ;
    
    private:
        static vector<double> * _temperature ;      // Temperature profile
        static vector<double> * _tssd ;             // Solubility profile
        static vector<double> * _totalContent ;     // Hydrogen profile
        static vector<double> * _hydrideContent ;   // Hydride profile
};

#endif /* Precipitation_hpp */
