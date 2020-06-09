#ifndef Precipitation_hpp
#define Precipitation_hpp

#include <stdio.h>
#include "Mechanism.hpp"
#include "PhysicsConstants.h"

using namespace std ;

class Precipitation : public Mechanism
{
    public:
        Precipitation(Sample* sample) ;
        static void defineEnergyPolynomial(double Eth0, double Eth1, double Eth2, double Eth3);
        static void computeCoeffs(Sample* sample) ;
        
    protected:
        static double _Eth0, _Eth1, _Eth2, _Eth3 ;
    
        static vector<double> _lever_rule ;
        static vector<double> _f_alpha ;
        static vector<double> _chi ;
        
        static double compute_xdelta(double T) ;
        static double convertToAtomFrac(double c) ;
        static double formation_energy(double T) ;
        
        static void compute_leverRule(Sample* sample) ;
        static void compute_f_alpha(Sample* sample) ;
        static void compute_chi(Sample* sample) ;
};

#endif /* Precipitation_hpp */
