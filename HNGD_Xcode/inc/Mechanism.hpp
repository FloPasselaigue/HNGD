#ifndef Mechanism_hpp
#define Mechanism_hpp

#include <stdio.h>
#include <vector>
#include "Sample.hpp"
#include "PhysicsConstants.h"

using namespace std ;

class Mechanism
{
    public:
        Mechanism(Sample* sample) ;
        static void defineEnergyPolynomial(double Eth0, double Eth1, double Eth2, double Eth3) ;
    
        vector<double>& computeRate() ;
        
        double timeStep() ;
    
        vector<double> returnKinetics();
        vector<double> returnDrivForce();
        vector<double> returnRate();
        
    protected:
        Sample * _sample ;
        const int _nbCells ;
    
        vector<double> _kinetic_factor ;
        vector<double> _driving_force ;
        vector<double> _rate ;

        virtual void computeKinetics() = 0 ;
        virtual void computeDrivForce()= 0 ;
    
        const double _gamma ;
        const double _v ;
        double volume_energy(double T) ;

        static double _Eth0, _Eth1, _Eth2, _Eth3 ;
        static double formation_energy(double T) ;
};

#endif /* Mechanism_hpp */
