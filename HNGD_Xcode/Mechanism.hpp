#ifndef Mechanism_hpp
#define Mechanism_hpp

#include <stdio.h>
#include <vector>
#include "Sample.hpp"

using namespace std ;

class Mechanism
{
    public:
        Mechanism(Sample* sample) ;
    
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
    
};

#endif /* Mechanism_hpp */
