#ifndef Dissolution_hpp
#define Dissolution_hpp

#include <stdio.h>
#include <vector>
#include <math.h>
#include "Mechanism.hpp"
#include "PhysicsConstants.h"

using namespace std ;

class Dissolution : public Mechanism
{
    
    public:
        Dissolution(Sample* sample, double Kd0, double Ed);
        void computeKinetics()  ;
        void computeDrivForce() ;
    
    private:
        const double _Kd0 ;
        const double _Ed ;
    
};


#endif /* Dissolution_hpp */
