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
        const vector<double> * _temperature ;
        const vector<double> * _tssd ;
        const vector<double> * _Css ;
        const vector<double> * _Ctot ;
};


#endif /* Dissolution_hpp */
