/**
This class implements hydride dissolution for the HNGD model
*/

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
    
        // Compute the kinetic coefficient at each position
        void computeKinetics()  ;
    
        // Compute the driving force at each position
        void computeDrivForce() ;
    
    private:
        const double _Kd0 ; // Preexponential factor for dissolution kinetics
        const double _Ed ;  // Activation energy for diffusion
        
        const vector<double> * _temperature ; // temperature profile
        const vector<double> * _tssd ;        // solubility profile
        const vector<double> * _Css ;         // solid solution profile
        const vector<double> * _Ctot ;        // hydrogen profile
};


#endif /* Dissolution_hpp */
