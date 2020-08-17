/**
    This class implements hydride growth for the HNGD model
*/

#ifndef Growth_hpp
#define Growth_hpp

#include <stdio.h>
#include "Precipitation.hpp"
#include "PhysicsConstants.h"

class Growth : public Precipitation
{
public:
    Growth(Sample* sample, double Kmob0, double Kth0, double Eg, double p);
    
    // Compute the kinetic coefficient at each position
    void computeKinetics()  ;
    
    // Compute the driving force at each position
    void computeDrivForce() ;
    
private:
    const double _Kth0 ;    // Preexponential factor for reaction-driven growth
    const double _Kmob0 ;   // Preexponential factor for diffusion-driven growth
    const double _Eg ;      // Activation energy for diffusion-driven growth
    const double _p ;       // Dimensionnality of growth (Avrami parameter)
    
    const vector<double> * _temperature ;   // Temperature profile
    const vector<double> * _tssd ;          // Solubility profile
    const vector<double> * _Cprec ;         // Hydride profile
    const vector<double> * _Ctot ;          // hydrogen profile
};


#endif /* Growth_hpp */
