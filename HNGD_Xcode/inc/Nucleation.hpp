/**
    This class implements hydride nucleation for the HNGD model.
 */

#ifndef Nucleation_hpp
#define Nucleation_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include <math.h>
#include "Precipitation.hpp"
#include "PhysicsConstants.h"

using namespace std ;

class Nucleation : public Precipitation
{
public:
    Nucleation(Sample* sample, double Kn0);
    
    // Redefinition of virtual functions from Mechanism superclass
    void computeKinetics()  ;
    void computeDrivForce() ;

private:
    const double _Kn0 ; // Preexponential factor for nucleation kinetics
    
    const vector<double> * _temperature ; // Temperature profile
    const vector<double> * _tssp ;        // Supersolubility profile
    const vector<double> * _Css ;         // Solid solution profile
};

#endif /* Nucleation_hpp */
