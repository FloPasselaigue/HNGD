#ifndef Diffusion_hpp
#define Diffusion_hpp

#include <stdio.h>
#include <vector>
#include "Sample.hpp"
#include "PhysicsConstants.h"

class Diffusion
{
public:
    Diffusion(Sample * sample, double D0, double Ed, double Q) ;
    
    vector<double> computeFlux() ;
    
    double timeStep() ;
    
    vector<double>& returnFlux() ;
    
private:
    Sample * _sample ;
    const int _nbCells ;
    const double _D0 ;
    const double _Ed ;
    const double _Q ;
    
    vector<double> _coeff_Fick ;
    vector<double> _flux ;
};

#endif /* Diffusion_hpp */
