#include "Dissolution.hpp"


Dissolution :: Dissolution(Sample* sample, double Kd0, double Ed) :
    Mechanism(sample),
    _Kd0(Kd0),
    _Ed(Ed)
{}

void Dissolution :: computeKinetics()
{
    vector<double> temperature = _sample->returnTemperature() ;
    
    for(int k=0; k<_nbCells; k++)
        _kinetic_factor[k] = (_Kd0 * exp(-_Ed / (kb * temperature[k]))) ;
}

void Dissolution :: computeDrivForce()
{
    vector<double> tssd = _sample->returnTSSp() ;
    vector<double> c_ss = _sample->returnSolutionContent() ;
    vector<double> c_tot= _sample->returnTotalContent() ;
    
    // Dissolution driving force cannot be negative
    for(int k=0; k<_nbCells; k++)
        _driving_force[k] = max(min(tssd[k]-c_ss[k], c_tot[k]-c_ss[k]), 0.) ;
    
}
