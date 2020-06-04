#include "Nucleation.hpp"

Nucleation :: Nucleation(Sample* sample, double Kn0):
    Precipitation(sample),
    _Kn0(Kn0)
{}

void Nucleation :: computeKinetics()
{
    vector<double> temperature = _sample->returnTemperature() ;
    
    for(int k=0; k<_nbCells; k++)
        _kinetic_factor[k] = _Kn0 * _f_alpha[k] * _chi[k]
        * exp(-formation_energy(temperature[k]) / (kb * temperature[k]) ) ;
}
    
void Nucleation :: computeDrivForce()
{
    vector<double> c_ss = _sample->returnSolutionContent() ;
    vector<double> tssp = _sample->returnTSSp() ;
    
    // Nucleation driving force cannot be positive
    for(int k=0; k<_nbCells; k++)
        _driving_force[k] = min( tssp[k] - c_ss[k] , 0.) ;
}


