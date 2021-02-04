#include "Nucleation.hpp"

Nucleation :: Nucleation(Sample* sample, double Kn0):
    Precipitation(sample),

    _Kn0         (Kn0),

    _temperature (& (sample->returnTemperature())),
    _tssp        (& (sample->returnTSSp())),
    _Css         (& (sample->returnSolutionContent()))
{}

void Nucleation :: computeKinetics()
{
    // Nucleation kinetics follow an Arrhenius law modified my a couple of factors
    // that include a hydrogen/hydride content dependency
    for(int k=0; k<_nbCells; k++)
        _kinetic_factor[k] = _Kn0 * _f_alpha[k]
        * exp( -formation_energy((*_temperature)[k]) / (kb * (*_temperature)[k]) ) ;

}
    
void Nucleation :: computeDrivForce()
{
    // Nucleation driving force cannot be positive
    for(int k=0; k<_nbCells; k++)
        _driving_force[k] = min( (*_tssp)[k] - (*_Css)[k] , 0.) ;
}

