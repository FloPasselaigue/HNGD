#include "Growth.hpp"
#include <iostream>

Growth :: Growth(Sample* sample, double Kmob0, double Kth0, double Eg, double p) :
    Precipitation   (sample),

    _Kmob0          (Kmob0),
    _Kth0           (Kth0),
    _Eg             (Eg),
    _p              (p),

    _temperature    (& (sample->returnTemperature())),
    _tssd           (& (sample->returnTSSd())),
    _Cprec          (& (sample->returnHydrideContent())),
    _Ctot           (& (sample->returnTotalContent()))
{}

void Growth :: computeKinetics()
{
    for(int k=0; k<_nbCells; k++)
    {
        // Each component of growth kinetics follow an Arrhenius law, modified by
        // a couple of factors depending on the hydrogen and hydride contents
        double Kmob = _Kmob0 * _f_alpha[k] * _lever_rule[k] * exp(-_Eg / (kb * (*_temperature)[k]));
        double Kth  = _Kth0  * _f_alpha[k] * _lever_rule[k] * exp(-formation_energy((*_temperature)[k]) / (kb * (*_temperature)[k]));
        
        _kinetic_factor[k] = 1. / (1./Kmob + 1./Kth) ;
    }
}

void Growth :: computeDrivForce()
{
    for(int k=0; k<_nbCells; k++)
    {
        // Advancement of the (1-x)H + xZr -> Zr_xH_(1-x)
        double x = (*_Cprec)[k] / ((*_Ctot)[k] - (*_tssd)[k]) ; 
        
        if(x>0. && x<1)
            _driving_force[k] = ((*_tssd)[k] - (*_Ctot)[k]) * _p * (1-x) * pow(-log(1-x), 1-1./_p) ;
        
        else
            _driving_force[k] = 0. ;
    }
}
