#include "Growth.hpp"

Growth :: Growth(Sample* sample, double Kmob0, double Kth0, double Eg, double p) :
    Precipitation(sample),
    _Kmob0(Kmob0),
    _Kth0(Kth0),
    _Eg(Eg),
    _p(p)
{}

void Growth :: computeKinetics()
{
    vector<double> temperature = _sample->returnTemperature();
    
    for(int k=0; k<_nbCells; k++)
    {
        double Kmob = _Kmob0 * _f_alpha[k] * _lever_rule[k] * exp(-_Eg / (kb * temperature[k]));
        double Kth  = _Kth0  * _f_alpha[k] * _lever_rule[k] * exp(-formation_energy(temperature[k]) / (kb * temperature[k]));
        
        _kinetic_factor[k] = 1. / (1./Kmob + 1./Kth) ;
    }
}

void Growth :: computeDrivForce()
{
    vector<double> c_prec = _sample->returnHydrideContent() ;
    vector<double> c_tot = _sample->returnTotalContent() ;
    vector<double> tssd = _sample->returnTSSd() ;
    
    for(int k=0; k<_nbCells; k++)
    {
        double x = c_prec[k] / (c_tot[k] - tssd[k]) ; // (c_prec-delta)/(c_tot-tssd)
        
        if(x>0. && x<1)
            _driving_force[k] = (tssd[k] - c_tot[k]) * _p * (1-x) * pow(-log(1-x), 1-1./_p) ;
        
        else
            _driving_force[k] = 0. ;
    }
}
