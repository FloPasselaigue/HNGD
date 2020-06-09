#include "Nucleation.hpp"

Nucleation :: Nucleation(Sample* sample, double Kn0):
    Precipitation(sample),
    _Kn0(Kn0)
//    _gamma(0.187),
//    _v(1.64e-5)
{}

void Nucleation :: computeKinetics()
{
    vector<double> temperature = _sample->returnTemperature() ;
    
    for(int k=0; k<_nbCells; k++)
        _kinetic_factor[k] = _Kn0 * _f_alpha[k] * _chi[k]
        * exp(-formation_energy(temperature[k]) / (kb * temperature[k]) ) ;
    
//    vector<double> c_ss = _sample->returnSolutionContent() ;
//    double a = 3.2e-10 ;
//
//    for(int k=0; k<_nbCells; k++)
//    {
//        double D = 1.08e-6 * exp(-0.46 / (kb*temperature[k])) ;
//        double x = convertToAtomFrac(c_ss[k]) ;
//        double Gv = volume_energy(temperature[k]) ;
//
//        _kinetic_factor[k] = 16 * 3.14 * pow(_gamma,2) * D * x / (pow(Gv, 2) * pow(a, 4)) ;
//    }
    
}
    
void Nucleation :: computeDrivForce()
{
    vector<double> c_ss = _sample->returnSolutionContent() ;
    vector<double> tssp = _sample->returnTSSp() ;
    
//    // Nucleation driving force cannot be positive
//    for(int k=0; k<_nbCells; k++)
//        _driving_force[k] = min( tssp[k] - c_ss[k] , 0.) ;
    
    
    vector<double> tssd = _sample->returnTSSd() ;
    double delta = 0 ;
    for(int k=0; k<_nbCells; k++)
    {
        double sigma = 0.15 * tssp[k] ; 
        double proba = 0.5 * erfc((delta + tssp[k] - c_ss[k]) / (sqrt(2) * sigma)) ;
        _driving_force[k] = min(0., proba * (tssd[k] - c_ss[k]) ) ;
    }
    
//    vector<double> temperature = _sample->returnTemperature() ;
//    for(int k=0; k<_nbCells; k++)
//    {
//        double Gv = volume_energy(temperature[k]) ;
//        double Gstar = 16 * 3.14 * pow(_gamma,3) / (3 * pow(Gv, 2)) / e ;
//        _driving_force[k] = 0.5 * c_ss[k] * exp(-Gstar / (kb * temperature[k])) ;
//    }
    
}

//double Nucleation :: volume_energy(double T)
//{
//    return formation_energy(T) * Na * e / _v ; //J/m3
//}

