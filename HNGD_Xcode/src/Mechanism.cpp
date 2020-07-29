#include "Mechanism.hpp"

// Static membres initialization
double Mechanism::_Eth0(0.) ;
double Mechanism::_Eth1(0.) ;
double Mechanism::_Eth2(0.) ;
double Mechanism::_Eth3(0.) ;

// Constructor
Mechanism :: Mechanism(Sample* sample):
    _nbCells(sample->returnNbCells()),
    _kinetic_factor(_nbCells),
    _driving_force(_nbCells),
    _rate(_nbCells),
    _gamma(0.187),
    _v(1.64e-5)
{}

// Formation energy fit
void Mechanism :: defineEnergyPolynomial(double Eth0, double Eth1, double Eth2, double Eth3)
{
    _Eth0 = Eth0 ;
    _Eth1 = Eth1 ;
    _Eth2 = Eth2 ;
    _Eth3 = Eth3 ;
}

// Compute the rate of the mechanism
vector<double>& Mechanism :: computeRate()
{
    computeKinetics();
    computeDrivForce();
    
    for(int k=0; k<_nbCells; k++)
        _rate[k] = _kinetic_factor[k] * _driving_force[k]  ;
    
    return _rate ;
}

// Compute the time step for this mechanism
double Mechanism :: timeStep()
{
    double max_K = *max_element(_kinetic_factor.begin(), _kinetic_factor.end()) ;
    
    if(max_K > 0)
        return .5 / max_K ;
    
    else
        return 1e6 ;
}

double Mechanism :: formation_energy(double T)
{
    return -_Eth0 + _Eth1 * T - _Eth2 * pow(T,2) + _Eth3 * pow(T,3);
}

double Mechanism :: volume_energy(double T)
{
    return formation_energy(T) * Na * e / _v ; //J/m3
}


// Getters
vector<double> & Mechanism :: returnKinetics()  {return _kinetic_factor;}
vector<double> & Mechanism :: returnDrivForce() {return _driving_force ;}

vector<double> & Mechanism :: returnRate() {return _rate ;}
