#include "Mechanism.hpp"

// Constructor
Mechanism :: Mechanism(Sample* sample):
    _nbCells(sample->returnNbCells()),
    _kinetic_factor(_nbCells),
    _driving_force(_nbCells),
    _rate(_nbCells)
{}

// Compute the rate of the mechanism as the product of kinetics and driving force
void Mechanism :: computeRate()
{
    computeKinetics();
    computeDrivForce();
    
    for(int k=0; k<_nbCells; k++)
        _rate[k] = _kinetic_factor[k] * _driving_force[k]  ;
    
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

// Getters
vector<double> & Mechanism :: returnKinetics()  {return _kinetic_factor;}
vector<double> & Mechanism :: returnDrivForce() {return _driving_force ;}

vector<double> & Mechanism :: returnRate() {return _rate ;}
