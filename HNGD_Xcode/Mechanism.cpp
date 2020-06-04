#include "Mechanism.hpp"


// Constructor
Mechanism :: Mechanism(Sample* sample):
    _nbCells(sample->returnNbCells()),
    _kinetic_factor(_nbCells),
    _driving_force(_nbCells),
    _rate(_nbCells)
    {_sample =  sample ;}


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
    double dt = 10. ;
    
    for(int k=0; k<_nbCells; k++)
        if(_driving_force[k] != 0.)
            dt = min(dt, 0.3/_kinetic_factor[k]);
    
    return dt ;
}


// Getters
vector<double> Mechanism :: returnKinetics()  {return _kinetic_factor;}
vector<double> Mechanism :: returnDrivForce() {return _driving_force ;}

vector<double> Mechanism :: returnRate() {return _rate ;}
