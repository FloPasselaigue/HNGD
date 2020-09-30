#include "Sample.hpp"
#include "PhysicsConstants.h"
#include <iostream>

Sample :: Sample(int nbCells, double bias, double tssp0, double Qp, double tssd0, double Qd) :

    _nbCells(nbCells),
    _bias   (bias),

    _tssp0  (tssp0),
    _Qp     (Qp),

    _tssd0  (tssd0),
    _Qd     (Qd)
{
    _position       = vector<double>(_nbCells) ;
    _temperature    = vector<double>(_nbCells) ;
    _solutionContent= vector<double>(_nbCells) ;
    _hydrideContent = vector<double>(_nbCells) ;
    _totalContent   = vector<double>(_nbCells) ;
    _tssd           = vector<double>(_nbCells) ;
    _tssp           = vector<double>(_nbCells) ;
    
    _t_since_T_changed = 0. ;
    _tau = 1e4 ;
}

// Compute the equilibrium for the initial conditions
void Sample :: computeEquilibrium()
{
    for(int k=0; k<_nbCells; k++)
    {
        _solutionContent[k] = min(_totalContent[k], _tssd[k]) ;
        _hydrideContent[k] = _totalContent[k] - _solutionContent[k] ;
    }
}


// Solubilities computation
void Sample :: computeTSS()
{
    double delta = .6 ;
    double g = 6  ;

    double c = 1. ;
    double b = g  ;
    double a = delta - g - 1 ;
    
    for(int k=0; k<_nbCells; k++)
    {
        _tssd[k] = _tssd0 * exp(-_Qd / (R * _temperature[k])) ;
        _tssd[k] *= a * pow(_hydrideContent[k]/17000, 2) + b * _hydrideContent[k]/17000 + c ;
        
        double tssp = _tssp0 * exp(-_Qp / (R * _temperature[k])) ;
        _tssp[k] = _tssd[k] + (tssp - _tssd[k]) * exp(-_t_since_T_changed / _tau) ;
    }
}

// Domain definition
void Sample :: computeLocations(double x0, double xEnd)
{
    // If the bias is negative, the mesh refinement is made on the "right" side
    double bias ;
    if(_bias < 0.)
        bias = - _bias ;
    else
        bias = _bias ;
        
    double  sum = 1. + bias ;
    for(int k=0; k<_nbCells-3; k++)
        sum = 1. + bias*sum ;

    const double initialLenght = (xEnd - x0)/sum ;

    _position[0] = x0 ;
    _position[1] = x0 + initialLenght ;
    
    for(int k=2; k<_nbCells-1; k++)
        _position[k] = _position[k-1] + bias*(_position[k-1] - _position[k-2]) ;
    
    _position[_nbCells-1] = xEnd  ;
    
    if(_bias < 0.)
    {
        reverse(_position.begin(), _position.end()) ;
        for(int k=0; k<_nbCells; k++)
            _position[k] = xEnd - _position[k] ;
    }
        
}

// Interpolation
void Sample :: spatialeInterpolation(vector<double>& refX, vector<double>& refY, vector<double>& vectorY)
{
    int k = 1 ;
    for(int i=0; i<_position.size(); i++) {
        if(_position[i] > refX[k])
            k++ ;
        vectorY[i] = refY[k-1] + (refY[k] - refY[k-1]) * (_position[i] - refX[k-1]) / (refX[k] - refX[k-1]);
    }
}
