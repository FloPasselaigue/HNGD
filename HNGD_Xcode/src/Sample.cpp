#include "Sample.hpp"
#include "PhysicsConstants.h"
#include <iostream>

Sample :: Sample(int nbCells, double bias, double tssp0, double Qp,
                 double tssd0, double Qd, double tau, double delta, double g):

    _nbCells(nbCells),
    _bias   (bias),

    _tssp0  (tssp0),
    _Qp     (Qp),

    _tssd0  (tssd0),
    _Qd     (Qd),

    _tau    (tau),

    _delta  (delta),
    _g      (g)
{
    _position       = vector<double>(_nbCells) ;
    _temperature    = vector<double>(_nbCells) ;
    _solutionContent= vector<double>(_nbCells) ;
    _hydrideContent = vector<double>(_nbCells) ;
    _totalContent   = vector<double>(_nbCells) ;
    _tssd           = vector<double>(_nbCells) ;
    _tssp           = vector<double>(_nbCells) ;

    _t_since_T_changed = 0. ;
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
    for(int k=0; k<_nbCells; k++)
    {

        // Default solubility fit
        double tssd = _tssd0 * exp(-_Qd / (R * _temperature[k])) ;

        // Hydride volume fraction
        double vf = vol_frac(k) ;

        // Modifying function
        double f_k = vf * (_g - ((1-_delta) * tssd + _g) * vf) ;

        // Modified solubility
        _tssd[k] = _tssd0 * exp(-_Qd / (R * _temperature[k])) + f_k ;


        // Default supersolubility
        double tssp = _tssp0 * exp(-_Qp / (R * _temperature[k])) ;

        // Modified supersolubility
        _tssp[k] = _tssd[k] + (tssp - _tssd[k]) * exp(-_t_since_T_changed / _tau) ;

    }
}

double Sample :: vol_frac(int k)
{
        // At. fraction of H at delta / delta+alpha boundary
    double T = _temperature[k] ;
    double x_delta = -9.93e-11*pow(T,3) + 8.48e-8*pow(T,2) - 5.73e-5*T + 0.623 ;

        // At. fraction solubility
    double x_alpha = _tssd[k]/Mh / (_tssd[k]/Mh + (1e6-_tssd[k])/Mzr) ;

        // Volume fraction
    double vf = _hydrideContent[k] / Mh / (_totalContent[k]/Mh + (1e6-_totalContent[k])/Mzr) / (x_delta-x_alpha) ;

    if (vf > 1.)
        vf = 1. ;

    return vf ;
}

// Domain definition
void Sample :: computeLocations(double x0, double xEnd, int _geometry)
{
    if (_geometry > 0)
    {
        // Polar
        // bias not implemented in polar geometry

        const double initialLenght = 2*M_PI/_nbCells;
        _position[0] = x0 ;
        for (int k=1; k<_nbCells; k++)
            _position[k] = _position[k-1] + initialLenght;
    }

    else
    {
        // Linear
        double  sum = 1. + _bias ;
        for(int k=0; k<_nbCells-3; k++)
            sum = 1. + _bias*sum ;

        const double initialLenght = (xEnd - x0)/sum ;

        _position[0] = x0 ;
        _position[1] = x0 + initialLenght ;

        for(int k=2; k<_nbCells-1; k++)
            _position[k] = _position[k-1] + _bias*(_position[k-1] - _position[k-2]) ;

        _position[_nbCells-1] = xEnd  ;
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
