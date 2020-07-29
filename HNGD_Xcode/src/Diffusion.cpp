#include <stdio.h>
#include "Diffusion.hpp"

Diffusion :: Diffusion(Sample * sample, double D0, double Ed, double Q):
    _nbCells(sample->returnNbCells()),
    _D0(D0),
    _Ed(Ed),
    _Q(Q),
    _coeff_Fick(vector<double>(_nbCells)),
    _flux(vector<double>(_nbCells)),
    _positions(& (sample->returnPosition())),
    _temperature(& (sample->returnTemperature())),
    _Css(& (sample->returnSolutionContent()))
{}

void Diffusion :: computeCoeff()
{
    for(int k=0; k<_nbCells; k++)
        _coeff_Fick[k] = _D0 * exp(-_Ed / (kb * (*_temperature)[k])) ;
}

vector<double> Diffusion :: computeFlux()
{
    double dC_dx(0.), dT_dx(0.) ;
    double flux_fick(0.), flux_soret(0.) ;
    
    for(int k=0; k<_nbCells-1; k++)
    {
        dC_dx = ((*_Css)[k+1] - (*_Css)[k]) / ((*_positions)[k+1] - (*_positions)[k]) ;
        flux_fick = - _coeff_Fick[k] * dC_dx ;
        
        dT_dx = ((*_temperature)[k+1] - (*_temperature)[k]) / ((*_positions)[k+1] - (*_positions)[k]) ;
        flux_soret = - _coeff_Fick[k] * _Q * (*_Css)[k] * dT_dx / (R * pow((*_temperature)[k], 2)) ;
        
        _flux[k] = flux_fick + flux_soret ;
    }
    
    return _flux ;
}

double  Diffusion :: timeStep()
{
    double dt = pow((*_positions)[1] - (*_positions)[0], 2) / (2 * _coeff_Fick[0]) ;
    for(int k=2; k<_nbCells; k++)
        dt = min(dt, pow((*_positions)[k] - (*_positions)[k-1], 2) / (2 * _coeff_Fick[k]));
    
    return dt ;
}

vector<double>&  Diffusion :: returnFlux() {return _flux;}
