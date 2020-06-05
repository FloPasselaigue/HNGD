#include <stdio.h>
#include "Diffusion.hpp"

Diffusion :: Diffusion(Sample * sample, double D0, double Ed, double Q):
    _nbCells(sample->returnNbCells()),
    _sample(sample),
    _D0(D0),
    _Ed(Ed),
    _Q(Q),
    _coeff_Fick(vector<double>(_nbCells)),
    _flux(vector<double>(_nbCells))
{}

vector<double> Diffusion :: computeFlux()
{
    vector<double> positions   = _sample->returnPosition() ;
    vector<double> temperature = _sample->returnTemperature() ;
    vector<double> c_ss = _sample->returnSolutionContent() ;
    
    double dC_dx(0.), dT_dx(0.) ;
    double flux_fick(0.), flux_soret(0.) ;
    
    for(int k=0; k<_nbCells-1; k++)
    {
        dC_dx = (c_ss[k+1] - c_ss[k]) / (positions[k+1] - positions[k]) ;
        flux_fick = - _coeff_Fick[k] * dC_dx ;
        
        dT_dx = (temperature[k+1] - temperature[k]) / (positions[k+1] - positions[k]) ;
        flux_soret = - _coeff_Fick[k] * _Q * c_ss[k] * dT_dx / (R * temperature[k]) ;
        
        _flux[k] = flux_fick + flux_soret ;
    }
    
    return _flux ;
}

double  Diffusion :: timeStep()
{
    vector<double> positions   = _sample->returnPosition() ;
    
    double dt = pow(positions[1] - positions[0], 2) / (2 * _coeff_Fick[0]) ;
    for(int k=2; k<_nbCells; k++)
        dt = min(dt, pow(positions[k] - positions[k-1], 2) / (2 * _coeff_Fick[k]));
    
    return dt ;
}

vector<double>&  Diffusion :: returnFlux() {return _flux;}
