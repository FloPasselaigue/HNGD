#include <stdio.h>
#include "Diffusion.hpp"

// TODO: geometry and radius can be taken from the sample
Diffusion :: Diffusion(Sample * sample, double D0, double Ed, double Q, int geometry, double radius):
    _nbCells    (sample->returnNbCells()),

    _D0         (D0),
    _Ed         (Ed),
    _Q          (Q),
    _coeff_Fick (vector<double>(_nbCells)),
    _flux       (vector<double>(_nbCells)),
    _dC_dx      (vector<double>(_nbCells)),
    _dT_dx      (vector<double>(_nbCells)),
    _geometry   (geometry),
    _radius     (radius),

    _positions  (& (sample->returnPosition())),
    _temperature(& (sample->returnTemperature())),
    _Css        (& (sample->returnSolutionContent())),
    _Cprec      (& (sample->returnSolutionContent()))
{}

void Diffusion :: computeCoeff()
{
    // The diffusion coefficient follows an Arrhenius law
    for(int k=0; k<_nbCells; k++)
        _coeff_Fick[k] = _D0 * exp(-_Ed / (kb * (*_temperature)[k])) ;
}

void Diffusion :: computeGradient()
{    // Find Temp and Solid Solution gradients at each node, depending on geometry type
    if (_geometry > 0)
        { // Polar Geometry
            for (int k=0; k<_nbCells-1; k++)
                {
                    _dC_dx[k] = 1/_radius * ((*_Css)[k+1] - (*_Css)[k]) / ((*_positions)[k+1] - (*_positions)[k]);
                    _dT_dx[k] = 1/_radius * ((*_temperature)[k+1] - (*_temperature)[k]) / ((*_positions)[k+1] - (*_positions)[k]);
                }

            // Last node's gradient is function of that node's and first node's Css
                _dC_dx[_nbCells-1] = 1/_radius * ((*_Css)[0] - (*_Css)[_nbCells-1]) / (2*M_PI - (*_positions)[_nbCells-1]);
                _dT_dx[_nbCells-1] = 1/_radius * ((*_temperature)[0] - (*_temperature)[_nbCells-1]) / (2*M_PI - (*_positions)[_nbCells-1]);

        }
    else
        { // Linear Geometry
        for (int k=0; k<_nbCells-1; k++)
            {
                _dC_dx[k] = ((*_Css)[k+1] - (*_Css)[k]) / ((*_positions)[k+1] - (*_positions)[k]) ;
                _dT_dx[k] = ((*_temperature)[k+1] - (*_temperature)[k]) / ((*_positions)[k+1] - (*_positions)[k]) ;
            }
        }
}

void Diffusion :: computeFlux()
{
    double flux_fick(0.), flux_soret(0.), _avgtemp(0.), _avgCss(0.) ;
    
    if (_geometry > 0)
        { // Polar Geometry
            for(int k=0; k<_nbCells; k++)
            {
                if (k==_nbCells-1) // Last node of sample
                    {
                    // Last node uses the first and last node's values
                    _avgtemp = ((*_temperature)[0]+(*_temperature)[_nbCells-1])/2 ;
                    _avgCss = ((*_Css)[0]+(*_Css)[_nbCells-1])/2 ;
                    }

                else
                    {
                    _avgtemp = ((*_temperature)[k]+(*_temperature)[k+1])/2 ;
                    _avgCss = ((*_Css)[k]+(*_Css)[k+1])/2 ;
                    }
                flux_fick = - _coeff_Fick[k] * _dC_dx[k] ;
                flux_soret = - _coeff_Fick[k] * _Q * _avgCss * _dT_dx[k] / (R * pow(_avgtemp, 2)) ;

                _flux[k] = flux_fick + flux_soret;
            }
        }
    else
        { // Linear Geometry
            for(int k=0; k<_nbCells-1; k++)
                {
                    flux_fick = - _coeff_Fick[k] * _dC_dx[k] ;
                    flux_soret = - _coeff_Fick[k] * _Q * (*_Css)[k] * _dT_dx[k] / (R * pow((*_temperature)[k], 2)) ;

                    _flux[k] = flux_fick + flux_soret ;
                }
        }
}

double  Diffusion :: timeStep()
{
    double dt ;
 
    // The time step associated with diffusion is computed using the convergence criterion from finite element theory
    if (_geometry>0)
    {
        // Polar, needs radius for cell size
        dt = pow(_radius*((*_positions)[1] - (*_positions)[0]), 2) / (2 * _coeff_Fick[0]) ;
        for(int k=2; k<_nbCells; k++)
            dt = min(dt, pow(_radius*((*_positions)[k] - (*_positions)[k-1]), 2) / (2 * _coeff_Fick[k]));
    }
    else
    {
        // Linear
        dt = pow((*_positions)[1] - (*_positions)[0], 2) / (2 * _coeff_Fick[0]) ;
        for(int k=2; k<_nbCells; k++)
            dt = min(dt, pow((*_positions)[k] - (*_positions)[k-1], 2) / (2 * _coeff_Fick[k]));
    }
    
    return dt ;
}

vector<double>&  Diffusion :: returnFlux() {return _flux;}
