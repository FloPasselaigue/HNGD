#include "HNGD.hpp"
#include <iostream>

HNGD :: HNGD(double* settings, double* physicalParameters):

    _sample(new Sample((int)settings[2],    // number of cells
                settings[7],                // bias
                physicalParameters[10],     // tssp0
                physicalParameters[11],     // Qp
                physicalParameters[12],     // tssd0
                physicalParameters[13])),   // Qd

    _dissolution(new Dissolution(_sample,
                 physicalParameters[0],     // Kd0
                 physicalParameters[1])),   // Ed

    _nucleation(new Nucleation(_sample,
                physicalParameters[2])),    // Kn0

    _growth(new Growth(_sample,
                physicalParameters[7],      // Kmob0
                physicalParameters[8],      // Kth0
                physicalParameters[9],      // Eg
                2.5)),                      // p

    _diffusion(new Diffusion(_sample,
                physicalParameters[14],     // D0
                physicalParameters[15],     // Ed
                physicalParameters[16])),   // Q

    _Css     (& (_sample->returnSolutionContent())),
    _Ctot    (& (_sample->returnTotalContent())),
    _Cprec   (& (_sample->returnHydrideContent())),
    _position(& (_sample->returnPosition())),
    _tssp    (& (_sample->returnTSSp())),
    _tssd    (& (_sample->returnTSSd())),

    _flux    (& (_diffusion->returnFlux())),
    _rateNuc (& (_nucleation->returnRate())),
    _rateGro (& (_growth->returnRate())),
    _rateDis (& (_dissolution->returnRate()))
{
    _NbCells = (int)settings[2] ;
    
    Mechanism :: defineEnergyPolynomial(
                    physicalParameters[3],  // Eth0
                    physicalParameters[4],  // Eth1
                    physicalParameters[5],  // Eth2
                    physicalParameters[6]); // Eth3
    
    // Create geometry
    _sample->computeLocations(0., settings[8]) ;
    
    // Time step management
    if(settings[9] < 0.)
        _auto_dt = true ;
    else {
        _auto_dt = false ;
        _dt = settings[9] ;
    }
    
}


void HNGD :: getInitialConditions(vector<double> pos_hyd, vector<double> hyd_inp,
                                  vector<double> pos_temp,vector<double> temp_inp)
{
    
    // Create initial hydrogen profile
    _sample->spatialeInterpolation(pos_hyd, hyd_inp, _sample->returnTotalContent()) ;
    _sample->spatialeInterpolation(pos_hyd, hyd_inp, _sample->returnSolutionContent()) ;
    _sample->setHydrideContent(vector<double>(_NbCells,1e-3));
    
    // Create initial temperature profile
    _sample->spatialeInterpolation(pos_temp, temp_inp, _sample->returnTemperature()) ;
    
    // Compute solubilities
    _sample->computeTSS() ;
    
    // Compute the initial equilibrium
    _sample->computeEquilibrium();
}

void HNGD :: compute()
{
    
    Precipitation :: computeCoeffs(_sample);
    _sample->computeTSS() ;
    _diffusion->computeCoeff() ;
    
    // Compute time step
    if(_auto_dt)
        computeTimeStep() ;
    
    // Compute hydrogen flux
    vector<double> flux = _diffusion->computeFlux() ;
    
    // Compute new hydrogen distribution
    vector<double> new_c_ss(_NbCells) ;
    new_c_ss[0] = (*_Css)[0] - _dt * (flux[0]) / ((*_position)[1] - (*_position)[0]);
    for(int k=1; k<_NbCells; k++)
        new_c_ss[k] = (*_Css)[k] - _dt * (flux[k] - flux[k-1]) / ((*_position)[k] - (*_position)[k-1]) ;
    
    _sample->setSolutionContent(new_c_ss) ;
    _sample->updateTotalContent() ;
    
    
    // Nucleation - Growth - Dissolution
    new_c_ss = *_Css ;
    vector<double> new_c_prec = *_Cprec ;
    
    _nucleation->computeRate() ;
    _growth->computeRate() ;
    _dissolution->computeRate() ;
    
    vector<double> rate(_NbCells, 0.);
    for(int k=0; k<_NbCells; k++)
    {
        if((*_Css)[k] > (*_tssp)[k])
            rate[k] += (*_rateNuc)[k] ;
        
        if((*_Css)[k] > (*_tssd)[k] && ((*_Cprec)[k] > 0 || rate[k] < 0))
            rate[k] += (*_rateGro)[k] ;
        
        if((*_Css)[k] < (*_tssd)[k] && ((*_Cprec)[k] > 0 || rate[k] < 0))
            rate[k] += (*_rateDis)[k] ;
    }
    
    for(int k=0; k<_NbCells; k++)
    {
        new_c_ss[k] += _dt * rate[k] ;
        new_c_prec[k] -= _dt * rate[k] ;
    }

    if(*min_element(new_c_ss.begin(), new_c_ss.end())<0 || *min_element(new_c_prec.begin(), new_c_prec.end())<0)
        std::cout << "/!\\ Negative Concentration /!\\ " << std::endl ;
    
    _sample->setSolutionContent(new_c_ss) ;
    _sample->setHydrideContent(new_c_prec);
    _sample->updateTotalContent() ;
}

void HNGD :: getInput(vector<double> pos_temp, vector<double> temp_inp)
{
    _sample->spatialeInterpolation(pos_temp, temp_inp, _sample->returnTemperature());
}


void HNGD :: computeTimeStep()
{
    double dt = 10 ;
    
    dt = min(dt, _diffusion->timeStep()) ;
    dt = min(dt, _dissolution->timeStep()) ;
    dt = min(dt, _nucleation->timeStep()) ;
    dt = min(dt, _growth->timeStep()) ;
    
    _dt = dt ;
}
