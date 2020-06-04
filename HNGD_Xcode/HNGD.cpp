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
                physicalParameters[2])),    // Kn

    _growth(new Growth(_sample,
                physicalParameters[7],      // Kmob0
                physicalParameters[8],      // Kth0
                physicalParameters[9],      // Eg
                2.5)),                      // p

    _diffusion(new Diffusion(_sample,
                physicalParameters[14],     // D0
                physicalParameters[15],     // Ed
                physicalParameters[16]))    // Q
{
    _NbCells = (int)settings[2] ;
    
    Precipitation :: defineEnergyPolynomial(
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
    
//    // Compute the initial equilibrium
//    _sample->computeEquilibrium();
}

void HNGD :: compute()
{
    
    Precipitation :: computeCoeffs(_sample);
    _sample->computeTSS() ;
    
    // Compute time step
    if(_auto_dt)
        computeTimeStep() ;
    
    // Compute hydrogen flux
    vector<double> flux = _diffusion->computeFlux() ;
    
    // Compute new hydrogen distribution
    vector<double> position = _sample->returnPosition() ;
    vector<double> c_ss = _sample->returnSolutionContent() ;
    vector<double> c_tot = _sample->returnTotalContent() ;
    vector<double> c_prec = _sample->returnHydrideContent() ;
    
    vector<double> new_c_ss(_NbCells) ;
    new_c_ss[0] = c_ss[0] - _dt * (flux[0]) / (position[1] - position[0]);
    for(int k=1; k<_NbCells; k++)
        new_c_ss[k] = c_ss[k] - _dt * (flux[k] - flux[k-1]) / (position[k] - position[k-1]) ;
    
    _sample->setSolutionContent(new_c_ss) ;
    
    
    // Nucleation - Growth - Dissolution
    
    vector<double> tssp = _sample->returnTSSp() ;
    vector<double> tssd = _sample->returnTSSd() ;
    
    vector<double> rateNuc = _nucleation->computeRate() ;
    vector<double> rateGro = _growth->computeRate() ;
    vector<double> rateDis = _dissolution->computeRate() ;
    
    vector<double> rate(_NbCells, 0.);
    for(int k=0; k<_NbCells; k++)
    {
        if(c_ss[k] > tssp[k])
            rate[k] += rateNuc[k] ;
        
        if(c_ss[k] > tssd[k] && (c_prec[k] > 0 || rate[k] < 0))
            rate[k] += rateGro[k] ;
        
        if(c_ss[k] < tssd[k] && (c_prec[k] > 0 || rate[k] < 0))
            rate[k] += rateDis[k] ;
    }
    
    
    for(int k=0; k<_NbCells; k++)
    {
        if(c_prec[k] + _dt * rate[k] < 0 || c_ss[k] + _dt * rate[k] < 0)
            std::cout << "Negative Concentration" << std::endl ;
        
        c_ss[k] += _dt * rate[k] ;
        c_prec[k] -= _dt * rate[k] ;
    }
    _sample->setSolutionContent(c_ss) ;
    _sample->setHydrideContent(c_prec);
    _sample->updateTotalContent() ;
}

void HNGD :: getInput(vector<double> pos_temp, vector<double> temp_inp)
{
    _sample->spatialeInterpolation(pos_temp, temp_inp, _sample->returnTemperature());
}


void HNGD :: computeTimeStep()
{
    double dt = _diffusion->timeStep() ;
    dt = min(dt, _dissolution->timeStep()) ;
    dt = min(dt, _nucleation->timeStep()) ;
    dt = min(dt, _growth->timeStep()) ;
    
    _dt = dt ;
    
    if(_dt < 0.)
        cout << 'p' ;
}
