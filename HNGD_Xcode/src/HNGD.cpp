#include "HNGD.hpp"
#include <iostream>

HNGD :: HNGD(double* settings, double* physicalParameters): 

    _sample(new Sample((int)settings[0],    // number of cells
                settings[1],                // bias
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
                physicalParameters[1],      // Ed
                physicalParameters[15],     // Q
                settings[6],                // Geometry Type
                settings[2])),              // Radius or Sample Length

    _Css     (& (_sample->returnSolutionContent())),// Solid solution profile
    _Ctot    (& (_sample->returnTotalContent())),   // Hydrogen profile
    _Cprec   (& (_sample->returnHydrideContent())), // Hydride profile
    _position(& (_sample->returnPosition())),       // Positions vector
    _tssp    (& (_sample->returnTSSp())),           // Supersolubility profile
    _tssd    (& (_sample->returnTSSd())),           // Solubility profile

    _flux    (& (_diffusion->returnFlux())),        // Hydrogen flux profile

    _rateNuc (& (_nucleation->returnRate())),       // Nucleation rate at each position
    _rateGro (& (_growth->returnRate())),           // Growth rate at each position
    _rateDis (& (_dissolution->returnRate()))       // Dissolution rate at each position
{
    _NbCells  = (int)settings[0] ;   // Number of nodes
    _geometry = (int)settings[6] ;   // Geometry Type
    _radius   = settings[2]; // Radius or Sample Length
    
    Precipitation :: defineEnergyPolynomial(
                    physicalParameters[3],  // Eth0
                    physicalParameters[4],  // Eth1
                    physicalParameters[5],  // Eth2
                    physicalParameters[6]); // Eth3
    
    // Create geometry
    double xEnd ;
    if(_geometry==1)
        xEnd = 2*M_PI ;
    else
        xEnd = settings[2] ;

    _sample->computeLocations(0., xEnd, _geometry) ;
    
    // Time step management
    if(settings[3] < 0.)
        _auto_dt = true ;
    else
    {
        _auto_dt = false ;
        _dt = settings[3] ;
    }
    
}


void HNGD :: getInitialConditions(vector<double> pos_hyd, vector<double> hyd_inp,
                                  vector<double> pos_temp,vector<double> temp_inp)
{
    
    // Create initial hydrogen profile
    _sample->spatialeInterpolation(pos_hyd, hyd_inp, _sample->returnTotalContent()) ;
    _sample->spatialeInterpolation(pos_hyd, hyd_inp, _sample->returnSolutionContent()) ;
    _sample->setHydrideContent(vector<double>(_NbCells,0.));
    
    // Create initial temperature profile
    _sample->spatialeInterpolation(pos_temp, temp_inp, _sample->returnTemperature()) ;
    
    // Compute solubilities
    _sample->computeTSS() ;
    
    // Compute the initial equilibrium
    _sample->computeEquilibrium();
}



void HNGD :: compute()
{
    //                  ---- COMPUTE THE PHYSICAL PARAMETERS ----
    
    // Compute the coefficients used to compute nucleation and growth kinetics
    Precipitation :: computeCoeffs(_sample);
    
    // Compute the supersolubility and solubility profiles
    _sample->computeTSS() ;
    
    // Compute the diffusion coefficient at each position
    _diffusion->computeCoeff() ;
    
    // Compute time step
    if(_auto_dt)
        computeTimeStep() ;
    
    
    
    //                ---- COMPUTE THE NEW SOLID SOLUTION PROFILE ----
    
    // Compute hydrogen flux
    _diffusion->computeGradient() ;
    _diffusion->computeFlux() ;
    
    // Compute new hydrogen distribution
    vector<double> new_c_ss(_NbCells) ;
    
    if (_geometry>0)
    {
        // Polar Geometry
        new_c_ss[0] = (*_Css)[0] - _dt * ((*_flux)[0] - (*_flux)[_NbCells-1]) / (_radius*(2*M_PI - (*_position)[_NbCells-1]));
        for (int k=1; k<_NbCells; k++)
            new_c_ss[k] = (*_Css)[k] - _dt * ((*_flux)[k] - (*_flux)[k-1]) / (_radius*((*_position)[k] - (*_position)[k-1])) ;
    }

    else
    {
        // Linear Geometry
        new_c_ss[0] = (*_Css)[0] - _dt * ((*_flux)[0]) / ((*_position)[1] - (*_position)[0]);
        for(int k=1; k<_NbCells; k++)
            new_c_ss[k] = (*_Css)[k] - _dt * ((*_flux)[k] - (*_flux)[k-1]) / ((*_position)[k] - (*_position)[k-1]) ;
    }
    
    _sample->setSolutionContent(new_c_ss) ;
    _sample->updateTotalContent() ;
    
    
    
    
    //                    ---- NUCLEATION  GROWTH  DISSOLUTION ----
    new_c_ss = *_Css ;
    vector<double> new_c_prec = *_Cprec ;
    
    // Compute the rate of each mechanism
    _nucleation->computeRate() ;
    _growth->computeRate() ;
    _dissolution->computeRate() ;
    
    // Precipitation/dissolution algorithm to compute the total rate at each position
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
    
    // Compute the new solid solution and hydride profiles
    for(int k=0; k<_NbCells; k++)
    {
        new_c_ss[k] += _dt * rate[k] ;
        new_c_prec[k] -= _dt * rate[k] ;
    }

    // Safeguard against negative concentration
    if(*min_element(new_c_ss.begin(), new_c_ss.end())<0 || *min_element(new_c_prec.begin(), new_c_prec.end())<0)
        std::cout << "/!\\ Negative Concentration /!\\ " << std::endl ;
    
    // Real 0 content if it goes below 1e-10 wt.ppm
    for(int k=0; k<_NbCells; k++)
        if(new_c_prec[k] < 1e-10)
            new_c_prec[k] = 0 ;
    
    // Update the hydrogen profiles
    _sample->setSolutionContent(new_c_ss) ;
    _sample->setHydrideContent(new_c_prec);
    _sample->updateTotalContent() ;
}



void HNGD :: getInput(vector<double> pos_temp, vector<double> temp_inp)
{
    // Update the state of the sample to compute the next state
    _sample->spatialeInterpolation(pos_temp, temp_inp, _sample->returnTemperature());
}




void HNGD :: computeTimeStep()
{
    // Compute the time step associated with each phenomenon and use the smallest one
    double dt = 10 ;
    
    dt = min(dt, _diffusion->timeStep()) ;
    dt = min(dt, _dissolution->timeStep()) ;
    dt = min(dt, _nucleation->timeStep()) ;
    dt = min(dt, _growth->timeStep()) ;
    
    _dt = dt ;
}
