/**
    This class implements the HNGD model. It uses secondary classes
    for each phenomenon: hydrogen Diffusion, hydride Nucleation, Growth and Dissolution.
    An additional class Sample is used to create the geometry, manage the temperature
    profile and compute the solubility/supersolubility profiles.
 */

#ifndef HNGD_hpp
#define HNGD_hpp

#include <stdio.h>
#include "Sample.hpp"
#include "Diffusion.hpp"
#include "Nucleation.hpp"
#include "Growth.hpp"
#include "Dissolution.hpp"
#include <algorithm>

class HNGD
{
    public:
        HNGD(double* settings, double* physicalParameters, double xEnd, int geometry); //TODO: geometry is part of settings
        
    // Use the information taken from the input files to create the initial state
        void getInitialConditions(vector<double> pos_hyd, vector<double> hyd_inp,
                                  vector<double> pos_temp,vector<double> temp_inp);
    
    // Get the current temperature from the main program
        void getInput(vector<double> pos_temp, vector<double> temp_inp);
    
    // Compute the evolution of the system during a time step
        void compute() ;
    
    // Compute the time step based on the kinetics of each penomenon
        void computeTimeStep();
    
    // Getters
        Sample*      returnSample   () {return _sample      ;} ;
        Diffusion*   returnDiff     () {return _diffusion   ;} ;
        Nucleation*  returnNuc      () {return _nucleation  ;} ;
        Growth *     returnGro      () {return _growth      ;} ;
        Dissolution* returnDiss     () {return _dissolution ;} ;
        
        double returnTimeStep() {return _dt;};
        
    private:
        Sample*      _sample      ; // Geometry, temperature and solubility management
        Diffusion*   _diffusion   ;
        Nucleation*  _nucleation  ;
        Growth*      _growth      ;
        Dissolution* _dissolution ;
    
        int _NbCells ;  // number of nodes in the geometry
        int _geometry ;    // Geometry type
        double _radius ;    // Radius or Sample Length
        bool _auto_dt ; // time step fixed by user or computed at each step
        double _dt ;    // time step value

        const vector<double> * _position ;  // List of positions
        const vector<double> * _Css ;       // Solid soltution profile
        const vector<double> * _Ctot ;      // Hydrogen profile
        const vector<double> * _Cprec ;     // Hydride profile
        const vector<double> * _tssp ;      // Supersolubility profile
        const vector<double> * _tssd ;      // Solubility profile

        const vector<double> * _flux ;      // Hydrogen flux
        
        const vector<double> * _rateNuc ;   // Rate of nucleation at each position
        const vector<double> * _rateGro ;   // Rate of growth at each position
        const vector<double> * _rateDis ;   // Rate of dissolution at each position
    
};



#endif /* HNGD_hpp */
