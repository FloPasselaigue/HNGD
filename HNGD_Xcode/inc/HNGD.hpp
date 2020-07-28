#ifndef HNGD_hpp
#define HNGD_hpp

#include <stdio.h>
#include "Sample.hpp"
#include "Diffusion.hpp"
#include "Nucleation.hpp"
#include "Growth.hpp"
#include "Dissolution.hpp"

class HNGD
{
    public:
        HNGD(double* settings, double* physicalParameters);
        
        void getInitialConditions(vector<double> pos_hyd, vector<double> hyd_inp,
                                  vector<double> pos_temp,vector<double> temp_inp);
    
        void getInput(vector<double> pos_temp, vector<double> temp_inp);
    
        void compute() ;
    
        void computeTimeStep();
    
    Sample*      returnSample   () {return _sample      ;} ;
    Diffusion*   returnDiff     () {return _diffusion   ;} ;
    Nucleation*  returnNuc      () {return _nucleation  ;} ;
    Growth *     returnGro      () {return _growth      ;} ;
    Dissolution* returnDiss     () {return _dissolution ;} ;
    
    double returnTimeStep() {return _dt;};
        
    private:
        Sample*      _sample      ;
        Diffusion*   _diffusion   ;
        Nucleation*  _nucleation  ;
        Growth*      _growth      ;
        Dissolution* _dissolution ;
    
        int _NbCells ;
        bool _auto_dt ;
        double _dt ;
    
};



#endif /* HNGD_hpp */
