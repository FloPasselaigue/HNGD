#ifndef Precipitation_hpp
#define Precipitation_hpp

#include <stdio.h>
#include "Mechanism.hpp"
#include "PhysicsConstants.h"

using namespace std ;

class Precipitation : public Mechanism
{
    public:
        Precipitation(Sample* sample) ;
        static void computeCoeffs(Sample* sample) ;
        
    protected:
    
        static vector<double> _lever_rule ;
        static vector<double> _f_alpha ;
        static vector<double> _chi ;
        
        static double compute_xdelta(double T) ;
        static double convertToAtomFrac(double c) ;
        
        static void compute_leverRule(Sample* sample) ;
        static void compute_f_alpha(Sample* sample) ;
        static void compute_chi(Sample* sample) ;
    
    private:
        static vector<double> * _totalContent ;
        static vector<double> * _hydrideContent ;
        static vector<double> * _tssd ;
        static vector<double> * _temperature ;
};

#endif /* Precipitation_hpp */
