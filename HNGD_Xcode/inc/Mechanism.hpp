/**
    This is a superclass for hydride nucleation/growth/dissolution.
    It defines the quantities (kinetics, driving force, rate) that exist
    in each of these phenomena. Daughter classes specifies how
    these quantities are computed.
 */

#ifndef Mechanism_hpp
#define Mechanism_hpp

#include <stdio.h>
#include <vector>
#include "Sample.hpp"
#include "PhysicsConstants.h"

using namespace std ;

class Mechanism
{
    public:
        Mechanism(Sample* sample) ;
    
    // Compute the rate of the mechanism based on the kinetics and driving force
        void computeRate() ;
        
    // Compute the time step associated with the mechanism based on the kinetics
        double timeStep() ;
    
    // Getters
        vector<double> & returnKinetics();
        vector<double> & returnDrivForce();
        vector<double> & returnRate();
        
    protected:
        const int _nbCells ; // Number of nodes
    
        vector<double> _kinetic_factor ;
        vector<double> _driving_force ;
        vector<double> _rate ;

    // Virtual functions to be redefined for each mechanism
    // to compute the kinetics and driving force
        virtual void computeKinetics() = 0 ;
        virtual void computeDrivForce()= 0 ;
    
};

#endif /* Mechanism_hpp */
