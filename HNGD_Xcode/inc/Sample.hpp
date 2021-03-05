/**
    This class creates the geometry, manages the temperature profile,
    and compute the solubility/supersolubility profiles for the HNGD model.
 */

#ifndef Sample_hpp
#define Sample_hpp

#include <stdio.h>
#include <vector>
#include <math.h>

using namespace std ;

class Sample
{
    
public:
    Sample(int nbCells, double bias, double tssp0, double Qp, double tssd0, double Qd);
    
    // Compute the solubility/supersolubility profiles
    void computeTSS();
    
    // Compute the equilibrium state for the initial condition
    void computeEquilibrium();
    
    // Create the geometry
    void computeLocations(double x0, double xEnd, int geometry);
    
    // Interpolate the input profile [refX; refY] on each point of the geometry
    void spatialeInterpolation(vector<double>& refX, vector<double>& refY, vector<double>& vectorY);
    
    // Getters
    const int returnNbCells() {return _nbCells;}
    vector<double>& returnPosition        () {return _position;}
    vector<double>& returnTemperature     () {return _temperature;}
    vector<double>& returnSolutionContent () {return _solutionContent;}
    vector<double>& returnHydrideContent  () {return _hydrideContent;}
    vector<double>& returnTotalContent    () {return _totalContent;}
    vector<double>& returnTSSp            () {return _tssp;}
    vector<double>& returnTSSd            () {return _tssd;}
    
    // Setters
    void setTemperature     (vector<double> temperature){_temperature     = vector<double>(temperature) ;}
    void setSolutionContent (vector<double> c_ss)       {_solutionContent = vector<double>(c_ss)        ;}
    void setHydrideContent  (vector<double> c_prec)     {_hydrideContent  = vector<double>(c_prec)      ;}
    
    void updateTotalContent () {for(int k=0; k<_nbCells; k++)
        _totalContent[k] = _solutionContent[k]+_hydrideContent[k] ;}
    
private:
    const int _nbCells ; // Number of nodes
    const double _bias ; // Bias of the geometry
    
    vector<double> _position ;          // Positions
    vector<double> _temperature ;       // Temperature profile
    vector<double> _solutionContent ;   // Solid solution profile
    vector<double> _hydrideContent ;    // Hydride profile
    vector<double> _totalContent ;      // Hydrogen profile
    
    const double _tssp0, _Qp ;  // Preexponential factor and activation energy for supersolubility
    const double _tssd0, _Qd ;  // Preexponential factor and activation energy for solubility
    
    vector<double> _tssp ;  // Supersolubility profile
    vector<double> _tssd ;  // Solubility profile
    
    vector<double> _dC_dx; // TODO: why are gradients needed here?
    vector<double> _dT_dx;
    
};



#endif /* Sample_hpp */
