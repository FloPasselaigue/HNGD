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
    
    void computeTSS();
    void computeEquilibrium();
    void computeLocations(double x0, double xEnd);
    void spatialeInterpolation(vector<double>& refX, vector<double>& refY, vector<double>& vectorY);
    
    const int returnNbCells() {return _nbCells;}

    vector<double>& returnPosition        () {return _position;}
    vector<double>& returnTemperature     () {return _temperature;}
    vector<double>& returnSolutionContent () {return _solutionContent;}
    vector<double>& returnHydrideContent  () {return _hydrideContent;}
    vector<double>& returnTotalContent    () {return _totalContent;}
    vector<double>& returnTSSp            () {return _tssp;}
    vector<double>& returnTSSd            () {return _tssd;}
    
    void setTemperature     (vector<double> temperature){_temperature     = vector<double>(temperature) ;}
    void setSolutionContent (vector<double> c_ss)       {_solutionContent = vector<double>(c_ss)        ;}
    void setHydrideContent  (vector<double> c_prec)     {_hydrideContent  = vector<double>(c_prec)      ;}
    
    void updateTotalContent () {for(int k=0; k<_nbCells; k++)
        _totalContent[k] = _solutionContent[k]+_hydrideContent[k] ;}
    
private:
    const int _nbCells ;
    const double _bias ;
    
    vector<double> _position ;
    vector<double> _temperature ;
    vector<double> _solutionContent ;
    vector<double> _hydrideContent ;
    vector<double> _totalContent ;
    
    const double _tssp0, _Qp ;
    const double _tssd0, _Qd ;
    vector<double> _tssp ;
    vector<double> _tssd ;
    
};



#endif /* Sample_hpp */
