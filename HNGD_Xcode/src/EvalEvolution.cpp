#include "EvalEvolution.hpp"
#include <math.h>

EvalEvolution::EvalEvolution()
{
    _mem_profile = vector<double>(0) ;
}

EvalEvolution::EvalEvolution(double criterion)
{
    _criterion = criterion ;
    _mem_profile = vector<double>(0) ;
}

bool EvalEvolution::evaluate(vector<double> profile)
{
    // Compute average relative difference between the
    // parameter profile and the memorized one
    double av = 0 ;
    
    for(int k=0; k<profile.size(); k++)
        av += abs(_mem_profile[k] - profile[k]) / _mem_profile[k] ;
    
    av /= profile.size() ;
    
    // If the average evolution is low, no printing
    if(av < _criterion)
        return false ;
    
    // Else, the current profile replaces the memorized
    else
    {
        _mem_profile = vector<double>(profile) ;
        return true ;
    }
    
}

void EvalEvolution::setCriterion(double value)
{
    _criterion = value ;
}

void EvalEvolution::setProfile(vector<double> profile)
{
    _mem_profile = vector<double>(profile);
}
