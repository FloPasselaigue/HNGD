#include "Precipitation.hpp"

vector<double> Precipitation::_lever_rule(0);
vector<double> Precipitation::_f_alpha(0);
vector<double> Precipitation::_chi(0);

vector<double> * Precipitation::_hydrideContent = new vector<double>(0) ;
vector<double> * Precipitation::_totalContent = new vector<double>(0) ;
vector<double> * Precipitation::_temperature = new vector<double>(0) ;
vector<double> * Precipitation::_tssd = new vector<double>(0) ;

// Constructor
Precipitation :: Precipitation(Sample* sample):
    Mechanism(sample)
{
    _lever_rule.resize(sample->returnNbCells()) ;
    _f_alpha.resize(sample->returnNbCells()) ;
    _chi.resize(sample->returnNbCells()) ;
    
    _hydrideContent = & (sample->returnHydrideContent()) ;
    _totalContent = & (sample->returnTotalContent()) ;
    _temperature = & (sample->returnTemperature()) ;
    _tssd = & (sample->returnTSSd()) ;
}

void Precipitation :: computeCoeffs(Sample* sample)
{
    compute_leverRule(sample) ;
    compute_f_alpha(sample) ;
    compute_chi(sample) ;
}

double Precipitation :: compute_xdelta(double T)
{
    // 3rd degree polynome fitting
    return -9.93e-11*pow(T,3) + 8.48e-8*pow(T,2) - 5.73e-5*T + 0.623 ;
}

double Precipitation :: convertToAtomFrac(double c)
{
    return c / (Mh * (c/Mh + (1e6 - c)/Mzr)) ;
}

void Precipitation :: compute_f_alpha(Sample* sample)
{
    for(int k=0; k<sample->returnNbCells(); k++)
    {
        // Hydride atom fraction
        double xhyd = convertToAtomFrac((*_hydrideContent)[k]) ;
        
        // alpha / alpha+delta boundaries
        double xdelta = compute_xdelta((*_temperature)[k]) ;
        
        _f_alpha[k] = 1. - xhyd / xdelta ;
    }
    
}


void Precipitation :: compute_leverRule(Sample* sample)
{
    for(int k=0; k<sample->returnNbCells(); k++)
    {
        // Total atomic fraction of hydrogen
        double xtot = convertToAtomFrac((*_totalContent)[k]);
        
        // TSSd atomic fraction
        double xalpha = convertToAtomFrac((*_tssd)[k]) ;
        
        // alpha / alpha+delta boundaries
        double xdelta = compute_xdelta((*_temperature)[k]) ;
        
        _lever_rule[k] = ((xtot - xalpha)/(xdelta - xalpha));
    }
    
}

void Precipitation :: compute_chi(Sample* sample)
{
    for(int k=0; k<sample->returnNbCells(); k++)
        _chi[k] = ((17000 - (*_hydrideContent)[k]) / 17000) ;
}

