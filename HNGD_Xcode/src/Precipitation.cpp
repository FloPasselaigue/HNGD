#include "Precipitation.hpp"

// Static members initialization

double Precipitation::_Eth0(0.) ;
double Precipitation::_Eth1(0.) ;
double Precipitation::_Eth2(0.) ;
double Precipitation::_Eth3(0.) ;

vector<double> Precipitation::_lever_rule(0);
vector<double> Precipitation::_f_alpha(0);

vector<double> * Precipitation::_hydrideContent = new vector<double>(0) ;
vector<double> * Precipitation::_totalContent   = new vector<double>(0) ;
vector<double> * Precipitation::_temperature    = new vector<double>(0) ;
vector<double> * Precipitation::_tssd           = new vector<double>(0) ;

// Constructor
Precipitation :: Precipitation(Sample* sample):
    Mechanism(sample)
{
    _lever_rule .resize(sample->returnNbCells()) ;
    _f_alpha    .resize(sample->returnNbCells()) ;
    
    _hydrideContent = & (sample->returnHydrideContent()) ;
    _totalContent   = & (sample->returnTotalContent()) ;
    _temperature    = & (sample->returnTemperature()) ;
    _tssd           = & (sample->returnTSSd()) ;
}

void Precipitation :: computeCoeffs(Sample* sample)
{
    compute_leverRule(sample) ;
    compute_f_alpha(sample) ;
}

// Formation energy fit
void Precipitation :: defineEnergyPolynomial(double Eth0, double Eth1, double Eth2, double Eth3)
{
    _Eth0 = Eth0 ;
    _Eth1 = Eth1 ;
    _Eth2 = Eth2 ;
    _Eth3 = Eth3 ;
}

double Precipitation :: compute_xdelta(double T)
{
    // 3rd degree polynomial fit
    return -9.93e-11*pow(T,3) + 8.48e-8*pow(T,2) - 5.73e-5*T + 0.623 ;
}

double Precipitation :: convertToAtomFrac(double c)
{
    return c / (Mh * (c/Mh + (1e6 - c)/Mzr)) ;
}

double Precipitation :: formation_energy(double T)
{
    return -_Eth0 + _Eth1 * T - _Eth2 * pow(T,2) + _Eth3 * pow(T,3);
}

void Precipitation :: compute_f_alpha(Sample* sample)
{
    for(int k=0; k<sample->returnNbCells(); k++)
    {
        // Hydride atom fraction
//        double xhyd = (*_hydrideContent)[k] / (Mh * ((*_totalContent)[k]/Mh + (1e6 - (*_totalContent)[k])/Mzr)) ;
        double xhyd = convertToAtomFrac((*_hydrideContent)[k]) ;
        
        // alpha / alpha+delta boundaries
        double xdelta = compute_xdelta((*_temperature)[k]) ;
        
        _f_alpha[k] = 1. - xhyd / (xdelta - convertToAtomFrac((*_tssd)[k])) ;
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
