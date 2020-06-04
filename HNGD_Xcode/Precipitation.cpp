#include "Precipitation.hpp"

// Static membres initialization
double Precipitation::_Eth0(0.) ;
double Precipitation::_Eth1(0.) ;
double Precipitation::_Eth2(0.) ;
double Precipitation::_Eth3(0.) ;

vector<double> Precipitation::_lever_rule(0);
vector<double> Precipitation::_f_alpha(0);
vector<double> Precipitation::_chi(0);

// Constructor
Precipitation :: Precipitation(Sample* sample):
    Mechanism(sample)
{
    _lever_rule.resize(sample->returnNbCells()) ;
    _f_alpha.resize(sample->returnNbCells()) ;
    _chi.resize(sample->returnNbCells()) ;
}

// Formation energy fit
void Precipitation :: defineEnergyPolynomial(double Eth0, double Eth1, double Eth2, double Eth3)
{
    _Eth0 = Eth0 ;
    _Eth1 = Eth1 ;
    _Eth2 = Eth2 ;
    _Eth3 = Eth3 ;
}

void Precipitation :: computeCoeffs(Sample* sample)
{
    compute_leverRule(sample) ;
    compute_f_alpha(sample) ;
    compute_chi(sample) ;
}

double Precipitation :: formation_energy(double T)
{
    return -_Eth0 + _Eth1 * T - _Eth2 * pow(T,2) + _Eth3 * pow(T,3);
}

double Precipitation :: compute_xdelta(double T)
{
    // 3rd degree polynome fitting
    return -9.93e-11*pow(T,3) + 8.48e-8*pow(T,2) - 5.73e-5*T + 0.623 ;
}

double Precipitation :: convertToAtomFrac(double c)
{
    return c / (Mh * (c/Mh + (1e6 - c/Mzr))) ;
}

void Precipitation :: compute_f_alpha(Sample* sample)
{
    vector<double> hydrideContent = sample->returnHydrideContent() ;
    vector<double> tssd = sample->returnTSSd() ;
    vector<double> temperature = sample->returnTemperature() ;
    
    
    for(int k=0; k<sample->returnNbCells(); k++)
    {
        // Hydride atom fraction
        double xhyd = convertToAtomFrac(hydrideContent[k]) ;
        
        // alpha / alpha+delta boundaries
        double xdelta = compute_xdelta(temperature[k]) ;
        
        // TSSd atom fraction
        double xalpha = convertToAtomFrac(tssd[k]) ;
        
        _f_alpha[k] = 1. - xhyd / xdelta ;
//        f[k] = 1. - xhyd/(xdelta - xalpha) ;
    }
    
}


void Precipitation :: compute_leverRule(Sample* sample)
{
    
    
    vector<double> totalContent = sample->returnTotalContent() ;
    vector<double> tssd = sample->returnTSSd() ;
    vector<double> temperature = sample->returnTemperature() ;
    
    for(int k=0; k<sample->returnNbCells(); k++)
    {
        // Total atomic fraction of hydrogen
        double xtot = convertToAtomFrac(totalContent[k]);
        
        // TSSd atomic fraction
        double xalpha = convertToAtomFrac(tssd[k]) ;
        
        // alpha / alpha+delta boundaries
        double xdelta = compute_xdelta(temperature[k]) ;
        
        _lever_rule[k] = ((xtot - xalpha)/(xdelta - xalpha));
    }
    
}

void Precipitation :: compute_chi(Sample* sample)
{
    vector<double> hydrideContent = sample->returnHydrideContent() ;
    
    for(int k=0; k<sample->returnNbCells(); k++)
        _chi[k] = ((17000 - hydrideContent[k]) / 17000) ;
    
}

