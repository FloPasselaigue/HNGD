#ifndef HydrogenBehaviorModel_H
#define HydrogenBehaviorModel_H

#include <vector>

const double kb     = 8.617e-5  ; //Boltzmann constant (eV/K)
const double R      = 8.314     ; //ideal gas constant (J/mol/K)
const double e      = 1.602e-19 ; //elementary charge (J/eV)
const double Mh     = 1.00784   ; //molar mass of Hydrogen (g/mol)
const double Mzr    = 91.224    ; //molar mass of zirconium (g/mol)
const double rhoZr  = 6.52      ; //volumic mass of zirconium (g/cm3)
const double rhoHyd = 5.66      ; //volumic mass of hydrides  (g/cm3)

using namespace std;

class HydrogenBehaviorModel
{
  public:
    HydrogenBehaviorModel         ();

    void getInitialConditions     (double* settings,double* physicalParameters,
      vector<double> pos_hyd, vector<double> hyd_inp, vector<double> pos_temp, vector<double> temp_inp);

    void getInput                 (double time, vector<double> pos_temp, vector<double> temp_inp, double dt);

    void computeProperties        ();
    void computeLocations         (double x0, double xEnd, double bias, std::vector<double>& x);

    double returnAdvancement      ();
    double returnTotalContent     ();
    double returnSolutionContent  ();
    double returnHydridesContent  ();
    double returnPosition         ();
    double returnTimeStep         ();

    std::vector<double>& returnTSSpVector             ();
    std::vector<double>& returnTSSdVector             ();
    std::vector<double>& returnTotalContentVector     ();
    std::vector<double>& returnSolutionContentVector  ();
    std::vector<double>& returnHydridesContentVector  ();
    std::vector<double>& returnKdVector               ();
    std::vector<double>& returnKgVector               ();
    std::vector<double>& returnKmobVector             ();
    std::vector<double>& returnKthVector              ();
    std::vector<double>& returnfVector                ();
    std::vector<double>& returnLeverVector            ();
    std::vector<double>& returnKnVector               ();
    std::vector<double>& returnPositionVector         ();
    std::vector<double>& returnTemperatureVector      ();
    std::vector<double>& returnCoeffFickVector        ();
    std::vector<double>& returnFlux                   ();

    void print();//Auxiliary method to find what's wrong

  private:
    double    _time                               ; // time (s)

    // Options
    short int _solubilities                       ; // Linear (1) or exponential (2) dependency
    short int _kinetics                           ; // Equilibrium calculation (0) ; contant kinetics (1) ; exponential dependency (3)
    short int _optionNGD                          ; // No precipitation/dissolution (0) ; only nucleation (1) ; nucleation + dissolution (2) ; nucleation + growth + dissolution (3)
    short int _fickLaw                            ; // No diffusion (0) ; constant diffusion coefficient (1) ; exponential dependency (2)
    short int _soretEffect                        ; // No Soret effect (0) ; with Soret effect (1)
    double    _bias                               ; // Bias for the spatial discretization
    double    _sampleLenght                       ; // Sample lenght in cm
    double    _dt                                 ; // Time step (s)
    int       _nbNodes                            ; // Number of nodes
    bool      _nuclIsHappening                    ; // Boolean to know if nucleation kinetics need to be taken into account to compute the time step
    bool      _growIsHappening                    ; // Boolean to know if growth kinetics need to be taken into account to compute the time step
    bool      _dissIsHappening                    ; // Boolean to know if dissolution kinetics need to be taken into account to compute the time step

    // Physical quantities
        // Space dependent
    std::vector<double> _positionsVector          ; // List of the node positions
    std::vector<double> _temperatureVector        ; // Temperature profile
    std::vector<double> _solutionContentVector    ; // Solid solution profile
    std::vector<double> _hydridesContentVector    ; // Hydrides profile
    std::vector<double> _totalContentVector       ; // Total hydrogen content profile
    std::vector<double> _flux                     ; // Hydrogen flux


    // Physical parameters
    double    _p                                  ; // Avrami parameter (dimentionnality)
    double    _Kn0, _Kn                           ; // Nucleation
    double    _Kd0, _Ediss, _Kd                   ; // Dissolution
    double    _tssp0, _Qtssp, _tssp               ; // Super solubility
    double    _tssd0, _Qtssd, _tssd               ; // Solubility
    double    _Kth0, _Kmob0, _Eg, _Kg             ; // Growth
    double    _Eth0, _Eth1, _Eth2, _Eth3          ; // Hydride formation energy
    double    _D0, _Ediff, _Qstar                 ; // Diffusion

    std::vector<double> _tsspVector               ; // Supersolubility
    std::vector<double> _tssdVector               ; // Solubility
    std::vector<double> _KdVector                 ; // Dissolution kinetics
    std::vector<double> _KnVector                 ; // Nucleation kinetics
    std::vector<double> _KmobVector               ; // Growth kinetics
    std::vector<double> _KthVector                ; // Growth kinetics
    std::vector<double> _f_alpha                  ; // Growth kinetics
    std::vector<double> _lever                    ; // Growth kinetics
    std::vector<double> _KgVector                 ; // Growth kinetics
    std::vector<double> _coeffFickVector          ; // Fick's law diffusion coefficient
    

    // Internal methods
    void      computeInitEq             () ; // Compute the initial equilibrium state

    void      computePhysicalParameters () ; // Compute the solubilities, kinetics, diffusion coefficient
    double    computeNucleationRate     () ;
    double    computeGrowthRate         (double delta_solutionContent) ;
    double    computeDissolutionRate    () ;
    double    computeNucleationEq       () ;
    double    computePrecipitationEq    () ;
    double    computeDissolutionEq      () ;
    double    leverRule                 () ;
    double    factor_f_alpha            () ;
    double    computeTimeStep           () ;

    void      computePhysicalParameters (int position) ;
    double    computeNucleationRate     (int position) ;
    double    computeGrowthRate         (double delta_solutionContent,int position) ;
    double    computeDissolutionRate    (int position) ;
    double    computeNucleationEq       (int position) ;
    double    computePrecipitationEq    (int position) ;
    double    computeDissolutionEq      (int position) ;
    double    leverRule                 (int position) ;
    double    factor_f_alpha            (int position) ;
    void      computeDiffusionFlux      (int position) ;
    void      spatialeInterpolation(std::vector<double>& refX, std::vector<double>& refY, std::vector<double>& vectorX, std::vector<double>& vectorY);

};

#endif
