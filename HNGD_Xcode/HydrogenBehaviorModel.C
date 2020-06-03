// Implementation of the class that contains the hydrogen behavior model.

#include "HydrogenBehaviorModel.h"
#include <iostream>
#include <math.h>
#include <vector>

// Constructor
HydrogenBehaviorModel :: HydrogenBehaviorModel()
{
  _solubilities    = 0  ;
  _kinetics        = 0  ;
  _optionNGD       = 0  ;
  _nbNodes         = 0  ;
  _fickLaw         = 0  ;
  _nuclIsHappening = false;
  _growIsHappening = false;
  _dissIsHappening = false;

  _positionsVector       =  std::vector<double>(0.) ;
  _solutionContentVector =  std::vector<double>(0.) ;
  _hydridesContentVector =  std::vector<double>(0.) ;
  _totalContentVector    =  std::vector<double>(0.) ;
  _temperatureVector     =  std::vector<double>(0.) ;
  _tsspVector            =  std::vector<double>(0.) ;
  _tssdVector            =  std::vector<double>(0.) ;
  _KdVector              =  std::vector<double>(0.) ;
  _KnVector              =  std::vector<double>(0.) ;
  _KgVector              =  std::vector<double>(0.) ;
  _f_alpha               =  std::vector<double>(0.) ;
  _lever                 =  std::vector<double>(0.) ;
  _coeffFickVector       =  std::vector<double>(0.) ;
  _flux                  =  std::vector<double>(0.) ;
  _probaVector           =  std::vector<double>(0.) ;

  _D0     = 0. ;
  _Ediff  = 0. ;
  _Qstar  = 0. ;
  _Kd0    = 0. ;
  _Ediss  = 0. ;
  _Kn0    = 0. ;
  _Eth0   = 0. ;
  _Eth1   = 0. ;
  _Eth2   = 0. ;
  _Eth3   = 0. ;
  _Kmob0  = 0. ;
  _Kth0   = 0. ;
  _Eg     = 0. ;
  _tssp0  = 0. ;
  _Qtssp  = 0. ;
  _tssd0  = 0. ;
  _Qtssd  = 0. ;
  _dt	    = 1. ;
  _p      = 2.5;

}

// Define the initial conditions
void HydrogenBehaviorModel :: getInitialConditions(double* settings,double* physicalParameters,
  vector<double> pos_hyd, vector<double> hyd_inp, vector<double> pos_temp, vector<double> temp_inp)
{
    _optionNGD    = (int)settings[1] ;
    _nbNodes      = (int)settings[2] ;
    _kinetics     = (int)settings[3] ;
    _solubilities = (int)settings[4] ;
    _fickLaw      = (int)settings[5] ;
    _soretEffect  = (int)settings[6] ;
    _bias         =      settings[7] ;
    _sampleLenght =      settings[8] ;
    _dt           =      settings[9] ;

    _Kd0    = physicalParameters[0]  ;
    _Ediss  = physicalParameters[1]  ;
    _Kn0    = physicalParameters[2]  ;
    _Eth0   = physicalParameters[3]  ;
    _Eth1   = physicalParameters[4]  ;
    _Eth2   = physicalParameters[5]  ;
    _Eth3   = physicalParameters[6]  ;
    _Kmob0  = physicalParameters[7]  ;
    _Kth0   = physicalParameters[8]  ;
    _Eg     = physicalParameters[9]  ;
    _tssp0  = physicalParameters[10] ;
    _Qtssp  = physicalParameters[11] ;
    _tssd0  = physicalParameters[12] ;
    _Qtssd  = physicalParameters[13] ;
    _D0     = physicalParameters[14] ;
    _Ediff  = physicalParameters[15] ;
    _Qstar  = physicalParameters[16] ;

    _positionsVector      .resize(_nbNodes);
    _totalContentVector   .resize(_nbNodes);
    _solutionContentVector.resize(_nbNodes);
    _temperatureVector    .resize(_nbNodes);
    _hydridesContentVector.resize(_nbNodes,0.);
    _tsspVector           .resize(_nbNodes,0.);
    _tssdVector           .resize(_nbNodes,0.);
    _KdVector             .resize(_nbNodes,0.);
    _KnVector             .resize(_nbNodes,0.);
    _KgVector             .resize(_nbNodes,0.);
    _KmobVector           .resize(_nbNodes,0.);
    _KthVector            .resize(_nbNodes,0.);
    _f_alpha              .resize(_nbNodes,0.);
    _lever                .resize(_nbNodes,0.);
    _coeffFickVector      .resize(_nbNodes,0.);
    _flux                 .resize(_nbNodes,0.);
    _probaVector          .resize(_nbNodes,1.);

    HydrogenBehaviorModel::computeLocations(0., _sampleLenght, _bias, _positionsVector);

    spatialeInterpolation(pos_hyd,  hyd_inp,  _positionsVector, _totalContentVector); // Hydrogen content profile
    spatialeInterpolation(pos_hyd,  hyd_inp,  _positionsVector, _solutionContentVector);
    spatialeInterpolation(pos_temp, temp_inp, _positionsVector, _temperatureVector); // Temperature profile

    if(settings[0]==1)
      computeInitEq();

    if(_dt<=0) // adaptative time step
      _dt = computeTimeStep() ;

}

// Overriding of the getInput method
void HydrogenBehaviorModel :: getInput(double time, vector<double> pos_temp, vector<double> temp_inp, double dt)
{
  _time = time;
  _dt   = dt ;

  spatialeInterpolation(pos_temp, temp_inp, _positionsVector, _temperatureVector);
}

// We assume the system is at equilibrium when the simulation starts
void HydrogenBehaviorModel :: computeInitEq()
{
  for(int k=0; k<_nbNodes; k++)
    computePhysicalParameters(k);

  if(_optionNGD==0)
    return;

  for(int k=0; k<_nbNodes; k++)
  {
    if(_optionNGD<3)
    {
        _solutionContentVector[k] = std::min(_totalContentVector[k], _tsspVector[k]);
        _hydridesContentVector[k] = _totalContentVector[k] - _solutionContentVector[k];
    }
    else
    {
        _solutionContentVector[k] = std::min(_totalContentVector[k], _tssdVector[k]);
        _hydridesContentVector[k] = _totalContentVector[k] - _solutionContentVector[k];
    }
  }
}

// Overriding of the computeProperties method
void HydrogenBehaviorModel :: computeProperties()
{
  
// Physical Parameters update (solubilities, kinetics, diffusion coefficient)
    for(int k=0; k<_nbNodes; k++)
        computePhysicalParameters(k) ;

// Time step update
    if(_dt<0)
      _dt = computeTimeStep() ;

// Solid solution distribution computation
    if(_fickLaw || _soretEffect)
    {
        // Hydrogen Flux
        for(int k=0; k<_nbNodes; k++)
              computeDiffusionFlux(k);

        // Diffusion
        std::vector<double> newSolutionContent(_nbNodes,0.);
        newSolutionContent[0] = _solutionContentVector[0] - _dt*(_flux[0])/(_positionsVector[1]-_positionsVector[0]) ;
        for(int k=1; k<_nbNodes; k++)
          newSolutionContent[k] = _solutionContentVector[k] - _dt*(_flux[k] - _flux[k-1])/(_positionsVector[k]-_positionsVector[k-1]);

        // Update of solid solution and total content
        for(int k=0; k<_nbNodes; k++)
        {
          _solutionContentVector[k] = newSolutionContent[k];
          _totalContentVector[k] = _solutionContentVector[k] + _hydridesContentVector[k];
        }

    }

// Nucleation/Growth/Dissolution
    if(_optionNGD>0)
    {
      for(int k=0; k<_nbNodes; k++)
      {
        double delta_solutionContent(0.) ;

        if(_kinetics==1) // Equilibrium calculation
        {
          //Nucleation if Css above TSSp
          if ((_optionNGD>0) && (_solutionContentVector[k] > _tsspVector[k]))
              delta_solutionContent = computeNucleationEq(k) ;

          //Dissolution if hydrides under TSSd
          if ((_optionNGD>1) && (_solutionContentVector[k] < _tssdVector[k] && _hydridesContentVector[k]>0.))
              delta_solutionContent = computeDissolutionEq(k) ;

          // Precipitation if Css above TSSp  OR  if hydrides>0 with Css above TSSd
          if ((_optionNGD>2) && (_solutionContentVector[k] > _tsspVector[k] || (_solutionContentVector[k]>_tssdVector[k] && _hydridesContentVector[k]>0.)))
              delta_solutionContent = computePrecipitationEq(k) ;

          _solutionContentVector[k] += delta_solutionContent ;
          _hydridesContentVector[k] -= delta_solutionContent ;
        }

        else
        {
//          // Nucleation over TSSp
//          if ((_optionNGD>0) && (_solutionContentVector[k] >= _tsspVector[k]))
//            delta_solutionContent += computeNucleationRate(k) ;
          
          // Nucleation always computed
          delta_solutionContent += computeNucleationRate(k) ;

          // Dissolution if hydrides under TSSd
          if ((_optionNGD>1) && (_solutionContentVector[k] < _tssdVector[k] && _hydridesContentVector[k]>0.))
            delta_solutionContent += computeDissolutionRate(k) ;

          // Growth if hydrides above TSSd
          if ((_optionNGD>2) && (_solutionContentVector[k] > _tssdVector[k] && (_hydridesContentVector[k]>0. || delta_solutionContent<0)))
            delta_solutionContent += computeGrowthRate(delta_solutionContent,k) ;

          if(_hydridesContentVector[k] - delta_solutionContent*_dt < 0)
            cout << 'p' ;
          
          _solutionContentVector[k] += delta_solutionContent*_dt ;
          _hydridesContentVector[k] -= delta_solutionContent*_dt ;

        }
      }
    }
} //end computeProperties


// ------------------------------------------ Space computation   ------------------------------------------
    // Diffusion
void HydrogenBehaviorModel :: computeDiffusionFlux(int position)
{
  double fluxFick(0.),fluxSoret(0.) ;
  double dC_dx(0.), dT_dx=(0.) ; ;
  
  if(_fickLaw)
  {
    if (!(position==_nbNodes-1))
      dC_dx = (_solutionContentVector[position+1]-_solutionContentVector[position])/(_positionsVector[position+1] - _positionsVector[position]);

    fluxFick = - _coeffFickVector[position] * dC_dx ;
  }

  if(_soretEffect)
  {
    if (!(position==_nbNodes-1))
      dT_dx = (_temperatureVector[position+1]-_temperatureVector[position])/(_positionsVector[position+1] - _positionsVector[position]);

    fluxSoret = - _coeffFickVector[position]*_Qstar*_solutionContentVector[position]* dT_dx/(R*pow(_temperatureVector[position],2))  ;
  }

  _flux[position] = fluxFick+fluxSoret ;
}

    // NGD equilibrium
double HydrogenBehaviorModel :: computeNucleationEq   (int position) {
  return -(_solutionContentVector[position] - _tsspVector[position]) ;}
double HydrogenBehaviorModel :: computePrecipitationEq(int position) {
  return -(_solutionContentVector[position] - _tssdVector[position]) ;}
double HydrogenBehaviorModel :: computeDissolutionEq  (int position) {
  return -(std::max(_solutionContentVector[position] - _tssdVector[position],
                    _solutionContentVector[position] - _totalContentVector[position])) ;}

    // NGD kinetics
double HydrogenBehaviorModel :: computeNucleationRate(int position)
{
  // if the timestep is too large it is re computed
  if(1./_KnVector[position] < _dt)
    std::cout << "Warning : time step too big for Nucleation: 1/Kn=" << 1./_KnVector[position] << " s " << "t="<<_time<<" dt="<<_dt<<'\n';

  //max to not precipitate more than equilibrium
  return -_KnVector[position] * _solutionContentVector[position] ;
  
//  return std::max(-_KnVector[position]*(_solutionContentVector[position] - _tsspVector[position]), -(1./_dt)*(_solutionContentVector[position] - _tsspVector[position]) );
}

double HydrogenBehaviorModel :: computeGrowthRate(double delta_solutionContent,int position)
{
  // if the timestep is too large it is re computed
  if(_KgVector[position]>0 && 1./_KgVector[position] < _dt)
    std::cout << "Warning : time step too big for Growth: 1/Kg=" << 1./_KgVector[position]<<" dt="<<_dt << " s " << "t="<<_time<<'\n';

  // Reaction advancement
  double x = (_hydridesContentVector[position]-delta_solutionContent) / (_totalContentVector[position]-_tssdVector[position]) ;
  if(x>=1.)
    return 0. ;
  else
  {
    double delta = (_totalContentVector[position] - _tssdVector[position])*_p*(1-x)*pow(-log(1-x),1-1./_p) ;
    return std::max(-_KgVector[position]*delta, -(1./_dt)*delta) ;
  }
}

double HydrogenBehaviorModel :: computeDissolutionRate(int position)
{
  // if the timestep is too large it is re computed
  if(1./_KdVector[position]<_dt)
    std::cout << "Warning : time step too big for Dissolution: 1/Kd=" << 1./_KdVector[position] << " s " << "t="<<_time<<" dt="<<_dt<<'\n';

  //min to not dissolve more than equilibrium
  double delta = std::max(_solutionContentVector[position] - _tssdVector[position], _solutionContentVector[position] - _totalContentVector[position] );
  return std::min(-_KdVector[position]*delta, -(1./_dt)*delta) ;
}

    // Physical parameters
void HydrogenBehaviorModel :: computePhysicalParameters(int position)
{
  // Solubilities
  switch(_solubilities)
  {
    case 0: // Linear dependence
      _tsspVector[position] = std::max(0.,_Qtssp*_temperatureVector[position] - _tssp0);
      _tssdVector[position] = std::max(0.,_Qtssd*_temperatureVector[position] - _tssd0);
      break ;

    case 1: // Exponential dependence
      _tsspVector[position] = _tssp0*exp(-_Qtssp/(R*_temperatureVector[position]));
      _tssdVector[position] = _tssd0*exp(-_Qtssd/(R*_temperatureVector[position]));
      break ;
  }

  //Kinetics
  switch(_kinetics)
  {
    case 1: // Nothing to do
      break ;

    case 2: // Constant kinetics
      _KdVector[position] = _Kd0  ;
      _KnVector[position] = _Kn0  ;
      _KgVector[position] = _Kmob0;
      break ;

    case 3: // Exponential dependence

      // Hydrides formation energy : cubic fit
      double Eth = -_Eth0 + _Eth1*_temperatureVector[position]
                   -_Eth2*pow(_temperatureVector[position],2)
                   +_Eth3*pow(_temperatureVector[position],3);

      // Nucleation kinetics
      if (_optionNGD>0) // && (_solutionContentVector[position] >= _tsspVector[position]))
      {
        double sigma = 1e-4 ; // .23 * _tsspVector[position] ;
        double delta = 0 ; //.30 * _tsspVector[position] ;
        _probaVector[position] = 0.5 * erfc((delta + _tsspVector[position] - _solutionContentVector[position]) / (1.4142 * sigma)) ;
        
        
        double Kn0 = _Kn0 * factor_f_alpha(position)  * _probaVector[position] ;
        _KnVector[position] = ((17000-_hydridesContentVector[position])/17000) * Kn0 * exp(-Eth/(kb*_temperatureVector[position])) ;

      }
      
      else
        _KnVector[position] = 0. ;


      // Dissolution kinetics
      if ((_optionNGD>1) && (_solutionContentVector[position] < _tssdVector[position] && _hydridesContentVector[position]>0.))
      {
        double modif = 1 ; //pow(_hydridesContentVector[position]/_totalContentVector[position], 2) ;
        
        _KdVector[position] = _Kd0 * exp(-_Ediss/(kb*_temperatureVector[position])) * modif ;
      }
      else
        _KdVector[position] = 0. ;

      // Growth kinetics
      if ((_optionNGD>2) && (_solutionContentVector[position] > _tssdVector[position] && _hydridesContentVector[position]>0.))
      {
        _KmobVector[position] = _Kmob0*leverRule(position)*factor_f_alpha(position)*exp(-_Eg/(kb*_temperatureVector[position]));
        _KthVector[position]  =  _Kth0*leverRule(position)*factor_f_alpha(position)*exp(-Eth/(kb*_temperatureVector[position]));
        _KgVector[position] = 1./(1./_KmobVector[position] + 1./_KthVector[position]) ;
      }
      else
        _KgVector[position] = 0. ;

    break;
  }

  // Diffusion
  if(_fickLaw || _soretEffect)
  {
    switch(std::max(_fickLaw,_soretEffect))
    {
      case 0:
        _coeffFickVector[position] = 0 ;
        break;

      case 1:
        _coeffFickVector[position] = _D0 ;
        break;

      case 2:
        _coeffFickVector[position] = _D0*exp(-_Ediff/(kb*_temperatureVector[position])) ;
        break;
    }
  }
}

double HydrogenBehaviorModel :: leverRule(int position)
{
  // Total atomic fraction of hydrogen
  double xtot   = _totalContentVector[position]   / (Mh*(  _totalContentVector[position]/Mh  +  (1e6 - _totalContentVector[position])/Mzr)) ;
  // TSSd atomic fraction
  double xalpha =     _tssdVector[position]       / (Mh*(    _tssdVector[position]/Mh        +  (1e6 -     _tssdVector[position]    )/Mzr)) ;
  // alpha / alpha+delta boundaries
  double xdelta = -9.93e-11*pow(_temperatureVector[position],3) + // 3rd degree polynome fitting
                  8.48e-8*pow(_temperatureVector[position],2) -
                  5.73e-5*_temperatureVector[position] + 0.623 ;

  _lever[position] = ((xtot - xalpha)/(xdelta - xalpha));
  // return the value of the lever rule
  return ((xtot - xalpha)/(xdelta - xalpha));
}

double HydrogenBehaviorModel :: factor_f_alpha(int position)
{
  double xhyd   = _hydridesContentVector[position]/ (Mh*(_hydridesContentVector[position]/Mh +(1e6 - _hydridesContentVector[position])/Mzr)) ;
  // alpha / alpha+delta boundaries
  double xdelta = -9.93e-11*pow(_temperatureVector[position],3) + // 3rd degree polynome fitting
                  8.48e-8*pow(_temperatureVector[position],2) -
                  5.73e-5*_temperatureVector[position] + 0.623 ;

 double xalpha =     _tssdVector[position]       / (Mh*(    _tssdVector[position]/Mh        +  (1e6 -     _tssdVector[position]    )/Mzr)) ;

  _f_alpha[position] = 1. - xhyd/(xdelta - xalpha) ;
  
  // return the value of the lever rule
  return 1. - xhyd/xdelta ;
}


//------------------------------------------ Auxiliary functions ------------------------------------------
    // Domain definition
void HydrogenBehaviorModel :: computeLocations(double x0, double xEnd, double bias, std::vector<double>& x)
{
  // number of cells for this domain
  const unsigned long n = x.size()  ;

  double  sum = 1. + bias ;
  for(int k=0; k<n-3; k++)
    sum = 1. + bias*sum ;

  const double initialLenght = (xEnd - x0)/sum ;

  x[0] = x0 ;
  x[1] = x0 + initialLenght ;
  for(int k=2; k<n-1; k++)
    x[k] = x[k-1] + bias*(x[k-1] - x[k-2]) ;
  x[n-1] = xEnd  ;
}

    // Interpolation
void HydrogenBehaviorModel :: spatialeInterpolation(std::vector<double>& refX, std::vector<double>& refY, std::vector<double>& vectorX, std::vector<double>& vectorY)
{
  int k = 1 ;
  for(int i=0; i<vectorX.size(); i++) {
    if(vectorX[i] > refX[k])
      k++ ;
    vectorY[i] = refY[k-1] + (refY[k] - refY[k-1]) * (vectorX[i] - refX[k-1]) / (refX[k] - refX[k-1]);
  }
}

    // Time step
double HydrogenBehaviorModel :: computeTimeStep()
{
  double dt = 10. ;
  double delta_t_fick  = dt ;
  double delta_t_soret = dt ;

  // We use the diffusion convergence criteria dt < dx^2/(2D) to find the minimum value of dt required to compute diffusion
  for(int k=1; k<_nbNodes; k++)
  {
    if(_fickLaw)
      delta_t_fick  = std::min(delta_t_fick,  pow(_positionsVector[k]-_positionsVector[k-1],2)/(2*_coeffFickVector[k]));

//    if(_soretEffect)
//    {
//      double dT_dx(0.) ;
//      if (!(k==_nbNodes-1))
//      {
//        dT_dx = (_temperatureVector[k+1]-_temperatureVector[k])/(_positionsVector[k+1] - _positionsVector[k]);
//
//        delta_t_soret = std::min(delta_t_soret, pow(_positionsVector[k]-_positionsVector[k-1],2)/(2*(_coeffFickVector[k]*_Qstar*dT_dx*_solutionContentVector[k]/(R*pow(_temperatureVector[k],2)))));
//      }
//    }
  }

  if(_fickLaw)
    dt = std::min(dt, delta_t_fick);

  if(_soretEffect)
    dt = std::min(dt, delta_t_soret);

  // NGD criteria
  _nuclIsHappening = false ;
  _dissIsHappening = false ;
  _growIsHappening = false ;
  if(_optionNGD>0)
  {
    for(int k=0; k<_nbNodes; k++)
    {
      if(_kinetics>1)
      {
          // Nucleation over TSSp
          if ((_optionNGD>0)) // && (_solutionContentVector[k] >= _tsspVector[k]))
            _nuclIsHappening = true ;

          // Dissolution if hydrides under TSSd
          if ((_optionNGD>1) && (_solutionContentVector[k] < _tssdVector[k] && _hydridesContentVector[k]>0.))
            _dissIsHappening = true ;

          // Growth if hydrides above TSSd
          if ((_optionNGD>2) && (_solutionContentVector[k] > _tssdVector[k] && _hydridesContentVector[k]>0.))
            _growIsHappening = true ;
      }
    }
  }
  double KnMax = _KnVector[0];
  double KgMax = _KgVector[0];
  double KdMax = _KdVector[0];
  for(int k=1; k<_nbNodes; k++)
  {
    KnMax = std::max(KnMax, _KnVector[k]);
    KgMax = std::max(KgMax, _KgVector[k]);
    KdMax = std::max(KdMax, _KdVector[k]);
  }
    // Nucleation
  if(_nuclIsHappening)
    dt = std::min(dt, 0.3/KnMax); //0.02
    // Growth
  if(_growIsHappening && KgMax>0)
    dt = std::min(dt, 0.3/KgMax); //2.5e-3
    // Dissolution
  if(_dissIsHappening)
    dt = std::min(dt, 0.3/KdMax); //0.03

  dt = std::max(dt,1e-6) ;

  // Security factor to be sure to converge
  double secuFactor = 1. ;
  return dt/secuFactor;

}


// ------------------------------------------   Getters   ------------------------------------------
    // Spacial profile
std::vector<double>& HydrogenBehaviorModel :: returnKdVector()              {return _KdVector ;}
std::vector<double>& HydrogenBehaviorModel :: returnKgVector()              {return _KgVector ;}
std::vector<double>& HydrogenBehaviorModel :: returnKmobVector()            {return _KmobVector ;}
std::vector<double>& HydrogenBehaviorModel :: returnKthVector()             {return _KthVector ;}
std::vector<double>& HydrogenBehaviorModel :: returnfVector()               {return _f_alpha ;}
std::vector<double>& HydrogenBehaviorModel :: returnLeverVector()           {return _lever ;}
std::vector<double>& HydrogenBehaviorModel :: returnKnVector()              {return _KnVector ;}
std::vector<double>& HydrogenBehaviorModel :: returnTSSpVector()            {return _tsspVector;}
std::vector<double>& HydrogenBehaviorModel :: returnTSSdVector()            {return _tssdVector;}
std::vector<double>& HydrogenBehaviorModel :: returnTotalContentVector()    {return _totalContentVector ;}
std::vector<double>& HydrogenBehaviorModel :: returnSolutionContentVector() {return _solutionContentVector ;}
std::vector<double>& HydrogenBehaviorModel :: returnHydridesContentVector() {return _hydridesContentVector;}
std::vector<double>& HydrogenBehaviorModel :: returnPositionVector()        {return _positionsVector;}
std::vector<double>& HydrogenBehaviorModel :: returnTemperatureVector()     {return _temperatureVector;}
std::vector<double>& HydrogenBehaviorModel :: returnCoeffFickVector()       {return _coeffFickVector;}
std::vector<double>& HydrogenBehaviorModel :: returnFlux()                  {return _flux ;}
std::vector<double>& HydrogenBehaviorModel :: returnProba()                 {return _probaVector ;}

double HydrogenBehaviorModel :: returnTimeStep() {return _dt ;}

double HydrogenBehaviorModel :: returnTotalContent()
{
  double sum(0.);
  for(int k=0; k<_nbNodes-2; k++)
    sum += (_totalContentVector[k] + 4*_totalContentVector[k+1] + _totalContentVector[k+2])*0.5*(_positionsVector[k+2]-_positionsVector[k])/3;
  return sum ;
}
double HydrogenBehaviorModel :: returnSolutionContent()
{
  double sum(0.);
  for(int k=0; k<_nbNodes-1; k++)
    sum += 0.5*(_solutionContentVector[k]+_solutionContentVector[k+1])*(_positionsVector[k+1]-_positionsVector[k]);
  return sum/_sampleLenght ;
}
double HydrogenBehaviorModel :: returnHydridesContent()
{
  double sum(0.);
  for(int k=0; k<_nbNodes-1; k++)
    sum += 0.5*(_hydridesContentVector[k]+_hydridesContentVector[k+1])*(_positionsVector[k+1]-_positionsVector[k]);
  return sum/_sampleLenght ;
}
double HydrogenBehaviorModel :: returnAdvancement()
{
  return (_totalContentVector[0]-_solutionContentVector[0])/(_totalContentVector[0]-_tssdVector[0]);
}

void HydrogenBehaviorModel::print(){ //Auxiliary method to find what's wrong

}
