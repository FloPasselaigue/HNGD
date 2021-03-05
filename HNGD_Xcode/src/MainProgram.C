/**
  This code is the implementation of the Hydride Nucleation-Growth-Dissolution (HNGD) model.
  It was developped by Florian Passelaigue, based on Evrard Lacroix's PhD work.
 
  The model development is described in
  E. Lacroix, P.-C. A. Simon, A. T. Motta, and J. Almer, “Zirconium hydride precipitation and
  dissolution kinetics in the hysteresis region in zirconium alloys” ASTM (submitted), 2019.
 
  The original implementation. verification and validation of this code is described in chapter 2 of
  F. Passelaigue, "Hydride Nucleation-Growth-Dissolution model: implementation in BISON",
  Master of Science thesis, The Pennsylvania State University, 2020
  available on Penn State library:
  https://etda.libraries.psu.edu/catalog/17572fpp8
 
  NB: the file structure was changed and new classes were created since the MS thesis was written.
    Refer to the README file for more details
 */

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>

#include "InOut.hpp"
#include "HNGD.hpp"
#include "EvalEvolution.hpp"

using namespace std;

int main(int argc, char* argv[])
{
  ofstream output;                  // output csv file

  // Interpolation function for temperature history
  // This function returns the temperature profile as a vector, based on the current time t
  // and the temperature steps specified in the temperature input file.
  vector<double> interpolate(double t, vector<double> time_stamps, vector<vector<double>> temp_stamps);
  
  
  // The following function and variables are used to
  // determine if the current state should be written in the output
  
  // This function checks if a time stamp specified in the temperature input file was passed.
  bool changeInterval(double t, double dt, vector<double> time_stamps);
  
  // The EvalEvolution objects determine if a given quantity
  // changed sufficiently to justify a print in the output
  EvalEvolution evalEvolTemp ; // for temperature
  EvalEvolution evalEvolHyd ;  // for hydrogen

  
  //----------------------- DEFINE EXECUTION FOLDER AND INPUT FILES NAMES --------------------
  
  // Path to the folder to use /*custom*/
  string path_exec = "/Users/fpp8/OneDrive - The Pennsylvania State University/Hydride_Modeling/Further HNGD/HNGD_Xcode/HNGD_Xcode/" ;
  
  // Default file names.
  // These files will be used if no argument is given to the program when launched.
  // These files must be in the folder defined by the path_exec variable
  string settings_name = "1_settings.txt" ;
  string treatment_name= "2_temperature.txt" ;
  string physics_name  = "3_physics.txt" ;
  string hydroIC_name  = "4_hydrogen.txt" ;
  string output_name   = "output" ;
  
  // Name of the folder containing the input files.
  // The input files specified by the argument given at
  // launch must be placed in this folder
  string input_folder = "input_files/" ;
  
  // Argument reading
  // The argument must be used as a prefix to all input files.
  if(argc > 1)
  {
    string name(argv[1]);
    
    settings_name  = input_folder + name + "_set.txt" ;
    treatment_name = input_folder + name + "_temp.txt";
    physics_name   = input_folder + name + "_phys.txt";
    hydroIC_name   = input_folder + name + "_hyd.txt" ;
    output_name    = name + "_out.csv" ;
  }
  
  output.open(path_exec + output_name, ios::out);
  if (output.fail())
  {
    cout  << "File opening error!\nProgram stopped.\n";
    exit(1);
  }

  
  
  //-------------------- INPUT READING ------------------------
  // More details on the format and information contained in each file in the README file

  // Simulation settings contained in the *_set.txt file
  short int nbSettings = 7 ;
  double settings[nbSettings];
  InOut::getSettings(nbSettings, settings, path_exec, settings_name);

  // Physical parameters contained in the *_phys.txt file
  int nbPhysicalParameters = 16 ;
  double physicalParameters[nbPhysicalParameters];
  InOut::getPhysics(nbPhysicalParameters, physicalParameters, path_exec, physics_name);

  // Thermal treatment contained in the *_temp.txt file
  vector<double>         time_temp(0), // input time stamps for temperature (s)
                         pos_temp(0);  // input positions for temperature (cm)
  vector<vector<double>> temp_inp;     // input temperature values (K)
  
  vector<vector<double>> thermal_treatment = InOut::getThermalTreatment(path_exec, treatment_name);
  pos_temp = thermal_treatment[0] ;
  time_temp= thermal_treatment[1] ;
  for(int k=0; k<thermal_treatment.size()-2; k++)
    temp_inp.push_back(thermal_treatment[k+2]) ;
  
  double t_end = time_temp[time_temp.size()-1];

  // Initial H profile contained in the *_hyd.txt file
  vector<vector<double>> hydrogenIC = InOut::getICHydrogen(path_exec, hydroIC_name);
  vector<double> pos_hyd = hydrogenIC[0] ;
  vector<double> hyd_inp = hydrogenIC[1] ;
  

  //---------------- SYSTEM INITIALIZATION ---------------------

  // The HNGD object collects the input information to build
  // a Sample and the objects associated with each phenomenon
//  HNGD hngd(settings, physicalParameters) ;
  unsigned int n = pos_temp.size() - 1;
  double xEnd = pos_temp[n];
  HNGD hngd(settings, physicalParameters, xEnd, settings[6]) ; // TODO: less parameters
  
  // Some of the settings are needed for the time loop
  int nbNodes        = settings[0];
  int nbPosPrint     = settings[0];
  double dtPrint     = settings[4];
  double critPrint   = settings[5] ;

  // Initialize the temperature and hydrogen profiles
  double t = 0. ;
  vector<double> temp = interpolate(t, time_temp, temp_inp);
  hngd.getInitialConditions(pos_hyd, hyd_inp, pos_temp, temp);
  
  // Associate the EvalEvolution objects to the profiles
  evalEvolHyd.setProfile(hngd.returnSample()->returnTotalContent()) ;
  evalEvolHyd.setCriterion(critPrint) ;
  
  evalEvolTemp.setProfile(hngd.returnSample()->returnTemperature()) ;
  evalEvolTemp.setCriterion(critPrint) ;

  // Initialize the output file
  const short int nbOutput = 5 ; /*custom*/
  int listPosPrint[nbPosPrint] ;
  InOut::writeInitialOutput(hngd, path_exec, output_name, nbNodes, nbOutput, nbPosPrint, listPosPrint, settings[6]); //TODO: less parameters
  InOut::writeOuput(hngd, path_exec, output_name, nbNodes, nbOutput, t, 0., nbPosPrint, listPosPrint);
  

  
  //-------------------- TIME LOOP --------------------------
  double printCountdown(0.);
  
  do
  {
  // Interpolation of input temperature using the function "interpolate" implemented below
    t += hngd.returnTimeStep() ;
    printCountdown += hngd.returnTimeStep() ;
    temp = interpolate(t, time_temp, temp_inp);

  // Compute the new system state
    hngd.getInput(pos_temp, temp);
    hngd.compute();

  // Write output if enough time has elapsed or if the temperature profile has changed
  // or if the hydrogen profile has changed of if a time stamp was reached
    if (printCountdown >= dtPrint ||
       evalEvolTemp.evaluate(hngd.returnSample()->returnTemperature()) ||
       evalEvolHyd.evaluate(hngd.returnSample()->returnTotalContent()) ||
        changeInterval(t, hngd.returnTimeStep(), time_temp))
    {
      InOut::writeOuput(hngd, path_exec, output_name, nbNodes, nbOutput, t, 0., nbPosPrint, listPosPrint);
      printCountdown = 0. ;
    }

  } while ( t < t_end );
  

  // ------------------- END OF COMPUTATION -----------------------

  cout  << "The calculation was performed!\n";

  // Sound notification at the end of simulation /*custom*/
  system("open \"/Users/fpp8/OneDrive - The Pennsylvania State University/Hydride_Modeling/Further HNGD/HNGD_Xcode/HNGD_Xcode/zelda.mp3\" -a VLC");
  
  return 0;
}


vector<double> interpolate(double t, vector<double> time_stamps, vector<vector<double>> temp_stamps)
{
  
  // If end of simulation, no interpolation
  if(t >= time_stamps[time_stamps.size()-1])
    return temp_stamps[temp_stamps.size()-1] ;
  
  else
  {
    int k=0 ;
    while(t >= time_stamps[k])
      k ++ ;

    vector<double> interpolated_temp_stamps(temp_stamps[0].size()) ;
    
    for(int i=0; i<temp_stamps[0].size(); i++)
      interpolated_temp_stamps[i] = temp_stamps[k-1][i] + (temp_stamps[k][i] - temp_stamps[k-1][i]) * (t - time_stamps[k-1]) / (time_stamps[k] - time_stamps[k-1]) ;

    return interpolated_temp_stamps ;
  }
}

bool changeInterval(double t, double dt, vector<double> time_stamps)
{
  // If end of simulation, no interpolation
  if(t >= time_stamps[time_stamps.size()-1] || t+dt >= time_stamps[time_stamps.size()-1])
    return true ;
  
  int k=0 ;
  while(t >= time_stamps[k])
    k ++ ;
  
  return t+dt > time_stamps[k] ;
}
