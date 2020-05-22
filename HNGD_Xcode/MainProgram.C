#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>

#include "InOut.hpp"
#include "HydrogenBehaviorModel.h"

using namespace std;

int main(int argc, char* argv[])
{
  double  t       = 0.,             // Time
          dt      = 0.,             // Time step
          dtPrint = 0.,             // Duration between 2 prints
          t_end   = 0.;             // Totale duration

  ofstream error,                   // error log file
           output;                  // output csv file

  // interpolation function for temperature history
  vector<double> interpolate(double t, vector<double> time_stamps, vector<vector<double>> temp_stamps);
  
  // HNGD object
  HydrogenBehaviorModel hydrogen_behavior;

  //----------------------- Define execution folder and file names --------------------
  string path_exec = "/Users/fpp8/OneDrive - The Pennsylvania State University/Hydride_Modeling/Further HNGD/HNGD_Xcode/HNGD_Xcode/" ;
  
  string input_folder = "input_files/" ;
  
  // Default file names
  string settings_name = "1_settings.txt" ;
  string treatment_name= "2_temperature.txt" ;
  string physics_name  = "3_physics.txt" ;
  string hydroIC_name  = "4_hydrogen.txt" ;
  string output_name   = "output" ;
  
  // Names can be given as parameters
  if(argc > 1)
  {
    string name(argv[1]);
    
    settings_name  = input_folder + name + "_set.txt" ;
    treatment_name = input_folder + name + "_temp.txt";
    physics_name   = input_folder + name + "_phys.txt";
    hydroIC_name   = input_folder + name + "_hyd.txt" ;
    output_name    = name + "_out.csv" ;
    
  }
  
  error.open   (path_exec + output_name + "_error.log", ios::out);
  output.open  (path_exec + output_name, ios::out);
  
  if (error.fail() || output.fail())
  {
    cout  << "File opening error!\nProgram stopped.\n";
    error << "File opening error!\nProgram stopped.\n";
    exit(1);
  }

  //-------------------- Input reading ------------------------

  // Simulation settings
  short int nbSettings = 12 ;
  double settings[nbSettings];
  InOut::getSettings(nbSettings, settings, path_exec, settings_name);

  // Physical parameters
  int nbPhysicalParameters = 18 ;
  double physicalParameters[nbPhysicalParameters];
  InOut::getPhysics(nbPhysicalParameters, physicalParameters, path_exec, physics_name);

  // Thermal treatment
  vector<double>         time_temp(0), // input time stamps for temperature (s)
                         pos_temp(0);  // input positions for temperature (cm)
  vector<vector<double>> temp_inp;     // input temperature values (K)
  
  vector<vector<double>> thermal_treatment = InOut::getThermalTreatment(path_exec, treatment_name);
  pos_temp = thermal_treatment[0] ;
  time_temp= thermal_treatment[1] ;
  for(int k=0; k<thermal_treatment.size()-2; k++)
    temp_inp.push_back(thermal_treatment[k+2]) ;
  
  t_end = time_temp[time_temp.size()-1];


  // Initial H profile
  vector<vector<double>> hydrogenIC = InOut::getICHydrogen(path_exec, hydroIC_name);
  vector<double> pos_hyd = hydrogenIC[0] ;
  vector<double> hyd_inp = hydrogenIC[1] ;
  

  //---------------- System Initialisation ---------------------

  // Settings
  short int typeSimu ;
  int nbNodes, nbPosPrint ;
  typeSimu    = (int)settings[0] ;
  nbNodes     = settings[2];
  dt          = settings[9];
  dtPrint     = settings[10];
  nbPosPrint  = min(settings[2],settings[11]);

  // Physics
  vector<double> temp = interpolate(t, time_temp, temp_inp);
  hydrogen_behavior.getInitialConditions(settings, physicalParameters, pos_hyd, hyd_inp, pos_temp, temp);

  // Output file
  const short int nbOutput = 5 ; /* HERE */
  int listPosPrint[nbPosPrint] ;
  if(typeSimu==1) // For a distribution simulation
  {
      InOut::writeInitialOutput(hydrogen_behavior, path_exec, output_name, nbNodes, nbOutput, nbPosPrint, listPosPrint);
      InOut::writeOuput(hydrogen_behavior, path_exec, output_name, nbNodes, nbOutput, t, 0., nbPosPrint, listPosPrint);
  }

  //-------------------- Time loop --------------------------
  double printCountdown(0.);
  if(typeSimu==1)  // For a distribution simulation
  {
    do
    {
    // Interpolation of input data using the function "interpolate" implemented below
      t += hydrogen_behavior.returnTimeStep() ;
      printCountdown += hydrogen_behavior.returnTimeStep() ;
      temp = interpolate(t, time_temp, temp_inp);

    // Computation
      hydrogen_behavior.getInput(t, pos_temp, temp, dt);
      hydrogen_behavior.computeProperties();

    // Write output
      if (printCountdown >= dtPrint)
      {
        InOut::writeOuput(hydrogen_behavior, path_exec, output_name, nbNodes, nbOutput, t, 0., nbPosPrint, listPosPrint);
        printCountdown = 0. ;
      }

    } while ( t < t_end );

    if(printCountdown>0.) // To be sure that the final state is printed
      InOut::writeOuput(hydrogen_behavior, path_exec, output_name, nbNodes, nbOutput, t, 0., nbPosPrint, listPosPrint);
  }

  // ------------------- End of computation -----------------------

  error << "The calculation was performed!\n";
  cout  << "The calculation was performed!\n";

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


// else  // for a TTT diagram
// {
//   // Output file
//   ofstream output ;
//   output.open("output.csv", std::ios_base::app);
//
//   // Completion parameter
//   double completionTarget = 0.99 ;
//   double completion = 0.;
//   for(int k=0; k<i-1; k++)
//   {
//     for(int j=0; j<time_temp[k+1]; j++)
//     {
//       temp = temp_inp[k] + (j/time_temp[k+1])*(temp_inp[k+1]-temp_inp[k]) ;
//       cout<<temp<<endl;
//       t = 0. ;
//       completion = 0.;
//       hydrogen_behavior = HydrogenBehaviorModel() ;
//       hydrogen_behavior.getInitialConditions(settings, physicalParameters, temp, 0.);
//       if(hydrogen_behavior.returnTSSdVector()[0] > hydrogen_behavior.returnTotalContentVector()[0])
//       {
//         cout<<"Error: at temperature "<<temp<<"K the solubility is higher than the hydrogen content\n";
//         continue ;
//       }
//
//       // Each simulation runs untill the completion reaches the target
//       while (completion < completionTarget)
//       {
//         hydrogen_behavior.getInput(t, temp, grad, dt);
//         hydrogen_behavior.computeProperties();
//         t += hydrogen_behavior.returnTimeStep() ;
//
//         printCountdown += hydrogen_behavior.returnTimeStep() ;
//         if (printCountdown >= dtPrint)
//         {
//           cout<<"t="<<t<<"s\n";
//           printCountdown = 0. ;
//         }
//
//         completion = hydrogen_behavior.returnAdvancement();
//       }
//       output << temp << ',' << t << ",\n" ;
//     }
//   }
//   output.close();
// }
