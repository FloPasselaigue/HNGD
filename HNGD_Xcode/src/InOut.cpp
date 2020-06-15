
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>

#include "InOut.hpp"

using namespace std;

// -------------------------------- Input reading --------------------------------
void InOut::getSettings(int nb, double* settings, string path_exec, string file_name)
{
    ifstream inpset ;
    inpset.open(path_exec + file_name);
    
    if(inpset.fail())
    {
        cout << "Setting file opening failed" << endl ;
        exit(1) ;
    }
    
    for(int k=0; k<nb; k++)
    {
      inpset >> settings[k] ;
      inpset.ignore(999,'\n');
      inpset.ignore(999,'\n');
    }
    
    InOut::writeSettingsInCheck(settings);
    
    inpset.close();
}

void InOut::getPhysics(int nb, double* physics, string path_exec, string file_name)
{
    ifstream inpphys ;
    inpphys.open(path_exec + file_name);
    
    if(inpphys.fail())
    {
        cout << "Physics file opening failed" << endl ;
        exit(1) ;
    }
    
    for(int k=0; k<nb; k++)
    {
      inpphys >> physics[k] ;
      inpphys.ignore(999,'\n');
    }
    
    InOut::writePhysicsInCheck(physics);
    
    inpphys.close();
}

vector<vector<double>> InOut::getThermalTreatment(string path_exec, string file_name)
{
    ifstream inptemp ;
    inptemp.open (path_exec + file_name);
    
    if(inptemp.fail())
    {
        cout << "Temperature file opening failed" << endl ;
        exit(1) ;
    }
    
    ofstream chkinp ;
    chkinp.open("input_check.txt", std::ios_base::app);
    chkinp << endl << endl << "time (s), temperature (K), gradient (K/cm)";
    
    vector<double> pos_temp(0) ;
    vector<double> time_temp(0) ;
    vector<vector<double>> temp_inp(0) ;
    
    // Positions
    string line_pos ;
    getline(inptemp, line_pos);
    istringstream stream_pos(line_pos) ;
    double pos ;
    while(stream_pos >> pos)
    pos_temp.push_back(pos) ;

    unsigned long nb_pos = pos_temp.size() ;

    unsigned short int i=0;
    do
    {
        double time,temperature ;
        inptemp >> time ;
        time_temp.push_back(time);

        if(i>0 && (time_temp[i] < time_temp[i-1])) // check that the input is consistent
        {
          cout  << "Input time not monotonically increasing!\nPlease check input.\n";
          exit(1);
        }

        vector<double> new_vec(0) ;
        for(int k=0; k<nb_pos; k++)
        {
          inptemp >> temperature ;
          new_vec.push_back(temperature);
        }
        temp_inp.push_back(new_vec) ;

        i++ ;

    } while(!(inptemp.eof()));

    temp_inp.pop_back() ;
    time_temp.pop_back();

//    for(int k=0; k<temp_inp.size(); k++)
//    {
//    chkinp << time_temp[k] ;
//    for(int j=0; j<temp_inp[0].size(); j++)
//      cout<< "\t" << temp_inp[k][j] ;
//
//    chkinp << endl ;
//    }

    cout << "Number of specified temperature steps = " << i-1 << endl;

    double t_end = time_temp[i-1];
    chkinp << endl << endl << "end time (s) = " << t_end;
    
    // Return the thermal treatment as a vector of vector<double>
    // First vector contains the positions
    // Second vector contains the time stamps
    // The others contain the temperature values
    vector<vector<double>> thermal_treatment(0) ;
    thermal_treatment.push_back(pos_temp) ;
    thermal_treatment.push_back(time_temp);
    for(int k=0; k<temp_inp.size(); k++)
        thermal_treatment.push_back(temp_inp[k]);
    
    return thermal_treatment ;
}

vector<vector<double>> InOut::getICHydrogen(string path_exec, string file_name)
{
    ifstream inphyd;
    inphyd.open(path_exec + file_name);
    if(inphyd.fail())
    {
        cout << "Hydrogen IC file opening failed" << endl ;
        exit(1) ;
    }
    
    vector<double> pos_hyd(0);
    vector<double> hyd_inp(0);
    
    string line_pos ;
    string line_hyd ;
    
    double pos ;
    double hyd ;
    
    getline(inphyd, line_pos);
    getline(inphyd, line_hyd);
    
    istringstream stream_pos_hyd(line_pos) ;
    istringstream stream_hyd(line_hyd) ;
    
    while(stream_pos_hyd >> pos){
      pos_hyd.push_back(pos) ;
      stream_hyd >> hyd ;
      hyd_inp.push_back(hyd) ;
    }
    
    vector<vector<double>> hydrogenIC(0) ;
    hydrogenIC.push_back(pos_hyd);
    hydrogenIC.push_back(hyd_inp);
    
    return hydrogenIC ;
}



// -------------------------------- Check file writing --------------------------------
void InOut::writeSettingsInCheck(double * settings)
{
  ofstream chkinp ;
  chkinp.open("input_check.txt", ios::out);

  //typeSimu
  switch((int)settings[0])
  {
    case 1: chkinp << "Computation of hydrogen distribution\n"; break;
    case 2: chkinp << "Computation of a TTT diagram\n" ; break ;
    default:cout   << "Invalid value for typeSimu in settings\nProgram has stopped\n";exit(1);break;
  }

  //optionNGD
  switch ((int)settings[1])
  {
    case 0: chkinp << "This simulation doesnt compute precipitation\n"; break;
    case 1: chkinp << "This simulation compute Nucleation\n" ; break ;
    case 2: chkinp << "This simulation compute Nucleation and Dissolution\n" ; break ;
    case 3: chkinp << "This simulation compute Nucleation Growth and Dissolution\n" ; break ;
    default:cout << "Invalid value for optionNGD in settings\nProgram has stopped\n";exit(1);break;
  }

  // nbNodes
  if((int)settings[2]>0)
    chkinp << "Space discretization: " << settings[2] << " cells\n";
  else
    chkinp << "No space discretization\n" ;

  // kinetics
  if((int)settings[1]>0){
    switch((int)settings[3])
    {
      case 1: chkinp << "1 -> This simulation compute the equilibrium\n" ; break;
      case 2: chkinp << "2 -> Kinetics are constant\n" ; break;
      case 3: chkinp << "3 -> Kinetics depend exponentially on the temperature\n" ; break;
      default:cout << "Invalid value for kinetics in settings\nProgram has stopped\n";exit(1);break;
    }
  }

  // solubility
  switch((int)settings[4])
  {
    case 0: chkinp << "0 -> Solubilities depend linearly on the temperature\n" ; break;
    case 1: chkinp << "1 -> Solubilities depend exponentially on the temperature\n" ; break;
    default:cout << "Invalid value for solubilities in settings\nProgram has stopped\n";exit(1);break;
  }

  // diffusion
  switch((int)settings[5])
  {
    case 0: chkinp << "0 -> No diffusion computation\n" ; break;
    case 1: chkinp << "1 -> Diffusion is computed with constant coefficient D\n" ; break;
    case 2: chkinp << "2 -> Diffusion is computed with exponential coefficient D\n" ; break;
    default:cout << "Invalid value for diffusion in settings\nProgram has stopped\n";exit(1);break;
  }

  // Soret effect
  switch((int)settings[6])
  {
    case 0: chkinp << "0 -> No Soret effect computation\n" ; break;
    case 1: chkinp << "1 -> Soret efffect is computed with constant coefficient D\n" ; break;
    case 2: chkinp << "2 -> Soret efffect is computed with exponential coefficient D\n" ; break;
    default:cout << "Invalid value for Soret effect in settings\nProgram has stopped\n";exit(1);break;
  }

  chkinp << "(If the two previous effet don't have the same option, D is computed according to the highest option)\n";

  // bias
  chkinp << "Bias: " << settings[7] << endl ;

  // Lenght
  chkinp << "Sample lenght: " << settings[8] << endl;

  // dt
  if(settings[9]<=0)
    chkinp << "Adaptative time step" << endl;
  else
    chkinp << "Time step: " << settings[9] << endl;

  // dtPrint
  chkinp << "Printing time step: " << settings[10] << endl << endl ;

  // nbPositionPrint
  chkinp << "The output is written for " << min(settings[2],settings[11]) << " locations" << endl << endl ;

  chkinp.close();

  if((int)settings[2]==0 && (int)settings[5]>0 && (int)settings[6]>0)
  {
    cout << "Error, you can't have diffusion without space discretization\n Program has stopped\n" ;
    exit(1);
  }
}

void InOut::writePhysicsInCheck(double * physicalParameters)
{
  ofstream chkinp ;
  chkinp.open("input_check.txt", std::ios_base::app);

  chkinp << "Kd0:\t"  << physicalParameters[0]  << "\ts-1\n" ;
  chkinp << "Ediss:\t"<< physicalParameters[1]  << "\teV/at\n" ;
  chkinp << "Kn0:\t"  << physicalParameters[2]  << "\ts-1\n" ;
  chkinp << "Eth0:\t" << physicalParameters[3]  << "\teV/mol\n" ;
  chkinp << "Eth1:\t" << physicalParameters[4]  << "\teV/mol/K\n" ;
  chkinp << "Eth2:\t" << physicalParameters[5]  << "\teV/mol/K2\n" ;
  chkinp << "Eth3:\t" << physicalParameters[6]  << "\teV/mol/K3\n" ;
  chkinp << "Kmob0:\t"<< physicalParameters[7]  << "\ts-1\n" ;
  chkinp << "Kth0:\t" << physicalParameters[8] << "\ts-1\n" ;
  chkinp << "Eg:\t"   << physicalParameters[9] << "\teV/at\n" ;
  chkinp << "TSSp0:\t"<< physicalParameters[10] << "\twt.ppm\n" ;
  chkinp << "Qtssp:\t"<< physicalParameters[11] << "\tJ/(mol.K)\n" ;
  chkinp << "TSSd0:\t"<< physicalParameters[12] << "\twt.ppm\n" ;
  chkinp << "Qtssd:\t"<< physicalParameters[13] << "\tJ/(mol.K)\n" ;
  chkinp << "D0:\t"   << physicalParameters[14] << "\tJ/(mol.K)\n" ;
  chkinp << "Ediff:\t"<< physicalParameters[15] << "\tJ/(mol.K)\n" ;
  chkinp << "Q*:\t"   << physicalParameters[16] << "\tJ/mol\n" ;
  chkinp << "tau:\t"  << physicalParameters[17] << "\ts\n" ;

  chkinp.close();
}


// -------------------------------- Output file writing --------------------------------

void InOut::writeOuput(HNGD hngd, string path_exec, string output_name, int nbNodes, int nbOutput, double t, double temp, int nbPosPrint, int* listPosPrint)
{
    ofstream output ;
    output.open(path_exec + output_name, std::ios_base::app);

    /* HERE */
    std::vector<double> listVector[nbOutput];
    listVector[0] = hngd.returnSample()->returnTotalContent();
    listVector[1] = hngd.returnSample()->returnSolutionContent();
    listVector[2] = hngd.returnSample()->returnHydrideContent();
    listVector[3] = hngd.returnDiss()->returnRate();
    listVector[4] = hngd.returnNuc()->returnRate();
    listVector[5] = hngd.returnGro()->returnRate();


    output << hngd.returnTimeStep() << "," << t << ","  ;
    for(int i=0; i<nbOutput; i++){
        for(int j=0; j<nbPosPrint; j++)
        {
            output << (listVector[i])[listPosPrint[j]] << ",";
        }
        output << ',' ;
    }
    output << "\n"   ;


    cout << "t= " << t << "s"<< endl;

    output.close();
}

void InOut :: writeInitialOutput(HNGD hngd, string path_exec, string output_name, int nbNodes, int nbOutput, int nbPosPrint, int* listPosPrint)
{
  ofstream output ;
  output.open(path_exec + output_name, std::ios_base::app);

  if(nbNodes>0)
  {
        output << " ,Positions:,";
        std::vector<double> positionVector = hngd.returnSample()->returnPosition();

        // Defining the positions where the values will be printed
        for(int i=0; i<nbNodes; i++)
          *(listPosPrint+i) = 0 ;

        if(nbPosPrint==nbNodes){
            for(int i=0; i<nbPosPrint; i++)
              *(listPosPrint+i) = i ;
        }
        else{
            for(int i=1; i<nbPosPrint-1; i++){
              while(positionVector[listPosPrint[i]] <= i*positionVector[nbNodes-1]/nbPosPrint)
                *(listPosPrint+i) += 1 ;
            }
            *(listPosPrint+nbPosPrint-1) = nbNodes-1;
        }
        // Print the positions as column head
        for(int j=0; j<nbOutput; j++){
          for(int i=0; i<nbPosPrint; i++){
            output << positionVector[listPosPrint[i]] << "," ;
          }
          output << ',' ;
        }
        output << "\n" ;

        output << "dt (s),Time (s)," ;
        string listOutputNames[nbOutput];
        /* HERE */
        listOutputNames[0] = "Ctot," ;
        listOutputNames[1] = "Css,"  ;
        listOutputNames[2] = "Cprec," ;
        listOutputNames[3] = "rD," ;
        listOutputNames[4] = "rN,"    ;
        listOutputNames[5] = "rG," ;
//        listOutputNames[5] = "TSSd," ;
//        listOutputNames[6] = "Kd," ;
//        listOutputNames[7] = "Kn," ;
//        listOutputNames[8] = "J," ;
        for(int j=0; j<nbOutput; j++){
          for(int i=0; i<nbPosPrint; i++) // Print the type of output as column "sub head"
            output << listOutputNames[j] ;
          output << ',' ;
        }
        output << '\n' ;

  }
  else //
  {
        output << "Time(s),Temp(K),TSSp,TSSd,Css,Chyd,Ctot" << endl;
  }

  output.close();
}
    
