#ifndef InOut_hpp
#define InOut_hpp

#include <stdio.h>
#include "HNGD.hpp"

using namespace std;

class InOut
{
    public:
        static void getSettings(int nb, double* settings, string path_exec, string file_name);
    
        static void getPhysics(int nb, double* physics, string path_exec, string file_name);
    
        static vector<vector<double>> getThermalTreatment(string path_exec, string file_name);
    
        static vector<vector<double>> getICHydrogen(string path_exec, string file_name);
    
        static void writeSettingsInCheck(double * settings);
        
        static void writePhysicsInCheck(double * physicalParameters);

        static void writeOuput(HNGD hngd, string path_exec, string output_name, int nbNodes, int nbOutput, double t, double temp, int nbPosPrint, int* listPosPrint);
            
        static void writeInitialOutput(HNGD hngd, string path_exec, string output_name, int nbNodes, int nbOutput, int nbPosPrint, int* listPosPrint) ;
};

#endif /* InOut_hpp */
