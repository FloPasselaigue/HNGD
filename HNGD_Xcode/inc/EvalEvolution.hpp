/**
    This class evaluates the average evolution of a given profile
    to determine if the current state should be written in the output
 */

#ifndef EvalEvolution_hpp
#define EvalEvolution_hpp

#include <stdio.h>
#include <vector>

using namespace std ;

class EvalEvolution
{
    public:
        EvalEvolution();
        EvalEvolution(double criterion);
    
        // Determine if the profile has evolved more than the criterion
        bool evaluate(vector<double> profile) ;
    
        void setCriterion(double value);
        void setProfile(vector<double> profile);
        
    private:
        double _criterion ; // Criterion of evolution
        vector<double> _mem_profile ; // Reference used to evaluate the evolution
};

#endif /* EvalEvolution_hpp */
