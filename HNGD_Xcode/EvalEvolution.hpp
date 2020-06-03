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
        bool evaluate(vector<double> profile) ;
        void setCriterion(double value);
        void setProfile(vector<double> profile);
        
    private:
        double _criterion ;
        vector<double> _mem_profile ;
};

#endif /* EvalEvolution_hpp */
