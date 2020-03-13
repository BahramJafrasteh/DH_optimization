#ifndef Loc_optimise_H
#define Loc_optimise_H
#include<vector>
#include <stdlib.h>
#include <sstream>
#include <stdio.h>
#include <iostream> 
#include <fstream>
#include <string>
#include <math.h>
#include <limits>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstring>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <climits>
#include <cfloat>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <cfloat>
#include "lapack.h"
#include<armadillo>
#include<random>

using namespace std;
using namespace arma;
inline int sgn(double val) {
    if (val < 0) return -1;
    if (val >= 0) return 1;
    return 1;
}
class Loc_optimise {
public :
    enum 
    {
        CompleteSearch,
        ABC
    };  
    Loc_optimise();
    ~Loc_optimise() {}
    
    void setDispLevel(int disp_level)
    {
        DispLevel = disp_level;
    }
    int getDispLevel() const
    {
        return(DispLevel);
    }
    void init( int Dimension)
    {
        
        //   NaND = numeric_limits<double>::quiet_NaN();
        CountFuncEval = 0;
        loop = 0;
        
    }
    
    string GetStdoutFromCommand(string cmd) {
        
        string data;
        FILE * stream;
        const int max_buffer = 256;
        char buffer[max_buffer];
        cmd.append(" 2>&1");
        
        stream = popen(cmd.c_str(), "r");
        if (stream) {
            while (!feof(stream))
                if (fgets(buffer, max_buffer, stream) != NULL) data.append(buffer);
                pclose(stream);
        }
        return data;
    }
    
    // tested-working
    virtual mat getOptParams() const = 0;
    
    virtual mat getmaxmindata() const = 0;
    virtual mat computeObjectiveVals() const=0;
    
    
    
    virtual void setOptParams(const mat param) const = 0;
    
    virtual int getOptNumParams() const = 0;
    
    double HyperVolume(int num_rep, int N, const mat fitness, const mat maxF);
    double calculateIGD(const mat fitness);
    double min_el(mat &X, int &R);
    
    bool Dominate(const mat FitX, int IndXi, 
                  const mat FitY, int indYi);
    
    
    void Complete_Search();
    void eval_FIS(const double Diversity, const double Err, 
                  double &C1,double &C2, double &W);
    // sort integers
    void Sort_integer(uvec &Sort, uvec &SortInd, const string mode = "A");
    // sort doulbes and floats
    void Sort_double( vec& Sort, uvec&SortInd, const string mode = "A");
    inline double rando()
    {
        return (double)rand() / ((double)RAND_MAX + 1);
    }
    
    mat XCandidate;
    void runDefaultOptimiser(mat CandidateDHs)
    {
        
        //FuncType = "DH";
        XCandidate = CandidateDHs;
        
        switch(DfOptmiser)
        {
            case CompleteSearch:
                Complete_Search();
                break;
            default:
                if (getDispLevel() > 0)
                    cout <<"Unknown optimisation algorithm. \n";
        }
    }
    
    
    void ChkBnd (double &A, double lb, double ub, int index_location)
    {
        cout << "Undefined \n";
        exit(1);
    }
    void ChkBnd (double &A, double lb, double ub)
    {
        if (A < lb)
            A = lb;
        if (A> ub)
            A = ub;
    }
    void ChkBnd (mat &A, mat lb, mat ub)
    {
        uvec indl = find(A < lb);
        A.elem(indl) = lb.elem(indl);
        uvec indu = find(A > ub);
        A.elem(indu) = lb.elem(indu);
        
    }
    
    void setFuncEvalTerminate(bool val)
    {
        funcEvalTerminate = val;
    }
    
    bool isFuncEvalTerminate() const
    {
        return funcEvalTerminate;
    }
    void setHybridLearning(bool val)
    {
        HybridLearn = val;
    }
    bool isHybridLearning() const
    {
        return HybridLearn;
    }
    void setObjectiveTolTerminate(bool val)
    {
        objectiveTolTerminate = val;
    }
    
    bool isObjectiveTolTerminate() const
    {
        return objectiveTolTerminate;
    }
    void setIterTerminate(bool val)
    {
        iterTerminate = val;
    }
    bool isIterTerminate() const
    {
        return iterTerminate;
    }
    void setMaxFuncEvals(unsigned int val)
    {
        maxFuncEval = val;
    }
    unsigned int getMaxFuncEvals() const
    {
        return maxFuncEval;
    }
    void setMaxIters(unsigned int val)
    {
        maxIters = val;
    }
    double getMaxIters()
    {
        return maxIters;
    }
    
    void setLb(mat& val)
    { 
        Lb = val;
    }
    mat getLb()
    {
        return Lb;
    }
    void setUb(mat& val)
    { 
        Ub = val;
    }
    mat getUb()
    {
        return Ub;
    } 
    
    void setObjectiveTol(double val)
    {
        objectiveTol = val;
    }
    double getObjectiveTol() const
    {
        return objectiveTol;
    }
    void setNPop(double val)
    {
        NPop = val;
    }
    double getNPop() const
    {
        return NPop;
    }
    void setParamTol(double val)
    {
        parameterTol = val;
    }
    double getParamTol() const
    {
        return parameterTol;
    }
    void setDefaultOptimiser(int val) const
    {
        DfOptmiser = val;
    }
    int getDefaultOptimiser() const
    {
        return DfOptmiser;
    }
    void setDefaultOptimiserStr(string val)
    {
        if(val == "CompleteSearch")
            DfOptmiser = CompleteSearch;
    }
    string getDefaultOptimiserStr() const
    {
        switch(DfOptmiser)
        {
            case CompleteSearch:
                return "CompleteSearch";
        }
    }
    
    void setNumObj(int Nh)
    {
        NumObjs = Nh;
    }
    int getNumObj() const
    {
        return(NumObjs);
    }
    void setInD(unsigned int dim)
    {
        InD = dim;
    }
    // Get the input dimension.
    inline unsigned int getInD() const
    {
        return InD;
    }
private:
    double objectiveTol;
    double parameterTol;
    
    int NPop;
    int NumObjs;
    int modelT;
    int InD;
    int DispLevel;
    
    mat Lb;
    mat Ub;
    mat TruePF;
    mat Final_EI;
    
    
    unsigned int maxIters;
    unsigned int OptNumParams;
    unsigned int maxFuncEval;
    
    
    mutable int DfOptmiser;
    mutable int loop;
    mutable int CountFuncEval;
    
    mutable mat ParetoF_fit;
    mutable mat ParetoF_pop;
    mutable mat ParetoF_vel;
    
    mutable double OldGlobalBest_fit;
    mutable double GlobalBest_fit;
    
    bool funcEvalTerminate;
    bool iterTerminate;
    bool objectiveTolTerminate;
    bool HybridLearn;
    
    string FuncType;
    
    void _init();
    
};
# endif
