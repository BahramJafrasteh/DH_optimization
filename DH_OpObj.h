#ifndef _DHO_
#define _DHO_
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif
#include "DH_OpAlg.h"
#include "lapack.h"
#include "ndlfortran.h"
#include <armadillo>
#include "CInterface.h"
using namespace std;
using namespace arma;
class DHO: public Loc_optimise, public CStreamInterface
{
public:
    enum PrepM
    {
        P1,
        P2,
        P3
    };
    enum FuncType
    {
        DH
    };
    
    DHO();
    DHO( mat &Xin, mat Xc, int numObj = 3,
         int Dimension = 4, int FuncT= DH,  double Dist_DHs = 20.0,
         double minAz = 0*M_PI,
         double maxAz = 360*M_PI/180,
         double dist_Az = 20*M_PI/180,
         double MinDip = 40*M_PI/180,
         double MaxDip = 90*M_PI/180,
         double dist_dip = 10*M_PI/180,
         double thrsh_ei = 2,
         string Mode = "Available Data", int *data_size = 0, int disp_level = 0);
    
    void init(); 
    void setDispLevel(int disp_level)
    {
        DispLevel = disp_level;
    }
    int getDispLevel() const
    {
        return(DispLevel);
    }
    
    void tDim(int Opt_Dim)
    {
        DimOpt = Opt_Dim;
    }
    int getOptDim() const
    {
        return(DimOpt);
    }
    
    void setThrsh_ei(double thrsh)
    {
        thrsh_ei = thrsh;
    }
    double getThrsh_ei() const
    {
        return(thrsh_ei);
    }    
    
    void setNumObj(int Nh)
    {
        NumObj = Nh;
    }
    int getNumObj() const
    {
        return(NumObj);
    }
    void setM(std::string val) 
    {
        mode = val;
    }
    inline std::string getM() const
    {
        return mode;
    }
    inline void setInpD(unsigned int dim)
    {
        inpD = dim;
    }
    // Get the input dimension.
    inline unsigned int getInpD() const
    {
        return inpD;
    }
    // Set the output dimension.
    inline void setOutD(unsigned int dim)
    {
        outD = dim;
    }
    // Get the output dimension.
    inline unsigned int getOutD() const
    {
        return outD;
    }
    void setPrepM(const int val)
    {
        PrepM = val;
    }
    
    
    void setFuncType(int val) const
    {
        FuncType = val;
    }
    int getFuncType() const
    {
        return FuncType;
    }
    string getFuncTypeStr() const
    {
        if ( getFuncType() == DH)
            return "DH";
    }
    //virtual void getM() const = 0;
    mat getOptParams() const;
    
    void setOptParams(const mat param) const;
    void setParetoPops(const mat param) const;
    void setParetoFits(const mat param) const;
    void setParetoLengthDH(const mat param) const;
    mat getmaxmindata() const;
    mat computeObjectiveVals() const;
    void WriteParams() const;
    void WriteToStream(ostream& out, const string FileName) const;
    void ReadFromStream(istream& in, const string FileName);
    
    //mat normpdf(const mat & x, const mat& u, const mat& sig) const;
    //mat normcdf(const mat & x, const mat& u, const mat& sig) const;
    //mat normcdf(double & x, const mat& u, const mat& sig) const;
    void StatisticsCalc()
    {
        MaxTotalin = Xtr.max();
        MinTotalin = Xtr.min();
        for (int i = 0; i < (getInpD() + getOutD()); i++)
        {
            MinData[i] = Xtr.col(i).min();
            MaxData[i] = Xtr.col(i).max();
            MeanData[i] = accu(Xtr.col(i))/getNumData();
            StData[i] = sqrt( accu( pow ( Xtr.col(i) - MeanData[i], 2) )/ (getNumData()-1) );
        }
    }
    int getOptNumParams() const
    {
        return NumParam;
    }
    int setNumData(const int val) 
    {
        NumData = val;
    }
    int getNumData() const
    {
        return NumData;
    }
    
    
    void setOptNumParams(unsigned int val)
    {
        NumParam = val;
    }
    void setModelType(const int val)
    {
        modelType = val;
    }
    int getModelType() const
    {
        return(modelType);
    }
    void setModelStr(const string val) 
    {
    }
    string getModelStr() const 
    {
    }
    void optimise_al();
    
    
    void none();
    void none(mat&X);
    void getNumHn();
    //vector<vector<double> >Xtr;
    mat Xtr;
    mat Xcollars;
    mat Xcollars_sub;
    mat TruePF;
    mat Rd;
    
    double maxRd;
    mutable mat onesmat;
    mutable mat CandidateDHs;
    mutable mat OptPars;
    mutable mat ParetoFits;
    mutable mat ParetoPops;
    mutable mat ParetoLengthDH;
    //vector<double>Ytr;
    //double* Ytr;
    int DispLevel;
    mutable double HV;
    mutable double IGD;
    mutable int modelType;
    int NumObj;
    mutable int inpD;
    mutable int outD;
    mutable int FuncType;
    mutable int PrepM;
    mutable bool YScale;
    mutable string mode;
    mutable int NumParam;
    mutable int NumData;
    
    double fact = 1.0;
    double thrsh_ei;
    
    mutable mat params;
    mutable mat MinData;
    mutable mat MaxData;
    mutable mat MeanData;
    mutable mat StData;
    mutable double MaxTotalin;
    mutable double MinTotalin;
    mutable double MaxTotalo;
    mutable double MinTotalo;
    int DimOpt;
    
private:
    void _init();
    
    
    
};
void WriteMToFile(const DHO& model, const string FileName);
DHO *ReadMFromFile(const string MFileName);
#endif // _DHO_
