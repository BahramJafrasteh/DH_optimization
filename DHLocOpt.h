#ifndef _NN_
#define _NN_
#include "DH_OpAlg.h"
#include <stdlib.h>
#include<time.h>
#include<armadillo>
#include "DH_OpObj.h"
using namespace arma;
int main(int argc, char* argv[]);
class CNN{
 public:
  CNN(int argc, char** argv): argc(argc), argv(argv)
  {   argNumber = 1;
  prepareM = 1; // preparation method
  disp_level = 0;
  setSeed();
  }
  void setSeed()
  {
    seed = time(NULL);
    srand(seed);
  }
  bool getline(std::istream& in, std::string& line)
{
  bool val = false;
  do
  {
    if(std::getline(in, line))
      val = true;
  }
  while(line[0]=='#');
  return val;
}

void Search();
void helpInfo();
void helpHeader();
void readDataFile(mat &X, int *data_size, const string fileName);
int* readDataSize(const string fileName);
    bool isCArgF() const
    {
      if(argv[argNumber][0]=='-')
	return true;
      else
	return false;
    }
      bool isCArg(std::string Name) const
    {
      return (argv[argNumber] == Name);
    }
      int getCArgNumber() const
    {
      return argNumber;
    }
  void setCArgNumber(int val)
    {
      argNumber = val;
    }
    
     std::string getCArg() const
    {
      return argv[argNumber];
    }
    int getIntCArg() const
    {
      return atol(argv[argNumber]);
    }
    int getDoubleCArg() const
    {
      return atof(argv[argNumber]);
    }
    int incArg()
    {
      argNumber++;
    }
    

    bool isFlgs() const
    {
      return flgs && getCArgNumber()<argc;
    }
  void setFlgs(bool val) 
    {
      flgs = val;
    }
        int getprepM() const
    {
      return prepareM;
    }
  void setprepM(int val)
    {
      prepareM = val;
    }
    int disp_level;
     private:
  bool flgs;
  unsigned long seed;
  int fileFormat;
  int argNumber;
  int prepareM;
  
  std::string mode;
  void isCurrentArgumentFlag();
 protected:
  int argc; 
  char** argv; 
};
#else /* _NN_ */
#endif /* _NN_ */
