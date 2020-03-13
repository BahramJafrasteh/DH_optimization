#ifndef CNdlInterfaces_H
#define CNdlInterfaces_H
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
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

const double MINVERSION = 0.1;
const double VERSION = 0.1;
inline string itoaa(long value);
 inline void tokenise(vector<string>& tokens,
			  const string& str,
			  const string& delimiters);
 
  
using namespace std;
class CStreamInterface 
{
public:
  virtual ~CStreamInterface() {}
  virtual void toStream(ostream& out) const
  {
    out << setiosflags(ios::fixed);
    out << setprecision(6);
    writeToStream(out, "version", getCurrentVersion());
    out << setiosflags(ios::scientific);
    out << setprecision(17);
//     WriteToStream(out);
  }
  
  virtual void toStream(ostream& out, const string fileName) const
  {
    
    WriteToStream(out, fileName);
  }  
  
  
  static double readVersionFromStream(istream& in) 
  {
    double ver = readDoubleFromStream(in, "version");
    if(ver<getMinCompatVersion())
      cout << "Error \n";
    return ver;

  }
  static int readIntFromStream(istream& in, const std::string fieldName)
  {
    string str = readStringFromStream(in, fieldName);
    return atol(str.c_str());
  }
  static double readDoubleFromStream(istream& in, const std::string fieldName)
  {
    string str = readStringFromStream(in, fieldName);
    return atof(str.c_str());
}
  static bool readBoolFromStream(istream& in, const std::string fieldName)
  {
    string str = readStringFromStream(in, fieldName);
    if(atol(str.c_str())!=0)
      return true;
    else
      return false;
  }
  static vector<int> readVectorIntFromStream(istream& in, const std::string fieldName)
  {
    vector<int> vals;
    string str = readStringFromStream(in, fieldName);
    vector<string> tokens;
    tokenise(tokens, str, " ");
    if(tokens.size()==0)
      cout <<  "Zero length vector<int>. \n";
    for(size_t i=0; i<tokens.size(); i++)
      vals.push_back(atol(tokens[i].c_str()));
    return vals;
  }
  static vector<unsigned int> readVectorUintFromStream(istream& in, const std::string fieldName) 
  {
    vector<unsigned int> vals;
    string str = readStringFromStream(in, fieldName);
    vector<string> tokens;
    tokenise(tokens, str, " ");
    if(tokens.size()==0)
      cout <<  "Zero length vector<int>. \n";
      for(size_t i=0; i<tokens.size(); i++)
      vals.push_back(atol(tokens[i].c_str()));
    return vals;
  }
  static string readStringFromStream(istream& in, const std::string fieldName)
  {
    string line;
    vector<string> tokens;
    getline(in, line);
    tokenise(tokens, line, "=");
    if(tokens.size()!=2 || tokens[0]!=fieldName)
      cout <<  "Error. \n";
    return tokens[1];
}
  static void writeToStream(ostream& out, const std::string fieldName, const vector<int> val)
  {
    out << fieldName << "=";
    for(size_t i=0; i<val.size()-1; i++)
      out << itoaa(val[i]) << " ";
    out << itoaa(val[val.size()-1]) << endl;
  }
  
  
  
  
  static void writeToStream(ostream& out, const std::string fieldName, const vector<unsigned int> val)
  {
    out << fieldName << "=";
    for(size_t i=0; i<val.size()-1; i++)
      out << itoaa(val[i]) << " ";
    out << itoaa(val[val.size()-1]) << endl;
  }
  static void writeToStream(ostream& out, const std::string fieldName, const int val)
  {
    out << fieldName << "=" << itoaa(val) << endl; 
  }
  static void writeToStream(ostream& out, const std::string fieldName, const unsigned int val)
  {
    out << fieldName << "=" << itoaa(val) << endl; 
  }
  static void writeToStream(ostream& out, const std::string fieldName, const double val)
  {
    out << fieldName << "=" << val << endl;
  }
  static void writeToStream(ostream& out, const std::string fieldName, const bool val)
  {
    out << fieldName << "=";
    if(val)
      out << "1" << endl;
    else
      out << "0" << endl;
  }
  static void writeToStream(ostream& out, const std::string fieldName, const std::string val)
  {
      out << fieldName << "=" << val << endl;
  }
  static void writeToStream(ostream& out, const std::string fieldName, const char* val)
  {
      out << fieldName << "=" << val << endl;
  }
  
  virtual void fromStream(istream& in, const string fileName)
  {
    readVersionFromStream(in);
    ReadFromStream(in, fileName);
  }

  double getCurrentVersion() const 
  {
    return VERSION;
  }
  static double getMinCompatVersion()
  {
    return MINVERSION;
  }
  virtual void WriteToStream(ostream& out, const string fileName) const = 0;
  virtual void ReadFromStream(istream& out, const string fileName) = 0;
  void toFile(const string fileName)
  {
    const string comment="";
    ofstream out(fileName.c_str());
    if(!out) cout << "Error \n" ;
    if(comment.size()>0)
      out << "# " << comment << endl;
    toStream(out);
    out.close();
  }
  void fromFile(const string fileName) 
  {
    ifstream in(fileName.c_str());
    if(!in.is_open())
      cout << "Error \n";
    try
    {
      fromStream(in, fileName);
    }
    catch(exception& e)
    {
    } 
    in.close();
  }
    
};


inline string itoaa(long value)

{
    enum { kMaxDigits = std::numeric_limits<long>::digits };
    std::string buf;
    buf.reserve( kMaxDigits ); // Pre-allocate enough space.
    int base = 10;
    long quotient = value;
    do 
    {
      buf += "0123456789abcdef"[ std::abs( quotient % base ) ];
      quotient /= base;
    } 
    while ( quotient );
    // Append the negative sign for base 10
    if ( value < 0 && base == 10) buf += '-';
    std::reverse( buf.begin(), buf.end() );
    return buf;
  }
 inline void tokenise(vector<string>& tokens,
			  const string& str,
			  const string& delimiters)
{
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  string::size_type pos = str.find_first_of(delimiters, lastPos);
  
  while (string::npos!=pos 
	 || string::npos!=lastPos)
  {
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);
    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);
  }
}

#endif
