#include "DH_OpAlg.h"

 using namespace arma;
Loc_optimise::Loc_optimise()
{
  _init();
}
void Loc_optimise::_init()
{
  setFuncEvalTerminate(true);
  setObjectiveTolTerminate(true);
  setIterTerminate(true);
  setHybridLearning(false);
  setMaxFuncEvals(20000);
  setMaxIters(300);
  setObjectiveTol(1e-3);
  setParamTol(1e-6);
  setNPop(20);
  setDefaultOptimiser(CompleteSearch);
}
/*Copyright (c) 2009, Yi Cao
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the distribution
*/


void Loc_optimise::Complete_Search()
{
  if (getDispLevel() > 0)
    cout << "// // // // // // // // // // // // // // //Multi-objective AgMOPSO optimisation// // // // // // // // // // // // // // //" << endl;
                   int ns = 0;
  cout << XCandidate(ns,0) << "\t" << XCandidate(ns,1) << "\t" << XCandidate(ns,2) << "\t" << XCandidate(ns,3) << "\t" << XCandidate(ns,4) << "\t Sample Number: "<< ns<< "\n";
  int length = XCandidate.n_rows;     
  cout << length << "\n";
  //Final_EI.resize(length,2);
  int mYh = XCandidate.n_cols-5;
  int mYst = XCandidate.n_cols-4;
  int pi = XCandidate.n_cols-3;
  int ei = XCandidate.n_cols-2;
  int len = XCandidate.n_cols-1;
  
  for (int i =0; i < length; i++)
  {
  setOptParams(XCandidate.row(i));
  mat Fit = computeObjectiveVals();
//Final_EI[i,0] = Fit[0];
  XCandidate(i,len) = Fit[0];
  XCandidate(i,pi) = Fit[1];
  XCandidate(i,ei) = Fit[2];
  XCandidate(i,mYh) = Fit[3];
  XCandidate(i,mYst) = Fit[4];
  
  //Final_EI[i,1] = Fit[1];
  cout << "Number : " << i<< endl;
  //cout << "Number : " << i <<" PI value : " <<XCandidate(i,pi) << ", EI value : " <<XCandidate(i,ei) <<" \n";
  }
//XCandidate.col( XCandidate.n_cols-2) = Final_EI.col(0);
  //XCandidate.col( XCandidate.n_cols-1) = Final_EI.col(1);
  XCandidate.save("Final_Results.txt", arma_ascii);
  
}

void Loc_optimise::Sort_integer(  uvec &Sort, uvec &SortInd, const string mode)
{
  uvec X = Sort;
  int sz= X.n_elem;
  int sizey = 0;
  int rind = 0;
  int s = 0;
  int indminx_tmp;
  double minx_tmp;
  
  vector<double> Y(sz, 0);
  vector<int> indY(sz, 0);
  vector<int> indX(sz, 0);
  vector<int> indselect(sz, 0);// index of selection for each value
  
  for (int i = 0; i <sz; i++)
  {
    Y[i] = X[i];
    indY[i] = i;
    indX[i] = indY[i];
  }
  
  sizey = Y.size();
  while ( sizey > 0 )
  {
    minx_tmp = Y[0];
    indminx_tmp = indY[0];
    for (int j = 0; j < sizey; j++)
    {
      if ( mode == "A")
      {
	if (Y[j] < minx_tmp)
	{
	  minx_tmp = Y[j];
	  indminx_tmp = indY[j];
	  //Sort[rind] = Y[j];
	}
      }
      else if (mode == "D")
      {
	if (Y[j] > minx_tmp)
	{
	  minx_tmp = Y[j];
	  indminx_tmp = indY[j];
	  //Sort[rind] = Y[j];
	}
      }	
      
    }
    SortInd[rind] = indminx_tmp;
    Sort[rind] = minx_tmp;
    //SortInd[indminx_tmp] = rind;
    //indselect[indminx_tmp] = 1;
    indselect[SortInd[rind]] = 1;
    s = 0;
    for (int i = 0; i < sz; i ++)
    {
      if (indselect[i] == 0)
      {
	Y[s] = X[i];
	indY[s] = indX[i];
	s++;
      }
    }
    sizey --;
    rind++;
  }
}


void Loc_optimise::Sort_double( vec& Sort, uvec&SortInd, const string mode)
{  
  vec X = Sort;
  int sz= X.n_elem;
  int sizey = 0;
  int rind = 0;
  int s = 0;
  int indminx_tmp;
  double minx_tmp;
  
  vector<double> Y(sz, 0);
  vector<int> indY(sz, 0);
  vector<int> indX(sz, 0);
  vector<int> indselect(sz, 0);// index of selection for each value
  
  for (int i = 0; i <sz; i++)
  {
    Y[i] = X[i];
    indY[i] = i;
    indX[i] = indY[i];
  }
  
  sizey = Y.size();
  while ( sizey > 0 )
  {
    minx_tmp = Y[0];
    indminx_tmp = indY[0];
    for (int j = 0; j < sizey; j++)
    {
      if ( mode == "A")// ascending
      {
	if (Y[j] < minx_tmp)
	{
	  minx_tmp = Y[j];
	  indminx_tmp = indY[j];
	}
      }
      else if (mode == "D")// descending
      {
	if (Y[j] > minx_tmp)
	{
	  minx_tmp = Y[j];
	  indminx_tmp = indY[j];
	}
      }	
      
    }
    SortInd[rind] = indminx_tmp;
    Sort[rind] = minx_tmp;
    indselect[SortInd[rind]] = 1;
    s = 0;
    for (int i = 0; i < sz; i ++)
    {
      if (indselect[i] == 0)
      {
	Y[s] = X[i];
	indY[s] = indX[i];
	s++;
      }
    }
    sizey --;
    rind++;
  }
  
}
