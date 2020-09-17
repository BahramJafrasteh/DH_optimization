#include "DH_OpObj.h"
#include "DH_OpAlg.h"
// #include "Backup/DHO.h"
// #include "Backup2/DHO.h"
const double ONE_OVER_SQRT_2PI = 0.39894228040143267793994605993438;
DHO::DHO(): Loc_optimise()
{
    _init();
}
void DHO::_init()
{
    
}

DHO::DHO( mat &Xin,
          mat Xc,
          int numObj,
          int Dimension,
          int FuncT,
          double Dist_DHs,
          double minAz,
          double maxAz,
          double dist_Az,
          double MinDip,
          double MaxDip,
          double dist_dip,
          double thrsh_ei,
          string Mode,
          int *data_size,
          int disp_level):Xtr(Xin),
          Xcollars(Xc), thrsh_ei(thrsh_ei)
          {
              setThrsh_ei(thrsh_ei);
              setNumData(data_size[0]);
              setInpD(data_size[1]);
              setOutD(1);
              setFuncType(FuncT);
              setNumObj(numObj);
              setOptNumParams(Dimension);
              setDispLevel(disp_level);
              
              setM(Mode);
              init();
              StatisticsCalc();
              //    prepareData();???????????????????????????????????
              mat d = max(Xtr, 0);
              Rd = sqrt( pow(abs(max(Xtr, 0)), 2) + pow(abs(min(Xtr, 0)), 2) );
              mat Rdpar = Rd.submat(0, 0, 0, 2);
              maxRd = Rdpar.max();
              Xcollars_sub = Xc.cols(0, 1);
              
              int indore = 3;
              int indstd = 4;
              MinYh = Xtr.col(indore).min();
              difMaxMinYh = Xtr.col(indore).max()- MinYh;
              
              MinYs = Xtr.col(indstd).min();
              difMaxMinYs = Xtr.col(indstd).max()- MinYs;

              
              
              
              
              int Nr = ((MaxData[0]-MinData[0])/Dist_DHs);
              int Nc = ((MaxData[1]-MinData[1])/Dist_DHs);
              int Naz = (maxAz - minAz)/dist_Az;
              int Ndp = (MaxDip - MinDip)/dist_dip + 1;
              int NCandidates = Nr*Nc*Naz*Ndp;
              
              int ns = 0;
              
              CandidateDHs.resize(NCandidates, 10);
              
              cout << "N Candidate DHs: "<< NCandidates << "\n";
              
              for (int i = 0; i < Nr; i++)
              {
                  
                  for (int j = 0; j < Nc; j++)
                  {
                      for (int k = 0; k < Naz; k++)
                      {
                          for (int l=0; l < Ndp; l++)
                          {
                              CandidateDHs(ns,0) = MinData[0] + Dist_DHs*i;
                              CandidateDHs(ns,1) = MinData[1] + Dist_DHs*j;
                              mat distancia = sum( pow(CandidateDHs(ns,0) - Xcollars.col(0), 2) + pow(CandidateDHs(ns,1) - Xcollars.col(1), 2), 1);
                              double ind_low = distancia.index_min();
                              CandidateDHs(ns,2) = Xcollars(ind_low, 2);
                              CandidateDHs(ns, 3) = minAz + dist_Az*k;
                              CandidateDHs(ns, 4) = MinDip + dist_dip*(double)l;
                              cout << CandidateDHs(ns,0) << "\t" << CandidateDHs(ns,1) << "\t" << Xcollars(ind_low, 2) << "\t" << CandidateDHs(ns,3) << "\t" << CandidateDHs(ns,4) << "\t Sample Number: "<< ns<< "\n";
                              ns += 1;
                          }
                      }
                  }
              }
              
              cout << "Finished\n";
          }
          
          
          
          void DHO::init()
          {
              OptPars.resize(1, getOptNumParams());
              if ( getFuncType() == DH){
                  params.resize(getInpD()+getOutD(), 2);
                  MinData.resize(getInpD()+getOutD(), 1);
                  MaxData.resize(getInpD()+getOutD(), 1);
                  StData.resize(getInpD()+getOutD(), 1);
                  MeanData.resize(getInpD()+getOutD(), 1);
              }
          }
          void DHO::WriteToStream(ostream& out, const string FileName) const
          {
              
              out << "Optimization algorithm : "<< getDefaultOptimiserStr()<<endl;
              out << "Dimension of optimization = " << getOptNumParams() << endl;
              out << "Number of objective functions = " << getNumObj() << endl;
              out << "optimization function : "<< getFuncTypeStr()<<endl;
              int inddot = FileName.find(".");
              string outp = FileName;
              if (inddot < 0)
              {
                  outp = FileName + "_pops.csv";
                  ParetoPops.save(outp, csv_ascii);
                  outp = FileName + "_fits.csv";
                  ParetoFits.save(outp, csv_ascii);
                  if ( getFuncType() == DH){
                      outp = FileName + "_lengths.csv";
                      ParetoLengthDH.save(outp, csv_ascii);}
              }
              else
              {
                  string FName = FileName.substr(0, inddot);
                  outp = FName + "_pops.csv";
                  ParetoPops.save(outp, csv_ascii);
                  outp = FName + "_fits.csv";
                  ParetoFits.save(outp, csv_ascii);
                  if ( getFuncType() == DH){
                      outp = FName + "_lengths.csv";
                      ParetoLengthDH.save(outp, csv_ascii);}
              }
              
          }
          
          void DHO::ReadFromStream(istream& in, const string FileName)
          {
              setModelType(readIntFromStream(in,  "NN Moldel :"));
              setDefaultOptimiser(readIntFromStream(in,  "Optimisation Algorithm :"));
              setNumObj(readIntFromStream(in,  "Number of Hidden Nodes :"));
              setInpD(readIntFromStream(in,  "Input Dimension :"));
              setOutD(readIntFromStream(in,  "Output Dimension :"));
              setPrepM( readIntFromStream(in,  "Preparation Method :") );
              MinTotalin = readDoubleFromStream(in, "MinTotalin :");
              MaxTotalin = readDoubleFromStream(in, "MaxTotalin :");
              MinTotalo = readDoubleFromStream(in, "MinTotalo :");
              MaxTotalo = readDoubleFromStream(in, "MaxTotalo :");
              setNumData(readIntFromStream(in , "Number of Data :"));
              setOptNumParams(readIntFromStream(in , "Number of Parameters :"));
              init();
              string inp = FileName+ "_OptPars.txt";
              OptPars.load(inp, csv_ascii);
              setOptParams(OptPars);
              Xtr.resize(getNumData(), getInpD());
              mat NXtr;
              mat Statistics;
              inp = FileName + "_train.txt";
              NXtr .load(inp , csv_ascii);
              Xtr = NXtr;
              inp = FileName + "_Statistics.txt";
              Statistics.load(inp ,csv_ascii);
              params = Statistics.cols(0, 1);
              MinData = Statistics.col(2);
              MaxData = Statistics.col(3);
              MeanData = Statistics.col(4);
              StData = Statistics.col(5);
          }
          DHO *ReadMFromFile(const string MFileName)
          {
              ifstream in(MFileName.c_str());
              if(!in)
              {cout << "Reading Model Aborted.\n" ;
                  exit(0);}
                  DHO* pmodel;
                  try
                  {
                      pmodel = new DHO();
                      pmodel->fromStream(in, MFileName);
                  }
                  catch (exception& er)
                  {
                  }
                  in.close();
                  return pmodel;
          }
          void WriteMToFile(const DHO& model, const string FileName)
          {
              int inddot = FileName.find(".");
              string FName = FileName;
              if (inddot > 0)
                  FName = FileName.substr(0, inddot) + "_info.txt";
              else
                  FName = FileName + "_info.txt";
              const string comment="";
              ofstream out(FName.c_str());
              if(!out)
              {
                  cout << "Writing Model Aborted.\n" ;
                  exit(1);
              }
              if(comment.size()>0)
                  out << "# " << comment << endl;
              model.toStream(out, FileName);
              out.close();
          }
          
          void DHO::optimise_al()
          {
              Loc_optimise::setDispLevel(getDispLevel());
              if ( getFuncType() == DH)
                  runDefaultOptimiser(CandidateDHs);
          }
          
          mat DHO::computeObjectiveVals() const
          {
              bool infeasible = false;
              int Dimension = getOptNumParams();
              
              double Length_drill = 0.0;
              mat Func_EI_L(1,5);
              switch ( getFuncType() )
              {
                  case DH:
                  {
                      int indore = 3;
                      int indstd = 4;
                      int indexes = 5;
                      int indx = 0;
                      int indy = 1;
                      int indz = 2;
                      uvec indloc = {0, 1, 2};
                      double md=20;//maximum distance from each block
                      double part = 5;
                      
                      double maxDepth = 200;
                      double minDepth = 100;
                      infeasible = false;
                      double thrsh = 5;
                      double azimuth = OptPars[3];
                      double diip=M_PI / 2- OptPars[4]; //pi/2-dip;
                      mat dist = pow(OptPars.submat(0,0,0,1) - Xcollars_sub.each_row(), 2);
                      mat distc = sqrt(sum(dist, 0));
                      int indc = distc.index_min();
                      
                      double     az = azimuth;
                      double dip = diip;
                      vec rad1;
                      vec rad2;
                      vec rad3;
                      if (dip>=0)
                      {
                          rad1 = regspace<vec>(0,  part, 5*maxRd);
                          //rad2 = regspace<vec>(0,  part,  Rd[1]);
                          //rad3 = regspace<vec>(0,  part,  Rd[2]);
                      }
                      else
                      {
                          rad1 = regspace<vec>(-5*maxRd,  -part,  0);
                          //rad2 = regspace<vec>(-Rd[1],  -part,  0);
                          //rad3 = regspace<vec>(-Rd[2],  -part,  0);
                      }
                      //effect of dip and azimuth
                      vec newx = OptPars[0] + rad1*sin(az) *cos(dip);
                      vec newy = OptPars[1] + rad1* cos(az) *cos(dip);
                      vec newz = OptPars[2] - rad1*sin(dip);
                      //split data
                      //bring the data into upperbound and lowerbound
                      uvec indn = newx <= Xtr.col(indx).max();
                      indn %= newx >= Xtr.col(indx).min();
                      indn %= newy <= Xtr.col(indy).max();
                      indn %= newy >= Xtr.col(indy).min();
                      indn %= newz <= Xtr.col(indz).max();
                      indn %= newz >= Xtr.col(indz).min();
                      
                      // return indices of non-zero elements
                      uvec inds = find(indn, 0);
                      
                      vec newx1 = newx(inds);
                      vec newy1 = newy(inds);
                      vec newz1 = newz(inds);
                      mat unz = unique<vec>(newz1);
                      if (unz.size()>1)
                      {
                          mat nexy = newx1;
                          nexy = join_horiz(nexy, newy1);
                          nexy = join_horiz(nexy, newz1);
                          
                          //maximum length of each drill hole
                          double maxl = sqrt ( accu( pow ( nexy.submat(0, 0, 0, indz)
                          - nexy.submat(nexy.n_rows-1, 0, nexy.n_rows-1, indz), 2 ) ) );
                          
                          
                          //( nexy.col(2).max() - nexy.col(2).min() ) * cos( abs(dip) );
                          /* while (maxl > maxDepth)
                           *    {
                           *      double indmax = nexy.col(indz).index_max();
                           *      nexy.shed_row(indmax);
                           *      maxl = sqrt ( accu( pow ( nexy.submat(0, 0, 0, indz)
                           *    - nexy.submat(nexy.n_rows-1, 0, nexy.n_rows-1, indz), 2 ) ) );
                      }*/
                          //     uvec maxl_remove = find(maxl < maxDepth);
                          mat Xyz = Xtr.col(indx);
                          Xyz = join_horiz( Xyz, Xtr.col(indy) );
                          Xyz = join_horiz( Xyz, Xtr.col(indz) );
                          
                          //     mat nexy1 = nexy.rows(maxl_remove);
                          mat nexy1;
                          int m = 0;
                          for (int i = 0; i < nexy.n_rows; i++)
                          {
                              mat distancia = sqrt( sum( pow(nexy.row(i) - Xyz.each_row(), 2), 1) );
                              //distance below a threshold
                              double ind_low = distancia.index_min();
                              mat temp = Xtr.row(ind_low);
                              nexy1.insert_rows(i, temp);
                              
                              /*uvec ind_low = find( distancia < md );
                               *      if (accu(ind_low) > 0){
                               *      mat temp = Xtr.rows(ind_low);
                               *      mat tt(1,1) ;
                               *      tt[0] = accu ( temp.col(4) ) / temp.n_rows;
                               *      nexy1.insert_rows(m, tt);
                               *      m++;
                          }
                          else {
                              nexy.shed_row(i);
                          }*/
                              /*if (i == 0)
                               *	nexy1 = accu ( temp.col(4) ) / temp.n_rows;
                               *      else
                               *	nexy1 = join_horiz( nexy1, accu ( temp.col(4) ) / temp.n_rows );*/
                              
                          }
                          // number of blocks
                          
                          uvec unq = find_unique(nexy1.col(indexes));
                          nexy1 = nexy1.rows(unq);
                          // fourth column is estimated value
                          
                          
                          int NOB = nexy1.n_rows;
                          if (NOB > 0)
                          {
                              // by default the forth column is ore grade value and the fifth one is its
                              
                              mat Yh_use = (nexy1.col(indore) - MinYh)/difMaxMinYh;
                              // fifth column is its standard deviation
                              mat Yst_use = (nexy1.col(indstd)-  MinYs)/ difMaxMinYs;
                              arma::rowvec smYh = median(Yh_use);
                              arma::rowvec smYst = median(Yst_use);
                              double mYh = smYh[0]; 
                              double mYst = smYst[0];
                              double thrsh_ei = getThrsh_ei();
                              
                              Length_drill = sqrt ( accu( pow ( nexy1.submat(0, 0, 0, indz)
                              - nexy1.submat(nexy1.n_rows-1, 0, nexy1.n_rows-1, indz), 2 ) ) );
                              Func_EI_L[0] = Length_drill;
                              
                              //PI computiation
                              Func_EI_L[1]   = normcdf(mYh, thrsh_ei, mYst );
                              
                              // EI computation
                              Func_EI_L[2]   = ( mYh - thrsh_ei ) * normcdf(mYh, thrsh_ei, mYst ) + mYst * normpdf(mYh, thrsh_ei, mYst );
                              Func_EI_L[3] = mYh;
                              Func_EI_L[4] = mYst;
                              
                              cout << "thrsh_ei: "<< thrsh_ei << "\t";
                              cout << "mYh: "<< mYh << "\t";
                              cout << "mYst: "<< mYst << "\t";
                              cout << "PI: "<< Func_EI_L[1] << "\t";
                              cout << "EI: "<< Func_EI_L[2] << "\n";
                              
                              
                              //if (Length_drill < minDepth || Length_drill > maxDepth)
                              //  infeasible = true;
                              //       nexy1.save("dd.txt", csv_ascii);
                              
                          }
                          else
                              infeasible = true;
                      }
                      else
                          infeasible = true;
                      if (infeasible)
                      {Func_EI_L[0] = INFINITY;
                          Func_EI_L[1] = INFINITY;
                      }
                      return Func_EI_L;
                  }
                  break;
              }
          }
          
          mat DHO::getmaxmindata() const {
              
              mat maxmindata = join_horiz<>(MinData, MaxData);
              return maxmindata;
          }
          
          void DHO::setOptParams(const mat param) const
          {
              OptPars = param;
          }
          
          
          
          mat DHO::getOptParams() const
          {
              mat param(1, getOptNumParams());
              param = OptPars;
              return (param);
          }
          
          
          
          
          
