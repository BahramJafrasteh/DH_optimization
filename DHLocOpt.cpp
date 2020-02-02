#include "DHLocOpt.h"
using namespace arma;
int main(int argc, char* argv[])
{
    CNN command(argc, argv);
    command.setFlgs(true);
    //command.setprepM(1);// 0: -1 and 1, 1: symmetric, 2: 0 and 1
    try
    {
        while (command.isFlgs())
        {
            std::string argument = command.getCArg();
            if (argv[command.getCArgNumber()][0] == '-')
            {
                if (command.isCArg("-?"))
                {
                    command.helpInfo();
                    exit(0);
                }
                else if (command.isCArg("-h"))
                {
                    command.helpInfo();
                    exit(0);
                }

                else

                    std::cout<<"Unrecognised flag: " << command.getCArg() << " under the provided command.";

            }
            else if (argument == "Search") // initialise with an old model.
                command.Search();
            command.incArg();
        }
        //if (disp_level > 0)
        //std::cout<< "Please provide a command \n";
    }
    catch (std::exception err)
    {
        return 0;
    }

}

void CNN::Search(){
    const clock_t begin_time = clock();
    // do something
    argNumber++;
    std::string Mode = "DHOptimize";
    disp_level = 2;
    int iters = 300;
    int Func_type = DHO::DH;
    string FuncTypeStr = "DH";
    std::string OptimizationAlgStr = "CompleteSearch";
    std::string modelName = "DH_Opt";
    double Dist_DHs = 20.0;
    double minAz = 0*M_PI;
    double maxAz = 360*M_PI/180;
    double dist_Az = 60*M_PI/180;//20*M_PI/180;
    double MinDip = 40*M_PI/180;
    double MaxDip = 90*M_PI/180;
    double dist_dip = 50*M_PI/180;
    double thrsh_ei = 26.9;
    int numObj = 1;
    int Dimension = 4;
    //string prepMStr = "P1";// P1 P2 P3
    std::string CollarFName;
    std::string outFName = "ParetoFronts.csv";
    int hybridstat;
    while (isFlgs()) {
        if (isCArgF())
        {
            if (isCArg("-?")) {
                helpInfo();
                exit(0);
            }
            else if (isCArg("-h")) {
                helpInfo();
                exit(0);
            }
            else if (isCArg("-dl"))
            {
                incArg();
                disp_level = getIntCArg();
            }
            else if (isCArg("-th")) {
                argNumber++;
                thrsh_ei = getDoubleCArg();
            }
            else if (isCArg("-daz")) {
                argNumber++;
                dist_Az = getDoubleCArg();
            }
            else if (isCArg("-ddh")) {
                argNumber++;
                Dist_DHs = getDoubleCArg();
            }
            else if (isCArg("-ddip")) {
                argNumber++;
                dist_dip = getDoubleCArg();
            }
            else if (isCArg("-mindip")) {
                argNumber++;
                MinDip = getDoubleCArg();
            }
            else if (isCArg("-maxdip")) {
                argNumber++;
                MaxDip = getDoubleCArg();
            }
            argNumber++;
        }
        else
            setFlgs(false);
    }
    string dataFName;
    mat X;
    int* data_size = new int[2];
    mat Xc;// collars data
    if ( FuncTypeStr == "DH")
    {
        string dataFName = argv[argNumber];
        if ((getCArgNumber() + 1) < argc)
            CollarFName = argv[getCArgNumber() + 1];
        data_size = readDataSize(dataFName);
        X.resize(data_size[0], data_size[1]);
        readDataFile(X, data_size, dataFName);
        //outFName = "ParetoFronts_DH.csv";
        if ((getCArgNumber() + 3) < argc)
            outFName = argv[getCArgNumber() + 2];
        data_size = readDataSize(CollarFName);
        Xc.resize(data_size[0], data_size[1]);
        readDataFile(Xc, data_size, CollarFName);



    }
    else{
        if (argNumber < argc)
            outFName = argv[argNumber];
        dataFName =   "ParetoFront/" + FuncTypeStr + ".txt";
        data_size = readDataSize(dataFName);
        X.resize(data_size[0], data_size[1]);
        readDataFile(X, data_size, dataFName);

    }



    mat lb(1, Dimension);
    mat ub(1, Dimension);
    if (FuncTypeStr == "DH")
    {
        Func_type = DHO::DH;
        Dimension = 5;
        // Lower Bound
        lb[0] = X.col(0).min();
        lb[1] = X.col(1).min();
        lb[2] = -180;
        lb[3] = 0;
        // Upper Bound
        ub[0] = X.col(0).max();
        ub[1] = X.col(1).max();
        ub[2] = 180;
        ub[3] = 50;
    }

    DHO* pmodel;

    if ( Func_type == DHO::DH)
        pmodel = new DHO(X, Xc, numObj, Dimension, Func_type, Dist_DHs , minAz, maxAz, dist_Az, MinDip, MaxDip, dist_dip,thrsh_ei, Mode, data_size, disp_level);

    if (OptimizationAlgStr == "CompleteSearch")
        pmodel->setDefaultOptimiser(Loc_optimise::CompleteSearch);
    else{
        cout << "Undefined Optimization Algorithm. \n";
        exit(1);}


        pmodel->setMaxIters(iters);
        pmodel->optimise_al();

        if (disp_level > 0)
            cout << "Total Execution Time : "<<float( clock () - begin_time ) /  CLOCKS_PER_SEC<< "s \n";
        WriteMToFile(*pmodel, outFName);
        exit(0);
}


void CNN::helpHeader()
{

}

void CNN::helpInfo()
{

}

void CNN::readDataFile(mat &X, int *data_size, const string fileName)
{
    ifstream in(fileName.c_str());
    if (!in.is_open())
    {cout<< "Reading Aborted. \n";
        exit(1);
    }
    string line;
    string token;
    int featNum;
    bool featureRead = false;
    int numData = data_size[0];
    int maxFeat = data_size[1];
    ifstream inToo(fileName.c_str());
    int pointNo = 0;
    while (getline(inToo, line))
    {
        if (line[line.size() - 1] == '\r')
            line.erase(line.size() - 1);
        featureRead = false;
        if (line[0] == '#')
            continue;
        else
        {
            int pos = 0;
            featNum = 0;
            while (pos < line.size())
            {
                token.erase();
                while (pos < line.size() && (line[pos] != '\t' && line[pos] != ',' && line[pos] != ' '))
                {
                    token += line[pos];
                    pos++;
                }
                pos++;
                if (token.size() > 0)
                {

                    // TODO Check that : is in the string.
                    string featStr = token.substr(0, pos);

                    //string featValStr = token.substr(ind + 1, token.size() - ind);
                    //int featNum = atoi(featStr.c_str());
                    if (featNum > maxFeat || pointNo < 0 || pointNo >= numData)
                    {

                        std::cout << "Error";
                    }

                    double featVal = atof(featStr.c_str());
                    //data_push[featNum] = featVal;
                    //     X[pointNo*maxFeat + featNum] = featVal;
                    X(pointNo, featNum) = featVal;
                    //X.setVal(featVal, pointNo, featNum);
                    //printf(" featNum %d, X %f\n", featNum, X.getVal(pointNo, featNum));
                    featNum ++;

                }
            }

        }
        //X.push_back(data_push);
        pointNo++;
    }

}


int* CNN::readDataSize(const string fileName)
{
    int *data_size = new int[2];
    ifstream in(fileName.c_str());
    if (!in.is_open())
        cout<< "File is open \n";
    string line;
    string token;
    int featNum;
    bool featureRead = false;
    int numData = 0;
    int maxFeat = 0;
    while (getline(in, line))
    {
        featureRead = false;
        if (line[line.size() - 1] == '\r')
        {      line.erase(line.size() - 1);
            //printf("line Size%d, r%d\n", line.size(), '\r');
        }
        if (line[0] == '#')
            continue;
        numData++;
        int pos = 0;
        featNum = 0;
        while (pos < line.size())
        {
            //     printf("pos %d\n", pos);
            token.erase();
            while ( pos < line.size() && ( line[pos] != '\t' && line[pos] != ',' && line[pos] != ' ' ) )
            {
                token += line[pos];

                pos++;
            }
            pos++;

            if (token.size() > 0)
            {
                //  fputs(featureRead ? "true" : "false", stdout);
                // deal with token.


                int ind = token.find('\t');
                if (ind == -1)
                    int ind = token.find(',');
                if (ind == -1)
                    int ind = token.find(' ');


                string featStr = token.substr(0, pos);

                //          int featNum = atoi(featStr.c_str());
                featNum ++;
                if (featNum > maxFeat)
                    maxFeat = featNum;

            }
        }

    }
    data_size[1] = maxFeat;
    data_size[0] = numData;
    if (disp_level > 0){
        cout << "Data number of features: " << maxFeat << endl;
        cout << "Number of data: " << numData << endl;}

        in.close();
        return (data_size);
}
