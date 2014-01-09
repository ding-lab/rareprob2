#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <string.h>
#include <algorithm>

#include "global.h"
#include "hmm.h"
#include "nrutil.h"

using namespace std;

//check X or R vector
bool isNot01(const char c)
{
    if(c != '0' && c != '1')
        return true;
    return false;
}

/**constructor:read the file stream and initial the class**/
genotype::genotype(string fname,int filetype)

{

	fstream ifs(fname.c_str());

	if(!ifs)
	{
		cerr << "can not open file " << fname << endl;
		exit(1);
	}

    filename = fname;
    int flag = filename.rfind("/");

	if(flag != string::npos)
	{
		filename =  filename.substr(flag+1);
	}

    /******count the matrix******/
    string str;
    int lineCount=1;
    string::iterator pos;

    switch(filetype)
    {
        case 1:
            ifs >> referenceR;
            ifs >> referenceX;

            pos = find_if(referenceR.begin(),referenceR.end(),isNot01);
            if(pos != referenceR.end())
            {
                //cout << *pos << endl;
                cerr << "genotype file has wrong format on line 1!" << endl;
                exit(1);
            }
            pos = find_if(referenceX.begin(),referenceX.end(),isNot01);
            if(pos != referenceX.end())
            {

                cerr << "genotype file has wrong format on line 2!" << endl;
                exit(1);
            }

            if( referenceR.length()!=referenceX.length())
            {
                cerr << "genotype file has wrong format on line 1 or 2!" << endl;
                exit(1);
            }

            numSite = referenceX.length();
            lineCount = 3;
            break;
        case 2:
            ifs >> referenceX;

            pos = find_if(referenceX.begin(),referenceX.end(),isNot01);
            if(pos != referenceX.end())
            {
                cerr << "genotype file has wrong format on line 1!" << endl;
                exit(1);
            }
            numSite = referenceX.length();

            lineCount = 2;
            break;
        case 3:
            ifs >> str;
            numSite = str.length() - 1;
            ifs.seekp(fstream::beg);

            lineCount = 1;

            break;
        default:
            cerr << "wrong file format!" << endl;
            exit(1);
            break;
    }

    /**initial all vectors*/
	countCase.assign(numSite,0);
	countControl.assign(numSite,0);
	countCaseControl.assign(numSite,0);
	mafCase.assign(numSite,0);
	mafControl.assign(numSite,0);
	maf.assign(numSite,0);


	int ind=0;
    int indC=0;
    int indA=0;

	/*----------------caculate countCausal,maf--------------*/
    while (ifs >> str)
    {

		if(str.length() != (unsigned)(numSite+1) && str.length() != (unsigned)numSite)
        {
            cerr << "the genotype file has wrong format on line " << lineCount + 1 << endl;
        	exit(1);
        }

        try
        {
            if(str[numSite]=='A')
            {
                indA++;
                for(int i = 0; i<numSite ;i++)
                {
                    if(str[i]=='1')
                    {
                        countCase[i]++;
                        countCaseControl[i]++;
                    }
                    else if(str[i]!='0')
                    {
                        cerr << "the genotype file has wrong format on line " << lineCount + 1 << endl;
                        exit(1);
                    }
                }
            }
            else if(str[numSite]=='C')
            {
                indC++;
                for(int i = 0; i<numSite ;i++)
                {
                    if(str[i]=='1')
                    {
                        countControl[i]++;
                        countCaseControl[i]++;
                    }
                    else if(str[i]!='0')
                    {
                        cerr << "the genotype file has wrong format on line " << lineCount + 1 << endl;
                        exit(1);
                    }
                }
            }
            else
            {
                cerr << "the genotype file has wrong format on line " << lineCount + 1 << endl;
                exit(1);
            }
            lineCount++;
        }
        catch(...)
        {
            cerr << "the genotype file has wrong format on line " << lineCount + 1 << endl;
            exit(1);
        }

	}
	ifs.close();
    ind = indA + indC;
    this->numCase = indA;
    this->numControl = indC;
    this->numCaseControl = ind;

	/*-----------caculate maf------------*/
	for(size_t i = 0;i<(unsigned)numSite;i++)
	{
		mafCase[i] = static_cast<double>(countCase[i])/static_cast<double>(numCase);
		mafControl[i] =static_cast<double>(countControl[i])/static_cast<double>(numControl);
		maf[i] = static_cast<double>(countCaseControl[i])/static_cast<double>(numCaseControl);
	}

}

/*--------output the data of class---------------*/
void genotype::print()
{
	cout<<"filename:"<<filename<<endl;
	cout<<"numCausal:"<<numCausal<<endl;
	cout<<"numNoncausal:"<<numNoncausal<<endl;
	cout<<"numSite:"<<numSite<<endl;
	cout<<"numCase:"<<numCase<<endl;
	cout<<"numControl:"<<numControl<<endl;
	cout<<"numCaseControl:"<<numCaseControl<<endl;

	cout<< "referenceX:" << this->referenceX << endl;
	cout<<"countCase:"<<endl;
	for(size_t i = 0; i<(unsigned)numSite; i++)
	{
		cout<<countCase[i]<<" ";
	}
	cout<<endl;
	cout<<"countControl:"<<endl;
	for(size_t i = 0; i<(unsigned)numSite; i++)
	{
		cout<<countControl[i]<<" ";
	}
	cout<<endl;
	cout<<"countCaseControl:"<<endl;
	for(size_t i = 0; i<(unsigned)numSite; i++)
	{
		cout<<countCaseControl[i]<<" ";
	}
	cout<<endl;
	cout<<"mafCase:"<<endl;
	for(size_t i = 0; i<(unsigned)numSite; i++)
	{
		cout<<mafCase[i]<<" ";
	}
	cout<<endl;
	cout<<"mafControl:"<<endl;
	for(size_t i = 0; i<(unsigned)numSite; i++)
	{
		cout<<mafControl[i]<<" ";
	}
	cout<<endl;
	cout<<"maf:"<<endl;
	for(size_t i = 0; i<(unsigned)numSite; i++)
	{
		cout<<maf[i]<<" ";
	}
	cout<<endl;
}




void initialR(string &R,const string &X,int num_site)
{

    char *chR=new char[X.length()+1];
    char *chX=new char[X.length()+1];
    strcpy(chX,X.c_str());

    int	*O=new int[num_site+1];
    for(int i=0;i<num_site;i++)
    {
        if(X[i]=='1')
            O[i+1]=1;
        else
            O[i+1]=2;

    }


	HMM  	hmm;
	int m = num_site;

	int	N=2;
	int	M=2;
	double 	**alpha;
	double	**beta;
	double	**gamma;

	int	seed=0; /* seed for random number generator */

	int	niter=0;
	double	logprobinit=0, logprobfinal=0;


    seed = hmmgetseed();
    InitHMM(&hmm, N, M, seed);

    alpha = dmatrix(1, m, 1, hmm.N);
	beta = dmatrix(1, m, 1, hmm.N);
    gamma = dmatrix(1, m, 1, hmm.N);

    //initial the hmm
    hmm.A[1][1] = 0.8;
    hmm.A[1][2] = 0.2;
    hmm.A[2][1] = 0.2;
    hmm.A[2][2] = 0.8;
    hmm.B[1][1] = 0.95;
    hmm.B[1][2] = 0.05;
    hmm.B[2][1] = 0.5;
    hmm.B[2][2] = 0.5;
    hmm.M = M;
    hmm.N = N;
    hmm.pi[1] = 0.5;
    hmm.pi[2] = 0.5;

    BaumWelch(&hmm, m, O, alpha, beta, gamma, &niter,
		&logprobinit, &logprobfinal);

	/* output the result hmm */
	/*
	FILE *tempFp=fopen("result.hmm","w");
	PrintHMM(tempFp,&hmm);
    fclose(tempFp);
    */

	int	**psi;
	double 	logproba;
	int *q = new int[num_site+1];

    double **delta = dmatrix(1, m+1, 1, hmm.N+1);
	psi = imatrix(1, m+1, 1, hmm.N+1);

	/* note: ViterbiLog() returns back with log(A[i][j]) instead
	** of leaving the A matrix alone. If you need the original A,
	** you can make a copy of hmm by calling CopyHMM */

	ViterbiLog(&hmm, m, O, delta, psi, q, &logproba);

    //R=q;
    for(int i = 1;i<=num_site;i++)
    {
		if(q[i]==1)
            R.push_back('0');
        else
            R.push_back('1');
    }

	free_imatrix(psi, 1, m, 1, hmm.N);
	free_dmatrix(delta, 1, m, 1, hmm.N);

	free_dmatrix(alpha, 1, m, 1, hmm.N);
	free_dmatrix(beta, 1, m, 1, hmm.N);
	free_dmatrix(gamma, 1, m, 1, hmm.N);
	FreeHMM(&hmm);

    delete []q;
    delete []O;
    delete []chX;
    delete []chR;
}


void arg_err()
{
    cerr << "-----------------------------------------------------------------" << endl;
    cerr << "Usage of RareProb:" << endl;
    cerr << "-----------------------------------------------------------------" << endl;
    cerr << "rareprob [-r|x] [-o omega-value] filename" << endl;
    cerr << "1. -r the genotype file has both region information and reference \
causal rare varients information" << endl;
    cerr << "2. -x the genotype file has only reference causal rare varients information" << endl;
    cerr << "3. default setting without -x or -r means the genotype file has \
neither reference causal rare varients nor region information" << endl;
    cerr << "4. -o omega-value: you can set the threshhold of omega-value through \
this option. And the legal value of omega-value is in (0,1). \
If you don't know how to set this value, you can ignore this option and use default value 0.5. " << endl;
    cerr << "5. filename: the genotype file. This argument is always been given." << endl;

    cerr << "e.g." << endl;
    cerr << "1. rareprob -r -o 0.7 genotypefile" << endl;
    cerr << "   It handles the first type genotype file with omega threshold 0.7." << endl;
    cerr << "1. rareprob -x -o 0.7 genotypefile" << endl;
    cerr << "   It handles the second type genotype file with omega threshold 0.7." << endl;
    cerr << "3. rareprob -r genotypefile" << endl;
    cerr << "   It handles the first type genotype file with default omega threshold 0.5." << endl;
    cerr << "4. rareprob -x genotypefile" << endl;
    cerr << "   It handles the second type genotype file with default omega threshold 0.5." << endl;
    cerr << "5. rareprob -o 0.7 genotypefile" << endl;
    cerr << "   It handles the third type genotype file with omega threshold 0.7." << endl;
    cerr << "6. rareprob genotypefile" << endl;
    cerr << "   It handles the third type genotype file with default omega threshold 0.5." << endl;

}


void rmtmpfile(const string filename)
{
    string filepath = "./tmp/"+filename ;

    remove((filepath+"_initR.seq").c_str());
    remove((filepath+"_initX.seq").c_str());
    remove((filepath+"_maf.seq").c_str());
    remove((filepath+"_newr.seq").c_str());
    remove((filepath+"_newx.seq").c_str());
    remove((filepath+"_Omega.seq").c_str());
    remove((filepath+"_parms").c_str());
    remove((filepath+"_P_VALUE").c_str());
    remove((filepath+"_referenceX.seq").c_str());
/*    remove((filepath+"_resultR.seq").c_str());*/
    remove((filepath+"_resultX.seq").c_str());
    remove((filepath+"_STATISTIC").c_str());
    remove((filepath+"_Wss.seq").c_str());
    remove((filepath+"_X.seq").c_str());
    remove((filepath+"_Zs.seq").c_str());
    remove((filepath+"_out_error").c_str());
    remove((filepath+"_out_temp").c_str());
}
