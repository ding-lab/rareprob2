#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <algorithm>

#include "hmm.h"
#include "nrutil.h"
#include "global.h"

using namespace std;


int main(int argc,char *argv[])
{


	/**---------get the parameters--------**/
	char *filepath = 0;
    const char *omega = "0.5";

	int num_site = 100;
    string tmpfile_path = "tmp/";
	string RFile ="";
	string RFileName="";
	
    /**
    filetype = 0; error;
    filetype = 1; data file hasn't R region or reference X information
    filetype = 2; data file has only reference X information
    filetype = 3; data file has both R and X information
    */
    int filetype = 0;

	if( argc == 2 )
	{
        if (argv[1][0] == '-')
        {
            arg_err();
            exit(1);
        }
        filepath = argv[1];
        filetype = 3;
	}
	else if(argc == 3)
	{
        if (argv[1][0] != '-')
        {
            arg_err();
            exit(1);
        }
        if (argv[1][1] != 'r' && argv[1][1] != 'x')
        {
            arg_err();
            exit(1);
        }
        if ( argv[1][1] == 'x')
        {

            filetype = 2;
        }
        if ( argv[1][1] == 'r' )
        {
            filetype = 1;
        }
        filepath = argv[2];
	}
	else if(argc == 4)
	{
        if ( argv[1][0] != '-' || argv[1][1] != 'o')
        {
            arg_err();
            exit(1);
        }

        omega = argv[2];

        double o = atof(omega);

        if( o <= 0 || o >= 1)
        {
            cerr << "omega " << omega << " is an illegal value! Please set omega in (0,1)." << endl;
            exit(1);
        }
        filetype = 3;
        filepath = argv[3];
	}
	else if(argc == 5)
	{
        if (argv[1][0] != '-')
        {
            arg_err();
            exit(1);
        }

        if ( argv[1][1] == 'x')
        {
            filetype = 2;
        }else if ( argv[1][1] == 'r' )
        {
            filetype = 1;
        }
        else
        {
            arg_err();
            exit(1);
        }

        if ( argv[2][0] != '-' || argv[2][1] != 'o')
        {
            arg_err();
            exit(1);
        }

        omega = argv[3];

        double o = atof(omega);

        if( o <= 0 || o >= 1)
        {
            cerr << "omega " << omega << " is an illegal value! Please set omega in (0,1)." << endl;
            exit(1);
        }

        filepath = argv[4];
	}
	else
	{

        arg_err();
        exit(1);
	}

    genotype g(string(filepath),filetype);

	num_site=g.getNumSite();
	ofstream os;

	/**save Maf vector And Count vector to file**/
	string maf_seqpath = tmpfile_path + g.getFilename() + "_maf.seq";

	os.open(maf_seqpath.c_str());
	if(!os)
	{
	    cerr << "creat file error1" << endl;
	    exit(1);
	}

	vector<double> mafCase = g.getMafCase();
	vector<double> mafControl = g.getMafControl();
	vector<double> maf = g.getMaf();
	vector<int> countCase = g.getCountCase();
	vector<int> countControl = g.getCountControl();
	for(vector<double>::const_iterator iter = mafCase.begin();iter!=mafCase.end();iter++)
	{
		os<<*iter<<" ";
	}
	os<<endl;

	for(vector<double>::const_iterator iter = mafControl.begin();iter!=mafControl.end();iter++)
	{
		os<<*iter<<" ";
	}
	os<<endl;

	for(vector<double>::const_iterator iter = maf.begin();iter!=maf.end();iter++)
	{
		os<<*iter<<" ";
	}
	os<<endl;

	for(vector<int>::const_iterator iter = countCase.begin();iter!=countCase.end();iter++)
	{
		os<<*iter<<" ";
	}
	os<<endl;

	for(vector<int>::const_iterator iter = countControl.begin();iter!=countControl.end();iter++)
	{
		os<<*iter<<" ";
	}
	os<<endl;

	os.clear();
	os.close();


	/**-----------save parameters to file-----------**/
	/*
	filename m n p_value mafseq_path
	m:means the number of site
	n:means the number of individuals
	na:means the number of cases
	nc:means the number of controls
	p_value:means the value when we used in intialing X, we assign p_value = 0.001
	mafseq_path:means the file path we save the maf and count vectors
	*/
	string parmpath = tmpfile_path + g.getFilename() + "_parms";
	os.open(parmpath.c_str());
    if(!os)
	{
	    cerr << "open file error1" << endl;
	    exit(1);
	}
	os <<"filename m n na nc p_value mafseq_path"<<endl;
	os  <<g.getFilename()<<" "<<g.getNumSite()<<" "
        <<g.getNumCaseControl()<<" "<<g.getNumCase()<<" "
        <<g.getNumControl()<<" "<<P_VALUE<<" "<<maf_seqpath<<endl;
	os.clear();
	os.close();

	string referenceXpath = tmpfile_path +g.getFilename() + "_referenceX.seq";
	string referencex=g.getReferenceX();
	string referencer=g.getReferenceR();
	os.open(referenceXpath.c_str());
	if(!os)
	{
	    cerr << "open file error2" << endl;
	    exit(1);
	}
	for(int i=0;i<num_site;i++)
	{
		os << referencex[i];
		os << ' ';
	}

	os << endl;
	os.clear();
	os.close();
	/**---execute r script to initial Xseq----**/

    string cmd = "Rscript -e \"tmpfile_path='"+ tmpfile_path +"'\" -e \"omega="+string(omega)+"\" -e \"path='"+ g.getFilename() +"'\" -e \"source('./r/initialX.r')\" >>"+ tmpfile_path + g.getFilename() + "_out_temp 2>>"+ tmpfile_path + g.getFilename() +"_out_error";
    system(cmd.c_str());

	/**-----initial R using HMM-----**/
//	string initX;

	/*read X seq from file*/
/*	ifstream ifs;
    ifs.open((tmpfile_path+g.getFilename()+"_initX.seq").c_str());
    if(!ifs)
	{
	    cerr << "open file error3" << endl;
	    exit(1);
	}

    char c;

    while (ifs.get(c))
    {
		if(c=='0'||c=='1')
			initX.push_back(c+1);
    }
	ifs.clear();
	ifs.close();

    string initR;
    initialR(initR,initX,num_site);

	//save R to file
	os.open((tmpfile_path+g.getFilename()+"_initR.seq").c_str());
	if(!os)
	{
	    cerr << "open file error4" << endl;
	    exit(1);
	}
	for(unsigned i=0;i<initR.length();i++)
	{
		os<<initR[i]<<" ";
	}
	os<<endl;
	os.clear();
	os.close();

*/
	/**----executing parameters estimating rscript----**/

    cmd = "Rscript -e \"tmpfile_path='" + tmpfile_path +"'\" -e \"path='"+ g.getFilename() +"'\" -e \"source('./r/estimate.r')\" >>"+  tmpfile_path + g.getFilename() + "_out_temp 2>>"+  tmpfile_path + g.getFilename() +"_out_error";
	system(cmd.c_str());

	/**--------statistic of result------------**/

	int countx=0;
	int count_result=0;
	int count_right=0;
	int countr=0;
	int countr_result=0;
	int countr_right=0;
	double p=0;
	double s=0;

    vector<int> selected_causal;

	ifstream ifsp;//file P_VALUE
	ifstream ifss;//file STATISTIC
	ifstream ifsx_result;//file _resultX.seq
	ifstream ifsr;


	ifsp.open((tmpfile_path+g.getFilename()+"_P_VALUE").c_str());
	ifss.open((tmpfile_path+g.getFilename()+"_STATISTIC").c_str());
	if(!ifsp)
	{
	    cerr << "open file error5" << endl;
	    exit(1);
	}
	if(!ifss)
	{
	    cerr << "open file error6" << endl;
	    exit(1);
	}
	ifsp >> p;
	ifss >> s;

	char cxresult;
	char cr;
	string resultr;
	string resultx;

	ifsx_result.open((tmpfile_path+g.getFilename()+"_resultX.seq").c_str());
	ifsr.open((tmpfile_path+g.getFilename()+"_resultR.seq").c_str());

	if(!ifsx_result)
	{
	    cerr << "open file error7" << endl;
	    exit(1);
	}
	if(!ifsr)
	{
	    cerr << "open file error8" << endl;
	    exit(1);
	}

	for(int j=0;j<num_site;j++)
	{
		ifsx_result >> cxresult;
		ifsr>>cr;
		resultr.push_back(cr);

		/**---X--*/
		if(referencex[j]=='1'&&cxresult=='1')
		{
			countx++;
			count_result++;
			count_right++;
			selected_causal.push_back(j);
		}
		if(referencex[j]=='1'&&cxresult!='1')
		{
			countx++;
		}
		if(referencex[j]!='1'&&cxresult=='1')
		{
			count_result++;
			selected_causal.push_back(j+1);
		}

		/**---R---*/
		if(referencer[j]=='1'&&cr=='1')
		{
			countr++;
			countr_result++;
			countr_right++;
		}
		if(referencer[j]=='1'&&cr!='1')
		{
			countr++;
		}
		if(referencer[j]!='1'&&cr=='1')
		{
			countr_result++;
		}
		resultx.push_back(cxresult);
	}
	ifsp >> p;
	ifss >> s;
	ifsx_result.clear();
	ifsx_result.close();
	ifsp.clear();
	ifsp.close();
	ifss.clear();
	ifss.close();
	ifsr.clear();
	ifsr.close();
	cout << "filename: " <<g.getFilename()<<endl;
	cout << "p_value: " << p<<endl;
	cout << "statistic: " << s<<endl;

    cout << "selected causal rare varients sites:" << endl;
    for(vector<int>::const_iterator iter = selected_causal.begin();iter!=selected_causal.end();++iter)
    {
        cout << *iter << " ";
    }
    cout << endl;

    if(filetype == 1 || filetype == 2)
    {
        cout << "count of x/selected_X/right_x: " << countx << " "
								<< count_result << " "
								<< count_right << endl;
    }

    if(filetype == 1)
    {
       	cout << "count of R/selected_R/right_R: " << countr << " "
								<< countr_result << " "
								<< countr_right << endl;
    }

	/**---clear the temp files-----**/
//	RFile = (tmpfile_path+g.getFilename()+"_resultR.seq");
//	RFileName = ("./RFile/"+g.getFilename()+"_resultR.seq");
	string cpcmd = "cp ./tmp/"+g.getFilename()+"_resultR.seq ./RFile/"+g.getFilename();	
	system(cpcmd.c_str());
    rmtmpfile(g.getFilename());

    return 0;
}
