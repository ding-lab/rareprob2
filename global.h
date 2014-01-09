#ifndef GLOBAL_H_INCLUDED
#define GLOBAL_H_INCLUDED
#include <fstream>
#include <string>
#include <vector>


#define P_VALUE 0.001


using namespace std;
/*

*/

class genotype
{
private:
    string filename;
    int numCausal;
    int numSite;
    int numNoncausal;


    int numCase;
    int numCaseControl;
    int numControl;



    vector<int> countCase;
    vector<int> countControl;
    vector<int> countCaseControl;

    vector<double> mafCase;
    vector<double> mafControl;
    vector<double> maf;

    string referenceX;
	string referenceR;
    string ex_referenceX;


    int filetype;


public:
    genotype(string fname,int filetype);
	//genotype(ifstream &ifs);
	void print();

	string getFilename(){return filename;}
	vector<double> getMafCase(){return mafCase;}
	vector<double> getMafControl(){return mafControl;}
	vector<double> getMaf(){return maf;}
	vector<int> getCountCase(){return countCase;}
	vector<int> getCountControl(){return countControl;}
	string getReferenceX(){return referenceX;}
	string getReferenceR(){return referenceR;}

	//int getNumCausal(){return numCausal;}
	int getNumSite(){return numSite;}
	int getNumCase(){return numCase;}
	int getNumControl(){return numControl;}
	int getNumCaseControl(){return numCaseControl;}



};


void initialR(string &R,const string &X,int num_site);
void arg_err();

void rmtmpfile(const string filename);



#endif // GLOBAL_H_INCLUDED
