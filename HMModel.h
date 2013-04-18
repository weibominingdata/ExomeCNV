#pragma once

#include <map>
#include <vector>
#include <string>
#include <iostream>

struct ReadDepthDataItem
{
	std::string chr;
	int startPos;
	int endPos;
	int windowSize;
	double count;        // in original scale
	double hFuntionGC;   // in original scale
	double logMap;       // in log scale
	double mprop;
	double gprop;
	double expected;
	double a1count;
	double a2count;
	double a1prop;
	double a2prop;
};

struct rdItem // storing info for read depth calculation
{
	rdItem():counts(0),num(0),avg(0){}
	void calAvg(){avg = counts/num;}
	double counts;
	int num;
	double avg;
};

struct rdvector
{
	// 0.1 as an interval
	std::vector<rdItem> v;
};

class rdmatrix
{
	// 0.1 as an interval
	// first map, second gc
public:
	rdmatrix()
	{
		m.resize(10);
		for(int i = 0; i < 10; ++i)
			m[i].v.resize(10);
	}

public:
	void setValue(double mprop, double gprop, double c)
	{
		int minx = getIndex(mprop);
		int ginx = getIndex(gprop);
		m[minx].v[ginx].counts += c;
		m[minx].v[ginx].num++;
	}
	void calAvg()
	{
		for(int i = 0; i < 10; ++i)
		{
			for(int j = 0; j < 10; ++j)
			{
				m[i].v[j].calAvg();
			}
		}
	}
	double getValue(double mprop, double gprop)
	{
		int minx = getIndex(mprop);
		int ginx = getIndex(gprop);
		return m[minx].v[ginx].avg;
	}

private:
	std::vector<rdvector> m;

private:
	int getIndex(double v)
	{
		//0<=v<=1.0
		if (v<0 || v > 1.0)
		{
			std::cerr << "error" << std::endl;
			return -1;
		}
		int i = 0;
		for(i=1; i < 10; ++i)
		{
			if (i/10.0 > v)
				break;
		}
		return i-1;
	}


};

struct ReadDepthData
{
	char * fileName;
	std::vector<ReadDepthDataItem> data;
	double largestReadCount;
	double medianReadCount;
	double medianlogmap;
	double medianhgc;
	double meanDis;
	rdmatrix rd;
	void loadData(char * f, bool allele=false);
	std::string chr;
};

class HMModel
{
public:
	HMModel(void);
	~HMModel(void);
	HMModel(const HMModel & m);


	// read depth data
	void loadReadDepthData(char * filename);
	void setReadDepthVariable();

	// read depth variable
	ReadDepthData inferData;


	void calculateMuAndPhiAllStatesCombined(bool init=false);
	void calculateMuAndPhiWithAutoRegressionAllStatesCombined();

	void startFromCoefficient();
    void setStates(int n){this->nSTATES=n;}

	double **mu;
	double *phi;

	double mixtureProportion;     // proportion of mixture component
	double mixtureProportionNormal; // proportion of mixture compoent of state 2


public:
	void writeResult(void);
	void inferAndEstimation(int rounds,bool wkv=true);

	bool USINGMAPPABILITY;
	bool USINGAUTOREGRESSION;
	bool USINGMIXTURECOMPONENT;
	bool REESTIMATETRANSITION;
	bool REESTIMATEINIT;
	bool HUMAN;
	bool POSTPROCESSING;
	bool GIVENSTATES;
	bool ALLELESPECIFICDATA;

public:
	void computAlpha(void);
	void computBeta(void);
	void computLikelihood(void);
	void computGamma(void);
	void reEstimation(bool transitionReestimate=true, bool initReestimation=true);//Baum-Welsh
	void findBestPath(bool viterbi=true);//Viterbi Algorithm
	void printVariable(void);
	void doOneRoundInference();
	void post_processing();
	void remove_smallCNV();
	void merge_CNV();
	int mostFrequentState(int l, int r);
	void getSegInfo(int l, int r, double &avemprop, double &avegprop, double &rd, double &expect);
	bool calculateAllelicConfiguration(int pos, int state, std::string &allelic_con, double &score);
	bool getAllelicConfiguration(int l_bound, int r_bound, int state, std::string &allelic_con, double &score, double &a1count, double &a2count, std::vector<std::string> &allelic_cons, std::vector<double> &scores);
	void getAllelicConfigTable(std::vector<int> &acntypes, std::map<int, std::vector<std::string> > &c, std::map<int, std::vector<double> > &v);


public:
	int nSTATES;
	int nLength;
	int normalStates; // index of normal states
	double **pAlpha; // now store the log value
	double **pBeta;  // now store the log value
	double **pGamma; // now store the log posterial value
	double *pPi;   // the initial probability
	double **pTranTbl; 
	double **pEmissTbl;
	double ***pTran;
	std::vector<double> cLikelihood;   // log likelihood
	int * inferenceResults;
	double largestReadCount;      // largest read count
	double medianReadCount;    // median of read count
	double medianLogmap;
	double medianHgc;
	double normalSelfTran;
	double otherSelfTran;
	std::vector<double> median;


	int nITRATION;


	void fillEmissTblItem(int site,int state);
	void fillEmissionTbl(void);
	void fillTranContinous(void); 
	void fillTranDiscrete(void);
	void setTranInitValue(double **pTranTbl);
	double lamda;   // state duration parameter. assume now it is not state specific


public:

	


//  SNP arrary data, not use for this time   
	std::map<int, double> cLRR;
	std::map<int, double> cBAF;
	std::map<int, double> cPFB;
	std::vector<std::string> cName;
	std::vector<int> cPos;

	char * snpdataname;
	char * infodataname;
	std::string path;
	//void setVariables(void);
	double ***pKexi;
	double **pGa;

public:
	std::string chrSymbol;
	void setFileName(char * sn, char * in);
public:
	double *pircn;
	double *murcn;
	double *sdrcn;
	double *pir;
	double *mur;
	double *sdr;
	double *pib;
	double **mub;
	double **sdb;

public:
	// for debug
	void writeKeyValue(int index);
	// for debug
	void writeKeyValueTable(char *filename, double **);

};
