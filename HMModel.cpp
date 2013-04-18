#include "HMModel.h"

// starting of my header files
#include <iostream>
#include <fstream>
#include <cmath>
#include <assert.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <algorithm>
#include <map>

// math tools, NB regression
#include "MathTools.h"

using namespace std;
// ending of my header files


void ReadDepthData::loadData(char *f, bool allele)
{
	int lastStartPos = 0;
	largestReadCount = 0;
	vector<double> readcounts;
	vector<double> logmap;
	vector<double> hgc;
	string str;
	ifstream fin(f);
	int lastEnd=-1;
	int totalDis=0;
	if (!allele){
		fin >> str >> str >> str >> str >> str >> str >> str
			>> str >> str >> str >> str >> str >> str >> str;
		string temp;

		while(!fin.eof())
		{
			ReadDepthDataItem item;
			fin >> item.chr;
			if (fin.eof())
				break;
			fin >> item.startPos >> item.endPos >> item.windowSize;
			fin >> temp >> item.mprop >> item.gprop;
			fin >> item.count;
			fin >> temp >> temp >> temp;
			fin	>> item.hFuntionGC >> item.logMap >> item.expected;
			if (item.count > largestReadCount)
				largestReadCount = item.count;
			readcounts.push_back(item.count);
			logmap.push_back(item.logMap);
			hgc.push_back(item.hFuntionGC);
			rd.setValue(item.mprop, item.gprop, item.count);
			data.push_back(item);
			if (lastEnd!=-1){
				totalDis+=item.startPos-lastEnd;
			}
			lastEnd=item.endPos;
		}
	} else {
		for(int i=0; i<18; ++i){
			fin >> str;
		}
		string temp;
		while(!fin.eof())
		{
			ReadDepthDataItem item;
			fin >> item.chr;
			if (fin.eof())
				break;
			fin >> item.startPos >> item.endPos >> item.windowSize;
			fin >> temp >> item.mprop >> item.gprop;
			fin >> item.count;
			fin >> item.a1count >> item.a2count;
			fin >> temp >> temp >> temp;
			fin	>> item.hFuntionGC >> item.logMap >> item.expected;
			fin >> item.a1prop >> item.a2prop;
			if (item.count > largestReadCount)
				largestReadCount = item.count;
			readcounts.push_back(item.count);
			logmap.push_back(item.logMap);
			hgc.push_back(item.hFuntionGC);
			rd.setValue(item.mprop, item.gprop, item.count);
			data.push_back(item);
			if (lastEnd!=-1){
				totalDis+=item.startPos-lastEnd;
			}
			lastEnd=item.endPos;
		}
	}

	chr = data[0].chr;
	rd.calAvg();
	
	cout << "Finish Reading Files" << endl;

	sort(readcounts.begin(),readcounts.end());
	medianReadCount = readcounts[int(readcounts.size()/2)];

	sort(logmap.begin(),logmap.end());
	medianlogmap = logmap[int(logmap.size()/2)];

	sort(hgc.begin(),hgc.end());
	medianhgc = hgc[int(hgc.size()/2)];

	meanDis = totalDis/int(data.size()-1);

	//int lastStartPos = 0;
	//string str;
	//ifstream fin(f);
	//fin >> str >> str >> str >> str >> str >> str >> str;
	//string temp;
	//
	//while(!fin.eof())
	//{
	//	ReadDepthDataItem item;
	//	fin >> item.chr;
	//	if (fin.eof())
	//		break;
	//	fin >> item.startPos >> item.endPos >> item.state >> item.logMap >> item.hFuntionGC
	//		>> item.count;
	//	item.windowSize = 500;
	//	data.push_back(item);
	//}

	//chr = data[0].chr;
}



////////////////////////////////////////////////////
// an auxiliary data structure for data processing 
// store LRR and BAF value
struct LRRBAF
{
	LRRBAF() {LRR=BAF=-10.0;}
	double LRR;
	double BAF;
};


HMModel::HMModel(void)
: nSTATES(0)
, nLength(0)
, pTranTbl(NULL)
, pEmissTbl(NULL)
, pAlpha(NULL)
, pBeta(NULL)
//, nOBSERVE(0)
, pGamma(NULL)
, pKexi(NULL)
, pGa(NULL)
, pPi(NULL)
, pircn(NULL)
, murcn(NULL)
, sdrcn(NULL)
, pir(0)
, mur(NULL)
, sdr(NULL)
, pib(NULL)
, mub(NULL)
, sdb(NULL)
, chrSymbol("1")
, lamda(0.00001)
, nITRATION(0)
, largestReadCount(0)
, medianReadCount(0)
, medianLogmap(0)
, medianHgc(0)
, mixtureProportion(0.1)
, mixtureProportionNormal(0.1)
, normalStates(2)
, normalSelfTran(0.995)
, otherSelfTran(0.95)
, USINGMAPPABILITY(true)
, USINGAUTOREGRESSION(true)
, USINGMIXTURECOMPONENT(true)
, REESTIMATETRANSITION(true)
, REESTIMATEINIT(true)
, HUMAN(true)
, POSTPROCESSING(false)
, GIVENSTATES(false)
, ALLELESPECIFICDATA(false)
{
}

void HMModel::loadReadDepthData(char * filename)
{
	setFileName(filename, filename);
	inferData.loadData(filename,ALLELESPECIFICDATA);
	chrSymbol = inferData.chr;
	largestReadCount = inferData.largestReadCount;
	medianReadCount = inferData.medianReadCount;
	setReadDepthVariable();

	startFromCoefficient();

	fillTranDiscrete();  // for read depth data, it is discrete time
	fillEmissionTbl();
}

void HMModel::startFromCoefficient()
{
	cout << "start with coefficient" << endl;
	cout << nSTATES << " " << normalStates << endl;


	double delta = 0.5;
	double coefficientforgc = 0.5;

	cout << delta << " " << coefficientforgc << " " << endl;

    double newintercept = log(medianReadCount)-log(normalStates)-medianLogmap-coefficientforgc*medianHgc;
    cout << newintercept << endl;


	for(int i = 0; i < nSTATES; ++i)
	{
		phi[i] = 1.0;
	}

	for(int i = 0; i < nLength; ++i)
	{
		for(int j = 0; j < nSTATES; ++j)
		{
			double offset = 0;
			if (j==0)
				offset = log(delta/2)+log(inferData.data[i].expected);
			else
				offset = log(j*1.0/2)+log(inferData.data[i].expected);
			mu[i][j] = exp(offset);
		}
	}

}

void HMModel::calculateMuAndPhiAllStatesCombined(bool init)
{
	// use glmNB to fill mu matrix and phi matrix 
	// for state 0, we only need to fit phi,
	// for other state, we need to get fitted mu and re-estimated phi

	// for state 0
	// load weights, load fitted value, in original scale, estimate phi
	int maxIt = 25;
	double convR = 1e-8;
	int nCovariate = 1;          // CN
	double *y = new double[nSTATES*nLength];
	double *fitted = new double[nSTATES*nLength];
	double *x = new double[nSTATES*nLength*nCovariate];
	double *prior = new double[nSTATES*nLength];
	double *weights = new double[nSTATES*nLength]; // will be changed in computation
	double *resid = new double[nSTATES*nLength]; 
	double *Xb = new double[nSTATES*nLength*nCovariate];  // used in the program
	double *offset = new double[nSTATES*nLength];         // offset, now log(state) is a offset
	double delta = 0.5; // delta for cn0

	// loaded x, y, z, fitted, weights
	// y x are not dependant with initial value
	int index = 0;
	for(int i = 0; i < nLength; ++i)
	{
		for(int j = 0; j < nSTATES; ++j)
			y[index++] = inferData.data[i].count+0.1;
	}
	index = 0;
	for(int i = 0; i < nLength; ++i)
	{
		x[index++] = log(delta/2);
		for(int j = 1; j < nSTATES; ++j)
		{
			x[index++] = log(j*0.5);
		}
	}



	// load offset
	index = 0;
	for(int i = 0; i < nLength; ++i)
	{
		for(int j = 1; j < nSTATES; ++j)
		{
			offset[index++] =log(inferData.data[i].expected+0.1);
		}
	}
	if (USINGMAPPABILITY)
	{
		index = 0;
		for(int i = 0; i < nLength; ++i)
		{
			offset[index++] += inferData.data[i].logMap/*+log(median[i]+0.01)*/;
			for(int j = 1; j < nSTATES; ++j)
			{
				offset[index++] += inferData.data[i].logMap/*+log(median[i]+0.01)*/;
			}
		}
	}


	// load weights
	index = 0;
	for(int i = 0; i < nLength; ++i)
	{
		for(int j = 0; j < nSTATES; ++j)
		{
			prior[index++] = exp(pGamma[i][j]);
		}
	}

	// update weights
	// given mu, readcount, overdispersion, have likelihood
	// calculate proportion
	if (!init && USINGMIXTURECOMPONENT)
	{
		index = 0;
		double w = 1;
		double l = 0;
		double proportion = 0;
		for(int i = 0; i < nLength; ++i)
		{
			for(int j = 0; j < nSTATES; ++j)
			{
				
				l = exp(MathTools::loglik_NB(1, phi[j], &mu[i][j], &inferData.data[i].count, &w));
				if (j==normalStates)
				{
					proportion = 
					(l*(1-mixtureProportionNormal))
					/(mixtureProportionNormal/largestReadCount+(1-mixtureProportionNormal)*l);
				}
				else
				{
					proportion = 
					(l*(1-mixtureProportion))
					/(mixtureProportion/largestReadCount+(1-mixtureProportion)*l);
				}
				
				prior[index++] *= proportion;
			}
		}
	}

	int dim[5];
	dim[0] = nLength*nSTATES;
	dim[1] = nCovariate;  
	dim[2] = maxIt; //
	dim[3] = 0; // false
	dim[4] = 1; // false, no offset

	int nIter = 0;
	int linkR = LOG;
	convR = 1e-8; 
	int rank = 0;
	double PHI = 0;  // not use
	double scale = 1.0;
	int de_resid = 0;
	int family = 0;
	double twologlik = 0;
	double scoreTestP = 0;
	int trace = 0;
	double beta = 0;



	int conv = MathTools::glmNB(dim, &nIter,y,prior, &linkR, offset, x, &convR, &rank,
		Xb, fitted, resid, weights, &PHI, &scale, &de_resid, &family,
		&twologlik, &scoreTestP, &trace, &beta);

    //if (beta>=0){
	    cout << "beta for CN content is " << beta << " " << endl;




	// save fitted
	    index = 0;
	    for(int i = 0; i < nLength; ++i)
	    {
		    for(int j = 0; j < nSTATES; ++j)
		    {
			    mu[i][j] = fitted[index++];
		    }
	    }
//    }
//    else{
//        cout << "beta becomes negative value, error will occur, just don't store the fitting value" << endl;
//    }




	// fix phi

	//delete []y; y=NULL;
	delete []x; x=NULL;
	//delete []prior; prior=NULL;
	//delete []fitted; fitted=NULL;
	delete []weights; weights=NULL;
	delete []resid; resid=NULL;
	delete []Xb; Xb=NULL;
	delete []offset;


		int cvPhi = MathTools::phi_ml(y, fitted, nLength*nSTATES, prior, maxIt,
			convR, &phi[0], 0, 0);
		//if (phi[0]>20)
		//	phi[0]=20;
		for(int j = 1; j < nSTATES; ++j)
		{
			phi[j] = phi[0];
		}


	delete []y; y=NULL;
	delete []fitted; fitted=NULL;
	delete []prior; prior=NULL;

}

void HMModel::calculateMuAndPhiWithAutoRegressionAllStatesCombined()
{
	// use glmNB to fill mu matrix and phi matrix 
	// for state 0, we only need to fit phi,
	// for other state, we need to get fitted mu and re-estimated phi

	// for state 0
	// load weights, load fitted value, in original scale, estimate phi
	int maxIt = 25;
	double convR = 1e-8;
	int nCovariate = 2;          // CN, autoRegression, //logmap
	double *y = new double[nSTATES*nLength];
	double *fitted = new double[nSTATES*nLength];
	double *x = new double[nSTATES*nLength*nCovariate];
	double *prior = new double[nSTATES*nLength];
	double *weights = new double[nSTATES*nLength]; // will be changed in computation
	double *resid = new double[nSTATES*nLength]; 
	double *Xb = new double[nSTATES*nLength*nCovariate];  // used in the program
	double *offset = new double[nSTATES*nLength];
	double delta = 0.5; // delta for cn0

	// loaded x, y, z, fitted, weights
	// y x are not dependant with initial value
	int index = 0;
	for(int i = 0; i < nLength; ++i)
	{
		for(int j = 0; j < nSTATES; ++j)
			y[index++] = inferData.data[i].count;
	}
	index = 0;
	for(int i = 0; i < nLength; ++i)
	{
		for(int j = 0; j < nSTATES; ++j)
		{
			if (i==0)
				x[index++] = 0;
			else
			{
				if (inferData.data[i-1].count == 0)
				{
					x[index++] = log(inferData.data[i-1].count+0.1)-log(mu[i-1][j]);
				}
				else
				{
					x[index++] = log(inferData.data[i-1].count)-log(mu[i-1][j]);
				}
			}   
		}
	}
	for(int i = 0; i < nLength; ++i)
	{
		for(int j = 0; j < nSTATES; ++j)
		{
			x[index++] = inferData.data[i].hFuntionGC;
		}
	}



	// load offset
	index = 0;
	for(int i = 0; i < nLength; ++i)
	{
		offset[index++] = log(delta);
		for(int j = 1; j < nSTATES; ++j)
		{
			offset[index++] = log(j*1.0);
		}
	}
	if (USINGMAPPABILITY)
	{
		index = 0;
		for(int i = 0; i < nLength; ++i)
		{
			offset[index++] += inferData.data[i].logMap/*+log(median[i]+0.01)*/;
			for(int j = 1; j < nSTATES; ++j)
			{
				offset[index++] += inferData.data[i].logMap/*+log(median[i]+0.01)*/;
			}
		}
	}



	// load weights
	index = 0;
	for(int i = 0; i < nLength; ++i)
	{
		for(int j = 0; j < nSTATES; ++j)
		{
			prior[index++] = exp(pGamma[i][j]);
		}
	}

	// update weights
	// given mu, readcount, overdispersion, have likelihood
	// calculate proportion
	if (USINGMIXTURECOMPONENT)
	{
		index = 0;
		double w = 1;
		double l = 0;
		double proportion = 0;
		for(int i = 0; i < nLength; ++i)
		{
			for(int j = 0; j < nSTATES; ++j)
			{
				
				l = exp(MathTools::loglik_NB(1, phi[j], &mu[i][j], &inferData.data[i].count, &w));
				if (j==normalStates)
				{
					proportion = 
					(l*(1-mixtureProportionNormal))
					/(mixtureProportionNormal/largestReadCount+(1-mixtureProportionNormal)*l);
				}
				else
				{
					proportion = 
					(l*(1-mixtureProportion))
					/(mixtureProportion/largestReadCount+(1-mixtureProportion)*l);
				}
				prior[index++] *= proportion;
			}
		}
	}

	int dim[5];
	dim[0] = nLength*nSTATES;
	dim[1] = nCovariate;  
	dim[2] = maxIt; //
	dim[3] = 0; // false
	dim[4] = 1; // false, no offset

	int nIter = 0;
	int linkR = LOG;
	convR = 1e-8; 
	int rank = 0;
	double PHI = 0;  // not use
	double scale = 1.0;
	int de_resid = 0;
	int family = 0;
	double twologlik = 0;
	double scoreTestP = 0;
	int trace = 0;
	double beta = 0;

	int conv = MathTools::glmNB(dim, &nIter,y,prior, &linkR, offset, x, &convR, &rank,
		Xb, fitted, resid, weights, &PHI, &scale, &de_resid, &family,
		&twologlik, &scoreTestP, &trace, &beta);


    if (beta >=0){
		cout << "beta for gc content is " << beta << " " << endl;
		//cout << "beta for mappability is " << beta << " " << endl;


		// save fitted
		// edited by Weibo, 04/27/2012
		// if fit is not good, basedline(state0) is much larger than the observation (>10), set state 0 fitted value as the observation
		index = 0;
		cout << "mu adjusted for state 0" << endl;
		int fitbadcondition = 10;
		for(int i = 0; i < nLength; ++i)
		{
			for(int j = 0; j < nSTATES; ++j)
			{
				mu[i][j] = fitted[index++];
			}
			if (mu[i][0]-inferData.data[i].count>= fitbadcondition)
			{
				mu[i][0] = inferData.data[i].count;
			}
		}
    }
    else{
        cout << "beta becomes negative value, error will occur, just don't store the fitting value" << endl;
    }



	//delete []y; y=NULL;
	delete []x; x=NULL;
	//delete []prior; prior=NULL;
	//delete []fitted; fitted=NULL;
	delete []weights; weights=NULL;
	delete []resid; resid=NULL;
	delete []Xb; Xb=NULL;
	delete []offset; offset = NULL;


		int cvPhi = MathTools::phi_ml(y, fitted, nLength*nSTATES, prior, maxIt,
			convR, &phi[0], 0, 0);
		if (phi[0]>20)
			phi[0]=20;
		for(int j = 1; j < nSTATES; ++j)
		{
			phi[j] = phi[0];
		}

	delete []y; y=NULL;
	delete []fitted; fitted=NULL;
	delete []prior; prior=NULL;
}

void HMModel::setTranInitValue(double **pTran)
{
	// load from the file, called transition_init.dat
	string tran_init="transition_init.dat";
	tran_init = path + tran_init;
    ifstream fin(tran_init.c_str());
    string comment;
    getline(fin, comment);
    while (comment.find("#ENDComments")!=0)
    {
    	getline(fin, comment);
    }
    for(int i = 0 ; i < nSTATES; ++i)
    {
    	for(int j = 0; j < nSTATES; ++j)
    	{
    		fin >> pTran[i][j];
    	}
    }
    fin.close();

}



void HMModel::setReadDepthVariable()
{
	if (!GIVENSTATES){
		// hard coding for the real human data
		if (HUMAN)
		{
			nSTATES = 7;
			normalStates = 2;
		}
		else
		{
			// mouse data modification
			nSTATES = 4;
			normalStates = 1;
		}
	}
	else{
		// nSTATES has been given
		normalStates=2;
	}

	nLength = inferData.data.size();

	mu = new double*[nLength];
	for(int i = 0; i < nLength; ++i)
	{
		mu[i] = new double[nSTATES];
	}
	// at the beginning, mu[:][i] = readcout 
	for(int i = 0; i < nLength; ++i)
	{
		for(int j = 0; j < nSTATES; ++j)
		{
			mu[i][j] = inferData.data[i].count;	
		}
	}

	phi = new double[nSTATES];

	inferenceResults = new int[nLength];

	pTranTbl = new double*[nSTATES];
	for(int i = 0; i < nSTATES; ++i)
	{
		pTranTbl[i] = new double[nSTATES];
	}

	pTran = new double **[nLength];
	for(int i = 0; i < nLength; ++i)
	{
		pTran[i] = new double *[nSTATES];
		for(int j = 0; j < nSTATES; ++j)
		{
			pTran[i][j] = new double[nSTATES];
		}
	}

	setTranInitValue(pTranTbl);
	cout << normalSelfTran << " " << otherSelfTran << endl;




	// revise pEmissTbl nSTATES * nLength
	pEmissTbl = new double*[nLength];
	for(int i = 0; i < nLength; ++i)
	{
		pEmissTbl[i] = new double[nSTATES];
	}

	// hard coding for real data
	pPi = new double[nSTATES];
	for(int i = 0; i < nSTATES; ++i)
	{
		pPi[i] = (1-normalSelfTran)/(nSTATES-1);
	}
	pPi[normalStates] = normalSelfTran;

	pAlpha = new double*[nLength];
	for(int i = 0; i < nLength; ++i)
	{
		pAlpha[i] = new double[nSTATES];
		for(int j = 0; j < nSTATES; ++j)
		{
			pAlpha[i][j] = 0;
		}
	}

	pBeta = new double*[nLength];
	for(int i = 0; i < nLength; ++i)
	{
		pBeta[i] = new double[nSTATES];
		for(int j = 0; j < nSTATES; ++j)
		{
			pBeta[i][j] = 0;
		}
	}

	pGamma = new double*[nLength];   // store the posterior probability
	for(int i = 0; i < nLength; ++i)
		pGamma[i] = new double[nSTATES];
	for(int i=0; i < nLength; ++i)
	{
		for(int j=0; j < nSTATES; ++j)
		{
			if (j==normalStates)
			{
				pGamma[i][j] = log(0.9);
			}
			else
			{
				pGamma[i][j] = log(0.1/(nSTATES-1));
			}
		}
	}

}



HMModel::~HMModel(void)
{

}

HMModel::HMModel(const HMModel & m)
{
}




void HMModel::inferAndEstimation(int rounds, bool wkv)
{
	if (wkv){
		writeKeyValue(0);
	}
	for(int i = 0; i < rounds; ++i)
	{
		doOneRoundInference();
		reEstimation(REESTIMATETRANSITION, REESTIMATEINIT);
		if (wkv){
			writeKeyValue(i+1);
		}
	}
	findBestPath(false);
	printVariable();
}

void HMModel::doOneRoundInference()
{
	nITRATION++;
	computAlpha();
	computBeta();
	computLikelihood();
	computGamma();
}

void HMModel::computAlpha(void)
{

	for(int i = 0; i < nSTATES; ++i)
	{
		pAlpha[0][i] = log(pPi[i]) + log(pEmissTbl[0][i]);
	}

	double *v = new double[nSTATES];
	for(int i = 1; i < nLength; ++i)
	{
		//if (i%100 == 0) cout << i << endl;
		for(int j = 0; j < nSTATES; ++j)
		{
			for(int k = 0; k < nSTATES; ++k)
			{
				v[k] = pAlpha[i-1][k]+log(pTran[i][k][j]);

			}
			pAlpha[i][j] = MathTools::logsumexp(v, nSTATES)+log(pEmissTbl[i][j]);
		}
	}
	delete []v;
}

void HMModel::computBeta(void)
{
	for(int i = 0;  i < nSTATES; ++i)
	{
		pBeta[nLength-1][i] = 0;
	}
	double *v = new double[nSTATES];
	for(int i = nLength-2; i >=0; --i)
	{
		//if (i%100 == 0) cout << i << endl;
		for(int j = 0; j < nSTATES; ++j)
		{
			for(int k = 0; k < nSTATES; ++k)
			{
				v[k] = pBeta[i+1][k] + log(pTran[i+1][j][k]) + log(pEmissTbl[i+1][k]);
			}
			pBeta[i][j] = MathTools::logsumexp(v, nSTATES);
		}
	}
	delete []v;
}

void HMModel::computGamma(void)
{
	// make sure comput Gamma is called after the corresponding computAlpha and computBeta
	for(int step = 0; step < nLength; ++step)
	{
		for(int i = 0; i < nSTATES; ++i)
		{
			pGamma[step][i] = pAlpha[step][i]+pBeta[step][i]-cLikelihood[nITRATION-1];
		}
	}
}

void HMModel::computLikelihood()
{
	double likelihood = 0;
	double *v = new double[nSTATES];
	for(int k = 0; k < nSTATES; ++k)
	{
		v[k] = pAlpha[nLength-1][k];
	}
	likelihood = MathTools::logsumexp(v, nSTATES);
	cLikelihood.push_back(likelihood);
	delete []v;
}


void HMModel::post_processing()
{
	// processing, including, remove <2 CNV, and merge CNV
	remove_smallCNV();
	merge_CNV();
}

void HMModel::remove_smallCNV()
{
	for(int i = 0; i < nLength; ++i)
	{
		if (inferenceResults[i] != normalStates)
		{
			int l = max(i-1,0);
			int r = min(i+1,nLength-1);
			if (inferenceResults[l] == normalStates && inferenceResults[r] == normalStates)
			{
				inferenceResults[i] = normalStates;
			}
		}
	}
}

void HMModel::merge_CNV()
{
	// when merging, only consider dup  and del
	// if consecutive cnvs appear, the boundaries will first be decided first, then the most frequent state will be assigned to this cnv
	int l = -1, r = -1;
	int cl = -1, cr=-1;
	bool incnv = false;
	for(int i = 0; i < nLength; ++i)
	{
		if (inferenceResults[i] != normalStates)
		{
			if (!incnv)
			{
				// entering a cnv
				incnv = true;
				cl = i;
			}
		}
		if (incnv && inferenceResults[min(i+1,nLength-1)] == normalStates)
		{
				incnv = false;
				cr = i;
				// find a cnv, check whether it can be merged with previous cnv
				if (l != -1)
				{
					// check whether these two can be merged
					if ((inferData.data[cr].endPos-inferData.data[cl].startPos+inferData.data[r].endPos-inferData.data[l].startPos) > 2*(inferData.data[cl].endPos-inferData.data[r].startPos)
							&& (inferData.data[cl].endPos-inferData.data[r].startPos) > min(inferData.data[cr].endPos-inferData.data[cl].startPos, inferData.data[r].endPos-inferData.data[l].startPos))
					{
						// can merge
						r = cr;
					}
					else
					{
						// can not
						// decide the state for l,r
						int s = mostFrequentState(l, r);
						for(int j = l; j <= r; ++j)
						{
							inferenceResults[j] = s;
						}
						l = cl;
						r = cr;
					}
				}
				else
				{
					l = cl;
					r = cr;
				}
		}

	}

}

int HMModel::mostFrequentState(int l, int r)
{
	vector<int> s(nSTATES);
	for(int i = l; i <=r; ++i)
	{
		s[inferenceResults[i]]++;
	}
	int maxCount = -1;
	int index = -1;
	for(int i = 0; i < nSTATES; ++i)
	{
		if (s[i] > maxCount)
		{
			maxCount = s[i];
			index = i;
		}
	}
	return index;
}



void HMModel::writeResult(void)
{
	
    //
	if (POSTPROCESSING)
	{
		cout << "post processing" << endl;
		post_processing();
	}

	int * cn = new int[nLength];
	//double max = -10;
	for(int i = 0; i < nLength; ++i)
	{
		int index = inferenceResults[i];
		cn[i] = index;
	}
	vector<string> allelic_cons(nLength,"");
	vector<double> allelic_scores(nLength,0);
	if (ALLELESPECIFICDATA){
		string a;
		double s;
		for(int i = 0; i< nLength; ++i){
			if (calculateAllelicConfiguration(i,cn[i],a,s)){
				allelic_cons[i]=a;
				allelic_scores[i]=s;
			}
		}
	}

	string fileName("chr");
	fileName += chrSymbol;
	int sPos =(string(snpdataname).find_last_of("/")==string::npos)?0:string(snpdataname).find_last_of("/")+1;
	fileName += string(snpdataname).substr(sPos,string(snpdataname).length());
	sPos =(string(infodataname).find_last_of("/")==string::npos)?0:string(infodataname).find_last_of("/")+1;
	fileName += string(infodataname).substr(sPos,string(infodataname).length());
	string snpName = path + "Jingerbread_"+fileName+"_SNP.dat";
	string segName = path + "Jingerbread_"+fileName+"_segment.dat";
	ofstream out(snpName.c_str());
	//ofstream out("JSNP.dat");
	ofstream out1(segName.c_str());
	//out << "name\t" << "state\t" << "stateP\t" << "CN\t" << endl;
	//out1 << "chr\t" << "start\t" << "end\t" << "state\t" << "cn\t" << "sample\t" << "snp1\t" << "snp2\t" << "score\t" << "n\t" << endl;
	if (ALLELESPECIFICDATA){
		out << "chr\t" << "str\t" << "end\t" <<"mprop\t" << "gprop\t" << "wincount\t" << "state\t" << "stateP\t" << "CN\t" << "Allelic.Configuration\t" << "Allelic.Score\t" << "a1count\t" << "a2count" << endl;
		out1 << "chr\t" << "start\t" << "end\t" << "state\t" << "cn\t" << "sample\t" << "score\t" << "n\t" << "mscore\t" << "ave.mprop\t" << "ave.gprop\t" << "winct.ratio\t" << "exp.winct\t" << "Allelic.Configuration\t" << "Allelic.TotalScore\t" << "Allelic.AvgScore\t" << "AvgA1count\t" << "AvgA2count" <<endl;
	}else{
		out << "chr\t" << "str\t" << "end\t" <<"mprop\t" << "gprop\t" << "wincount\t" << "state\t" << "stateP\t" << "CN" << endl;
		out1 << "chr\t" << "start\t" << "end\t" << "state\t" << "cn\t" << "sample\t" << "score\t" << "n\t" << "mscore\t" << "ave.mprop\t" << "ave.gprop\t" << "winct.ratio\t" << "exp.winct" << endl;
	}
	
	int start = 0;
	int end = 0;
	bool cnv = false;
	int cnvtype = -1;

	//
	start = 0;
	end = 0;
	cnv =false;
	cnvtype = -1;

	for(int i = 0; i < nLength; ++i)
	{
		if (cn[i] != normalStates)
		{
			if (!cnv)
			{
				cnv = true;
				start = i;
				cnvtype = cn[i];
			}
			else
			{
				if (cnvtype != cn[i])
				{
					end = i-1;
					double score = 0;
					for(int j = start; j <= end; ++j)
						score += exp(pGamma[j][inferenceResults[j]]);
					double mscore = score/(end-start+1);
					double amprop,agprop,rd,expect;
					getSegInfo(start,end,amprop,agprop,rd,expect);
					if (ALLELESPECIFICDATA){
						string a="";
						double s=0;
						double a1count=0,a2count=0;
						if (!getAllelicConfiguration(start, end, cn[start], a, s, a1count, a2count,allelic_cons, allelic_scores)){
							a="";
							s=0;
						}
						out1 << chrSymbol << "\t" << inferData.data[start].startPos << "\t" << inferData.data[end].endPos << "\t"
								<< inferenceResults[start] << "\t" << cn[start] << "\t" << "test\t"
								<< score << "\t" << end-start+1 << "\t" << mscore << "\t" << amprop << "\t" << agprop << "\t" << rd << "\t" << expect << "\t"
								<< a << "\t" << s << "\t" << s/(end-start+1) <<  "\t" << a1count << "\t" << a2count << endl;

					}else{
						out1 << chrSymbol << "\t" << inferData.data[start].startPos << "\t" << inferData.data[end].endPos << "\t"
								<< inferenceResults[start] << "\t" << cn[start] << "\t" << "test\t"
								<< score << "\t" << end-start+1 << "\t" << mscore << "\t" << amprop << "\t" << agprop << "\t" << rd << "\t" << expect << endl;
					}
					// a start of a new CNV
					start = i;
					cnvtype = cn[i];
				}
			}
		}
		else
		{
			if (cnv)
			{
				end = i-1;
				// out a cnv segment
				double score = 0;
				for(int j = start; j <= end; ++j)
					score += exp(pGamma[j][inferenceResults[j]]);
				double mscore = score/(end-start+1);
				double amprop,agprop,rd,expect;
				getSegInfo(start,end,amprop,agprop,rd,expect);
				if (ALLELESPECIFICDATA){
					string a="";
					double s=0;
					double a1count=0, a2count=0;
					if (!getAllelicConfiguration(start, end, cn[start], a, s, a1count, a2count, allelic_cons, allelic_scores)){
						a="";
						s=0;
					}
					out1 << chrSymbol << "\t" << inferData.data[start].startPos << "\t" << inferData.data[end].endPos << "\t"
							<< inferenceResults[start] << "\t" << cn[start] << "\t" << "test\t"
							<< score << "\t" << end-start+1 << "\t" << mscore << "\t" << amprop << "\t" << agprop << "\t" << rd << "\t" << expect << "\t"
							<< a << "\t" << s << "\t" << s/(end-start+1) <<  "\t" << a1count << "\t" << a2count << endl;

				}else{
					out1 << chrSymbol << "\t" << inferData.data[start].startPos << "\t" << inferData.data[end].endPos << "\t"
							<< inferenceResults[start] << "\t" << cn[start] << "\t" << "test\t"
							<< score << "\t" << end-start+1 << "\t" << mscore << "\t" << amprop << "\t" << agprop << "\t" << rd << "\t" << expect << endl;
				}

			}
			cnv = false;
		}
		if (ALLELESPECIFICDATA){
			out << chrSymbol << "\t" << inferData.data[i].startPos << "\t" << inferData.data[i].endPos <<"\t" << inferData.data[i].mprop <<"\t" << inferData.data[i].gprop
					<<"\t" << inferData.data[i].count <<"\t" <<  inferenceResults[i] << "\t" << exp(pGamma[i][inferenceResults[i]]) <<"\t"
				<< cn[i] << "\t" << allelic_cons[i] << "\t" << allelic_scores[i] << "\t" << inferData.data[i].a1count << "\t" << inferData.data[i].a2count << endl;
		}else{
		out << chrSymbol << "\t" << inferData.data[i].startPos << "\t" << inferData.data[i].endPos <<"\t" << inferData.data[i].mprop <<"\t" << inferData.data[i].gprop
				<<"\t" << inferData.data[i].count <<"\t" <<  inferenceResults[i] << "\t" << exp(pGamma[i][inferenceResults[i]]) <<"\t"
			<< cn[i] << endl;
		}
	}

	if (cnv)
	{
		end = nLength-1;
		// out a cnv segment
		double score = 0;
		for(int j = start; j <= end; ++j)
			score += exp(pGamma[j][inferenceResults[j]]);
		double mscore = score/(end-start+1);
		double amprop,agprop,rd,expect;
		getSegInfo(start,end,amprop,agprop,rd,expect);
		if (ALLELESPECIFICDATA){
			string a="";
			double s=0;
			double a1count=0,a2count=0;
			if (!getAllelicConfiguration(start, end, cn[start], a, s, a1count, a2count, allelic_cons, allelic_scores)){
				a="";
				s=0;
			}
			out1 << chrSymbol << "\t" << inferData.data[start].startPos << "\t" << inferData.data[end].endPos << "\t"
					<< inferenceResults[start] << "\t" << cn[start] << "\t" << "test\t"
					<< score << "\t" << end-start+1 << "\t" << mscore << "\t" << amprop << "\t" << agprop << "\t" << rd << "\t" << expect << "\t"
					<< a << "\t" << s << "\t" << s/(end-start+1) <<  "\t" << a1count << "\t" << a2count << endl;

		}else{
			out1 << chrSymbol << "\t" << inferData.data[start].startPos << "\t" << inferData.data[end].endPos << "\t"
					<< inferenceResults[start] << "\t" << cn[start] << "\t" << "test\t"
					<< score << "\t" << end-start+1 << "\t" << mscore << "\t" << amprop << "\t" << agprop << "\t" << rd << "\t" << expect << endl;
		}

	}


	out.close();
	out1.close();
	delete []cn;
	cn = NULL;

}

void HMModel::getSegInfo(int l, int r, double &avemprop, double &avegprop, double &rd, double &expect)
{
	// calculate the average mprop, gprop, win count
	int s = r-l+1;
	double mprop=0, gprop=0, wc=0;
	for(int i = l; i<=r; ++i)
	{
		mprop += inferData.data[i].mprop;
		gprop += inferData.data[i].gprop;
		wc += inferData.data[i].count;
	}
	avemprop = mprop/s;
	avegprop = gprop/s;
	wc /= s;
	double expw = inferData.rd.getValue(avemprop,avegprop);
	rd = wc / (expw+0.00001);
	expect=expw;
}



void HMModel::reEstimation(bool transitionReestimate, bool initReestimation)
{
    if (initReestimation)
    {
		cout << "initial probability re-estimated" << endl;
		// update initial probability
		for(int i = 0; i < nSTATES; ++i)
		{
			pPi[i] = exp(pAlpha[0][i] + pBeta[0][i] - cLikelihood[nITRATION-1]);
		}
    }

	if (transitionReestimate)
	{
		cout << "transition re-estimated" << endl;
		// update transition probability
		// first we need to create a temp transition matrix to store new values
		double **newTran = new double*[nSTATES];
		for(int i = 0; i < nSTATES; ++i)
		{
			newTran[i] = new double[nSTATES];
			for(int j = 0; j < nSTATES; ++j)
				newTran[i][j] = 0;
		}
		// then we create a temp cjk function
		double **c = new double *[nSTATES];
		for(int i = 0; i < nSTATES; ++i)
		{
			c[i] = new double[nSTATES];
			for(int j = 0; j < nSTATES; ++j)
				c[i][j] = 0;
		}
		double *v = new double[nLength-1];
		for(int j = 0; j < nSTATES; ++j)
		{
			for(int k = 0; k < nSTATES; ++k)
			{
				for(int i = 1; i < nLength; ++i)
				{
					v[i-1] = pAlpha[i-1][j]+log(pEmissTbl[i][k])+pBeta[i][k];
				}
				c[j][k] = MathTools::logsumexp(v, nLength-1);
			}
		}
		delete []v;
		v = new double[nSTATES];
		for(int j = 0; j < nSTATES; ++j)
		{
			for(int k = 0; k < nSTATES; ++k)
			{
				newTran[j][k] = log(pTranTbl[j][k]);
				newTran[j][k] += c[j][k];
				for(int l = 0; l < nSTATES; ++l)
				{
					v[l] = log(pTranTbl[j][l])+c[j][l];
				}
				double vsum = MathTools::logsumexp(v, nSTATES);
				newTran[j][k] -= vsum;
				newTran[j][k] = exp(newTran[j][k]);
			}
		}
		delete []v;
		// put the new value to Tran Table
		for(int i = 0; i < nSTATES; ++i)
		{
			for(int j = 0; j < nSTATES; ++j)
				pTranTbl[i][j] = newTran[i][j];
		}
		fillTranDiscrete();
		for(int i = 0; i < nSTATES; ++i)
		{
			delete []newTran[i];
			newTran[i] = NULL;
		}
		delete []newTran;
		newTran = NULL;
		for(int i = 0; i < nSTATES; ++i)
		{
			delete []c[i];
			c[i] = NULL;
		}
		delete []c;
		c = NULL;
	}

	calculateMuAndPhiAllStatesCombined();
	if (USINGAUTOREGRESSION)
		calculateMuAndPhiWithAutoRegressionAllStatesCombined();


	fillEmissionTbl();

}



void HMModel::fillEmissionTbl(void)
{

	for(int i = 0; i < nLength; ++i)
	{
		for(int j = 0; j < nSTATES; ++j)
		{
			fillEmissTblItem(i,j);
		}
	}
}

void HMModel::fillEmissTblItem(int site, int state)
{

	double m = mu[site][state];
	double y = inferData.data[site].count;
	double w = 1; // weight is not useful when calculating the emission probability
	double e = exp(MathTools::loglik_NB(1, phi[state], &m, &y, &w));
	//pEmissTbl[site][state] = e;
	if (state==2)
	{
		pEmissTbl[site][state] = mixtureProportionNormal/largestReadCount + (1-mixtureProportionNormal)*e;
	}
	else
	{
		pEmissTbl[site][state] = mixtureProportion/largestReadCount + (1-mixtureProportion)*e;
	}


	double alleleeffect=1.0;
	double alleleCount=inferData.data[site].a1count + inferData.data[site].a2count;
	if (ALLELESPECIFICDATA && alleleCount>0){

		switch(state){
		case 0:
			//alleleeffect=MathTools::betaBinomial(inferData.data[site].count,inferData.data[site].a1count,0.99,0.1);
			alleleeffect=1/std::max(1.0,alleleCount);
			break;
		case 1:
			alleleeffect=0.5*MathTools::betaBinomial(alleleCount,inferData.data[site].a1count,0.01,0.1)
			                     +0.5*MathTools::betaBinomial(alleleCount,inferData.data[site].a1count,0.99,0.1);
			break;
		case 2:
			alleleeffect=MathTools::betaBinomial(alleleCount,inferData.data[site].a1count,0.5,0.1);
			break;
		case 3:
			alleleeffect=0.5*MathTools::betaBinomial(alleleCount,inferData.data[site].a1count,0.33,0.1)
						                     +0.5*MathTools::betaBinomial(alleleCount,inferData.data[site].a1count,0.67,0.1);
			break;
		case 4:
			alleleeffect=0.33*MathTools::betaBinomial(alleleCount,inferData.data[site].a1count,0.25,0.1)
			                               +0.33*MathTools::betaBinomial(alleleCount,inferData.data[site].a1count,0.5,0.1)
			                               +0.33*MathTools::betaBinomial(alleleCount,inferData.data[site].a1count,0.75,0.1);
			break;
		case 5:
			alleleeffect=0.25*MathTools::betaBinomial(alleleCount,inferData.data[site].a1count,0.2,0.1)
			                               +0.25*MathTools::betaBinomial(alleleCount,inferData.data[site].a1count,0.4,0.1)
			                               +0.25*MathTools::betaBinomial(alleleCount,inferData.data[site].a1count,0.6,0.1)
			                               +0.25*MathTools::betaBinomial(alleleCount,inferData.data[site].a1count,0.8,0.1);
			break;
		case 6:
			alleleeffect=0.17*MathTools::betaBinomial(alleleCount,inferData.data[site].a1count,0.17,0.1)
			                               +0.17*MathTools::betaBinomial(alleleCount,inferData.data[site].a1count,0.33,0.1)
			                               +0.17*MathTools::betaBinomial(alleleCount,inferData.data[site].a1count,0.5,0.1)
										   +0.17*MathTools::betaBinomial(alleleCount,inferData.data[site].a1count,0.67,0.1)
			                               +0.17*MathTools::betaBinomial(alleleCount,inferData.data[site].a1count,0.83,0.1);
			break;
		default:
			alleleeffect=1.0;
			break;
		}
	}

	pEmissTbl[site][state] *= alleleeffect;

//	if (site>=505400 && site<=505410){
//		cout << "site: " << site << " state: "<< state << " starting pos " << inferData.data[site].startPos << " phi: " << phi[state] << " mu: " << m << " count:  " << y << "emission: " << log(pEmissTbl[site][state]) << "allele effect is :" << alleleeffect << endl;
//	}

	if (pEmissTbl[site][state] <1e-12)
		pEmissTbl[site][state] = 1e-12;

}

void HMModel::setFileName(char * sn, char * in)
{
	snpdataname = sn;
	infodataname = in;

	int sPos =(string(snpdataname).find_last_of("/")==string::npos)?0:string(snpdataname).find_last_of("/")+1;
    path = "";
    if (sPos==0)
    	path="./";
    else
    	path=string(snpdataname).substr(0,sPos);

}

void HMModel::printVariable(void)
{
	for(int i = 0; i < nITRATION; ++i)
		cout << "The log likelihood value for the " << i+1 << " round inference is " << cLikelihood[i] << endl;
}


void HMModel::fillTranContinous()
{
	for(int i = 0; i < nSTATES; ++i)
	{
		for(int j = 0; j < nSTATES; ++j)
		{
			pTran[0][i][j] = pTranTbl[i][j];
		}
	}
	int step = 1;
	map<int, double>::iterator it = cLRR.begin();
	int pre = (*it).first;
	++it;
	for(; it != cLRR.end(); ++it, ++step)
	{
		int cur = (*it).first;
		for(int i = 0; i < nSTATES; ++i)
		{
			double s = 0;
			for(int j = 0; j < nSTATES; ++j)
			{
				if (i != j)
				{
					double temp = pTranTbl[i][j]*(1-exp(-lamda*(cur-pre)));
					if (temp > 1e-10)
						pTran[step][i][j] = temp;
					else
						pTran[step][i][j] = 1e-10;
					s+=pTran[step][i][j];
				}
				else
				{
					pTran[step][i][j] = 0;
				}
			}
			pTran[step][i][i] = 1-s;
		}
		pre = cur;
	}
}

void HMModel::fillTranDiscrete()
{
	for(int step = 0; step < nLength; ++step)
		for(int i = 0; i < nSTATES; ++i)
			for(int j = 0; j < nSTATES; ++j)
			{
				if (pTranTbl[i][j] > 1e-10)
					pTran[step][i][j] = pTranTbl[i][j];
				else
					pTran[step][i][j] = 1e-10;
				if (step<nLength-1){
					int dis=inferData.data[step+1].startPos-inferData.data[step].endPos;
					double f=exp(-dis/inferData.meanDis);
					pTran[step][i][j]=f*(pTranTbl[i][j])+(1-f)*(pTranTbl[normalStates][j]);
				}
			}
}

void HMModel::findBestPath(bool viterbi)
{
	if (viterbi)
	{
		// viterbi
		double **v = new double*[nLength];  // in log scale
		int **pathm = new int*[nLength];
		for(int i = 0; i < nLength; ++i)
		{
			v[i] = new double[nSTATES];
			pathm[i] = new int[nSTATES];
		}

		for(int i = 0; i < nSTATES; ++i)
		{
			v[0][i] = log(pPi[i])+log(pEmissTbl[0][i]);
		}

		for(int i = 1; i < nLength; ++i)
		{
			for(int z = 0; z < nSTATES; ++z)
			{
				double maxValue = -1e10;
				int index = -1;
				for(int j = 0; j < nSTATES; ++j)
				{
					double value = v[i-1][j]+log(pTran[i][j][z]);
					if (value > maxValue)
					{
						maxValue = value;
						index = j;
					}
				}
				v[i][z] = maxValue + log(pEmissTbl[i][z]);
				pathm[i-1][z] = index;
			}
		}

		ofstream fout("v.log");
		for(int i = 0; i < nLength; ++i)
		{
			for(int j = 0; j < nSTATES; ++j)
				fout << v[i][j] << " ";
			fout << endl;
		}
		fout.close();	

		fout.open("path.log");
		for(int i = 0; i < nLength; ++i)
		{
			for(int j = 0; j < nSTATES; ++j)
				fout << pathm[i][j] << " ";
			fout << endl;
		}
		fout.close();



		double maxValue = -1e10;
		int index = -1;
		for(int z = 0; z < nSTATES; ++z)
		{
			if (v[nLength-1][z] > maxValue)
			{
				maxValue = v[nLength-1][z];
				index = z;
			}
		}

		inferenceResults[nLength-1] = index;
		cout << "The log probability of the most likely path is " << maxValue << endl;

		for(int i = nLength -2; i >= 0; --i)
			inferenceResults[i] = pathm[i][inferenceResults[i+1]];



		for(int i = 0; i < nLength; ++i)
		{
			delete []v[i];
			v[i] = NULL;
			delete []pathm[i];
			pathm[i] = NULL;
		}
		delete []v;
		v = NULL;
		delete []pathm;
		pathm = NULL;
	}
	else
	{

		for(int i = 0; i < nLength; ++i)
		{
			double delP = 0;
			double norP = exp(pGamma[i][normalStates]);
			double dupP = 0;
			int index = -1;
			double max = -1;
			for(int j = 0; j < nSTATES; ++j)
			{

				if (exp(pGamma[i][j]) > max)
				{
					index = j;
					max = exp(pGamma[i][j]);
				}
			}
			inferenceResults[i] = index;
		}
	}

}

void HMModel::writeKeyValueTable(char *filename, double ** data)
{
	ofstream fout(filename);
	fout.precision(6);
	fout << "chr  " << "str  " << "end  ";
	for(int i = 0; i < nSTATES; ++i)
		fout << "column  ";
	fout << endl;


	for(int i = 0; i < nLength; ++i)
	{
		fout << inferData.chr << " " << inferData.data[i].startPos << " " << inferData.data[i].endPos << " ";
		for(int j = 0; j < nSTATES; ++j)
		{
			fout << data[i][j] << " ";
		}
	}

	fout.close();
}

void HMModel::writeKeyValue(int index)
{
	// 
	stringstream s;
	s << index;
	string postfix = s.str()+".log";
	// alpha,beta,gamma,transition,emission
	string filename("alpha");
	filename += postfix;
	filename = path+filename;
//	writeKeyValueTable(filename.c_str(), pAlpha);
	ofstream fout(filename.c_str());
	fout.precision(8);
	fout << "chr  " << "str  " << "end  ";
    for(int i = 0; i < nSTATES; ++i)
		fout << "column  ";
	fout << endl;
	for(int i = 0; i < nLength; ++i)
	{
		fout << inferData.chr << " " << inferData.data[i].startPos << " " << inferData.data[i].endPos << " ";
		for(int j = 0; j < nSTATES; ++j)
			fout << pAlpha[i][j] << " ";
		fout << endl;
	}
	fout.close();

	filename = "beta";
	filename += postfix;
	filename = path+filename;
//	writeKeyValueTable(filename.c_str(), pBeta);
	fout.open(filename.c_str());
	fout << "chr  " << "str  " << "end  ";
    for(int i = 0; i < nSTATES; ++i)
		fout << "column  ";
	fout << endl;
	for(int i = 0; i < nLength; ++i)
	{
		fout << inferData.chr << " " << inferData.data[i].startPos << " " << inferData.data[i].endPos << " ";
		for(int j = 0; j < nSTATES; ++j)
			fout << pBeta[i][j] << " ";
		fout << endl;
	}
	fout.close();

	filename = "gamma";
	filename += postfix;
	filename = path+filename;
//	writeKeyValueTable(filename.c_str(), pGamma);
	fout.open(filename.c_str());
	fout << "chr  " << "str  " << "end  ";
	for(int i = 0; i < nSTATES; ++i)
		fout << "column  ";
	fout << "total ";
	fout << endl;
	for(int i = 0; i < nLength; ++i)
	{
		double sum = 0;
		fout << inferData.chr << " " << inferData.data[i].startPos << " " << inferData.data[i].endPos << " ";
		for(int j = 0; j < nSTATES; ++j)
		{
			fout << exp(pGamma[i][j]) << " ";
			sum += exp(pGamma[i][j]);
		}
		fout << sum;
		fout << endl;
	}
	fout.close();

	filename = "transition";
	filename += postfix;
	filename = path+filename;
	fout.open(filename.c_str());
	for(int step=0; step<nLength;++step){
		for(int i = 0; i < nSTATES; ++i)
		{
			for(int j = 0; j < nSTATES; ++j)
				fout << pTran[step][i][j] << " ";
			fout << endl;
		}
		fout << "-------------------------------------------------------" << endl;
	}
	fout.close();

	filename = "emission";
	filename += postfix;
	filename = path+filename;
	fout.open(filename.c_str());
	fout << "chr  " << "str  " << "end  ";
	for(int i = 0; i < nSTATES; ++i)
		fout << "column  ";
	fout << endl;
	for(int i = 0; i < nLength; ++i)
	{
		fout << inferData.chr << " " << inferData.data[i].startPos << " " << inferData.data[i].endPos << " ";
		for(int j = 0; j < nSTATES; ++j)
			fout << log(pEmissTbl[i][j]) << " ";
		fout << endl;
	}
	fout.close();

	filename = "mu";
	filename += postfix;
	filename = path+filename;
	fout.open(filename.c_str());
	fout << "chr  " << "str  " << "end  ";
	for(int i = 0; i < nSTATES; ++i)
		fout << "column  ";
	fout << endl;
	for(int i = 0; i < nLength; ++i)
	{
		fout << inferData.chr << " " << inferData.data[i].startPos << " " << inferData.data[i].endPos << " ";
		for(int j = 0; j < nSTATES; ++j)
			fout << mu[i][j] << " ";
		fout << endl;
	}
	fout.close();

	filename = "phi";
	filename += postfix;
	filename = path+filename;
	fout.open(filename.c_str());
	for(int j = 0; j < nSTATES; ++j)
		fout << phi[j] << " ";
	fout.close();

	filename = "init";
	filename += postfix;
	filename = path+filename;
	fout.open(filename.c_str());
	for(int j = 0; j < nSTATES; ++j)
		fout << pPi[j] << " ";
	fout.close();
}

//find the most likely allelic configuration of a window, according to its most likely state

bool HMModel::calculateAllelicConfiguration(int pos, int state, std::string &allelic_con, double &score){
	allelic_con="";
	score=0;
	vector<int> acntypes; // number of allelic configurations of each state
	map<int, vector<string> > c; // allelic configurations of each state
	map<int, vector<double> > v; // proporation of selecting one allele of each state
	getAllelicConfigTable(acntypes, c, v);


	if (inferData.data[pos].a1count>0 || inferData.data[pos].a2count>0){
		double acount=inferData.data[pos].a1count+inferData.data[pos].a2count;
		vector<double> posteriori(7,0.0);
		switch(state){
		case 0:
			allelic_con=c[0][0];
			score=exp(pGamma[pos][0]);
			break;
		case 2:
			allelic_con=c[2][0];
			score=exp(pGamma[pos][2]);
			break;
		default:
			int n=acntypes[state];
			vector<double> val(n,0);
			for(int i=0; i<n; ++i){
				val[i]=MathTools::betaBinomial(acount,inferData.data[pos].a1count,v[state][i],0.1);
			}
			int mcon=-1;
			double maxScore=0;
			double totalScore=0;
			for(int i=0; i<n; ++i){
				if (mcon==-1){
					mcon=i;
					maxScore=val[i];
				}
				else{
					if (val[i]>maxScore){
						mcon=i;
						maxScore=val[i];
					}
				}
				totalScore+=val[i];
			}
			allelic_con=c[state][mcon];
			score=exp(pGamma[pos][state])*maxScore/totalScore;
//			if (pos==186){
//				cout << "wincount is " << inferData.data[pos].count << "a1count is " << inferData.data[pos].a1count << "a2count is " << inferData.data[pos].a2count << endl;
//				cout << "state is: " << state << endl;
//				cout << "n: " << n << endl;
//				cout << "acount: " << acount << endl;
//				for(int i=0; i<n; ++i){
//					cout << "v[state][i] " << state << " " << i << ": is " << v[state][i] << " and the value[i] is " << val[i] << endl;
//				}
//                cout << "the detected max con is " << mcon << endl;
//                cout << "the maximum value is " << maxScore << endl;
//                cout << "the choosen allelic_con is " << allelic_con << endl;
//			}
			break;
		}
		return true;
	}else{
		return false;
	}
}


// find the majority allelic configuration of the CNV bounded by l_bound, r_bound
// calculate the aveconfidence score of the majority allelic configuration
// return whether the CNV has allelic configuration,
// and if yes, return its allelic configuration, and the total confidence score, avg a1count avg a2count,

bool HMModel::getAllelicConfiguration(int l_bound, int r_bound, int state, string &allelic_con, double &score, double &a1count, double &a2count, vector<string> &allelic_cons, vector<double> &scores){
	if (r_bound<l_bound){
		return false;
	}
	int n=r_bound-l_bound+1;
	score=0;
	allelic_con="";
	a1count=0;
	a2count=0;
	// to find, first needs to count, 1) how many with allelic configurations, 2) the majority allelic configurations
	vector<int> acntypes; // number of allelic configurations of each state
	map<int, vector<string> > c; // allelic configurations of each state
	map<int, vector<double> > v; // proporation of selecting one allele of each state
	getAllelicConfigTable(acntypes, c, v);
	int types=acntypes[state];
	vector<int> counts(types+1,0); // the last element will be the number of no allelic specific windows
	for(int i=0; i<n; ++i){
		a1count+=inferData.data[l_bound+i].a1count;
		a2count+=inferData.data[l_bound+i].a2count;
		if (allelic_cons[l_bound+i]==""){
			counts[types]++;
		}else{
			int idx=-1;
			for(int j=0; j<types; ++j){
				if (allelic_cons[l_bound+i]==c[state][j]){
					idx=j;
					break;
				}
			}
			if (idx==-1){
				cout << "error in getAllelicConfiguration, position: " <<l_bound+i << "state: " << state << "allelic configuration: " << allelic_cons[l_bound+i] << endl;
				return false;
			}else{
				counts[idx]++;
			}
		}
	}
	a1count/=n;
	a2count/=n;
	// which one is the majority?
	int idx=-1;
	int maxCount=-1;
	for(int i=0; i<types+1;++i){
		if (idx==-1){
			idx=i;
			maxCount=counts[idx];
		}else{
			if (counts[i]>maxCount){
				idx=i;
				maxCount=counts[i];
			}
		}
	}
	if (idx==types){
		return false;
	}else{
		allelic_con=c[state][idx];
		for(int i=0; i<n; ++i){
			if (allelic_cons[l_bound+i]==allelic_con){
				score+=scores[l_bound+i];
			}
		}
		return true;
	}

}

void HMModel::getAllelicConfigTable(vector<int> &acntypes, map<int, vector<string> > &c, map<int, vector<double> > &v){
	acntypes.clear();
	acntypes.push_back(1); // state 0
	acntypes.push_back(2); // state 1
	acntypes.push_back(1); // state 2
	acntypes.push_back(2); // state 3
	acntypes.push_back(3); // state 4
	acntypes.push_back(4); // state 5
	acntypes.push_back(5); // state 6
	c.clear(); // allelic configurations of each state
	vector<string> t;
	t.push_back("-");
	c[0]=t;
	t.clear();
	t.push_back("AB");
	c[2]=t;
	for(int i=1; i<=6; ++i){
		if (i==1){
			t.clear();
			t.push_back("B");
			t.push_back("A");
			c[i]=t;
			continue;
		}
		if (i==2) continue;
		t.clear();
		int n=acntypes[i];
		for(int j=0; j<n; ++j){
			string config="";
			int nA=(j+1)>i?0:(j+1);
			int nB=i-nA;
			string Apart(nA,'A');
			string Bpart(nB,'B');
			config=Apart+Bpart;
			t.push_back(config);
		}
		c[i]=t;
	}
	v.clear(); // proporation of selecting one allele of each state
	vector<double> tt;
	tt.push_back(0.99);
	v[0]=tt; // actually, for state 0, we use uniform distribution rather than beta binormial distribution. Here just for the sake of coding
	tt.clear();
	tt.push_back(0.01);
	tt.push_back(0.99);
	v[1]=tt;
	tt.clear();
	tt.push_back(0.5);
	v[2]=tt;
	tt.clear();
	tt.push_back(0.33);
	tt.push_back(0.67);
	v[3]=tt;
	tt.clear();
	tt.push_back(0.25);
	tt.push_back(0.5);
	tt.push_back(0.75);
	v[4]=tt;
	tt.clear();
	tt.push_back(0.2);
	tt.push_back(0.4);
	tt.push_back(0.6);
	tt.push_back(0.8);
	v[5]=tt;
	tt.clear();
	tt.push_back(0.17);
	tt.push_back(0.33);
	tt.push_back(0.5);
	tt.push_back(0.67);
	tt.push_back(0.83);
	v[6]=tt;
	tt.clear();
}

////////////////////////////////////// codes that are not active ///////////////////

