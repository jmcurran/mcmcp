#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include "Fastnorm.h"
#include <time.h>

#include <string>
#include <vector>

using namespace std;

double rgamma(double a, double scale);
double lgamma(double arg);


double rbeta(double a, double b)
{
	double g1 = rgamma(a, 1);
	double g2 = rgamma(b, 1);

	return g1/(g1+g2);
}

double logdgamma(double x, double a, double shape)
{
	return (a-1)*log(x)-x/shape;
}

double logdbeta(double x, double a, double b)
{
	return (a-1)*log(x)+(b-1)*log(1-x)+lgamma(a+b)-lgamma(a)-lgamma(b);
}

double fmax(double a, double b)
{
	if(a>b)
		return a;
	return b;
}

static unsigned long next1, next2;


int G2[][4]={
	{1,2,1,1},
	{2,2,1,1},
	{1,1,1,2},
	{1,2,1,2},
	{2,2,1,2},
	{1,1,2,2},
	{1,2,2,2}};


int G2X[][4]={
	{1,2,1,1},
	{1,1,1,2},
	{1,2,1,2}};




int G3[][4]={
	{2,3,1,1},
	{1,3,1,2},
	{2,3,1,2},
	{3,3,1,2},
	{1,2,1,3},
	{2,2,1,3},
	{2,3,1,3},
	{1,3,2,2},
	{1,1,2,3},
	{1,2,2,3},
	{1,3,2,3},
	{1,2,3,3}};

int G4[][4]={
	{3,4,1,2},
	{2,4,1,3},
	{2,3,1,4},
	{1,4,2,3},
	{1,3,2,4},
	{1,2,3,4}};

double logdchisq(double x, double df)
{
	return (0.5*df-1)*log(x)-0.5*x;
}



static int ix, iy, iz; //  Wichman Hill seeds
//static unsigned long next1, next2;

void init_generator(int nx, int ny, int nz/*int n1, int n2*/)
{
	/*next1=n1;
	next2=n2;*/
	ix = nx;
	iy = ny;
	iz = nz;
}


double runif()
{
 /*next1=30903*(next1 & 65535L)+(next1 >> 16);
 next2=18000*(next2 & 65535L)+(next2 >> 16);
 return((next1<<16)+next2)/4294967295.0;*/

/*Wichmann & Hill random number generator */

  ix=(171*ix)%30269;
  iy=(172*iy)%30307;
  iz=(170*iz)%30323;
  return fmod((ix/30269.0)+(iy/30307.0)+(iz/30323.0),1.0);
}




class Locus{
	string m_strName;
	int m_nPeaks;
	double m_dArea[4];
	double m_dTotal;
	vector<string> m_vAlleles;
	int m_nGenotype[4];

public:
	Locus(){
		m_nPeaks=0;
		
		int i;

		for(i=0;i<4;i++){
			m_dArea[i]=0;
			m_nGenotype[i]=0;
		}
		m_dTotal=0;
	};

	Locus(const string& strName, int nPeaks, const string& strAreas, const string& strAlleles){
		m_strName = strName;
		m_nPeaks=nPeaks;
		
		int i;
		for(i=0;i<4;i++){
			m_dArea[i]=0;
			m_nGenotype[i]=-1;
		}

		string strTmp = strAreas;
		string strVal;
		string::size_type nPos=strTmp.find(",");
		
		i=0;
		m_dTotal=0;
		while(nPos!=string::npos){
			strVal = strTmp.substr(0, nPos);
			strTmp = strTmp.substr(nPos+1);
			m_dArea[i]=atof(strVal.c_str());
			m_dTotal+=m_dArea[i];
			i++;
			nPos=strTmp.find(",");
		}
		m_dArea[i]=atof(strTmp.c_str());
		m_dTotal+=m_dArea[i];


		strTmp = strAlleles;
		nPos=strTmp.find(",");
		
		while(nPos!=string::npos){
			strVal = strTmp.substr(0, nPos);
			strTmp = strTmp.substr(nPos+1);
			m_vAlleles.push_back(strVal);
			nPos=strTmp.find(",");
		}
		m_vAlleles.push_back(strTmp);

	};

	Locus(const Locus& l){
		m_strName = l.m_strName;
		m_nPeaks = l.m_nPeaks;

		int i;
		for(i=0;i<4;i++){
			m_dArea[i]=l.m_dArea[i];
			m_nGenotype[i]=l.m_nGenotype[i];
		}

		m_vAlleles = l.m_vAlleles;
		m_dTotal = l.m_dTotal;
	};

	const Locus& operator=(const Locus& l){
		m_strName = l.m_strName;
		m_nPeaks = l.m_nPeaks;

		int i;
		for(i=0;i<4;i++){
			m_dArea[i]=l.m_dArea[i];
			m_nGenotype[i]=l.m_nGenotype[i];
		}

		m_vAlleles = l.m_vAlleles;
		m_dTotal = l.m_dTotal;

		return *this;
	};

	void CopyGenotypes(const Locus& l){
		int i;
		for(i=0;i<4;i++)
			m_nGenotype[i]=l.m_nGenotype[i];
	};


	double loglik(double mx){
		double dExpected[]={0,0,0,0};
		
		dExpected[m_nGenotype[0]]=mx*m_dTotal*0.5;
		dExpected[m_nGenotype[1]]+=mx*m_dTotal*0.5;
		dExpected[m_nGenotype[2]]+=(1-mx)*m_dTotal*0.5;
		dExpected[m_nGenotype[3]]+=(1-mx)*m_dTotal*0.5;

		double dSum=0;
		int i;
		double o, e;

		for(i=0;i<m_nPeaks;i++){
			o=m_dArea[i];
			e=dExpected[i];
			dSum+=(o-e)*(o-e)/e;
		}

		return dSum;
	}

	void rand(){
		int g;

		switch(m_nPeaks){
		case 1:
			m_nGenotype[0]=0;
			m_nGenotype[1]=0;
			m_nGenotype[2]=0;
			m_nGenotype[3]=0;
			break;
		case 2:
			if(m_strName=="AM"){
				g=floor(3*runif());
				m_nGenotype[0]=G2X[g][0]-1;
				m_nGenotype[1]=G2X[g][1]-1;
				m_nGenotype[2]=G2X[g][2]-1;
				m_nGenotype[3]=G2X[g][3]-1;	
			}else{
				g=floor(7*runif());
					
				m_nGenotype[0]=G2[g][0]-1;
				m_nGenotype[1]=G2[g][1]-1;
				m_nGenotype[2]=G2[g][2]-1;
				m_nGenotype[3]=G2[g][3]-1;
			}
			break;
		case 3:
			g=floor(12*runif());
				
			m_nGenotype[0]=G3[g][0]-1;
			m_nGenotype[1]=G3[g][1]-1;
			m_nGenotype[2]=G3[g][2]-1;
			m_nGenotype[3]=G3[g][3]-1;	
			break;
		case 4:
			g=floor(6*runif());			
			m_nGenotype[0]=G4[g][0]-1;
			m_nGenotype[1]=G4[g][1]-1;
			m_nGenotype[2]=G4[g][2]-1;
			m_nGenotype[3]=G4[g][3]-1;
	
		/*	m_nGenotype[0]=3-1;
			m_nGenotype[1]=4-1;
			m_nGenotype[2]=1-1;
			m_nGenotype[3]=2-1; */	
			break;
		}
	}

	void Scale(){
		int i;

		for(i=0;i<m_nPeaks;i++)
			m_dArea[i]/=m_dTotal;
		m_dTotal=1;
	};


	string Genotype(double dMx, bool bMinor=true){
		ostringstream oss;
		string strResult;
		string strA1, strA2;

		if(bMinor){
			strA1=m_vAlleles[m_nGenotype[0]];
			strA2=m_vAlleles[m_nGenotype[1]];
		}
		/*else if(dMx>0.5&&bMinor){
			strA1=m_vAlleles[m_nGenotype[2]];
			strA2=m_vAlleles[m_nGenotype[3]];
		}else if(dMx<0.5&&!bMinor){
			strA1=m_vAlleles[m_nGenotype[0]];
			strA2=m_vAlleles[m_nGenotype[1]];
		}*/
		else{
			strA1=m_vAlleles[m_nGenotype[2]];
			strA2=m_vAlleles[m_nGenotype[3]];
		}

		if(m_strName=="AM"){
			if(strA1=="Y"&&strA2=="X"){
				string strTemp = strA1;
				strA1 = strA2;
				strA2 = strTemp;
			}
		}else{
			if(atof(strA1.c_str())>atof(strA2.c_str())){
				string strTemp = strA1;
				strA1 = strA2;
				strA2 = strTemp;
			}
		}


		oss << "\"" << strA1 << "/" << strA2 << "\"";
		strResult= oss.str();
		return strResult;
	};
};

class Profile{
public:
	int m_nLoci;
	Locus *m_pLoci;

public:
	Profile(){
		m_nLoci=0;
		m_pLoci=NULL;
	};

	Profile(int nLoci){
		m_nLoci=nLoci;
		m_pLoci=new Locus[nLoci];
	};

	~Profile(){
		m_nLoci=0;
		delete [] m_pLoci;
	};

	Profile(const Profile& p){
		m_nLoci=p.m_nLoci;
		m_pLoci = new Locus[m_nLoci];

		int i;
		for(i=0;i<m_nLoci;i++)
			m_pLoci[i] = Locus(p.m_pLoci[i]);
	};

	const Profile &operator=(const Profile& p){
		m_nLoci=p.m_nLoci;
		if(m_pLoci!=NULL)
			delete [] m_pLoci;
		m_pLoci = new Locus[m_nLoci];

		int i;
		for(i=0;i<m_nLoci;i++)
			m_pLoci[i] = Locus(p.m_pLoci[i]);
		return *this;
	};

	void CopyGenotypes(const Profile& p){
		int i;
		for(i=0;i<m_nLoci;i++)
			m_pLoci[i].CopyGenotypes(p.m_pLoci[i]);
	};

	Locus& ElementAt(int nIndex){
		return m_pLoci[nIndex];
	};

	Locus GetAt(int nIndex) const{
		//map<int, Combination> &lresult = const_cast<map<int, Combination>&>(m_mapComb);
		
		return m_pLoci[nIndex];
	}

	Locus operator[](int nIndex) const{
		return GetAt(nIndex);
	}

	Locus& operator[](int nIndex){
		return ElementAt(nIndex);
	}

	double loglik(double mx, double *mxl, double *alpha, double a, double b){
		double dSum1=0;
		double dSum2=logdbeta(mx/0.48,alpha[0],alpha[1]);
		int i;

		for(i=0;i<m_nLoci;i++){
			dSum1+=m_pLoci[i].loglik(mxl[i]);
			dSum2+=logdbeta(mxl[i],a, b);
		}

		return logdchisq(dSum1,4*m_nLoci-1)+dSum2+logdgamma(alpha[0], 1, 1000)+logdgamma(alpha[1], 1, 1000);


	}

	void Read(const string& strFileName){
		ifstream f1(strFileName.c_str());
		if(f1.is_open()){
		/*	vector<string> vLoci(11);
			vLoci[0] = "D3";
			vLoci[1] = "vWA";
			vLoci[2] = "D16";
			vLoci[3] = "D2";
			vLoci[4] = "AM";
			vLoci[5] = "D8";
			vLoci[6] = "D21";
			vLoci[7] = "D18";
			vLoci[8] = "D19";
			vLoci[9] = "TH01";
			vLoci[10] = "FGA";*/

		/*	vector<string> vLoci(16);
			
				
			vLoci[0] = "D8";
			vLoci[1] = "D21";
			vLoci[2] = "D7";
			vLoci[3] = "CSF";
			vLoci[4] = "D3";
			vLoci[5] = "TH01";
			vLoci[6] = "D13";
			vLoci[7] = "D16";
			vLoci[8] = "D2";
			vLoci[9] = "D19";
			vLoci[10] = "vWA";
			vLoci[11] = "TPOX";
			vLoci[12] = "D18";
			vLoci[13] = "AM";
			vLoci[14] = "D5";
			vLoci[15] = "FGA"; */

			vector<string> vLoci(13);
			
				
			vLoci[0] = "D3";
			vLoci[1] = "vWA";
			vLoci[2] = "FGA";
			vLoci[3] = "D8";
			vLoci[4] = "D21";
			vLoci[5] = "D18";
			vLoci[6] = "D5";
			vLoci[7] = "D13";
			vLoci[8] = "D7";
			vLoci[9] = "D16";
			vLoci[10] = "TH01";
			vLoci[11] = "TPOX";
			vLoci[12] = "CSF";
			

			int i, j;
			char szLine[100];
			string strLine, strAllele, strPeak;
			ostringstream ossAlleles, ossPeaks;
			int nPeaks,nP;
			int nPos;

			for(i=0;i<m_nLoci;i++){
				nPeaks=0;
				ossAlleles.str("");
				ossPeaks.str("");

				for(j=0;j<4;j++){
					f1.getline(szLine, 100);
					strLine = szLine;
					nPos = strLine.find(",");
					strAllele = strLine.substr(0, nPos);
					strPeak = strLine.substr(nPos+1);
					nP = atoi(strPeak.c_str());
					if(nP!=0){
						ossAlleles << strAllele << ',';
						ossPeaks << strPeak <<  ',';
						nPeaks++;
					}
				}

				strAllele = ossAlleles.str();
				strPeak = ossPeaks.str();

				strAllele = strAllele.substr(0, strAllele.size()-1);
				strPeak = strPeak.substr(0, strPeak.size()-1);

				(*this)[i] = Locus(vLoci[i], nPeaks, strPeak, strAllele);
			}
		}
	};

	void Scale(){
		int i;
		
		for(i=0;i<m_nLoci;i++)
			m_pLoci[i].Scale();
	};

};


double fmin(double a, double b)
{
	if(a<b)
		return a;
	return b;
}

void fitbeta(double m, double dUpperQuantile, double *a, double *b, double dP = 0.99, double dLowerLim = 1, double dUpperLim = 500);

main()
{
	int nLoci = 13;
	Profile g0(nLoci), g1(nLoci);
	double alpha0[2], alpha1[2];
	double a0, b0, a1, b1;
	double mx0, mx1;
	double **mxlocus;

	int i;

	double **ppdLookupTable = new double*[46000];
	ifstream fin("lookuptable.csv");
	string strA0, strA1, strD;

	for(i=0;i<46000;i++)
	{
		ppdLookupTable[i] = new double[2];
		getline(fin, strD, ',');
		getline(fin, strA0, ',');
		getline(fin, strA1);

		ppdLookupTable[i][0] = atof(strA0.c_str());
		ppdLookupTable[i][1] = atof(strA1.c_str());
	}
	fin.close();

	cout << "Read lookup table" << endl;



	time_t s=time(NULL);
	srand(s);

/*	g0[0] = Locus("D3",3,"1910,953,398","16,17,18");
	g0[1] = Locus("vWA",3,"3207,2583,521","17,18,15");
	g0[2] = Locus("D16",4,"1906,1603,874,475","11,9,13,12");
	g0[3] = Locus("D2",4,"2863,1630,1024,1005","17,20,22,23");
	g0[4] = Locus("AM",2,"5609,552","X,Y");
	g0[5] = Locus("D8",3,"2630,2103,272","10,14,13");
	g0[6] = Locus("D21",2,"4525,378","30,29");
	g0[7] = Locus("D18",4,"3042,1860,557,367","17,12,19,15");
	g0[8] = Locus("D19",3,"2212,1753,857","14,13,15");
	g0[9] = Locus("TH01",4,"1514,1410,651,512","6,8,7,9.3");
	g0[10] = Locus("FGA",2,"3929,738","20,25");*/


/*	g0[0] = Locus("D3",3,"4950,3153,1766","17,16,18");
	g0[1] = Locus("vWA",3,"4661,3334,855","18,17,15");
	g0[2] = Locus("D16",4,"3532,2749,691,634","11,9,12,13");
	g0[3] = Locus("D2",4,"4029,3713,898,622","17,20,23,22");
	g0[4] = Locus("AM",2,"9871,902","X,Y");
	g0[5] = Locus("D8",3,"3269,2949,1259","10,14,13");
	g0[6] = Locus("D21",2,"875,6703","29,30");
	g0[7] = Locus("D18",4,"3806,2719,610,584","12,17,15,19");
	g0[8] = Locus("D19",3,"3638,3225,973","14,13,15");
	g0[9] = Locus("TH01",4,"5496,2943,897,614","8,6,9.3,7");
	g0[10] = Locus("FGA",3,"8705,764,501","20,25,22"); */

	g0.Read("WangTW.csv");
//	g0.Scale();

	g1=g0;

	
//	init_generator(rand(),rand());
//	initnorm(rand());
	//init_generator(1,2,3);
	init_generator(4,5,6);
	initnorm(2);

	time_t start, finish;

	time(&start);

	mxlocus = new double*[2];

	mxlocus[0] = new double[nLoci];
	mxlocus[1] = new double[nLoci];

	alpha0[0] = rgamma(1,1000);
	alpha0[1] = rgamma(1,1000);
	
	mx0= 0.02+(rbeta(alpha0[0], alpha0[1])*0.46);
	
	double mode = mx0;
	double sd = 0.035;
	double qupper = mx0+ 2.32*sd;
	double prob = 0.99;

	//fitbeta(mode, qupper, &a0, &b0, prob);

	int nIdx = (int)floor((mx0 - 0.02)*1e5);
	nIdx = nIdx >= 46000 ? 45999 : nIdx;

	a0 = ppdLookupTable[nIdx][0];
	b0 = ppdLookupTable[nIdx][1];

	for(i=0;i<nLoci;i++)
	{
		g0[i].rand();
		mxlocus[0][i] = rbeta(a0,b0);
	}

	double l0 = g0.loglik(mx0, mxlocus[0] ,alpha0, a0, b0);
	double l1,p,u;

	long N=0;
	long nAccept=0, nAcceptOld=0, nRestart=0;
	ofstream f1("results.10k.csv");

	long nIterations = 10000000;
	long nBurnTime = 1000;
	int nSample = 1000;

	for(N=0;N<nIterations+nBurnTime;N++){
		
		g1.CopyGenotypes(g0);
		for(i=0;i<nLoci;i++)
			mxlocus[1][i] = mxlocus[0][i];

		i=(int)(floor(nLoci*runif()));
		
		g1[i].rand();
		alpha1[0] = rgamma(1,1000);
		alpha1[1] = rgamma(1,1000);

		mx1 = 0.02+(rbeta(alpha1[0], alpha1[1])*0.46);


		mode = mx1;
		qupper = mx1+ 2.32*sd;
		
		//fitbeta(mode, qupper, &a1, &b1, prob);

		int nIdx1 = (int)floor((mx1 - 0.02)*1e5);
		nIdx1 = nIdx1 >= 46000 ? 45999 : nIdx1;

		a1 = ppdLookupTable[nIdx1][0];
		b1 = ppdLookupTable[nIdx1][1];

/*		if(fabs(mode-(a1-1)/(a1+b1-2))>1e-2)
			cout << "Error " << a1 << " " << b1 << " " << mode<< " " << N << endl; */

		for(i=0;i<nLoci;i++)
			mxlocus[1][i] = rbeta(a1,b1);

		l1=g1.loglik(mx1, mxlocus[1], alpha1, a1, b1);
		p=fmin(0,l1-l0);
		u=log(runif());

		if(u<=p)
		{
			g0.CopyGenotypes(g1);
			mx0=mx1;
			a0=a1;
			b0=b1;

			for(i=0;i<nLoci;i++)
				mxlocus[0][i]=mxlocus[1][i];
			alpha0[0] = alpha1[0];
			alpha0[1] = alpha1[1];

			l0=l1;
			nAcceptOld = nAccept;
			nAccept++;
		}
		
		if(N>=nBurnTime && N%nSample==0){
			for(i=0;i<nLoci;i++)
				f1 << g1[i].Genotype(1-mxlocus[1][i],false) << " ";
			for(i=0;i<nLoci;i++)
				f1 << g1[i].Genotype(mxlocus[1][i],true) << " ";
			for(i=0;i<nLoci;i++)
				f1 << mxlocus[1][i] << " ";

			f1 << mx1 << " " << alpha1[0] << " " << alpha1[1] << " " << a1 << " " << b1 << endl;
		}

		if(N>0 && N%10000==0)
		{
			ofstream f2("seeds.dat");
			f2 << next1 << endl;
			f2 << next2 << endl;
			f2.close();
			
			cout << N << " " << (double)nAccept/(double)N << " " << nAccept << endl;
		}
	}

	cout << (double)nAccept/(double)N << endl;

	time(&finish);

	cout << "Elapsed time is " << difftime(finish,start) << " seconds" << endl;


	return 0;
}