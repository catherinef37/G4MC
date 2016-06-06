// HRSRand.cpp: implementation of the HRSRand class.
//
//////////////////////////////////////////////////////////////////////

#include "HRSRand.hh"
#include "TTimeStamp.h"   //from root, can provide ns

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

HRSRand::HRSRand(unsigned seed)
{
	//initial class variables
	//ensure the seed is different even though begin at the same time

	if (seed==0)
	{
		TTimeStamp pTS;
		Seed=unsigned(time(0))+pTS.GetNanoSec();
	}
	else Seed=seed;

	//printf("\n***HRSRand::HRSRand(): set the random seed to %d ***\n",Seed);
	srand(Seed);

	mSeedLR=rand();
	mSeedDLR=rand();
}

HRSRand::~HRSRand()
{
}

//boxmuller gauss number generator
//input: mean m, standard deviation s 
double HRSRand::fGaus(double m, double s)	
{				        
	if(s==0.0) return m;

	double x1, x2, w, y1;
	static double y2;
	static int use_last = 0;

	if (use_last)		        /* use value from previous call */
	{
		y1 = y2;
		use_last = 0;
	}
	else
	{
		do {
			x1 = 2.0 * double(rand())/double(RAND_MAX) - 1.0;
			x2 = 2.0 * double(rand())/double(RAND_MAX) - 1.0;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );

		w = sqrt( (-2.0 * log( w ) ) / w );
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}

	return( m + y1 * s );
}

//return an unsigned form o to RAND_MAX
unsigned HRSRand::iRand()
{
	return rand();
}

unsigned HRSRand::iRand(unsigned low, unsigned high)
{
	unsigned range = high - low ;
	double temp=double(rand())/double(RAND_MAX);
	return low + unsigned(range * temp);
}

//return a random value from 0 to 1
double HRSRand::fRand()
{
	return double(rand())/double(RAND_MAX);
}

double HRSRand::fRand(double low, double high)
{
	double range = high - low ;
	double temp=double(rand())/double(RAND_MAX);
	return low + range * temp;
}


//return an unsigned from 0 to m,
//which following a linear distribution a*x+c
unsigned HRSRand::LinearRand(double a,unsigned c,unsigned m)
{
	if (m<=0)  m=RAND_MAX;
	unsigned x=mSeedLR;
	if (x==0) x=rand();
	mSeedLR=unsigned(a*x+c)%m;	
	return mSeedLR;
}

//return an unsigned from low to high
//which following a linear distribution a*x+c
unsigned HRSRand::LinearRand(double a,unsigned c,unsigned low,unsigned high)
{
	unsigned m=high-low;
	return low+LinearRand(a,c,m);
}


//generate a non-uniform x using the constrain(distribution) of f(x)
double HRSRand::VNReject(double xlow, double xhigh, double ylow,double yhigh,func f)
{
	double x,y;
	int counter=0;
	while(counter<100000)
	{
		x=this->fRand(xlow,xhigh);
		y=this->fRand(ylow,yhigh);
		if(f(x)>y) return x;
		counter++;
	};
	
	printf("VNReject() can not generate a value following the given distribution. Quit after 100000 trial...\n");
	return -counter;
}

//generate a non-uniform x using the constrain(distribution) of f(x,npar,par)
//this function f should be in this form double f(double,int,double *)
double HRSRand::VNReject(double xlow, double xhigh, double ylow,double yhigh,double *par,func2 f)
{
	double x,y;
	int counter=0;
	while(counter<100000)
	{
		x=this->fRand(xlow,xhigh);
		y=this->fRand(ylow,yhigh);
		if(f(x,par)>y) return x;
		counter++;
	};
	
	printf("VNReject() can not generate a value following the given distribution. Quit after 100000 trial...\n");
	return -counter;
}

//generate a sequence of x[i] and y[i] using the constrain f(x)
void HRSRand::VNRejectSeq(double xlow, double xhigh, double ylow,double yhigh,func f,
						  unsigned n, double xo[], double yo[])
{
	//input:
	//	n:	  the number of points
	//  x[] y[]: output x and y array 
	double x,y;
	unsigned counter=0;
	while(counter<n)
	{
		x=this->fRand(xlow,xhigh);
		y=this->fRand(ylow,yhigh);
		if(f(x)>y)
		{
			xo[counter]=x;
			yo[counter]=y;
			counter++;
		}
	}
	return;
}
