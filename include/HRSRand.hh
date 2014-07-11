// HRSRand.h: unsignederface for the HRSRand class.
//
//////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#ifndef _HRSRAND_H_
#define _HRSRAND_H_

typedef double (*func)(double);
typedef double (*func2)(double,double*);

class HRSRand
{
public:
	HRSRand(unsigned seed=0);
	virtual ~HRSRand();

	//boxmuller gauss number generator
	//input: mean m, standard deviation s 
	static double   fGaus(double mean=0.0, double sigma=1.0);

	//return an unsigned from 0 to RAND_MAX
	static unsigned iRand();
	//return an unsigned from low to high
	static unsigned iRand(unsigned low, unsigned high);

	//return an double from 0 to 1
	static double   fRand();
	//return a double from low to high 
	static double   fRand(double low, double high);


	//return an unsigned from low to high
	//which following a linear distribution a*x+c
	unsigned LinearRand(double a,unsigned c,unsigned low,unsigned high);  
	//return an unsigned from 0 to m
	//which following a linear distribution a*x+c
	unsigned LinearRand(double a,unsigned c,unsigned m);   


	//Use the von Neumann rejection technique to generate an abnormal distribution f(x)
	//where f(x) is define as double f(double)
	double VNReject(double xlow, double xhigh, double ylow,double yhigh,func f);


	//generate a non-uniform x using the constrain(distribution) of f(x,npar,par)
	//this function f should be in this form double f(double,double *)
	double VNReject(double xlow, double xhigh, double ylow,double yhigh,double *par,func2 f);

	//generate a sequence of x[i] and y[i] using the constrain f(x)
	void VNRejectSeq(double xlow, double xhigh, double ylow,double yhigh,func f,
		unsigned n, double x[], double y[]);


public:
	unsigned Seed;

private:
	unsigned mSeedLR,mSeedDLR;   //used by linear random model
};

#endif //_HRSRAND_H_
