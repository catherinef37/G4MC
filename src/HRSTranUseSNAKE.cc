#include "HRSTranUseSNAKE.hh"
#include "GlobalDebuger.hh"

#include <iostream>
#include <iomanip>
#include <math.h>
#include "stdio.h"
#include "stdlib.h"


using namespace std;

//#define DEBUG_HRS_TRANSPORT 2
//#define DEBUG_HRS_RECONSTRUCTION 2

//input of forward:
//pVector_tg[0] = x_tg;
//pVector_tg[1] = theta_tg;
//pVector_tg[2] = y_tg;
//pVector_tg[3] = phi_tg;
//pVector_tg[4] = delta;
//output
//pVector_fp[0] = x_fp;
//pVector_fp[1] = theta_fp;
//pVector_fp[2] = y_fp;
//pVector_fp[3] = phi_fp;
//pVector_fp[4] = delta;		// = delta @ tg; no change
//input of backward:
//pVector_fp[0] = x_fp;
//pVector_fp[1] = theta_fp;
//pVector_fp[2] = y_fp;
//pVector_fp[3] = phi_fp;
//pVector_fp[4] = x_tg;
//output:
//pVector_rec_tg[0] = x_rec ;	// = x_tg; no change
//pVector_rec_tg[1] = theta_rec;
//pVector_rec_tg[2] = y_rec;
//pVector_rec_tg[3] = phi_rec;
//pVector_rec_tg[4] = delta;

namespace SNAKE
{
	//flat random number generator between [0,1)
	double fRand()
	{
		return double(rand())/double(RAND_MAX);
	}

	//boxmuller gauss number generator
	double fRandGaus(double m=0, double s=1.0)	/* normal random variate generator */
	{				        /* mean m, standard deviation s */
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
			do {;
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

	//return random number in [low,High) following a*x+c prob density
	double fLinearRand(double a,double c,double low,double high)
	{
		double m=double(RAND_MAX);
		double x,y;

		double ylow=a*low+c,yhigh=a*high+c;
		if(ylow<yhigh) 
		{
			double tmp=ylow;
			ylow=yhigh;
			yhigh=tmp;
		}

		do{
			x=low+(high-low)*((double)rand())/m;
			y=ylow+(yhigh-ylow)*((double)rand())/m;

			if((a*x+c)>=y) break;
		}while(true);
		return x;
	}


	bool SNAKEThruHRS(HRSTransport* pSnakeModel,double pX0_tg_m,
		const double *V5_tg,double* pV5_fp, double* pV5_rec)
	{

		//pV5_tg is the vector to be used in forward transportation functions 
		//in unit of meter, rad

		double pV5_tg[5]={V5_tg[0],V5_tg[1],V5_tg[2],V5_tg[3],V5_tg[4]};
		//Note that the 'anlge' defined in snake is tan(angle)
		pV5_tg[1]=tan(pV5_tg[1]);
		pV5_tg[3]=tan(pV5_tg[3]);

#if defined DEBUG_HRS_RECONSTRUCTION && (DEBUG_HRS_RECONSTRUCTION>=2)
		char str[1024];
		sprintf(str,"IN:  %8.4f %8.4f %8.4f %8.4f %8.4f",
			pV5_tg[0],pV5_tg[1],pV5_tg[2],pV5_tg[3],pV5_tg[4]);
		cout<<str<<endl;
#endif

#if defined DEBUG_HRS_RECONSTRUCTION && (DEBUG_HRS_RECONSTRUCTION>=3)
		//convert to meter and rad
		double pV5[5]={pV5_tg[0],pV5_tg[1],pV5_tg[2],pV5_tg[3],pV5_tg[4]};
		CMP_HRS_TRANSPORTATION(pV5_tg,pV5);
#endif


		bool bGoodParticle=false;


		///******************************************************************** //
		///*                        Forward Transportation                    * //
		///******************************************************************** //
		bGoodParticle=pSnakeModel->Forward(pV5_tg,pV5_fp); 
		pV5_fp[4]=pX0_tg_m;

		///********************************************************************* //
		///*                       Wire chamber smearing                       * //
		///********************************************************************* //
		bool bApplyVDCSmearing=true;
#if defined DEBUG_HRS_RECONSTRUCTION && (DEBUG_HRS_RECONSTRUCTION>=4)
		bApplyVDCSmearing=false;
#endif
		if(bApplyVDCSmearing)
		{
			// Resolutions 
			//double mHRSres = 2.9e-4;
			//double mBeamRes = 3.0e-5;
			double mWireChamberRes_x = 0.0013; //m;
			double mWireChamberRes_y = 0.0013; //m;
			double mWireChamberRes_theta = 0.0003; //rad;
			double mWireChamberRes_phi = 0.0003; //rad;
			pV5_fp[0] += mWireChamberRes_x * fRandGaus(0,1.0);
			pV5_fp[2] += mWireChamberRes_y * fRandGaus(0,1.0);
			pV5_fp[1] += mWireChamberRes_theta * fRandGaus(0,1.0);
			pV5_fp[3] += mWireChamberRes_phi * fRandGaus(0,1.0);
		}

		///************************************************************************* //
		///*                        Particle Reconstruction                        * //
		///************************************************************************* //

		if(bGoodParticle)
		{
			pSnakeModel->Backward(pV5_fp,pV5_rec); 

#if defined DEBUG_HRS_RECONSTRUCTION && (DEBUG_HRS_RECONSTRUCTION>=2)
			sprintf(str,"FP:  %8.4f %8.4f %8.4f %8.4f %8.4f \n",
				pV5_fp[0],pV5_fp[1],pV5_fp[2],pV5_fp[3],pV5_fp[4]);
			cout<<str;
			sprintf(str,"OUT: %8.4f %8.4f %8.4f %8.4f %8.4f \n",
				pV5_rec[0],pV5_rec[1],pV5_rec[2],pV5_rec[3],pV5_rec[4]);
			cout<<str<<endl;
#endif
		}

		//Note that the 'anlge' defined in snake is tan(angle)
		pV5_rec[1]=atan(pV5_rec[1]);
		pV5_rec[3]=atan(pV5_rec[3]);
		pV5_fp[1]=atan(pV5_fp[1]);
		pV5_fp[3]=atan(pV5_fp[3]);

		return bGoodParticle;
	}


	bool SNAKEForward(HRSTransport* pSnakeModel,const double *V5_tg,double* pV5_fp)
	{

		//pV5_tg is the vector to be used in forward transportation functions 
		//in unit of meter, rad

		double pV5_tg[5]={V5_tg[0],V5_tg[1],V5_tg[2],V5_tg[3],V5_tg[4]};
		//Note that the 'anlge' defined in snake is tan(angle)
		pV5_tg[1]=tan(pV5_tg[1]);
		pV5_tg[3]=tan(pV5_tg[3]);

#if defined DEBUG_HRS_RECONSTRUCTION && (DEBUG_HRS_RECONSTRUCTION>=2)
		char str[1024];
		sprintf(str,"IN:  %8.4f %8.4f %8.4f %8.4f %8.4f",
			pV5_tg[0],pV5_tg[1],pV5_tg[2],pV5_tg[3],pV5_tg[4]);
		cout<<str<<endl;
#endif

#if defined DEBUG_HRS_RECONSTRUCTION && (DEBUG_HRS_RECONSTRUCTION>=3)
		//convert to meter and rad
		double pV5[5]={pV5_tg[0],pV5_tg[1],pV5_tg[2],pV5_tg[3],pV5_tg[4]};
		CMP_HRS_TRANSPORTATION(pV5_tg,pV5);
#endif


		bool bGoodParticle=false;


		///******************************************************************** //
		///*                        Forward Transportation                    * //
		///******************************************************************** //
		bGoodParticle=pSnakeModel->Forward(pV5_tg,pV5_fp); 

		///********************************************************************* //
		///*                       Wire chamber smearing                       * //
		///********************************************************************* //
		bool bApplyVDCSmearing=true;
#if defined DEBUG_HRS_RECONSTRUCTION && (DEBUG_HRS_RECONSTRUCTION>=4)
		bApplyVDCSmearing=false;
#endif
		if(bApplyVDCSmearing)
		{
			// Resolutions 
			//double mHRSres = 2.9e-4;
			//double mBeamRes = 3.0e-5;
			double mWireChamberRes_x = 0.0013; //m;
			double mWireChamberRes_y = 0.0013; //m;
			double mWireChamberRes_theta = 0.0003; //rad;
			double mWireChamberRes_phi = 0.0003; //rad;
			pV5_fp[0] += mWireChamberRes_x * fRandGaus(0,1.0);
			pV5_fp[2] += mWireChamberRes_y * fRandGaus(0,1.0);
			pV5_fp[1] += mWireChamberRes_theta * fRandGaus(0,1.0);
			pV5_fp[3] += mWireChamberRes_phi * fRandGaus(0,1.0);
		}

		//Note that the 'anlge' defined in snake is tan(angle)
		pV5_fp[1]=atan(pV5_fp[1]);
		pV5_fp[3]=atan(pV5_fp[3]);

		return bGoodParticle;
	}

	void SNAKEBackward(HRSTransport* pSnakeModel,const double* V5_fp, double* pV5_rec)
	{

		//pV5_tg is the vector to be used in forward transportation functions 
		//in unit of meter, rad

		double pV5_fp[5]={V5_fp[0],V5_fp[1],V5_fp[2],V5_fp[3],V5_fp[4]};
		//Note that the 'anlge' defined in snake is tan(angle)
		pV5_fp[1]=tan(pV5_fp[1]);
		pV5_fp[3]=tan(pV5_fp[3]);

		pSnakeModel->Backward(pV5_fp,pV5_rec); 

#if defined DEBUG_HRS_RECONSTRUCTION && (DEBUG_HRS_RECONSTRUCTION>=2)
		sprintf(str,"FP:  %8.4f %8.4f %8.4f %8.4f %8.4f \n",
			pV5_fp[0],pV5_fp[1],pV5_fp[2],pV5_fp[3],pV5_fp[4]);
		cout<<str;
		sprintf(str,"OUT: %8.4f %8.4f %8.4f %8.4f %8.4f \n",
			pV5_rec[0],pV5_rec[1],pV5_rec[2],pV5_rec[3],pV5_rec[4]);
		cout<<str<<endl;
#endif

		//Note that the 'anlge' defined in snake is tan(angle)
		pV5_rec[1]=atan(pV5_rec[1]);
		pV5_rec[3]=atan(pV5_rec[3]);

	}


} //end of namespace 