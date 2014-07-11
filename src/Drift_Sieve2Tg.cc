//Drift particle in the field

#include "BField_Helm.hh"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TObject.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TVector3.h"

#include "time.h"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <ios>
#include <vector>

using namespace std;


#define USEG4 1
#ifdef USEG4
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4Field.hh"
#endif

//Use_NystromRK4 just cost 60% of the time that used by the normal RK4
//It can only be used in NO E field situation
#define Use_NystromRK4 1

//DRIFT_VERBOSE: level 1 just print first and last step
//level 2: level 1 + write into a txt file for each interval of printstep(about 1mm)
//level 3: level 1 + print each step on the screen
//#define DRIFT_VERBOSE 1

//#define DRIFT_BENCHMARK 1		//will bentch mark for Drift() and DriftPath()
//#define DRIFT_BENCHMARK2 1    //will benchmark for Comp2Traj()
#if defined DRIFT_BENCHMARK || defined DRIFT_BENCHMARK2
static double gTotalTime=0;
static int gTotalN=0;
static double gTotalTime1=0;
static int gTotalN1=0;
static double gTotalTime2=0;
static int gTotalN2=0;
#endif

#include "Drift_Sieve2Tg.hh"
#include "HRSTransform_TCSNHCS.hh"

namespace DriftSieve2Tg
{
	static const double printstep = 0.001;      //print one step every 1 mm
	static const double c = 2.99796458E+08;		//[m/s]
	static const double q_e = 1.602E-19;		//[coulombs]=[C]
	static const double GeV2Kg = 1.783E-27;		//1 GeV/c2 = 1.783E-27 kg

	static const int kNVar=6;                   // number of first-order equations 

	static BField_Helm *gBField_Helm=0;
	static double kTargetOffset=-0.87693;

	//This variable will be set be Caller
	static int gUseSeptumInDrift=0;

#ifdef USEG4
	static const G4Field *gField = 0;
#endif

	////////////////////////////////////////////
	void  Init()
	{
		gBField_Helm=BField_Helm::GetInstance(); 
		if(gBField_Helm) 
		{
			double pOrigin[3];
			gBField_Helm->GetOrigin(pOrigin); 
			kTargetOffset=pOrigin[2]*0.01; //converted to meter
		}
#ifdef USEG4
		gField = G4TransportationManager::GetTransportationManager()->GetFieldManager()->GetDetectorField();
#endif
	}

	void SetUseSeptumInDrift(int val)
	{
		gUseSeptumInDrift=val;
	}

	void GetBField(const double pos[], double b[])
	{ 
#ifdef USEG4
		if(!gUseSeptumInDrift)
		{
			//for No septum field situation
			//this part is fast, it cost 10 ms for each call to Drift()  
			double x[]={pos[0]*100,pos[1]*100,pos[2]*100,0};  //convert to cm
			gBField_Helm->GetBField(x,b);
		}
		else
		{		
			//use G4 field manager, all position must be in Hall system and in G4 units
			//only need to do this when VB locate behind septum, which can be tell
			//by VB's location, rotation, or both.
			//it cost 12 ms for each call to Drift()
			double tmpPoint[4]={pos[0]*m,pos[1]*m,pos[2]*m,0};  //convert to Geant4 unit
			double tmpField[6];
			gField->GetFieldValue(tmpPoint,tmpField);
			for(int i=0;i<3;i++) {b[i]=tmpField[i]/tesla;}  //G4 return field in tesla,
		}
#else
		double x[]={pos[0]*100,pos[1]*100,pos[2]*100,0};  //convert to cm
		gBField_Helm->GetBField(x,b);
#endif
	}

		
	double DistChord(double *StartPoint,double *MidPoint,double *EndPoint)
	{
		// Calculate the deviation of the trajectory: dist2chord

		double ax = EndPoint[0] - StartPoint[0];
		double ay = EndPoint[1] - StartPoint[1];
		double az = EndPoint[2] - StartPoint[2];
		double dx = MidPoint[0] - StartPoint[0];
		double dy = MidPoint[1] - StartPoint[1];
		double dz = MidPoint[2] - StartPoint[2];
		double d2 = ax * ax + ay * ay + az*az;

		if (d2 != 0.) 
		{
			double s = (ax * dx + ay * dy + az * dz) / d2;
			dx -= (s * ax);
			dy -= (s * ay);
			dz -= (s * az);
		}

		return sqrt(dx * dx + dy * dy + dz * dz);
	}

	////////////////////////////////////////////////////////////////////////////////////
	//By Jixie: currently we never use the err[]
	double NystromRK4( const double* dxdt, double dt, const double* x, double* xo, 
		 double m0, double q, double &dist2chord /*, double* err*/)
	{
		double S = dt;
		double S5 = dt*0.5;
		double S4 = dt*0.25;
		double S6 = dt/6.0;

		const double *v = &x[3];
		double velocity2=v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
		double gamma = 1.0/sqrt(1.0-(velocity2)/c/c);
		const double m = m0*gamma;
		double coff=q/m;

		double R[3] = { x[0], x[1], x[2] };
		double A[3] = { dxdt[0], dxdt[1], dxdt[2] };

		double pStartPoint[3],pMidPoint[3],pEndPoint[3];
		pStartPoint[0] = R[0];
		pStartPoint[1] = R[1];
		pStartPoint[2] = R[2];

		// Point 1
		double K1[3] = { dxdt[3], dxdt[4], dxdt[5] };

		// Point 2
		double p[3] = { R[0]+S5*(A[0]+S4*K1[0]),
			R[1]+S5*(A[1]+S4*K1[1]),
			R[2]+S5*(A[2]+S4*K1[2])};

		double b[3];
		GetBField(p, b);

		double A2[3] = { A[0]+S5*K1[0], A[1]+S5*K1[1], A[2]+S5*K1[2] };
		double K2[3] = { (A2[1]*b[2]-A2[2]*b[1])*coff,
			(A2[2]*b[0]-A2[0]*b[2])*coff,
			(A2[0]*b[1]-A2[1]*b[0])*coff};

		pMidPoint[0] = p[0];
		pMidPoint[1] = p[1];
		pMidPoint[2] = p[2];

		// Point 3
		double A3[3] = { A[0]+S5*K2[0], A[1]+S5*K2[1], A[2]+S5*K2[2] };
		double K3[3] = { (A3[1]*b[2]-A3[2]*b[1])*coff,
			(A3[2]*b[0]-A3[0]*b[2])*coff,
			(A3[0]*b[1]-A3[1]*b[0])*coff};

		// Point 4
		p[0] = R[0]+S*(A[0]+S5*K3[0]);
		p[1] = R[1]+S*(A[1]+S5*K3[1]);
		p[2] = R[2]+S*(A[2]+S5*K3[2]);

		GetBField(p, b);

		double A4[3] = { A[0]+S*K3[0], A[1]+S*K3[1], A[2]+S*K3[2] };
		double K4[3] = { (A4[1]*b[2]-A4[2]*b[1])*coff,
			(A4[2]*b[0]-A4[0]*b[2])*coff,
			(A4[0]*b[1]-A4[1]*b[0])*coff};

		// New position
		xo[0] = R[0]+S*(A[0]+S6*(K1[0]+K2[0]+K3[0]));
		xo[1] = R[1]+S*(A[1]+S6*(K1[1]+K2[1]+K3[1]));
		xo[2] = R[2]+S*(A[2]+S6*(K1[2]+K2[2]+K3[2]));

		pEndPoint[0] = xo[0];
		pEndPoint[1] = xo[1];
		pEndPoint[2] = xo[2];

		// New direction
		xo[3] = A[0]+S6*(K1[0]+K4[0]+2.*(K2[0]+K3[0]));
		xo[4] = A[1]+S6*(K1[1]+K4[1]+2.*(K2[1]+K3[1]));
		xo[5] = A[2]+S6*(K1[2]+K4[2]+2.*(K2[2]+K3[2]));
 
		//distance to the chord, will be used to adjust step length
		dist2chord = DistChord( pStartPoint, pMidPoint, pEndPoint);
	
		//we do not use it now, comment it out for speed
		// Errors
		//err[3] = S*fabs(K1[0]-K2[0]-K3[0]+K4[0]);
		//err[4] = S*fabs(K1[1]-K2[1]-K3[1]+K4[1]);
		//err[5] = S*fabs(K1[2]-K2[2]-K3[2]+K4[2]);
		//err[0] = S*err[3];
		//err[1] = S*err[4];
		//err[2] = S*err[5];

		//do not do this if E field is on
		//do normalization such that the magnitue do not change
		double normF = sqrt(velocity2/(xo[3]*xo[3]+xo[4]*xo[4]+xo[5]*xo[5]));
		xo[3]*=normF;
		xo[4]*=normF;
		xo[5]*=normF;

		double tl2 = pow(xo[0]-x[0],2.0) + pow(xo[1]-x[1],2.0) + pow(xo[2]-x[2],2.0);
		return sqrt(tl2);
	}

	////////////////////////////////////////////////////////////////////////////////////
	double rk4(void(dnx)(double,const double [],double [],double,double), 
		double ti, double tf, double xi[], double xf[],double m0, double q)
	{
		const int n=kNVar;
		double h, t, x[n], dx[n];
		double k1[n],k2[n],k3[n],k4[n];
		int j;

		h = tf-ti;
		t = ti;
		//k1
		dnx(t, xi, dx, m0, q);
		for (j=0; j<n; j++)
		{
			k1[j] = h*dx[j];
			x[j]  = xi[j] + k1[j]/2.0;  
		}      
		//k2
		dnx(t+h/2.0, x, dx, m0, q);
		for (j=0; j<n; j++)
		{
			k2[j] = h*dx[j];
			x[j]  = xi[j] + k2[j]/2.0;  
		}
		//k3
		dnx(t+h/2.0, x, dx, m0, q);
		for (j=0; j<n; j++)
		{
			k3[j] = h*dx[j];
			x[j]  = xi[j] + k3[j];  
		}      
		//k4 and result      
		dnx(t+h, x, dx, m0, q);
		for (j=0; j<n; j++)
		{
			k4[j] = h*dx[j];
			xf[j] = xi[j] + k1[j]/6.0+k2[j]/3.0+k3[j]/3.0+k4[j]/6.0;
		}     

		//do not do this if E field exist
		//do normalization such that the magnitue do not change
		double velocity2 = xi[3]*xi[3]+xi[4]*xi[4]+xi[5]*xi[5];
		double normF = sqrt(velocity2/(xf[3]*xf[3]+xf[4]*xf[4]+xf[5]*xf[5]));
		xf[3]*=normF;
		xf[4]*=normF;
		xf[5]*=normF;

		double tl2 = pow(xf[0]-xi[0],2.0) + pow(xf[1]-xi[1],2.0) + pow(xf[2]-xi[2],2.0);
		return sqrt(tl2);
	}


	//For a system of n first-order ODEs
	//x [] array - x values
	//dx[] array - dx/dt values
	void dxdt(double t, const double x[], double dx[], double m0, double q)
	{
		double b[3];
		GetBField(x,b);
		const double *v = &x[3];
		const double gamma = 1.0/sqrt(1.0-(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])/c/c);
		const double m = m0*gamma;
		const double coff=q/m;	//divided is much slower than multiply, so calculate it once
		/* first order */
		dx[0] = v[0];
		dx[1] = v[1];
		dx[2] = v[2];
		/* second order */
		dx[3] = (v[1]*b[2]-b[1]*v[2])*coff;
		dx[4] = (v[2]*b[0]-b[2]*v[0])*coff;
		dx[5] = (v[0]*b[1]-b[0]*v[1])*coff;
	}    

	/////////////////////////////////////////////////////////////////////////////
	//return negative if error
	//This routine is used to drift particle in only B field, will stop once z across z_tr_limit 
	//or track_length>tracklengthlimit, which comes earlier 
	//definition of backward and forward: If Pz<0 go backward, otherwise forward
	int Drift2Ztr_Fast(TVector3& VX, TVector3& VP, double angle_rad, double m_gev, double q,
		double z_tr_limit, double tracklengthlimit, double steplength, double ztarget)
	{
#ifdef DRIFT_BENCHMARK 
		clock_t start = clock();
#endif
		if(!gBField_Helm) Init();

		const int n=6;                  // number of first-order equations 
		double dt=0.0, dl=0.0;
		double ti = 0.0, xi[n], xf[n];
		double dnx[n];
		int i;
		double tracklength=0.0, beta, ptot;

		double E0=sqrt(m_gev*m_gev+VP.Mag2());
		TVector3 VV = VP * (c/E0);		//in [m/s]
		double m0 = m_gev*GeV2Kg;		//turn into Kg [Kg]
		q *= q_e;						//turn in coulombs [C]

		double V0=VV.Mag();				//in [m/s]


		////////////////////////////////////////////////
		// initial information 
		xi[0] = VX.x();         // initial position in x (m)
		xi[1] = VX.y();         // initial position in y (m)
		xi[2] = VX.z();
		xi[3] = VV.x();			// initial speed in x direction (m.s)
		xi[4] = VV.y();			// initial speed in y direction (m/s)
		xi[5] = VV.z();  	    // initial speed in z direction (m/s)
		for (i=0; i<n; i++)  xf[i] = xi[i]; 

		// step size for integration (s)
		//in order to get accurate result, I set the time interval limit as the time to travel 1 mm
		const double kStepLimitHigh = 1.0E-4;
		const double kStepLimitLow = (steplength > 1.0E-6) ? 1.0E-6 : steplength;
		const double kDtLimit = kStepLimitHigh/V0;

		const double kDist2ChordHigh = 1.0E-6;
		const double kDist2ChordLow  = 1.0E-7;
		double dist2chord = (kDist2ChordHigh + kDist2ChordLow) / 2.0;

		dt = steplength/V0;
		// end of initial information 
		////////////////////////////////////////////////

		double sinHRS=sin(angle_rad);
		double cosHRS=cos(angle_rad);
		double new_z_tr = xi[0]*sinHRS+(xi[2]-ztarget)*cosHRS;	
		bool bIsBackward=(VP.z()<0.0)?true:false;

		if( bIsBackward && new_z_tr<z_tr_limit ) 
		{
			cout<<"\n***Error: Drift2Ztr_Fast() try to drift backward, but initial z_tr("
				<<new_z_tr<<") < z_tr_limit("<<z_tr_limit<<")\n\n";
			return -1;
		}
		else if( !bIsBackward && new_z_tr>z_tr_limit ) 
		{
			cout<<"\n***Error: Drift2Ztr_Fast() try to drift forward, but initial z_tr("
				<<new_z_tr<<") > z_tr_limit("<<z_tr_limit<<")\n\n";
			return -1;
		}

#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=1)
		cout<<"\nDrift2Ztr_Fast(): Swing the "<<((q>0)?"positron":"electron")
			<<(bIsBackward?" backward":" forward") <<" with Z_tr_limit="<<z_tr_limit
			<<" Ztarget="<<ztarget<<" and Tracklengthlimit="<<tracklengthlimit<<endl;
#endif

#ifdef DRIFT_VERBOSE
		printf("START: x        y        z       Px       Py       Pz     Ptot   Theta     Phi\n");
		printf("%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %7.2f\n",VX.x(),VX.y(),VX.z(),
			VP.x(),VP.y(),VP.z(),VP.Mag(),VP.Theta()*180./3.14159, VP.Phi()*180./3.14159); 
#endif

#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=2)
		int count=0,printinterval=(printstep>steplength)?int(printstep/steplength):1;
		const double cm=0.01;
		const double deg=acos(-1.0)/180.;

		// output: file and formats 
#if (DRIFT_VERBOSE>=3)
		ostream &file=cout;
#else
		ofstream file;
		file.open ("rk4_drift2ztr_fast.txt",fstream::app);		// write results to this file
#endif

		// output format 
		file.precision(5);
		//file.setf(ios::scientific | ios::showpoint);
		//cout.precision(6);
		//cout.setf(ios::fixed | ios::showpoint);
		file<< setw(12) << "t(ns)"  << setw(12) << "x(cm)" << setw(12) << "y(cm)"  << setw(12) << "z(cm)"
			<< setw(12) << "Px(GeV)" << setw(12) << "Py(GeV)" << setw(12) << "Pz(GeV)" << setw(12) 
			<< "Ptot(GeV)"<< setw(12) << " TL(cm)"<< setw(12) << " Theta(deg)"<< setw(12) << " Phi(deg)"
			<<endl; 
#endif
		// integration of ODE 
		// drift particle in various step length till it is overshoot
		// then reduce the step length for approaching
		while (tracklength < tracklengthlimit)
		{		
			/////////////////////////////////////
			if ( (dist2chord > kDist2ChordHigh) && (dt > kDtLimit) ) 
			{
				//rerun this step with a different step length
				dt *= 0.5;
			} 
			else 
			{
				//place trigger to stop integration
				//stop drift when z_tr across z_tr_limit
				new_z_tr = xf[0]*sinHRS+(xf[2]-ztarget)*cosHRS;
				//if( (bIsBackward && new_z_tr<=z_tr_limit) || (!bIsBackward && new_z_tr>=z_tr_limit) )
				if( bIsBackward ^ (new_z_tr>z_tr_limit) )
				{
#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=2)			
					cout<<"Overshoot at end plane: (z_tr_limit="<<z_tr_limit
						<<"), current z_tr="<<new_z_tr<<"\n";
#endif
					break;
				}

				// prepare for the next step
				for (i=0; i<n; i++)  xi[i] = xf[i]; 
				tracklength += dl;
				ti += dt;

				if (dist2chord < kDist2ChordLow) dt *= 2.0;

#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=2)
				if(!(count++%printinterval))
				{
					beta = sqrt(xi[3]*xi[3]+xi[4]*xi[4]+xi[5]*xi[5])/c;
					ptot = beta * ( m_gev / sqrt(1.0 - beta*beta) ); //in [GeV/c]		
					VV.SetXYZ(xf[3],xf[4],xf[5]);
					file<< setw(12) << ti*1.0E+09 
						<< setw(12) << xi[0]/cm << setw(12) << xi[1]/cm << setw(12) << xi[2]/cm
						<< setw(12) << xi[3]/c*ptot << setw(12) << xi[4]/c*ptot << setw(12) << xi[5]/c*ptot 
						<< setw(12) << ptot << setw(12) << tracklength/cm 
						<< setw(12) << VV.Theta()/deg << setw(12) << VV.Phi()/deg << endl;
				}
#endif  //#ifdef DRIFT_VERBOSE

			}

			dxdt(ti,xi,dnx,m0,q);
			dl = NystromRK4(dnx, dt, xi, xf, m0, q, dist2chord);
		}

		//already overshoot, now do approching
		if (tracklength < tracklengthlimit)
		{
			double pre_z_tr = xi[0]*sinHRS+(xi[2]-ztarget)*cosHRS;
			double dist2limit = fabs(new_z_tr-z_tr_limit);
			while( dist2limit > kStepLimitLow  )
			{		
				new_z_tr = xf[0]*sinHRS+(xf[2]-ztarget)*cosHRS;
				//if( (bIsBackward && new_z_tr<=z_tr_limit) || (!bIsBackward && new_z_tr>=z_tr_limit) )
				if( bIsBackward ^ (new_z_tr>z_tr_limit) )
				{
					//dt *= 0.5;
					dt = 0.9*fabs(pre_z_tr-z_tr_limit)/V0;
#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=2)
					cout<<"Approching: dist2limit="<<dist2limit<<"  dl="<<dl
						<<"  steplength="<<steplength<<"  new_z_tr="<<new_z_tr
						<<"  NewStepLimit="<< dt*V0<<endl;
#endif
				} 
				else 
				{
					// prepare for the next step
					for (i=0; i<n; i++)  xi[i] = xf[i]; 
					tracklength += dl;
					ti += dt;	
					pre_z_tr = new_z_tr;

#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=2)
					if(!(count++%printinterval))
					{
						beta = sqrt(xi[3]*xi[3]+xi[4]*xi[4]+xi[5]*xi[5])/c;
						ptot = beta * ( m_gev / sqrt(1.0 - beta*beta) ); //in [GeV/c]		
						VV.SetXYZ(xf[3],xf[4],xf[5]);
						file<< setw(12) << ti*1.0E+09 
							<< setw(12) << xi[0]/cm << setw(12) << xi[1]/cm << setw(12) << xi[2]/cm
							<< setw(12) << xi[3]/c*ptot << setw(12) << xi[4]/c*ptot << setw(12) << xi[5]/c*ptot 
							<< setw(12) << ptot << setw(12) << tracklength/cm 
							<< setw(12) << VV.Theta()/deg << setw(12) << VV.Phi()/deg << endl;
					}
#endif  //#ifdef DRIFT_VERBOSE
				}

				dxdt(ti,xi,dnx,m0,q);
				dl = NystromRK4(dnx, dt, xi, xf, m0, q, dist2chord);
				dist2limit = fabs(new_z_tr-z_tr_limit);	
			}

#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=1)			
			cout<<"Stop at end plane: (z_tr_limit="<<z_tr_limit
				<<"), current z_tr="<<new_z_tr<<"\n";
#endif
		}

		////////////////////////////////////////////////////////////////

		if(tracklength<1.0E-8) 
		{
			cout<<"***DriftSieve2Tg::Drift2Ztr_Fast() warning: wrong argument, did not move at all ***\n";
			return -1;
		}

		//prepare for return
		VX.SetXYZ(xf[0],xf[1],xf[2]); 
		VP.SetXYZ(xf[3],xf[4],xf[5]);
		beta = sqrt(xf[3]*xf[3]+xf[4]*xf[4]+xf[5]*xf[5])/c;
		ptot = beta * ( m_gev / sqrt(1.0 - beta*beta) ); //in [GeV/c]
		VP.SetMag(ptot);

#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=2)
		//print last step, need this when printstep is larger than actuall step
		if(count%printinterval)
		{
			file<< setw(12) << ti*1.0E+09 
				<< setw(12) << xi[0]/cm << setw(12) << xi[1]/cm << setw(12) << xi[2]/cm
				<< setw(12) << xi[3]/c*ptot << setw(12) << xi[4]/c*ptot << setw(12) << xi[5]/c*ptot 
				<< setw(12) << ptot << setw(12) << tracklength/cm 
				<< setw(12) << VP.Theta()/deg  << setw(12) << VP.Phi()/deg<< endl;
			file<<endl;
		}

#if (DRIFT_VERBOSE==2) 
		file.close();
#endif

#endif

#ifdef DRIFT_VERBOSE
		printf("END:   x        y        z       Px       Py       Pz     Ptot   Theta     Phi\n");
		printf("%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %7.2f\n",VX.x(),VX.y(),VX.z(),
			VP.x(),VP.y(),VP.z(),VP.Mag(),VP.Theta()*180./3.14159, VP.Phi()*180./3.14159); 
#endif

#ifdef DRIFT_BENCHMARK 		
		//cout<<"For this OS, CLOCKS_PER_SEC="<<CLOCKS_PER_SEC<<endl;
		clock_t end = clock();
		gTotalTime+=(double)(end-start)*1000.0/(double)CLOCKS_PER_SEC;
		gTotalN++;
		if(DRIFT_BENCHMARK>=2)
		{
			printf("Drift2Ztr() Elapsed time: %6.3f milisecond\r",
				(double)(end-start)*1000.0/(double)CLOCKS_PER_SEC);
		}
		if(!(gTotalN%2000))
		{
			printf("Drift_SL2TG():Drift2Ztr() Average elapsed time: %6.3f milisecond\n",
				gTotalTime/double(gTotalN));
			printf("Drift_SL2TG():Drift2Ztr() Total elapsed time: %6.0f milisecond over %d calls.\n",
				gTotalTime,gTotalN);
		}
#endif
		return 0;
	}

	/////////////////////////////////////////////////////////////////////////////
	//return negative if error
	//This routine is used to drift particle in the field, will stop once z_tr across z_tr_limit 
	//or track_length>tracklengthlimit, which comes earlier 
	//definition of backward and forward: If Pz<0 go backward, otherwise forward
	int Drift2Ztr(TVector3& VX, TVector3& VP, double angle_rad, double m_gev, double q,
		double z_tr_limit, double tracklengthlimit, double steplength, double ztarget)
	{
#ifdef DRIFT_BENCHMARK 
		clock_t start = clock();
#endif
		if(!gBField_Helm) Init();

		const int n=6;                  // number of first-order equations 
		double ti=0.0, dt=0.0, dl=0.0;
		double xi[n], xf[n];
		double dnx[n];
		int i;
		double tracklength=0.0, beta, ptot;

		double E0=sqrt(m_gev*m_gev+VP.Mag2());
		TVector3 VV = VP * (c/E0);		//in [m/s]
		double m0 = m_gev*GeV2Kg;		//turn into Kg [Kg]
		q *= q_e;						//turn in coulombs [C]

		double V0=VV.Mag();				//in [m/s]


		////////////////////////////////////////////////
		// initial information 
		ti = 0.0;               // initial value for variable t
		xi[0] = VX.x();         // initial position in x (m)
		xi[1] = VX.y();         // initial position in y (m)
		xi[2] = VX.z();
		xi[3] = VV.x();			// initial speed in x direction (m.s)
		xi[4] = VV.y();			// initial speed in y direction (m/s)
		xi[5] = VV.z();  	    // initial speed in z direction (m/s)
		for (i=0; i<n; i++)  xf[i] = xi[i]; 

		// step size for integration (s)
		//in order to get accurate result, I set the time interval limit as the time to travel 1 mm
		const double kStepLimitHigh = 1.0E-4;
		const double kStepLimitLow = (steplength > 1.0E-6) ? 1.0E-6 : steplength;
		const double kDtLimit = kStepLimitHigh/V0;

		const double kDist2ChordHigh = 1.0E-6;
		const double kDist2ChordLow  = 1.0E-7;
		double dist2chord = (kDist2ChordHigh + kDist2ChordLow) / 2.0;

		dt = steplength/V0;
		// end of initial information 
		////////////////////////////////////////////////

		double sinHRS=sin(angle_rad);
		double cosHRS=cos(angle_rad);
		double new_z_tr = xi[0]*sinHRS+(xi[2]-ztarget)*cosHRS;	
		bool bIsBackward=(VP.z()<0.0)?true:false;

		if( bIsBackward && new_z_tr<z_tr_limit ) 
		{
			cout<<"\n***Error: Drift2Ztr() try to drift backward, but initial z_tr("
				<<new_z_tr<<") < z_tr_limit("<<z_tr_limit<<")\n\n";
			return -1;
		}
		else if( !bIsBackward && new_z_tr>z_tr_limit ) 
		{
			cout<<"\n***Error: Drift2Ztr() try to drift forward, but initial z_tr("
				<<new_z_tr<<") > z_tr_limit("<<z_tr_limit<<")\n\n";
			return -1;
		}

#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=1)
		cout<<"\nDrift2Ztr(): Swing the "<<((q>0)?"positron":"electron")
			<<(bIsBackward?" backward":" forward") <<" with Z_tr_limit="<<z_tr_limit
			<<" Ztarget="<<ztarget<<" and Tracklengthlimit="<<tracklengthlimit<<endl;
#endif

#ifdef DRIFT_VERBOSE
		printf("START: x        y        z       Px       Py       Pz     Ptot   Theta     Phi\n");
		printf("%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %7.2f\n",VX.x(),VX.y(),VX.z(),
			VP.x(),VP.y(),VP.z(),VP.Mag(),VP.Theta()*180./3.14159, VP.Phi()*180./3.14159); 
#endif

#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=2)
		int count=0,printinterval=(printstep>steplength)?int(printstep/steplength):1;
		const double cm=0.01;
		const double deg=acos(-1.0)/180.;

		// output: file and formats 
#if (DRIFT_VERBOSE>=3)
		ostream &file=cout;
#else
		ofstream file;
		file.open ("rk4_drift2ztr.txt",fstream::app);		// write results to this file
#endif

		// output format 
		file.precision(5);
		//file.setf(ios::scientific | ios::showpoint);
		//cout.precision(6);
		//cout.setf(ios::fixed | ios::showpoint);
		file<< setw(12) << "t(ns)"  << setw(12) << "x(cm)" << setw(12) << "y(cm)"  << setw(12) << "z(cm)"
			<< setw(12) << "Px(GeV)" << setw(12) << "Py(GeV)" << setw(12) << "Pz(GeV)" << setw(12) 
			<< "Ptot(GeV)"<< setw(12) << " TL(cm)"<< setw(12) << " Theta(deg)"<< setw(12) << " Phi(deg)"
			<<endl; 
#endif

		// integration of ODE 
		// drift particle in various step length till it is overshoot
		// then reduce the step length for approaching
		while (tracklength < tracklengthlimit)
		{		
			/////////////////////////////////////
			if ( (dist2chord > kDist2ChordHigh) && (dt > kDtLimit) ) 
			{
				//rerun this step with a different step length
				dt *= 0.5;
			} 
			else 
			{
				//place trigger to stop integration
				//stop drift when z_tr across z_tr_limit
				new_z_tr = xf[0]*sinHRS+(xf[2]-ztarget)*cosHRS;
				//if( (bIsBackward && new_z_tr<=z_tr_limit) || (!bIsBackward && new_z_tr>=z_tr_limit) )
				if( bIsBackward ^ (new_z_tr>z_tr_limit) )
				{
#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=2)			
					cout<<"Overshoot at end plane: (z_tr_limit="<<z_tr_limit
						<<"), current z_tr="<<new_z_tr<<"\n";
#endif
					break;
				}

				// prepare for the next step
				for (i=0; i<n; i++)  xi[i] = xf[i]; 
				tracklength += dl;
				ti += dt;

				if (dist2chord < kDist2ChordLow) dt *= 2.0;

#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=2)
				if(!(count++%printinterval))
				{
					beta = sqrt(xi[3]*xi[3]+xi[4]*xi[4]+xi[5]*xi[5])/c;
					ptot = beta * ( m_gev / sqrt(1.0 - beta*beta) ); //in [GeV/c]		
					VV.SetXYZ(xf[3],xf[4],xf[5]);
					file<< setw(12) << ti*1.0E+09 
						<< setw(12) << xi[0]/cm << setw(12) << xi[1]/cm << setw(12) << xi[2]/cm
						<< setw(12) << xi[3]/c*ptot << setw(12) << xi[4]/c*ptot << setw(12) << xi[5]/c*ptot 
						<< setw(12) << ptot << setw(12) << tracklength/cm 
						<< setw(12) << VV.Theta()/deg << setw(12) << VV.Phi()/deg << endl;
				}
#endif  //#ifdef DRIFT_VERBOSE

			}

#ifdef Use_NystromRK4
			dxdt(ti,xi,dnx,m0,q);
			//=================================//
			dl = NystromRK4(dnx, dt, xi, xf, m0, q, dist2chord);
			//=================================//
#else
			//=================================//
			dl = rk4(dxdt, ti, ti+dt, xi, xf, m0, q);
			//=================================//
#endif

		}

		//already overshoot, now do approching
		if (tracklength < tracklengthlimit)
		{
			double pre_z_tr = xi[0]*sinHRS+(xi[2]-ztarget)*cosHRS;
			double dist2limit = fabs(new_z_tr-z_tr_limit);
			while( dist2limit > kStepLimitLow  )
			{		
				new_z_tr = xf[0]*sinHRS+(xf[2]-ztarget)*cosHRS;	
				//if( (bIsBackward && new_z_tr<=z_tr_limit) || (!bIsBackward && new_z_tr>=z_tr_limit) )
				if( bIsBackward ^ (new_z_tr>z_tr_limit) )
				{
					//dt *= 0.5;
					dt = 0.9*fabs(pre_z_tr-z_tr_limit)/V0;
#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=2)
					cout<<"Approching: dist2limit="<<dist2limit<<"  dl="<<dl
						<<"  steplength="<<steplength<<"  new_z_tr="<<new_z_tr
						<<"  NewStepLimit="<< dt*V0<<endl;
#endif
				} 
				else 
				{
					// prepare for the next step
					for (i=0; i<n; i++)  xi[i] = xf[i]; 
					tracklength += dl;
					ti += dt;	
					pre_z_tr = new_z_tr;

#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=2)
					if(!(count++%printinterval))
					{
						beta = sqrt(xi[3]*xi[3]+xi[4]*xi[4]+xi[5]*xi[5])/c;
						ptot = beta * ( m_gev / sqrt(1.0 - beta*beta) ); //in [GeV/c]		
						VV.SetXYZ(xf[3],xf[4],xf[5]);
						file<< setw(12) << ti*1.0E+09 
							<< setw(12) << xi[0]/cm << setw(12) << xi[1]/cm << setw(12) << xi[2]/cm
							<< setw(12) << xi[3]/c*ptot << setw(12) << xi[4]/c*ptot << setw(12) << xi[5]/c*ptot 
							<< setw(12) << ptot << setw(12) << tracklength/cm 
							<< setw(12) << VV.Theta()/deg << setw(12) << VV.Phi()/deg << endl;
					}
#endif  //#ifdef DRIFT_VERBOSE
				}

#ifdef Use_NystromRK4
				double dnx[n];
				dxdt(ti,xi,dnx,m0,q);
				//=================================//
				dl = NystromRK4(dnx, dt, xi, xf, m0, q, dist2chord);
				//=================================//
#else
				tf = ti + dt;
				//=================================//
				dl = rk4(dxdt, ti, tf, xi, xf, m0, q);
				//=================================//
#endif
				dist2limit = fabs(new_z_tr-z_tr_limit);	
			}//end of while

#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=1)			
			cout<<"Stop at end plane: (z_tr_limit="<<z_tr_limit
				<<"), current z_tr="<<new_z_tr<<"\n";
#endif
		}

		////////////////////////////////////////////////////////////////
		if(tracklength<1.0E-8) 
		{
			cout<<"***DriftSieve2Tg::Drift2Ztr() warning: wrong input argument, did not move at all ***\n";			
			return -1;
		}

		//prepare for return
		VX.SetXYZ(xf[0],xf[1],xf[2]); 
		VP.SetXYZ(xf[3],xf[4],xf[5]);
		beta = sqrt(xf[3]*xf[3]+xf[4]*xf[4]+xf[5]*xf[5])/c;
		ptot = beta * ( m_gev / sqrt(1.0 - beta*beta) ); //in [GeV/c]
		VP.SetMag(ptot);

#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=2)
		//print last step, need this when printstep is larger than actuall step
		if(count%printinterval)
		{
			file<< setw(12) << ti*1.0E+09 
				<< setw(12) << xi[0]/cm << setw(12) << xi[1]/cm << setw(12) << xi[2]/cm
				<< setw(12) << xi[3]/c*ptot << setw(12) << xi[4]/c*ptot << setw(12) << xi[5]/c*ptot 
				<< setw(12) << ptot << setw(12) << tracklength/cm 
				<< setw(12) << VP.Theta()/deg  << setw(12) << VP.Phi()/deg<< endl;
			file<<endl;
		}

#if (DRIFT_VERBOSE==2) 
		file.close();
#endif

#endif

#ifdef DRIFT_VERBOSE
		printf("END:   x        y        z       Px       Py       Pz     Ptot   Theta     Phi\n");
		printf("%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %7.2f\n",VX.x(),VX.y(),VX.z(),
			VP.x(),VP.y(),VP.z(),VP.Mag(),VP.Theta()*180./3.14159, VP.Phi()*180./3.14159); 
#endif

#ifdef DRIFT_BENCHMARK 		
		//cout<<"For this OS, CLOCKS_PER_SEC="<<CLOCKS_PER_SEC<<endl;
		clock_t end = clock();
		gTotalTime+=(double)(end-start)*1000.0/(double)CLOCKS_PER_SEC;
		gTotalN++;
		if(DRIFT_BENCHMARK>=2)
		{
			printf("Drift2Ztr() Elapsed time: %6.3f milisecond\r",
				(double)(end-start)*1000.0/(double)CLOCKS_PER_SEC);
		}
		if(!(gTotalN%2000))
		{
			printf("Drift_SL2TG():Drift2Ztr() Average elapsed time: %6.3f milisecond\n",
				gTotalTime/double(gTotalN));
			printf("Drift_SL2TG():Drift2Ztr() Total elapsed time: %6.0f milisecond over %d calls.\n",
				gTotalTime,gTotalN);
		}
#endif
		return 0;
	}

	//////////////////////////////////////////////////////////////////////////////////////

	//Drift in TCS,
	int Drift2Ztr(const double *V5In, double *V5Out, double& z_tr, double P0_gev, double angle_rad, double m_gev, double q,
		double z_tr_limit, double tracklengthlimit, double steplength, double ztarget)
	{
		//convert to HCS
		double x,y,z,ptot,theta,phi;
		Transform::X_TCS2HCS(V5In[0],V5In[2],z_tr,angle_rad,x,y,z);
		z += ztarget;
		Transform::P_TCS2HCS(V5In[1],V5In[3],angle_rad,theta,phi);
		ptot = (1.0+V5In[4])*P0_gev; 
		
		TVector3 VX(x,y,z);
		TVector3 VP;
		VP.SetMagThetaPhi(ptot,theta,phi); 

		bool bIsBackWard = (z_tr_limit < z_tr) ? true : false;
		
		if(bIsBackWard) {VP *= -1.0;q *= -1.0;}		
		int ret = Drift2Ztr(VX,VP,angle_rad,m_gev,q,z_tr_limit,tracklengthlimit,steplength,ztarget);	
		if(bIsBackWard) {VP *= -1.0;}

		//comvert back to TCS
		Transform::X_HCS2TCS(VX.X(),VX.Y(),VX.Z()-ztarget,angle_rad,V5Out[0],V5Out[2],z_tr);
		Transform::P_HCS2TCS(VP.Theta(),VP.Phi(),angle_rad,V5Out[1],V5Out[3]);
		//delta does not change, if it does, use the next line
		V5Out[4] = (VP.Mag()-P0_gev)/P0_gev;

		return ret;
	}


	//////////////////////////////////////////////////////////////////////////////////////
	//This routine is used to drift particle in the field, will stop once z across zlimit 
	//or track_length>tracklengthlimit 
	//definition of backward and forward: If Pz<0 go backward, otherwise forward
	int Drift2Z(TVector3& VX, TVector3& VP, double m_gev, double q, double zlimit,
		double tracklengthlimit, double steplength)
	{
#ifdef DRIFT_BENCHMARK 
		clock_t start = clock();
#endif
		if(!gBField_Helm) Init();

		const int n=6;                  // number of first-order equations 
		double ti=0.0, tf=0.0, dt=0.0, dl=0.0;
		double xi[n], xf[n];
		double dnx[n];
		int i;
		double tracklength=0.0, beta, ptot;

		double E0=sqrt(m_gev*m_gev+VP.Mag2());
		TVector3 VV = VP * (c/E0);		//in [m/s]
		double m0 = m_gev*GeV2Kg;		//turn into Kg [Kg]
		q *= q_e;						//turn in coulombs [C]

		double V0=VV.Mag();				//in [m/s]


		//////////////////////////////////////////////////
		// initial information 
		ti = tf = 0.0;               // initial value for variable t
		xi[0] = VX.x();         // initial position in x (m)
		xi[1] = VX.y();         // initial position in y (m)
		xi[2] = VX.z();
		xi[3] = VV.x();			// initial speed in x direction (m.s)
		xi[4] = VV.y();			// initial speed in y direction (m/s)
		xi[5] = VV.z();  	    // initial speed in z direction (m/s)
		for (i=0; i<n; i++)  xf[i] = xi[i]; 

		// step size for integration (s)
		//in order to get accurate result, I set the time interval limit as the time to travel 1 mm
		const double kStepLimitHigh = 1.0E-3;
		const double kStepLimitLow = (steplength > 1.0E-6) ? 1.0E-6 : steplength;
		const double kDtLimit = kStepLimitHigh/V0;

		const double kDist2ChordHigh = 1.0E-6;
		const double kDist2ChordLow  = 1.0E-7;
		double dist2chord = (kDist2ChordHigh + kDist2ChordLow) / 2.0;

		dt = steplength/V0;
		// end of initial information 
		//////////////////////////////////////////////////

		bool bIsBackward=(VP.z()<0.0)?true:false;

		if( bIsBackward && xi[2]<zlimit ) 
		{
			cout<<"\n***Error: Drift2Z() try to drift backward, but initial z("
				<<xi[2]<<") < zlimit("<<zlimit<<")\n\n";
			return -1;
		}
		else if( !bIsBackward && xi[2]>zlimit ) 
		{
			cout<<"\n***Error: Drift2Z() try to drift forward, but initial z("
				<<xi[2]<<") > zlimit("<<zlimit<<")\n\n";
			return -1;
		}

#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=1)
		cout<<"\nDrift2Z(): Swing the "<<((q>0)?"positron":"electron")
			<<(bIsBackward?" backward":" forward") <<" with Zlimit="<<zlimit
			<<" and Tracklengthlimit="<<tracklengthlimit<<endl;
#endif

#ifdef DRIFT_VERBOSE
		printf("START: x        y        z       Px       Py       Pz     Ptot   Theta     Phi\n");
		printf("%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %7.2f\n",VX.x(),VX.y(),VX.z(),
			VP.x(),VP.y(),VP.z(),VP.Mag(),VP.Theta()*180./3.14159, VP.Phi()*180./3.14159); 
#endif

		//////////////////////////////////////////////////


#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=2)
		int count=0,printinterval=(printstep>steplength)?int(printstep/steplength):1;
		const double cm=0.01;
		const double deg=acos(-1.0)/180.;

		// output: file and formats 
#if (DRIFT_VERBOSE>=3)
		ostream &file=cout;
#else
		ofstream file;
		file.open ("rk4_drift2z.txt",fstream::app);		// write results to this file
#endif

		// output format
		file.precision(5);
		//file.setf(ios::scientific | ios::showpoint);
		//cout.precision(6);
		//cout.setf(ios::fixed | ios::showpoint);
		file<< setw(12) << "t(ns)"  << setw(12) << "x(cm)" << setw(12) << "y(cm)"  << setw(12) << "z(cm)"
			<< setw(12) << "Px(GeV)" << setw(12) << "Py(GeV)" << setw(12) << "Pz(GeV)" << setw(12) 
			<< "Ptot(GeV)"<< setw(12) << " TL(cm)"<< setw(12) << " Theta(deg)"<< setw(12) << " Phi(deg)"
			<<endl; 
#endif

		// integration of ODE 
		// drift particle in various step length till it is overshoot
		// then reduce the step length for approaching
		while (tracklength < tracklengthlimit)
		{		
			/////////////////////////////////////
			if ( (dist2chord > kDist2ChordHigh) && (dt > kDtLimit) ) 
			{
				//rerun this step with a different step length
				dt *= 0.5;
			} 
			else 
			{
				//place trigger to stop integration
				//stop drift when z across zlimit
				//if( (bIsBackward && xf[2]<=zlimit) || (!bIsBackward && xf[2]>=zlimit) )
				if( bIsBackward ^ (xf[2]>zlimit) )
				{
	#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=2)
					cout<<"Overshoot at end plane: zlmit="<<zlimit<<" current z="<<xf[2]<<"\n";
	#endif
					break;
				}	

				// prepare for the next step
				for (i=0; i<n; i++)  xi[i] = xf[i]; 
				tracklength += dl;
				ti += dt;
				
				if (dist2chord < kDist2ChordLow) dt *= 2.0;

#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=2)
				if(!(count++%printinterval))
				{
					beta = sqrt(xi[3]*xi[3]+xi[4]*xi[4]+xi[5]*xi[5])/c;
					ptot = beta * ( m_gev / sqrt(1.0 - beta*beta) ); //in [GeV/c]		
					VV.SetXYZ(xf[3],xf[4],xf[5]);
					file<< setw(12) << ti*1.0E+09 
						<< setw(12) << xi[0]/cm << setw(12) << xi[1]/cm << setw(12) << xi[2]/cm
						<< setw(12) << xi[3]/c*ptot << setw(12) << xi[4]/c*ptot << setw(12) << xi[5]/c*ptot 
						<< setw(12) << ptot << setw(12) << tracklength/cm 
						<< setw(12) << VV.Theta()/deg << setw(12) << VV.Phi()/deg << endl;
				}
#endif  //#ifdef DRIFT_VERBOSE

			}

			/////////////////////////////////////

#ifdef Use_NystromRK4
			dxdt(ti,xi,dnx,m0,q);
			//=================================//
			dl = NystromRK4(dnx, dt, xi, xf, m0, q, dist2chord);
			//=================================//
#else
			//=================================//
			dl = rk4(dxdt, ti, ti+dt, xi, xf, m0, q);
			//=================================//
#endif

			/////////////////////////////////////

		}

		//already overshoot, now do approching
		if (tracklength < tracklengthlimit)
		{		
			double dist2limit = fabs(xf[2]-zlimit);
			while( dist2limit > kStepLimitLow  )
			{		
				if( bIsBackward ^ (xf[2]>zlimit) )
				{
					//dt *= 0.5;
					dt = 0.9*fabs(xi[2]-zlimit)/V0;
#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=2)
					cout<<"Approching: dist2limit="<<dist2limit<<"  dl="<<dl
						<<"  steplength="<<steplength<<"  new_z="<<xf[2]
					<<"  NewStepLimit="<< dt*V0<<endl;
#endif
				} 
				else 
				{
					// prepare for the next step
					for (i=0; i<n; i++)  xi[i] = xf[i]; 
					tracklength += dl;
					ti += dt;

					if (dist2chord < kDist2ChordLow) dt *= 2.0;

#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=2)
					if(!(count++%printinterval))
					{
						beta = sqrt(xi[3]*xi[3]+xi[4]*xi[4]+xi[5]*xi[5])/c;
						ptot = beta * ( m_gev / sqrt(1.0 - beta*beta) ); //in [GeV/c]		
						VV.SetXYZ(xf[3],xf[4],xf[5]);
						file<< setw(12) << ti*1.0E+09 
							<< setw(12) << xi[0]/cm << setw(12) << xi[1]/cm << setw(12) << xi[2]/cm
							<< setw(12) << xi[3]/c*ptot << setw(12) << xi[4]/c*ptot << setw(12) << xi[5]/c*ptot 
							<< setw(12) << ptot << setw(12) << tracklength/cm 
							<< setw(12) << VV.Theta()/deg << setw(12) << VV.Phi()/deg << endl;
					}
#endif  //#ifdef DRIFT_VERBOSE
				}

#ifdef Use_NystromRK4
				dxdt(ti,xi,dnx,m0,q);
				//=================================//
				dl = NystromRK4(dnx, dt, xi, xf, m0, q, dist2chord);
				//=================================//
#else
				//=================================//
				dl = rk4(dxdt, ti, ti+dt, xi, xf, m0, q);
				//=================================//
#endif
				dist2limit = fabs(xf[2]-zlimit);
			}//end of while
#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=1)			
			cout<<"Stop at end plane: (zlimit="<<zlimit
				<<"), current z="<<xf[2]<<"\n";
#endif
		}

		////////////////////////////////////////////////////////////////
		if(tracklength<1.0E-8) 
		{
			cout<<"***DriftSieve2Tg::Drift2Z() warning: wrong argument, did not move at all ***\n";
			return -1;
		}

		//prepare for return
		VX.SetXYZ(xf[0],xf[1],xf[2]); 
		VP.SetXYZ(xf[3],xf[4],xf[5]);
		beta = sqrt(xf[3]*xf[3]+xf[4]*xf[4]+xf[5]*xf[5])/c;
		ptot = beta * ( m_gev / sqrt(1.0 - beta*beta) ); //in [GeV/c]
		VP.SetMag(ptot);

#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=2)
		//print last step, need this when printstep is larger than actuall step
		if(count%printinterval)
		{
			file<< setw(12) << ti*1.0E+09 
				<< setw(12) << xi[0]/cm << setw(12) << xi[1]/cm << setw(12) << xi[2]/cm
				<< setw(12) << xi[3]/c*ptot << setw(12) << xi[4]/c*ptot << setw(12) << xi[5]/c*ptot 
				<< setw(12) << ptot << setw(12) << tracklength/cm 
				<< setw(12) << VP.Theta()/deg  << setw(12) << VP.Phi()/deg<< endl;
			file<<endl;
		}

#if (DRIFT_VERBOSE==2) 
		file.close();
#endif

#endif

#ifdef DRIFT_VERBOSE
		printf("END:   x        y        z       Px       Py       Pz     Ptot   Theta     Phi\n");
		printf("%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %7.2f\n",VX.x(),VX.y(),VX.z(),
			VP.x(),VP.y(),VP.z(),VP.Mag(),VP.Theta()*180./3.14159, VP.Phi()*180./3.14159); 
#endif

#ifdef DRIFT_BENCHMARK 		
		//cout<<"For this OS, CLOCKS_PER_SEC="<<CLOCKS_PER_SEC<<endl;
		clock_t end = clock();
		gTotalTime2+=(double)(end-start)*1000.0/(double)CLOCKS_PER_SEC;
		gTotalN2++;
		if(DRIFT_BENCHMARK>=2)
		{
			printf("Drift2Z() Elapsed time: %6.3f milisecond\r",
				(double)(end-start)*1000.0/(double)CLOCKS_PER_SEC);
		}
		if(!(gTotalN2%2000))
		{
			printf("Drift_SL2TG():Drift2Z() Average elapsed time: %6.3f milisecond\n",
				gTotalTime2/double(gTotalN2));
			printf("Drift_SL2TG():Drift2Z() Total elapsed time: %6.0f milisecond over %d calls.\n",
				gTotalTime2,gTotalN2);
		}
#endif
		return 0;
	}


	//////////////////////////////////////////////////////////////////////////////////////
	//This routine is used to drift particle in the field, will stop once z across zlimit 
	//or track_length>tracklengthlimit 
	//definition of backward and forward: If Pz<0 go backward, otherwise forward
	int DriftPath(TVector3& VX, TVector3& VP, double m_gev, double q, double zlimit,
		double tracklengthlimit, double steplength, vector <vector <double> > &traj)
	{
#ifdef DRIFT_BENCHMARK 
		clock_t start = clock();
#endif
		if(!gBField_Helm) Init();

		bool bIsBackward=(VP.z()<0.0)?true:false;

#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=1)
		cout<<"\nDriftPath(): Swing the "<<((q>0)?"positron":"electron")
			<<(bIsBackward?" backward":" forward") <<" with Zlimit="<<zlimit
			<<" and Tracklengthlimit="<<tracklengthlimit<<endl;
#endif

		const int n=6;                  // number of first-order equations 
		double ti=0.0, tf=0.0, dt=0.0, dl=0.0;
		double xi[n], xf[n];
		int i;
		double tracklength=0.0, beta, ptot;

		double E0=sqrt(m_gev*m_gev+VP.Mag2());
		TVector3 VV = VP * (c/E0);		//in [m/s]
		double m0 = m_gev*GeV2Kg;		//turn into Kg [Kg]
		q *= q_e;						//turn in coulombs [C]

		double V0=VV.Mag();				//in [m/s]

		//////////////////////////////////////////////
		// initial information 
		ti = tf = 0.0;               // initial value for variable t
		xi[0] = VX.x();         // initial position in x (m)
		xi[1] = VX.y();         // initial position in y (m)
		xi[2] = VX.z();
		xi[3] = VV.x();			// initial speed in x direction (m.s)
		xi[4] = VV.y();			// initial speed in y direction (m/s)
		xi[5] = VV.z();  	    // initial speed in z direction (m/s)
		for (i=0; i<n; i++)  xf[i] = xi[i]; 

		// step size for integration (s)
		//in order to get accurate result, I set the time interval limit as the time to travel 1 mm
		const double kStepLimitHigh = 1.0E-3;
		const double kStepLimitLow = (steplength > 1.0E-6) ? 1.0E-6 : steplength;
		const double kDtLimit = kStepLimitHigh/V0;

		const double kDist2ChordHigh = 1.0E-6;
		const double kDist2ChordLow  = 1.0E-7;
		double dist2chord = (kDist2ChordHigh + kDist2ChordLow) / 2.0;

		dt = steplength/V0;
		// end of initial information 
		//////////////////////////////////////////////

#ifdef DRIFT_VERBOSE
		printf("START: x        y        z       Px       Py       Pz     Ptot   Theta     Phi\n");
		printf("%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %7.2f\n",VX.x(),VX.y(),VX.z(),
			VP.x(),VP.y(),VP.z(),VP.Mag(),VP.Theta()*180./3.14159, VP.Phi()*180./3.14159); 
#endif

		//////////////////////////////////////////////////
		traj.clear(); 
		//record the initial step
		vector<double> point(6,0.0);
		//point.push_back(xi[0]);
		//point.push_back(xi[1]);
		//point.push_back(xi[2]);
		//point.push_back(xi[3]);
		//point.push_back(xi[4]);
		//point.push_back(xi[5]);
		//traj.push_back(point); 
		//////////////////////////////////////////////////


#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=2)
		int count=0,printinterval=(printstep>steplength)?int(printstep/steplength):1;
		const double cm=0.01;
		const double deg=acos(-1.0)/180.;

		// output: file and formats 
#if (DRIFT_VERBOSE>=3)
		ostream &file=cout;
#else
		ofstream file;
		file.open ("rk4_driftpath.txt",fstream::app);		// write results to this file
#endif

		// output format
		file.precision(5);
		//file.setf(ios::scientific | ios::showpoint);
		//cout.precision(6);
		//cout.setf(ios::fixed | ios::showpoint);
		file<< setw(12) << "t(ns)"  << setw(12) << "x(cm)" << setw(12) << "y(cm)"  << setw(12) << "z(cm)"
			<< setw(12) << "Px(GeV)" << setw(12) << "Py(GeV)" << setw(12) << "Pz(GeV)" << setw(12) 
			<< "Ptot(GeV)"<< setw(12) << " TL(cm)"<< setw(12) << " Theta(deg)"<< setw(12) << " Phi(deg)"
			<<endl; 
#endif

		// integration of ODE 
		while (tracklength < tracklengthlimit)
		{		
			//////////////////////////////////////////////////
			//adjust step length if needed
			if ( (dist2chord > kDist2ChordHigh) && (dt > kDtLimit) ) 
			{
				//rerun this step with a different step length
				dt *= 0.5;
			} 
			else 
			{
				// prepare for the next step
				ti = tf;
				for (i=0; i<n; i++)  xi[i] = xf[i]; 
				tracklength += dl;
			
#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=2)
				if(!(count++%printinterval))
				{
					beta = sqrt(xi[3]*xi[3]+xi[4]*xi[4]+xi[5]*xi[5])/c;
					ptot = beta * ( m_gev / sqrt(1.0 - beta*beta) ); //in [GeV/c]		
					VV.SetXYZ(xf[3],xf[4],xf[5]);
					file<< setw(12) << ti*1.0E+09 
						<< setw(12) << xi[0]/cm << setw(12) << xi[1]/cm << setw(12) << xi[2]/cm
						<< setw(12) << xi[3]/c*ptot << setw(12) << xi[4]/c*ptot << setw(12) << xi[5]/c*ptot 
						<< setw(12) << ptot << setw(12) << tracklength/cm 
						<< setw(12) << VV.Theta()/deg << setw(12) << VV.Phi()/deg << endl;
				}
#endif  //#ifdef DRIFT_VERBOSE

				//approching is defined as distance to limit is smaller than last step length
				//or smaller than given steplength or smaller than 1 mm
				double dist2limit = fabs(xi[2]-zlimit);
				if( dist2limit < kStepLimitLow || (dist2limit < dl && dl > kStepLimitLow) )
				{			
					//approching the limit
					double pApprochStepLimit = (dist2limit > 3.0*kStepLimitLow) ? 0.3*dist2limit : kStepLimitLow;
					dt = pApprochStepLimit/V0;
#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=2)
					cout<<"Approching: dist2limit="<<dist2limit<<"  dl="<<dl
						<<"  steplength="<<steplength<<"  z="<<xi[2]
						<<"  NewStepLimit="<< pApprochStepLimit<<endl;
#endif
				}
				else if (dist2chord < kDist2ChordLow) dt *= 2.0;

				//record this step, place here so it can also record the 1st step
				for(int kk=0;kk<6;kk++) point[kk]=xi[kk];
				traj.push_back(point); 

			}
			//////////////////////////////////////////////////

			/////////////////////////////////////
			//place trigger to stop integration
			//stop drift when z across zlimit

			//if( (bIsBackward && xi[2]<=zlimit) || (!bIsBackward && xi[2]>=zlimit) )
			if( bIsBackward ^ (xi[2]>zlimit) )
			{
#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=1)
				cout<<"Stop at zlmit="<<zlimit<<" current z="<<xi[2]<<"\n";
#endif
				break;
			}	
			/////////////////////////////////////

			tf = ti + dt;
#ifdef Use_NystromRK4
			double dnx[n];
			dxdt(ti,xi,dnx,m0,q);
			//=================================//
			dl = NystromRK4(dnx, dt, xi, xf, m0, q, dist2chord);
			//=================================//
#else
			//=================================//
			dl = rk4(dxdt, ti, tf, xi, xf, m0, q);
			//=================================//
#endif

			//////////////////////////////////////////////////

		}

		if(tracklength<1.0E-8) 
		{
			cout<<"***DriftSieve2Tg::DriftPath() warning: wrong argument, did not move at all ***\n";
			return -1;
		}

		//prepare for return
		VX.SetXYZ(xf[0],xf[1],xf[2]); 
		VP.SetXYZ(xf[3],xf[4],xf[5]);
		beta = sqrt(xf[3]*xf[3]+xf[4]*xf[4]+xf[5]*xf[5])/c;
		ptot = beta * ( m_gev / sqrt(1.0 - beta*beta) ); //in [GeV/c]
		VP.SetMag(ptot);

#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=2)
		//print last step, need this when printstep is larger than actuall step
		if(count%printinterval)
		{
			file<< setw(12) << ti*1.0E+09 
				<< setw(12) << xi[0]/cm << setw(12) << xi[1]/cm << setw(12) << xi[2]/cm
				<< setw(12) << xi[3]/c*ptot << setw(12) << xi[4]/c*ptot << setw(12) << xi[5]/c*ptot 
				<< setw(12) << ptot << setw(12) << tracklength/cm 
				<< setw(12) << VP.Theta()/deg  << setw(12) << VP.Phi()/deg<< endl;
			file<<endl;
		}
#if (DRIFT_VERBOSE==2) 
		file.close();
#endif

#endif

#ifdef DRIFT_VERBOSE
		printf("END:   x        y        z       Px       Py       Pz     Ptot   Theta     Phi\n");
		printf("%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %7.2f\n",VX.x(),VX.y(),VX.z(),
			VP.x(),VP.y(),VP.z(),VP.Mag(),VP.Theta()*180./3.14159, VP.Phi()*180./3.14159); 
#endif

#ifdef DRIFT_BENCHMARK 		
		//cout<<"For this OS, CLOCKS_PER_SEC="<<CLOCKS_PER_SEC<<endl;
		clock_t end = clock();
		gTotalTime1+=(double)(end-start)*1000.0/(double)CLOCKS_PER_SEC;
		gTotalN1++;
		if(DRIFT_BENCHMARK>=2)
		{
			printf("DriftPath() Elapsed time: %6.3f milisecond\r",
				(double)(end-start)*1000.0/(double)CLOCKS_PER_SEC);
		}
		if(!(gTotalN1%4000))
		{
			printf("Drift_SL2TG():DriftPath() Average elapsed time: %6.3f milisecond\n",
				gTotalTime1/double(gTotalN1));
			//printf("Drift_SL2TG():DriftPath() Total elapsed time: %6.0f milisecond over %d calls.\n",
			//	gTotalTime1,gTotalN1);
		}
#endif
		return 0;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////
	bool Cmp2Traj(vector < vector <double> > &trajrec, vector < vector <double> > &trajbpm, 
		size_t &indexI, size_t &indexJ, double &dist)
	{
#ifdef DRIFT_BENCHMARK2
		clock_t start = clock();
#endif
		size_t s_rec=trajrec.size();
		size_t s_bpm=trajbpm.size();

		bool found=false;
		dist=10000.0; indexI=10000;indexJ=10000;
		double tmpdist;
		for (size_t i=0;i<s_bpm;i++)
		{
			//take one point at bpm traj, compare it to 9 points at rec traj, 
			//record the most close distance and index at bpm traj
			for(int j=int(i)-4;j<=int(i)+4;j++)
			{
				if(j>=0 && j<int(s_rec)) 
				{
					//only compare those points that z are within 1mm
					if(fabs(trajrec[j][2]-trajbpm[i][2])<0.001) 
					{
						tmpdist = sqrt( pow(trajrec[j][0]-trajbpm[i][0],2.0) + 
							pow(trajrec[j][1]-trajbpm[i][1],2.0) +
							pow(trajrec[j][2]-trajbpm[i][2],2.0) );
						if(tmpdist<dist) 
						{
							found=true;
							dist = tmpdist; 
							indexI = i;
							indexJ = j;
						}
					}
				} //end if(j>=0)
			} //end for j
		}  //end for i

#ifdef DRIFT_BENCHMARK2 		
		//cout<<"For this OS, CLOCKS_PER_SEC="<<CLOCKS_PER_SEC<<endl;
		clock_t end = clock();
		gTotalTime+=(double)(end-start)*1000.0/(double)CLOCKS_PER_SEC;
		gTotalN++;
		if(DRIFT_BENCHMARK2>=2)
		{
			printf("Cmp2Traj() Elapsed time: %6.3f milisecond\r",
				(double)(end-start)*1000.0/(double)CLOCKS_PER_SEC);
		}
		if(!(gTotalN%2000))
		{
			printf("Drift_SL2TG():Cmp2Traj() Average elapsed time: %6.3f milisecond\n",
				gTotalTime/gTotalN);
			//printf("Drift_SL2TG():Cmp2Traj() Total elapsed time: %6.0f milisecond over %d calls.\n",
			//	gTotalTime,gTotalN);
		}
#endif
		//cout<<"Most closest distance("<<dist<<" m): trajbpm["<<indexI<<"] to trajrec["<<indexJ<<"]\n";
		return found;
	}


	/////////////////////////////////////////////////////////////////////////////
	//return negative if error
	//This routine is used to drift particle in the field, will stop once z_tr across z_tr_limit 
	//or track_length>tracklengthlimit, which comes earlier 
	//definition of backward and forward: If Pz<0 go backward, otherwise forward
	int Drift2Ztr_EqualStep(TVector3& VX, TVector3& VP, double angle_rad, double m_gev, double q,
		double z_tr_limit, double tracklengthlimit, double steplength, double ztarget)
	{
#ifdef DRIFT_BENCHMARK 
		clock_t start = clock();
#endif
		if(!gBField_Helm) Init();

		bool bIsBackward=(VP.z()<0.0)?true:false;

#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=1)
		cout<<"\nDrift2Ztr_EqualStep(): Swing the "<<((q>0)?"positron":"electron")
			<<(bIsBackward?" backward":" forward") <<" with z_tr_limit="<<z_tr_limit
			<<" Ztarget="<<ztarget<<" and Tracklengthlimit="<<tracklengthlimit<<endl;
#endif

		const int n=6;                  // number of first-order equations 
		double ti, tf, dt;
		double xi[n], xf[n];
		int i;
		double tracklength=0.0, beta, ptot;

		double new_z_tr;
		double sinHRS=sin(angle_rad);
		double cosHRS=cos(angle_rad);


		double E0=sqrt(m_gev*m_gev+VP.Mag2());
		TVector3 VV = VP * (c/E0);		//in [m/s]
		double m0 = m_gev*GeV2Kg;		//turn into Kg [Kg]
		q *= q_e;						//turn in coulombs [C]

		double V0=VV.Mag();				//in [m/s]


		////////////////////////////////////////////////
		// initial information 
		ti = tf = 0.0;               // initial value for variable t
		xi[0] = VX.x();         // initial position in x (m)
		xi[1] = VX.y();         // initial position in y (m)
		xi[2] = VX.z();
		xi[3] = VV.x();			// initial speed in x direction (m.s)
		xi[4] = VV.y();			// initial speed in y direction (m/s)
		xi[5] = VV.z();  	    // initial speed in z direction (m/s)

		// step size for integration (s)
		//in order to get accurate result, I set the time interval as the time to travel 100 um
		double dist2chord;
		dt = steplength/VV.Mag();
		// end of initial information 
		////////////////////////////////////////////////

#ifdef DRIFT_VERBOSE
		printf("START: x        y        z       Px       Py       Pz     Ptot   Theta     Phi\n");
		printf("%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %7.2f\n",VX.x(),VX.y(),VX.z(),
			VP.x(),VP.y(),VP.z(),VP.Mag(),VP.Theta()*180./3.14159, VP.Phi()*180./3.14159); 
#endif

#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=2)
		int count=0,printinterval=(printstep>steplength)?int(printstep/steplength):1;
		const double cm=0.01;
		const double deg=acos(-1.0)/180.;

		// output: file and formats 
#if (DRIFT_VERBOSE>=3)
		ostream &file=cout;
#else
		ofstream file;
		file.open ("rk4_drift2ztr.txt",fstream::app);		// write results to this file
#endif

		// output format 
		file.precision(5);
		//file.setf(ios::scientific | ios::showpoint);
		//cout.precision(6);
		//cout.setf(ios::fixed | ios::showpoint);
		file<< setw(12) << "t(ns)"  << setw(12) << "x(cm)" << setw(12) << "y(cm)"  << setw(12) << "z(cm)"
			<< setw(12) << "Px(GeV)" << setw(12) << "Py(GeV)" << setw(12) << "Pz(GeV)" << setw(12) 
			<< "Ptot(GeV)"<< setw(12) << " TL(cm)"<< setw(12) << " Theta(deg)"<< setw(12) << " Phi(deg)"
			<<endl; 
#endif

		// integration of ODE 
		while (tracklength < tracklengthlimit)
		{		
#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=2)
			if(!(count++%printinterval))
			{
				beta = sqrt(xi[3]*xi[3]+xi[4]*xi[4]+xi[5]*xi[5])/c;
				ptot = beta * ( m_gev / sqrt(1.0 - beta*beta) ); //in [GeV/c]		
				VV.SetXYZ(xf[3],xf[4],xf[5]);
				file<< setw(12) << ti*1.0E+09 
					<< setw(12) << xi[0]/cm << setw(12) << xi[1]/cm << setw(12) << xi[2]/cm
					<< setw(12) << xi[3]/c*ptot << setw(12) << xi[4]/c*ptot << setw(12) << xi[5]/c*ptot 
					<< setw(12) << ptot << setw(12) << tracklength/cm 
					<< setw(12) << VV.Theta()/deg << setw(12) << VV.Phi()/deg << endl;
			}
#endif  //#ifdef DRIFT_VERBOSE

			new_z_tr = xi[0]*sinHRS+(xi[2]-ztarget)*cosHRS;

			//stop drift when z_tr across z_tr_limit
			//if( (bIsBackward && new_z_tr<=z_tr_limit) || (!bIsBackward && new_z_tr>=z_tr_limit) )
			if( bIsBackward ^ (new_z_tr>z_tr_limit) )
			{
#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=1)			
				cout<<"stop at end plane(Z_tr="<<z_tr_limit
					<<"), current z_tr="<<new_z_tr<<"\n";
#endif
				break;
			}


			tf = ti + dt;
#ifdef Use_NystromRK4
			double dnx[n];
			dxdt(ti,xi,dnx,m0,q);
			//=================================//
			tracklength += NystromRK4(dnx, dt, xi, xf, m0, q, dist2chord);
			//=================================//
#else
			//=================================//
			tracklength += rk4(dxdt, ti, tf, xi, xf, m0, q);
			//=================================//
#endif


			//the particle will always gain energy from this routine 
			//need to reset the total momentum	
			double normalizefactor=V0/sqrt(xf[3]*xf[3]+xf[4]*xf[4]+xf[5]*xf[5]);		 
			xf[3]*=normalizefactor;xf[4]*=normalizefactor;xf[5]*=normalizefactor;


			// prepare for the next step
			ti = tf;
			for (i = 0; i<n; i++)  xi[i] = xf[i]; 
		}

		if(tracklength<1.0E-8) 
		{
			cout<<"***DriftSieve2Tg::Drift2Ztr_EqualStep() warning: wrong argument, did not move at all ***\n";
			return -1;
		}

		//prepare for return
		VX.SetXYZ(xf[0],xf[1],xf[2]); 
		VP.SetXYZ(xf[3],xf[4],xf[5]);
		beta = sqrt(xf[3]*xf[3]+xf[4]*xf[4]+xf[5]*xf[5])/c;
		ptot = beta * ( m_gev / sqrt(1.0 - beta*beta) ); //in [GeV/c]
		VP.SetMag(ptot);

#if defined DRIFT_VERBOSE && (DRIFT_VERBOSE>=2)
		//print last step, need this when printstep is larger than actuall step
		if(count%printinterval)
		{
			file<< setw(12) << ti*1.0E+09 
				<< setw(12) << xi[0]/cm << setw(12) << xi[1]/cm << setw(12) << xi[2]/cm
				<< setw(12) << xi[3]/c*ptot << setw(12) << xi[4]/c*ptot << setw(12) << xi[5]/c*ptot 
				<< setw(12) << ptot << setw(12) << tracklength/cm 
				<< setw(12) << VP.Theta()/deg  << setw(12) << VP.Phi()/deg<< endl;
			file<<endl;
		}

#if (DRIFT_VERBOSE==2) 
		file.close();
#endif

#endif

#ifdef DRIFT_VERBOSE
		printf("END:   x        y        z       Px       Py       Pz     Ptot   Theta     Phi\n");
		printf("%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %7.2f\n",VX.x(),VX.y(),VX.z(),
			VP.x(),VP.y(),VP.z(),VP.Mag(),VP.Theta()*180./3.14159, VP.Phi()*180./3.14159); 
#endif

#ifdef DRIFT_BENCHMARK 		
		//cout<<"For this OS, CLOCKS_PER_SEC="<<CLOCKS_PER_SEC<<endl;
		clock_t end = clock();
		gTotalTime+=(double)(end-start)*1000.0/(double)CLOCKS_PER_SEC;
		gTotalN++;
		if(DRIFT_BENCHMARK>=2)
		{
			printf("Drift2Ztr_EqualStep() Elapsed time: %6.3f milisecond\r",
				(double)(end-start)*1000.0/(double)CLOCKS_PER_SEC);
		}
		if(!(gTotalN%2000))
		{
			printf("Drift_SL2TG():Drift2Ztr_EqualStep() Average elapsed time: %6.3f milisecond\n",
				gTotalTime/double(gTotalN));
			printf("Drift_SL2TG():Drift2Ztr_EqualStep() Total elapsed time: %6.0f milisecond over %d calls.\n",
				gTotalTime,gTotalN);
		}
#endif
		return 0;
	}

	//////////////////////////////////////////////////////////////////////////////////////

}  //end of namespace
