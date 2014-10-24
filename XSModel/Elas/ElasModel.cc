#include <cstdio>
#include <cmath>
#include <vector>

using namespace std;
#include "ElasModel.hh"
#include "gauss_legendre.h"

namespace  ElasModel
{

	extern "C"
	{
		void carbon_(double* Ei, double* ang, double* xs);
	}

	static const double kDEG = 3.14159265358979323846/180.0;
	static const double kMEV = 1.0e-3;
	static const double kAlpha = 1 / 137.035999679;
	static const double kPi = 3.1415926535897932384626;
	static const double kC = 299792458;
	static const double kQe = 1.602176565e-19;
	static const double kHbar = 1.054571726e-34;
	static const double kFm = kHbar*kC / 1e-15 / 1e9 / kQe;

	static const int QUADRATURE_ORDER = 128;
	static const double kIntegralMin = 0.0;
	static const double kIntegralMax = 10.0;
	static double kRho0GE_He4, kRho0GE_N14, kRho0GM_N14;
	bool bIsInit = false;


	double rhoGE_He4(double x, void* data);
	double rhoGE_N14(double x, void* data);
	double rhoGM_N14(double x, void* data);
	void Init()
	{
		double temp = gauss_legendre(QUADRATURE_ORDER, rhoGE_He4, NULL, kIntegralMin, kIntegralMax);
		kRho0GE_He4 = 1. / temp;
		temp = gauss_legendre(QUADRATURE_ORDER, rhoGE_N14, NULL, kIntegralMin, kIntegralMax);
		kRho0GE_N14 = 1. / temp;
		temp = gauss_legendre(QUADRATURE_ORDER, rhoGM_N14, NULL, kIntegralMin, kIntegralMax);
		kRho0GM_N14 = 1. / temp;
		bIsInit = true;
	}
	double FF_C12(double Q2)
	{
		double XA;
		double Fc12;
		/*
		C     K. Slifer 10/04/02
		C     12C Form Factor
		C     See K. C. Stansfield et al. PRC 3, 1448 (1971)
		*/
		double HBARC = 0.197327053; // GEV-FM
		double ZT = 6.0;            // ATOMIC NUMBER|CHARGE

		double XALPHA = (ZT-2.0)/3.0;
		double Q2FM = Q2/(HBARC*HBARC); // fm^-2 // Corrected on Feb 17, 2011

		if (Q2FM<3.2)
			XA = 1.64;              // fm
		else if(Q2FM>3.5)
			XA = 1.68;              // fm
		else
			XA = 0;                 // fm

		Fc12 = 1.0-XALPHA/(2.0*(2.+3.*XALPHA))*Q2FM*XA*XA;
		Fc12 = Fc12*exp(-(Q2FM*XA*XA)/4.0);

		if ((Q2FM>3.2)&&(Q2FM<3.5)) // DIFFRACTION MINIMUM
			Fc12 = 1.0/100000.0;
		return Fc12;
	}

	double FF_All(double Q2, double Z)
	{
		//input: Q2 in GeV2, Z is atomic number

		double XA = 1.64; //unit in Fm
		double FF;

		double HBARC = 0.197327053; // GEV-FM
		double XALPHA = (Z-2.0)/3.0;
		double Q2FM = Q2/(HBARC*HBARC); // fm^-2 // Corrected on Feb 17, 2011

		FF = 1.0-XALPHA/(2.0*(2.+3.*XALPHA))*Q2FM*XA*XA;
		FF = FF*exp(-(Q2FM*XA*XA)/4.0); 
		if (FF<1.0e-6) FF=1.0e-6;

		return FF;
	}

	//wrapper for Larry's elas C12 XS
	double Carbon(double Ei, double theta)
	{
		double e = Ei*1000;
		double t = theta/kDEG;
		double XS;
		//larry's code take beam energy in MeV and angle in degree
		//return xs in fm^-2
		carbon_(&e, &t, &XS);

		XS *=1.0e4; // fm^2 (1e-30m^2) to microbarn (1e-34m^2)
		return XS;
	}


	// FF parametrization from Arrington Phys. Rev. C 69, 022201 (2004)
	void  FF_H(double qsq, double &gep, double &gmp)
	{
		double fmup = 5.586/2.;
		gep=1 + qsq*2.94+qsq*qsq*3.04 - pow(qsq,3.0)*2.255 + 
			pow(qsq,4.0)*2.002 - 0.5338*pow(qsq,5.0)+
			4.875e-2*pow(qsq,6.0);
		gep=1./gep;

		gmp=1 + qsq*3.0 + 1.39*pow(qsq,2.0) + 0.122*pow(qsq,3.0) - 
			8.34e-3*pow(qsq,4.0) + 4.25e-4*pow(qsq,5) - 
			7.79e-6*pow(qsq,6);
		gmp=1./gmp*fmup;
	}


	//Hydrogen elastic cross section	
	// FF parametrization from Arrington Phys. Rev. C 69, 022201 (2004)
	double XS_H(double Ei, double Theta)
	{
		/* Formula are from Alexandre Deur thesis */

		double Mtg=0.93889;
		double Ef, nu, tau, Q2, Q2fm;
		double XSmott, recoil, Gep, Gmp, W1, W2, XS;

		/* mili barn value in GeV */
		double hbarc2 = 0.38938;

		/* Scattered energy */
		Ef = Ei / (1.0 + 2.0 * Ei * sin(Theta/2.0) * sin(Theta/2.0) / Mtg);

		nu  = Ei - Ef;
		tau = nu / (2.0 * Mtg);   //this is Q^2/(4M^2)
		Q2  = nu * 2.0 * Mtg;
		Q2fm = Q2 / (0.1973269631 * 0.1973269631);

		double alpha = 1.0 / 137.035999679;
		/* Mott cross section */
		XSmott = alpha * cos(Theta/2.0) / (2.0 * Ei * sin(Theta/2.0) * sin(Theta/2.0));
		XSmott = XSmott * XSmott * hbarc2;  // in mili barn

		/* Recoil */
		recoil = Ef / Ei;

		/* Get form factors */
		/* Beware : The current FF_H function needs Q2 in m, not in fm */
		FF_H(Q2, Gep, Gmp);

		W1 = Gmp * Gmp * tau;
		W2 = (Gep * Gep + Gmp * Gmp * tau) / (1.0 + tau);

		/* Cross section in hbarc2*/
		XS = recoil * XSmott * (W2 + 2.0 * tan(Theta/2.0) * tan(Theta/2.0) * W1);

		return XS*1000.0;  //in microbarn
	}

	//Xiaohui's H elas parametrization
	double XS_H1(double Ei, double theta)
	{
		// reference: Xiaohui thesis, eq. 1.39

		int Z = 1;
		double M = 0.93889; // GeV
		double eP = Ei / (1 + Ei / M * (1 - cos(theta))); // scattered electron energy
		double Q2 = 4 * Ei * eP * sin(theta / 2) * sin(theta / 2); // units GeV2
		double tau = Q2 / 4.0 / M / M;

		double x1 = Q2;
		double x2 = x1*Q2;
		double x3 = x2*Q2;
		double x4 = x3*Q2;
		double x5 = x4*Q2;
		double x6 = x5*Q2;  
		double GE = 1.0 / (1.0 + 2.94 * x1 + 3.04 * x2 - 2.255 * x3 + 2.002 * x4 - 0.5338 * x5 + 0.04875 * x6);
		double GM = 2.793 / (1 + 3.0 * x1 + 1.39 * x2 + 0.122 * x3 - 0.00834 * x4 + 4.25e-4 * x5 - 7.79e-6 * x6);

		double sigma = Z * Z * kAlpha * kAlpha / Q2 * (eP / Ei)*(eP / Ei)*(2 * tau * GM * GM + 1 / tan(theta / 2) / tan(theta / 2) / (1 + tau)*(GE * GE + tau * GM * GM));

		return sigma * kFm * kFm * 1e4; // microbarn
	}

	//He3 form factor
	/* Beware : The current FF_He3 function needs Q2 in fm, not in meter */
	void FF_He3(double Q2FM, double & Fc, double & Fm)
	{
		/* Function to compute a SOG form factor
		* Q2FM : momentum transfer squared (fm-2)
		* Nr : number of Gaussians
		* Ga : Gaussians rms (usually 0.8)
		* R2 : Gaussians central positions
		* Q2_c, Q2_m : Gaussians amplitudes
		********************************************************/
		int Nr = 12;
		double Ga, Ga2, a, b, tmp, sum_c, sum_m, Qr, Q;

		Ga = 0.65;

		double R2[]  = { 0.1, 0.5, 0.9, 1.3, 1.6, 2.0,
			2.4, 2.9, 3.4, 4.0, 4.6, 5.2
		};
		double Q2_c[] = { 0.027614, 0.170847, 0.219805, 0.170486,
			0.134453, 0.100953, 0.074310, 0.053970,
			0.023689, 0.017502, 0.002034, 0.004338
		};
		double Q2_m[] = { 0.059785, 0.138368, 0.281326, 0.000037,
			0.289808, 0.019056, 0.114825, 0.042296,
			0.028345, 0.018312, 0.007843, 0.000000
		};

		Q = sqrt(Q2FM);
		Ga2 = Ga * Ga;
		a = exp(-Q2FM * Ga2 / 4.0);

		sum_c = 0.0;
		sum_m = 0.0;

		for(int i=0; i<Nr; i++)
		{
			b = 2.0 * R2[i] * R2[i] / Ga2;
			Qr = Q * R2[i];

			if(Qr>-1.0e-9 && Qr<1.0e-9)
				tmp = 1.0 + b;
			else
				tmp = cos(Qr) + b * sin(Qr) / Qr;

			sum_c = sum_c + tmp * Q2_c[i] / (1.0 + b);
			sum_m = sum_m + tmp * Q2_m[i] / (1.0 + b);
		}

		Fc = a * sum_c;
		Fm = a * sum_m;
	}


	//Helium3 cross section
	double XS_He3(double Ei, double Theta)
	{
		/* Formula are from Alexandre Deur thesis */
		double Mtg=2.9438;
		double Ef, nu, tau, Q2, Q2fm;
		double XSmott, recoil, Fc, Fm, Ge, Gm, W1, W2, XS;

		/* mili barn value in GeV */
		double hbarc2 = 0.38938;

		/* Scattered energy */
		Ef = Ei / (1.0 + 2.0 * Ei * sin(Theta/2.0) * sin(Theta/2.0) / Mtg);

		nu  = Ei - Ef;
		tau = nu / (2.0 * Mtg);
		Q2  = nu * 2.0 * Mtg;
		Q2fm = Q2 / (0.1973269631 * 0.1973269631);

		double alpha = 1.0 / 137.035999679;
		/* Mott cross section */
		XSmott = alpha * cos(Theta/2.0) / (2.0 * Ei * sin(Theta/2.0) * sin(Theta/2.0));
		XSmott = XSmott * XSmott * hbarc2;  //in mb

		/* Recoil */
		recoil = Ef / Ei;

		/* Get form factors */
		/* Beware : The current FF_He3 function needs Q2 in fm, not in meter */
		FF_He3(Q2fm, Fc, Fm);
		Ge = 2.0 * Fc;  // Z * Fc
		Gm = -6.3663 * Fm;  // mu_A * Fm

		W1 = Gm * Gm * tau;
		W2 = (Ge * Ge + Gm * Gm * tau) / (1.0 + tau);

		/* Cross section in hbarc2*/
		XS = recoil * XSmott * (W2 + 2.0 * tan(Theta/2.0) * tan(Theta/2.0) * W1);

		return XS*1000.0;  //in microbarn
	}

	double rhoGE_He4(double x, void* data)
	{
		return 4 * kPi * x * x * (1 + 0.445 * x * x / 1.008 / 1.008) / (1 + exp((x - 1.008) / 0.327));
	}


	double funcGE_He4(double x, void* data)
	{
		double qfm = *(double*) data;
		return kRho0GE_He4 * rhoGE_He4(x, NULL) * sin(qfm * x) / (qfm * x);
	}


	double GE_He4(double q2)
	{
		// calculate GE from density function
		double q = sqrt(q2);
		double qfm = q / kFm;

		double params[1] = {qfm};
		return fabs(gauss_legendre(QUADRATURE_ORDER, funcGE_He4, params, kIntegralMin, kIntegralMax));
	}

	double XS_He4(double Ei, double theta) 
	{  
		if(!bIsInit) Init();
		// reference: xiaohui thesis, eq. 1.39
		const double kAlpha = 1 / 137.035999679;
		int Z = 2;
		double M = 3.7284; // GeV
		double eP = Ei / (1 + Ei / M * (1 - cos(theta))); // scattered electron energy
		double Q2 = 4 * Ei * eP * sin(theta / 2) * sin(theta / 2); // units GeV2
		// printf("theta = %3.1f Q2 = %3.1f\n", theta / kDEG, Q2 / kFm / kFm);
		double tau = Q2 / 4.0 / M / M;
		double GE = GE_He4(Q2);
		double GM = 0.0;
		double sigma = Z*Z * kAlpha*kAlpha / Q2 * (eP/Ei)*(eP/Ei) *
			( 2 * tau * GM*GM + 1 / pow(tan(theta/2),2.0) / (1 + tau) * (GE*GE + tau * GM*GM) );

		return sigma*kFm*kFm*1e4; // microbarn
	}


	//Tantalum Form Factor.
	//input Q2 in GeV^2, theta in rad
	double FF_Ta(double Q2)
	{
		// Note that here F.F(Q2=0)=1, not Z
		// MT  = 180947.88              ! MeV

		double Z,A,Q2FM,XB,XC,XF,FF;
		Z   = 73.;
		A   = 180.94788;

		Q2FM= Q2/0.1973/0.1973;      //! Convert fom GeV^2 to fm^-2
		XB  = 2.4;                   //! fm
		XC  = 1.07*pow(A,1./3.);     //! fm

		XF  = 1.0/(1.0 + 1.0/6.0 * Q2FM * XC*XC);
		XF  = XF*exp(-1./6.*Q2FM*XB*XB);

		FF  = XF*XF;

		return FF;
	}


	double XS_Ta184(double Ei, double theta)
	{
		double Mtg=183.84*0.9315;
		double Recoil, Ef, Q2;

		Recoil = 1.0/(1.0+Ei/Mtg*(1.-cos(theta)));
		Ef = Recoil*Ei;
		Q2 = 2.*Mtg*(Ei-Ef);

		double Z = 74;
		double alpha = 1.0/137.035989561;
		// Calculated cross sections are often written in units of hbar2c2/GeV2 
		//(approximately 0.3894 mb).
		double hbc2 = 0.38938;  // turn GeV to mbar

		double C = cos(theta/2.0);
		double S = Z*alpha*C/(2.*Ei*(1.-C*C));
		double Mott = S*S*hbc2*1000.0; // microbarn

		double FF = FF_All(Q2, Z);

		return Recoil*Mott*FF*FF; // microbarn
	}



	//wrapper for Stansfield's elas C12 XS
	double XS_C12(double Ei, double theta)
	{
		double Recoil, Ef, Q2;
		double Mtg = 11.178;

		/*
		C     Ei        - INCIDENT ENERGY (GEV)
		C     theta     - SCATTERING ANGLE (RADIANS)
		C     Mtg       - TARGET MASS (GEV)
		C     Recoil    - RECOIL FACTOR
		C     Ef        - SCATTERED ENERGY (GEV)
		C     Q2        - MOMENTUM TRANSFER SQUARED (GEV**2)
		C     MOTT      - MOTT CROSS SECTION
		C     FF_C12    - NUCLEAR FORM FACTOR
		*/

		Recoil = 1.0/(1.0+Ei/Mtg*(1.-cos(theta)));
		Ef = Recoil*Ei;
		Q2 = 2.*Mtg*(Ei-Ef);

		double Z = 6;
		double alpha = 1.0/137.035989561;
		// Calculated cross sections are often written in units of hbar2c2/GeV2 (approximately 0.3894 mb).
		double hbc2 = 0.38938;  // turn GeV to mbar

		double C = cos(theta/2.0);
		double S = Z*alpha*C/(2.*Ei*(1.-C*C));
		double Mott = S*S*hbc2*1000.0; // microbarn

		double FF = FF_C12(Q2);

		return Recoil*Mott*FF*FF; // microbarn
	}




	double rhoGM_N14(double x, void* data)
	{
		return 4 * kPi * x * x * (1 + 0.25193 * x * x / 6.751) * exp(-x * x / 6.751);
	}

	double funcGM_N14(double x, void* data)
	{
		double qfm = *(double*) data;
		return kRho0GM_N14 * rhoGM_N14(x, NULL) * sin(qfm * x) / (qfm * x);
	}


	double rhoGE_N14(double x, void* data)
	{
		return 4 * kPi * x * x * (1 + 1.234 * x * x / 3.0976) * exp(-x * x / 3.0976);
	}


	double funcGE_N14(double x, void* data)
	{
		double qfm = *(double*) data;
		return kRho0GE_N14 * rhoGE_N14(x, NULL) * sin(qfm * x) / (qfm * x);
	}

	double GE_N14(double q2)
	{
		// calculate GE from density function

		double q = sqrt(q2);
		double qfm = q / kFm;

		double params[1] = {qfm};
		return fabs(gauss_legendre(QUADRATURE_ORDER, funcGE_N14, params, kIntegralMin, kIntegralMax));
	}


	double GM_N14(double q2)
	{
		// calculate GM from density function
		double q = sqrt(q2);
		double qfm = q / kFm;

		double params[1] = {qfm};
		return fabs(gauss_legendre(QUADRATURE_ORDER, funcGM_N14, params, kIntegralMin, kIntegralMax));
	}

	double XS_N14(double Ei, double theta)
	{
		// reference: Xiaohui thesis, eq. 1.39
		if(!bIsInit) Init();
		int Z = 7;
		double M = 13.04378; // GeV
		double eP = Ei / (1 + Ei / M * (1 - cos(theta))); // scattered electron energy
		double Q2 = 4 * Ei * eP * sin(theta / 2) * sin(theta / 2); // units GeV2
		double tau = Q2 / 4.0 / M / M;
		double GE = GE_N14(Q2);
		double GM = GM_N14(Q2);
		double sigma = Z * Z * kAlpha * kAlpha / Q2 * (eP / Ei)*(eP / Ei)*(2 * tau * GM * GM + 1 / tan(theta / 2) / tan(theta / 2) / (1 + tau)*(GE * GE + tau * GM * GM));

		return sigma * kFm * kFm * 1e4; // microbarn
	}


	double GetXS(int Z, int N, double Ei, double Theta, double Mtg_GeV, int iUseLarry)
	{
		if ((Z==6)&&(N==6)) 
		{
			//Larry Cardman's model is too slow 
			if (iUseLarry==1) return Carbon(Ei, Theta);
			else return XS_C12(Ei, Theta);
		}
		else if ((Z==1)&&(N==0))  return XS_H1(Ei,Theta); 
		else if ((Z==2)&&(N==1))  return XS_He3(Ei,Theta); 
		else if ((Z==2)&&(N==2))  return XS_He4(Ei,Theta); 
		else if ((Z==7)&&(N==7))  return XS_N14(Ei,Theta); 
		else if ((Z==74)&&(N==110))  return XS_Ta184(Ei,Theta); 
		else
		{
			const double AMU = 0.9314941;
			if(Mtg_GeV<=1.0E-9) Mtg_GeV=(Z+N)*AMU;
			return GetXS_All(Ei, Theta, Z, Mtg_GeV);
		}
		return 0.0;
	}

	double GetXS_All(double Ei, double theta, int iZ, double Mtg)
	{
		double Recoil, Ef, Q2;

		Recoil = 1.0/(1.0+Ei/Mtg*(1.-cos(theta)));
		Ef = Recoil*Ei;
		Q2 = 2.*Mtg*(Ei-Ef);

		double Z = iZ;
		double alpha = 1.0/137.035989561;
		// Calculated cross sections are often written in units of hbar2c2/GeV2 (approximately 0.3894 mb).
		double hbc2 = 0.38938;  // turn GeV to mbar

		double C = cos(theta/2.0);
		double S = Z*alpha*C/(2.*Ei*(1.-C*C));
		double Mott = S*S*hbc2*1000.0; // microbarn

		double FF = FF_All(Q2, iZ);

		return Recoil*Mott*FF*FF; // microbarn
	}

}

