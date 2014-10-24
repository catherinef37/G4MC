///////////////////////////////////////////////////////////////////
//By Jixie Zhang, 20131205  jixie@jlab.org
//return Jerry Miller's RCS (photon-proton) XS in nb/GeV^2/sr between 2-12 GeV
//Use 2nd order interpolation to get logXS for each incident energy 
//
//details:
//I fit logXS vs cos(Theta_cm) with pol7 for 16 incident energies (E[0]-E[15])
//to create 16 pol7 curves between 2 to 12 GeV.
//For each given incident energy Ei and Theta_CM, I will find 
//E0,E1,E2 from the energy table such that E0 < Ei <= E1 < E2
//If Ei<E[0], I will use  E0=E[0], E1=E[1], E2=E[2]
//If Ei>E[14], I will use  E0=E[13], E1=E[14], E2=E[15]
//Then calculate logXS0, logXS1, and logXS2 using E0, E1, and E2, 
//Finally interpolate Ei for logXS, return exp(logXS)
//Note that there is no guarantee for photon incident energy beyond 2-12 GeV.   
///////////////////////////////////////////////////////////////////

#ifndef  RCSXS_Miller_H
#define  RCSXS_Miller_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
using namespace std;

namespace Compton
{
  #include "RCS_Miller_Para_Final.inc"  
  //#include "RCS_Miller_Para_test.inc"  
  //#include "RCS_Miller_Para.inc"  

  //input: photon incident energy in GeV and 
  //       photon scattering angle in CM frame in radian 
  //return RCS XS in nb/GeV^2/sr
  double GetRCSXS_old(double Ei, double Theta_cm_g)
  {
    //Pol7 fitted parameters for J. Miller's RCS XS at E=2.000 GeV
    const double dPara_Pol7_E2p0[8] = {
      1.156613E+00, 2.759328E+00, 1.433561E+00, -5.086756E-01, -2.868878E-01, 
      1.215566E+00, 6.199179E-01, -8.282572E-01 
    };

    //Pol7 fitted parameters for J. Miller's RCS XS at E=2.500 GeV
    const double dPara_Pol7_E2p5[8] = {
      1.638564E-02, 3.007689E+00, 1.608685E+00, -4.353289E-01, -3.059095E-01, 
      1.308890E+00, 7.543261E-01, -8.389771E-01 
    };

    //Pol7 fitted parameters for J. Miller's RCS XS at E=3.000 GeV
    const double dPara_Pol7_E3p0[8] = {
      -9.584140E-01, 3.207695E+00, 1.759580E+00, -3.762341E-01, -3.403237E-01, 
      1.386289E+00, 8.932510E-01, -8.282480E-01 
    };

    //Pol7 fitted parameters for J. Miller's RCS XS at E=3.500 GeV
    const double dPara_Pol7_E3p5[8] = {
      -1.810878E+00, 3.367359E+00, 1.891875E+00, -3.039620E-01, -4.017965E-01, 
      1.412227E+00, 1.052525E+00, -7.764302E-01 
    };

    //Pol7 fitted parameters for J. Miller's RCS XS at E=4.000 GeV
    const double dPara_Pol7_E4p0[8] = {
      -2.569402E+00, 3.501014E+00, 2.015439E+00, -2.120783E-01, -4.998494E-01, 
      1.374132E+00, 1.232652E+00, -6.839922E-01 
    };

    //Pol7 fitted parameters for J. Miller's RCS XS at E=4.500 GeV
    const double dPara_Pol7_E4p5[8] = {
      -3.252587E+00, 3.610476E+00, 2.116546E+00, -9.473633E-02, -5.999195E-01, 
      1.275553E+00, 1.417999E+00, -5.569097E-01 
    };

    //Pol7 fitted parameters for J. Miller's RCS XS at E=5.000 GeV
    const double dPara_Pol7_E5p0[8] = {
      -3.875227E+00, 3.714619E+00, 2.200749E+00, -1.963378E-02, -6.922464E-01, 
      1.211639E+00, 1.600043E+00, -4.405773E-01 
    };

    //Pol7 fitted parameters for J. Miller's RCS XS at E=5.500 GeV
    const double dPara_Pol7_E5p5[8] = {
      -4.447882E+00, 3.804472E+00, 2.282810E+00, 8.831010E-02, -7.973278E-01, 
      1.063101E+00, 1.785220E+00, -2.794365E-01 
    };

    //Pol7 fitted parameters for J. Miller's RCS XS at E=6.000 GeV
    const double dPara_Pol7_E6p0[8] = {
      -4.978414E+00, 3.893853E+00, 2.358090E+00, 1.285199E-01, -9.056032E-01, 
      1.013393E+00, 1.973237E+00, -1.668049E-01 
    };

    //Pol7 fitted parameters for J. Miller's RCS XS at E=6.600 GeV
    const double dPara_Pol7_E6p6[8] = {
      -5.568443E+00, 3.997300E+00, 2.450282E+00, 1.214555E-01, -1.017850E+00,
      1.024395E+00, 2.165315E+00, -6.227354E-02
    };

    //Pol7 fitted parameters for J. Miller's RCS XS at E=7.300 GeV
    const double dPara_Pol7_E7p3[8] = {
      -6.202416E+00, 4.092473E+00, 2.583395E+00, 1.907419E-01, -1.239637E+00, 
      8.812899E-01, 2.442635E+00, 1.473698E-01 
    };

    //Pol7 fitted parameters for J. Miller's RCS XS at E=8.000 GeV
    const double dPara_Pol7_E8p0[8] = {
      -6.786660E+00, 4.166332E+00, 2.718320E+00, 2.842456E-01, -1.481214E+00, 
      7.049851E-01, 2.718752E+00, 3.708147E-01 
    };

    //Pol7 fitted parameters for J. Miller's RCS XS at E=8.800 GeV
    const double dPara_Pol7_E8p8[8] = {
      -7.401229E+00, 4.238400E+00, 2.871654E+00, 3.977364E-01, -1.792098E+00,
      5.179595E-01, 3.049175E+00, 6.015992E-01
    };

    //Pol7 fitted parameters for J. Miller's RCS XS at E=10.000 GeV
    const double dPara_Pol7_E10p0[8] = {
      -8.233906E+00, 4.309203E+00, 3.054716E+00, 7.426125E-01, -2.214239E+00, 
      -1.143753E-01, 3.514689E+00, 1.138638E+00 
    };

    //Pol7 fitted parameters for J. Miller's RCS XS at E=11.000 GeV
    const double dPara_Pol7_E11p0[8] = {
      -8.862408E+00, 4.360690E+00, 3.174455E+00, 1.014378E+00, -2.521065E+00,
      -5.998640E-01, 3.869595E+00, 1.541981E+00
    };

    //Pol7 fitted parameters for J. Miller's RCS XS at E=12.000 GeV
    const double dPara_Pol7_E12p0[8] = {
      -9.439849E+00, 4.410439E+00, 3.248512E+00, 1.263725E+00, -2.726996E+00, 
      -1.083364E+00, 4.152567E+00, 1.946779E+00 
    };

    const int N = 16;
    double E[]={2.0,2.5,3.0,3.5,4.0,   
		4.5,5.0,5.5,6.0,6.6,
		7.3,8.0,8.8,10.0,11.0,  
		12.0};
    const double *dPara_Pol7[]={
      dPara_Pol7_E2p0,dPara_Pol7_E2p5,dPara_Pol7_E3p0,dPara_Pol7_E3p5,dPara_Pol7_E4p0,
      dPara_Pol7_E4p5,dPara_Pol7_E5p0,dPara_Pol7_E5p5,dPara_Pol7_E6p0,dPara_Pol7_E6p6,
      dPara_Pol7_E7p3,dPara_Pol7_E8p0,dPara_Pol7_E8p8,dPara_Pol7_E10p0,dPara_Pol7_E11p0,
      dPara_Pol7_E12p0
    };
    const double *par0 = dPara_Pol7[0];
    const double *par1 = dPara_Pol7[1];
    const double *par2 = dPara_Pol7[2];
    double E0=E[0], E1=E[1], E2=E[2];
    if(Ei<E[0])
      {
	//do nothing, using default values
	;
      }
    else if(Ei>E[N-2])   
      {
	par0 = dPara_Pol7[N-3]; 
	par1 = dPara_Pol7[N-2];
	par2 = dPara_Pol7[N-1];
	E0=E[N-3]; E1=E[N-2]; E2=E[N-1];
      }
    else
      {
	for(int i=0; i<N-2; i++)
	  {
	    if(Ei>E[i] && Ei<=E[i+1]) 
	      {
		par0 = dPara_Pol7[i]; 
		par1 = dPara_Pol7[i+1];
		par2 = dPara_Pol7[i+2];
		E0=E[i]; E1=E[i+1]; E2=E[i+2];
		break;
	      }
	  }  
      }
		
		
    //2nd order interpolation
    double cosTheta = cos(Theta_cm_g);
    double XS=0, logXS0=0, logXS1=0, logXS2=0;
    for(int i=0;i<=7;i++) logXS0 += par0[i] * pow(cosTheta,double(i)); 
    for(int i=0;i<=7;i++) logXS1 += par1[i] * pow(cosTheta,double(i)); 
    for(int i=0;i<=7;i++) logXS2 += par2[i] * pow(cosTheta,double(i)); 
		
    double factor0 = (Ei-E1)*(Ei-E2)/((E0-E1)*(E0-E2));
    double factor1 = (Ei-E0)*(Ei-E2)/((E1-E0)*(E1-E2));
    double factor2 = (Ei-E0)*(Ei-E1)/((E2-E0)*(E2-E1));
    double logXS = factor0*logXS0 + factor1*logXS1 + factor2*logXS2;
    XS = exp(logXS);
    //cout<<"E0="<<E0<<"  E1="<<E1<<"  E2="<<E2<<endl;
    //cout<<"par0="<<par0[0]<<"  par1="<<par1[0]<<"  par2="<<par2[0]<<endl;
    //cout<<"logXS0="<<logXS0<<"  logXS1="<<logXS1<<"  logXS2="<<logXS2<<endl;
    return XS;
  }

  double GetRCSXS_Miller(double Ei, double Theta_cm_g)
  {
    const int N = kEnergyN;
    const double *E = kIncidentEnergy;
    const double **dPara_Pol7 = dPara_Pol7_XS;

    const double *par0 = dPara_Pol7[0];
    const double *par1 = dPara_Pol7[1];
    const double *par2 = dPara_Pol7[2];
    double E0=E[0], E1=E[1], E2=E[2];
    if(Ei<E[0])
      {
	//do nothing, using default values
	;
      }
    else if(Ei>E[N-2])   
      {
	par0 = dPara_Pol7[N-3]; 
	par1 = dPara_Pol7[N-2];
	par2 = dPara_Pol7[N-1];
	E0=E[N-3]; E1=E[N-2]; E2=E[N-1];
      }
    else
      {
	for(int i=0; i<N-2; i++)
	  {
	    if(Ei>E[i] && Ei<=E[i+1]) 
	      {
		par0 = dPara_Pol7[i]; 
		par1 = dPara_Pol7[i+1];
		par2 = dPara_Pol7[i+2];
		E0=E[i]; E1=E[i+1]; E2=E[i+2];
		break;
	      }
	  }  
      }


    //2nd order interpolation
    double cosTheta = cos(Theta_cm_g);
    double XS=0, logXS0=0, logXS1=0, logXS2=0;
    for(int i=0;i<=7;i++) logXS0 += par0[i] * pow(cosTheta,double(i)); 
    for(int i=0;i<=7;i++) logXS1 += par1[i] * pow(cosTheta,double(i)); 
    for(int i=0;i<=7;i++) logXS2 += par2[i] * pow(cosTheta,double(i)); 

    double factor0 = (Ei-E1)*(Ei-E2)/((E0-E1)*(E0-E2));
    double factor1 = (Ei-E0)*(Ei-E2)/((E1-E0)*(E1-E2));
    double factor2 = (Ei-E0)*(Ei-E1)/((E2-E0)*(E2-E1));
    double logXS = factor0*logXS0 + factor1*logXS1 + factor2*logXS2;
    XS = exp(logXS);

    //cout<<"E0="<<E0<<"  E1="<<E1<<"  E2="<<E2<<endl;
    //cout<<"par0="<<par0[0]<<"  par1="<<par1[0]<<"  par2="<<par2[0]<<endl;
    //cout<<"logXS0="<<logXS0<<"  logXS1="<<logXS1<<"  logXS2="<<logXS2<<endl;
    return XS;

  }

  double GetRCSXS(double Ei, double Theta_CM)
  {
    const double E0=3.15, E1=4.275, E2=5.35;
    const double kPara0[]={1.11, 0.712918, 0.57};
    
    double pPara[6]={0.712918,-0.2226,0.4792,-11.11,24.87,23.78};

    if(Ei<E0-0.2)  pPara[0] = kPara0[0];
    else if (Ei>E2+0.2)  pPara[0] = kPara0[2];
    else
      {
	//interpolate for pPara[0]
	double factor0 = (Ei-E1)*(Ei-E2)/((E0-E1)*(E0-E2));
	double factor1 = (Ei-E0)*(Ei-E2)/((E1-E0)*(E1-E2));
	double factor2 = (Ei-E0)*(Ei-E1)/((E2-E0)*(E2-E1));
	pPara[0] = factor0*kPara0[0] + factor1*kPara0[1] + factor2*kPara0[2];
      }

    double x = cos(Theta_CM);
    double Data2Miller = pPara[0];
    for(int i=1;i<=5;i++) Data2Miller += pPara[i]*pow(x,double(i));
    double XSMiller = GetRCSXS_Miller(Ei,Theta_CM);  //in nb
    return XSMiller*Data2Miller;
  }

  double GetRCSKLL(double Ei, double Theta_cm_g)
  {
    const int N = kEnergyN;
    const double *E = kIncidentEnergy;
    const double **dPara_Pol7 = dPara_Pol7_KLL;

    const double *par0 = dPara_Pol7[0];
    const double *par1 = dPara_Pol7[1];
    const double *par2 = dPara_Pol7[2];
    double E0=E[0], E1=E[1], E2=E[2];
    if(Ei<E[0])
      {
	//do nothing, using default values
	;
      }
    else if(Ei>E[N-2])   
      {
	par0 = dPara_Pol7[N-3]; 
	par1 = dPara_Pol7[N-2];
	par2 = dPara_Pol7[N-1];
	E0=E[N-3]; E1=E[N-2]; E2=E[N-1];
      }
    else
      {
	for(int i=0; i<N-2; i++)
	  {
	    if(Ei>E[i] && Ei<=E[i+1]) 
	      {
		par0 = dPara_Pol7[i]; 
		par1 = dPara_Pol7[i+1];
		par2 = dPara_Pol7[i+2];
		E0=E[i]; E1=E[i+1]; E2=E[i+2];
		break;
	      }
	  }  
      }

    //2nd order interpolation
    double cosTheta = cos(Theta_cm_g);
    double Var=0, Var0=0, Var1=0, Var2=0;
    for(int i=0;i<=7;i++) Var0 += par0[i] * pow(cosTheta,double(i)); 
    for(int i=0;i<=7;i++) Var1 += par1[i] * pow(cosTheta,double(i)); 
    for(int i=0;i<=7;i++) Var2 += par2[i] * pow(cosTheta,double(i)); 

    double factor0 = (Ei-E1)*(Ei-E2)/((E0-E1)*(E0-E2));
    double factor1 = (Ei-E0)*(Ei-E2)/((E1-E0)*(E1-E2));
    double factor2 = (Ei-E0)*(Ei-E1)/((E2-E0)*(E2-E1));
    Var = factor0*Var0 + factor1*Var1 + factor2*Var2;
			
    return Var;

  }

  double GetRCSALL(double Ei, double Theta_cm_g)
  {
    const int N = kEnergyN;
    const double *E = kIncidentEnergy;
    const double **dPara_Pol7 = dPara_Pol7_ALL;

    const double *par0 = dPara_Pol7[0];
    const double *par1 = dPara_Pol7[1];
    const double *par2 = dPara_Pol7[2];
    double E0=E[0], E1=E[1], E2=E[2];
    if(Ei<E[0])
      {
	//do nothing, using default values
	;
      }
    else if(Ei>E[N-2])   
      {
	par0 = dPara_Pol7[N-3]; 
	par1 = dPara_Pol7[N-2];
	par2 = dPara_Pol7[N-1];
	E0=E[N-3]; E1=E[N-2]; E2=E[N-1];
      }
    else
      {
	for(int i=0; i<N-2; i++)
	  {
	    if(Ei>E[i] && Ei<=E[i+1]) 
	      {
		par0 = dPara_Pol7[i]; 
		par1 = dPara_Pol7[i+1];
		par2 = dPara_Pol7[i+2];
		E0=E[i]; E1=E[i+1]; E2=E[i+2];
		break;
	      }
	  }  
      }

    //2nd order interpolation
    double cosTheta = cos(Theta_cm_g);
    double Var=0, Var0=0, Var1=0, Var2=0;
    for(int i=0;i<=7;i++) Var0 += par0[i] * pow(cosTheta,double(i)); 
    for(int i=0;i<=7;i++) Var1 += par1[i] * pow(cosTheta,double(i)); 
    for(int i=0;i<=7;i++) Var2 += par2[i] * pow(cosTheta,double(i)); 

    double factor0 = (Ei-E1)*(Ei-E2)/((E0-E1)*(E0-E2));
    double factor1 = (Ei-E0)*(Ei-E2)/((E1-E0)*(E1-E2));
    double factor2 = (Ei-E0)*(Ei-E1)/((E2-E0)*(E2-E1));
    Var = factor0*Var0 + factor1*Var1 + factor2*Var2;
			
    return Var;
  }

  //par: Theta_cm_deg
  double XSVsE_FixedThetacm_deg(double *x, double *par)
  {
    const double deg = atan(1.0)/45.;
    double Ei = x[0];
    double Theta_cm = par[0]*deg;
    //cout<<"XSVsE_FixedThetacm_deg(): Theta_cm="<<Theta_cm/deg<<endl;
    return GetRCSXS(Ei, Theta_cm);
  }

  double GetTheta_cm(double Ei, double Theta_lab_rad)
  {
    const double M = 0.9383;
    double cosTheta = cos(Theta_lab_rad);
    double Ef = Ei/(1+Ei/M*(1.0-cosTheta));
    double t  = -2.0*Ei*Ef*(1.0-cosTheta);
    double s  = M*M + 2.0*M*Ei;
    double sinHalfTheta_cm = sqrt(-t*s)/(s-M*M);
    double Theta_cm = 2.0*asin(sinHalfTheta_cm);
    return Theta_cm;
  }

  //par: Theta_lab_deg
  double XSVsE_FixedThetalab_deg(double *x, double *par)
  {
    const double deg = atan(1.0)/45.;
    double Ei = x[0];
    double Theta_lab = par[0]*deg;
    double Theta_cm = GetTheta_cm(Ei,Theta_lab);
    return GetRCSXS(Ei, Theta_cm);
  }

  //return RCS XS in nbar
  double GetRCSXS_lab(double Ei, double Theta_lab_rad)
  {
    double Theta_cm = GetTheta_cm(Ei,Theta_lab_rad);
    return GetRCSXS(Ei, Theta_cm);
  }

  //return RCS KLL
  double GetRCSKLL_lab(double Ei, double Theta_lab_rad)
  {
    double Theta_cm = GetTheta_cm(Ei,Theta_lab_rad);
    return GetRCSKLL(Ei, Theta_cm);
  }
  //return RCS ALL
  double GetRCSALL_lab(double Ei, double Theta_lab_rad)
  {
    double Theta_cm = GetTheta_cm(Ei,Theta_lab_rad);
    return GetRCSALL(Ei, Theta_cm);
  }

  ////////////////////////////////////////////////////////
  //x is cos(Theta_cm)
  double RCSXSVsCosTh_cm(double *x, double *par)
  {
    double E=par[0];
    double Theta_cm_g=acos(x[0]);
    return GetRCSXS(E,Theta_cm_g);
  }

  //x is cos(Theta_cm)
  double RCSKLLVsCosTh_cm(double *x, double *par)
  {
    double E=par[0];
    double Theta_cm_g=acos(x[0]);
    return GetRCSKLL(E,Theta_cm_g);
  }

  //x is cos(Theta_cm)
  double RCSALLVsCosTh_cm(double *x, double *par)
  {
    double E=par[0];
    double Theta_cm_g=acos(x[0]);
    return GetRCSALL(E,Theta_cm_g);
  }

  //x is Theta_cm
  double RCSXSVsTh_cm(double *x, double *par)
  {
    double E=par[0];
    double Theta_cm_g=x[0];
    return GetRCSXS(E,Theta_cm_g);
  }

  //x is Theta_cm
  double RCSKLLVsTh_cm(double *x, double *par)
  {
    double E=par[0];
    double Theta_cm_g=x[0];
    return GetRCSKLL(E,Theta_cm_g);
  }

  //x is Theta_cm
  double RCSALLVsTh_cm(double *x, double *par)
  {
    double E=par[0];
    double Theta_cm_g=x[0];
    return GetRCSALL(E,Theta_cm_g);
  }

}  //end of namespace
#endif

