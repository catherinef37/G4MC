////////////////////////////////////////////////////////////////////////
// Standard HRS transport functions FOR PREX
// 
////////////////////////////////////////////////////////////////////////

#include <cmath>

//#include "prex_forward_tuneX.h"
#include "prex_forward.hh"
#include "hamcPREXTrans.hh"
#include <iostream>

//#define NICKIETEST 1
//#define ACCTEST 1
//using namespace SNoSepta;
using namespace std;

const float m2cm = 100.0;
const double kDEG = 3.14159265358979323846/180.0;

hamcPREXTrans::hamcPREXTrans()//:cModelAngle(12.5*kDEG)
:cModelAngle(5.0*kDEG)
{
    // Nothing to do
}

hamcPREXTrans::~hamcPREXTrans()
{
    // Nothing to do
}

bool hamcPREXTrans::TransLeftHRS(double* pV5)
{
  float vector_jjl[]={(float)pV5[0],(float)pV5[1],(float)pV5[2],(float)pV5[3],(float)pV5[4]};
	int iii = 5; int *ii = &iii;

#if defined NICKIETEST
	for( int i = 0; i < 5; i++){
	  cout << pV5[i] << " ";
	}
	cout << endl;
#endif
	float x_test, y_test;
	
	//Nickie adds this to take into account acceptance from septum
	//x_test = x_sp_q1en_(vector_jjl, ii)*m2cm;
        //y_test = y_sp_q1en_(vector_jjl, ii)*m2cm;
        //x_test = x_test - 0.9;
        //if( (x_test*x_test + y_test*y_test) > (14.92*14.92) )
	//return false;
	
	//Target to Q1 exit, circle of radius 14.92 cm
	x_test = x_sp_q1ex_(vector_jjl, ii)*m2cm;
	y_test = y_sp_q1ex_(vector_jjl, ii)*m2cm;
	//x_test = x_test - 0.9;
	if( (x_test*x_test + y_test*y_test) > (14.92*14.92) )
		return false;

	//Nickie: this used to be commented
	//Target to dipole entrance, trapezoid 522.0cm>x>498.1cm  |y| < -0.1924*x-19.24
	//x_test = x_sp_den_(vector_jjl, ii)*m2cm;
	//y_test = y_sp_den_(vector_jjl, ii)*m2cm;
	//if( (x_test>522.0) || (x_test<498.1) || fabs(y_test) > fabs(-0.1924*x_test-19.24) )
	//return false;

	//Target to dipole exit, trapezoid -46.19cm<x<46.19cm  |y| < -0.0161*x+12.5
	x_test = x_sp_dex_(vector_jjl, ii)*m2cm;
	y_test = y_sp_dex_(vector_jjl, ii)*m2cm;
	//cout<<"dipole_exit:(x,y)=\t"<<x_test<<"\t "<<y_test<<endl;
	if( fabs(x_test)>46.19 || fabs(y_test) > fabs(-0.0161*x_test+12.5) )
		return false;

	//Target to Q3 entrance, circle of radius 30.0 cm
	x_test = x_sp_q3en_(vector_jjl, ii)*m2cm;
	y_test = y_sp_q3en_(vector_jjl, ii)*m2cm;
	if( (x_test*x_test + y_test*y_test) > (30.0*30.0) )
		return false;

	//Target to Q3 exit, circle of radius 30.0 cm  -> 28.0cm
	x_test = x_sp_q3ex_(vector_jjl, ii)*m2cm;
	y_test = y_sp_q3ex_(vector_jjl, ii)*m2cm;
	//x_test = (x_test + 1.0) / (28.0);
	//y_test = y_test / (30.0);
	if( (x_test*x_test + y_test*y_test) > (30.0*30.0) )
		return false;

	/////////////////////////////////////////////////////////////
	/* If we reach this point, it means the test was succesful */
	float x_fp     = x_sp_fp_(vector_jjl, ii);
	float theta_fp = t_sp_fp_(vector_jjl, ii);
	float y_fp     = y_sp_fp_(vector_jjl, ii);
	float phi_fp   = p_sp_fp_(vector_jjl, ii);

	//reset the vector and return it back to the caller
	pV5[0] = (double)x_fp;
	pV5[1] = (double)theta_fp;
	pV5[2] = (double)y_fp;
	pV5[3] = (double)phi_fp;
	//pV5[4] = (double)delta_fp;  // delta is not change

#if defined NICKIETEST
	for( int i = 0; i < 5; i++){
	  cout << pV5[i] << " ";
	}
	cout << endl;
#endif
	return true;
}

//bool hamcPREXTrans::TransRightHRS(double* pV5){
//return false;
//}

bool hamcPREXTrans::TransRightHRS(double* pV5)
{
  //cout << "Transporting.\n";
	float vector_jjl[]={(float)pV5[0],(float)pV5[1],(float)pV5[2],(float)pV5[3],(float)pV5[4]};
	int iii = 5; int *ii = &iii;
#if defined NICKIETEST
	for( int i = 0; i < 5; i++){
	  cout << pV5[i] << " ";
	}
	cout << endl;
#endif
	float x_test, y_test;

	//Nickie adds the collimator:
	x_test = x_sp_col_(vector_jjl, ii);
        y_test = y_sp_col_(vector_jjl, ii);
	//cout << "!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	//cout << x_test << " " << y_test << endl;
	float atxlo = 0.032;
	float atyhi = 0.041;
	float atr   = 0.20;
	float atc   = 0.229;
	float rad1  = 0.205;
	float ycent1= 0.145;
	float rad2  = 0.20;
	float ycent2= 0.229;
	float vtop  = 0.117;
	float vright= 0.04;
	float yline1= 0.1474;
	float slope1= -1.88;
	float ysign = 1.0;

	float xdif = x_test;
	float ydif = ysign*(y_test);
	float xtrial,ytrial,xtmp;
	
	ytrial = ydif-atc;
	xtrial = atr*atr - ytrial*ytrial;
	if (xtrial > 0) {
	  xtmp = sqrt(xtrial);
	  if ((xdif > atxlo) && (ydif < atyhi) &&
	      (x_test < xtmp)) return false;     // in A_T hole (a triangle with arc)
	}
	// Now check if in main acceptance
	// Remember x is vertical (transport), y is horizontal
	ytrial = ydif-ycent1;
	xtrial = rad1*rad1 - ytrial*ytrial;

	if (xtrial < 0) return false;   // outside outer circle
	xtrial = sqrt(xtrial);

	if (x_test > xtrial || x_test < -1.0*xtrial) return false;  // outside outer circle
	ytrial = ydif-ycent2;
	xtrial = rad2*rad2 - ytrial*ytrial;
	if (xtrial > 0) {
	  xtrial = sqrt(xtrial);

	  if ((x_test > 0 && x_test < xtrial) ||
	      (x_test < 0 && x_test > -1.0*xtrial)) return false; // outside innner circle
	}
	if (x_test > vtop) return false;       // above upper line
	if (x_test < -1.0*vtop) return false;  // beneath lower line
	if (ydif > vright) return false;     // beyond right border
	xtrial = slope1*ydif + yline1; // Champhor line
	if (x_test > xtrial) return false;     // above upper line
	if (x_test < -1.0*xtrial) return false;// below lower line

	//return 1;  // inside main acceptance.


	//Target to Q1 exit, circle of radius 14.92 cm
	x_test = x_sp_q1ex_(vector_jjl, ii)*m2cm;
	y_test = y_sp_q1ex_(vector_jjl, ii)*m2cm;
	//x_test = x_test + 0.9;
	//cout << "Man, it didn't even pass the first cut!!" << endl;
	//cout << x_test*x_test << " " <<  y_test*y_test << " " << 14.92*14.92 << endl;
	if( (x_test*x_test + y_test*y_test) > (14.92*14.92) )
		return false;
	//cout << "q1ex" << endl;
	//Target to dipole entrance, trapezoid -522.0cm<x<-498.1cm  |y| < -0.1924*x-19.24
	//x_test = x_sp_dent_(vector_jjl, ii)*m2cm;
	//y_test = y_sp_dent_(vector_jjl, ii)*m2cm;
	//if( (x_test<-522.0) || (x_test>-498.1) || fabs(y_test) > fabs(-0.1924*x_test-19.24) )
	//	return false;

	//Target to dipole exit, trapezoid -46.19cm<x<46.19cm  |y| < -0.0161*x+12.5
	x_test = x_sp_dex_(vector_jjl, ii)*m2cm;
	y_test = y_sp_dex_(vector_jjl, ii)*m2cm;
	if( (x_test<-46.19) || (x_test>46.19) || fabs(y_test) > fabs(-0.0161*x_test+12.5) )
		return false;
	//cout << "dex" << endl;
	//Target to Q3 entrance, circle of radius 30.0 cm
	x_test = x_sp_q3en_(vector_jjl, ii)*m2cm;
	y_test = y_sp_q3en_(vector_jjl, ii)*m2cm;
	if( (x_test*x_test + y_test*y_test) > (30.0*30.0) )
		return false;
	//cout << "q3en" << endl;
	//Target to Q3 exit, circle of radius 30.0 cm  -> 28.0cm
	x_test = x_sp_q3ex_(vector_jjl, ii)*m2cm;
	y_test = y_sp_q3ex_(vector_jjl, ii)*m2cm;
	//x_test = (x_test - 1.0) / (28.0);
	//y_test = y_test / (30.0);
	if( (x_test*x_test + y_test*y_test) > (30.0*30.0) )
		return false;
	//cout << "q3ex" << endl;
	/////////////////////////////////////////////////////////////
	/* If we reach this point, it means the test was succesful */

	float x_fp     = x_sp_fp_(vector_jjl, ii);
	float theta_fp = t_sp_fp_(vector_jjl, ii);
	float y_fp     = y_sp_fp_(vector_jjl, ii);
	float phi_fp   = p_sp_fp_(vector_jjl, ii);
#if defined NICKIETEST
	cout << x_fp << " " << theta_fp << " " << y_fp << " " << phi_fp << endl;
#endif

	//reset the vector and return it back to the caller
	pV5[0] = (double)x_fp;
	pV5[1] = (double)theta_fp;
	pV5[2] = (double)y_fp;
	pV5[3] = (double)phi_fp;
	//pV5[4] = (double)delta_fp;  // delta is not change
#if defined NICKIETEST
	for( int i = 0; i < 5; i++){
	  cout << pV5[i] << " ";
	}
	cout << endl;
#endif
	return true;
}


bool hamcPREXTrans::TransRightHRS_C(double* pV5)
{
  //cout << "Transporting, but C mode\n";
	float vector_jjl[]={(float)pV5[0],(float)pV5[1],(float)pV5[2],(float)pV5[3],(float)pV5[4]};
	int iii = 5; int *ii = &iii;
#if defined NICKIETEST
	for( int i = 0; i < 5; i++){
	  cout << pV5[i] << " ";
	}
	cout << endl;
#endif
	float x_test, y_test;
	
	//Nickie adds the collimator:
	x_test = vector_jjl[0];
        y_test = vector_jjl[2];
	//cout << "!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	float atxlo = 0.032;
	float atyhi = 0.041;
	float atr   = 0.20;
	float atc   = 0.229;
	float rad1  = 0.205;
	float ycent1= 0.145;
	float rad2  = 0.20;
	float ycent2= 0.229;
	float vtop  = 0.117;
	float vright= 0.04;
	float yline1= 0.1474;
	float slope1= -1.88;
	float ysign = 1.0;

	float xdif = x_test;
	float ydif = ysign*(y_test);
	float xtrial,ytrial,xtmp;
	
	ytrial = ydif-atc;
	xtrial = atr*atr - ytrial*ytrial;
	if (xtrial > 0) {
	  xtmp = sqrt(xtrial);
	  if ((xdif > atxlo) && (ydif < atyhi) &&
	      (x_test < xtmp)) return false;     // in A_T hole (a triangle with arc)
	}
	// Now check if in main acceptance
	// Remember x is vertical (transport), y is horizontal
	ytrial = ydif-ycent1;
	xtrial = rad1*rad1 - ytrial*ytrial;

	if (xtrial < 0) return false;   // outside outer circle
	xtrial = sqrt(xtrial);

	if (x_test > xtrial || x_test < -1.0*xtrial) return false;  // outside outer circle
	ytrial = ydif-ycent2;
	xtrial = rad2*rad2 - ytrial*ytrial;
	if (xtrial > 0) {
	  xtrial = sqrt(xtrial);

	  if ((x_test > 0 && x_test < xtrial) ||
	      (x_test < 0 && x_test > -1.0*xtrial)) return false; // outside innner circle
	}
	if (x_test > vtop) return false;       // above upper line
	if (x_test < -1.0*vtop) return false;  // beneath lower line
	if (ydif > vright) return false;     // beyond right border
	xtrial = slope1*ydif + yline1; // Champhor line
	if (x_test > xtrial) return false;     // above upper line
	if (x_test < -1.0*xtrial) return false;// below lower line

	//return 1;  // inside main acceptance.
	//cout << "Passed in transfer" << endl;
	//Target to Q1 exit, circle of radius 14.92 cm
	//cout << "ALRIGHT, Let's see how it goes!! " << endl;
	x_test = x_sp_cq1x_(vector_jjl, ii)*m2cm;
	y_test = y_sp_cq1x_(vector_jjl, ii)*m2cm;
	//x_test = x_test + 0.9;
	//cout << "Man, it didn't even pass the first cut!!" << endl;
	//cout << x_test*x_test << " " <<  y_test*y_test << " " << 14.92*14.92 << endl;
	if( (x_test*x_test + y_test*y_test) > (14.92*14.92) )
		return false;
	//cout << "q1ex" << endl;
	//Target to dipole entrance, trapezoid -522.0cm<x<-498.1cm  |y| < -0.1924*x-19.24
	//x_test = x_sp_dent_(vector_jjl, ii)*m2cm;
	//y_test = y_sp_dent_(vector_jjl, ii)*m2cm;
	//if( (x_test<-522.0) || (x_test>-498.1) || fabs(y_test) > fabs(-0.1924*x_test-19.24) )
	//	return false;

	//Target to dipole exit, trapezoid -46.19cm<x<46.19cm  |y| < -0.0161*x+12.5
	x_test = x_sp_cdex_(vector_jjl, ii)*m2cm;
	y_test = y_sp_cdex_(vector_jjl, ii)*m2cm;
	if( (x_test<-46.19) || (x_test>46.19) || fabs(y_test) > fabs(-0.0161*x_test+12.5) )
		return false;
	//cout << "dex" << endl;
	//Target to Q3 entrance, circle of radius 30.0 cm
	x_test = x_sp_cq3e_(vector_jjl, ii)*m2cm;
	y_test = y_sp_cq3e_(vector_jjl, ii)*m2cm;
	if( (x_test*x_test + y_test*y_test) > (30.0*30.0) )
		return false;
	//cout << "q3en" << endl;
	//Target to Q3 exit, circle of radius 30.0 cm  -> 28.0cm
	x_test = x_sp_cq3x_(vector_jjl, ii)*m2cm;
	y_test = y_sp_cq3x_(vector_jjl, ii)*m2cm;
	//x_test = (x_test - 1.0) / (28.0);
	//y_test = y_test / (30.0);
	if( (x_test*x_test + y_test*y_test) > (30.0*30.0) )
		return false;
	//cout << "q3ex" << endl;
	/////////////////////////////////////////////////////////////
	/* If we reach this point, it means the test was succesful */

	float x_fp     = x_sp_cfp_(vector_jjl, ii);
	float theta_fp = t_sp_cfp_(vector_jjl, ii);
	float y_fp     = y_sp_cfp_(vector_jjl, ii);
	float phi_fp   = p_sp_cfp_(vector_jjl, ii);
#if defined NICKIETEST
	cout << x_fp << " " << theta_fp << " " << y_fp << " " << phi_fp << endl;
#endif

	//reset the vector and return it back to the caller
	pV5[0] = (double)x_fp;
	pV5[1] = (double)theta_fp;
	pV5[2] = (double)y_fp;
	pV5[3] = (double)phi_fp;
	//pV5[4] = (double)delta_fp;  // delta is not change
#if defined NICKIETEST
	for( int i = 0; i < 5; i++){
	  cout << pV5[i] << " ";
	}
	cout << endl;
#endif
	return true;
}


void hamcPREXTrans::ReconLeftHRS(double* pV5){
  return;
}
/*
{   
	float vector_jjl[]={pV5[0],pV5[1],pV5[2],pV5[3],pV5[4]};
	int iii = 5; int *ii = &iii;
    int jjj = 1; int *jj = &jjj;
    
	vector_jjl[1] = vector_jjl[1] - txfit_sp_(vector_jjl, jj);
	float x_or      = vector_jjl[4];
	float delta_rec = delta_sp_(vector_jjl, ii);
	float theta_rec = theta_sp_(vector_jjl, ii);
	float phi_rec   = phi_sp_(vector_jjl, ii); 
	float y_rec     = y00_sp_(vector_jjl, ii); 

	//reset the vector and return it back to the caller
	pV5[0] = (double)x_or;
	pV5[1] = (double)theta_rec;
	pV5[2] = (double)y_rec;
	pV5[3] = (double)phi_rec;
	pV5[4] = (double)delta_rec;
}
*/

void hamcPREXTrans::ReconRightHRS(double* pV5){
  return;
}
/*
{
	float vector_jjl[]={pV5[0],pV5[1],pV5[2],pV5[3],pV5[4]};
	int iii = 5; int *ii = &iii;
    int jjj = 1; int *jj = &jjj;
*/  
	/* Orthogonalize theta as JJL asks*/
/*	vector_jjl[1]   = vector_jjl[1] - txfit_r12p5_(vector_jjl, jj);
	float x_or      = vector_jjl[4];
	float delta_rec = delta_r12p5_(vector_jjl, ii);
	float theta_rec = theta_r12p5_(vector_jjl, ii);
	float phi_rec   = phi_r12p5_(vector_jjl, ii); 
	float y_rec     = y00_r12p5_(vector_jjl, ii); 

	//reset the vector and return it back to the caller
	pV5[0] = (double)x_or;
	pV5[1] = (double)theta_rec;
	pV5[2] = (double)y_rec;
	pV5[3] = (double)phi_rec;
	pV5[4] = (double)delta_rec;

} 
*/

void hamcPREXTrans::Acceptance(double* pV5, double* x_check, double* theta_check, double* y_check, double* phi_check, int* acc_bool){
  //cout << "Checking the acceptance..." << endl;
  float vector_jjl[]={(float)pV5[0],(float)pV5[1],(float)pV5[2],(float)pV5[3],(float)pV5[4]};
  int iii = 5; int *ii = &iii;

  x_check[0]      = x_sp_q1ex_(vector_jjl, ii);  y_check[0]    = y_sp_q1ex_(vector_jjl, ii);
  theta_check[0]  = t_sp_q1ex_(vector_jjl, ii);  phi_check[0]  = p_sp_q1ex_(vector_jjl, ii);
  x_check[1]      = x_sp_q2ex_(vector_jjl, ii);  y_check[1]    = y_sp_q2ex_(vector_jjl, ii);
  theta_check[1]  = t_sp_q2ex_(vector_jjl, ii);  phi_check[1]  = p_sp_q2ex_(vector_jjl, ii);
  x_check[2]      = x_sp_den_(vector_jjl, ii);   y_check[2]    = y_sp_den_(vector_jjl, ii);
  theta_check[2]  = t_sp_den_(vector_jjl, ii);   phi_check[2]  = p_sp_den_(vector_jjl, ii);
  x_check[3]      = x_sp_dex_(vector_jjl, ii);   y_check[3]    = y_sp_dex_(vector_jjl, ii);
  theta_check[3]  = t_sp_dex_(vector_jjl, ii);   phi_check[3]  = p_sp_dex_(vector_jjl, ii);
  x_check[4]      = x_sp_q3en_(vector_jjl, ii);  y_check[4]    = y_sp_q3en_(vector_jjl, ii);
  theta_check[4]  = t_sp_q3en_(vector_jjl, ii);  phi_check[4]  = p_sp_q3en_(vector_jjl, ii);
  x_check[5]      = x_sp_q3ex_(vector_jjl, ii);  y_check[5]    = y_sp_q3ex_(vector_jjl, ii);
  theta_check[5]  = t_sp_q3ex_(vector_jjl, ii);  phi_check[5]  = p_sp_q3ex_(vector_jjl, ii);
  x_check[6]      = x_sp_sen_(vector_jjl, ii);   y_check[6]    = y_sp_sen_(vector_jjl, ii);
  theta_check[6]  = t_sp_sen_(vector_jjl, ii);   phi_check[6]  = p_sp_sen_(vector_jjl, ii);
  x_check[7]      = x_sp_sm_(vector_jjl, ii);    y_check[7]    = y_sp_sm_(vector_jjl, ii);
  theta_check[7]  = t_sp_sm_(vector_jjl, ii);    phi_check[7]  = p_sp_sm_(vector_jjl, ii);
  x_check[8]      = x_sp_sex_(vector_jjl, ii);   y_check[8]    = y_sp_sex_(vector_jjl, ii);
  theta_check[8]  = t_sp_sex_(vector_jjl, ii);   phi_check[8]  = p_sp_sex_(vector_jjl, ii);
  x_check[9]      = x_sp_col_(vector_jjl, ii);   y_check[9]    = y_sp_col_(vector_jjl, ii);
  theta_check[9]  = t_sp_col_(vector_jjl, ii);   phi_check[9]  = p_sp_col_(vector_jjl, ii);
  x_check[10]     = 0.;                          y_check[10]   = 0.;
  theta_check[10] = 0.;                          phi_check[10] = 0.;

  float x_test, y_test;

  //Nickie adds the septum cuts
  //sen: +0.088 < x < +0.382
  //-0.120 < y < 0.120
  //sm: +0.088 < x < +0.382
  //-0.120 < y < 0.120
  //sex: +0.088 < x < +0.382
  //-0.120 < y < 0.120  

  x_test = x_sp_sen_(vector_jjl, ii);
  y_test = y_sp_sen_(vector_jjl, ii);
  if( x_test < 0.088 || x_test > 0.382 || y_test < -0.120 || y_test > 0.120 )
    return;

#ifdef ACCTEST
  cout << "sen ";
#endif 

  acc_bool[6] = 1;

  x_test = x_sp_sm_(vector_jjl, ii);
  y_test = y_sp_sm_(vector_jjl, ii);
  if( x_test < 0.088 || x_test > 0.382 || y_test < -0.120 || y_test > 0.120 )
    return;

#ifdef ACCTEST
  cout << "sm ";
#endif 

  acc_bool[7] = 1;

  x_test = x_sp_sex_(vector_jjl, ii);
  y_test = y_sp_sex_(vector_jjl, ii);
  if( x_test < 0.088 || x_test > 0.382 || y_test < -0.120 || y_test > 0.120 )
    return;

#ifdef ACCTEST
  cout << "sex ";
#endif 

  acc_bool[8] = 1;
  
  //Nickie adds the collimator:
  x_test = x_sp_col_(vector_jjl, ii);
  y_test = y_sp_col_(vector_jjl, ii);
  float atxlo = 0.032;
  float atyhi = 0.041;
  float atr   = 0.20;
  float atc   = 0.229;
  float rad1  = 0.205;
  float ycent1= 0.145;
  float rad2  = 0.20;
  float ycent2= 0.229;
  float vtop  = 0.117;
  float vright= 0.04;
  float yline1= 0.1474;
  float slope1= -1.88;
  float ysign = 1.0;
  
  float xdif = x_test;
  float ydif = ysign*(y_test);
  float xtrial,ytrial,xtmp;
  
  ytrial = ydif-atc;
  xtrial = atr*atr - ytrial*ytrial;
  if (xtrial > 0) {
    xtmp = sqrt(xtrial);
    if ((xdif > atxlo) && (ydif < atyhi) &&
	(x_test < xtmp)) return;     // in A_T hole (a triangle with arc)
  }
  
 // Now check if in main acceptance
  // Remember x is vertical (transport), y is horizontal
  ytrial = ydif-ycent1;
  xtrial = rad1*rad1 - ytrial*ytrial;
  if (xtrial < 0) return;   // outside outer circle
  xtrial = sqrt(xtrial);
  if (x_test > xtrial || x_test < -1.0*xtrial) return;  // outside outer circle
  ytrial = ydif-ycent2;
  xtrial = rad2*rad2 - ytrial*ytrial;
  if (xtrial > 0) {
    xtrial = sqrt(xtrial);
    
    if ((x_test > 0 && x_test < xtrial) ||
	(x_test < 0 && x_test > -1.0*xtrial)) return; // outside innner circle
  }
  if (x_test > vtop) return;       // above upper line
  if (x_test < -1.0*vtop) return;  // beneath lower line
  if (ydif > vright) return;     // beyond right border
  xtrial = slope1*ydif + yline1; // Champhor line
  if (x_test > xtrial) return;     // above upper line
  if (x_test < -1.0*xtrial) return;// below lower line
  //return 1;  // inside main acceptance.

#ifdef ACCTEST
  cout << "col ";
#endif 

  acc_bool[9] = 1;

  //Target to Q1 exit, circle of radius 14.92 cm
  x_test = x_sp_q1ex_(vector_jjl, ii)*m2cm;
  y_test = y_sp_q1ex_(vector_jjl, ii)*m2cm;
  //x_test = x_test + 0.9;
  //cout << x_test*x_test << " " <<  y_test*y_test << " " << 14.92*14.92 << endl;
  if( (x_test*x_test + y_test*y_test) > (14.92*14.92) )
    return;

#ifdef ACCTEST
  cout << "q1ex ";
#endif 

  acc_bool[0] = 1;

  //cout << "q1ex" << endl;
  //Target to dipole entrance, trapezoid -522.0cm<x<-498.1cm  |y| < -0.1924*x-19.24
  x_test = x_sp_den_(vector_jjl, ii);
  y_test = y_sp_den_(vector_jjl, ii);
  if( (x_test<-5.220) || (x_test>-4.981) || ( y_test < -(-0.1924*x_test-.1924) ) || ( y_test > (-0.1924*x_test-.1924) ) )
    return;

#ifdef ACCTEST
  cout << "den ";
#endif 

  acc_bool[2] = 1;
  
  //Target to dipole exit, trapezoid -46.19cm<x<46.19cm  |y| < -0.0161*x+12.5
  x_test = x_sp_dex_(vector_jjl, ii)*m2cm;
  y_test = y_sp_dex_(vector_jjl, ii)*m2cm;
  if( (x_test<-46.19) || (x_test>46.19) || fabs(y_test) > fabs(-0.0161*x_test+12.5) )
    return;
  
#ifdef ACCTEST
  cout << "dex ";
#endif 

  acc_bool[3] = 1;

  //cout << "The fox ";
  x_test = x_sp_q2ex_(vector_jjl, ii)*m2cm;
  y_test = y_sp_q2ex_(vector_jjl, ii)*m2cm;
  if( (x_test*x_test + y_test*y_test) > (30.0*30.0) )
    return;

#ifdef ACCTEST
  cout << "q2ex ";
#endif 

  acc_bool[1] = 1;
  //cout << "is mad!" << endl;

  //cout << "dex" << endl;
  //Target to Q3 entrance, circle of radius 30.0 cm
  x_test = x_sp_q3en_(vector_jjl, ii)*m2cm;
  y_test = y_sp_q3en_(vector_jjl, ii)*m2cm;
  if( (x_test*x_test + y_test*y_test) > (30.0*30.0) )
    return;

#ifdef ACCTEST
  cout << "q3en ";
#endif 

  acc_bool[4] = 1;
  //cout << "q3en" << endl;
  //Target to Q3 exit, circle of radius 30.0 cm  -> 28.0cm
  x_test = x_sp_q3ex_(vector_jjl, ii)*m2cm;
  y_test = y_sp_q3ex_(vector_jjl, ii)*m2cm;
  //x_test = (x_test - 1.0) / (28.0);
  //y_test = y_test / (30.0);
  if( (x_test*x_test + y_test*y_test) > (30.0*30.0) )
    return;

#ifdef ACCTEST
  cout << "q3ex ";
#endif 
  
  acc_bool[5] = 1;
  //cout << "q3ex" << endl;
  /////////////////////////////////////////////////////////////
  /* If we reach this point, it means the test was succesful */

#ifdef ACCTEST
  cout << "made it!" << endl;
#endif 
  
  return;
}

void hamcPREXTrans::Acceptance_C(double* pV5, double* x_check, double* theta_check, double* y_check, double* phi_check, int* acc_bool){
  //cout << "Checking the acceptance..., but C mode" << endl;
  float vector_jjl[]={(float)pV5[0],(float)pV5[1],(float)pV5[2],(float)pV5[3],(float)pV5[4]};
  int iii = 5; int *ii = &iii;

  x_check[0]      = x_sp_cq1x_(vector_jjl, ii);   y_check[0]    = y_sp_cq1x_(vector_jjl, ii);
  theta_check[0]  = t_sp_cq1x_(vector_jjl, ii);   phi_check[0]  = p_sp_cq1x_(vector_jjl, ii);
  x_check[1]      = 0.;                           y_check[1]    = 0.;
  theta_check[1]  = 0.;                           phi_check[1]  = 0.;
  x_check[2]      = x_sp_cden_(vector_jjl, ii);   y_check[2]    = y_sp_cden_(vector_jjl, ii);
  theta_check[2]  = t_sp_cden_(vector_jjl, ii);   phi_check[2]  = p_sp_cden_(vector_jjl, ii);
  x_check[3]      = x_sp_cdex_(vector_jjl, ii);   y_check[3]    = y_sp_cdex_(vector_jjl, ii);
  theta_check[3]  = t_sp_cdex_(vector_jjl, ii);   phi_check[3]  = p_sp_cdex_(vector_jjl, ii);
  x_check[4]      = x_sp_cq3e_(vector_jjl, ii);   y_check[4]    = y_sp_cq3e_(vector_jjl, ii);
  theta_check[4]  = t_sp_cq3e_(vector_jjl, ii);   phi_check[4]  = p_sp_cq3e_(vector_jjl, ii);
  x_check[5]      = x_sp_cq3x_(vector_jjl, ii);   y_check[5]    = y_sp_cq3x_(vector_jjl, ii);
  theta_check[5]  = t_sp_cq3x_(vector_jjl, ii);   phi_check[5]  = p_sp_cq3x_(vector_jjl, ii);
  x_check[6]      = 0.;                           y_check[6]    = 0.;
  theta_check[6]  = 0.;                           phi_check[6]  = 0.;
  x_check[7]      = 0.;                           y_check[7]    = 0.;
  theta_check[7]  = 0.;                           phi_check[7]  = 0.;
  x_check[8]      = 0.;                           y_check[8]    = 0.;
  theta_check[8]  = 0.;                           phi_check[8]  = 0.;
  x_check[9]      = 0.;                           y_check[9]    = 0.;
  theta_check[9]  = 0.;                           phi_check[9]  = 0.;
  x_check[10]     = 0.;                           y_check[10]   = 0.;
  theta_check[10] = 0.;                           phi_check[10] = 0.;

  float x_test, y_test;

  //Nickie adds the septum cuts
  //sen: +0.088 < x < +0.382
  //-0.120 < y < 0.120
  //sm: +0.088 < x < +0.382
  //-0.120 < y < 0.120
  //sex: +0.088 < x < +0.382
  //-0.120 < y < 0.120  

  //x_test = x_sp_sen_(vector_jjl, ii);
  //y_test = y_sp_sen_(vector_jjl, ii);
  //if( x_test < 0.088 || x_test > 0.382 || y_test < -0.120 || y_test > 0.120 )
  //return;

#ifdef ACCTEST
  cout << "sen ";
#endif 

  acc_bool[6] = 1;

  //x_test = x_sp_sm_(vector_jjl, ii);
  //y_test = y_sp_sm_(vector_jjl, ii);
  //if( x_test < 0.088 || x_test > 0.382 || y_test < -0.120 || y_test > 0.120 )
  //return;

#ifdef ACCTEST
  cout << "sm ";
#endif 

  acc_bool[7] = 1;

  //x_test = x_sp_sex_(vector_jjl, ii);
  //y_test = y_sp_sex_(vector_jjl, ii);
  //if( x_test < 0.088 || x_test > 0.382 || y_test < -0.120 || y_test > 0.120 )
  //return;

#ifdef ACCTEST
  cout << "sex ";
#endif 

  acc_bool[8] = 1;

  //Nickie adds the collimator:
  x_test = vector_jjl[0];
  y_test = vector_jjl[2];
  float atxlo = 0.032;
  float atyhi = 0.041;
  float atr   = 0.20;
  float atc   = 0.229;
  float rad1  = 0.205;
  float ycent1= 0.145;
  float rad2  = 0.20;
  float ycent2= 0.229;
  float vtop  = 0.117;
  float vright= 0.04;
  float yline1= 0.1474;
  float slope1= -1.88;
  float ysign = 1.0;
  
  float xdif = x_test;
  float ydif = ysign*(y_test);
  float xtrial,ytrial,xtmp;
  
  ytrial = ydif-atc;
  xtrial = atr*atr - ytrial*ytrial;
  if (xtrial > 0) {
    xtmp = sqrt(xtrial);
    if ((xdif > atxlo) && (ydif < atyhi) &&
	(x_test < xtmp)) return;     // in A_T hole (a triangle with arc)
  }
  
 // Now check if in main acceptance
  // Remember x is vertical (transport), y is horizontal
  ytrial = ydif-ycent1;
  xtrial = rad1*rad1 - ytrial*ytrial;
  if (xtrial < 0) return;   // outside outer circle
  xtrial = sqrt(xtrial);
  if (x_test > xtrial || x_test < -1.0*xtrial) return;  // outside outer circle
  ytrial = ydif-ycent2;
  xtrial = rad2*rad2 - ytrial*ytrial;
  if (xtrial > 0) {
    xtrial = sqrt(xtrial);
    
    if ((x_test > 0 && x_test < xtrial) ||
	(x_test < 0 && x_test > -1.0*xtrial)) return; // outside innner circle
  }
  if (x_test > vtop) return;       // above upper line
  if (x_test < -1.0*vtop) return;  // beneath lower line
  if (ydif > vright) return;     // beyond right border
  xtrial = slope1*ydif + yline1; // Champhor line
  if (x_test > xtrial) return;     // above upper line
  if (x_test < -1.0*xtrial) return;// below lower line
  //return 1;  // inside main acceptance.
  
#ifdef ACCTEST
  cout << "col ";
#endif 

  acc_bool[9] = 1;

  //Target to Q1 exit, circle of radius 14.92 cm
  x_test = x_sp_cq1x_(vector_jjl, ii)*m2cm;
  y_test = y_sp_cq1x_(vector_jjl, ii)*m2cm;
  //x_test = x_test + 0.9;
  //cout << x_test*x_test << " " <<  y_test*y_test << " " << 14.92*14.92 << endl;
  if( (x_test*x_test + y_test*y_test) > (14.92*14.92) )
    return;

#ifdef ACCTEST
  cout << "q1ex ";
#endif 

  acc_bool[0] = 1;

  //cout << "q1ex" << endl;
  //Target to dipole entrance, trapezoid -522.0cm<x<-498.1cm  |y| < -0.1924*x-19.24
  x_test = x_sp_cden_(vector_jjl, ii);
  y_test = y_sp_cden_(vector_jjl, ii);
  if( (x_test<-5.220) || (x_test>-4.981) || ( y_test < -(-0.1924*x_test-.1924) ) || ( y_test > (-0.1924*x_test-.1924) ) )
    return;

#ifdef ACCTEST
  cout << "den ";
#endif 

  acc_bool[2] = 1;
  
  //Target to dipole exit, trapezoid -46.19cm<x<46.19cm  |y| < -0.0161*x+12.5
  x_test = x_sp_cdex_(vector_jjl, ii)*m2cm;
  y_test = y_sp_cdex_(vector_jjl, ii)*m2cm;
  if( (x_test<-46.19) || (x_test>46.19) || fabs(y_test) > fabs(-0.0161*x_test+12.5) )
    return;
  
#ifdef ACCTEST
  cout << "dex ";
#endif 

  acc_bool[3] = 1;

  //cout << "The fox ";
  //x_test = x_sp_cq2x_(vector_jjl, ii)*m2cm;
  //y_test = y_sp_cq2x_(vector_jjl, ii)*m2cm;
  //if( (x_test*x_test + y_test*y_test) > (30.0*30.0) )
  //return;

#ifdef ACCTEST
  cout << "q2ex ";
#endif 

  acc_bool[1] = 1;
  //cout << "is mad!" << endl;

  //cout << "dex" << endl;
  //Target to Q3 entrance, circle of radius 30.0 cm
  x_test = x_sp_cq3e_(vector_jjl, ii)*m2cm;
  y_test = y_sp_cq3e_(vector_jjl, ii)*m2cm;
  if( (x_test*x_test + y_test*y_test) > (30.0*30.0) )
    return;

#ifdef ACCTEST
  cout << "q3en ";
#endif 

  acc_bool[4] = 1;
  //cout << "q3en" << endl;
  //Target to Q3 exit, circle of radius 30.0 cm  -> 28.0cm
  x_test = x_sp_cq3x_(vector_jjl, ii)*m2cm;
  y_test = y_sp_cq3x_(vector_jjl, ii)*m2cm;
  //x_test = (x_test - 1.0) / (28.0);
  //y_test = y_test / (30.0);
  if( (x_test*x_test + y_test*y_test) > (30.0*30.0) )
    return;

#ifdef ACCTEST
  cout << "q3ex ";
#endif 
  
  acc_bool[5] = 1;
  //cout << "q3ex" << endl;
  /////////////////////////////////////////////////////////////
  /* If we reach this point, it means the test was succesful */

#ifdef ACCTEST
  cout << "made it!" << endl;
#endif 
  
  return;
}
