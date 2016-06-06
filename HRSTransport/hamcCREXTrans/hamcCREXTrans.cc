////////////////////////////////////////////////////////////////////////
// Standard HRS transport functions FOR CREX
// 
////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "crex_4degr.hh"
#include "hamcCREXTrans.hh"
#include <iostream>

//#define NICKIETEST 1
//#define ACCTEST 1
//using namespace SNoSepta;
using namespace std;

const float m2cm = 100.0;
const double kDEG = 3.14159265358979323846/180.0;

hamcCREXTrans::hamcCREXTrans()//:cModelAngle(12.5*kDEG)
:cModelAngle(5.0*kDEG)
{
    // Nothing to do
}

hamcCREXTrans::~hamcCREXTrans()
{
    // Nothing to do
}

bool hamcCREXTrans::TransLeftHRS(double* pV5)
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
	//x_test = x_s4_q1en_(vector_jjl, ii)*m2cm;
        //y_test = y_s4_q1en_(vector_jjl, ii)*m2cm;
        //x_test = x_test - 0.9;
        //if( (x_test*x_test + y_test*y_test) > (14.92*14.92) )
	//return false;
	
	//Target to Q1 exit, circle of radius 14.92 cm
	x_test = x_s4_q1ex_(vector_jjl, ii)*m2cm;
	y_test = y_s4_q1ex_(vector_jjl, ii)*m2cm;
	//x_test = x_test - 0.9;
	if( (x_test*x_test + y_test*y_test) > (14.92*14.92) )
		return false;

	//Nickie: this used to be commented
	//Target to dipole entrance, trapezoid 522.0cm>x>498.1cm  |y| < -0.1924*x-19.24
	//x_test = x_s4_den_(vector_jjl, ii)*m2cm;
	//y_test = y_s4_den_(vector_jjl, ii)*m2cm;
	//if( (x_test>522.0) || (x_test<498.1) || fabs(y_test) > fabs(-0.1924*x_test-19.24) )
	//return false;

	//Target to dipole exit, trapezoid -46.19cm<x<46.19cm  |y| < -0.0161*x+12.5
	x_test = x_s4_dex_(vector_jjl, ii)*m2cm;
	y_test = y_s4_dex_(vector_jjl, ii)*m2cm;
	//cout<<"dipole_exit:(x,y)=\t"<<x_test<<"\t "<<y_test<<endl;
	if( fabs(x_test)>46.19 || fabs(y_test) > fabs(-0.0161*x_test+12.5) )
		return false;

	//Target to Q3 entrance, circle of radius 30.0 cm
	x_test = x_s4_q3en_(vector_jjl, ii)*m2cm;
	y_test = y_s4_q3en_(vector_jjl, ii)*m2cm;
	if( (x_test*x_test + y_test*y_test) > (30.0*30.0) )
		return false;

	//Target to Q3 exit, circle of radius 30.0 cm  -> 28.0cm
	x_test = x_s4_q3ex_(vector_jjl, ii)*m2cm;
	y_test = y_s4_q3ex_(vector_jjl, ii)*m2cm;
	//x_test = (x_test + 1.0) / (28.0);
	//y_test = y_test / (30.0);
	if( (x_test*x_test + y_test*y_test) > (30.0*30.0) )
		return false;

	/////////////////////////////////////////////////////////////
	/* If we reach this point, it means the test was succesful */
	float x_fp     = x_s4_fp_(vector_jjl, ii);
	float theta_fp = t_s4_fp_(vector_jjl, ii);
	float y_fp     = y_s4_fp_(vector_jjl, ii);
	float phi_fp   = p_s4_fp_(vector_jjl, ii);

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

//bool hamcCREXTrans::TransRightHRS(double* pV5){
//return false;
//}

bool hamcCREXTrans::TransRightHRS(double* pV5)
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
	x_test = x_s4_sext_(vector_jjl, ii);
	y_test = y_s4_sext_(vector_jjl, ii);
	if( x_test < 0.088 || x_test > 0.382 || y_test < -0.120 || y_test > 0.120 )
	  return false;

	x_test = x_s4_q1en_(vector_jjl, ii)*m2cm;
	y_test = y_s4_q1en_(vector_jjl, ii)*m2cm;
	if( ( x_test < -6.09 ) || ( x_test > 6.09 ) || ( y_test < -3.145 ) || ( y_test > 3.145 ) )
	  return false;//collimator 1?
	
	//Target to Q1 exit, circle of radius 14.92 cm
	x_test = x_s4_q1ex_(vector_jjl, ii)*m2cm;
	y_test = y_s4_q1ex_(vector_jjl, ii)*m2cm;
	//x_test = x_test + 0.9;
	//cout << "Man, it didn't even pass the first cut!!" << endl;
	//cout << x_test*x_test << " " <<  y_test*y_test << " " << 14.92*14.92 << endl;
	if( (x_test*x_test + y_test*y_test) > (14.92*14.92) )
	  return false;


	//cout << "q1ex" << endl;
	//Target to dipole entrance, trapezoid -522.0cm<x<-498.1cm  |y| < -0.1924*x-19.24
	x_test = x_s4_den_(vector_jjl, ii);
	y_test = y_s4_den_(vector_jjl, ii);
	//cout << "x and y: " << x_test << " " << y_test << endl;
	//cout << -6.19 << " < x_fail < " << -5.9 << endl;
	//cout << -(-0.1924*x_test-.1924) << " > y_fail > " << (-0.1924*x_test-.1924) << endl;
	if( (x_test<-6.19) || (x_test>-5.9) || ( y_test < -(-0.1924*x_test-.161) ) || ( y_test > (-0.1924*x_test-.161) ) )
	  return false;
	//cout << "den" << endl;
	//Target to dipole exit, trapezoid -46.19cm<x<46.19cm  |y| < -0.0161*x+12.5
	x_test = x_s4_dex_(vector_jjl, ii)*m2cm;
	y_test = y_s4_dex_(vector_jjl, ii)*m2cm;
	if( (x_test<-46.19) || (x_test>46.19) || fabs(y_test) > fabs(-0.0161*x_test+12.5) )
		return false;
	//cout << "dex" << endl;
	//Target to Q3 entrance, circle of radius 30.0 cm
	x_test = x_s4_q3en_(vector_jjl, ii)*m2cm;
	y_test = y_s4_q3en_(vector_jjl, ii)*m2cm;
	if( (x_test*x_test + y_test*y_test) > (30.0*30.0) )
		return false;
	//cout << "q3en" << endl;
	//Target to Q3 exit, circle of radius 30.0 cm  -> 28.0cm
	x_test = x_s4_q3ex_(vector_jjl, ii)*m2cm;
	y_test = y_s4_q3ex_(vector_jjl, ii)*m2cm;
	//x_test = (x_test - 1.0) / (28.0);
	//y_test = y_test / (30.0);
	if( (x_test*x_test + y_test*y_test) > (30.0*30.0) )
		return false;
	//cout << "q3ex" << endl;
	/////////////////////////////////////////////////////////////
	/* If we reach this point, it means the test was succesful */

	float x_fp     = x_s4_fp_(vector_jjl, ii);
	float theta_fp = t_s4_fp_(vector_jjl, ii);
	float y_fp     = y_s4_fp_(vector_jjl, ii);
	float phi_fp   = p_s4_fp_(vector_jjl, ii);
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

bool hamcCREXTrans::TransRightHRS_C(double* pV5)//NOT FINISHED
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
	x_test = x_s4_sext_(vector_jjl, ii);
	y_test = y_s4_sext_(vector_jjl, ii);
	if( x_test < 0.088 || x_test > 0.382 || y_test < -0.120 || y_test > 0.120 )
	  return false;

	x_test = x_s4_q1en_(vector_jjl, ii)*m2cm;
	y_test = y_s4_q1en_(vector_jjl, ii)*m2cm;
	if( ( x_test < -6.09 ) || ( x_test > 6.09 ) || ( y_test < -3.145 ) || ( y_test > 3.145 ) )
	  return false;//collimator 1?
	
	//Target to Q1 exit, circle of radius 14.92 cm
	x_test = x_s4_q1ex_(vector_jjl, ii)*m2cm;
	y_test = y_s4_q1ex_(vector_jjl, ii)*m2cm;
	//x_test = x_test + 0.9;
	//cout << "Man, it didn't even pass the first cut!!" << endl;
	//cout << x_test*x_test << " " <<  y_test*y_test << " " << 14.92*14.92 << endl;
	if( (x_test*x_test + y_test*y_test) > (14.92*14.92) )
	  return false;


	//cout << "q1ex" << endl;
	//Target to dipole entrance, trapezoid -522.0cm<x<-498.1cm  |y| < -0.1924*x-19.24
	x_test = x_s4_den_(vector_jjl, ii);
	y_test = y_s4_den_(vector_jjl, ii);
	//cout << "x and y: " << x_test << " " << y_test << endl;
	//cout << -6.19 << " < x_fail < " << -5.9 << endl;
	//cout << -(-0.1924*x_test-.1924) << " > y_fail > " << (-0.1924*x_test-.1924) << endl;
	if( (x_test<-6.19) || (x_test>-5.9) || ( y_test < -(-0.1924*x_test-.161) ) || ( y_test > (-0.1924*x_test-.161) ) )
	  return false;
	//cout << "den" << endl;
	//Target to dipole exit, trapezoid -46.19cm<x<46.19cm  |y| < -0.0161*x+12.5
	x_test = x_s4_dex_(vector_jjl, ii)*m2cm;
	y_test = y_s4_dex_(vector_jjl, ii)*m2cm;
	if( (x_test<-46.19) || (x_test>46.19) || fabs(y_test) > fabs(-0.0161*x_test+12.5) )
		return false;
	//cout << "dex" << endl;
	//Target to Q3 entrance, circle of radius 30.0 cm
	x_test = x_s4_q3en_(vector_jjl, ii)*m2cm;
	y_test = y_s4_q3en_(vector_jjl, ii)*m2cm;
	if( (x_test*x_test + y_test*y_test) > (30.0*30.0) )
		return false;
	//cout << "q3en" << endl;
	//Target to Q3 exit, circle of radius 30.0 cm  -> 28.0cm
	x_test = x_s4_q3ex_(vector_jjl, ii)*m2cm;
	y_test = y_s4_q3ex_(vector_jjl, ii)*m2cm;
	//x_test = (x_test - 1.0) / (28.0);
	//y_test = y_test / (30.0);
	if( (x_test*x_test + y_test*y_test) > (30.0*30.0) )
		return false;
	//cout << "q3ex" << endl;
	/////////////////////////////////////////////////////////////
	/* If we reach this point, it means the test was succesful */

	float x_fp     = x_s4_fp_(vector_jjl, ii);
	float theta_fp = t_s4_fp_(vector_jjl, ii);
	float y_fp     = y_s4_fp_(vector_jjl, ii);
	float phi_fp   = p_s4_fp_(vector_jjl, ii);
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


void hamcCREXTrans::ReconLeftHRS(double* pV5){
  return;
}
/*
{   
	float vector_jjl[]={pV5[0],pV5[1],pV5[2],pV5[3],pV5[4]};
	int iii = 5; int *ii = &iii;
    int jjj = 1; int *jj = &jjj;
    
	vector_jjl[1] = vector_jjl[1] - txfit_s4_(vector_jjl, jj);
	float x_or      = vector_jjl[4];
	float delta_rec = delta_s4_(vector_jjl, ii);
	float theta_rec = theta_s4_(vector_jjl, ii);
	float phi_rec   = phi_s4_(vector_jjl, ii); 
	float y_rec     = y00_s4_(vector_jjl, ii); 

	//reset the vector and return it back to the caller
	pV5[0] = (double)x_or;
	pV5[1] = (double)theta_rec;
	pV5[2] = (double)y_rec;
	pV5[3] = (double)phi_rec;
	pV5[4] = (double)delta_rec;
}
*/

void hamcCREXTrans::ReconRightHRS(double* pV5){
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

void hamcCREXTrans::Acceptance(double* pV5, double* x_check, double* theta_check, double* y_check, double* phi_check, int* acc_bool){
  //cout << "Checking the acceptance..." << endl;
  float vector_jjl[]={(float)pV5[0],(float)pV5[1],(float)pV5[2],(float)pV5[3],(float)pV5[4]};
  int iii = 5; int *ii = &iii;

  x_check[0]  = x_s4_q1ex_(vector_jjl, ii);  theta_check[0]  = t_s4_q1ex_(vector_jjl, ii);
  y_check[0]  = p_s4_q1ex_(vector_jjl, ii);  phi_check[0]    = p_s4_q1ex_(vector_jjl, ii);
  x_check[1]  = 0.;                          theta_check[1]  = 0.;
  y_check[1]  = 0.;                          phi_check[1]    = 0.;
  x_check[2]  = x_s4_den_(vector_jjl, ii);   theta_check[2]  = t_s4_den_(vector_jjl, ii);
  y_check[2]  = y_s4_den_(vector_jjl, ii);   phi_check[2]    = p_s4_den_(vector_jjl, ii);
  x_check[3]  = x_s4_dex_(vector_jjl, ii);   theta_check[3]  = t_s4_dex_(vector_jjl, ii);
  y_check[3]  = y_s4_dex_(vector_jjl, ii);   phi_check[3]    = p_s4_dex_(vector_jjl, ii);
  x_check[4]  = x_s4_q3en_(vector_jjl, ii);  theta_check[4]  = t_s4_q3en_(vector_jjl, ii);
  y_check[4]  = y_s4_q3en_(vector_jjl, ii);  phi_check[4]    = p_s4_q3en_(vector_jjl, ii);
  x_check[5]  = x_s4_q3ex_(vector_jjl, ii);  theta_check[5]  = t_s4_q3ex_(vector_jjl, ii);
  y_check[5]  = y_s4_q3ex_(vector_jjl, ii);  phi_check[5]    = p_s4_q3ex_(vector_jjl, ii);
  x_check[6]  = 0.;                          theta_check[6]  = 0.;
  y_check[6]  = 0.;                          phi_check[6]    = 0.;
  x_check[7]  = 0.;                          theta_check[7]  = 0.;
  y_check[7]  = 0.;                          phi_check[7]    = 0.;
  x_check[8]  = x_s4_sext_(vector_jjl, ii);  theta_check[8]  = t_s4_sext_(vector_jjl, ii);
  y_check[8]  = y_s4_sext_(vector_jjl, ii);  phi_check[8]    = p_s4_sext_(vector_jjl, ii);
  x_check[9]  = 0.;                          theta_check[9]  = 0.;
  y_check[9]  = 0.;                          phi_check[9]    = 0.;
  x_check[10] = x_s4_q1en_(vector_jjl, ii);  theta_check[10] = t_s4_q1en_(vector_jjl, ii);
  y_check[10] = y_s4_q1en_(vector_jjl, ii);  phi_check[10]   = p_s4_q1en_(vector_jjl, ii);

  float x_test, y_test;

  //Nickie adds the septum cuts
  //sen: +0.088 < x < +0.382
  //-0.120 < y < 0.120
  //sm: +0.088 < x < +0.382
  //-0.120 < y < 0.120
  //sex: +0.088 < x < +0.382
  //-0.120 < y < 0.120  

  //x_test = x_s4_sen_(vector_jjl, ii);
  //y_test = y_s4_sen_(vector_jjl, ii);
  //if( x_test < 0.088 || x_test > 0.382 || y_test < -0.120 || y_test > 0.120 )
  //return;

#ifdef ACCTEST
  cout << "sen ";
#endif 

  acc_bool[6] = 1;

  //x_test = x_s4_sm_(vector_jjl, ii);
  //y_test = y_s4_sm_(vector_jjl, ii);
  //if( x_test < 0.088 || x_test > 0.382 || y_test < -0.120 || y_test > 0.120 )
  //return;

#ifdef ACCTEST
  cout << "sm ";
#endif 

  acc_bool[7] = 1;

  x_test = x_s4_sext_(vector_jjl, ii);
  y_test = y_s4_sext_(vector_jjl, ii);
  if( x_test < 0.088 || x_test > 0.382 || y_test < -0.120 || y_test > 0.120 )
    return;

#ifdef ACCTEST
  cout << "sex ";
#endif 

  acc_bool[8] = 1;

#ifdef ACCTEST
  cout << "col ";
#endif 
  
  acc_bool[9] = 1;

  x_test = x_s4_q1en_(vector_jjl, ii)*m2cm;
  y_test = y_s4_q1en_(vector_jjl, ii)*m2cm;
  if( ( x_test < -6.09 ) || ( x_test > 6.09 ) || ( y_test < -3.145 ) || ( y_test > 3.145 ) )
    return;//collimator 1?
  
  //Target to Q1 exit, circle of radius 14.92 cm
  x_test = x_s4_q1ex_(vector_jjl, ii)*m2cm;
  y_test = y_s4_q1ex_(vector_jjl, ii)*m2cm;
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
  x_test = x_s4_den_(vector_jjl, ii);
  y_test = y_s4_den_(vector_jjl, ii);
  if( (x_test<-6.19) || (x_test>-5.9) || ( y_test < -(-0.1924*x_test-.161) ) || ( y_test > (-0.1924*x_test-.161) ) )
    return;

#ifdef ACCTEST
  cout << "den ";
#endif 

  acc_bool[2] = 1;
  
  //Target to dipole exit, trapezoid -46.19cm<x<46.19cm  |y| < -0.0161*x+12.5
  x_test = x_s4_dex_(vector_jjl, ii)*m2cm;
  y_test = y_s4_dex_(vector_jjl, ii)*m2cm;
  if( (x_test<-46.19) || (x_test>46.19) || fabs(y_test) > fabs(-0.0161*x_test+12.5) )
    return;
  
#ifdef ACCTEST
  cout << "dex ";
#endif 

  acc_bool[3] = 1;

  //cout << "The fox ";
  //x_test = x_s4_q2ex_(vector_jjl, ii)*m2cm;
  //y_test = y_s4_q2ex_(vector_jjl, ii)*m2cm;
  //if( (x_test*x_test + y_test*y_test) > (30.0*30.0) )
  //return;

#ifdef ACCTEST
  cout << "q2ex ";
#endif 

  acc_bool[1] = 1;
  //cout << "is mad!" << endl;

  //cout << "dex" << endl;
  //Target to Q3 entrance, circle of radius 30.0 cm
  x_test = x_s4_q3en_(vector_jjl, ii)*m2cm;
  y_test = y_s4_q3en_(vector_jjl, ii)*m2cm;
  if( (x_test*x_test + y_test*y_test) > (30.0*30.0) )
    return;

#ifdef ACCTEST
  cout << "q3en ";
#endif 

  acc_bool[4] = 1;
  //cout << "q3en" << endl;
  //Target to Q3 exit, circle of radius 30.0 cm  -> 28.0cm
  x_test = x_s4_q3ex_(vector_jjl, ii)*m2cm;
  y_test = y_s4_q3ex_(vector_jjl, ii)*m2cm;
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

void hamcCREXTrans::Acceptance_C(double* pV5, double* x_check, double* theta_check, double* y_check, double* phi_check, int* acc_bool){
  //cout << "Checking the acceptance..." << endl;
  float vector_jjl[]={(float)pV5[0],(float)pV5[1],(float)pV5[2],(float)pV5[3],(float)pV5[4]};
  int iii = 5; int *ii = &iii;

  x_check[0]  = x_s4_q1ex_(vector_jjl, ii);  theta_check[0]  = t_s4_q1ex_(vector_jjl, ii);
  y_check[0]  = p_s4_q1ex_(vector_jjl, ii);  phi_check[0]    = p_s4_q1ex_(vector_jjl, ii);
  x_check[1]  = 0.;                          theta_check[1]  = 0.;
  y_check[1]  = 0.;                          phi_check[1]    = 0.;
  x_check[2]  = x_s4_den_(vector_jjl, ii);   theta_check[2]  = t_s4_den_(vector_jjl, ii);
  y_check[2]  = y_s4_den_(vector_jjl, ii);   phi_check[2]    = p_s4_den_(vector_jjl, ii);
  x_check[3]  = x_s4_dex_(vector_jjl, ii);   theta_check[3]  = t_s4_dex_(vector_jjl, ii);
  y_check[3]  = y_s4_dex_(vector_jjl, ii);   phi_check[3]    = p_s4_dex_(vector_jjl, ii);
  x_check[4]  = x_s4_q3en_(vector_jjl, ii);  theta_check[4]  = t_s4_q3en_(vector_jjl, ii);
  y_check[4]  = y_s4_q3en_(vector_jjl, ii);  phi_check[4]    = p_s4_q3en_(vector_jjl, ii);
  x_check[5]  = x_s4_q3ex_(vector_jjl, ii);  theta_check[5]  = t_s4_q3ex_(vector_jjl, ii);
  y_check[5]  = y_s4_q3ex_(vector_jjl, ii);  phi_check[5]    = p_s4_q3ex_(vector_jjl, ii);
  x_check[6]  = 0.;                          theta_check[6]  = 0.;
  y_check[6]  = 0.;                          phi_check[6]    = 0.;
  x_check[7]  = 0.;                          theta_check[7]  = 0.;
  y_check[7]  = 0.;                          phi_check[7]    = 0.;
  x_check[8]  = x_s4_sext_(vector_jjl, ii);  theta_check[8]  = t_s4_sext_(vector_jjl, ii);
  y_check[8]  = y_s4_sext_(vector_jjl, ii);  phi_check[8]    = p_s4_sext_(vector_jjl, ii);
  x_check[9]  = 0.;                          theta_check[9]  = 0.;
  y_check[9]  = 0.;                          phi_check[9]    = 0.;
  x_check[10] = x_s4_q1en_(vector_jjl, ii);  theta_check[10] = t_s4_q1en_(vector_jjl, ii);
  y_check[10] = y_s4_q1en_(vector_jjl, ii);  phi_check[10]   = p_s4_q1en_(vector_jjl, ii);


  float x_test, y_test;

  //Nickie adds the septum cuts
  //sen: +0.088 < x < +0.382
  //-0.120 < y < 0.120
  //sm: +0.088 < x < +0.382
  //-0.120 < y < 0.120
  //sex: +0.088 < x < +0.382
  //-0.120 < y < 0.120  

  //x_test = x_s4_sen_(vector_jjl, ii);
  //y_test = y_s4_sen_(vector_jjl, ii);
  //if( x_test < 0.088 || x_test > 0.382 || y_test < -0.120 || y_test > 0.120 )
  //return;

#ifdef ACCTEST
  cout << "sen ";
#endif 

  acc_bool[6] = 1;

  //x_test = x_s4_sm_(vector_jjl, ii);
  //y_test = y_s4_sm_(vector_jjl, ii);
  //if( x_test < 0.088 || x_test > 0.382 || y_test < -0.120 || y_test > 0.120 )
  //return;

#ifdef ACCTEST
  cout << "sm ";
#endif 

  acc_bool[7] = 1;

  x_test = x_s4_sext_(vector_jjl, ii);
  y_test = y_s4_sext_(vector_jjl, ii);
  if( x_test < 0.088 || x_test > 0.382 || y_test < -0.120 || y_test > 0.120 )
    return;

#ifdef ACCTEST
  cout << "sex ";
#endif 

  acc_bool[8] = 1;

#ifdef ACCTEST
  cout << "col ";
#endif 
  
  acc_bool[9] = 1;

  x_test = x_s4_q1en_(vector_jjl, ii)*m2cm;
  y_test = y_s4_q1en_(vector_jjl, ii)*m2cm;
  if( ( x_test < -6.09 ) || ( x_test > 6.09 ) || ( y_test < -3.145 ) || ( y_test > 3.145 ) )
    return;//collimator 1?
  
  //Target to Q1 exit, circle of radius 14.92 cm
  x_test = x_s4_q1ex_(vector_jjl, ii)*m2cm;
  y_test = y_s4_q1ex_(vector_jjl, ii)*m2cm;
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
  x_test = x_s4_den_(vector_jjl, ii);
  y_test = y_s4_den_(vector_jjl, ii);
  if( (x_test<-6.19) || (x_test>-5.9) || ( y_test < -(-0.1924*x_test-.161) ) || ( y_test > (-0.1924*x_test-.161) ) )
    return;

#ifdef ACCTEST
  cout << "den ";
#endif 

  acc_bool[2] = 1;
  
  //Target to dipole exit, trapezoid -46.19cm<x<46.19cm  |y| < -0.0161*x+12.5
  x_test = x_s4_dex_(vector_jjl, ii)*m2cm;
  y_test = y_s4_dex_(vector_jjl, ii)*m2cm;
  if( (x_test<-46.19) || (x_test>46.19) || fabs(y_test) > fabs(-0.0161*x_test+12.5) )
    return;
  
#ifdef ACCTEST
  cout << "dex ";
#endif 

  acc_bool[3] = 1;

  //cout << "The fox ";
  //x_test = x_s4_q2ex_(vector_jjl, ii)*m2cm;
  //y_test = y_s4_q2ex_(vector_jjl, ii)*m2cm;
  //if( (x_test*x_test + y_test*y_test) > (30.0*30.0) )
  //return;

#ifdef ACCTEST
  cout << "q2ex ";
#endif 

  acc_bool[1] = 1;
  //cout << "is mad!" << endl;

  //cout << "dex" << endl;
  //Target to Q3 entrance, circle of radius 30.0 cm
  x_test = x_s4_q3en_(vector_jjl, ii)*m2cm;
  y_test = y_s4_q3en_(vector_jjl, ii)*m2cm;
  if( (x_test*x_test + y_test*y_test) > (30.0*30.0) )
    return;

#ifdef ACCTEST
  cout << "q3en ";
#endif 

  acc_bool[4] = 1;
  //cout << "q3en" << endl;
  //Target to Q3 exit, circle of radius 30.0 cm  -> 28.0cm
  x_test = x_s4_q3ex_(vector_jjl, ii)*m2cm;
  y_test = y_s4_q3ex_(vector_jjl, ii)*m2cm;
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
