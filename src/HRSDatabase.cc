#include "HRSDatabase.hh"
#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
using namespace std;

HRSDatabase::HRSDatabase(int mymodel) : model(mymodel) {
  //I'm not giving the user the option of loading the table herself.
  if( model == 0 ) {
    ene_min = 0.6;
    ene_step = 0.05;
    imax = 14;
    jmax = 66;
    kmax = 6;
    LoadTable("horpb.dat" , 0);
    LoadTable("horpb1.dat", 1);
  } else if( model == 1 ) {
    ene_min = 0.5;
    ene_step = 0.05;
    imax = 62;
    jmax = 141;
    kmax = 3;
    LoadTable("ca48_fsu.dat" , 0);
    LoadTable("ca48_fsu_stretched.dat", 1);
  }
}

/////////////////////////////////////////////////////////////////////////

HRSDatabase::~HRSDatabase(){
}

////////////////////////////////////////////////////////////////////////

void HRSDatabase::LoadTable(string filename, int index){
  //indices:
  //1. Pb Horowitz
  //2. Pb Horowitz, stretched by 1%

  char  str1[50], str2[50], str3[50];
  float f1, f2, f3, f4, f5, f6;

  FILE * fp;

  fp = fopen (filename.data(), "r");
  //cout << imax << " " << jmax << " " << kmax << endl;
  for( int i = 0; i < imax; i++ ){

    fscanf(fp, "%s", str1);

    for( int j = 0; j < jmax; j++ ){
      //cout << "j: " << j << " " << model << endl;
      if( model == 0 ){
	fscanf(fp, "%f %f %f %f %f %f", &f1, &f2, &f3, &f4, &f5, &f6);
	//format: &angle, &cross_section, &ignore1, &asymmetry, &ignore2, &ignore3
	/*
	  printf("|%f|"  , f1 );
	  printf("|%f|"  , f2 );
	  printf("|%f|"  , f3 );
	  printf("|%f|"  , f4 );
	  printf("|%f|"  , f5 );
	  printf("|%f|\n", f6 );
      */
	table[index][i][j][0] = f1;
	table[index][i][j][1] = f2;
	table[index][i][j][2] = f3;//ignore
	table[index][i][j][3] = f4;
	table[index][i][j][4] = f5;//ignore
	table[index][i][j][5] = f6;//ignore
      }else if ( model == 1 ){
	fscanf(fp, "%f %f %f", &f1, &f2, &f3);
	//format: &angle, &cross_section, &ignore1, &asymmetry, &ignore2, &ignore3
	
	//printf("|%f|"  , f1 );
	//printf("|%f|"  , f2 );
	//printf("|%f|\n"  , f3 );
	
	table[index][i][j][0] = f1;
	table[index][i][j][1] = f2;
	table[index][i][j][2] = 0.;//ignore
	table[index][i][j][3] = f3;
	table[index][i][j][4] = 0.;//ignore
	table[index][i][j][5] = 0.;//ignore
      }
    }
  }
  fclose(fp);

  cout << "Table loaded" << endl;
}

void HRSDatabase::PrintTable(int index){
  //indices:
  //1. Pb Horowitz
  //2. Pb Horowitz, stretched by 1%
  for( int i = 0; i < imax; i++ ){
    for( int j = 0; j < jmax; j++ ){
      for( int k = 0; k < kmax; k++ ){
	cout << table[index][i][j][k] << " ";
      }
      cout << endl;
    }
  }

}

float HRSDatabase::Interpolate(float energy, float theta, int index, int which){
  //cout << "Preparing to interpolate." << endl;
  //cout << energy << " " << ene_min << " " << ene_step << endl;
  //index is which model
  //indices:
  //0. Pb Horowitz
  //1. Pb Horowitz, stretched by 1%
  //which is whether you want:
  //0. Cross section
  //1. Asymmetry
  float last_theta = 0.0; //Last: theta,
  float curr_theta = 0.0; //Current: theta,
  float last_energy = 0.0; //Last: theta,
  float curr_energy = 0.0; //Current: theta,

  int   last_i     = 0;
  int   curr_i     = 0;
  int   last_j     = 0;
  int   curr_j     = 0;

  //int i = float( ( energy - ene_min * 1000. ) / ( ene_step * 1000. ) + 0.5);
  
  //cout << energy << " " << i << " " << ( energy - ene_min * 1000. ) / ( ene_step * 1000. ) << endl;

  //for( int i = 0; i < imax; i++ ){
  for( int j = 0; j < jmax; j++ ){
    if( (float)curr_theta < (float)theta){
      last_theta = curr_theta;
      curr_theta = table[index][0][j][0];
      last_j     = curr_j;
      curr_j     = j;
      //cout << last_theta << " " << theta << " " << curr_theta << 
      //", and j = " << j << "/" << jmax << endl;
    }else{
      //cout << "FOUND!" << " i = " << i << "/" << imax << ", and j = " << j << "/" << jmax << endl;
      j = jmax; //exit loops
    }
  }

  for( int i = 0; i < imax; i++ ){
    if( (float)curr_energy < (float)energy){
      last_energy = curr_energy;
      curr_energy = ene_min + ene_step * i;
      curr_energy *= 1000.;
      last_i     = curr_i;
      curr_i     = i;
      //cout << (float)last_energy << " " << energy << " " << (float)curr_energy << 
      //", and i = " << i << "/" << imax << endl;
    }else{
      //cout << "FOUND!" << " i = " << i << "/" << imax << ", and j = " << j << "/" << jmax << endl;
      i = imax; //exit loops
    }
  }
  
  float xs_t1_e1    = table[index][last_i][last_j][1];
  float a_t1_e1     = table[index][last_i][last_j][3];

  float xs_t2_e1    = table[index][last_i][curr_j][1];
  float a_t2_e1     = table[index][last_i][curr_j][3];

  float xs_t1_e2    = table[index][curr_i][last_j][1];
  float a_t1_e2     = table[index][curr_i][last_j][3];

  float xs_t2_e2    = table[index][curr_i][curr_j][1];
  float a_t2_e2     = table[index][curr_i][curr_j][3];
  

  float answer;
  if( which == 0 ){
    float m_e1 = ( xs_t2_e1 - xs_t1_e1 ) / ( curr_theta - last_theta );
    float m_e2 = ( xs_t2_e2 - xs_t1_e2 ) / ( curr_theta - last_theta );
    float xs1  = xs_t1_e1 + m_e1 * ( theta - last_theta );
    float xs2  = xs_t1_e2 + m_e2 * ( theta - last_theta );

    float m_t  = ( xs2 - xs1 ) / ( curr_energy - last_energy );
    answer = xs1 + m_t * ( energy - last_energy );

  }else if( which == 1 ){
    float m_e1 = ( a_t2_e1 - a_t1_e1 ) / ( curr_theta - last_theta );
    float m_e2 = ( a_t2_e2 - a_t1_e2 ) / ( curr_theta - last_theta );
    float a1  = a_t1_e1 + m_e1 * ( theta - last_theta );
    float a2  = a_t1_e2 + m_e2 * ( theta - last_theta );

    float m_t  = ( a2 - a1 ) / ( curr_energy - last_energy );
    answer = a1 + m_t * ( energy - last_energy );

  }else{
    cout << "Not an interpolation option: " << which << endl;
  }

  //cout << answer << " " << curr_theta << " " << last_theta << " " << curr_energy << " " << last_energy << endl;

  /*
  if( which == 0 ){
    float m = ( curr_xs - last_xs ) / ( curr_theta - last_theta );
    answer = last_xs + m * ( theta - last_theta );
  }else if( which == 1 ){
    float m = ( curr_a - last_a ) / ( curr_theta - last_theta );
    answer = last_a  + m * ( theta - last_theta );
  }else{
    cout << "Not an interpolation option: " << which << endl;
  }
  */
  //cout << last_theta << " " << theta << " " << curr_theta << endl;
  return answer;

}

float HRSDatabase::TotalCrossSection(float energy, float theta1, float theta2, float phi1, float phi2){//must be input in degrees
  float th1 = theta1 * 3.141592654 / 180.;
  float th2 = theta2 * 3.141592654 / 180.;
  float dth;
  float ph1 = phi1   * 3.141592654 / 180.;
  float ph2 = phi2   * 3.141592654 / 180.;
  float dph = ph2 - ph1;

  float curr_theta = 0.0; //Current: theta,
  float curr_xs    = 0.0; //Cross section, and

  //int i = ( ( energy - 550 ) / 50);
  int i = float( ( energy - ene_min * 1000.) / ( ene_step * 1000. ) + 0.5);
  //cout << energy << " " << i << endl;

  float answer;

  //for( int i = 0; i < imax; i++ ){
  for( int j = 0; j < jmax; j++ ){
    curr_theta = table[0][i][j][0] * 3.141592654 / 180.;
    curr_xs    = table[0][i][j][1];
    float theta_step = 0.2 * 3.141592654 / 180.;
    //dth        = -cos( curr_theta + theta_step / 2. ) + cos( curr_theta - theta_step / 2. );
    dth        = sin( curr_theta + theta_step / 2. ) - sin( curr_theta - theta_step / 2. );
    if( th1 < curr_theta && curr_theta < th2 )
      answer += curr_xs * dph * dth;
    //cout << answer << endl;
  }
  
  //cout << last_theta << " " << theta << " " << curr_theta << endl;
  return answer;

}
