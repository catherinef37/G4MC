#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
using namespace std;

extern"C" {
  float x_sp_fp_(float *input, int *dim);
  float t_sp_fp_(float *input, int *dim);
  float y_sp_fp_(float *input, int *dim);
  float p_sp_fp_(float *input, int *dim);
  float l_sp_fp_(float *input, int *dim);
  float x_sp_cden_(float *input, int *dim);
  float t_sp_cden_(float *input, int *dim);
  float y_sp_cden_(float *input, int *dim);
  float p_sp_cden_(float *input, int *dim);
  float l_sp_cden_(float *input, int *dim);
}

int main(){
  
  //ofstream OUTFILE1;
  //ofstream OUTFILE2;

  int dim = 5;
  
  //0 -0.00129336 -0.000263345 0.00377103 -0.00492643 
  //-7.40116090E-02 This is a sample to show we are in tune X
  //-0.0740116 -0.0126819 -0.00901139 -0.00489415
  
  //float input[] = {0, -0.00129336, -0.000263345, 0.00377103, -0.00492643 };
  //float input[] = {0.01, 0.01, 0.01, 0.01, 0.01 };
  float input [5] = {0.0, 0.0, 0.0, 0.0, 0.0 };
  float output[5] = {0.0, 0.0, 0.0, 0.0, 0.0 };
  output[0] = x_sp_cden_(input, &dim);
  output[1] = t_sp_cden_(input, &dim);
  output[2] = y_sp_cden_(input, &dim);
  output[3] = p_sp_cden_(input, &dim);
  output[4] = l_sp_cden_(input, &dim);
  for(int i = 0; i < 5; i++){
    //cout << input [i] << " ";
    cout << output[i] << " ";
  }

  /*
  //OUTFILE1.open("phiy.dat");
  //for(int i = 0; i < 200; i++ ){
  float tanphi = ( -10. + 10. / 100. * i ) * 3.141592654 / 180.;
  float input[] = {0., 0., 0., tanphi, 0.};
  output[0] = x_sp_cfp_(input, &dim);
  output[1] = t_sp_cfp_(input, &dim);
  output[2] = y_sp_cfp_(input, &dim);
  output[3] = p_sp_cfp_(input, &dim);
  output[4] = l_sp_cfp_(input, &dim);
  //for(int i = 0; i < 5; i++){
  //cout << input[i] << " ";
  //cout << output[i] << " ";
  //}
  OUTFILE1 << atan(tanphi) << " " << output[2];
  OUTFILE1 << endl;
}
  OUTFILE1.close();
  //cout << endl;
  float  pi         = 3.14159265358979;  
  //float z           = 0.05; //m                                     #how far up the target are we looking?
  float  planeangle = 85.* pi / 180.;
  int    index      = 0;
  for(float z = -0.1; z < 0.15; z+=0.05){
    string name       = "phiyprime";
    string which[]      = {"1", "2", "3", "4", "5"};
    name.append(which[index]);
    name.append(".dat");
    cout << name;
    OUTFILE2.open(name.data());
    index++;
    for (int i = 0; i < 200; i++){
      float ztanphi   = ( -10. + 10. / 100. * i ) * pi / 180.;//               #phi at z
      float zphi      = atan(ztanphi);//                                       #tan of phi, in other words, transport phi
      float drift     = z * sin(planeangle) / sin( pi -planeangle - zphi );
      float y         = z * cos(planeangle)    ;
      y              -= z * ztanphi;
      float input2[] = {0., 0., y, ztanphi, 0.};
      float vby = y_sp_cfp_(input2, &dim);
      OUTFILE2 << zphi << " " << vby << endl;
    }
    OUTFILE2.close();
  }
*/
  return 0;
}
