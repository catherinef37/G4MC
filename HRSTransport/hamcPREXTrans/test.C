#include <iostream>
#include <math.h>
using namespace std;

extern"C" {
  float x_sp_fp_(float *input, int *dim);
  float t_sp_fp_(float *input, int *dim);
  float y_sp_fp_(float *input, int *dim);
  float p_sp_fp_(float *input, int *dim);
  float l_sp_fp_(float *input, int *dim);
}

int main(){
  
  int dim = 5;
  
  //0 -0.00129336 -0.000263345 0.00377103 -0.00492643 
  //-7.40116090E-02 This is a sample to show we are in tune X
  //-0.0740116 -0.0126819 -0.00901139 -0.00489415

  float input[] = {0, -0.00129336, -0.000263345, 0.00377103, -0.00492643 };
  float output[]= {x_sp_fp_(input, &dim), t_sp_fp_(input, &dim),
		   y_sp_fp_(input, &dim), p_sp_fp_(input, &dim),
		   l_sp_fp_(input, &dim)};
  for(int i = 0; i < 5; i++){
    cout << output[i] << " ";
  }
  cout << endl;
  return 0;
}
