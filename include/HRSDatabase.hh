#ifndef HRSDATABASE_HH
#define HRSDATABASE_HH

#include <string>
#include <iostream>
using namespace std;

class HRSDatabase{

public: // with description
  
  HRSDatabase(int);  
  ~HRSDatabase();  
  void  PrintTable(int);
  float Interpolate(float, float, int, int);
  float TotalCrossSection(float energy, float theta1, float theta2, float phi1, float phi2);

private:
  //not giving the user the ability to load tables herself - must do it through class
  void  LoadTable(string, int);
  float table[6][200][200][200];
  int model;
  int imax, jmax, kmax;
  float ene_min, ene_step;
  //float table;
  
};

#endif

