//calculate eleastic scatered energy
#include <iostream>
#include "math.h"
#include <stdio.h>
using namespace std;

//inout:
//beam in gev, theta in rad, target mass in GeV
double GetElasE(double beam_gev, double theta_rad, double M_gev=0.938)
{
  double Ei=beam_gev, t=theta_rad, M=M_gev;
  double Ef = Ei / (1+2*Ei/M*pow(sin(t/2),2.0));
  cout<< "Ei="<<Ei<< " Theta="<<t<<" M="<<M<<" ==> Ef =" <<Ef<<endl;
  return Ef; 
}


int main(int argc, char** argv)
{
  if(argc<4)
    {
      cout<<"\t Usage: GetElasE <beam_gev>  <theta_deg>  <targetmass_gev> \n"
	  <<"\t        HRS Angle is 5.65 deg, targetmass = 0.9315 * atomic_mass_number\n" 
	  <<"\t        Some useful number: M_H=0.938, M_He4=3.726, M_C12=11.178\n";
      return -1;
    }
  double beam=2.257, theta=5.65/57.3, M=0.938;
  beam=atof(argv[1]);
  theta=atof(argv[2])/57.3;
  M=atof(argv[3]);
  GetElasE(beam,theta,M);
  return 0;
}
