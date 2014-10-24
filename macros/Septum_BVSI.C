#include <iostream>
#include <stdio.h>
#include <math.h>
using namespace std;

//from S2
double Septum_Right_BVSI(double I)
{
	//B_gauss=(-3.32e-09)*x4+(2.99e-06)x3+(-6.46e-04)x2+(14.31)x+(6.38);
	double B_gauss =  -3.32e-09*pow(I,4.0) + 2.99e-06*pow(I,3.0) -
	6.46e-04*pow(I,2.0) + 14.31*I + 6.38;
	return -B_gauss; //the measurement is in the oppsite direction
}


double Septum_Left_BVSI(double I)
{
	//B_gauss=(-2.87e-09)*x4+(2.59e-06)x3+(-5.64e-04)x2+(14.52)x+(-3.14);
	double B_gauss =  -2.87e-09*pow(I,4.0) + 2.59e-06*pow(I,3.0) -
	5.64e-04*pow(I,2.0) + 14.52*I - 3.14;
	return B_gauss;
}

//Reference point:  
//Left(S3):  400A -->   5801.0 gauss --> 1.2876 GeV/c  
//Right(S2): 400A -->  -5719.0 gauss --> 1.2869 GeV/c  
//Left(S3):  950A -->  13167.0 gauss --> 2.8088 GeV/c  
//Right(S2): 950A --> -12873.0 gauss --> 2.8094 GeV/c  


//use 950I reference point to predict the current
// B0/P0 = B1/P1 ==> P1 = P0 * B1/B0;  or B1 = B0 * P1/P0 ;
void GetSeptumCurrent(double P/*, double &Ileft=0, double &Iright=0*/)
{
	double Ileft=0, Iright=0;
//left
    double Iref=950.0, Pref=2.809; //950A ==> 2.809 GeV/c
	
	double Bleft = Septum_Left_BVSI(Iref) * P/Pref;
	
	double tmpI = Iref*P/Pref;
	double tmpB = Septum_Left_BVSI(tmpI);
	for(int i=0;i<50;i++)
	{
		if(fabs((Bleft-tmpB)/Bleft)>1.0E-5)
		{
			tmpI += (Bleft-tmpB)/14.52;	
			tmpB = Septum_Left_BVSI(tmpI);	
			//cout<<"B="<<tmpB<<"  Bexpected="<<Bleft<<"  I="<<tmpI<<"  deltaB="<<Bleft-tmpB<<endl; 
		}
		else break;
	}
	Ileft = tmpI;
	//cout<<"P="<<P<<"  Ileft="<<Ileft<<"  Bleft="<<tmpB<<endl<<endl; 
	
	
	double Bright = Septum_Right_BVSI(Iref) * P/Pref;
	tmpB = Septum_Right_BVSI(tmpI);
	for(int i=0;i<50;i++)
	{
		if(fabs((Bright-tmpB)/Bleft)>1.0E-5)
		{
			tmpI -= (Bright-tmpB)/14.31;	
			tmpB = Septum_Right_BVSI(tmpI);	
			//cout<<"B="<<tmpB<<"  Bexpected="<<Bright<<"  I="<<tmpI<<"  deltaB="<<Bright-tmpB<<endl; 
		}
		else break;
	}
	Iright = tmpI;
	cout<<"P="<<P<<"  Ileft="<<Ileft<<"  Iright="<<Iright <<"  B="<<tmpB<<endl<<endl; 
}


int main(int argc, char** argv)
{
  if(argc<2)
    {
      cout<<"\t Usage: GetSeptaCurrent <HRSMomentum_gev>";
      return -1;
    }
  double P=2.254;
  P=atof(argv[1]);

  GetSeptumCurrent(P);
  return 0;
}
