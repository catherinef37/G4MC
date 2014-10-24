
#include "BigBiteMagField.hh"
#include "math.h"

//------------------------------
BigBiteMagField::BigBiteMagField(double armangle, double fieldval, double field2tg,
								 double pivotxoffset, double pivotyoffset, double pivotzoffset)
{
	detangle=armangle;		
	//make sure it is between [0,2Pi]
	dist2tg=field2tg; //150*cm;
	tgx=pivotxoffset;
	tgy=pivotyoffset;
	tgz=pivotzoffset;
	this->SetMagField(fieldval); 
}

//------------------------------

BigBiteMagField::~BigBiteMagField()
{
	;
}

//------------------------------
void BigBiteMagField::GetFieldValue(const double pos[3],double *B) const
{
	//This field is a local field, it will be placed into the magnetical volumn only
	//therefore no need to check the boundary of this volumn
	B[0] = Bx; B[1] = By; B[2] = Bz;
	//just for debugging
	//G4cout<<"x="<<pos[0]/cm<<" y="<<pos[1]/cm<<" z="<<pos[2]/cm<<"  Rxz="<<Rxz<<"(cm)\n\t"
	//	<<" dphi="<<(phi-detangle)/deg<<" dphi_limit="<<atan2(14*cm,dist2tg)/deg<<"(deg) \n\t"
	//	<<"==>  Bx="<<B[0]/tesla<<" Bz="<<B[2]/tesla<<" (Tesla)"<<G4endl;
	
	return;
}

//------------------------------
void BigBiteMagField::SetMagField(double fieldval)
{
	Btot = fieldval;  //Bx=-0.92*tesla if detangle=0;

	//?  The field manager is placed inside a container
	//should I put the field according to the global (hall) coor.
	//or the local coor.?
	Bx = fieldval*cos(detangle);
	By = 0; 
	Bz = -fieldval*sin(detangle); 
}
