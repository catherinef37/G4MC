#ifndef ndMF_VAR
#define ndMF_VAR 1

#include "globals.hh"
#include "G4MagneticField.hh"


class BigBiteDetectorConstruction;

class BigBiteMagField : public G4MagneticField
{
public:
	//armangle is the phi angle relative to the beam line 
	BigBiteMagField(double armangle, double fieldval=-0.92*tesla, double field2tg=1.50*m,
		double pivotxoffset=0.0, double pivotyoffset=0.0, double pivotzoffset=0.0);
	~BigBiteMagField();

public:
	void GetFieldValue(const double point[3], double *field ) const;
	void SetTargetPosition(double x,double y,double z) {tgx=x;tgy=y;tgz=z;};
	void SetDetectorAngle(double v) {detangle=v;};
	void SetMagField(double);

private:
	double Bx,By,Bz,Btot;
	double dist2tg;				//the distance of field front face to target center, 
								//usually 0.4m+bgbite_container_front_face
	double tgx,tgy,tgz;			//target positon
	double detangle;			//bigbite arm angle
};

#endif

