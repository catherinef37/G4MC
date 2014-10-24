// ********************************************************************
//
// $Id: HRSFZMagField.hh,v 1.0, 2010/12/26   HRS Exp $
// GEANT4 tag $Name: geant4-09-04 $
//

#ifndef HRSFZMagField_H
#define HRSFZMagField_H

#include "G4ThreeVector.hh"
#include "G4MagneticField.hh"

//class HRSFZMagFieldMessenger;
class HRSFZMagField : public G4MagneticField
{
public:

	HRSFZMagField() ;                
	~HRSFZMagField() ;  

	inline void GetFieldValue(const G4double pos[4], G4double *Bfield ) const;
	//  Point[4] x,y,z,time
	//  Return as Bfield[0], [1], [2] the magnetic field x, y & z components

	inline void SetFZB1Pos3V(G4ThreeVector v) { mFZB1BPos3V = v; }
	inline G4ThreeVector GetFZB1Pos3V() const { return mFZB1BPos3V; }

	inline void SetFZB2Pos3V(G4ThreeVector v) { mFZB2BPos3V = v; }
	inline G4ThreeVector GetFZB2Pos3V() const { return mFZB2BPos3V; }

	inline void SetFZB1Field3V(G4ThreeVector v) { mFZB1Field3V = v; bIsFZB1Uniform=true;}
	inline G4ThreeVector GetFZB1Field3V() const { return mFZB1Field3V; }

	inline void SetFZB2Field3V(G4ThreeVector v) { mFZB2Field3V = v; bIsFZB2Uniform=true;}
	inline G4ThreeVector GetFZB2Field3V() const { return mFZB2Field3V; }
	
	inline void SetFZB1TiltedAngle(G4double v) { mFZB1TiltedAngle = v; }
	inline G4double GetFZB1TiltedAngle() { return mFZB1TiltedAngle; }

	inline void SetFZB2TiltedAngle(G4double v) { mFZB2TiltedAngle = v; }
	inline G4double GetFZB2TiltedAngle() { return mFZB2TiltedAngle; }

private:
	//HRSFZMagFieldMessenger* messenger;
	//If the field is Uniform I use this flag to speed up the program
	bool bIsFZB1Uniform;
	bool bIsFZB2Uniform;
	G4ThreeVector mFZB1BPos3V,mFZB1Field3V;
	G4ThreeVector mFZB2BPos3V,mFZB2Field3V;

	double mPivotXOffset,mPivotYOffset,mPivotZOffset;
	double mFZB1TiltedAngle,mFZB1PosX,mFZB1PosY,mFZB1PosZ;
	double mFZB1Bx,mFZB1By,mFZB1Bz;
	double mFZB2TiltedAngle,mFZB2PosX,mFZB2PosY,mFZB2PosZ;
	double mFZB2Bx,mFZB2By,mFZB2Bz;

	//the following is the MAX after rotation
	double mFZB1MaxX,mFZB1MaxY,mFZB1MaxZ,mFZB2MaxX,mFZB2MaxY,mFZB2MaxZ;

	//the half after rotation
	double mFZB1HalfY,mFZB2HalfY;
};

#endif
