/*
Based on G4 example.
Used to create the class for detector construction
*/

#ifndef ndDC_VAR
#define ndDC_VAR 1

#include "BigBiteMagField.hh"
#include "BigBiteHit.hh"
#include "HRSVisAttribute.hh"

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

#include "HRSStdSD.hh"

class BigBiteDetectorConstruction : public G4VUserDetectorConstruction, public HRSVisAttribute
{
public:
	BigBiteDetectorConstruction(G4LogicalVolume *mother=0);
	~BigBiteDetectorConstruction();

	// the function which builds everything
	G4VPhysicalVolume* Construct();

private:
	void GetConfig();

private:

	G4LogicalVolume* logMotherVol;
	int    mSetupBigBite, mSetupBigBiteSieve;

	double mPivotXOffset,mPivotYOffset,mPivotZOffset;
	double mPivot2BigBiteFace;   // distance from target center to BigBite front face 
	
	double mBigBiteAngle,mBigBiteTiltAngle;
	
	int    mSetupHAND;
	double mPivot2HANDLeadWall;

};

#endif

