// ********************************************************************
//
// $Id: HMSDetectorConstruction.hh,v 1.00, 2013/10/06 HRS Exp $
// --------------------------------------------------------------
//
//Construct HMS detector
//This is not an independent class, it has to be invoke by class HRSDetectorConstruction

#ifndef HMSDetectorConstruction_H_
#define HMSDetectorConstruction_H_ 1

#include "globals.hh"	//for units and g4types and g4io
#include "HRSVisAttribute.hh"
#include "G4VUserDetectorConstruction.hh"
#include "HRSMaterial.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;

class HMSDetectorConstruction : public G4VUserDetectorConstruction, public HRSVisAttribute
{
public:
	HMSDetectorConstruction(G4LogicalVolume *mother=0);
	~HMSDetectorConstruction();

	G4VPhysicalVolume* Construct();
	G4VPhysicalVolume* ConstructHMS(G4LogicalVolume* motherLogical);


private:
	void GetConfig();
	void ConstructMaterials();

private:
	G4LogicalVolume* mMotherLogVol;
	HRSMaterialManager* mMaterialManager;

private:
	//config block
	int    mSetupHMS;
	double mPivotXOffset,mPivotYOffset,mPivotZOffset;
	double mPivot2HMSFace;
	double mHMSAngle;
	//config block end

};

#endif  //HMSDetectorConstruction_H_

