// ********************************************************************
//
// $Id: SBSDetectorConstruction.hh,v 1.00, 2013/10/06 HRS Exp $
// --------------------------------------------------------------
//
//Construct RTPC detector
//This is not an independent class, it has to be invoke by class HRSDetectorConstruction

#ifndef SBSDetectorConstruction_H_
#define SBSDetectorConstruction_H_ 1

#include "globals.hh"	//for units and g4types and g4io
#include "HRSVisAttribute.hh"
#include "G4VUserDetectorConstruction.hh"
#include "HRSMaterial.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;

class SBSDetectorConstruction : public G4VUserDetectorConstruction, public HRSVisAttribute
{
public:
	SBSDetectorConstruction(G4LogicalVolume *mother=0);
	~SBSDetectorConstruction();

	G4VPhysicalVolume* Construct();
	G4VPhysicalVolume* ConstructSBSCCPrototype(G4LogicalVolume* motherLogical);


private:
	void GetConfig();
	void ConstructMaterials();

private:
	G4LogicalVolume* mMotherLogVol;
	HRSMaterialManager* mMaterialManager;

private:
	//config block
	double mPivotXOffset,mPivotYOffset,mPivotZOffset;
	double mPivot2SuperBigBiteFace;
	double mSuperBigBiteAngle;
	//config block end

};

#endif  //SBSDetectorConstruction_H_

