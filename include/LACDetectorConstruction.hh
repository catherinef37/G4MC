// ********************************************************************
//
// $Id: LACDetectorConstruction.hh,v 1.00, 2013/10/06 HRS Exp $
// --------------------------------------------------------------
//
//This is the template to construct a detector
//This is not an independent class, it has to be invoke by class HRSDetectorConstruction

#ifndef LACDetectorConstruction_H_
#define LACDetectorConstruction_H_ 1

#include "globals.hh"	//for units and g4types and g4io
#include "HRSVisAttribute.hh"
#include "G4VUserDetectorConstruction.hh"
#include "HRSMaterial.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;

class LACDetectorConstruction : public G4VUserDetectorConstruction, public HRSVisAttribute
{
public:
	LACDetectorConstruction(G4LogicalVolume *mother=0);
	~LACDetectorConstruction();

	G4VPhysicalVolume* Construct();


private:
	void GetConfig();
	void ConstructMaterials();

private:
	G4LogicalVolume* mMotherLogVol;
	HRSMaterialManager* mMaterialManager;

private:
	//config block
	int    mSetupLAC;
	double mPivotXOffset,mPivotYOffset,mPivotZOffset;
	double mPivot2LACFace;
	double mLACAngle,mLACTiltAngle;
	double mLACYOffset;
	//config block end

};

#endif  //LACDetectorConstruction_H_

