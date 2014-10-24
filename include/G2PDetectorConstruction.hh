// ********************************************************************
// $Id: G2PDetectorConstruction.hh,v 1.0, 2010/12/26 HRS Exp $
// --------------------------------------------------------------
//

#ifndef G2PDetectorConstruction_h
#define G2PDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "HRSVisAttribute.hh"
#include "HRSMaterial.hh"


class HRSEMFieldSetup;

class G4VPhysicalVolume;
class G4VSensitiveDetector;
class G4LogicalVolume;

class G2PDetectorConstruction : public G4VUserDetectorConstruction, public HRSVisAttribute
{
public:
	G2PDetectorConstruction(G4LogicalVolume *mother=0);
	virtual ~G2PDetectorConstruction();

public:
	virtual G4VPhysicalVolume* Construct();
	G4VPhysicalVolume* ConstructG2PScatChamber(G4LogicalVolume* motherLogical);
	G4VPhysicalVolume* ConstructG2PTarget(G4LogicalVolume* motherLogical);
	G4VPhysicalVolume* ConstructG2PPVCTarget(G4LogicalVolume* motherLogical);
	G4VPhysicalVolume* ConstructG2PSeptumNSieve(G4LogicalVolume* motherLogical);
	G4VPhysicalVolume* ConstructG2PBeamDump(G4LogicalVolume* motherLogical);
	G4VPhysicalVolume* ConstructG2PThirdArm(G4LogicalVolume* motherLogical);
	G4VPhysicalVolume* ConstructG2PChicane(G4LogicalVolume* motherLogical);
	G4VPhysicalVolume* ConstructG2PPlatform(G4LogicalVolume* motherLogical);
	
	G4VPhysicalVolume* ConstructHRS(G4LogicalVolume* motherLogical);
	G4VPhysicalVolume* ConstructG2PHRS(G4LogicalVolume* motherLogical);

private:
	void ConstructMaterial();
	void GetConfig();

private:
	
	HRSMaterial* mMaterialManager;
	G4LogicalVolume *mMotherLogVol;

	G4Material *G2P_NH3He;


private:

	//////////////////////////
	//the following can be found in Detector_G2P.ini
	double mScatChamberRin,mScatChamberRout,mScatChamberL;
	double mScatChamberExitWindowThick;

	double mShieldLN2Rin,mShieldLN2Rout,mShieldLN2WindowThick,mShieldLN2L;
	double mShieldLHeRin,mShieldLHeRout,mShieldLHeL;

	int    mTargetType,mSetupG2PTarget;
	double mTargetL;

	int    mSetupG2PScatChamber,mSetupCoil;

	int    mSetupBeamDump;
	double mBeamDumpWidth,mBeamDumpHeight,mBeamDumpThick;
	double mPivot2BeamDumpFace;
	
	int    mSetupThirdArm;
	double mThirdArmAngle,mSetupThirdArmVD,mThirdArmRotZAngle,mPivot2ThirdArmFace;

	int    mSetupChicane,mSetupChicaneVD;

	int    mSetupPlatform;

	
	//////////////////////////
	//the following can be found in Detector_G2P.ini or Detecotor.ini
	double mPivotXOffset,mPivotYOffset,mPivotZOffset;
	double mScatChamberXOffset,mScatChamberYOffset,mScatChamberZOffset;
	double mTargetXOffset,mTargetYOffset,mTargetZOffset;

	int    mSetupLHRS,mSetupRHRS,mSetupLSieveSlit,mSetupRSieveSlit;
	double mLHRSAngle,mLSeptumAngle,mRHRSAngle,mRSeptumAngle;

	double mPivot2LSieveFace,mPivot2RSieveFace;

	double mPivot2LHRSVBFace,mPivot2RHRSVBFace;
	double mHRSVBWidth,mHRSVBHeight,mHRSVBThick;
};

#endif

