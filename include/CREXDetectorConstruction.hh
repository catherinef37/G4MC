/*
By Jixie Zhang @20121112
This is the class to built the target platform for CREX experiment
In order to optimize the design, a lot of variables will be read from
an input file 'Detector_CREX.ini'
*/

#ifndef CREX_Detector
#define CREX_Detector 1

#include "HRSVisAttribute.hh"

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

#include "HRSMaterial.hh"
#include "HRSStdSD.hh"

class CREXDetectorConstruction : public G4VUserDetectorConstruction, public HRSVisAttribute
{
public:
	// constructors and destructor
	CREXDetectorConstruction(G4LogicalVolume *mother=0);
	~CREXDetectorConstruction();

	// the function which builds everything
	G4VPhysicalVolume* Construct();
	
private:
	void ConstructMaterial();
	G4VPhysicalVolume* ConstructTargetChamber(G4LogicalVolume *pMotherLogVol);
	G4VPhysicalVolume* ConstructPREXTarget(G4LogicalVolume *pMotherLogVol);
	G4VPhysicalVolume* ConstructCREXTarget(G4LogicalVolume *pMotherLogVol);
	G4VPhysicalVolume* ConstructSeptumNSieve(G4LogicalVolume *pMotherLogVol);
	G4VPhysicalVolume* ConstructDVCSSolenoid(G4LogicalVolume *pMotherLogVol);

	void GetConfig();

private:

	HRSMaterial* mMaterialManager;
	//Declaration of all materials 
  G4Material *calcium,*air,*vacuum,*stainlesssteel,*aluminum,*heliumGas,*lead,*diamond,*lead208;
	G4Material *theTargetMaterial;

	G4LogicalVolume* mMotherLogVol;
	

private:

	//////////////////////////
	//the following can be found in Detector_CREX.ini
	int    mSetupStdScatChamber;
	double mScatChamberRin,mScatChamberRout,mScatChamberL;
	double mScatChamberEntranceWindowThick,mScatChamberExitWindowThick;

        int    mSetupPREXTarget,mSetupCREXTarget,mTargetType;
	double mTargetW,mTargetH,mTargetL;
	
	double mUpBlockRin,mUpBlockThick;
	double mDownBlockRin;
	double mDownBlockThickAt0cm,mDownBlockThickAt5cm;
	double mUpBlockLength,mDownBlockLength;

	//thickness of target window, not for scattering chamber
	double mUpCapThick,mDownCapThick;

	//vacuum chamber R,thickness and up|down cap thickness
	double mVCRin,mVCThick,mVCUpCapThick,mVCDownCapThick,mVCLength;

	//VacuumChamber to the center of scattering chamber at Z
	double mVC2SCZOffset;


	//////////////////////////
	//the following can be found in Detector_CREX.ini or Detecotor.ini
	double mPivotXOffset,mPivotYOffset,mPivotZOffset;
	double mScatChamberXOffset,mScatChamberYOffset,mScatChamberZOffset;
	double mTargetXOffset,mTargetYOffset,mTargetZOffset;

	int    mSetupLHRS,mSetupRHRS,mSetupLSieveSlit,mSetupRSieveSlit;
	double mLHRSAngle,mLSeptumAngle,mRHRSAngle,mRSeptumAngle;

	double mPivot2LSieveFace,mPivot2RSieveFace;

        double mPivot2LHRSVBFace,mPivot2RHRSVBFace,mPivot2CHRSVBFace;
	double mHRSVBWidth,mHRSVBHeight,mHRSVBThick;

};

#endif

