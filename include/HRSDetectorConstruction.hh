// ********************************************************************
// $Id: HRSDetectorConstruction.hh,v 1.0, 2010/12/26 HRS Exp $
//
//
// ********************************************************************

#ifndef HRSDetectorConstruction_h
#define HRSDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "HRSVisAttribute.hh"
#include "HRSMaterial.hh"

#include <iostream>
#include <map>
using namespace std;

class HRSEMFieldSetup;

class G4VPhysicalVolume;
class G4VSensitiveDetector;
class G4LogicalVolume;

class HRSDetectorConstruction : public G4VUserDetectorConstruction, public HRSVisAttribute
{
public:
	HRSDetectorConstruction();
	virtual ~HRSDetectorConstruction();

public:
	virtual G4VPhysicalVolume* Construct();

	G4VPhysicalVolume* ConstructHRS(G4LogicalVolume* motherLogical);
	G4VPhysicalVolume* ConstructRadiator(G4LogicalVolume* motherLogical);

public:

	G4VPhysicalVolume* GetHallPhysVol(){return mHallPhysVol;};

private:

	void DumpGeometricalTree(G4VPhysicalVolume* aPhysVol,G4int depth=0,ostream &fout=G4cout);
	void GenerateGeometryMacro(G4VPhysicalVolume* aPhysVol,G4int depth=0);
	void GetConfig();

private:
	G4VPhysicalVolume* mHallPhysVol;

	//keep the information if the log vol has been printed or not
	map<string, int> mIsLogVolPrinted;   

	HRSMaterial* mMaterialManager;
	HRSEMFieldSetup* mEMFieldSetup;

private:

	//the following will be read from configuration file Detector.ini
	double mHallX,mHallY,mHallZ;
	double mFieldX,mFieldY,mFieldZ;
	int    mBMaterialType;

	int    mSetupAirAxis;

	double mPivotXOffset,mPivotYOffset,mPivotZOffset;
	double mScatChamberXOffset,mScatChamberYOffset,mScatChamberZOffset;
	double mTargetXOffset,mTargetYOffset,mTargetZOffset;

	int	   mSetupVirtualDetector;
	double mVirtualDetectorWidth,mVirtualDetectorHeight,mVirtualDetectorThick;
	double mVDRotYAngle,mVDRotXAngle;
	double mPivot2VDFace;
	string mVDPhysVolName;

	int    mSetupLHRS,mSetupRHRS;
	int    mSetupG2PGeometry;
	int    mSetupCREXGeometry; 
	int    mSetupRTPCGeometry;
	int    mSetupBigBite;
	int    mSetupSuperBigBite;
	int    mSetupHMS;
	int    mSetupLAC;

	//for HRS
	int    mSetupLSieveSlit,mSetupRSieveSlit;
	double mLHRSAngle,mLSeptumAngle,mRHRSAngle,mRSeptumAngle;

	//radiator
	int    mSetupRadiator;
	int    mSetupRadiatorVD;  //1 just set it up as VD, 2 set it up as virtual boundary
	//the material, 1 iron, 2 for copper, 3 for tantalum, 4 for tungsten
	int    mRadiatorType; 
	double mRaditorThickness, mRadiator2Pivot; 
};

#endif

