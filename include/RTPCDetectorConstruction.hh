// ********************************************************************
//
// $Id: RTPCDetectorConstruction.hh,v 1.00, 2013/10/06 HRS Exp $
// --------------------------------------------------------------
//
//Construct RTPC detector
//This is not an independent class, it has to be invoke by class HRSDetectorConstruction

#ifndef RTPCDetectorConstruction_H_
#define RTPCDetectorConstruction_H_ 1

#include "globals.hh"	//for units and g4types and g4io
#include "HRSVisAttribute.hh"
#include "G4VUserDetectorConstruction.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;

class RTPCDetectorConstruction : public G4VUserDetectorConstruction, public HRSVisAttribute
{
public:
	RTPCDetectorConstruction(G4LogicalVolume *mother=0);
	~RTPCDetectorConstruction();

	// the function which builds everything
	G4VPhysicalVolume* Construct();
	G4VPhysicalVolume* ConstructSolenoid(G4LogicalVolume *pMotherLogVol);

private:
	void GetConfig();
	void ConstructMaterials();

private:
	G4LogicalVolume* mMotherLogVol;
	double mTargetXOffset,mTargetYOffset,mTargetZOffset;

private:
	G4Material* vaccum;
	G4Material* air;

	G4Material* heliumGas;
	G4Material* aluminum;
	G4Material* bonusGas;
	G4Material* mixHeGas;
	G4Material* DMEGas;
	G4Material* mylar;
	G4Material* kapton;
	G4Material* deuteriumGas;
	G4Material* copper;
	G4Material* stainlesssteel;

	G4Material* H2TgGas;
	G4Material* D2TgGas;
	G4Material* He3TgGas;
	G4Material* He4TgGas;

	G4Material* ultem;
	G4Material* SiO2;
	G4Material* epoxy;
	G4Material* G10FR4;
	G4Material* rohacel71;
	G4Material* pcbNchip;
	G4Material* cable;

	G4Material* targetMaterial;
	G4Material* targetWallMaterial;

	G4double mD2GasD,mHeGasD,mMixHeGasD,mMixDMEGasD,mMixDMEHeD;

	G4double mD2GasL,mD2GasR,mD2GasT,mD2GasP,mHeGasT,mHeGasP;
	G4double mMixGasT,mMixGasP,mRatioHe2DME;

	G4int mTargetType;  //1=H2, 2=D2,3=He3,4=He4
	G4int mTgWallMaterialType;  //1 is kapton, 2 is aluminum
	G4double mRTPCLength,mTgWallThick;
	G4double m1stMylarR,m1stAlThick,m1stMylarThick;
	G4double m2ndMylarR,m2ndAlThick,m2ndMylarThick;
	
	G4double mGEM1R,mGEM2R,mGEM3R,mPCBReadOutR;

	G4double mBStepLimit,mDCStepLimit;

	G4double mBedPlateThick,mBedPlateLowEdge,mBedPlateHighEdge;
	G4double mInnerGapSpThick,mG10FR4Thick;
	G4double mGEM1SpThick,mGEM2SpThick,mGEM3SpThick,mReadOutSpThick;

	G4int mSetupSolenoid;
	G4double mSolenoidPosX,mSolenoidPosY,mSolenoidPosZ;

	G4int mSetupEntranceNExitCap,mSetupEndPlateNCover, mSetupCableNChip;

	G4double mLGEMV,mRGEMV,mLCathodeV,mRCathodeV;
};

#endif  //RTPCDetectorConstruction_H_

