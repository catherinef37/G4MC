// ********************************************************************
//
// $Id: SBSDetectorConstruction.cc,v 1.00, 2013/10/06 HRS Exp $
// --------------------------------------------------------------
//
// ********************************************************************
//
#include "SBSDetectorConstruction.hh"

#include "G4FieldManager.hh"
#include "G4MagneticField.hh"
#include "G4ChordFinder.hh"
#include "G4TransportationManager.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4UniformMagField.hh"

#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trap.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

#include "G4PVReplica.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4Polycone.hh"
#include "G4AssemblyVolume.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4UserLimits.hh"
#include "G4ios.hh"

#include "HRSStdSD.hh"
#include "HRSDCSD.hh"
#include "UsageManager.hh"
#include "HRSMaterial.hh"

#include "BField_SBS.hh"

extern UsageManager* gConfig;

////////////////////////////////////////////////////////////////////////////////////////////////
SBSDetectorConstruction::SBSDetectorConstruction(G4LogicalVolume *mother): 
mMotherLogVol(mother) 
{
	GetConfig();

	mMaterialManager = HRSMaterial::GetHRSMaterialManager();
	//construct the material manager, this line should behind GetConfig()
	//since it need to access the buffer of gConfig
	ConstructMaterials();

	G4cout<<"Contruct SBS geometry ... done! "<<G4endl;
}

SBSDetectorConstruction::~SBSDetectorConstruction()
{
	//I might need to delete the materials
	//But it does not hurt if I don't, since this class will have only one instance
	//in the program
	G4cout<<"Delete RTPC geometry ... done! "<<G4endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////
void SBSDetectorConstruction::GetConfig()
{
	gConfig->ReadFile("Detector_SBS.ini");
	//////////////////////////////////////////////////////////////////
	//global variables
	gConfig->GetParameter("PivotXOffset",mPivotXOffset);
	mPivotXOffset*=mm;
	gConfig->GetParameter("PivotYOffset",mPivotYOffset);
	mPivotYOffset*=mm;
	gConfig->GetParameter("PivotZOffset",mPivotZOffset);
	mPivotZOffset*=mm;

	//////////////////////////////////////////////////////////////////
	//module variables
	gConfig->GetParameter("SuperBigBiteAngle",mSuperBigBiteAngle);
	mSuperBigBiteAngle*=deg;
	gConfig->GetParameter("Pivot2SuperBigBiteFace",mPivot2SuperBigBiteFace);
	mPivot2SuperBigBiteFace*=mm;
}


////////////////////////////////////////////////////////////////////////////////////////////////
void SBSDetectorConstruction::ConstructMaterials()
{
	//a lot of materials have been built in the HRSMaterial
	//you can get the HRSMaterialManager then use them
	//or you can build your own material
	//Note that you have to delete whatever you have built
	//during the deconstruction
	//mMaterialManager = HRSMaterial::GetHRSMaterialManager();
}


////////////////////////////////////////////////////////////////////////////////////////////////


G4VPhysicalVolume* SBSDetectorConstruction::Construct()
{
	const double inch=2.54*cm;
	// All G4VSensitiveDetector will be managed (deleted) by SDManager
	//therefore no need to delete it at the end of this subroutine
	G4String SDname;
	G4SDManager* SDman = G4SDManager::GetSDMpointer();

	G4VSensitiveDetector* SBSVBSD=new HRSStdSD(SDname="SBSVBSD");
	SDman->AddNewDetector(SBSVBSD); 

	G4RotationMatrix *pRotSBS=new G4RotationMatrix();
	pRotSBS->rotateY(-mSuperBigBiteAngle); 

	G4VPhysicalVolume* SBSContainerPhys = 0;
	/////////////////////////
	// SBS Container
	/////////////////////////
	//This container is 300 inches long, its front face is the 
	//upstream end of the magnet, 63 inches to the target center.

	double pSBSX=95*inch, pSBSY=180*inch, pSBSZ=300*inch;

	double pSBSPos_X=(mPivot2SuperBigBiteFace+pSBSZ/2.0)*sin(mSuperBigBiteAngle)+mPivotXOffset;
	double pSBSPos_Y=mPivotYOffset;
	double pSBSPos_Z=(mPivot2SuperBigBiteFace+pSBSZ/2.0)*cos(mSuperBigBiteAngle)+mPivotZOffset;

	G4VSolid* SBSContainerSolid = new G4Box("SBSContainerBox",pSBSX/2.0+1.0*cm,
		pSBSY/2.0+1.0*cm,pSBSZ/2.0+1.0*cm);
	G4LogicalVolume* SBSContainerLogical = new G4LogicalVolume(SBSContainerSolid,
		mMaterialManager->air,"SBSContainerLogical",0,0,0);
	SBSContainerLogical->SetVisAttributes(HallVisAtt);  

	SBSContainerPhys=new G4PVPlacement(pRotSBS,
		G4ThreeVector(pSBSPos_X,pSBSPos_Y,pSBSPos_Z),
		SBSContainerLogical,"SBSContainerPhys",mMotherLogVol,0,0,0);

	/////////////////////////////////////
	//The magnet block
	/////////////////////////////////////
	//the whole box is 92"(width) x 173"(height) x 50" (length)
	//the aperture is 47(w) x 121.9(h) x 122(length) cm
	double pSBSMagX=92*inch, pSBSMagY=173*inch, pSBSMagZ=122*cm;
	G4VSolid* SBSMagWholeSolid = new G4Box("SBSMagWholeBox",pSBSMagX/2.0,
		pSBSMagY/2.0,pSBSMagZ/2.0);
	//the aperture that will be cut out, need to make it 2 mm longer 
	double pSBSMagAperX=47*cm, pSBSMagAperY=121.9*cm, pSBSMagAperZ=pSBSMagZ+1*mm;
	G4VSolid* SBSMagAperSolid = new G4Box("SBSMagAperBox",pSBSMagAperX/2.0,
		pSBSMagAperY/2.0,pSBSMagAperZ/2.0);

	G4SubtractionSolid* SBSMagBlockSolid = new G4SubtractionSolid("SBSMagBlockSolid",
		SBSMagWholeSolid,SBSMagAperSolid);

	G4LogicalVolume* SBSMagBlockLogical = new G4LogicalVolume(SBSMagBlockSolid,
		mMaterialManager->siliconsteel,"SBSMagBlockLogical",0,0,0);
	SBSMagBlockLogical->SetVisAttributes(DarkBlueVisAtt);  

	double pSBSMagPos_X=0,pSBSMagPos_Y=0;
	double pSBSMagPos_Z=-pSBSZ/2+pSBSMagZ/2;
	new G4PVPlacement(0,G4ThreeVector(pSBSMagPos_X,pSBSMagPos_Y,pSBSMagPos_Z),
		SBSMagBlockLogical,"SBSMagBlockPhys",SBSContainerLogical,0,0,0);


	/////////////////////////////////////
	//The magnet field volumn 
	/////////////////////////////////////
	//Put an uniform magnetic field into the aperture, this is a local field 
	//Do not use the global field manager. What you need to do is 
	//pass this field manager to the logical volumn during construction

	//I use a field map now, do not need this local field any more
	//if you really want to use this local field, you should disable 
	//the SBS field first (from the ini file)
	bool UseUniformField = false;
	if(UseUniformField)
	{
		G4FieldManager* SBSFieldMan = new G4FieldManager();
		//G4FieldManager* fieldMgr
		//	= G4TransportationManager::GetTransportationManager()->GetFieldManager();
		double pFieldVal = -1.8*tesla;
		G4ThreeVector SBSFieldV3( pFieldVal*cos(mSuperBigBiteAngle), 0.0,
			-pFieldVal*sin(mSuperBigBiteAngle)); 
		G4UniformMagField* SBSMagField = new G4UniformMagField(SBSFieldV3);
		SBSFieldMan->SetDetectorField(SBSMagField);
		SBSFieldMan->CreateChordFinder(SBSMagField);

		G4VSolid* SBSMagSolid = new G4Box("SBSMagBox",pSBSMagAperX/2.0,
			pSBSMagAperY/2.0,pSBSMagZ/2.0);

		G4LogicalVolume* SBSMagLogical = new G4LogicalVolume(SBSMagSolid,
			mMaterialManager->heliumGas,"SBSMagLogical",0,0,0,0);
		SBSMagLogical->SetFieldManager(SBSFieldMan,true);
		SBSMagLogical->SetVisAttributes(HallVisAtt);  

		new G4PVPlacement(0,G4ThreeVector(pSBSMagPos_X,pSBSMagPos_Y,pSBSMagPos_Z),
			SBSMagLogical,"SBSMagPhys",SBSContainerLogical,0,0,0);
	}
	else
	{	
		//By Jixie: 
		//I have to check the position and anlge for SBS field in case the user
		//forgot to change them in the ini file
		double pSBS_OriginX = 0;
		gConfig->GetParameter("Tosca_SBS_OriginX",pSBS_OriginX);
		pSBS_OriginX *= cm;
		double pSBS_OriginY = 0;
		gConfig->GetParameter("Tosca_SBS_OriginY",pSBS_OriginY);
		pSBS_OriginY *= cm;
		double pSBS_OriginZ = 0;
		gConfig->GetParameter("Tosca_SBS_OriginZ",pSBS_OriginZ);
		pSBS_OriginZ *= cm;

		double pPivot2SBSFieldCenter = mPivot2SuperBigBiteFace+61*cm;
		double pSBSFieldX = pPivot2SBSFieldCenter * sin(mSuperBigBiteAngle);
		double pSBSFieldY = pSBS_OriginY;
		double pSBSFieldZ = pPivot2SBSFieldCenter * cos(mSuperBigBiteAngle);
		

		double pSBSRotAngle1 = 0.0;
		gConfig->GetParameter("Tosca_SBS_RotAngle1",pSBSRotAngle1); 
		pSBSRotAngle1 *= deg;
		double pSBSRotAngle2 = 0.0;
		gConfig->GetParameter("Tosca_SBS_RotAngle2",pSBSRotAngle2);
		pSBSRotAngle2 *= deg;
		double pSBSRotAngle3 = 0.0;
		gConfig->GetParameter("Tosca_SBS_RotAngle3",pSBSRotAngle3);
		pSBSRotAngle3 *= deg;

		if( fabs(pSBSFieldX-pSBS_OriginX)>1*mm || fabs(pSBSFieldZ-pSBS_OriginZ)>1*mm || 
			fabs(pSBSRotAngle1 - (360*deg-mSuperBigBiteAngle))>1*deg )
		{
			cout<<"\n*** Warning: Your SBS field origin and angle does not match to SBS geometry,     ***\n"
				<<  "*** I have corrected them for you this time, please correct BField_SBS.ini ASAP. ***\n";
			gConfig->SetParameter("Tosca_SBS_OriginX",pSBSFieldX/cm);
			gConfig->SetParameter("Tosca_SBS_OriginZ",pSBSFieldZ/cm);
			gConfig->SetParameter("Tosca_SBS_RotAngle1",360.0-mSuperBigBiteAngle/deg);
			//now I have to update these 3 variable in pBField_SBS
			BField_SBS* pBField_SBS = BField_SBS::GetInstance();
			pBField_SBS->SetOrigin(pSBSFieldX/cm,pSBSFieldY/cm,pSBSFieldZ/cm); 
			pBField_SBS->SetRotAngle(360*deg-mSuperBigBiteAngle,pSBSRotAngle2,pSBSRotAngle3); 
		}
	}

	/////////////////////////////////////
	//virtual detector
	/////////////////////////////////////
	int mSetupSBSVD=1;
	if(mSetupSBSVD)
	{
		double pSBSVDX=pSBSMagAperX;
		double pSBSVDY=pSBSMagAperY;
		double SBSVDThick=1.0*mm;
		G4VSolid* SBSVDSolid = new G4Box("SBSVDBox",
			pSBSVDX/2.0, pSBSVDY/2.0, SBSVDThick/2.0);
		G4LogicalVolume* SBSVDLogical = new G4LogicalVolume(SBSVDSolid, 
			mMaterialManager->air,"SBSVDLogical",0,0,0);

		//5 cm behind the magnect 
		double pSBSVDPos_Z=-pSBSZ/2+pSBSMagZ+5*cm;

		new G4PVPlacement(0,G4ThreeVector(0,0,pSBSVDPos_Z),
			SBSVDLogical,"virtualBoundaryPhys_SBS",SBSContainerLogical,0,0,0);
		SBSVDLogical->SetSensitiveDetector(SBSVBSD);
		SBSVDLogical->SetVisAttributes(LightYellowVisAtt); 
	}

	return SBSContainerPhys;
}


//this routine is built for test SBS CC prototye
//do not need any longer
G4VPhysicalVolume* SBSDetectorConstruction::ConstructSBSCCPrototype(G4LogicalVolume* motherLogical)
{
	const double inch=2.54*cm;
	// All G4VSensitiveDetector will be managed (deleted) by SDManager
	//therefore no need to delete it at the end of this subroutine
	G4String SDname;
	G4SDManager* SDman = G4SDManager::GetSDMpointer();

	G4VSensitiveDetector* SBSCCSD=new HRSStdSD(SDname="SBSCCSD");
	SDman->AddNewDetector(SBSCCSD); 

	G4RotationMatrix *pRotSBS=new G4RotationMatrix();
	pRotSBS->rotateY(-mSuperBigBiteAngle); 

	/////////////////////////
	// SBS CC prototype Container
	/////////////////////////

	double pSBSX=20*inch, pSBSY=20*inch, pSBSZ=40*inch;

	double pSBSPos_X=(mPivot2SuperBigBiteFace+pSBSZ/2.0)*sin(mSuperBigBiteAngle)+mPivotXOffset;
	double pSBSPos_Y=mPivotYOffset;
	double pSBSPos_Z=(mPivot2SuperBigBiteFace+pSBSZ/2.0)*cos(mSuperBigBiteAngle)+mPivotZOffset;

	G4VSolid* SBSContainerSolid = new G4Box("SBSContainerBox",pSBSX/2.0+1.0*cm,
		pSBSY/2.0+1.0*cm,pSBSZ/2.0+5.0*cm);
	G4LogicalVolume* SBSContainerLogical = new G4LogicalVolume(SBSContainerSolid,
		mMaterialManager->air,"SBSContainerLogical",0,0,0);
	SBSContainerLogical->SetVisAttributes(DarkBlueVisAtt);  

	G4VPhysicalVolume* SBSContainerPhys=new G4PVPlacement(pRotSBS,
		G4ThreeVector(pSBSPos_X,pSBSPos_Y,pSBSPos_Z),
		SBSContainerLogical,"SBSContainerPhys",motherLogical,0,0,0);


	/////////////////////////////////////
	//virtual detector
	/////////////////////////////////////
	int mSetupSBSVD=1;
	if(mSetupSBSVD)
	{
		double pSBSVDX=pSBSX;
		double pSBSVDY=pSBSY;
		double SBSVDThick=1.0*mm;
		G4VSolid* SBSVDSolid = new G4Box("SBSVDBox",
			pSBSVDX/2.0, pSBSVDY/2.0, SBSVDThick/2.0);
		G4LogicalVolume* SBSVDLogical = new G4LogicalVolume(SBSVDSolid, 
			mMaterialManager->air,"SBSVDLogical",0,0,0);

		G4ThreeVector pVDCenter;
		if(mSetupSBSVD==1)
		{
			//0.5 cm in front of the box
			pVDCenter.set(0,0,-pSBSZ/2.-SBSVDThick/2.-0.5*cm);
		}
		else if(mSetupSBSVD==2)
		{
			//0.5 cm behind the box
			pVDCenter.set(0,0,pSBSZ/2.+SBSVDThick/2.+0.5*cm);
		}
		else
		{
			//in the center
			pVDCenter.set(0,0,0);
		}

		new G4PVPlacement(0,pVDCenter,
			SBSVDLogical,"virtualBoundaryPhys_SBS",SBSContainerLogical,0,0,0);

		SBSVDLogical->SetSensitiveDetector(SBSCCSD); 
		SBSVDLogical->SetVisAttributes(LightYellowVisAtt); 
	}

	return SBSContainerPhys;
}
