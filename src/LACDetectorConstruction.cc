// ********************************************************************
//
// $Id: LACDetectorConstruction.cc,v 1.00, 2013/10/06 HRS Exp $
// --------------------------------------------------------------
//
// ********************************************************************
//
#include "LACDetectorConstruction.hh"

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
#include "HRSCalorimeterSD.hh"
#include "UsageManager.hh"
#include "HRSMaterial.hh"


extern UsageManager* gConfig;

////////////////////////////////////////////////////////////////////////////////////////////////
LACDetectorConstruction::LACDetectorConstruction(G4LogicalVolume *mother): 
mMotherLogVol(mother) 
{
	GetConfig();

	mMaterialManager = HRSMaterial::GetHRSMaterialManager();
	//construct the material manager, this line should behind GetConfig()
	//since it need to access the buffer of gConfig
	ConstructMaterials();
	
	G4cout<<"Contruct LAC geometry ... done! "<<G4endl;
}

LACDetectorConstruction::~LACDetectorConstruction()
{
	//I might need to delete the materials
	//But it does not hurt if I don't, since this class will have only one instance
	//in the program
	G4cout<<"Delete LAC geometry ... done! "<<G4endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////
void LACDetectorConstruction::GetConfig()
{
	gConfig->ReadFile("Detector_LAC.ini");
	//////////////////////////////////////////////////////////////////
	//global variables
	gConfig->GetParameter("PivotXOffset",mPivotXOffset);
	mPivotXOffset*=mm;
	gConfig->GetParameter("PivotYOffset",mPivotYOffset);
	mPivotYOffset*=mm;
	gConfig->GetParameter("PivotZOffset",mPivotZOffset);
	mPivotZOffset*=mm;

	//////////////////////////////////////////////////////////////////
	
	gConfig->GetParameter("SetupLAC",mSetupLAC);

	//module variables
	gConfig->GetParameter("LACAngle",mLACAngle);
	mLACAngle*=deg;
	gConfig->GetParameter("LACTiltAngle",mLACTiltAngle);
	mLACTiltAngle*=deg;
	gConfig->GetParameter("Pivot2LACFace",mPivot2LACFace);
	mPivot2LACFace*=mm;
	gConfig->GetParameter("LACYOffset",mLACYOffset);
	mLACYOffset*=mm;
}


////////////////////////////////////////////////////////////////////////////////////////////////
void LACDetectorConstruction::ConstructMaterials()
{
	//a lot of materials have been built in the HRSMaterial
	//you can get the HRSMaterialManager then use them
	//or you can build your own material
	//Note that you have to delete whatever you have built
	//during the deconstruction
	//mMaterialManager = HRSMaterial::GetHRSMaterialManager();
}


////////////////////////////////////////////////////////////////////////////////////////////////


G4VPhysicalVolume* LACDetectorConstruction::Construct()
{
	// All G4VSensitiveDetector will be managed (deleted) by SDManager
	//therefore no need to delete it at the end of this subroutine
	G4String SDname;
	G4SDManager* SDman = G4SDManager::GetSDMpointer();

	G4VSensitiveDetector* LACSD=new HRSCalorimeterSD(SDname="LACSD");
	SDman->AddNewDetector(LACSD); 
	G4VSensitiveDetector* LACVBSD=new HRSStdSD(SDname="LACVBSD");
	SDman->AddNewDetector(LACVBSD); 

	G4RotationMatrix *pRotLAC=new G4RotationMatrix();
	pRotLAC->rotateY(-mLACAngle); 
	if(fabs(mLACTiltAngle/deg)>1.0E-5) pRotLAC->rotateY(-mLACTiltAngle); 

	double pLACWidth=217.0*cm;
	double pLACHeight=400.0*cm;
	double pLACLeadThick=2.0*mm;
	double pLACSCThick=15.0*mm;

	G4VPhysicalVolume* LACContainerPhys = 0;
	/////////////////////////
	// LAC Container
	/////////////////////////
	//This container is 220(W) x 420(H) x  1.7x33+20 cm, 
	double pLACX=pLACWidth, pLACY=pLACHeight;
	double pLACZ=(pLACSCThick+pLACLeadThick)*33;

	double pLACPos_X=(mPivot2LACFace+pLACZ/2.0)*sin(mLACAngle)+mPivotXOffset;
	double pLACPos_Y=mPivotYOffset+mLACYOffset;
	double pLACPos_Z=(mPivot2LACFace+pLACZ/2.0)*cos(mLACAngle)+mPivotZOffset;

	G4VSolid* LACContainerSolid = new G4Box("LACContainerBox",pLACX/2.0+5.0*cm,
		pLACY/2.0+5.0*cm,pLACZ/2.0+5.0*cm);
	G4LogicalVolume* LACContainerLogical = new G4LogicalVolume(LACContainerSolid,
		mMaterialManager->air,"LACContainerLogical",0,0,0);
	LACContainerLogical->SetVisAttributes(HallVisAtt);  

	LACContainerPhys=new G4PVPlacement(pRotLAC,
		G4ThreeVector(pLACPos_X,pLACPos_Y,pLACPos_Z),
		LACContainerLogical,"LACContainerPhys",mMotherLogVol,0,0,0);


	/////////////////////////////////////
	//the lead and SC layers
	/////////////////////////////////////
	
	if(mSetupLAC)
	{	

		G4VSolid* LACLeadLayerSolid = new G4Box("LACLeadLayerBox",pLACWidth/2.0,
			pLACHeight/2.0,pLACLeadThick/2.0);
		G4LogicalVolume* LACLeadLayerLogical = new G4LogicalVolume(LACLeadLayerSolid,
			mMaterialManager->lead,"LACLeadLayerLogical",0,0,0);
		LACLeadLayerLogical->SetVisAttributes(LeadVisAtt);  

		G4VSolid* LACSCLayerSolid = new G4Box("LACSCLayerBox",pLACWidth/2.0,
			pLACHeight/2.0,pLACSCThick/2.0);
		G4LogicalVolume* LACSCLayerLogical = new G4LogicalVolume(LACSCLayerSolid,
			mMaterialManager->scintillator,"LACSCLayerLogical",0,0,0);

		LACSCLayerLogical->SetVisAttributes(LightYellowVisAtt); 
		LACSCLayerLogical->SetSensitiveDetector(LACSD);


		double Zface = -pLACZ/2;
		double tmpZ = 0;
		for(int ii=0;ii<33;ii++)
		{
			//place lead layer
			tmpZ = Zface + (ii+0.5)*pLACLeadThick + ii*pLACSCThick;
			new G4PVPlacement(0,G4ThreeVector(0,0,tmpZ),
				LACLeadLayerLogical,"LACLeadLayerPhys",LACContainerLogical,true,ii+1);

			//place SC layer			
			tmpZ = Zface + (ii+1.0)*pLACLeadThick + (ii+0.5)*pLACSCThick;
			new G4PVPlacement(0,G4ThreeVector(0,0,tmpZ),
				LACSCLayerLogical,"LACSCLayerPhys",LACContainerLogical,true,ii+1);
		}

	}

	/////////////////////////////////////
	//virtual detector
	/////////////////////////////////////
	int mSetupLACVD=1;
	if(mSetupLACVD)
	{
		double pLACVDX=pLACX;
		double pLACVDY=pLACY;
		double LACVDThick=2.0*mm;
		G4VSolid* LACVDSolid = new G4Box("LACVDBox",
			pLACVDX/2.0, pLACVDY/2.0, LACVDThick/2.0);
		G4LogicalVolume* LACVDLogical = new G4LogicalVolume(LACVDSolid, 
			mMaterialManager->air,"LACVDLogical",0,0,0);

			//2 cm to the back of the layer
		double pLACVDPos_Z=pLACZ/2+2.0*cm;

		new G4PVPlacement(0,G4ThreeVector(0,0,pLACVDPos_Z),
			LACVDLogical,"virtualBoundaryPhys_LAC",LACContainerLogical,0,0,0);
		LACVDLogical->SetSensitiveDetector(LACVBSD);
		//LACVDLogical->SetVisAttributes(LightYellowVisAtt); 
		LACVDLogical->SetVisAttributes(DarkBlueVisAtt);
	}

	return LACContainerPhys;
}

