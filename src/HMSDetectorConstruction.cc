// ********************************************************************
//
// $Id: HMSDetectorConstruction.cc,v 1.00, 2013/10/06 HRS Exp $
// --------------------------------------------------------------
//
// ********************************************************************
//
#include "HMSDetectorConstruction.hh"

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
#include "G4Trd.hh"

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


extern UsageManager* gConfig;

////////////////////////////////////////////////////////////////////////////////////////////////
HMSDetectorConstruction::HMSDetectorConstruction(G4LogicalVolume *mother): 
mMotherLogVol(mother) 
{
	GetConfig();

	mMaterialManager = HRSMaterial::GetHRSMaterialManager();
	//construct the material manager, this line should behind GetConfig()
	//since it need to access the buffer of gConfig
	ConstructMaterials();

	G4cout<<"Contruct HMS geometry ... done! "<<G4endl;
}

HMSDetectorConstruction::~HMSDetectorConstruction()
{
	//I might need to delete the materials
	//But it does not hurt if I don't, since this class will have only one instance
	//in the program
	G4cout<<"Delete HMS geometry ... done! "<<G4endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////
void HMSDetectorConstruction::GetConfig()
{
	gConfig->ReadFile("Detector_HMS.ini");
	//////////////////////////////////////////////////////////////////
	//global variables
	gConfig->GetParameter("PivotXOffset",mPivotXOffset);
	mPivotXOffset*=mm;
	gConfig->GetParameter("PivotYOffset",mPivotYOffset);
	mPivotYOffset*=mm;
	gConfig->GetParameter("PivotZOffset",mPivotZOffset);
	mPivotZOffset*=mm;

	gConfig->GetParameter("SetupHMS",mSetupHMS);
	//////////////////////////////////////////////////////////////////
	//module variables
	gConfig->GetParameter("HMSAngle",mHMSAngle);
	mHMSAngle*=deg;
	gConfig->GetParameter("Pivot2HMSFace",mPivot2HMSFace);
	mPivot2HMSFace*=mm;
}


////////////////////////////////////////////////////////////////////////////////////////////////
void HMSDetectorConstruction::ConstructMaterials()
{
	//a lot of materials have been built in the HRSMaterial
	//you can get the HRSMaterialManager then use them
	//or you can build your own material
	//Note that you have to delete whatever you have built
	//during the deconstruction
	//mMaterialManager = HRSMaterial::GetHRSMaterialManager();
}


////////////////////////////////////////////////////////////////////////////////////////////////


G4VPhysicalVolume* HMSDetectorConstruction::Construct()
{
	const double inch=2.54*cm;
	// All G4VSensitiveDetector will be managed (deleted) by SDManager
	//therefore no need to delete it at the end of this subroutine
	G4String SDname;
	G4SDManager* SDman = G4SDManager::GetSDMpointer();

	G4VSensitiveDetector* HMSVBSD=new HRSStdSD(SDname="HMSVBSD");
	SDman->AddNewDetector(HMSVBSD); 

	G4RotationMatrix *pRotHMS=new G4RotationMatrix();
	pRotHMS->rotateY(-mHMSAngle); 

	G4VPhysicalVolume* HMSAperContainerPhys = 0;
	/////////////////////////
	// HMS Container
	/////////////////////////
	//This container is 8.5 cm long, its front face is the 
	//upstream end of the collimater,166.37 cm to the target center.

	double pHMSX=45*inch, pHMSY=70*inch, pHMSZ=8.5*cm;

	double pHMSPos_X=(mPivot2HMSFace+pHMSZ/2.0)*sin(mHMSAngle)+mPivotXOffset;
	double pHMSPos_Y=mPivotYOffset;
	double pHMSPos_Z=(mPivot2HMSFace+pHMSZ/2.0)*cos(mHMSAngle)+mPivotZOffset;

	G4VSolid* HMSAperContainerSolid = new G4Box("HMSAperContainerBox",pHMSX/2.0+1.0*cm,
		pHMSY/2.0+1.0*cm,pHMSZ/2.0+0.1*cm);
	G4LogicalVolume* HMSAperContainerLogical = new G4LogicalVolume(HMSAperContainerSolid,
		mMaterialManager->air,"HMSAperContainerLogical",0,0,0);
	HMSAperContainerLogical->SetVisAttributes(HallVisAtt);  

	HMSAperContainerPhys=new G4PVPlacement(pRotHMS,
		G4ThreeVector(pHMSPos_X,pHMSPos_Y,pHMSPos_Z),
		HMSAperContainerLogical,"HMSAperContainerPhys",mMotherLogVol,0,0,0);

	/////////////////////////////////////
	//The collimater block
	/////////////////////////////////////

	//////////////////
	//This is the shape of the aperture, BG=2*AH,  AD=2*BC
	//Entrance: AH=4.575*cm, BC=11.646*cm AD=23.292*cm, BG=9.150*cm  
	//Exit:     AH=4.759*cm, BC=12.114*cm AD=24.228*cm, BG=9.518*cm
	//  A ______ H
	//   /      \                                            1
	//  /        \                                           2
	//B|          | G
	// |          |
	// |          |
	// |          |
	// |          |
	//C\          / F
	//  \        /
	//  D ------ E
	///////////////////


	//the whole box
	double pHMSCollimaterX=30*inch, pHMSCollimaterY=60*inch, pHMSCollimaterZ=6.3*cm;
	G4VSolid* HMSCollimaterWholeSolid = new G4Box("HMSCollimaterWholeBox",
		pHMSCollimaterX/2.0,pHMSCollimaterY/2.0,pHMSCollimaterZ/2.0);
	//the aperture that will be cut out, need to make it 2 mm longer 
	double pColliAper_dh_entr=4.575*cm, pColliAper_dv_entr=11.646*cm;
	double pColliAper_dh_exit=4.759*cm, pColliAper_dv_exit=12.114*cm;
	double pColliAper_dz=pHMSCollimaterZ/2 + +1*mm; 
	G4VSolid* HMSCollimaterAperSolid = new G4Trd("HMSCollimaterAperTrd",
		pColliAper_dh_exit, pColliAper_dh_entr,
		pColliAper_dv_exit, pColliAper_dv_entr,
		pColliAper_dz);

	G4SubtractionSolid* HMSCollimaterSolid = new G4SubtractionSolid("HMSCollimaterSolid",
		HMSCollimaterWholeSolid,HMSCollimaterAperSolid);

	G4LogicalVolume* HMSCollimaterLogical = new G4LogicalVolume(HMSCollimaterSolid,
		mMaterialManager->absorber,"HMSCollimaterLogical",0,0,0);
	HMSCollimaterLogical->SetVisAttributes(IronVisAtt);  

	double pHMSCollimaterPos_X=0,pHMSCollimaterPos_Y=0;
	double pHMSCollimaterPos_Z=-pHMSZ/2+pHMSCollimaterZ/2;
	new G4PVPlacement(0,
		G4ThreeVector(pHMSCollimaterPos_X,pHMSCollimaterPos_Y,pHMSCollimaterPos_Z),
		HMSCollimaterLogical,"HMSCollimaterPhys",HMSAperContainerLogical,0,0,0);

	/////////////////////////////////////
	//virtual detector
	/////////////////////////////////////
	int mSetupHMSVD=1;
	if(mSetupHMSVD)
	{
		double pHMSVDX=pColliAper_dh_exit*2.0;
		double pHMSVDY=pColliAper_dv_exit*2.0;
		double HMSVDThick=1.0*mm;
		G4VSolid* HMSVDSolid = new G4Box("HMSVDBox",
			pHMSVDX/2.0, pHMSVDY/2.0, HMSVDThick/2.0);
		G4LogicalVolume* HMSVDLogical = new G4LogicalVolume(HMSVDSolid, 
			mMaterialManager->air,"HMSVDLogical",0,0,0);
		HMSVDLogical->SetSensitiveDetector(HMSVBSD);
		HMSVDLogical->SetVisAttributes(LightYellowVisAtt); 

		//2 mm behind the magnect 
		double pHMSVDPos_Z=-pHMSZ/2+pHMSCollimaterZ+2*mm;
		new G4PVPlacement(0,G4ThreeVector(0,0,pHMSVDPos_Z),
			HMSVDLogical,"virtualBoundaryPhys_HMS",HMSAperContainerLogical,0,0,0);
	}

	if(mSetupHMS>=2) ConstructHMS(mMotherLogVol);

	return HMSAperContainerPhys;
}


G4VPhysicalVolume* HMSDetectorConstruction::ConstructHMS(G4LogicalVolume* motherLogical)
{
	G4VPhysicalVolume* theHMSPhys=0;

	G4RotationMatrix *pRotX90deg=new G4RotationMatrix();
	pRotX90deg->rotateX(-90*deg); 

	//Build the HMS, .....
	//parameters need to be verified
	/////////////////////////
	// HMS QQDQ Containner
	/////////////////////////
	//Build a container using polycone, covering 270+/-12 degrees, 1.66m to 15.76m
	//3.05m below the beam line is the ground, need to subtract everything beneath that.
	//looks like this:
	/*
	HMS container:covering 24 degrees in X-Z plane, 10.05 m height.
	Stuff inside this containner will also use the pivot as the origin, but do not 
	need to worry about the rotation of the HMS.
	//                         7.0m above beam line 
	//                       -----------------------------| 15.76m from pivot,
	//                      /                             |
	//                     /                              |
	//                    /                               |
	//                   /                           E    |
	//                  /                         L       |
	//          --------                       O          |
	//         /                            P             |
	//    -----                          I                |
	//----|----- Q1 --- Q2 --- Q3 --- D   ----------------|------ beam line -----
	//    |                                               |
	//    ------------------------------------------------|
	//    1.66m from pivot, 3.05 m below beam line

	*/

	double pHMSContainerRin=mPivot2HMSFace+15*cm; //1.66*m,
	double pHMSContainerRout=15.76*m;
	double pBeamLine2Ground;
	pBeamLine2Ground=-3.05*m;
	//build the container with polycone

	const int kNPlane_HMSContainer=7;
	double rInner_HMSContainer[] = {pHMSContainerRin,pHMSContainerRin,2.5*m,
		3.7*m,9.0*m,pHMSContainerRout-3.0*m,pHMSContainerRout};
	double rOuter_HMSContainer[] = {pHMSContainerRout,pHMSContainerRout,pHMSContainerRout,
		pHMSContainerRout,pHMSContainerRout,pHMSContainerRout,pHMSContainerRout};
	double zPlane_HMSContainer[] = {-2.0*m,1.0*m,1.0*m,
		2.0*m,2.0*m,7.0*m,7.0*m};
	G4Polycone* HMSContainerSolid = new G4Polycone("HMSContainer",258.0*deg,24.0*deg,
		kNPlane_HMSContainer,zPlane_HMSContainer,rInner_HMSContainer,rOuter_HMSContainer);

	G4LogicalVolume* HMSContainerLogical = new G4LogicalVolume(HMSContainerSolid,
		mMaterialManager->vacuum,"HMSContainerLogical",0,0,0);

	HMSContainerLogical->SetVisAttributes(HallVisAtt); 
 
	G4RotationMatrix *pRotHMSContainer=new G4RotationMatrix();
	pRotHMSContainer->rotateX(-270*deg);
	pRotHMSContainer->rotateZ(-mHMSAngle);  

	if(mSetupHMS>=2)
	{
		new G4PVPlacement(pRotHMSContainer,G4ThreeVector(0,0,0),
			HMSContainerLogical,"HMSContainerPhys",motherLogical,0,0,0);
	}


	/////////////////////////
	// HMS Q1 
	/////////////////////////
	double pHallCenter2Q1Face=mPivot2HMSFace+25*cm;
	double pQ1Rin=25.0*cm;
	double pQ1Rout=65.0*cm;
	double pQ1Length=189*cm;

	G4VSolid* Q1Solid = new G4Tubs("Q1Tub",pQ1Rin,pQ1Rout,pQ1Length/2.0,0.0,360.0*deg);

	G4LogicalVolume* Q1Logical = new G4LogicalVolume(Q1Solid,
		mMaterialManager->siliconsteel,"Q1Logical",0,0,0);

	Q1Logical->SetVisAttributes(IronVisAtt); 

	if(mSetupHMS>=2)
	{
		double pQ1Pos_Z=(pHallCenter2Q1Face+pQ1Length/2.0);
		new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pQ1Pos_Z,0),
			Q1Logical,"Q1Phys",HMSContainerLogical,0,0,0);
	}


	/////////////////////////
	// HMS Q2 
	/////////////////////////
	double pHallCenter2Q2Face=4.0*m;
	double pQ2Rin=35.0*cm;
	double pQ2Rout=85.0*cm;
	double pQ2Length=210*cm;

	//By Jixie: HepRApp.jar has a bug to view this tubs. I am sure this tube is oriented correctly
	G4VSolid* Q2Solid = new G4Tubs("Q2Tub",pQ2Rin,pQ2Rout,pQ2Length/2.0,0.0,360.0*deg);

	G4LogicalVolume* Q2Logical = new G4LogicalVolume(Q2Solid,
		mMaterialManager->siliconsteel,"Q2Logical",0,0,0);

	Q2Logical->SetVisAttributes(IronVisAtt); 

	if(mSetupHMS>=3)
	{
		double pQ2Pos_Z=(pHallCenter2Q2Face+pQ2Length/2.0);
		new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pQ2Pos_Z,0),
			Q2Logical,"Q2Phys",HMSContainerLogical,0,0,0);
	}


	/////////////////////////
	// HMS Q3 
	/////////////////////////

	double pHallCenter2Q3Face=6.5*m;
	double pQ3Rin=35.0*cm;
	double pQ3Rout=85.0*cm;
	double pQ3Length=210*cm;

	G4VSolid* Q3Solid = new G4Tubs("Q3Tub",pQ3Rin,pQ3Rout,pQ3Length/2.0,0.0,360.0*deg);

	G4LogicalVolume* Q3Logical = new G4LogicalVolume(Q3Solid,
		mMaterialManager->siliconsteel,"Q3Logical",0,0,0);

	Q3Logical->SetVisAttributes(IronVisAtt); 

	if(mSetupHMS>=3)
	{
		double pQ3Pos_Z=(pHallCenter2Q3Face+pQ3Length/2.0);
		new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pQ3Pos_Z,0),
			Q3Logical,"Q3Phys",HMSContainerLogical,0,0,0);
	}

	/////////////////////////
	// HMS Dipole 
	/////////////////////////
	//The dipole is built as a disc subtraced subtract another disc to get the side 
	//then subtract by the tunnel disc
	double pDipoleBendAngle=25*deg, pDipoleFaceAngle=30*deg;
	double pDipoleR=8.4*m;
	double pDipoleRprime=pDipoleR*sin(pDipoleBendAngle/2)/
		sin(180*deg-pDipoleBendAngle/2-pDipoleFaceAngle);
	//double pDipoleRprime=4.0518*m;

	double pDipoleRCenterY=pDipoleR;
	double pDipoleRCenterZ=9.34*m;


	//the center of Rprime relative to R 
	double pRprime2R_Y=pDipoleRprime*cos(pDipoleFaceAngle)-pDipoleR;	// =-5.535*m;
	double pRprime2R_Z=pDipoleRprime*sin(pDipoleFaceAngle);				// =1.865*m;

	//the original disc
	G4VSolid* DipoleWholeTub = new G4Tubs("DipoleWholeTub",
		pDipoleR-0.8*m,pDipoleR+0.8*m,0.4*m,
		172*deg,pDipoleBendAngle+16*deg);
	//the disc to be subtracted, musu be thicker 
	G4VSolid* DipolePrimeTub = new G4Tubs("DipolePrimeTub",
		0,pDipoleR,0.5*m,
		180*deg+pDipoleFaceAngle+pDipoleBendAngle,360*deg-pDipoleFaceAngle*2-pDipoleBendAngle);
	//subtract the small tube to form the shape of the sides
	G4SubtractionSolid* DipoleWithSides = new G4SubtractionSolid("DipoleWithSides",
		DipoleWholeTub,DipolePrimeTub,
		0,G4ThreeVector(pRprime2R_Y,-pRprime2R_Z,0));

	//the tunnel disc, I use a rectangle shape here
	//G4VSolid* DipoleTunnelTub = new G4Tubs("DipoleTunnelTub",
	//	pDipoleR-0.4*m,pDipoleR+0.4*m,0.125*m,
	//	170*deg,pDipoleBendAngle+20*deg);
	//The shape of the dipole tunnel is a trapzoid, x=+/-0.4; y=+/-(0.125*(1-(1.25*x/8.40))
	//To build this shape, I have to use polycone
	double dy=0.125*(1.25*0.40/8.40)*m;
	const int kNPlane_DipoleTunnel=4;
	double rInner_DipoleTunnel[] = {pDipoleR-0.4*m,pDipoleR-0.4*m,pDipoleR-0.4*m,pDipoleR-0.4*m};
	double rOuter_DipoleTunnel[] = {pDipoleR-0.4*m,pDipoleR+0.4*m,pDipoleR+0.4*m,pDipoleR-0.4*m};
	double zPlane_DipoleTunnel[] = {-0.125*m-dy,-0.125*m+dy,0.125*m-dy,0.125*m+dy};
	G4Polycone* DipoleTunnelCone = new G4Polycone("DipoleTunnelCone",
		170.0*deg,pDipoleBendAngle+20*deg,
		kNPlane_DipoleTunnel,zPlane_DipoleTunnel,rInner_DipoleTunnel,rOuter_DipoleTunnel);

	//subtract the tunnel disc
	G4SubtractionSolid* DipoleSolid = new G4SubtractionSolid("DipoleSolid",
		DipoleWithSides,DipoleTunnelCone);

	G4LogicalVolume* DipoleLogical = new G4LogicalVolume(DipoleSolid,
		mMaterialManager->siliconsteel,"DipoleLogical",0,0,0);

	DipoleLogical->SetVisAttributes(OrangeVisAtt); 

	G4RotationMatrix *pRotDipoleInContainer=new G4RotationMatrix();
	pRotDipoleInContainer->rotateY(90*deg); 

	if(mSetupHMS>=4)
	{
		new G4PVPlacement(pRotDipoleInContainer,
			G4ThreeVector(0,-pDipoleRCenterZ,pDipoleRCenterY),
			DipoleLogical,"DipolePhys",HMSContainerLogical,0,0,0);
	}


	//////////////////////////////////////////////////////////

	return theHMSPhys;
}


/////////////////////////////////////////////////////////////////////

