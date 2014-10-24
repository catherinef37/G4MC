// ********************************************************************
//
// $Id: BigBiteDetectorConstruction.cc,v 1.02, 2012/09/26 HRS Exp $
// --------------------------------------------------------------
//
// ********************************************************************
//
#include "BigBiteDetectorConstruction.hh"

#include "G4FieldManager.hh"
#include "G4MagneticField.hh"
#include "G4ChordFinder.hh"
#include "G4TransportationManager.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

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
#include "UsageManager.hh"
#include "HRSMaterial.hh"

//The G4RotationMatrix rotates the whole coordinate system, 
//Looking from top, it always rotate clockwise  

//for bigbite detector, to save time, not build the individual bar, 
//just build one layer of SC
#define DoNotBuildSCBar 1

extern UsageManager* gConfig;

////////////////////////////////////////////////////////////////////////////////////////////////
BigBiteDetectorConstruction::BigBiteDetectorConstruction(G4LogicalVolume *mother): 
logMotherVol(mother) 
{
	GetConfig();
}

BigBiteDetectorConstruction::~BigBiteDetectorConstruction()
{
	;
}

////////////////////////////////////////////////////////////////////////////////////////////////
void BigBiteDetectorConstruction::GetConfig()
{
	gConfig->ReadFile("Detector_BigBite.ini");

	gConfig->GetParameter("PivotXOffset",mPivotXOffset);
	mPivotXOffset*=mm;
	gConfig->GetParameter("PivotYOffset",mPivotYOffset);
	mPivotYOffset*=mm;
	gConfig->GetParameter("PivotZOffset",mPivotZOffset);
	mPivotZOffset*=mm;

	mSetupBigBite=0;
	gConfig->GetParameter("SetupBigBite",mSetupBigBite);

	gConfig->GetParameter("BigBiteAngle",mBigBiteAngle);
	mBigBiteAngle*=deg;
	gConfig->GetParameter("BigBiteTiltAngle",mBigBiteTiltAngle);
	mBigBiteTiltAngle*=deg;
	gConfig->GetParameter("Pivot2BigBiteFace",mPivot2BigBiteFace);
	mPivot2BigBiteFace*=mm;

	mSetupBigBiteSieve=0;
	gConfig->GetParameter("SetupBigBiteSieve",mSetupBigBiteSieve);

	mSetupHAND=0;
	gConfig->GetParameter("SetupHAND",mSetupHAND);
	gConfig->GetParameter("Pivot2HANDLeadWall",mPivot2HANDLeadWall);
	mPivot2HANDLeadWall*=mm;
}


////////////////////////////////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* BigBiteDetectorConstruction::Construct()
{

	G4SDManager* SDMan = G4SDManager::GetSDMpointer();

	G4VSensitiveDetector* SDBBSC1 = new HRSStdSD("BBSC1");
	SDMan->AddNewDetector(SDBBSC1);
	G4VSensitiveDetector* SDBBSC2 = new HRSStdSD("BBSC2");
	SDMan->AddNewDetector(SDBBSC2);

	G4VSensitiveDetector* SDVeto  = new HRSStdSD("VETO");
	SDMan->AddNewDetector(SDVeto);
	G4VSensitiveDetector* SDNDSC1 = new HRSStdSD("NDSC1");
	SDMan->AddNewDetector(SDNDSC1);
	G4VSensitiveDetector* SDNDSC2 = new HRSStdSD("NDSC2");
	SDMan->AddNewDetector(SDNDSC2);
	G4VSensitiveDetector* SDNDSC3 = new HRSStdSD("NDSC3");
	SDMan->AddNewDetector(SDNDSC3);
	G4VSensitiveDetector* SDNDSC4 = new HRSStdSD("NDSC4");
	SDMan->AddNewDetector(SDNDSC4);

	//############################################
	//Declaration and building of all materials used in the BigBite
	//  - vaccum, vide
	//  - hydrogen (for the target and the scintillator and ethane)
	//  - carbon (for scintillator of the scintillators and ethane)
	//  - nitrogen and oxygen (for air)
	//  - helium (for wire chambers)
	//  - iron, lead 
	G4Material *scintillator,*EthaneHe,*air,*vide,*iron,*lead,*copper,*aluminum,*tungsten;

	HRSMaterial* fHRSMaterialManager=HRSMaterial::GetHRSMaterialManager();
	scintillator=fHRSMaterialManager->scintillator;
	EthaneHe=fHRSMaterialManager->EthaneHe;
	air=fHRSMaterialManager->air;
	vide=fHRSMaterialManager->vacuum;
	iron=fHRSMaterialManager->iron;
	lead=fHRSMaterialManager->lead;
	copper=fHRSMaterialManager->copper;
	aluminum=fHRSMaterialManager->aluminum;
	tungsten=fHRSMaterialManager->tungsten;


	// Then ; building of the detection apparatus...
	//############################################

	//the mother box was built along x axis, just need to rotate 180 deg to set it alone -x
	// Since this box is the mother Volume, no translation & rotation.
	double BBArmRotAngle=mBigBiteAngle-90*deg;
	G4RotationMatrix *RotBB=new G4RotationMatrix();
	RotBB->rotateY(-BBArmRotAngle); 


	//############################################

	//////////////////////////
	//The mother box
	//////////////////////////

	//In this way we can place everything into the box as if placed them in the hall 
	double BBContainerX = 7.0*m;  //
	double BBContainerY = 4.0*m;  // Size of this box, large enough but no too big
	double BBContainerZ = 2.0*m;  //  

	G4Box* BBContainer = new G4Box("BBContainer",
		BBContainerX/2+5*cm,BBContainerY/2,BBContainerZ/2);

	G4LogicalVolume* logBBContainer = new G4LogicalVolume(BBContainer, 
		vide, "logBBContainer", 0, 0, 0);
	logBBContainer->SetVisAttributes(HallVisAtt);


	//the position at the hall
	double BBContainerPosX=(mPivot2BigBiteFace+BBContainerX/2)*cos(BBArmRotAngle)+mPivotXOffset;
	double BBContainerPosY=mPivotYOffset;
	double BBContainerPosZ=-(mPivot2BigBiteFace+BBContainerX/2)*sin(BBArmRotAngle)+mPivotZOffset;
	G4PVPlacement* phyBBContainer= new G4PVPlacement(RotBB,
		G4ThreeVector(BBContainerPosX,BBContainerPosY,BBContainerPosZ),
		logBBContainer,"BigBiteContainner",logMotherVol,0,0);


	//////////////////////////
	//The Aperture 
	//////////////////////////
	//the sieve slit is a 4cm thick tungsten block, 0.695m high and 0.24m wide

	double BBSieveX = 0.04*m;
	double BBSieveY = 0.6985*m; 
	double BBSieveZ = 0.24*m+8*cm; 
	G4VSolid* BBAperturebox = new G4Box("BBAperturebox",BBSieveX/2-0.5*mm,1.6*m/2,0.7*m/2);
	G4VSolid* BBSievebox = new G4Box("BBSievebox",BBSieveX/2,BBSieveY/2,BBSieveZ/2);
	G4SubtractionSolid* BBApertureVol = new G4SubtractionSolid("BBApertureVol",
		BBAperturebox,BBSievebox);

	G4LogicalVolume* logBBApertureVol = new G4LogicalVolume(BBApertureVol, 
		aluminum, "logBBApertureVol", 0, 0, 0);
	logBBApertureVol->SetVisAttributes(WhiteVisAtt);

	//xpos relative to its containner(3 cm away from the front face);
	const double sieve2MBFrontFace=0.03*m;
	double BBSievePosX=-BBContainerX/2+sieve2MBFrontFace+BBSieveX/2;
	new G4PVPlacement(0,G4ThreeVector(BBSievePosX,0,0),
		logBBApertureVol,"phyBBApertureVol",logBBContainer,0,0);


	//////////////////////////
	//The sieve slit
	//////////////////////////
	//do we really want to place the sieve, if yes, need to built the sieve holes(d=1.91cm)

	if(mSetupBigBiteSieve)
	{
		//sieve holes uniform distributed, 13(row) x 7 col, edge collums located at +/- 22.86/2 cm
		//gap between rows is 5cm, gap between col is 22.86/6
		//there are 4 enlongated holes, will take care later (see NIM paper) 
		const double dy=5.0*cm,dz=22.86*cm/6;
		const double BBSieveHoleR=1.91*cm/2;

		//subtract 91 holes
		G4VSolid* BBSieveVol[96];
		char VolName[100];
		double holePosZ,holePosY;		
		int idx=0;

		G4Tubs* BBSievehole = new G4Tubs("BBSievehole",
			0,BBSieveHoleR,BBSieveX/2+1*mm,0.0,360.0*deg);
		G4RotationMatrix* RotY90deg=new G4RotationMatrix();
		RotY90deg->rotateY(-90*deg);

		//the holes might be tilted.....need to check
		for(int i=0;i<7;i++)
		{
			holePosZ=(i-3)*dz;
			for(int j=0;j<13;j++)
			{
				holePosY=(j-6)*dy;
				sprintf(VolName,"sieve_subhole_%02d",idx+1);
				if(idx==0)
				{
					BBSieveVol[idx] = new G4SubtractionSolid(VolName,BBSievebox,
						BBSievehole,RotY90deg,G4ThreeVector(0,holePosY,holePosZ));
				}
				else
				{
					BBSieveVol[idx] = new G4SubtractionSolid(VolName,BBSieveVol[idx-1],
						BBSievehole,RotY90deg,G4ThreeVector(0,holePosY,holePosZ));
				}
				idx++;
			}
		}


		//subtract longated holes,
		G4Box* Sieve2HoleBox = new G4Box("Sieve2HoleBox",BBSieveX/2+1*mm,BBSieveHoleR,dz/2);
		//subtract 5 times at the following location...
		//j=3, i=1-2  and i=4-5
		//j=9, i=2-3-4
		//j=11, i=1-2
		int jj[5]={3,3,9,9,11};	//row index
		int ii[5]={1,4,2,3,1};	//col index
		for(idx=91;idx<96;idx++)
		{
			holePosY=(jj[idx-91]-6)*dy;
			holePosZ=(ii[idx-91]-3+0.5)*dz;
			sprintf(VolName,"sieve_subhole_%02d",idx+1);
			BBSieveVol[idx] = new G4SubtractionSolid(VolName,BBSieveVol[idx-1],
				Sieve2HoleBox,0,G4ThreeVector(0,holePosY,holePosZ));
		}

		G4LogicalVolume* logBBSieveVol = new G4LogicalVolume(BBSieveVol[95], 
			tungsten, "logBBSieveVol", 0, 0, 0);
		logBBSieveVol->SetVisAttributes(YellowVisAtt);

		new G4PVPlacement(0,G4ThreeVector(BBSievePosX,0,0),
			logBBSieveVol,"phyBBSieveVol",logBBContainer,0,0);
	}


	//############################################	
	//some constants
	const double BBCoilAngle=25*deg;    // angle of trapzoid side for the coil
	G4RotationMatrix *RotBBCoil=new G4RotationMatrix();
	RotBBCoil->rotateZ(-BBCoilAngle);   // (-) in mother box coordinate system

	const double kCoilThickness=28*cm;
	const double kCoilWidth=20*cm;
	const double kCoilGap=36*cm;		//gap between 2 coils
	const double kFieldGap=28*cm;	//gap of the field


	//////////////////////////
	//the inner iron of the magnet
	//////////////////////////
	//This is the size of the inner iron, which is a normal trap union a triangle
	//this triangle can be part of a rotated rectangle
	//this solid should be subtracted from the coil solid to form the coil geometry
	//The coil shape can be achieved by just extending the size here by 2xcoilthickness 
	//The effective field Volumn, similarly, can be achieved by extending the size 
	//by 1xcoilthickness  
	//size of inner iron
	double BBInnerIronX = 752.9*mm;         // length along the track (bottom)
	double BBInnerIronZ = kCoilWidth+(kCoilGap-kFieldGap)/2;  // transverse size
	double BBInnerIronY = 817.4*mm;         // vertical seize of opening apperture
	double BBInnerIronX2 = BBInnerIronX-BBInnerIronY*tan(BBCoilAngle); // length (top)//BBInnerIronX2=371.74*mm
	//the front face of the inner iron is 0.5m from the front face of the mother box
	const double BBField2MBFrontFace=0.5*m+(BBInnerIronX+BBInnerIronX2)/2/2;
	const double BBFiledPosX = -BBContainerX/2+BBField2MBFrontFace;

	// built the inner iron geometry
	//First, the big trapezoid ;
	G4Trap* BBInnerIrontrap = new G4Trap("BBInnerIrontrap",
		BBInnerIronZ, BBInnerIronY, BBInnerIronX, BBInnerIronX2);
	//Second, the missing part (tilted box)
	G4Box* BBInnerIronbox = new G4Box("BBInnerIronbox",
		BBInnerIronX2*cos(BBCoilAngle)/2,BBInnerIronX2*sin(BBCoilAngle)/2,BBInnerIronZ/2);
	//Then, we have to move & rotate the box when merge the compounds
	G4UnionSolid* BBInnerIronVol = new G4UnionSolid("BBInnerIronVol",
		BBInnerIrontrap, BBInnerIronbox,
		RotBBCoil,G4ThreeVector((BBInnerIronX2-BBInnerIronX)/2/2,BBInnerIronY/2,0.));

	//Then we create the logical Volume
	G4LogicalVolume* logBBInnerIronVol = new G4LogicalVolume(BBInnerIronVol,        
		iron, "logBBInnerIronVol",0,0,0);
	logBBInnerIronVol->SetVisAttributes(DarkBlueVisAtt);

	//place the iron 
	new G4PVPlacement(0,G4ThreeVector(BBFiledPosX, 0., kFieldGap/2+BBInnerIronZ/2),
		logBBInnerIronVol,"physBBInnerIronVol",logBBContainer,1,0);
	new G4PVPlacement(0,G4ThreeVector(BBFiledPosX, 0., -(kFieldGap/2+BBInnerIronZ/2)),
		logBBInnerIronVol,"physBBInnerIronVol",logBBContainer,1,1);


	//////////////////////////
	//The field Volumn
	//////////////////////////
	//size of effective field
	double BBFieldX = BBInnerIronX + kCoilThickness;  // length along the track (bottom)
	double BBFieldZ = kFieldGap;             // transverse size
	double BBFieldY = BBInnerIronY + kCoilThickness;  // vertical seize of opening apperture
	double BBFieldX2 = BBFieldX-BBFieldY*tan(BBCoilAngle);	 // length (top)


	// built the effective field geometry
	//First, the big trapezoid ;
	G4Trap* BBFieldtrap = new G4Trap("BBFieldtrap",BBFieldZ, BBFieldY, BBFieldX, BBFieldX2);
	//Second, the missing part (tilted box)
	G4Box* BBFieldbox = new G4Box("BBFieldbox",
		BBFieldX2*cos(BBCoilAngle)/2,BBFieldX2*sin(BBCoilAngle)/2,BBFieldZ/2);
	//Then, we have to move & rotate the box when merge the compounds
	G4UnionSolid* BBFieldVol = new G4UnionSolid("BBFieldVol",
		BBFieldtrap, BBFieldbox,
		RotBBCoil,G4ThreeVector((BBFieldX2-BBFieldX)/2/2,BBFieldY/2,0.));

	//############################################	
	//Here we create a magnetic field inside this Volume, do not use the global field manager
	//otherwise this field will be applied to the whole world
	double BBFieldValue = -0.92*tesla;
	double BBField2Pivot = mPivot2BigBiteFace + 0.4*m;
	BigBiteMagField *BBMField = new BigBiteMagField(mBigBiteAngle,BBFieldValue,BBField2Pivot,
		mPivotXOffset,mPivotYOffset,mPivotZOffset);
	//G4FieldManager* BBMFieldMgr = 
	//	G4TransportationManager::GetTransportationManager()->GetFieldManager();
	G4FieldManager* BBMFieldMgr = new G4FieldManager();
	BBMFieldMgr->SetDetectorField(BBMField);
	//For the trajectory calculation
	BBMFieldMgr->CreateChordFinder(BBMField);
	BBMFieldMgr->GetChordFinder()->SetDeltaChord(1.e-2);
	//############################################	

	//should we use air or helium here?
	G4LogicalVolume* logBBFieldVol = new G4LogicalVolume(BBFieldVol,        
		air, "logBBFieldVol",BBMFieldMgr,0,0);
	logBBFieldVol->SetVisAttributes(PurpleVisAtt);  

	//Then we make the Physical Volume & put it in the mother box
	new G4PVPlacement(0,G4ThreeVector(BBFiledPosX, 0., 0.),
		logBBFieldVol,"physBBFieldVol",logBBContainer,0,0);


	//////////////////////////
	//The coil Volumn
	//////////////////////////
	//size of coil
	double BBCoilX = BBInnerIronX + 2*kCoilThickness;    // length along the track (bottom)
	double BBCoilZ = kCoilWidth;                         // transverse size
	double BBCoilY = BBInnerIronY + 2*kCoilThickness;    // vertical seize of opening apperture
	double BBCoilX2 = BBCoilX-BBCoilY*tan(BBCoilAngle);	 // length (top)

	//built the coil geometry
	//First, the big trapezoid ;
	G4Trap* BBCoiltrap = new G4Trap("BBCoiltrap",BBCoilZ, BBCoilY, BBCoilX, BBCoilX2);
	//Second, the missing part (tilted box)
	G4Box* BBCoilbox = new G4Box("BBCoilbox",
		BBCoilX2*cos(BBCoilAngle)/2,BBCoilX2*sin(BBCoilAngle)/2,BBCoilZ/2);
	//Then, we have to move & rotate the box when merge the compounds
	G4UnionSolid* BBCoilshape = new G4UnionSolid("BBCoilshape",
		BBCoiltrap, BBCoilbox,
		RotBBCoil,G4ThreeVector((BBCoilX2-BBCoilX)/2/2,BBCoilY/2,0.));
	//then subtract the inner iron from the coil
	G4SubtractionSolid* BBCoilVol = new G4SubtractionSolid("BBCoilVol",
		BBCoilshape,BBInnerIrontrap);

	G4LogicalVolume* logBBCoilVol = new G4LogicalVolume(BBCoilVol,        
		copper, "logBBCoilVol",0,0,0);
	logBBCoilVol->SetVisAttributes(CuBrownVisAtt);

	new G4PVPlacement(0,G4ThreeVector(BBFiledPosX, 0., kCoilGap/2+BBCoilZ/2),
		logBBCoilVol,"physBBCoilVol",logBBContainer,0,0);
	new G4PVPlacement(0,G4ThreeVector(BBFiledPosX, 0., -(kCoilGap/2+BBCoilZ/2)),
		logBBCoilVol,"physBBCoilVol",logBBContainer,0,1);


	//////////////////////////
	//the outter iron of the magnet
	//////////////////////////
	double BBOuterIronX=(BBInnerIronX+BBInnerIronX2)/2;
	double BBOuterIronY=BBCoilY+BBCoilX2*cos(BBCoilAngle)*sin(BBCoilAngle)+60*cm;
	double BBOuterIronZ=BBInnerIronZ*2+kCoilGap+60*cm;

	G4Box* BBOuterIronshape = new G4Box("BBOuterIronshape",
		BBOuterIronX/2,BBOuterIronY/2,BBOuterIronZ/2);

	//built the larger coil geometry then subtract it
	//It has to be a little bit larger ...
	double BBLargeCoilX = BBCoilX + 1*cm;        // length along the track (bottom)
	double BBLargeCoilZ = kCoilGap + 2*kCoilWidth; // transverse size
	double BBLargeCoilY = BBCoilY + 1*cm;        // vertical seize of opening apperture
	double BBLargeCoilX2 = BBLargeCoilX-BBLargeCoilY*tan(BBCoilAngle);	 // length (top)

	G4Trap* BBLargeCoiltrap = new G4Trap("BBLargeCoiltrap_sub",
		BBLargeCoilZ, BBLargeCoilY, BBLargeCoilX, BBLargeCoilX2);
	G4Box* BBLargeCoilbox = new G4Box("BBLargeCoilbox",
		BBLargeCoilX2*cos(BBCoilAngle)/2,BBLargeCoilX2*sin(BBCoilAngle)/2,BBLargeCoilZ/2);
	G4UnionSolid* BBLargeCoilshape = new G4UnionSolid("BBLargeCoilshape",
		BBLargeCoiltrap, BBLargeCoilbox,
		RotBBCoil,G4ThreeVector((BBLargeCoilX2-BBLargeCoilX)/2/2,BBLargeCoilY/2,0.));


	//Have to subtract part of back since the coil shape is thinner than the iron
	double BBIronBoxsubx=BBInnerIronX2;
	double BBIronBoxsuby=BBCoilY/2+BBCoilX2*cos(BBCoilAngle)*sin(BBCoilAngle);
	double BBIronBoxsubz=BBLargeCoilZ;
	G4Box *BBIronBoxsub=new G4Box("BBIronBoxsub",
		BBIronBoxsubx/2,BBIronBoxsuby/2,BBIronBoxsubz/2);

	double BBIronBoxsub_PosX=BBIronBoxsubx/2 + 
		(BBCoilX2*cos(BBCoilAngle)*cos(BBCoilAngle)-(BBCoilX+BBCoilX2)/2/2);  //relative to the field center
	G4SubtractionSolid* BBOuterIron_subbox = new G4SubtractionSolid(
		"BBOuterIron_subbox",BBOuterIronshape,BBIronBoxsub,
		0,G4ThreeVector(BBIronBoxsub_PosX,BBIronBoxsuby/2,0));

	///////////////debug this box
	////Then we create the logical Volume
	//G4LogicalVolume* logBBIronBoxVol = new G4LogicalVolume(BBIronBoxsub,        
	//	iron, "logBBIronBoxVol",0,0,0);
	//logBBIronBoxVol->SetVisAttributes(DarkBlueVisAtt);
	////place the iron 
	//new G4PVPlacement(0,G4ThreeVector(BBFiledPosX+BBIronBoxsub_PosX,BBIronBoxy/2,0),
	//	logBBIronBoxVol,"physBBIronBoxVol",logBBContainer,0,0);
	/////////////////debug this box

	//subtract the coil shape
	G4SubtractionSolid* BBOuterIron_subboxcoil = new G4SubtractionSolid(
		"BBOuterIron_subboxcoil",BBOuterIron_subbox,BBLargeCoilshape,
		0,G4ThreeVector(0,0,0));

	//for some unknown reasons, If I union the inner iron here, 
	//the iron can not build ......I gave this way up
	//union the inner iron
	//G4UnionSolid* BBIron_1 = new G4UnionSolid("BBIron_1",
	//	BBOuterIron_subboxcoil,BBInnerIronVol,0,G4ThreeVector(0,0,kFieldGap/2+BBInnerIronZ/2));
	//G4UnionSolid* BBIronVol = new G4UnionSolid("BBIronVol",
	//	BBIron_1,BBInnerIronVol,0,G4ThreeVector(0,0,-(kFieldGap/2+BBInnerIronZ/2)));


	//Then we create the logical Volume
	G4LogicalVolume* logBBIronVol = new G4LogicalVolume(BBOuterIron_subboxcoil,        
		iron, "logBBIronVol",0,0,0);
	logBBIronVol->SetVisAttributes(DarkBlueVisAtt);

	//place the iron 
	new G4PVPlacement(0,G4ThreeVector(BBFiledPosX, 0., 0),
		logBBIronVol,"physBBIronVol",logBBContainer,0,0);

	//////////////////////////////////////////////////////
	//here we built the bigbite detector box, which contains 4 layers, MWDC1, MQDC2, dE and E
	//the first two is drift chamber, which can be ignored
	//dE and E layer are scintillators, both contains 24 bars of 50cm long and 8.6cm wide
	//dE layer is 0.3cm thick and E layer is 3cm thick

	//relative position information;
	//MWDC1 front face is 0.27m/cos(mBigBiteTiltAngle)=63.9cm to the magnetic field center(BBFiledPosX)
	//where 0.27m is the height of MWDC1 center at its front face
	//MWDC2's front face is 0.757m from the front face of MWDC1
	//dE layer front face is 0.973m from the front face of MWDC1
	//E layer front face is 1.023m from the front face of MWDC1
	//////////////////////////////////////////////////////

	//The rotation of the detector could be 25, 16, or 10 deg
	//const double mBigBiteTiltAngle=16*deg;    // angle of BigBite planes

	G4RotationMatrix *RotBBDet=new G4RotationMatrix();
	RotBBDet->rotateZ(-mBigBiteTiltAngle); 

	//////////////////////////
	//the detector container
	//////////////////////////
	const double BBField2MWDC1=0.27*m/cos(mBigBiteTiltAngle)+40*cm;  //MWDC1 front face to field center
	const double BBField2MWDC2=BBField2MWDC1+0.757*m;
	const double BBField2dEPlane=BBField2MWDC1+0.973*m;
	const double BBField2EPlane=BBField2MWDC1+1.023*m;	
	//put the detector box's front face just 1 cm before MWDC1 such that this container will
	//not overlap with the iron
	const double BBField2DetContainer=BBField2MWDC1-1*cm;
	const double MWDCThick=7*cm;

	double BBDetContainerX=1.2*m;
	double BBDetContainerY=2.4*m;
	double BBDetContainerZ=1.2*m;
	double BBDetContainerPosX=BBFiledPosX+
		(BBField2DetContainer+BBDetContainerX/2)*cos(mBigBiteTiltAngle);
	double BBDetContainerPosY=(BBField2DetContainer+BBDetContainerX/2)*sin(mBigBiteTiltAngle);
	G4Box *BBDetContainer=new G4Box("BBDetContainer",
		BBDetContainerX/2,BBDetContainerY/2,BBDetContainerZ/2);

	//Then we create the logical Volume
	G4LogicalVolume* logBBDetContainer = new G4LogicalVolume(BBDetContainer,        
		air, "logBBDetContainer",0,0,0);
	logBBDetContainer->SetVisAttributes(HallVisAtt);

	//place the containner
	new G4PVPlacement(RotBBDet,G4ThreeVector(BBDetContainerPosX, BBDetContainerPosY, 0),
		logBBDetContainer,"physBBDetContainer",logBBContainer,0,0);


	//////////////////////////
	//the MWDC1
	//////////////////////////

	double BBMWDC1X=MWDCThick;
	double BBMWDC1Y=1.4*m;
	double BBMWDC1Z=0.35*m;
	double BBMWDC1PosX=BBField2MWDC1-(BBField2DetContainer+BBDetContainerX/2)
		+BBMWDC1X/2;  
	//relative to its container
	G4Box *BBMWDC1Vol=new G4Box("BBMWDC1Vol",BBMWDC1X/2,BBMWDC1Y/2,BBMWDC1Z/2);

	//Then we create the logical Volume
	G4LogicalVolume* logBBMWDC1Vol = new G4LogicalVolume(BBMWDC1Vol,        
		EthaneHe, "logBBMWDC1Vol",0,0,0);
	logBBMWDC1Vol->SetVisAttributes(PurpleVisAtt);

	//place this Volume
	new G4PVPlacement(0,G4ThreeVector(BBMWDC1PosX, 0, 0),
		logBBMWDC1Vol,"physBBMWDC1Vol",logBBDetContainer,0,0);


	//////////////////////////
	//the MWDC2
	//////////////////////////

	double BBMWDC2X=MWDCThick;
	double BBMWDC2Y=0.086*m*24;
	double BBMWDC2Z=0.5*m;
	double BBMWDC2PosX=BBField2MWDC2-(BBField2DetContainer+BBDetContainerX/2)
		+BBMWDC2X/2;		
	//relative to its container
	G4Box *BBMWDC2Vol=new G4Box("BBMWDC2Vol",BBMWDC2X/2,BBMWDC2Y/2,BBMWDC2Z/2);

	//Then we create the logical Volume
	G4LogicalVolume* logBBMWDC2Vol = new G4LogicalVolume(BBMWDC2Vol,        
		EthaneHe, "logBBMWDC2Vol",0,0,0);
	logBBMWDC2Vol->SetVisAttributes(PurpleVisAtt);

	//place this Volume
	new G4PVPlacement(0,G4ThreeVector(BBMWDC2PosX, 0, 0),
		logBBMWDC2Vol,"physBBMWDC2Vol",logBBDetContainer,0,0);



	//////////////////////////
	//the dE plane (SC1 plane)
	//////////////////////////
	//24 bars of 0.3 x 8.6 x 50 cm

	//We can place 24 SC bars here, which will reduce the work on digitization
	//but make this program runs slow, especially the visulization
	//Or we can place one big layer here

#ifdef DoNotBuildSCBar

	double BBdEPlaneX=3.0*cm;
	double BBdEPlaneY=8.6*cm*24;
	double BBdEPlaneZ=50*cm;
	G4Box *BBdEPlaneVol=new G4Box("BBdEPlaneVol",BBdEPlaneX/2,BBdEPlaneY/2,BBdEPlaneZ/2);

	//Then we create the logical Volume
	G4LogicalVolume* logBBdEPlaneVol = new G4LogicalVolume(BBdEPlaneVol,        
		scintillator, "logBBdEPlaneVol",0,0,0);
	logBBdEPlaneVol->SetVisAttributes(LightYellowVisAtt);

	//place this Volume
	//position relative to its container
	double BBdEPlanePosX=BBField2dEPlane-(BBField2DetContainer+BBDetContainerX/2)
		+BBdEPlaneX/2;		
	new G4PVPlacement(0,G4ThreeVector(BBdEPlanePosX, 0, 0),
		logBBdEPlaneVol,"physBBdEPlaneVol",logBBDetContainer,0,0);

#else

	bool bBuildSCBar=true;

	const int BBSC1BarNum=24;
	double BBSC1BarX=0.3*cm, BBSC1BarY=8.6*cm, BBSC1BarZ=50*cm;
	G4Box *BBSC1Bar=new G4Box("BBSC1Bar",BBSC1BarX/2,BBSC1BarY/2,BBSC1BarZ/2);

	double BBSC1VolX=BBSC1BarX, BBSC1VolY=BBSC1BarY*BBSC1BarNum, BBSC1VolZ=BBSC1BarZ;
	G4Box *BBSC1Vol=new G4Box("BBSC1Vol",BBSC1VolX/2,BBSC1VolY/2,BBSC1VolZ/2);

	G4LogicalVolume* logBBSC1Vol;
	if(bBuildSCBar)
	{
		logBBSC1Vol= new G4LogicalVolume(BBSC1Vol, 
			vide, "logBBSC1Vol", 0, 0, 0);
		logBBSC1Vol->SetVisAttributes(HallVisAtt);

		G4LogicalVolume* logBBSC1Bar = new G4LogicalVolume(BBSC1Bar,
			scintillator, "logBBSC1Bar",0,0,0);
		logBBSC1Bar->SetVisAttributes(LightGreenVisAtt);
		logBBSC1Bar->SetSensitiveDetector(SDBBSC1);

		// Placement of the 24 bars vith replication
		new G4PVReplica("replicaBBSC1Vol",logBBSC1Bar,logBBSC1Vol,
			kYAxis,BBSC1BarNum,BBSC1BarY,-BBSC1VolY/2+0.5*BBSC1BarY);
	}
	else
	{
		logBBSC1Vol = new G4LogicalVolume(BBSC1Vol, 
			scintillator, "logBBSC1Vol", 0, 0, 0);
		logBBSC1Vol->SetVisAttributes(LightGreenVisAtt);
		logBBSC1Vol->SetSensitiveDetector(SDBBSC1);
	}

	//position relative to its container
	double BBSC1VolPosX=BBField2dEPlane-(BBField2DetContainer+BBDetContainerX/2)
		+BBSC1VolX/2;
	new G4PVPlacement(0,G4ThreeVector(BBSC1VolPosX, 0, 0),
		logBBSC1Vol,"physBBSC1Vol",logBBDetContainer,0,0);

#endif


	//////////////////////////
	//the E plane (SC2 plane)
	//////////////////////////
	//24 bars of 3 x 8.6 x 50 cm

	//We can place 24 SC bars here, which will reduce the work on digitization
	//but make this program runs slow, especially the visulization
	//Or we can place one big layer here

#ifdef DoNotBuildSCBar

	double BBEPlaneX=3.0*cm;
	double BBEPlaneY=8.6*cm*24;
	double BBEPlaneZ=50*cm;
	G4Box *BBEPlaneVol=new G4Box("BBEPlaneVol",BBEPlaneX/2,BBEPlaneY/2,BBEPlaneZ/2);

	//Then we create the logical Volume
	G4LogicalVolume* logBBEPlaneVol = new G4LogicalVolume(BBEPlaneVol,        
		scintillator, "logBBEPlaneVol",0,0,0);
	logBBEPlaneVol->SetVisAttributes(LightYellowVisAtt);

	//place this Volume
	//position relative to its container
	double BBEPlanePosX=BBField2EPlane-(BBField2DetContainer+BBDetContainerX/2)
		+BBEPlaneX/2;		
	new G4PVPlacement(0,G4ThreeVector(BBEPlanePosX, 0, 0),
		logBBEPlaneVol,"physBBEPlaneVol",logBBDetContainer,0,0);

#else

	const int BBSC2BarNum=24;
	double BBSC2BarX=3*cm, BBSC2BarY=8.6*cm, BBSC2BarZ=50*cm;
	G4Box *BBSC2Bar=new G4Box("BBSC2Bar",BBSC2BarX/2,BBSC2BarY/2,BBSC2BarZ/2);

	double BBSC2VolX=BBSC2BarX, BBSC2VolY=BBSC2BarY*BBSC2BarNum, BBSC2VolZ=BBSC2BarZ;
	G4Box *BBSC2Vol=new G4Box("BBSC2Vol",BBSC2VolX/2,BBSC2VolY/2,BBSC2VolZ/2);

	G4LogicalVolume* logBBSC2Vol;
	if(bBuildSCBar)
	{
		logBBSC2Vol= new G4LogicalVolume(BBSC2Vol, 
			vide, "logBBSC2Vol", 0, 0, 0);
		logBBSC2Vol->SetVisAttributes(HallVisAtt);

		G4LogicalVolume* logBBSC2Bar = new G4LogicalVolume(BBSC2Bar,
			scintillator, "logBBSC2Bar",0,0,0);
		logBBSC2Bar->SetVisAttributes(LightGreenVisAtt);
		logBBSC2Bar->SetSensitiveDetector(SDBBSC2);

		// Placement of the 24 bars vith replication
		new G4PVReplica("replicaBBSC2Vol",logBBSC2Bar,logBBSC2Vol,
			kYAxis,BBSC2BarNum,BBSC2BarY,-BBSC2VolY/2+0.5*BBSC2BarY);
	}
	else
	{
		logBBSC2Vol = new G4LogicalVolume(BBSC2Vol, 
			scintillator, "logBBSC2Vol", 0, 0, 0);
		logBBSC2Vol->SetVisAttributes(LightGreenVisAtt);
		logBBSC2Vol->SetSensitiveDetector(SDBBSC2);
	}

	//position relative to its container
	double BBSC2VolPosX=BBField2EPlane-(BBField2DetContainer+BBDetContainerX/2)
		+BBSC2VolX/2;
	new G4PVPlacement(0,G4ThreeVector(BBSC2VolPosX, 0, 0),
		logBBSC2Vol,"physBBSC2Vol",logBBDetContainer,0,0);

#endif

	//////////////////////////////////////////
	//From here we built the Hall A Neutron Detector
	//HAND is made of 88 main detecting bars arranged in four layers. The thickness 
	//of each bar in these layers is 10 cm, the length is 100 cm, and the height 
	//varies with the smaller bars placed in front of the larger bars. There is also 
	//a thinner "Veto" layer that contains 64 bars with dimensions of 2 x 11 x 70 cm.
	//It will run very slow if we build all the bars here

	//the best way is just built one layer then digitize it according to its position
	//at EndOfEventAction

	//////////////////////////////////////////
	//The front face of the lead wall to the target is 4.4m, or 3.3m to the mother box front face
	//Here the container of HAND is 10 cm in front of the lead wall

	//const double mPivot2HANDLeadWall=3.3*m+mPivot2BigBiteFace;  //4.4*m;

	const double LeadWall2NDContainer=10*cm;
	const double Veto2NDContainer=40*cm;   //Veto1 front face to the NDContainer front face
	const double NDContainer2MBFrontFace=mPivot2HANDLeadWall-mPivot2BigBiteFace-LeadWall2NDContainer;   
	const double BBField2NDContainer = NDContainer2MBFrontFace-BBField2MBFrontFace; // usefull later

	//////////////////////////////////////////

	if(mSetupHAND)
	{
		//////////////////////////
		//the neutron detector container
		//////////////////////////
		double NDContainerX=1.2*m;
		double NDContainerY=3.6*m;
		double NDContainerZ=2.6*m;
		G4Box *NDContainer=new G4Box("NDContainer",
			NDContainerX/2,NDContainerY/2,NDContainerZ/2);

		//Then we create the logical Volume
		G4LogicalVolume* logNDContainer = new G4LogicalVolume(NDContainer,        
			air, "logNDContainer",0,0,0);
		logNDContainer->SetVisAttributes(HallVisAtt);

		//place the HAND containner
		double NDContainerPosX=BBFiledPosX+BBField2NDContainer+NDContainerX/2;
		new G4PVPlacement(0,G4ThreeVector(NDContainerPosX, 0, 0),
			logNDContainer,"physNDContainer",logBBContainer,0,0);


		//////////////////////////
		//the iron+lead+iron wall 
		//////////////////////////
		//a 9.08 cm thick wall, made up of 4 cm of iron casing surrounding the 5.08 cm thick lead
		//I am going to place here 2cm iron + 5.08cm lead + 2cm iron 


		double HalfLeadWallX=25.4*mm;  //Half thickness of the lead wall 
		double HalfLeadWallY=1775.*mm; //Half height of wall
		double HalfLeadWallZ=1200.*mm; //Half width of wall

		G4Box* LeadWallVol = new G4Box("LeadWallVol",
			HalfLeadWallX,HalfLeadWallY,HalfLeadWallZ);

		G4LogicalVolume* logLeadWallVol = new G4LogicalVolume(LeadWallVol,
			lead, "logLeadWallVol",  0, 0, 0);
		logLeadWallVol->SetVisAttributes(GrayVisAtt);

		double HalfIronWallX=10.0*mm;  //Half thickness of iron wall  
		G4Box* IronWallVol = new G4Box("IronWallVol",
			HalfIronWallX,HalfLeadWallY,HalfLeadWallZ);

		G4LogicalVolume* logIronWallVol = new G4LogicalVolume(IronWallVol,
			iron, "logIronWallVol",  0, 0, 0);
		logIronWallVol->SetVisAttributes(GrayVisAtt);

		double IronWall1PosX=-NDContainerX/2+LeadWall2NDContainer+HalfIronWallX;
		double LeadWallPosX=IronWall1PosX+HalfIronWallX+HalfLeadWallX;
		double IronWall2PosX=LeadWallPosX+HalfLeadWallX+HalfIronWallX;

		//place 2cm iron + 5.08cm lead + 2cm iron 
		new G4PVPlacement(0,G4ThreeVector(IronWall1PosX,0,0),
			logIronWallVol,"phyIronWall1Vol",logNDContainer,0,0);
		new G4PVPlacement(0,G4ThreeVector(LeadWallPosX,0,0),
			logLeadWallVol,"phyLeadWallVol",logNDContainer,0,0);
		new G4PVPlacement(0,G4ThreeVector(IronWall2PosX,0,0),
			logIronWallVol,"phyIronWall2Vol",logNDContainer,0,0);


#ifdef DoNotBuildSCBar

		//////////////////////////
		//the Veto plane  
		//////////////////////////

		double VetoBarX=2*cm, VetoBarY=11*cm, VetoBarZ=70*cm;
		double VetoVolX=VetoBarX, VetoVolZ=VetoBarZ*2;
		double  Veto1VolY=12*VetoBarY, Veto2VolY=10*VetoBarY;
		G4Box* Veto1Vol = new G4Box("Veto1Vol",VetoVolX/2,Veto1VolY/2,VetoVolZ/2);
		G4Box* Veto2Vol = new G4Box("Veto2Vol",VetoVolX/2,Veto2VolY/2,VetoVolZ/2);

		G4LogicalVolume* logVeto1Vol = new G4LogicalVolume(Veto1Vol,
			scintillator, "logVeto1Vol",  0, 0, 0);
		logVeto1Vol->SetVisAttributes(PurpleVisAtt);
		logVeto1Vol->SetSensitiveDetector(SDVeto);

		G4LogicalVolume* logVeto2Vol = new G4LogicalVolume(Veto2Vol,
			scintillator, "logVeto2Vol",  0, 0, 0);
		logVeto2Vol->SetVisAttributes(PurpleVisAtt);
		logVeto2Vol->SetSensitiveDetector(SDVeto);

		//place the Veto1 and Veto2 into the container	
		double Veto1PosX=-NDContainerX/2+Veto2NDContainer+VetoVolX/2;	
		double Veto2PosX=Veto1PosX+VetoVolX/2+2*mm+VetoVolX/2;
		double Veto2PosY=1.5*m-5*VetoBarY;
		new G4PVPlacement( 0,G4ThreeVector(Veto1PosX,0,0),
			logVeto1Vol,"phyVeto1Vol",logNDContainer,0,0);
		new G4PVPlacement(0,G4ThreeVector(Veto2PosX,Veto2PosY,0),
			logVeto2Vol,"phyVeto2UpVol",logNDContainer,0,0);
		new G4PVPlacement(0,G4ThreeVector(Veto2PosX,-Veto2PosY,0),
			logVeto2Vol,"phyVeto2DownVol",logNDContainer,0,0);

		//////////////////////////
		//the 4 sc planes
		//////////////////////////

		double NDSCVolX = 10*cm, NDSCVolY = 300*cm, NDSCVolZ = 100*cm;
		G4Box* NDSCVol = new G4Box("NDSCVol",NDSCVolX/2,NDSCVolY/2,NDSCVolZ/2);

		G4LogicalVolume* logNDSC1Vol = new G4LogicalVolume(NDSCVol, 
			scintillator, "logNDSC1Vol", 0, 0, 0);
		logNDSC1Vol->SetVisAttributes(DarkBlueVisAtt);
		logNDSC1Vol->SetSensitiveDetector(SDNDSC1);

		double NDSC1PosX=Veto2PosX+VetoVolX/2+NDSCVolX/2+2*mm;	
		new G4PVPlacement(0,G4ThreeVector(NDSC1PosX,0,0),
			logNDSC1Vol,"phyNDSC1Vol",logNDContainer,0,0);


		G4LogicalVolume* logNDSC2Vol = new G4LogicalVolume(NDSCVol, 
			scintillator, "logNDSC2Vol", 0, 0, 0);
		logNDSC2Vol->SetVisAttributes(LightYellowVisAtt);
		logNDSC2Vol->SetSensitiveDetector(SDNDSC2);

		double NDSC2PosX=NDSC1PosX+NDSCVolX+2*mm;	
		new G4PVPlacement(0,G4ThreeVector(NDSC2PosX,0,0),
			logNDSC2Vol,"phyNDSC2Vol",logNDContainer,0,0);

		G4LogicalVolume* logNDSC3Vol = new G4LogicalVolume(NDSCVol, 
			scintillator, "logNDSC3Vol", 0, 0, 0);
		logNDSC3Vol->SetVisAttributes(RedVisAtt);
		logNDSC3Vol->SetSensitiveDetector(SDNDSC3);

		double NDSC3PosX=NDSC2PosX+NDSCVolX+2*mm;	
		new G4PVPlacement(0,G4ThreeVector(NDSC3PosX,0,0),
			logNDSC3Vol,"phyNDSC3Vol",logNDContainer,0,0);

		G4LogicalVolume* logNDSC4Vol = new G4LogicalVolume(NDSCVol, 
			scintillator, "logNDSC4Vol", 0, 0, 0);
		logNDSC4Vol->SetVisAttributes(LightGreenVisAtt);
		logNDSC4Vol->SetSensitiveDetector(SDNDSC4);

		double NDSC4PosX=NDSC3PosX+NDSCVolX+2*mm;	
		new G4PVPlacement(0,G4ThreeVector(NDSC4PosX,0,0),
			logNDSC4Vol,"phyNDSC4Vol",logNDContainer,0,0);

#else

		//////////////////////////
		//the Veto plane  
		//////////////////////////
		//2 thinner "Veto" layer that contains 2x32 bars (side by side) of 2 x 11 x 70 cm.
		//2x12 bars in the first layer and evenly distributed in center y
		//2x20 bars in the 2nd layer and 10 bars distributed at both end of y

		double VetoBarX=2*cm, VetoBarY=11*cm, VetoBarZ=70*cm;
		G4Box* VetoBar = new G4Box("VetoBar",VetoBarX/2,VetoBarY/2,VetoBarZ/2);

		double VetoVolX=2*cm+0.4*cm, VetoVolY=300*cm, VetoVolZ=70*cm*2;
		G4Box* VetoVol = new G4Box("VetoVol",VetoVolX/2,VetoVolY/2,VetoVolZ/2);

		G4LogicalVolume* logVetoBar = new G4LogicalVolume(VetoBar,
			scintillator, "logVetoBar",  0, 0, 0);
		logVetoBar->SetVisAttributes(PurpleVisAtt);

		logVetoBar->SetSensitiveDetector(SDVeto);

		G4LogicalVolume* logVeto1Vol = new G4LogicalVolume(VetoVol,
			vide, "logVeto1Vol",  0, 0, 0);
		logVeto1Vol->SetVisAttributes(HallVisAtt);
		G4LogicalVolume* logVeto2Vol = new G4LogicalVolume(VetoVol,
			vide, "logVeto2Vol",  0, 0, 0);
		logVeto2Vol->SetVisAttributes(HallVisAtt);

		//place the Veto1 and Veto2 into the container	
		double Veto1PosX=-NDContainerX/2+Veto2NDContainer+VetoVolX/2;	
		double Veto2PosX=Veto1PosX+VetoVolX/2+VetoVolX/2;
		new G4PVPlacement( 0,G4ThreeVector(Veto1PosX,0,0),
			logVeto1Vol,"phyVeto1Vol",logNDContainer,0,0);
		new G4PVPlacement(0,G4ThreeVector(Veto2PosX,0,0),
			logVeto2Vol,"phyVeto2Vol",logNDContainer,0,0);

		for(int i=0;i<12;i++)
		{
			double thePosY=(i-6+0.5)*VetoBarY;
			// Placement of the 12 bars in the nagtive z half
			new G4PVPlacement(0,G4ThreeVector(0,thePosY,-VetoBarZ/2),
				logVetoBar,"phyVetoBar",logVeto1Vol,true,i);
			// Placement of the 12 bars in positive z half
			new G4PVPlacement(0,G4ThreeVector(0,thePosY,VetoBarZ/2),
				logVetoBar,"phyVetoBar",logVeto1Vol,true,23-i);
		}

		for(int i=0;i<10;i++)
		{
			double thePosY=-VetoVolY/2+(i+0.5)*VetoBarY;
			// Placement of the 12 bars in the nagtive z half
			new G4PVPlacement(0,G4ThreeVector(0,thePosY,-VetoBarZ/2),
				logVetoBar,"phyVetoBar",logVeto2Vol,true,i);
			new G4PVPlacement(0,G4ThreeVector(0,-thePosY,-VetoBarZ/2),
				logVetoBar,"phyVetoBar",logVeto2Vol,true,19-i);
			// Placement of the 12 bars in positive z half
			new G4PVPlacement(0,G4ThreeVector(0,thePosY,VetoBarZ/2),
				logVetoBar,"phyVetoBar",logVeto2Vol,true,20+i);
			new G4PVPlacement(0,G4ThreeVector(0,-thePosY,VetoBarZ/2),
				logVetoBar,"phyVetoBar",logVeto2Vol,true,39-i);
		}


		//////////////////////////
		//the first auxilliary scintillator plane  
		//////////////////////////
		//30 bars of 10x10x100 cm

		const int NDSC1BarNum=30;
		double NDSC1BarX = 10*cm, NDSC1BarY = 10*cm, NDSC1BarZ = 100*cm; 
		G4Box* NDSC1Bar = new G4Box("NDSC1Bar",NDSC1BarX/2,NDSC1BarY/2,NDSC1BarZ/2);

		// Creation of the associated logical Volume
		G4LogicalVolume* logNDSC1Bar = new G4LogicalVolume(NDSC1Bar, 
			scintillator, "logNDSC1Bar", 0, 0, 0);
		logNDSC1Bar->SetVisAttributes(DarkBlueVisAtt);

		logNDSC1Bar->SetSensitiveDetector(SDNDSC1);

		double NDSC1VolX = 10*cm+0.2*cm, NDSC1VolY = 300*cm, NDSC1VolZ = 100*cm;
		G4Box* NDSC1Vol = new G4Box("NDSC1Vol",NDSC1VolX/2,NDSC1VolY/2,NDSC1VolZ/2);

		G4LogicalVolume* logNDSC1Vol = new G4LogicalVolume(NDSC1Vol, 
			vide, "logNDSC1Vol", 0, 0, 0);
		logNDSC1Vol->SetVisAttributes(HallVisAtt);

		// Placement of the 30 bars vith replication
		new G4PVReplica("replicaNDSC1Vol",logNDSC1Bar,logNDSC1Vol,
			kYAxis,NDSC1BarNum,NDSC1BarY,-NDSC1VolY/2+0.5*NDSC1BarY);

		double NDSC1PosX=Veto2PosX+VetoVolX/2+NDSC1VolX/2;	
		new G4PVPlacement(0,G4ThreeVector(NDSC1PosX,0,0),
			logNDSC1Vol,"phyNDSC1Vol",logNDContainer,0,0);

		//////////////////////////
		//the 2nd auxilliary scintillator plane  
		//////////////////////////
		//24 bars of 10x12.5x100 cm

		const int NDSC2BarNum=24;
		double NDSC2BarX = 10*cm, NDSC2BarY = 12.5*cm, NDSC2BarZ = 100*cm; 
		G4Box* NDSC2Bar = new G4Box("NDSC2Bar",NDSC2BarX/2,NDSC2BarY/2,NDSC2BarZ/2);

		// Creation of the associated logical Volume
		G4LogicalVolume* logNDSC2Bar = new G4LogicalVolume(NDSC2Bar, 
			scintillator, "logNDSC2Bar", 0, 0, 0);
		logNDSC2Bar->SetVisAttributes(LightYellowVisAtt);

		logNDSC2Bar->SetSensitiveDetector(SDNDSC2);

		double NDSC2VolX = 10*cm+0.4*cm, NDSC2VolY = 300*cm, NDSC2VolZ = 100*cm;
		G4Box* NDSC2Vol = new G4Box("NDSC2Vol",NDSC2VolX/2,NDSC2VolY/2,NDSC2VolZ/2);

		G4LogicalVolume* logNDSC2Vol = new G4LogicalVolume(NDSC2Vol, 
			vide, "logNDSC2Vol", 0, 0, 0);
		logNDSC2Vol->SetVisAttributes(HallVisAtt);

		// Placement of the 24 bars vith replication
		new G4PVReplica("replicaNDSC2Vol",logNDSC2Bar,logNDSC2Vol,
			kYAxis,NDSC2BarNum,NDSC2BarY,-NDSC2VolY/2+0.5*NDSC2BarY);

		double NDSC2PosX=NDSC1PosX+NDSC1VolX/2+NDSC2VolX/2;	
		new G4PVPlacement(0,G4ThreeVector(NDSC2PosX,0,0),
			logNDSC2Vol,"phyNDSC1Vol",logNDContainer,0,0);


		//////////////////////////
		//the 3rd auxilliary scintillator plane  
		//////////////////////////
		//three kinds of bars, 22 bars in total 
		//12 NDSC3 bars of 10x15x100 cm, 8 NDSC2 bars and 2NDSC1 bars
		//6 NDSC3 + 4 NDSC2 + 2 NDSC1 + 4 NDSC2 + 6NDSC3

		double NDSC3BarX = 10*cm, NDSC3BarY = 15*cm, NDSC3BarZ = 100*cm; 
		G4Box* NDSC3Bar = new G4Box("NDSC3Bar",NDSC3BarX/2,NDSC3BarY/2,NDSC3BarZ/2);

		// Creation of the associated logical Volume
		G4LogicalVolume* logNDSC3Bar = new G4LogicalVolume(NDSC3Bar, 
			scintillator, "logNDSC3Bar", 0, 0, 0);
		logNDSC3Bar->SetVisAttributes(RedVisAtt);

		logNDSC3Bar->SetSensitiveDetector(SDNDSC3);

		double NDSC3VolX = 10*cm+0.4*cm, NDSC3VolY = 300*cm, NDSC3VolZ = 100*cm;
		G4Box* NDSC3Vol = new G4Box("NDSC3Vol",NDSC3VolX/2,NDSC3VolY/2,NDSC3VolZ/2);

		G4LogicalVolume* logNDSC3Vol = new G4LogicalVolume(NDSC3Vol, 
			vide, "logNDSC3Vol", 0, 0, 0);
		logNDSC3Vol->SetVisAttributes(HallVisAtt);


		double NDSC3PosX=NDSC2PosX+NDSC2VolX/2+NDSC3VolX/2;	

		new G4PVPlacement(0,G4ThreeVector(NDSC3PosX,0,0),
			logNDSC3Vol,"phyNDSC3Vol",logNDContainer,0,0);

		for(int i=0;i<6;i++)
		{
			double thePosY=-NDSC3VolY/2+(i+0.5)*NDSC3BarY;
			// Placement of the 6 NDSC3 bars in the bottom half
			new G4PVPlacement(0,G4ThreeVector(0,thePosY,0),
				logNDSC3Bar,"phyNDSC3Bar",logNDSC3Vol,true,i);
			// Placement of the 6 NDSC3 bars in the top half
			new G4PVPlacement(0,G4ThreeVector(0,-thePosY,0),
				logNDSC3Bar,"phyNDSC3Bar",logNDSC3Vol,true,21-i);
		}
		for(int i=6;i<10;i++)
		{
			double thePosY=-NDSC3VolY/2+6*NDSC3BarY+(i-6+0.5)*NDSC2BarY;
			// Placement of the 4 NDSC2 bars in the bottom half
			new G4PVPlacement(0,G4ThreeVector(0,thePosY,0),
				logNDSC2Bar,"phyNDSC3Bar",logNDSC3Vol,true,i);
			// Placement of the 4 NDSC2 bars in the top half
			new G4PVPlacement(0,G4ThreeVector(0,-thePosY,0),
				logNDSC2Bar,"phyNDSC3Bar",logNDSC3Vol,true,15-i);
		}

		// Placement of the 2 NDSC1 bars in the center
		new G4PVPlacement(0,G4ThreeVector(0,-NDSC1BarY/2,0),
			logNDSC1Bar,"phyNDSC3Bar",logNDSC3Vol,true,10);
		new G4PVPlacement(0,G4ThreeVector(0,NDSC1BarY/2,0),
			logNDSC1Bar,"phyNDSC3Bar",logNDSC3Vol,true,11);


		//////////////////////////
		//the 4th auxilliary scintillator plane  
		//////////////////////////
		//12 bars of 10x25x100 cm

		const int NDSC4BarNum=12;
		double NDSC4BarX = 10*cm, NDSC4BarY = 25*cm, NDSC4BarZ = 100*cm; 
		G4Box* NDSC4Bar = new G4Box("NDSC4Bar",NDSC4BarX/2,NDSC4BarY/2,NDSC4BarZ/2);

		// Creation of the associated logical Volume
		G4LogicalVolume* logNDSC4Bar = new G4LogicalVolume(NDSC4Bar, 
			scintillator, "logNDSC4Bar", 0, 0, 0);
		logNDSC4Bar->SetVisAttributes(LightGreenVisAtt);

		logNDSC4Bar->SetSensitiveDetector(SDNDSC4);

		double NDSC4VolX = 10*cm+0.4*cm, NDSC4VolY = 300*cm, NDSC4VolZ = 100*cm;
		G4Box* NDSC4Vol = new G4Box("NDSC4Vol",NDSC4VolX/2,NDSC4VolY/2,NDSC4VolZ/2);

		G4LogicalVolume* logNDSC4Vol = new G4LogicalVolume(NDSC4Vol, 
			vide, "logNDSC4Vol", 0, 0, 0);
		logNDSC4Vol->SetVisAttributes(HallVisAtt);

		// Placement of the 12 bars vith replication
		new G4PVReplica("replicaNDSC4Vol",logNDSC4Bar,logNDSC4Vol,
			kYAxis,NDSC4BarNum,NDSC4BarY,-NDSC4VolY/2+0.5*NDSC4BarY);

		double NDSC4PosX=NDSC3PosX+NDSC3VolX/2+NDSC4VolX/2;	
		new G4PVPlacement(0,G4ThreeVector(NDSC4PosX,0,0),
			logNDSC4Vol,"phyNDSC4Vol",logNDContainer,0,0);

#endif
	} //end of HAND
	//////////////////////

	return phyBBContainer;
}


/////////////////////////////////////////////////////////////////////

