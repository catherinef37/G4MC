/*///////////////////////////////////////////////////////////////////////////////////////////////
By Jixie Zhang @ 20121112
This class is used to build the geometry for CREX. In order to optimized the design for the
target, I put all variables into the configuration file 'Detector_CREX.ini'.   
Note that Detector.ini, is also needed, but Detector_CREX.ini will overwrite any variable 
in Detector.ini, it will not take effects.
/*///////////////////////////////////////////////////////////////////////////////////////////////
#include "CREXDetectorConstruction.hh"

#include "G4FieldManager.hh"
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

extern UsageManager* gConfig;
////////////////////////////////////////////////////////////////////////////////////////////////
// Declaration of defaults constructeurs & destructeurs
CREXDetectorConstruction::CREXDetectorConstruction(G4LogicalVolume *mother) : 
mMotherLogVol(mother) 
{
	GetConfig();	
	mMaterialManager=HRSMaterial::GetHRSMaterialManager();
	ConstructMaterial();
	G4cout<<"Contrstruct CREX geometry ... done! "<<G4endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////
CREXDetectorConstruction::~CREXDetectorConstruction()
{
	G4cout<<"Delete CREX geometry ... done! "<<G4endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////
void CREXDetectorConstruction::GetConfig()
{
	gConfig->ReadFile("Detector_CREX.ini");
	//////////////////////////////////////////////////////////////////
	//the following is just for crex

	gConfig->GetParameter("ScatChamberRin",mScatChamberRin);
	mScatChamberRin*=mm;
	gConfig->GetParameter("ScatChamberRout",mScatChamberRout);
	mScatChamberRout*=mm;
	gConfig->GetParameter("ScatChamberL",mScatChamberL);
	mScatChamberL*=mm;

	gConfig->GetParameter("ScatChamberEntranceWindowThick",mScatChamberEntranceWindowThick);
	mScatChamberEntranceWindowThick*=mm;
	gConfig->GetParameter("ScatChamberExitWindowThick",mScatChamberExitWindowThick);
	mScatChamberExitWindowThick*=mm;

	
	gConfig->GetParameter("SetupCREXTarget",mSetupCREXTarget);
	gConfig->GetParameter("TargetType",mTargetType);
	gConfig->GetParameter("TargetW",mTargetW); mTargetW*=mm;
	gConfig->GetParameter("TargetH",mTargetH); mTargetH*=mm;
	gConfig->GetParameter("TargetL",mTargetL); mTargetL*=mm;

	gConfig->GetParameter("SetupStdScatChamber",mSetupStdScatChamber);

	gConfig->GetParameter("UpBlockRin",mUpBlockRin); 
	mUpBlockRin*=mm;
	gConfig->GetParameter("UpBlockThick",mUpBlockThick); 
	mUpBlockThick*=mm;
	gConfig->GetParameter("DownBlockRin",mDownBlockRin); 
	mDownBlockRin*=mm;
	gConfig->GetParameter("DownBlockThickAt0cm",mDownBlockThickAt0cm); 
	mDownBlockThickAt0cm*=mm;
	gConfig->GetParameter("DownBlockThickAt5cm",mDownBlockThickAt5cm); 
	mDownBlockThickAt5cm*=mm;
	gConfig->GetParameter("UpBlockLength",mUpBlockLength); 
	mUpBlockLength*=mm;
	gConfig->GetParameter("DownBlockLength",mDownBlockLength); 
	mDownBlockLength*=mm;

	gConfig->GetParameter("UpCapThick",mUpCapThick); 
	mUpCapThick*=mm;
	gConfig->GetParameter("DownCapThick",mDownCapThick); 
	mDownCapThick*=mm;

	gConfig->GetParameter("VCRin",mVCRin); 
	mVCRin*=mm;
	gConfig->GetParameter("VCThick",mVCThick); 
	mVCThick*=mm;
	gConfig->GetParameter("VCUpCapThick",mVCUpCapThick); 
	mVCUpCapThick*=mm;
	gConfig->GetParameter("VCDownCapThick",mVCDownCapThick); 
	mVCDownCapThick*=mm;
	gConfig->GetParameter("VCLength",mVCLength); 
	mVCLength*=mm;

	gConfig->GetParameter("VC2SCZOffset",mVC2SCZOffset); 
	mVC2SCZOffset*=mm;

	////////////////////////////////////////////////////////////////
	gConfig->GetParameter("Pivot2LSieveFace",mPivot2LSieveFace);
	mPivot2LSieveFace*=mm;
	gConfig->GetParameter("Pivot2RSieveFace",mPivot2RSieveFace);
	mPivot2RSieveFace*=mm;

	gConfig->GetParameter("Pivot2LHRSVBFace",mPivot2LHRSVBFace);
	mPivot2LHRSVBFace*=mm;
	gConfig->GetParameter("Pivot2RHRSVBFace",mPivot2RHRSVBFace);
	mPivot2RHRSVBFace*=mm;

	gConfig->GetParameter("HRSVBWidth",mHRSVBWidth);
	mHRSVBWidth*=mm;
	gConfig->GetParameter("HRSVBHeight",mHRSVBHeight);
	mHRSVBHeight*=mm;
	gConfig->GetParameter("HRSVBThick",mHRSVBThick);
	mHRSVBThick*=mm;

	////////////////////////////////////////////////////////////
	//the following is for global (HRS)
	
	gConfig->GetParameter("SetupLHRS",mSetupLHRS);
	gConfig->GetParameter("SetupRHRS",mSetupRHRS);
	gConfig->GetParameter("SetupLSieveSlit",mSetupLSieveSlit);
	gConfig->GetParameter("SetupRSieveSlit",mSetupRSieveSlit);

	gConfig->GetParameter("LHRSAngle",mLHRSAngle);
	mLHRSAngle*=deg;
	gConfig->GetParameter("LSeptumAngle",mLSeptumAngle);
	mLSeptumAngle*=deg;
	gConfig->GetParameter("RHRSAngle",mRHRSAngle);
	mRHRSAngle*=deg;
	gConfig->GetParameter("RSeptumAngle",mRSeptumAngle);
	mRSeptumAngle*=deg;
	

	gConfig->GetParameter("PivotXOffset",mPivotXOffset);
	mPivotXOffset*=mm;
	gConfig->GetParameter("PivotYOffset",mPivotYOffset);
	mPivotYOffset*=mm;
	gConfig->GetParameter("PivotZOffset",mPivotZOffset);
	mPivotZOffset*=mm;

	gConfig->GetParameter("ScatChamberXOffset",mScatChamberXOffset);
	mScatChamberXOffset*=mm;
	gConfig->GetParameter("ScatChamberYOffset",mScatChamberYOffset);
	mScatChamberYOffset*=mm;
	gConfig->GetParameter("ScatChamberZOffset",mScatChamberZOffset);
	mScatChamberZOffset*=mm;
	gConfig->GetParameter("TargetXOffset",mTargetXOffset);
	mTargetXOffset*=mm;
	gConfig->GetParameter("TargetYOffset",mTargetYOffset);
	mTargetYOffset*=mm;
	gConfig->GetParameter("TargetZOffset",mTargetZOffset);
	mTargetZOffset*=mm;
}

////////////////////////////////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* CREXDetectorConstruction::Construct()
{
	G4VPhysicalVolume *theCREXPhys=0;
	
	/////////////////////////
	//Hall A scatter chamber 
	/////////////////////////
	if(mSetupStdScatChamber)  ConstructTargetChamber(mMotherLogVol);

	/////////////////////////
	//g2p target container and target 
	/////////////////////////	

	if(mSetupCREXTarget)  ConstructTarget(mMotherLogVol);


	//DVCS Solenoid
	//ConstructDVCSSolenoid(mMotherLogVol);

	/////////////////////////
	// Sieve slit, septum window and HRS VirtualBoundary
	/////////////////////////
	if(mSetupLHRS || mSetupRHRS) ConstructSeptumNSieve(mMotherLogVol);


	return theCREXPhys;
}

////////////////////////////////////////////////////////////////////////////////////////////////
void CREXDetectorConstruction::ConstructMaterial()
{
	calcium=mMaterialManager->calcium;
	air=mMaterialManager->air;
	vacuum=mMaterialManager->vacuum;
	stainlesssteel=mMaterialManager->stainlesssteel;
	aluminum=mMaterialManager->aluminum;
	heliumGas=mMaterialManager->heliumGas;
	
	if (mTargetType==0) theTargetMaterial=mMaterialManager->vacuum;
	else if (mTargetType==1) theTargetMaterial=mMaterialManager->NH3He;
	else if (mTargetType==2) theTargetMaterial=mMaterialManager->CH2;
	else if (mTargetType==3) theTargetMaterial=mMaterialManager->carbon;
	else if (mTargetType==4) theTargetMaterial=mMaterialManager->tantalum;
	else if (mTargetType==5) theTargetMaterial=mMaterialManager->liquidH2;
	else if (mTargetType==6) theTargetMaterial=mMaterialManager->liquidD2;
	else if (mTargetType==7) theTargetMaterial=mMaterialManager->liquidHe3;
	else if (mTargetType==8) theTargetMaterial=mMaterialManager->liquidHe;
	else if (mTargetType==9) theTargetMaterial=mMaterialManager->aluminum;
	else if (mTargetType==10) theTargetMaterial=mMaterialManager->copper;
	else if (mTargetType==11) theTargetMaterial=mMaterialManager->lead;
	else if (mTargetType==12) theTargetMaterial=mMaterialManager->tungsten;
	else if (mTargetType==13) theTargetMaterial=mMaterialManager->stainlesssteel;
	else if (mTargetType==14) theTargetMaterial=mMaterialManager->kapton;
	else if (mTargetType==15) theTargetMaterial=mMaterialManager->air;       //room temperature
	else if (mTargetType==16) theTargetMaterial=mMaterialManager->heliumGas; //room temperature
	else if (mTargetType==17) theTargetMaterial=mMaterialManager->calcium;  
	else theTargetMaterial=mMaterialManager->calcium;

}

////////////////////////////////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* CREXDetectorConstruction::ConstructTargetChamber(G4LogicalVolume *pMotherLogVol)
{
	const double inch=2.54*cm;

	G4RotationMatrix *pRotScatInHall=new G4RotationMatrix();
	pRotScatInHall->rotateX(90.*deg);

	double startphi=0.*deg, deltaphi=360.*deg;
	/////////////////////////////////////////////////////////
	//scattering chamber container
	/////////////////////////////////////////////////////////
	//This is just a container to enclose the taraget chamber, 10 mm larger than the 
	//scattering chamber itself.
	//With this container, all stuff inside do not need a rotation
	
	//The scattering chamber containner is made of helium gas, 

	double pScatChamberContainerRin=mScatChamberRin-1*cm;      
	double pScatChamberContainerRout=mScatChamberRout+1*cm;
	double pScatChamberContainerL=mScatChamberL+(3.50+17.0+1.25)*inch*2+10.0*mm;
	G4VSolid* scatChamberContainerExtendedSolid = new G4Tubs("scatChamberContainerExtendedTubs",
		pScatChamberContainerRin,pScatChamberContainerRout,
		pScatChamberContainerL/2.0,0.,360.*deg);
	G4VSolid* scatChamberContainerExtraSolid = new G4Tubs("scatChamberContainerExtraTubs",
		0,pScatChamberContainerRout+1*mm,17.25*inch/2.0,0.,360.*deg);
	G4SubtractionSolid* scatChamberContainerSolid=new G4SubtractionSolid("scatChamberContainerSolid",
		scatChamberContainerExtendedSolid,scatChamberContainerExtraSolid,
		0,G4ThreeVector(0,0,-mScatChamberL/2-17.25*inch/2.0));

	G4LogicalVolume* scatChamberContainerLogical = new G4LogicalVolume(scatChamberContainerSolid,
		heliumGas,"scatChamberContainerLogical",0,0,0);
	scatChamberContainerLogical->SetVisAttributes(HallVisAtt); 

	G4VPhysicalVolume* scatChamberContainerPhys=new G4PVPlacement(pRotScatInHall,
		G4ThreeVector(mScatChamberXOffset,mScatChamberYOffset,mScatChamberZOffset),
		scatChamberContainerLogical,"scatChamberContainerPhys",pMotherLogVol,0,0);

	//////////////////////////
	// build the target chamber.
	//////////////////////////
	//Build the simplified scatter chamber,it contains 2 windows of rectangles 
	//The following already defined in the config file
	//double mScatChamberRin=17.875*inch,mScatChamberRout=18.875*inch,mScatChamberL=27.25*inch;

	G4VSolid* scatChamberWholeSolid=0;
	//If mSetupScatChamber==1, setup the body only, 
	//If mSetupScatChamber==2, setup the body plus top flange and bottom flange, this
	//will make the program slower
	if(mSetupStdScatChamber==1)
	{
		scatChamberWholeSolid = new G4Tubs("scatChamberWholeTubs",
			mScatChamberRin,mScatChamberRout,mScatChamberL/2.0,0.,360.*deg);
	}
	else if(mSetupStdScatChamber>=2)
	{
		startphi=0.*deg; deltaphi=360.*deg;
		const int kNPlane_SC=11;
		double rInner_SC[] = {0,0,mScatChamberRin,
			mScatChamberRin,mScatChamberRin,mScatChamberRin,
			mScatChamberRin,mScatChamberRin,mScatChamberRin,
			0,0};
		double rOuter_SC[] = {
			mScatChamberRout+1.0*inch,mScatChamberRout+1.0*inch,mScatChamberRout+1.0*inch,
			mScatChamberRout,mScatChamberRout,mScatChamberRout+1.0*inch,
			mScatChamberRout+1.0*inch,mScatChamberRout,mScatChamberRout,
			mScatChamberRout+1*inch,mScatChamberRout+1*inch
		};
		double zPlane_SC[] = {
			-mScatChamberL/2-4.50*inch,-mScatChamberL/2-3.25*inch,-mScatChamberL/2-1.0*inch,
			-mScatChamberL/2,mScatChamberL/2+0.25*inch,mScatChamberL/2+1.25*inch,
			mScatChamberL/2+3.50*inch,mScatChamberL/2+3.50*inch,mScatChamberL/2+20.5*inch,
			mScatChamberL/2+20.5*inch,mScatChamberL/2+21.75*inch
		};

		G4Polycone* SCWholeSolid = new G4Polycone("SCPolycone",startphi,deltaphi,
			kNPlane_SC,zPlane_SC,rInner_SC,rOuter_SC);

		scatChamberWholeSolid = SCWholeSolid;
	}


	//these are the subtraction part, not the scatter chamber itself
	double pSCWindowRin=mScatChamberRin-1*mm;
	double pSCWindowRout=mScatChamberRout+1*mm;
	double pSCEntranceWindowH=6.44*inch;
	double pSCDownCapH=15.0*inch;
	
	//rectangle EntranceWindow covering 80 to 100 degrees
	startphi=80*deg;deltaphi=20*deg;
	G4VSolid* SCEntranceWindowSolid = new G4Tubs("SCEntranceWindowTubs",
		pSCWindowRin,pSCWindowRout,pSCEntranceWindowH/2.0,startphi,deltaphi);

	//rectangle DownCap covering -225 to 45 degrees
	startphi=-225*deg;deltaphi=270*deg;
	G4VSolid* SCDownCapSolid = new G4Tubs("SCDownCapTubs",
		pSCWindowRin,pSCWindowRout,pSCDownCapH/2.0,startphi,deltaphi);


	// subtract the Entrance window 
	G4SubtractionSolid* SCSubtractEntranceSolid=new G4SubtractionSolid(
		"SCSubtractEntrance",scatChamberWholeSolid,SCEntranceWindowSolid);

	// subtract the Exit window 
	G4SubtractionSolid* SCSubtractEntranceNExitSolid=new G4SubtractionSolid(
		"SCSubtractEntranceNExit",SCSubtractEntranceSolid,SCDownCapSolid);

	//setup the scatter chamber
	G4LogicalVolume* scatChamberLogical = new G4LogicalVolume(
		SCSubtractEntranceNExitSolid,aluminum,"scatChamberLogical",0,0,0);
	scatChamberLogical->SetVisAttributes(WhiteVisAtt); 

	new G4PVPlacement(0,G4ThreeVector(),
		scatChamberLogical,"scatChamberPhys",scatChamberContainerLogical,0,0);

	
	/////////////////////////
	// target chamber window covers 
	/////////////////////////
	//Covers for EntranceWindow 

	//EntranceWindowCover
	double pSCEntranceWindowCoverH=pSCEntranceWindowH+0.8*inch;
	double pSCEntranceWindowCoverRin=mScatChamberRout;
	double pSCEntranceWindowCoverRout=pSCEntranceWindowCoverRin+mScatChamberEntranceWindowThick;

	startphi=78*deg;deltaphi=24*deg;
	G4VSolid* SCEntranceWindowCoverSolid = new G4Tubs("SCEntranceWindowCoverTubs",
		pSCEntranceWindowCoverRin,pSCEntranceWindowCoverRout,
		pSCEntranceWindowCoverH/2.0,startphi,deltaphi);
	G4LogicalVolume* SCEntranceWindowCoverLogical = new G4LogicalVolume(
		SCEntranceWindowCoverSolid,aluminum,"SCEntranceWindowCoverLogical",0,0,0);
	SCEntranceWindowCoverLogical->SetVisAttributes(LightYellowVisAtt); 

	new G4PVPlacement(0,G4ThreeVector(),SCEntranceWindowCoverLogical,
		"SCEntranceWindowCoverPhys",scatChamberContainerLogical,false,0);

	//DownCapCover
	double pSCDownCapCoverH=pSCDownCapH+0.8*inch;
	double pSCDownCapCoverRin=mScatChamberRout;
	double pSCDownCapCoverRout=mScatChamberRout+mScatChamberExitWindowThick;

	startphi=-227*deg;deltaphi=274*deg;
	G4VSolid* SCDownCapCoverSolid = new G4Tubs("SCDownCapCoverTubs",
		pSCDownCapCoverRin,pSCDownCapCoverRout,pSCDownCapCoverH/2.0,startphi,deltaphi);
	G4LogicalVolume* SCDownCapCoverLogical = new G4LogicalVolume(
		SCDownCapCoverSolid,aluminum,"SCDownCapCoverLogical",0,0,0);
	SCDownCapCoverLogical->SetVisAttributes(LightYellowVisAtt); 

	new G4PVPlacement(0,G4ThreeVector(),SCDownCapCoverLogical,
		"SCDownCapCoverPhys",scatChamberContainerLogical,false,0);


	return scatChamberContainerPhys;
}


////////////////////////////////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* CREXDetectorConstruction::ConstructTarget(G4LogicalVolume *pMotherLogVol)
{
	G4SDManager* SDMan = G4SDManager::GetSDMpointer();
	
	G4VSensitiveDetector* SDUpBlock = new HRSStdSD("CREXUpBlock");
	SDMan->AddNewDetector(SDUpBlock);
	G4VSensitiveDetector* SDDownBlock = new HRSStdSD("CREXDownBlock");
	SDMan->AddNewDetector(SDDownBlock);

	G4VSensitiveDetector* SDTarget  = new HRSStdSD("CREXTarget");
	SDMan->AddNewDetector(SDTarget);

	double startphi=0, deltaphi=360*deg;
	//############################################
	
	G4RotationMatrix *pRotX90deg=new G4RotationMatrix();
	pRotX90deg->rotateX(90.*deg);
	G4RotationMatrix *pRotX270deg=new G4RotationMatrix();
	pRotX270deg->rotateX(270.*deg);

	/////////////////////////////////////////////////////////
	//inner of the scattering chamber container
	/////////////////////////////////////////////////////////
	//A container to enclose everything inside the taraget chamber, 
	//10 mm larger than the scattering chamber itself.
	//With this container, all stuff inside do not need a rotation
	

	double targetContainerRin=0;      
	double targetContainerRout=mScatChamberRin-1*cm;
	double targetContainerL=mScatChamberL-1*cm;
	G4VSolid* targetContainerSolid = new G4Tubs("targetContainerTubs",
		targetContainerRin,targetContainerRout,
		targetContainerL/2.0,0.,360.*deg);

	G4LogicalVolume* targetContainerLogical = new G4LogicalVolume(targetContainerSolid,
		vacuum,"targetContainerLogical",0,0,0);
	//By Jixie: Add this step limit 
	double pScatChamberStepLimit=10;
	gConfig->GetArgument("ScatChamberStepLimit",pScatChamberStepLimit);
	pScatChamberStepLimit*=mm;
	G4UserLimits* uSCStepLimits = new G4UserLimits(pScatChamberStepLimit);
	targetContainerLogical->SetUserLimits(uSCStepLimits);
	targetContainerLogical->SetVisAttributes(HallVisAtt); 

	G4VPhysicalVolume* targetContainerPhys=new G4PVPlacement(pRotX90deg,
		G4ThreeVector(mScatChamberXOffset,mScatChamberYOffset,mScatChamberZOffset),
		targetContainerLogical,"targetContainerPhys",pMotherLogVol,0,0);

	/////////////////////////////////////////////////////////////////
	//since target containner has been rotated around X aixs by 90 degree (clockwisely), 
	//new z' = old y and new y' = old -z
	//In order to get original lab coordinate layout, we have to rotate by X an extra 270 degree.

	//////////////////////////
	//the target
	//////////////////////////

	double pTg2SC_X=mTargetXOffset-mScatChamberXOffset;
	double pTg2SC_Y=mTargetYOffset-mScatChamberYOffset;
	double pTg2SC_Z=mTargetZOffset-mScatChamberZOffset;

	G4VSolid* targetSolid=new G4Box("targetBox",mTargetW/2.0,
		mTargetH/2.0,mTargetL/2.0);

	G4LogicalVolume* targetLogical = new G4LogicalVolume(
		targetSolid,calcium,"targetLogical",0,0,0);
	//By Jixie: Add this step limit in the target
	double pTargetStepLimit=10;
	gConfig->GetArgument("TargetStepLimit",pTargetStepLimit);
	pTargetStepLimit*=mm;
	G4UserLimits* uTgtepLimits = new G4UserLimits(pTargetStepLimit);
	targetLogical->SetUserLimits(uTgtepLimits);
	targetLogical->SetVisAttributes(LightPurpleVisAtt); 

	targetLogical->SetSensitiveDetector(SDTarget);

	new G4PVPlacement(pRotX270deg,G4ThreeVector(pTg2SC_X,-pTg2SC_Z,pTg2SC_Y),
		targetLogical,"targetPhys",targetContainerLogical,0,0);


	//////////////////////////
	//the vacuum Chamber tube
	//////////////////////////

	//the tube without caps
	double pVCLengthWithoutCaps=mVCLength-mVCUpCapThick-mVCDownCapThick;
	double pVCRout=mVCRin+mVCThick;

	G4VSolid* VCSolid = new G4Tubs("VCTubs",
		mVCRin,pVCRout,pVCLengthWithoutCaps/2.0,0,360*deg);

	G4LogicalVolume* VCLogical = new G4LogicalVolume(
		VCSolid,stainlesssteel,"VCLogical",0,0,0);
	VCLogical->SetVisAttributes(SteelVisAtt); 

	new G4PVPlacement(pRotX270deg,G4ThreeVector(0,-mVC2SCZOffset,0),
		VCLogical,"VCPhys",targetContainerLogical,0,0);

	
	//////////////////////////
	//the vacuum Chamber up cap
	//////////////////////////
	G4VSolid* VCUpCapSolid = new G4Tubs("VCUpCapTubs",
		mUpBlockRin,pVCRout,mVCUpCapThick/2.0,0,360*deg);

	G4LogicalVolume* VCUpCapLogical = new G4LogicalVolume(
		VCUpCapSolid,stainlesssteel,"VCUpCapLogical",0,0,0);
	VCUpCapLogical->SetVisAttributes(SteelVisAtt); 

	double pVCUpCap2SC_Z=mVC2SCZOffset-(mVCLength/2-mVCUpCapThick/2);
	new G4PVPlacement(pRotX270deg,G4ThreeVector(0,-pVCUpCap2SC_Z,0),
		VCUpCapLogical,"VCUpCapPhys",targetContainerLogical,0,0);


	//////////////////////////
	//the vacuum Chamber down cap
	//////////////////////////
	G4VSolid* VCDownCapSolid = new G4Tubs("VCDownCapTubs",
		mDownBlockRin,pVCRout,mVCDownCapThick/2.0,0,360*deg);

	G4LogicalVolume* VCDownCapLogical = new G4LogicalVolume(
		VCDownCapSolid,stainlesssteel,"VCDownCapLogical",0,0,0);
	VCDownCapLogical->SetVisAttributes(SteelVisAtt); 

	double pVCDownCap2SC_Z=mVC2SCZOffset+(mVCLength/2-mVCDownCapThick/2);
	new G4PVPlacement(pRotX270deg,G4ThreeVector(0,-pVCDownCap2SC_Z,0),
		VCDownCapLogical,"VCDownCapPhys",targetContainerLogical,0,0);


	//////////////////////////
	//the upstream blocker
	//////////////////////////
	G4VSolid* UpBlockSolid = new G4Tubs("UpBlockTubs",
		mUpBlockRin,mUpBlockRin+mUpBlockThick,mUpBlockLength/2.0,0,360*deg);

	G4LogicalVolume* UpBlockLogical = new G4LogicalVolume(
		UpBlockSolid,aluminum,"UpBlockLogical",0,0,0);
	UpBlockLogical->SetVisAttributes(WhiteVisAtt); 
	UpBlockLogical->SetSensitiveDetector(SDUpBlock);

	double pUpBlock2SC_Z=mVC2SCZOffset-(mVCLength/2+mUpBlockLength/2);
	new G4PVPlacement(pRotX270deg,G4ThreeVector(0,-pUpBlock2SC_Z,0),
		UpBlockLogical,"UpBlockPhys",targetContainerLogical,0,0);


	//////////////////////////
	//the upstream end cap
	//////////////////////////
	G4VSolid* UpCapSolid = new G4Tubs("UpCapTubs",
		0,mUpBlockRin,mUpCapThick/2.0,0,360*deg);

	G4LogicalVolume* UpCapLogical = new G4LogicalVolume(
		UpCapSolid,aluminum,"UpCapLogical",0,0,0);
	UpCapLogical->SetVisAttributes(LightYellowVisAtt); 

	double pUpCap2SC_Z=mVC2SCZOffset-(mVCLength/2+mUpBlockLength-mUpCapThick/2);
	new G4PVPlacement(pRotX270deg,G4ThreeVector(0,-pUpCap2SC_Z,0),
		UpCapLogical,"UpCapPhys",targetContainerLogical,0,0);

	
	//////////////////////////
	//the downstream end cap
	//////////////////////////
	G4VSolid* DownCapSolid = new G4Tubs("DownCapTubs",
		0,mDownBlockRin,mDownCapThick/2.0,0,360*deg);

	G4LogicalVolume* DownCapLogical = new G4LogicalVolume(
		DownCapSolid,aluminum,"DownCapLogical",0,0,0);
	DownCapLogical->SetVisAttributes(LightYellowVisAtt); 

	double pDownCap2SC_Z=mVC2SCZOffset+mVCLength/2-mDownCapThick/2;
	new G4PVPlacement(pRotX270deg,G4ThreeVector(0,-pDownCap2SC_Z,0),
		DownCapLogical,"DownCapPhys",targetContainerLogical,0,0);

	//////////////////////////
	//the downstream blocker
	//////////////////////////
	//the downstream block can be taped, starting at a thickness of 2 mm nearest the target 
	//and increasing to 4 mm thick 5 cm further downstream.

	startphi=0.*deg; deltaphi=360.*deg;
	const int kNPlane_DownBlock=3;
	double pDownBlockRout=mDownBlockRin+mDownBlockThickAt5cm;
	double rInner_DownBlock[] = {mDownBlockRin,mDownBlockRin,mDownBlockRin};
	double rOuter_DownBlock[] = {mDownBlockRin+mDownBlockThickAt0cm,pDownBlockRout,pDownBlockRout};
	double zPlane_DownBlock[] = {0,5*cm,mDownBlockLength};

	G4Polycone* DownBlockSolid = new G4Polycone("DownBlockPolycone",startphi,deltaphi,
		kNPlane_DownBlock,zPlane_DownBlock,rInner_DownBlock,rOuter_DownBlock);

	G4LogicalVolume* DownBlockLogical = new G4LogicalVolume(
		DownBlockSolid,aluminum,"DownBlockLogical",0,0,0);
	DownBlockLogical->SetVisAttributes(WhiteVisAtt);
	DownBlockLogical->SetSensitiveDetector(SDDownBlock); 

	double pDownBlockFace2VC=mVC2SCZOffset+mVCLength/2;
	new G4PVPlacement(pRotX270deg,G4ThreeVector(0,-pDownBlockFace2VC,0),
		DownBlockLogical,"DownBlockPhys",targetContainerLogical,0,0);


	////////////////////////////////////////////////
	return targetContainerPhys;

}

/////////////////////////////////////////////////////////////////////
//this is for g2p, need to update this routine later
G4VPhysicalVolume* CREXDetectorConstruction::ConstructSeptumNSieve(G4LogicalVolume* motherLogical)
{
	const double inch=2.54*cm;
	double startphi,deltaphi;
	G4String SDname;

	// All G4VSensitiveDetector will be managed (deleted) by SDManager
	//therefore no need to delete it at the end of this subroutine
	G4VSensitiveDetector* sieveSlitSD=new HRSStdSD(SDname="sieveSlit");
	G4VSensitiveDetector* septumWindowSD=new HRSStdSD(SDname="septumWindow");

	G4SDManager* SDman = G4SDManager::GetSDMpointer();

	int pSeptum_UseUniformB=0;
	gConfig->GetParameter("Septum_UseUniformB",pSeptum_UseUniformB); 
	double pSeptumCurrentRatioL=1.0, pSeptumCurrentRatioR=1.0;  
	gConfig->GetParameter("Septum_CurrentRatioL",pSeptumCurrentRatioL);
	gConfig->GetParameter("Septum_CurrentRatioR",pSeptumCurrentRatioR);
	int pUseSeptumField=(pSeptum_UseUniformB==0 && 
		(fabs(pSeptumCurrentRatioL)>1.0E-08 || fabs(pSeptumCurrentRatioR)>1.0E-08) )?1:0;
	//set up the septum only if there is an angle difference	
	bool mSetupSeptumBlock=((mLHRSAngle-mLSeptumAngle)/deg>0.5)?true:false;

	/////////////////////////////////////////////////
	//From Hall A NIM Paper, the standard sieve slit
	//Each spectrometer is equipped with a set of collimators, positioned 1:109 +/- 0.005 
	//and 1:101+/-0.005 m from the target on the left and right spectrometers, respectively. 
	//There is a large collimator, made of 80 mm thick tungsten, with a 121.8 mm vertical 
	//opening and a 62:9 mm horizontal opening at the entrance face. The opening in this 
	//collimator expands to 129.7 by 66.8 mm at the exit face. A second smaller collimator, 
	//made of the same material, is 50.0 by 21.3 mm at the entrance face and 53.2 by 22:6 mm 
	//at the exit face. The third collimator is the sieve slit, which is used to study the 
	//optical properties of the spectro- meters. The sieve is a 5 mm thick stainless steel 
	//sheet with a pattern of 49 holes 7 x 7, spaced 25 mm apart vertically and 12:5 mm apart 
	//horizontally. Two of the holes, one in the center and one displaced two rows vertically 
	//and one horizontally, are 4 mm in diameter. The rest of the holes are 2 mm in diameter. 
	//The sieve slits are positioned 75 mm further from the target than the other collimators.

	//double mLHRSAngle=12.5*deg,mRHRSAngle=-12.5*deg;
	//double mLSeptumAngle=5.69*deg,mRSeptumAngle=-5.69*deg;
	G4RotationMatrix *pRotLHRS=new G4RotationMatrix();
	pRotLHRS->rotateY(-mLHRSAngle); 
	G4RotationMatrix *pRotRHRS=new G4RotationMatrix();
	pRotRHRS->rotateY(-mRHRSAngle); 
	G4RotationMatrix *pRotLSeptum=new G4RotationMatrix();
	pRotLSeptum->rotateY(-mLSeptumAngle); 
	G4RotationMatrix *pRotRSeptum=new G4RotationMatrix();
	pRotRSeptum->rotateY(-mRSeptumAngle); 

	G4RotationMatrix *pRotX90deg=new G4RotationMatrix();
	pRotX90deg->rotateX(-90*deg); 
	/////////////////////////////////////////////////
	G4VPhysicalVolume* theG2PSeptumPhys=0;

	///////////////////////////////////////
	//Sieve Slit for HRS-Angle=5.65
	///////////////////////////////////////

	////////////////the following is for 6 degree sieve//////////////

	//these 2 will be read from ini file, it is 80.0cm for 6 deg and 120 cm for 12.5 deg
	//double mPivot2LSieveFace=79.96*cm;
	//double mPivot2RSieveFace=79.96*cm;

	double pSieveSlitX=2.205*inch; //33.13*mm
	double pSieveSlitY=5.134*inch; //130.40*mm
	double pSieveSlitZ=0.2*inch;

	double pSieveSlitHoleR=0.6985*mm;           //radius of small hole 0.055/2 inch
	double pSieveSlitLargeHoleR=1.3462*mm;      //radius of large hole 0.106/2 inch
	double pSieveSlitHoldL=pSieveSlitZ+0.1*mm;  //need to make it longer to avoid round off in the subtraction

	//the big center hole horizontal and vertical offset, From servey
	double pSieveSlitLargeHoleH=0*mm;    //positive means shift away from the beam line
	double pSieveSlitLargeHoleV=0*mm;

	//please note that the following constants are for the sieve in the right arm only
	//need to mirror(flip) it in order to put in the left arm  
	double pSieveSlitDeltaH[8]={0.537*inch, 0.188*inch, 0.188*inch, 0.188*inch,
		0.241*inch, 0.241*inch, 0.241*inch, 0.381*inch}; //in inch, from left to right
	double pSieveSlitDeltaV[8]={0.496*inch, 0.524*inch, 0.524*inch, 0.524*inch,
		0.524*inch, 0.524*inch, 0.524*inch, 1.494*inch}; //in inch, from top to buttom

	//the whole position relative to the slit center 
	double pSieveSlitHolePosH[7], pSieveSlitHolePosV[7];
	for(int ii=0;ii<7;ii++)
	{
		pSieveSlitHolePosH[ii] = (ii==0)?pSieveSlitX/2.0:pSieveSlitHolePosH[ii-1];
		pSieveSlitHolePosH[ii] -= pSieveSlitDeltaH[ii];

		pSieveSlitHolePosV[ii] = (ii==0)?pSieveSlitY/2.0:pSieveSlitHolePosV[ii-1];
		pSieveSlitHolePosV[ii] -= pSieveSlitDeltaV[ii];
	}

	////////////////the following is for 12.5 degree sieve//////////////	
	
	if(!mSetupSeptumBlock)
	{
		//these 2 will be read from ini file, it is 80.0cm for 6 deg and 120 cm for 12.5 deg
		//mPivot2LSieveFace=115.6*cm;
		//mPivot2RSieveFace=116.39*cm;

		pSieveSlitX=5.687*inch; //144.45*mm
		pSieveSlitY=7.750*inch; //196.85*mm
		pSieveSlitZ=0.25*inch;

		pSieveSlitHoleR=1.0033*mm;           //radius of small hole 0.079/2 inch
		pSieveSlitLargeHoleR=1.994*mm;       //radius of large hole 0.157/2 inch
		pSieveSlitHoldL=pSieveSlitZ+0.1*mm;  //need to make it longer to avoid round off in the subtraction

		double pFirstHole2Top=0.799*inch, pFirstHole2Left=1.306*inch;
		double pHoleSpanH=0.513*inch, pHoleSpanV=1.025*inch;
		for(int ii=0;ii<7;ii++)
		{
			pSieveSlitHolePosH[ii] = (ii==0)?pSieveSlitX/2.0-pFirstHole2Left:pSieveSlitHolePosH[ii-1]-pHoleSpanH;
			pSieveSlitHolePosV[ii] = (ii==0)?pSieveSlitY/2.0-pFirstHole2Top:pSieveSlitHolePosV[ii-1]-pHoleSpanV;
		}

		//the big center hole horizontal and vertical offset, From Alan
		pSieveSlitLargeHoleH=0.0;
		pSieveSlitLargeHoleV=0.0;
	}


	//now start to build box then subtract 49 holes 
	G4VSolid* sieveSlitWholeSolid=new G4Box("sieveSlitWholeBox",pSieveSlitX/2.0,
		pSieveSlitY/2.0,pSieveSlitZ/2.0);
	startphi=0.0*deg;deltaphi=360.0*deg;
	G4VSolid* sieveSlitHoleSolid=new G4Tubs("sieveSlitHoleTubs",0,pSieveSlitHoleR,
		pSieveSlitHoldL/2.0,startphi,deltaphi); 
	G4VSolid* sieveSlitLargeHoleSolid=new G4Tubs("sieveSlitLargeHoleTubs",0,
		pSieveSlitLargeHoleR,pSieveSlitHoldL/2.0,startphi,deltaphi); 

	G4SubtractionSolid* sieveSlitSolid=(G4SubtractionSolid*)sieveSlitWholeSolid;
	char strName[100];
	for(int ih=0;ih<7;ih++)
	{
		for(int iv=0;iv<7;iv++)
		{
			sprintf(strName,"sieveSlitHole_H%d_V%d",ih,iv);
			if((ih==3 && iv==3) || (ih==4 && iv==1)) 
			{
				//now dig large holes in the block
				sieveSlitSolid=new G4SubtractionSolid(strName,sieveSlitSolid,
					sieveSlitLargeHoleSolid,0,
					G4ThreeVector(pSieveSlitHolePosH[ih],pSieveSlitHolePosV[iv],0));
			}
			else
			{
				//now dig small holes in the block
				sieveSlitSolid=new G4SubtractionSolid(strName,sieveSlitSolid,
					sieveSlitHoleSolid,0,
					G4ThreeVector(pSieveSlitHolePosH[ih],pSieveSlitHolePosV[iv],0));
			}
		}
	}
	sieveSlitSolid->SetName("sieveSlitSolid");

	G4LogicalVolume* sieveSlitLogical = new G4LogicalVolume(sieveSlitSolid,
		mMaterialManager->tungsten,"sieveSlitLogical",0,0,0);
	sieveSlitLogical->SetVisAttributes(LeadVisAtt); 

	SDman->AddNewDetector(sieveSlitSD);
	sieveSlitLogical->SetSensitiveDetector(sieveSlitSD);


	//calculate the center position in the Lab frame
	double pSieveSlitCenterHOffset=pSieveSlitLargeHoleH-pSieveSlitHolePosH[3];
	double pSieveSlitCenterVOffset=pSieveSlitLargeHoleV-pSieveSlitHolePosV[3];


	//place the sieve slits in the hall
	double pLSieveSlitPos_X=(mPivot2LSieveFace+pSieveSlitZ/2.0)*sin(mLSeptumAngle)+
		pSieveSlitCenterHOffset+mPivotXOffset;
	double pLSieveSlitPos_Y=pSieveSlitCenterVOffset+mPivotYOffset;
	double pLSieveSlitPos_Z=(mPivot2LSieveFace+pSieveSlitZ/2.0)*cos(mLSeptumAngle)+mPivotZOffset;
	double pRSieveSlitPos_X=(mPivot2RSieveFace+pSieveSlitZ/2.0)*sin(mRSeptumAngle)-
		pSieveSlitCenterHOffset+mPivotXOffset;
	double pRSieveSlitPos_Y=pSieveSlitCenterVOffset+mPivotYOffset;
	double pRSieveSlitPos_Z=(mPivot2RSieveFace+pSieveSlitZ/2.0)*cos(mRSeptumAngle)+mPivotZOffset;

	if(mSetupLSieveSlit)
	{ 
		G4RotationMatrix *pRotLSieve=new G4RotationMatrix();
		pRotLSieve->rotateY(-mLSeptumAngle-180*deg);
		new G4PVPlacement(pRotLSieve,
			G4ThreeVector(pLSieveSlitPos_X,pLSieveSlitPos_Y,pLSieveSlitPos_Z),
			sieveSlitLogical,"leftSieveSlitPhys",motherLogical,0,0);
	}
	if(mSetupRSieveSlit)
	{
		new G4PVPlacement(pRotRSeptum,
			G4ThreeVector(pRSieveSlitPos_X,pRSieveSlitPos_Y,pRSieveSlitPos_Z),
			sieveSlitLogical,"rightSieveSlitPhys",motherLogical,0,0);
	}

	/////////////////////////
	// Septum block 
	/////////////////////////
	//by Jixie: Allow one to setup septum without HRS

	if(mSetupSeptumBlock)
	{
		////////////////////////////////////////////////////////////////////
		//Septum block, 140 cm width, 84.4 cm height and 74 cm in length, silicon steel
		//Tunnel size: started from x=8.4cm, 30.4cm wide and 24.4 cm in height
		//located at z=700 mm, no rotation
		double pSeptumX=140.0*cm;
		double pSeptumY=84.4*cm;
		double pSeptumZ=74.0*cm;
		double pSeptumTunnelX=30.4*cm;
		double pSeptumTunnelY=24.4*cm-2.0*inch;  //By Jixie @20120205: Add 2 inches of iron
		double pSeptumBeamHoleX=7.8*cm;
		double pSeptumBeamHoleY=8.0*cm;

		double pSeptumTunnelPos_X=8.4*cm+pSeptumTunnelX/2.0;
		double pSeptumPos_Z=70.0*cm;		
		gConfig->GetParameter("Septum_OriginZ",pSeptumPos_Z); 
		pSeptumPos_Z*=cm;

		G4VSolid* septumBlockSolid = new G4Box("septumBlockBox",pSeptumX/2.0,
			pSeptumY/2.0,pSeptumZ/2.0);

		//Left and right tunnels, treat the cu coils as part of the block
		//By Jixie: I reduced this by 0.2cm for the Helium bag
		G4VSolid* septumTunnelSolid = new G4Box("septumTunnelBox",pSeptumTunnelX/2.0-0.2*cm,
			pSeptumTunnelY/2.0-0.2*cm,pSeptumZ/2.0+1.0*mm);

		//beam pine hole
		G4VSolid* septumBeamHoleSolid = new G4Box("septumBeamHoleBox",pSeptumBeamHoleX/2.0,
			pSeptumBeamHoleY/2.0,pSeptumZ/2.0+1.0*mm);

		//dig 3 holes, left, right tunnel and beam hole
		G4SubtractionSolid* septumBlockSubLSolid=new G4SubtractionSolid("septumBlockSubL",
			septumBlockSolid,septumTunnelSolid,0,G4ThreeVector(pSeptumTunnelPos_X,0,0));
		G4SubtractionSolid* septumBlockSubLRSolid=new G4SubtractionSolid("septumBlockSubLR",
			septumBlockSubLSolid,septumTunnelSolid,0,G4ThreeVector(-pSeptumTunnelPos_X,0,0));
		G4SubtractionSolid* septumBlockSubLRCSolid=new G4SubtractionSolid("septumBlockSubLRC",
			septumBlockSubLRSolid,septumBeamHoleSolid);

		G4LogicalVolume* septumLogical = new G4LogicalVolume(septumBlockSubLRCSolid,
			mMaterialManager->siliconsteel,"septumLogical",0,0,0);
		septumLogical->SetVisAttributes(IronVisAtt);

		//put it in the hall, no rotation
		new G4PVPlacement(0,G4ThreeVector(0,0,pSeptumPos_Z),
			septumLogical,"septumPhys",motherLogical,0,0,0);


		double pSeptumCoilRadius = 22.06*cm;
		double pSeptumCoilCenterX = 2.83*cm;
		double pSeptumCoilCenterY = -10.25*cm;
		double pSeptumCoilThickness = 4.5*cm;
		G4VSolid* septumCoilCylinderSolid = new G4Tubs("septumCoilCylinderTub",
			0,pSeptumCoilRadius,pSeptumCoilThickness/2.0,0,360*deg);

		double septumCoilRecX = 25.0*cm;
		double septumCoilRecY = 15.0*cm;
		G4VSolid* septumCoilRecSolid = new G4Box("septumCoilRecBox",
			septumCoilRecX/2.0,septumCoilRecY/2.0,pSeptumCoilThickness/2.0);

		G4IntersectionSolid* septumCoilSolid = new G4IntersectionSolid("septumCoilSolid",
			septumCoilCylinderSolid,septumCoilRecSolid,0,
			G4ThreeVector(-pSeptumCoilCenterX-septumCoilRecX/2.0,-pSeptumCoilCenterY+septumCoilRecY/2.0,0));

		G4LogicalVolume* septumCoilLogical = new G4LogicalVolume(septumCoilSolid,
			mMaterialManager->copper,"septumCoilLogical",0,0,0);
		septumCoilLogical->SetVisAttributes(CuBrownVisAtt);

		//place 16 copies into the septum container
		double pSeptumCoilPos_X_in   = pSeptumTunnelPos_X-pSeptumTunnelX/2.0-pSeptumCoilThickness/2.0;
		double pSeptumCoilPos_X_out  = pSeptumTunnelPos_X+pSeptumTunnelX/2.0+pSeptumCoilThickness/2.0;
		double pSeptumCoilPos_Y      = pSeptumCoilCenterY-pSeptumTunnelY/2.0;
		double pSeptumCoilPos_Z_up   = pSeptumPos_Z-pSeptumZ/2.0+pSeptumCoilCenterX;
		double pSeptumCoilPos_Z_down = pSeptumPos_Z+pSeptumZ/2.0-pSeptumCoilCenterX;


		//place up stream coils in the following order (looking downstream)
		//#####2######1###7######6#####
		//######      |###|      |#####
		//######      |###|      |#####
		//#####3######0###4######5#####
		G4RotationMatrix* pSeptumCoilRotFrontDown = new G4RotationMatrix();
		pSeptumCoilRotFrontDown->rotateY(90*deg);
		G4RotationMatrix* pSeptumCoilRotFrontUp = new G4RotationMatrix();
		pSeptumCoilRotFrontUp->rotateY(90*deg);
		pSeptumCoilRotFrontUp->rotateX(180*deg);

		new G4PVPlacement(pSeptumCoilRotFrontDown,
			G4ThreeVector(pSeptumCoilPos_X_in,pSeptumCoilPos_Y,pSeptumCoilPos_Z_up),
			septumCoilLogical,"septumCoilPhys",motherLogical,true,0,0);
		new G4PVPlacement(pSeptumCoilRotFrontUp,
			G4ThreeVector(pSeptumCoilPos_X_in,-pSeptumCoilPos_Y,pSeptumCoilPos_Z_up),
			septumCoilLogical,"septumCoilPhys",motherLogical,true,1,0);
		new G4PVPlacement(pSeptumCoilRotFrontUp,
			G4ThreeVector(pSeptumCoilPos_X_out,-pSeptumCoilPos_Y,pSeptumCoilPos_Z_up),
			septumCoilLogical,"septumCoilPhys",motherLogical,true,2,0);
		new G4PVPlacement(pSeptumCoilRotFrontDown,
			G4ThreeVector(pSeptumCoilPos_X_out,pSeptumCoilPos_Y,pSeptumCoilPos_Z_up),
			septumCoilLogical,"septumCoilPhys",motherLogical,true,3,0);

		new G4PVPlacement(pSeptumCoilRotFrontDown,
			G4ThreeVector(-pSeptumCoilPos_X_out,pSeptumCoilPos_Y,pSeptumCoilPos_Z_up),
			septumCoilLogical,"septumCoilPhys",motherLogical,true,4,0);
		new G4PVPlacement(pSeptumCoilRotFrontUp,
			G4ThreeVector(-pSeptumCoilPos_X_out,-pSeptumCoilPos_Y,pSeptumCoilPos_Z_up),
			septumCoilLogical,"septumCoilPhys",motherLogical,true,5,0);
		new G4PVPlacement(pSeptumCoilRotFrontUp,
			G4ThreeVector(-pSeptumCoilPos_X_in,-pSeptumCoilPos_Y,pSeptumCoilPos_Z_up),
			septumCoilLogical,"septumCoilPhys",motherLogical,true,6,0);
		new G4PVPlacement(pSeptumCoilRotFrontDown,
			G4ThreeVector(-pSeptumCoilPos_X_in,pSeptumCoilPos_Y,pSeptumCoilPos_Z_up),
			septumCoilLogical,"septumCoilPhys",motherLogical,true,7,0);


		//place down stream coils in the following order (looking downstream)
		//####10######9###15#####14####
		//######      |###|      |#####
		//######      |###|      |#####
		//####11######8###12#####13####
		G4RotationMatrix* pSeptumCoilRotBackDown = new G4RotationMatrix();
		pSeptumCoilRotBackDown->rotateY(270*deg);
		G4RotationMatrix* pSeptumCoilRotBackUp = new G4RotationMatrix();
		pSeptumCoilRotBackUp->rotateY(270*deg);
		pSeptumCoilRotBackUp->rotateX(180*deg);

		new G4PVPlacement(pSeptumCoilRotBackDown,
			G4ThreeVector(pSeptumCoilPos_X_in,pSeptumCoilPos_Y,pSeptumCoilPos_Z_down),
			septumCoilLogical,"septumCoilPhys",motherLogical,true,8,0);
		new G4PVPlacement(pSeptumCoilRotBackUp,
			G4ThreeVector(pSeptumCoilPos_X_in,-pSeptumCoilPos_Y,pSeptumCoilPos_Z_down),
			septumCoilLogical,"septumCoilPhys",motherLogical,true,9,0);
		new G4PVPlacement(pSeptumCoilRotBackUp,
			G4ThreeVector(pSeptumCoilPos_X_out,-pSeptumCoilPos_Y,pSeptumCoilPos_Z_down),
			septumCoilLogical,"septumCoilPhys",motherLogical,true,10,0);
		new G4PVPlacement(pSeptumCoilRotBackDown,
			G4ThreeVector(pSeptumCoilPos_X_out,pSeptumCoilPos_Y,pSeptumCoilPos_Z_down),
			septumCoilLogical,"septumCoilPhys",motherLogical,true,11,0);

		new G4PVPlacement(pSeptumCoilRotBackDown,
			G4ThreeVector(-pSeptumCoilPos_X_out,pSeptumCoilPos_Y,pSeptumCoilPos_Z_down),
			septumCoilLogical,"septumCoilPhys",motherLogical,true,12,0);
		new G4PVPlacement(pSeptumCoilRotBackUp,
			G4ThreeVector(-pSeptumCoilPos_X_out,-pSeptumCoilPos_Y,pSeptumCoilPos_Z_down),
			septumCoilLogical,"septumCoilPhys",motherLogical,true,13,0);
		new G4PVPlacement(pSeptumCoilRotBackUp,
			G4ThreeVector(-pSeptumCoilPos_X_in,-pSeptumCoilPos_Y,pSeptumCoilPos_Z_down),
			septumCoilLogical,"septumCoilPhys",motherLogical,true,14,0);
		new G4PVPlacement(pSeptumCoilRotBackDown,
			G4ThreeVector(-pSeptumCoilPos_X_in,pSeptumCoilPos_Y,pSeptumCoilPos_Z_down),
			septumCoilLogical,"septumCoilPhys",motherLogical,true,15,0);
	}
	
	/////////////////////////
	// HRS Virtual Boundary, 
	/////////////////////////
	//For 6 degrees, if septum field is valid and argument UseSeptumPlusStdHRS==1, will 
	//place the virtual boundary 1 cm before HRSContainer (pHRSContainerRin-6*cm = 1.40m)
	//otherwise will place the VB at the position given by the Detector.ini
	//For 12.5 degrees, always place VB 1.40m away from the hall center

	int pUseSeptumPlusStdHRS=0;
	gConfig->GetArgument("UseSeptumPlusStdHRS",pUseSeptumPlusStdHRS);
	if((pUseSeptumField && pUseSeptumPlusStdHRS) || !mSetupSeptumBlock) 
	{
		//this part is trying to place a virtual boundary 1 cm in front of the HRSContainer
		//It is for the case that we use the septum field to propogate electrons to Q1 entrance 
		//and then use the STD HRS transportation other than use 5.69 degrees HRS transportation

		//Place both left and right VB for HRS, which is pHRSContainerRin-6.0*cm away from the 
		//hall center(1.40m). This aperture is a round disk of 30 cm diameter
		//The real Q1 vacumn entrance to hall center is 1.312m, 

		G4VSolid* HRSVBSolid = new G4Tubs("HRSVBTub",0.0,15*cm,
			mHRSVBThick/2.0,0.0,360.0*deg);
		G4LogicalVolume* HRSVBLogical = new G4LogicalVolume(HRSVBSolid,
			mMaterialManager->mylar,"HRSVBLogical",0,0,0);
		SDman->AddNewDetector(septumWindowSD);
		HRSVBLogical->SetSensitiveDetector(septumWindowSD);
		HRSVBLogical->SetVisAttributes(LightYellowVisAtt); 

		double pHallCenter2VB=1.40*m;
		gConfig->SetParameter("Pivot2LHRSVBFace",pHallCenter2VB-mPivotZOffset*cos(mLHRSAngle));
		gConfig->SetParameter("Pivot2RHRSVBFace",pHallCenter2VB-mPivotZOffset*cos(mRHRSAngle)); 

		if(mSetupLHRS)
		{
			double pLHRSVBPos_X=(pHallCenter2VB-mHRSVBThick/2)*sin(mLHRSAngle)+mPivotXOffset;
			double pLHRSVBPos_Y=mPivotYOffset;
			//no need to correct for pivot since the distance is from the hall center
			double pLHRSVBPos_Z=(pHallCenter2VB-mHRSVBThick/2.0)*cos(mLHRSAngle);
			new G4PVPlacement(pRotLHRS,G4ThreeVector(pLHRSVBPos_X,pLHRSVBPos_Y,pLHRSVBPos_Z),
				HRSVBLogical,"virtualBoundaryPhys_LHRS",motherLogical,0,0,0);
		}
		if(mSetupRHRS)
		{
			double pRHRSVBPos_X=(pHallCenter2VB-mHRSVBThick/2)*sin(mRHRSAngle)+mPivotXOffset;
			double pRHRSVBPos_Y=mPivotYOffset;
			//no need to correct for pivot since the distance is from the hall center
			double pRHRSVBPos_Z=(pHallCenter2VB-mHRSVBThick/2)*cos(mRHRSAngle); 
			new G4PVPlacement(pRotRHRS,G4ThreeVector(pRHRSVBPos_X,pRHRSVBPos_Y,pRHRSVBPos_Z),
				HRSVBLogical,"virtualBoundaryPhys_RHRS",motherLogical,0,0,0);
		}
	}
	else
	{
		//place VB @ Septum entrance window, 10.4cm width and 24.4cm height, 
		//The following declared as module variales already
		//double mHRSVBWidth=104*mm;
		//double mHRSVBHeight=244*mm;
		//double mHRSVBThick=0.0508*mm; 
		//acceptance is 20 mrad
		G4VSolid* septumWindowSolid = new G4Box("septumWindowBox",mHRSVBWidth/2.0,
			mHRSVBHeight/2.0,mHRSVBThick/2.0);
		G4LogicalVolume* septumWindowLogical = new G4LogicalVolume(septumWindowSolid,
			mMaterialManager->mylar,"septumWindowLogical",0,0,0);
		SDman->AddNewDetector(septumWindowSD);
		septumWindowLogical->SetSensitiveDetector(septumWindowSD);
		septumWindowLogical->SetVisAttributes(LightYellowVisAtt); 


		double pTunnel2Beam_X=8.4*cm;
		//put both left and right septum entrance window, which should not less than pTunnel2Beam_X
		if(mSetupLHRS)
		{
			double pLSeptumWindowPos_X=(mPivot2LHRSVBFace+mHRSVBThick/2.0)*
				sin(mLSeptumAngle)+mPivotXOffset;
			if(mPivot2LHRSVBFace>1190*mm)
			{
				//need to shift it in X to make it barely touch the septum tunnel		
				pLSeptumWindowPos_X=mHRSVBWidth/2*cos(mLSeptumAngle)+pTunnel2Beam_X;
			}
			double pLSeptumWindowPos_Y=mPivotYOffset;
			double pLSeptumWindowPos_Z=(mPivot2LHRSVBFace+mHRSVBThick/2.0)*
				cos(mLSeptumAngle)+mPivotZOffset;
			new G4PVPlacement(pRotLSeptum,
				G4ThreeVector(pLSeptumWindowPos_X,pLSeptumWindowPos_Y,pLSeptumWindowPos_Z),
				septumWindowLogical,"virtualBoundaryPhys_LHRS",motherLogical,0,0,0);
		}
		if(mSetupRHRS)
		{
			double pRSeptumWindowPos_X=(mPivot2RHRSVBFace+mHRSVBThick/2.)*
				sin(mRSeptumAngle)+mPivotXOffset;			
			if(mPivot2LHRSVBFace>1190*mm)
			{
				//need to shift it in X to make it barely touch the septum tunnel		
				pRSeptumWindowPos_X=-mHRSVBWidth/2*cos(mRSeptumAngle)-pTunnel2Beam_X;
			}
			double pRSeptumWindowPos_Y=mPivotYOffset;
			double pRSeptumWindowPos_Z=(mPivot2RHRSVBFace+mHRSVBThick/2.0)*
				cos(mRSeptumAngle)+mPivotZOffset;
			new G4PVPlacement(pRotRSeptum,
				G4ThreeVector(pRSeptumWindowPos_X,pRSeptumWindowPos_Y,pRSeptumWindowPos_Z),
				septumWindowLogical,"virtualBoundaryPhys_RHRS",motherLogical,0,0,0);
		}
	}


	return theG2PSeptumPhys;
}


//build DVCS solenoid, will move to another file later
G4VPhysicalVolume* CREXDetectorConstruction::ConstructDVCSSolenoid(G4LogicalVolume *pMotherLogVol)
{
	double startphi=0.*deg, deltaphi=360.*deg;

	G4VPhysicalVolume* thePhysVol=0;
	int mSetupDVCSSolenoid=1;
	//build the DVCS solenoid, Rin=230mm and Rout=912 mm
	//I build the whole body, not just the coil. Those numbers are measured from
	//pdf file: http://wwwold.jlab.org/Hall-B/secure/e1-dvcs/michel/sol/drawings/Ensemble.pdf
	//
	//|<--------- 510 mm -------->|<--228.8-->|	      
	//-----------------------------------------        --
	//|                                       |         132.7 mm in vertical
	//|                                       |        __
	//|                                      /                             p
	//|                                    /                               p
	//|                                  /                                 p 
	//-----------------------------------
	//|
	//|                           +(center)
	//|                        -->|     |<--106.2 mm 
	//-----------------------------------
	//|                                  \                                 p
	//|                                    \                               p 
	//|                                      \                             p
	//|                                       |
	//|                                       |  
	//-----------------------------------------


	if(mSetupDVCSSolenoid>=1)
	{
		startphi=0.*deg; deltaphi=360.*deg;
		const int kNPlane_DVCSCoil=3;
		double rInner_DVCSCoil[] = {115*mm,115*mm,313.3*mm};
		double rOuter_DVCSCoil[] = {456*mm,456*mm,456*mm};
		double zPlane_DVCSCoil[] = {-510*mm,106.2*mm,228.8*mm};

		G4Polycone* DVCSCoilSolid = new G4Polycone("DVCSCoilPolycone",startphi,deltaphi,
			kNPlane_DVCSCoil,zPlane_DVCSCoil,rInner_DVCSCoil,rOuter_DVCSCoil);
	
		G4LogicalVolume* DVCSCoilLogical = new G4LogicalVolume(DVCSCoilSolid,
			mMaterialManager->stainlesssteel,"DVCSCoilLogical",0,0,0);

		DVCSCoilLogical->SetVisAttributes(SteelVisAtt); 
		
		//define the center position, these number should be available in the ini file
		double DVCSCoilPos_x=0*mm;
		double DVCSCoilPos_y=0*mm;
		double DVCSCoilPos_z=-100*mm;
		thePhysVol=new G4PVPlacement(0,G4ThreeVector(DVCSCoilPos_x,DVCSCoilPos_y,DVCSCoilPos_z),
			DVCSCoilLogical,"DVCSCoilPhys",pMotherLogVol,false,0);
	}

	return thePhysVol;
}


