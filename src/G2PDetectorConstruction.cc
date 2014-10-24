// ********************************************************************
// $Id: G2PDetectorConstruction.cc,v 1.02, 2010/12/26 HRS Exp $
//
// ********************************************************************
//
#include <stdio.h>
#include <math.h>
#include <iostream>
using namespace std;

#include "G2PDetectorConstruction.hh"

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

#include "HRSStdSD.hh"
#include "HRSDCSD.hh"
#include "UsageManager.hh"
#include "HRSMaterial.hh"
#include "HRSEMFieldSetup.hh"

//To verify some geometry, I use this flag to place some unit
//just for debugging
//#define G4DEBUG_GEOMETRY 0

//The G4RotationMatrix rotates the whole coordinate system, 
//Looking from top, it always rotate clockwise  

extern UsageManager* gConfig;	

/////////////////////////////////////////////////////////////////////
G2PDetectorConstruction::G2PDetectorConstruction(G4LogicalVolume *mother) : 
mMotherLogVol(mother) 
{
	GetConfig();
	mMaterialManager=HRSMaterial::GetHRSMaterialManager();
	
	ConstructMaterial();
	G4cout<<"Contrstruct G2P geometry ... done! "<<G4endl;
}

G2PDetectorConstruction::~G2PDetectorConstruction()
{
	if (G2P_NH3He)  delete G2P_NH3He;
	G4cout<<"Delete G2P geometry ... done! "<<G4endl;
}

/////////////////////////////////////////////////////////////////////
void G2PDetectorConstruction::GetConfig()
{
	gConfig->ReadFile("Detector_G2P.ini");
	//////////////////////////////////////////////////////////////////
	//the following is just for g2p

	gConfig->GetParameter("ScatChamberRin",mScatChamberRin);
	mScatChamberRin*=mm;
	gConfig->GetParameter("ScatChamberRout",mScatChamberRout);
	mScatChamberRout*=mm;
	gConfig->GetParameter("ScatChamberL",mScatChamberL);
	mScatChamberL*=mm;

	gConfig->GetParameter("ScatChamberExitWindowThick",mScatChamberExitWindowThick);
	mScatChamberExitWindowThick*=mm;
	gConfig->GetParameter("ShieldLN2WindowThick",mShieldLN2WindowThick);
	mShieldLN2WindowThick*=mm;

	gConfig->GetParameter("ShieldLN2Rin",mShieldLN2Rin);
	mShieldLN2Rin*=mm;
	gConfig->GetParameter("ShieldLN2Rout",mShieldLN2Rout);
	mShieldLN2Rout*=mm;
	gConfig->GetParameter("ShieldLN2L",mShieldLN2L);
	mShieldLN2L*=mm;
	gConfig->GetParameter("ShieldLHeRin",mShieldLHeRin);
	mShieldLHeRin*=mm;
	gConfig->GetParameter("ShieldLHeRout",mShieldLHeRout);
	mShieldLHeRout*=mm;
	gConfig->GetParameter("ShieldLHeL",mShieldLHeL);
	mShieldLHeL*=mm;

	gConfig->GetParameter("TargetType",mTargetType);
	mSetupG2PTarget=1;
	gConfig->GetParameter("SetupG2PTarget",mSetupG2PTarget);

	gConfig->GetParameter("TargetL",mTargetL);
	mTargetL*=mm;

	gConfig->GetParameter("SetupG2PScatChamber",mSetupG2PScatChamber);
	gConfig->GetParameter("SetupCoil",mSetupCoil);

	gConfig->GetParameter("SetupBeamDump",mSetupBeamDump);
	gConfig->GetParameter("BeamDumpWidth",mBeamDumpWidth);
	mBeamDumpWidth*=mm;
	gConfig->GetParameter("BeamDumpHeight",mBeamDumpHeight);
	mBeamDumpHeight*=mm;
	gConfig->GetParameter("BeamDumpThick",mBeamDumpThick);
	mBeamDumpThick*=mm;
	gConfig->GetParameter("Pivot2BeamDumpFace",mPivot2BeamDumpFace);
	mPivot2BeamDumpFace*=mm;


	gConfig->GetParameter("SetupThirdArm",mSetupThirdArm);
	gConfig->GetParameter("ThirdArmAngle",mThirdArmAngle);
	mThirdArmAngle*=deg;
	gConfig->GetParameter("SetupThirdArmVD",mSetupThirdArmVD);
	gConfig->GetParameter("ThirdArmRotZAngle",mThirdArmRotZAngle);
	mThirdArmRotZAngle*=deg;
	gConfig->GetParameter("Pivot2ThirdArmFace",mPivot2ThirdArmFace);
	mPivot2ThirdArmFace*=mm;

	gConfig->GetParameter("SetupChicane",mSetupChicane);
	gConfig->GetParameter("SetupChicaneVD",mSetupChicaneVD);

	mSetupPlatform=0;
	gConfig->GetParameter("SetupPlatform",mSetupPlatform);
	
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

	///////////////////////////////////////////////////////////////
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

	G4cout<<"\n****Load G2P detector config parameters done!***"<<G4endl;

	return ;
}

void G2PDetectorConstruction::ConstructMaterial()
{
	//Due to the fact that the HRSMaterial is built before reading Detector_G2P.ini
	//The mMaterialManager->NH3He was built with 55% NH3 volumn ratio. To allow user 
	//to change the packing ratio, I have to create the material here.

	double pSolidNH3D;	
	gConfig->GetParameter("SolidNH3D",pSolidNH3D);
	pSolidNH3D*=mg/cm3;

	double pLiquidHeD;	
	gConfig->GetParameter("LiquidHeD",pLiquidHeD);
	pLiquidHeD*=mg/cm3;

	double pNH3WeightRatio;	
	gConfig->GetParameter("NH3WeightRatio",pNH3WeightRatio);
	
	double density;
	density = 1.0/((1.0-pNH3WeightRatio)/pLiquidHeD+pNH3WeightRatio/pSolidNH3D);
	G2P_NH3He = new G4Material("G2P_NH3He", density, 2);
	G2P_NH3He->AddMaterial(mMaterialManager->solidNH3, pNH3WeightRatio);
	G2P_NH3He->AddMaterial(mMaterialManager->liquidHe, 1.0-pNH3WeightRatio);
}

/////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* G2PDetectorConstruction::Construct()
{	
	G4VPhysicalVolume* thePhysVol=0;

	/////////////////////////
	// G2P chicane magnets
	/////////////////////////
	if(mSetupChicane) this->ConstructG2PChicane(mMotherLogVol);

	/////////////////////////
	// G2P local beam dump
	/////////////////////////
	if(mSetupBeamDump) this->ConstructG2PBeamDump(mMotherLogVol);

	/////////////////////////
	//g2p scatter chamber 
	/////////////////////////
	//since the coil will always rotate with the chamber, I use a container to hold them.
	//the ScatChamberContainer start from mShieldLHeRin to 
	//mScatChamberRout+pScatChamberEntranceWindow2UnionL+10*mm
	//The advantage is that all staff inside the containner have the same rotation 
	if(mSetupG2PScatChamber) this->ConstructG2PScatChamber(mMotherLogVol);


	/////////////////////////
	//g2p target container and target vessel
	/////////////////////////	
	if(mSetupG2PTarget==99) this->ConstructG2PPVCTarget(mMotherLogVol);
	else if(mSetupG2PTarget) this->ConstructG2PTarget(mMotherLogVol);

	
	/////////////////////////
	// Sieve slit, septum window and HRS VirtualBoundary
	/////////////////////////
	if(mSetupLHRS || mSetupRHRS)  
	{
		this->ConstructG2PSeptumNSieve(mMotherLogVol);
		//this->ConstructHRS(mMotherLogVol);
	}

	/////////////////////////
	// 3rd arm
	/////////////////////////
	if(mSetupThirdArm)  this->ConstructG2PThirdArm(mMotherLogVol);


	/////////////////////////
	// the platform
	/////////////////////////
	if(mSetupPlatform)  this->ConstructG2PPlatform(mMotherLogVol);

	return thePhysVol;
}

/////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* G2PDetectorConstruction::ConstructG2PScatChamber(G4LogicalVolume* motherLogical)
{
	const double inch=2.54*cm;
	double startphi,deltaphi;

	G4RotationMatrix* pRotX90deg = new G4RotationMatrix();
	pRotX90deg->rotateX(90.*deg);
	G4RotationMatrix* pRotX180deg = new G4RotationMatrix();
	pRotX180deg->rotateX(180.*deg);
	G4RotationMatrix* pRotX270deg = new G4RotationMatrix();
	pRotX270deg->rotateX(270.*deg);

	G4RotationMatrix* pRotY90deg = new G4RotationMatrix();
	pRotY90deg->rotateY(90.*deg);
	G4RotationMatrix* pRotY180deg = new G4RotationMatrix();
	pRotY180deg->rotateY(180.*deg);
	G4RotationMatrix* pRotY270deg = new G4RotationMatrix();
	pRotY270deg->rotateY(270.*deg);

	/////////////////////////////////////////////////////////
	//By default, the scatter chamber and field are in longitudinal status,
	//For transverse configuration, we must rotate the scatter chamber and coil
	//about Y axis by -90 degrees, minus here means anti-clockwise in the overlook view 

	//Note that the scatter chamber is put in a vertical orientation
	//By Jixie: I have verified that rotate by matrix A then by Matrix B should be written
	//as Rotate=A*B, always put the first rotation at the left

#ifdef G4DEBUG_GEOMETRY
	if(G4DEBUG_GEOMETRY>=1)
	{
		G4RotationMatrix *pRotBField=new G4RotationMatrix();
		BField_Helm *pHelmField=BField_Helm::GetInstance();
		pRotBField=(G4RotationMatrix*)(pHelmField->GetRotation_L2F());

		G4RotationMatrix* pRotX90deg = new G4RotationMatrix();
		pRotX90deg->rotateX(90.*deg);
		G4RotationMatrix *pRotScatChamber=new G4RotationMatrix();
		*pRotScatChamber = (*pRotX90deg) * (*pRotBField) ; //rotate with the coil	

		cout<<"\n pRotBField ==> EulerAngles: phi="<<pRotBField->getPhi()/deg
			<<"  theta="<<pRotBField->getTheta()/deg
			<<"  psi="<<pRotBField->getPsi()/deg<<endl;

		cout<<"\n pRotScatChamber ==> EulerAngles: phi="<<pRotScatChamber->getPhi()/deg
			<<"  theta="<<pRotScatChamber->getTheta()/deg
			<<"  psi="<<pRotScatChamber->getPsi()/deg<<endl;
	}
#endif

	// Magnetic field Rotation----------------------------------------------------------
	BField_Helm *pHelmField=BField_Helm::GetInstance();
	G4RotationMatrix *pRotBField=(G4RotationMatrix*)(pHelmField->GetRotation_L2F());
	G4RotationMatrix *pRotScatInHall=new G4RotationMatrix();
	*pRotScatInHall = (*pRotX90deg) * (*pRotBField) ; //rotate with the coil

	/////////////////////////////////////////////////////////
#ifdef G4DEBUG_GEOMETRY
	//this part is used to compare with pRotScatInHall, this is the way I figure out
	//how pRotScatInHall should be constructed

	G4RotationMatrix *pRotX90Z270deg=new G4RotationMatrix();
	pRotX90Z270deg->rotateX(90*deg);
	pRotX90Z270deg->rotateZ(270*deg); 

	G4RotationMatrix *pRotX90Z353deg=new G4RotationMatrix();
	pRotX90Z353deg->rotateX(90*deg);
	pRotX90Z353deg->rotateZ(353*deg); 

	if(G4DEBUG_GEOMETRY>=1)
	{
		cout<<"\n pRotBField ==> EulerAngles: phi="<<pRotBField->getPhi()/deg
			<<"  theta="<<pRotBField->getTheta()/deg
			<<"  psi="<<pRotBField->getPsi()/deg<<endl;

		cout<<"\n pRotX90deg ==> EulerAngles: phi="<<pRotX90deg->getPhi()/deg
			<<"  theta="<<pRotX90deg->getTheta()/deg
			<<"  psi="<<pRotX90deg->getPsi()/deg<<endl;

		cout<<"\n pRotX90Z270deg ==> EulerAngles: phi="<<pRotX90Z270deg->getPhi()/deg
			<<"  theta="<<pRotX90Z270deg->getTheta()/deg
			<<"  psi="<<pRotX90Z270deg->getPsi()/deg<<endl;

		cout<<"\n pRotX90Z353deg ==> EulerAngles: phi="<<pRotX90Z353deg->getPhi()/deg
			<<"  theta="<<pRotX90Z353deg->getTheta()/deg
			<<"  psi="<<pRotX90Z353deg->getPsi()/deg<<endl;

		//The following is used to verify that rotate by matrix A then rotate by Matrix B should 
		//be written as Rotate=A*B, always put the first rotation at the left
		cout<<"\n (*pRotX90deg) * (*pRotBField) ==> EulerAngles: phi="<<pRotScatInHall->getPhi()/deg
			<<"  theta="<<pRotScatInHall->getTheta()/deg
			<<"  psi="<<pRotScatInHall->getPsi()/deg<<endl;
		//to compare AB with BA	
		G4RotationMatrix pRot_BA=(*pRotBField) * (*pRotX90deg);
		cout<<"\n (*pRotBField) * (*pRotX90deg) ==> EulerAngles: phi="<<pRot_BA.getPhi()/deg
			<<"  theta="<<pRot_BA.getTheta()/deg
			<<"  psi="<<pRot_BA.getPsi()/deg<<endl;

		cout<<"\n pRotScatInHall ==> EulerAngles: phi="<<pRotScatInHall->getPhi()/deg
			<<"  theta="<<pRotScatInHall->getTheta()/deg
			<<"  psi="<<pRotScatInHall->getPsi()/deg<<endl;
	}
#endif
	/////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////
	//scattering chamber container
	/////////////////////////////////////////////////////////
	//This is just a container to enclose the taraget chamber, 10 mm larger than the 
	//scattering chamber itself.
	//The scattering chamber containner is made of helium gas, 
	//will put vacuum inside using innerScatChamber
	//double pScatChamberContainerR=(479.425+0.508+40+10.0)*mm, mScatChamberL=53.5*inch;
	double pSCEntranceWindowLongFlangeL=40.0*mm;		//the thickness of the flange
	double pScatChamberContainerRin=mShieldLHeRin;      //double mShieldLHeRin=38.1*mm
	double pScatChamberContainerRout=mScatChamberRout+mScatChamberExitWindowThick+
		pSCEntranceWindowLongFlangeL+10.0*mm;
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
		mMaterialManager->heliumGas,"scatChamberContainerLogical",0,0,0);
	G4VPhysicalVolume* scatChamberContainerPhys=new G4PVPlacement(pRotScatInHall,
		G4ThreeVector(mScatChamberXOffset,mScatChamberYOffset,mScatChamberZOffset),
		scatChamberContainerLogical,"scatChamberContainerPhys",motherLogical,0,0);
	scatChamberContainerLogical->SetVisAttributes(HallVisAtt); 

#ifdef G4DEBUG_GEOMETRY
	if(G4DEBUG_GEOMETRY==1)
	{	
		//debug only
		new G4PVPlacement(pRotX90Z270deg,G4ThreeVector(),scatChamberContainerLogical,
			"scatChamberContainerPhys_Tran",motherLogical,0,0);
	}
	else if(G4DEBUG_GEOMETRY==2)
	{	
		new G4PVPlacement(pRotX90deg,G4ThreeVector(),scatChamberContainerLogical,
			"scatChamberContainerPhys_Long",motherLogical,0,0);
	}
	else if(G4DEBUG_GEOMETRY==3)
	{	
		new G4PVPlacement(pRotX90Z353deg,G4ThreeVector(),scatChamberContainerLogical,
			"scatChamberContainerPhys_Gep",motherLogical,0,0);
	}
#endif


	/////////////////////////
	// scatter chamber 
	/////////////////////////
	//Build the scatter chamber,it contains 5 windows,  4 rectangles and 1 circle
	//target chamber container || scattering chamber container
	//double mScatChamberRin=17.875*inch,mScatChamberRout=18.875*inch,mScatChamberL=27.25*inch;

	G4VSolid* scatChamberWholeSolid=0;
	//If mSetupG2PScatChamber==1, setup the body only, 
	//If mSetupG2PScatChamber==2, setup the body plus top flange and bottom flange, this
	//will make the program slower
	if(mSetupG2PScatChamber==1)
	{
		scatChamberWholeSolid = new G4Tubs("scatChamberWholeTubs",
			mScatChamberRin,mScatChamberRout,mScatChamberL/2.0,0.,360.*deg);
	}
	else if(mSetupG2PScatChamber>=2)
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
	else if(mSetupG2PScatChamber>=30)
	{
		//build the body, then union with top and bottom flanges
		//this is the same as above
		startphi=0.*deg; deltaphi=360.*deg;
		G4VSolid* SCBodySolid = new G4Tubs("SCBodyTubs",
			mScatChamberRin,mScatChamberRout,mScatChamberL/2.0,startphi,deltaphi);

		const int kNPlane_SCBottom=4;
		double rInner_SCBottom[] = {0,0,mScatChamberRin,mScatChamberRin};
		double rOuter_SCBottom[] = {mScatChamberRout+1.0*inch,mScatChamberRout+1.0*inch,
			mScatChamberRout+1.0*inch,mScatChamberRout};
		double zPlane_SCBottom[] = {-mScatChamberL/2-4.50*inch,-mScatChamberL/2-3.25*inch,
			-mScatChamberL/2-1.0*inch,-mScatChamberL/2};
		G4Polycone* SCBottomFlangeSolid = new G4Polycone("SCBottomFlangePolycone",startphi,deltaphi,
			kNPlane_SCBottom,zPlane_SCBottom,rInner_SCBottom,rOuter_SCBottom);

		const int kNPlane_SCTop=8;
		double rInner_SCTop[] = {
			mScatChamberRin,mScatChamberRin,mScatChamberRin,
			mScatChamberRin,mScatChamberRin,mScatChamberRin,
			0,0
		};
		double rOuter_SCTop[] = {
			mScatChamberRout,mScatChamberRout,mScatChamberRout+1.0*inch,
			mScatChamberRout+1.0*inch,mScatChamberRout,mScatChamberRout,
			mScatChamberRout+1*inch,mScatChamberRout+1*inch 
		};
		double zPlane_SCTop[] = { 
			mScatChamberL/2,mScatChamberL/2+0.25*inch,mScatChamberL/2+1.25*inch,
			mScatChamberL/2+3.50*inch,mScatChamberL/2+3.50*inch,mScatChamberL/2+20.5*inch,
			mScatChamberL/2+20.5*inch,mScatChamberL/2+21.75*inch
		};

		G4Polycone* SCTopFlangeSolid = new G4Polycone("SCTopFlangePolycone",startphi,deltaphi,
			kNPlane_SCTop,zPlane_SCTop,rInner_SCTop,rOuter_SCTop);

		G4UnionSolid* SCBodyUnionBSolid = new G4UnionSolid("SCBodyUnionBSolid",
			SCBodySolid,SCBottomFlangeSolid);
		G4UnionSolid* SCBodyUnionBTSolid = new G4UnionSolid("SCBodyUnionBTSolid",
			SCBodyUnionBSolid,SCTopFlangeSolid);

		scatChamberWholeSolid = SCBodyUnionBTSolid;
	}


	//these are for the subtraction part, not the scatter chamber itself
	double pSCExitWindowH=15.0*inch;
	double pSCEntranceWindowH=6.44*inch;
	double pSCEntranceWindowVoffset=-1.766*inch;
	double pSCWindowRin=mShieldLN2Rin-1.0*cm;
	double pSCWindowRout=mScatChamberRout+1.0*cm;

	//In a bollean solid, if SolidA=SolidB-SolidC, it require solidC thicker than solidB 
	//the following angles are listed in the drawing
	//EntranceWindowTran (phi=72+36) for transverse, rec
	//EntranceWindowLong (d=5.75',phi=180) for logitidinal,  circle
	//EntranceWindowGEP (phi=156.45+7.1) for 20 degree GEP, rec
	//ExitWindowTran (phi=252+36), ExitWindowLong(phi=309+60), ExitWindowLeftSide (phi=19+32), 
	//Viewwindow1(d=3.75'), viewwindow2(d=3.75') are not included 

	//the following angles come from the angles in the drawing subtract 90 degrees, this is 
	//the G4 Lab convention
	//EntranceWindowTran (phi=342+36) for transverse, reduce to 7 degree for g2p 
	//EntranceWindowLong (d=5.75',phi=90) for logitidinal,  circle
	//EntranceWindowGEP (phi=66.5+7) for 20 degree GEP, rec
	//ExitWindowTran (phi=162+36), ExitWindowLong(phi=219+60), ExitWindowLeftSide (phi=289+32), 
	//Viewwindow1(d=3.75'), viewwindow2(d=3.75') are not included 

	//rectangle ExitWindowTran
	startphi=162.0*deg;deltaphi=36*deg;
	G4VSolid* SCExitWindowTranSolid = new G4Tubs("SCExitWindowTranTubs",
		pSCWindowRin,pSCWindowRout,pSCExitWindowH/2.0,startphi,deltaphi);

	//rectangle ExitWindowLong
	startphi=219.0*deg;deltaphi=60*deg;
	G4VSolid* SCExitWindowLongSolid = new G4Tubs("SCExitWindowLongTubs",
		pSCWindowRin,pSCWindowRout,pSCExitWindowH/2.0,startphi,deltaphi);

	//rectangle ExitWindowLeftSide
	startphi=289.0*deg;deltaphi=32*deg;
	G4VSolid* SCExitWindowLeftSideSolid = new G4Tubs("SCExitWindowLeftSideTubs",
		pSCWindowRin,pSCWindowRout,pSCExitWindowH/2.0,startphi,deltaphi);


	//////////////////////////////////////////////////////////////////////////
	//rectangle EntranceWindowTran
	//rectangle EntranceWindowTran, for SANE transverse, LN2 shielding will subtract this solid
	startphi=342.0*deg;deltaphi=36*deg;
	G4VSolid* SCEntranceWindowTranSANESolid = new G4Tubs("SCEntranceWindowTranSANETubs",
		pSCWindowRin,pSCWindowRout,pSCExitWindowH/2.0,startphi,deltaphi);

	//rectangle EntranceWindowTran, for G2P transverse
	//In G2P, Al patch this original window and dig a small window on the patch to use it as
	//the transverse entrance. The size of this entrance window is identical to the GEP entrace
	//in this program, I dig this small window from the scattering chamber
	//But I still dig the big window in LN2 shieding since Al did not patch it
	startphi=356.25*deg;deltaphi=7.5*deg;
	G4VSolid* SCEntranceWindowTranSolid = new G4Tubs("SCEntranceWindowTranTubs",
		pSCWindowRin,pSCWindowRout,pSCEntranceWindowH/2.0,startphi,deltaphi);
	//////////////////////////////////////////////////////////////////////////

	//circle, located at phi=180 degree, EntranceWindowLong subtracted part
	double pSCEntranceWindowLongR=5.75*inch/2.0;
	double pSCEntranceWindowLongL_Half=pSCWindowRout-pSCWindowRin;
	double pSCEntranceWindowLongSubY=(pSCWindowRout+pSCWindowRin)/2.0;
	startphi=0.0*deg;deltaphi=360.0*deg;
	G4VSolid* SCEntranceWindowLongSubSolid = new G4Tubs("SCEntranceWindowLongSubTubs",
		0.0,pSCEntranceWindowLongR,pSCEntranceWindowLongL_Half,startphi,deltaphi);

	//circle, located at phi=180 degree, EntranceWindowLong unioned part, 
	//the flange, outer diameter is 8 inch, opening is 5.6 inch,
	double pSCEntranceWindowLongFlangeRout=8.0*inch/2;
	double pSCEntranceWindowLongFlangeRin=5.6*inch/2;
	//double pSCEntranceWindowLongFlangeL=40.0*mm;  //the thickness of the flange, has been delared 
	//this flange attached to the openning, started position is
	//sqrt(SCRout^2 - Ropen^2) 
	double pSCEntranceWindowLongFlangeY = sqrt(pow(mScatChamberRout,2)-pow(pSCEntranceWindowLongR,2)) + 
		pSCEntranceWindowLongFlangeL/2.0;
	startphi=0.0*deg;deltaphi=360.0*deg;
	G4VSolid* SCEntranceWindowLongFlangeSolid = new G4Tubs("SCEntranceWindowLongFlangeTubs",
		pSCEntranceWindowLongFlangeRin,pSCEntranceWindowLongFlangeRout,
		pSCEntranceWindowLongFlangeL/2.0,startphi,deltaphi);

	//rectangle EntranceWindowGEP 
	startphi=66.25*deg;deltaphi=7.5*deg;
	G4VSolid* SCEntranceWindowGEPSolid = new G4Tubs("SCEntranceWindowGEPTubs",
		pSCWindowRin,pSCWindowRout,pSCEntranceWindowH/2.0,startphi,deltaphi);

	//////////////////////////////////////////////////////////////////////
	//build the boolean geometry

	// subtract the transverse exit 
	G4SubtractionSolid* SCSubtractExitTSolid=new G4SubtractionSolid(
		"SCSubtractExitT",scatChamberWholeSolid,SCExitWindowTranSolid);

	// subtract the longitudinal exit 
	G4SubtractionSolid* SCSubtractExitTNLSolid=new G4SubtractionSolid(
		"SCSubtractExitTNL",SCSubtractExitTSolid,SCExitWindowLongSolid);

	// subtract the left side exit, the orignal size 
	G4SubtractionSolid* SCSubtractExitTNLNSSolid=new G4SubtractionSolid(
		"SCSubtractExitTNLNS",SCSubtractExitTNLSolid,SCExitWindowLeftSideSolid);

	// subtract the transverse entrance 
	G4SubtractionSolid* SCSubtractExitNEntranceTSolid=new G4SubtractionSolid(
		"SCSubtractExitNEntranceT",SCSubtractExitTNLNSSolid,SCEntranceWindowTranSolid,
		0,G4ThreeVector(0,0,pSCEntranceWindowVoffset));

	// subtract the logitudinal entrance, which is a circle 
	//need to rotate about X (in container Coor, in Hall Coor it is Y) by 90 degrees then 
	//subtract it to make a circle hole
	G4SubtractionSolid* SCSubtractExitNEntranceTNLSolid=new G4SubtractionSolid(
		"SCSubtractExitNEntranceTNL",SCSubtractExitNEntranceTSolid,
		SCEntranceWindowLongSubSolid,pRotX90deg,G4ThreeVector(0,pSCEntranceWindowLongSubY,0));

	//subtract GEP entrance windows, 1.766 inch below
	G4SubtractionSolid* SCSubtractExitNEntranceTNLNGSolid=new G4SubtractionSolid(
		"SCSubtractExitNEntranceTNLNG",SCSubtractExitNEntranceTNLSolid,
		SCEntranceWindowGEPSolid,0,G4ThreeVector(0,0,pSCEntranceWindowVoffset));

	//union with the flange(the protrude part) of logitudinal entrance
	G4UnionSolid *SCSubtractExitNEntranceTNLNGUnionLSolid=new G4UnionSolid(
		"SCSubtractExitNEntranceTNLNGUnionL",SCSubtractExitNEntranceTNLNGSolid,
		SCEntranceWindowLongFlangeSolid,pRotX90deg,G4ThreeVector(0,pSCEntranceWindowLongFlangeY,0));

	//setup the scatter chamber
	G4LogicalVolume* scatChamberLogical = new G4LogicalVolume(
		SCSubtractExitNEntranceTNLNGUnionLSolid,mMaterialManager->aluminum,"scatChamberLogical",0,0,0);
	scatChamberLogical->SetVisAttributes(WhiteVisAtt); 

	new G4PVPlacement(0,G4ThreeVector(),scatChamberLogical,"scatChamberPhys",
		scatChamberContainerLogical,0,0);

	/////////////////////////
	// target chamber window covers 
	/////////////////////////
	//Covers for EntranceWindow Tran,Long,Side ,ExitWindowTran, ExitWindowLong,  
	//ExitWindowLeftSide, viewwindow1(d=3.75'), viewwindow2(d=3.75') are not included 
	/////////////////////////////////////////////////////////////////////
	//Al told that the exit window for 0, 20 and 90 degrees target field are all 20 mil of aluminum
	//double mScatChamberExitWindowThick=0.02*inch;
	double pSCExitWindowCoverH=pSCExitWindowH+0.8*inch;
	double pSCEntranceWindowCoverH=pSCEntranceWindowH+0.8*inch;
	double pSCExitWindowCoverRin=mScatChamberRout;
	double pSCExitWindowCoverRout=mScatChamberRout+mScatChamberExitWindowThick;

	//ExitWindowTranCover
	startphi=162.0*deg;deltaphi=36*deg;
	G4VSolid* SCExitWindowTranCoverSolid = new G4Tubs("SCExitWindowTranTubs",
		pSCExitWindowCoverRin,pSCExitWindowCoverRout,pSCExitWindowCoverH/2.0,startphi,deltaphi);
	G4LogicalVolume* SCExitWindowTranCoverLogical = new G4LogicalVolume(
		SCExitWindowTranCoverSolid,mMaterialManager->aluminum,"SCExitWindowTranCoverLogical",0,0,0);
	new G4PVPlacement(0,G4ThreeVector(),SCExitWindowTranCoverLogical,
		"SCExitWindowTranCoverPhys",scatChamberContainerLogical,false,0);
	SCExitWindowTranCoverLogical->SetVisAttributes(LightYellowVisAtt); 

	//ExitWindowLongCover
	startphi=219.0*deg;deltaphi=60*deg;
	G4VSolid* SCExitWindowLongCoverSolid = new G4Tubs("SCExitWindowLongTubs",
		pSCExitWindowCoverRin,pSCExitWindowCoverRout,pSCExitWindowCoverH/2.0,startphi,deltaphi);
	G4LogicalVolume* SCExitWindowLongCoverLogical = new G4LogicalVolume(
		SCExitWindowLongCoverSolid,mMaterialManager->aluminum,"SCExitWindowLongCoverLogical",0,0,0);
	SCExitWindowLongCoverLogical->SetVisAttributes(LightYellowVisAtt); 

	new G4PVPlacement(0,G4ThreeVector(),SCExitWindowLongCoverLogical,
		"SCExitWindowLongCoverPhys",scatChamberContainerLogical,false,0);

	//ExitWindowLeftSideCover
	startphi=289.0*deg;deltaphi=32*deg;
	G4VSolid* SCExitWindowLeftSideCoverSolid = new G4Tubs("SCExitWindowLeftSideTubs",
		pSCExitWindowCoverRin,pSCExitWindowCoverRout,pSCExitWindowCoverH/2.0,startphi,deltaphi);
	G4LogicalVolume* SCExitWindowLeftSideCoverLogical = new G4LogicalVolume(
		SCExitWindowLeftSideCoverSolid,mMaterialManager->aluminum,
		"SCExitWindowLeftSideCoverLogical",0,0,0);
	SCExitWindowLeftSideCoverLogical->SetVisAttributes(LightYellowVisAtt); 

	new G4PVPlacement(0,G4ThreeVector(),SCExitWindowLeftSideCoverLogical,
		"SCExitWindowLeftSideCoverPhys",scatChamberContainerLogical,false,0);

	/////////////////////////////////////////////////////////////////////
	//Al told that the entrance window for 20 and 90 degrees target field are 7 mil of aluminum
	double pSCEntranceWindowCoverRin=mScatChamberRout;
	double pSCEntranceWindowCoverRout=mScatChamberRout+0.007*inch;

	//EntranceWindowTranCover, G2P window
	startphi=356.25*deg;deltaphi=7.5*deg;
	G4VSolid* SCEntranceWindowTranCoverSolid = new G4Tubs("SCEntranceWindowTranCoverTubs",
		pSCEntranceWindowCoverRin,pSCEntranceWindowCoverRout,
		pSCEntranceWindowCoverH/2.0,startphi,deltaphi);

	G4LogicalVolume* SCEntranceWindowTranCoverLogical = new G4LogicalVolume(
		SCEntranceWindowTranCoverSolid,mMaterialManager->aluminum,
		"SCEntranceWindowTranCoverLogical",0,0,0);
	SCEntranceWindowTranCoverLogical->SetVisAttributes(LightYellowVisAtt); 

	new G4PVPlacement(0,G4ThreeVector(0,0,pSCEntranceWindowVoffset),
		SCEntranceWindowTranCoverLogical,"SCEntranceWindowTranCoverPhys",
		scatChamberContainerLogical,false,0);


	////EntranceWindowLongCover, this rotation is very strange, but it works!!!
	//!!!Note, only work if the BField rotation is about Y axis
	////circle, located at 180 phi=degree, this entrance window is for 0 degree target field
	//Al told me that the material is 20 mil Beryllium
	double pBeThick=0.02*inch;
	double pSCEntranceWindowLongCoverY=pSCEntranceWindowLongFlangeY
		+pSCEntranceWindowLongFlangeL/2+pBeThick/2.0;
	startphi=0.0*deg;deltaphi=360.0*deg;
	G4VSolid* SCEntranceWindowLongCoverSolid = new G4Tubs("SCEntranceWindowLongCoverTubs",
		0.0,pSCEntranceWindowLongFlangeRout,pBeThick/2.0,startphi,deltaphi);
	G4LogicalVolume* SCEntranceWindowLongCoverLogical = new G4LogicalVolume(
		SCEntranceWindowLongCoverSolid,mMaterialManager->beryllium,
		"SCEntranceWindowLongCoverLogical",0,0,0);
	SCEntranceWindowLongCoverLogical->SetVisAttributes(LightYellowVisAtt);

	G4RotationMatrix* pRotEntranceWindowLongCover=new G4RotationMatrix();
	*pRotEntranceWindowLongCover=(*pRotBField)*(*pRotBField)*(*pRotX90deg);
	/////////////////////////////////////////////////////////
#ifdef G4DEBUG_GEOMETRY
	if(G4DEBUG_GEOMETRY>=2)
	{
		cout<<"\n (*pRotBField)*(*pRotBField)*(*pRotX90deg) ==> EulerAngles: phi="
			<<pRotEntranceWindowLongCover->getPhi()/deg
			<<"  theta="<<pRotEntranceWindowLongCover->getTheta()/deg
			<<"  psi="<<pRotEntranceWindowLongCover->getPsi()/deg<<endl;

		G4RotationMatrix* pRotEn2Cover=new G4RotationMatrix();
		*pRotEn2Cover=(*pRotX270deg)*(*pRotBField)*(*pRotBField);
		cout<<"\n (*pRotX270deg)*(*pRotBField)*(*pRotBField) ==> EulerAngles: phi="
			<<pRotEn2Cover->getPhi()/deg
			<<"  theta="<<pRotEn2Cover->getTheta()/deg
			<<"  psi="<<pRotEn2Cover->getPsi()/deg<<endl;

		*pRotEn2Cover=(*pRotY90deg)*(*pRotBField);
		cout<<"\n (*pRotY90deg)*(*pRotBField) ==> EulerAngles: phi="<<pRotEn2Cover->getPhi()/deg
			<<"  theta="<<pRotEn2Cover->getTheta()/deg
			<<"  psi="<<pRotEn2Cover->getPsi()/deg<<endl;
	}
#endif
	/////////////////////////////////////////////////////////
	new G4PVPlacement(pRotEntranceWindowLongCover,G4ThreeVector(0,pSCEntranceWindowLongCoverY,0),
		SCEntranceWindowLongCoverLogical,"SCEntranceWindowLongCoverPhys",
		scatChamberContainerLogical,false,0); 

	//EntranceWindowGEPCover, for GEP
	startphi=66.25*deg;deltaphi=7.5*deg;
	G4VSolid* SCEntranceWindowGEPCoverSolid = new G4Tubs("SCEntranceWindowGEPCoverTubs",
		pSCEntranceWindowCoverRin,pSCEntranceWindowCoverRout,
		pSCEntranceWindowCoverH/2.0,startphi,deltaphi);
	G4LogicalVolume* SCEntranceWindowGEPCoverLogical = new G4LogicalVolume(
		SCEntranceWindowGEPCoverSolid,mMaterialManager->aluminum,
		"SCEntranceWindowGEPCoverLogical",0,0,0);
	SCEntranceWindowGEPCoverLogical->SetVisAttributes(LightYellowVisAtt); 

	new G4PVPlacement(0,G4ThreeVector(0,0,pSCEntranceWindowVoffset),SCEntranceWindowGEPCoverLogical,
		"SCEntranceWindowGEPCoverPhys",scatChamberContainerLogical,false,0);
	/////////////////////////////////////////////////////////////////////

	/////////////////////////
	// innerOfChamber 
	/////////////////////////
	//inside the scatter chamber, it is vacuum, all stuff inside the chamber will use 
	//the local coordination of this mother volumn, 
	startphi=0.0*deg;deltaphi=360.0*deg;
	G4VSolid* innerOfChamberSolid = new G4Tubs("innerOfChamberTubs",mShieldLHeRin,mScatChamberRin,
		mScatChamberL/2.0,startphi,deltaphi);
	G4LogicalVolume* innerOfChamberLogical = new G4LogicalVolume(innerOfChamberSolid,
		mMaterialManager->vacuum,"innerOfChamberLogical",0,0,0);
	//By Jixie: Add this step limit can help in calculating the integrated BdL
	double pScatChamberStepLimit=10;
	gConfig->GetArgument("ScatChamberStepLimit",pScatChamberStepLimit);
	pScatChamberStepLimit*=mm;
	G4UserLimits* uSCStepLimits = new G4UserLimits(pScatChamberStepLimit);
	innerOfChamberLogical->SetUserLimits(uSCStepLimits);
	innerOfChamberLogical->SetVisAttributes(HallVisAtt);  //white invisible

	new G4PVPlacement(0,G4ThreeVector(),innerOfChamberLogical,"innerOfChamberPhys",
		scatChamberContainerLogical,0,0);


	/////////////////////////
	// LN2 shielding 
	/////////////////////////

	//LN2 Shield cylinder itself is 1/4 inch, but have windows similar to the scattering chamber
	//there are also covers for each opening. In order to simplify the geometry, I union a
	//thin aluminum foil (1.5 mil) onto the outer surface of it as the cover
	//LN2 Shield cynlider with Din=32.5", Dout=33" and L=33.25", the cover is 0.0015'=0.0381mm,
	//Drawing# 67504-E-0015 sheet 2
	//double mShieldLN2Rin=16.15*inch,mShieldLN2Rout=16.50*inch, mShieldLN2L=33.25*inch;
	startphi=0.0*deg;deltaphi=360.0*deg;
	G4VSolid* shieldLN2WholeSolid = new G4Tubs("shieldLN2WholeTubs",mShieldLN2Rin,mShieldLN2Rout,
		mShieldLN2L/2.0,startphi,deltaphi);

	if(mShieldLN2WindowThick/mm < 1.0E-05) mShieldLN2WindowThick=1.0E-06*mm;
	G4VSolid* shieldLN2CoverSolid = new G4Tubs("shieldLN2CoverTubs",mShieldLN2Rout,
		mShieldLN2Rout+mShieldLN2WindowThick,mShieldLN2L/2.0,startphi,deltaphi);

	//////////////////////////////////////////////////////////////////////
	//build the boolean geometry

	// subtract the transverse exit 
	G4SubtractionSolid* ShieldLN2SubtractExitTSolid=new G4SubtractionSolid(
		"ShieldLN2SubtractExitT",shieldLN2WholeSolid,SCExitWindowTranSolid);

	// subtract the longitudinal exit 
	G4SubtractionSolid* ShieldLN2SubtractExitTNLSolid=new G4SubtractionSolid(
		"ShieldLN2SubtractExitTNL",ShieldLN2SubtractExitTSolid,SCExitWindowLongSolid);

	// subtract the left side exit 
	G4SubtractionSolid* ShieldLN2SubtractExitTNLNSSolid=new G4SubtractionSolid(
		"ShieldLN2SubtractExitTNLNS",ShieldLN2SubtractExitTNLSolid,SCExitWindowLeftSideSolid);

	// subtract the SANE transverse entrance 
	G4SubtractionSolid* ShieldLN2SubtractExitNEntranceTSolid=new G4SubtractionSolid(
		"ShieldLN2SubtractExitNEntranceT",
		ShieldLN2SubtractExitTNLNSSolid,SCEntranceWindowTranSANESolid);

	// subtract the logitudinal entrance, which is a circle 
	//need to rotate about X (in container Coor, in Hall Coor it is Y) by 90 degrees then 
	//subtract it to make a circle hole
	G4SubtractionSolid* ShieldLN2SubtractExitNEntranceTNLSolid=new G4SubtractionSolid(
		"ShieldLN2SubtractExitNEntranceTNL",ShieldLN2SubtractExitNEntranceTSolid,
		SCEntranceWindowLongSubSolid,pRotX90deg,G4ThreeVector(0,pSCEntranceWindowLongSubY,0));

	//subtract GEP entrance windows, 1.766 inch below
	// Al make it a circle, since GEP angle is not finalize yet, I do not want to change it now
	// need to change it later
	G4SubtractionSolid* ShieldLN2SubtractExitNEntranceTNLNGSolid=new G4SubtractionSolid(
		"ShieldLN2SubtractExitNEntranceTNLNG",ShieldLN2SubtractExitNEntranceTNLSolid,
		SCEntranceWindowGEPSolid,0,G4ThreeVector(0,0,pSCEntranceWindowVoffset));

	//Union with the cover 
	G4VSolid* ShieldLN2Solid=(G4VSolid*)ShieldLN2SubtractExitNEntranceTNLNGSolid;
	if(mShieldLN2WindowThick/mm > 1.0E-05)
	{
		ShieldLN2Solid=new G4UnionSolid("shieldLN2Soild",
			ShieldLN2SubtractExitNEntranceTNLNGSolid,shieldLN2CoverSolid);
	}
	/////////////////////////////////////////////////////////////////////
	//now build and place the geometry

	G4LogicalVolume* shieldLN2Logical = new G4LogicalVolume(ShieldLN2Solid,
		mMaterialManager->aluminum,"shieldLN2Logical",0,0,0);
	shieldLN2Logical->SetVisAttributes(WhiteVisAtt);  //white

	new G4PVPlacement(0,G4ThreeVector(),shieldLN2Logical,"shieldLN2Phys",
		innerOfChamberLogical,0,0);


	/////////////////////////
	// target field coils 
	/////////////////////////
	//
	if (mSetupCoil==1)
	{
		// just 2 helm coils, downstream and upstream, connection part are not included
		//Downtream and upstream coil
		//this number here is for Bz along Z axis
		startphi=0.*deg; deltaphi=360.*deg;
		const int kNPlane_Coil=6;
		const double kInnerLineTheta=39.65*deg; //The inner surface is alone a line with this angle
		double rInner_Coil[] = {4.04*inch,4.04*inch,3.82*inch/tan(kInnerLineTheta),
			11.60*inch,12.06*inch,12.06*inch};
		double rOuter_Coil[] = {5.14*inch,11.60*inch,12.56*inch,12.56*inch,12.56*inch,12.56*inch};
		double zPlane_DownCoil[] = {1.57*inch,3.37*inch,3.82*inch,9.61*inch,9.611*inch,10.21*inch};
		double zPlane_UpCoil[kNPlane_Coil];
		for(int ii=0;ii<kNPlane_Coil;ii++) zPlane_UpCoil[ii]=-1*zPlane_DownCoil[ii];

		G4VSolid* downCoilSolid = new G4Polycone("DownCoilPcon",startphi,deltaphi,
			kNPlane_Coil,zPlane_DownCoil,rInner_Coil,rOuter_Coil);
		G4VSolid* upCoilSolid = new G4Polycone("UpCoilPcon",startphi,deltaphi,
			kNPlane_Coil,zPlane_UpCoil,rInner_Coil,rOuter_Coil);
		G4LogicalVolume* downCoilLogical = new G4LogicalVolume(downCoilSolid,
			mMaterialManager->stainlesssteel,"downCoilLogical",0,0,0);
		downCoilLogical->SetVisAttributes(DarkBlueVisAtt);
		G4LogicalVolume* upCoilLogical = new G4LogicalVolume(upCoilSolid,
			mMaterialManager->stainlesssteel,"upCoilLogical",0,0,0);
		upCoilLogical->SetVisAttributes(DarkBlueVisAtt);  //dark blue

		new G4PVPlacement(pRotX270deg,G4ThreeVector(),downCoilLogical,
			"DownCoilPhys",innerOfChamberLogical,false,0);
		new G4PVPlacement(pRotX270deg,G4ThreeVector(),upCoilLogical,
			"UpCoilPhys",innerOfChamberLogical,false,0);

	}
	else if (mSetupCoil>=2)
	{
		//Hall B Coil
		startphi=0.*deg;deltaphi=360.*deg;
		G4RotationMatrix* pRotDCoil = new G4RotationMatrix();
		pRotDCoil->rotateY(90.*deg);
		pRotDCoil->rotateX(-90.*deg);
		G4RotationMatrix* pRotUCoil = new G4RotationMatrix();
		pRotUCoil->rotateY(90.*deg);
		pRotUCoil->rotateX(90.*deg);

		//Coils
		const int kNTubs=2;
		double rInnerTubs[] = {116.5*mm,116.5*mm};
		double rOuterTubs[] = {153.0*mm,153.0*mm};
		double zPlaneTubs[] = {-16.0*mm,16.0*mm};
		G4Polycone* InnerCoilSolid = new G4Polycone("InnerCoilSolid",startphi,deltaphi,
			kNTubs,zPlaneTubs,rInnerTubs,rOuterTubs);
		//By Jixie: The tub is not shown correctly in the visulization, Chao use a polycone solid instead
		//This is a bug of HEPREP viewer itself
		//G4Tubs* InnerCoilSolid = new G4Tubs("InnerCoilSolid",116.5*mm,153.0*mm,32.0*mm/2.0,startphi,deltaphi);
		G4LogicalVolume* InnerCoilLogical = new G4LogicalVolume(InnerCoilSolid,
			mMaterialManager->copper,"InnerCoilLogical",0,0,0);
		InnerCoilLogical->SetVisAttributes(CuBrownVisAtt);
		new G4PVPlacement(pRotDCoil,G4ThreeVector(0,-68.0*mm,0),InnerCoilLogical,
			"InnerCoilPhys",innerOfChamberLogical,true,0,0);
		new G4PVPlacement(pRotUCoil,G4ThreeVector(0,68.0*mm,0),InnerCoilLogical,
			"InnerCoilPhys",innerOfChamberLogical,true,1,0);
		rInnerTubs[0]=rInnerTubs[1]=165.5*mm;
		rOuterTubs[0]=rOuterTubs[1]=222.5*mm;
		zPlaneTubs[0]=-24.0*mm;zPlaneTubs[1]=24.0*mm;
		G4Polycone* OuterCoilSolid = new G4Polycone("OuterCoilSolid",startphi,deltaphi,
			kNTubs,zPlaneTubs,rInnerTubs,rOuterTubs);
		//By Jixie: The tub is not shown correctly in the visulization, chao use a polycone solid instead        
		//G4Tubs* OuterCoilSolid = new G4Tubs("OuterCoilSolid",165.5*mm,222.5*mm,48.0*mm/2.0,startphi,deltaphi);
		G4LogicalVolume* OuterCoilLogical = new G4LogicalVolume(OuterCoilSolid,
			mMaterialManager->copper,"OuterCoilLogical",0,0,0);
		OuterCoilLogical->SetVisAttributes(CuBrownVisAtt);
		new G4PVPlacement(pRotDCoil,G4ThreeVector(0,-99.0*mm,0),OuterCoilLogical,
			"OuterCoilPhys",innerOfChamberLogical,true,2,0);
		new G4PVPlacement(pRotUCoil,G4ThreeVector(0,99.0*mm,0),OuterCoilLogical,
			"OuterCoilPhys",innerOfChamberLogical,true,3,0);

		//Coil Former
		const int kNPlane_Coil_F=9;
		double rInner_Coil_F[] = {105.0*mm,105.0*mm,153.0*mm,
			153.0*mm,153.0*mm,153.0*mm,153.0*mm,222.5*mm,222.5*mm};
		double rOuter_Coil_F[] = {130.0*mm,130.0+34.5/tan(17.0*deg)*mm,130.0+34.5/tan(17.0*deg)*mm,
			250.0*mm,250.0*mm,244.5*mm,244.5*mm,244.5*mm,244.5*mm};
		double zPlane_Coil_F[] = {0.0*mm,34.5*mm,34.5*mm,
			120*tan(17.0*deg)*mm,62.0*mm,62.0*mm,74.0*mm,74.0*mm,76.0*mm};
		G4Polycone* CoilFPolycone = new G4Polycone("CoilFPcon",startphi,deltaphi,
			kNPlane_Coil_F,zPlane_Coil_F,rInner_Coil_F,rOuter_Coil_F);
		G4SubtractionSolid* CoilFSub1 = new G4SubtractionSolid("CoilFSub1",
			CoilFPolycone,InnerCoilSolid,0,G4ThreeVector(0,0,27.5*mm));
		G4SubtractionSolid* CoilFSolid = new G4SubtractionSolid("CoilFSolid",
			CoilFSub1,OuterCoilSolid,0,G4ThreeVector(0,0,58.5*mm));
		G4LogicalVolume* CoilFLogical = new G4LogicalVolume(CoilFSolid,
			mMaterialManager->stainlesssteel,"CoilFLogical",0,0,0);
		CoilFLogical->SetVisAttributes(DarkBlueVisAtt);
		new G4PVPlacement(pRotDCoil,G4ThreeVector(0,-40.5*mm,0),
			CoilFLogical,"CoilFormerPhys",innerOfChamberLogical,true,0,0);
		new G4PVPlacement(pRotUCoil,G4ThreeVector(0,40.5*mm,0),
			CoilFLogical,"CoilFormerPhys",innerOfChamberLogical,true,1,0);

		//Coil Retainer
		const int kNPlane_Coil_R=7;
		double rInner_Coil_R[] = {105.0*mm,105.0*mm,105.0*mm,105.0*mm,
			(150.0-9.0*tan(50.0*deg))*mm,(150.0-9.0*tan(50.0*deg))*mm,150.0*mm};
		double rOuter_Coil_R[] = {116.5*mm,116.5*mm,153.0*mm,153.0*mm,
			153.0*mm,165.5*mm,165.5*mm};
		double zPlane_Coil_R[] = {0.0*mm,9.0*mm,9.0*mm,(48.0-45.0/tan(50.0*deg))*mm,
			39.0*mm,39.0*mm,48.0*mm};
		G4Polycone* CoilRSolid = new G4Polycone("CoilRSolid",startphi,deltaphi,
			kNPlane_Coil_R,zPlane_Coil_R,rInner_Coil_R,rOuter_Coil_R);
		G4LogicalVolume* CoilRLogical = new G4LogicalVolume(CoilRSolid,
			mMaterialManager->stainlesssteel,"CoilRLogical",0,0,0);
		CoilRLogical->SetVisAttributes(DarkBlueVisAtt);
		new G4PVPlacement(pRotDCoil,G4ThreeVector(0,-75.0*mm,0),
			CoilRLogical,"CoilRetainerPhys",innerOfChamberLogical,true,0,0);
		new G4PVPlacement(pRotUCoil,G4ThreeVector(0,75.0*mm,0),
			CoilRLogical,"CoilRetainerPhys",innerOfChamberLogical,true,1,0);

		//Joint Flange
		const int kNPlane_Joint_F=5;
		double rInner_Joint_F[] = {170.0*mm,170.0*mm,170.0*mm,170.0*mm,177.0*mm};
		double rOuter_Joint_F[] = {221.0*mm,221.0*mm,240.0*mm,240.0*mm,240.0*mm};
		double zPlane_Joint_F[] = {0.0*mm,6.0*mm,6.0*mm,15.0*mm,22.0*mm};
		G4Polycone* JointFSolid = new G4Polycone("JointFSolid",startphi,deltaphi,
			kNPlane_Joint_F,zPlane_Joint_F,rInner_Joint_F,rOuter_Joint_F);
		G4LogicalVolume* JointFLogical = new G4LogicalVolume(JointFSolid,
			mMaterialManager->plastic,"ClampPLogical",0,0,0);
		JointFLogical->SetVisAttributes(YellowVisAtt);
		new G4PVPlacement(pRotDCoil,G4ThreeVector(0,-123.0*mm,0),
			JointFLogical,"JointFlangePhys",innerOfChamberLogical,true,0,0);
		new G4PVPlacement(pRotUCoil,G4ThreeVector(0,123.0*mm,0),
			JointFLogical,"JointFlangePhys",innerOfChamberLogical,true,1,0);

		//Left Side Surface
		const int kNPlane_Spinning=3;
		double rInner_Spinning[] = {102.0*mm,102.0*mm,241.5*mm};
		double rOuter_Spinning[] = {105.0*mm,105.0*mm,244.5*mm};
		double zPlane_Spinning[] = {40.5*mm,(123.0-45.0/tan(50.0*deg))*mm,205.7*mm};
		G4Polycone* SpinningSolid = new G4Polycone("SpinningSolid",startphi,deltaphi,
			kNPlane_Spinning,zPlane_Spinning,rInner_Spinning,rOuter_Spinning);
		rInnerTubs[0]=rInnerTubs[1]=244.5*mm;
		rOuterTubs[0]=rOuterTubs[1]=250.0*mm;
		zPlaneTubs[0]=-103.2/2.0*mm;zPlaneTubs[1]=103.2/2.0*mm;
		G4Polycone* LCTubSolid = new G4Polycone("LCTubSolid",startphi,deltaphi,
			kNTubs,zPlaneTubs,rInnerTubs,rOuterTubs);
		//G4Tubs* LCTubSolid = new G4Tubs("LCTubSolid",244.5*mm,250.0*mm,103.2/2.0*mm,startphi,deltaphi);
		G4LogicalVolume* SpinningLogical = new G4LogicalVolume(SpinningSolid,
			mMaterialManager->stainlesssteel,"SpinningLogical",0,0,0);
		G4LogicalVolume* LCTubLogical = new G4LogicalVolume(LCTubSolid,
			mMaterialManager->stainlesssteel,"LeftCoverTubLogical",0,0,0);
		SpinningLogical->SetVisAttributes(DarkBlueVisAtt);
		LCTubLogical->SetVisAttributes(DarkBlueVisAtt);
		new G4PVPlacement(pRotDCoil,G4ThreeVector(),
			SpinningLogical,"LeftSurfPhys",innerOfChamberLogical,true,0,0);
		new G4PVPlacement(pRotDCoil,G4ThreeVector(0,-154.1*mm,0),
			LCTubLogical,"LeftSurfPhys",innerOfChamberLogical,true,1,0);

		//Right Side Surface
		G4Tubs* FridgeTubsOuter = new G4Tubs("FridgeTubsOuter",0.0*mm,22.25*mm,500.0/2.0*mm,
			startphi,deltaphi);
		G4RotationMatrix* pRotTub = new G4RotationMatrix();
		pRotTub->rotateY(25.*deg);
		G4Tubs* HeBaseFlgTubs1 = new G4Tubs("HeBaseFlgTubs1",38.1*mm,250.0*mm,19.5/2.0*mm,
			startphi,deltaphi);
		G4Tubs* HeBaseFlgTubs2 = new G4Tubs("HeBaseFlgTubs2",244.5*mm,250.0*mm,99.4/2.0*mm,
			startphi,deltaphi);
		G4UnionSolid* HeBaseFlgU1 = new G4UnionSolid("HeBaseFlgU1",
			HeBaseFlgTubs1,HeBaseFlgTubs2,0,G4ThreeVector(0,0,-39.95*mm));
		G4SubtractionSolid* HeBaseFlgSolid = new G4SubtractionSolid("HeBaseFlgSolid",
			HeBaseFlgU1,FridgeTubsOuter,
			pRotTub,G4ThreeVector(-(192.15+37.0/sin(25.0*deg))*tan(25.0*deg)*mm,0,0));
		G4LogicalVolume *HeBaseFlgLogical = new G4LogicalVolume(HeBaseFlgSolid,
			mMaterialManager->stainlesssteel,"HeBaseFlgLogical",0,0,0);
		HeBaseFlgLogical->SetVisAttributes(DarkBlueVisAtt);
		G4Tubs* HeCanBulkheadTubs1 = new G4Tubs("HeCanBulkheadTubs1",38.1*mm,105.0*mm,9.5/2.0*mm,
			startphi,deltaphi);
		G4Tubs* HeCanBulkheadTubs2 = new G4Tubs("HeCanBulkheadTubs2",102.0*mm,105.0*mm,24.4/2.0*mm,
			startphi,deltaphi);
		G4UnionSolid* HeCanBulkheadU1 = new G4UnionSolid("HeCBU1",
			HeCanBulkheadTubs1,HeCanBulkheadTubs2,0,G4ThreeVector(0,0,-7.45*mm));
		G4SubtractionSolid* HeCBSolid = new G4SubtractionSolid("HeCBSolid",
			HeCanBulkheadU1,FridgeTubsOuter,
			pRotTub,G4ThreeVector(-(60.15+37.0/sin(25.0*deg))*tan(25.0*deg)*mm,0,0));
		G4LogicalVolume* HeCBLogical = new G4LogicalVolume(HeCBSolid,
			mMaterialManager->stainlesssteel,"HeCBLogical",0,0,0);
		HeCBLogical->SetVisAttributes(DarkBlueVisAtt);
		//By Jixie: The tub is not shown correctly in the HEPREP visulization, 
		// Chao use an intersection solid to replace it, this is a bug of HEPREP viewer
		//G4Tubs* BeamTubSolid = new G4Tubs("BeamTubSolid",36.5*mm,36.6*mm,146.5/2.0*mm,startphi,deltaphi);
		G4Tubs* BeamTubs = new G4Tubs("BeamTubs",36.5*mm,36.6*mm,300.0/2.0*mm,startphi,deltaphi);
		G4Box* BeamBox = new G4Box("BeamBox",200.0/2.0*mm,200.0/2.0*mm,146.5/2.0*mm);
		G4IntersectionSolid* BeamTubSolid = new G4IntersectionSolid("BeamTubSolid",BeamBox,BeamTubs);
		G4LogicalVolume* BeamTubLogical = new G4LogicalVolume(BeamTubSolid,
			mMaterialManager->stainlesssteel,"BeamTubLogical",0,0,0);
		BeamTubLogical->SetVisAttributes(DarkBlueVisAtt);
		G4Tubs* FridgeTubs = new G4Tubs("FridgeTubs",21.85*mm,22.25*mm,300.0/2.0*mm,startphi,deltaphi);
		G4Box* FridgeBox = new G4Box("FridgeBox",200.0/2.0*mm,200.0/2.0*mm,146.5/2.0*mm);
		G4IntersectionSolid* FridgeTubSolid = new G4IntersectionSolid("FridgeTubSolid",
			FridgeBox,FridgeTubs,pRotTub,G4ThreeVector());
		G4LogicalVolume* FridgeTubLogical = new G4LogicalVolume(FridgeTubSolid,
			mMaterialManager->stainlesssteel,"FridgeTubLogical",0,0,0);
		FridgeTubLogical->SetVisAttributes(DarkBlueVisAtt);
		new G4PVPlacement(pRotUCoil,G4ThreeVector(0,192.15*mm,0),
			HeBaseFlgLogical,"RightSurfPhys",innerOfChamberLogical,true,0,0);
		new G4PVPlacement(pRotUCoil,G4ThreeVector(0,60.15*mm,0),
			HeCBLogical,"RightSurfPhys",innerOfChamberLogical,true,1,0);
		new G4PVPlacement(pRotUCoil,G4ThreeVector(0,128.65*mm,0),
			BeamTubLogical,"RightSurfPhys",innerOfChamberLogical,true,2,0);
		new G4PVPlacement(pRotUCoil,
			G4ThreeVector(0,128.65*mm,-(128.65+37.0/sin(25.0*deg))*tan(25.0*deg)*mm),
			FridgeTubLogical,"RightSurfPhys",innerOfChamberLogical,true,3,0);

		//Left Side LHe Container
		const int kNPlane_LeftHe=4;
		double rInner_LeftHe[] = {222.5*mm,222.5*mm,150.0*mm,244.5*mm};
		double rOuter_LeftHe[] = {244.5*mm,244.5*mm,244.5*mm,244.5*mm};
		double zPlane_LeftHe[] = {116.5*mm,123.0*mm,123.0*mm,205.7*mm};
		double zPlane_LeftHe_S[] = {-1.0*mm,6.0*mm,6.0*mm,15.0*mm,22.0*mm};
		G4Polycone* LeftHeSPolycone = new G4Polycone("LeftHeSPcon",startphi,deltaphi,
			kNPlane_Joint_F,zPlane_LeftHe_S,rInner_Joint_F,rOuter_Joint_F);
		G4Polycone* LeftHePolycone = new G4Polycone("LeftHePcon",startphi,deltaphi,
			kNPlane_LeftHe,zPlane_LeftHe,rInner_LeftHe,rOuter_LeftHe);
		G4SubtractionSolid* LeftHeSolid = new G4SubtractionSolid("LeftHeSolid",
			LeftHePolycone,LeftHeSPolycone,0,G4ThreeVector(0,0,123.0*mm));
		G4LogicalVolume* LeftHeLogical = new G4LogicalVolume(LeftHeSolid,
			mMaterialManager->liquidHe,"LeftHeLogical",0,0,0);
		LeftHeLogical->SetVisAttributes(WhiteVisAtt);
		new G4PVPlacement(pRotDCoil,G4ThreeVector(),
			LeftHeLogical,"LeftHePhys",innerOfChamberLogical,false,0);

		//Right Side LHe Container
		const int kNPlane_RightHe=5;
		double rInner_RightHe[] = {38.1*mm,38.1*mm,38.1*mm,38.1*mm,38.1*mm};
		double rOuter_RightHe[] = {105.0*mm,105.0*mm,150.0*mm,244.5*mm,244.5*mm};
		double zPlane_RightHe[] = {64.9*mm,(123.0-45.0/tan(50.0*deg))*mm,123.0*mm,123.0*mm,182.4*mm};
		G4Polycone* RightHePolycone = new G4Polycone("RightHePcon",startphi,deltaphi,
			kNPlane_RightHe,zPlane_RightHe,rInner_RightHe,rOuter_RightHe);
		G4SubtractionSolid* RightHeS1 = new G4SubtractionSolid("RightHeS1",
			RightHePolycone,FridgeTubsOuter,pRotTub,G4ThreeVector(-37.0/cos(25.0*deg)*mm,0,0));
		G4SubtractionSolid* RightHeS2 = new G4SubtractionSolid("RightHeS2",
			RightHeS1,LeftHeSPolycone,0,G4ThreeVector(0,0,123.0*mm));
		G4Tubs* RightHeTubs = new G4Tubs("RightHeTubs",222.5*mm,244.5*mm,(6.5+1.0)/2.0*mm,
			startphi,deltaphi);
		G4UnionSolid* RightHeSolid = new G4UnionSolid("RightHeSolid",
			RightHeS2,RightHeTubs,0,G4ThreeVector(0,0,120.25));
		G4LogicalVolume* RightHeLogical = new G4LogicalVolume(RightHeSolid,
			mMaterialManager->liquidHe,"RightHeLogical",0,0,0);
		RightHeLogical->SetVisAttributes(WhiteVisAtt);
		new G4PVPlacement(pRotUCoil,G4ThreeVector(),
			RightHeLogical,"RightHePhys",innerOfChamberLogical,false,0);

		//Large Wedge
		const int kNPlane_LWedge=10;
		double rInner_LWedge[] = {250.0*mm,(130.0+9.5/tan(17.0*deg))*mm,(130.0+9.5/tan(17.0*deg))*mm,
			130.0*mm,112.0*mm,112.0*mm,130.0*mm,(130.0+9.5/tan(17.0*deg))*mm,
			(130.0+9.5/tan(17.0*deg))*mm,250.0*mm};
		double rOuter_LWedge[] = {250.0*mm,250.0*mm,220.0*mm,
			220.0*mm,220.0*mm,220.0*mm,220.0*mm,220.0*mm,
			250.0*mm,250.0*mm};
		double zPlane_LWedge[] = {-(40.5+120*tan(17.0*deg))*mm,-50.0*mm,-50.0*mm,
			-40.5*mm,-40.5*mm,40.5*mm,40.5*mm,50.0*mm,
			50.0*mm,(40.5+120*tan(17.0*deg))*mm};
		G4Polycone* LWedgePolycone = new G4Polycone("LWedgePcon",-20.0*deg,40.0*deg,
			kNPlane_LWedge,zPlane_LWedge,rInner_LWedge,rOuter_LWedge);
		G4Tubs* LWedgeTubs = new G4Tubs("LWedgeTubs",0.0*mm,38.2*mm,1000/2.0*mm,
			startphi,deltaphi);//The orininal radius of this hole is 31.5 mm, change it to 38.2 mm
		G4Tubs* LWedgeTubs2 = new G4Tubs("LWedgeTubs2",0.0*mm,130.0*mm,40.0/2.0*mm,startphi,deltaphi);
		G4RotationMatrix* pRotWedgeTubs = new G4RotationMatrix();
		pRotWedgeTubs->rotateY(90.*deg);
		G4SubtractionSolid* LWedgeS1 = new G4SubtractionSolid("LWedgeS1",
			LWedgePolycone,LWedgeTubs,pRotWedgeTubs,G4ThreeVector());
		G4SubtractionSolid* LWedgeSolid = new G4SubtractionSolid("LWedgeSolid",LWedgeS1,LWedgeTubs2);
		G4LogicalVolume* LWedgeLogical = new G4LogicalVolume(LWedgeSolid,
			mMaterialManager->aluminum,"LWedgeLogical",0,0,0);
		LWedgeLogical->SetVisAttributes(SilverVisAtt);
		G4RotationMatrix* pRotWedge1 = new G4RotationMatrix();
		pRotWedge1->rotateY(90.*deg);
		pRotWedge1->rotateX(90.*deg);
		new G4PVPlacement(pRotWedge1,G4ThreeVector(),
			LWedgeLogical,"LargeWedgePhys",innerOfChamberLogical,true,0,0);
		G4RotationMatrix* pRotWedge2 = new G4RotationMatrix();
		pRotWedge2->rotateY(-90.*deg);
		pRotWedge2->rotateX(90.*deg);
		new G4PVPlacement(pRotWedge2,G4ThreeVector(),
			LWedgeLogical,"LargeWedgePhys",innerOfChamberLogical,true,1,0);

		//Small Wedge
		const int kNPlane_SWedge=10;
		double rInner_SWedge[] = {250.0*mm,130.0*mm,112.0*mm,112.0*mm,
			130.0*mm,130.0*mm,112.0*mm,112.0*mm,130.0*mm,250.0*mm};
		double rOuter_SWedge[] = {250.0*mm,250.0*mm,250.0*mm,250.0*mm,
			250.0*mm,250.0*mm,250.0*mm,250.0*mm,250.0*mm,250.0*mm};
		double zPlane_SWedge[] = {-(40.5+120*tan(17.0*deg))*mm,-40.5*mm,-40.5*mm,-20.0*mm,
			-20.0*mm,20.0*mm,20.0*mm,40.5*mm,40.5*mm,(40.5+120*tan(17.0*deg))*mm};
		G4Polycone* SWedgeSolid = new G4Polycone("SWedgeSolid",-7.5/2.0*deg,7.5*deg,
			kNPlane_SWedge,zPlane_SWedge,rInner_SWedge,rOuter_SWedge);
		G4LogicalVolume* SWedgeLogical= new G4LogicalVolume(SWedgeSolid,
			mMaterialManager->stainlesssteel,"SWedgeLogical",0,0,0);
		SWedgeLogical->SetVisAttributes(GrayVisAtt);
		G4RotationMatrix* pRotWedge3 = new G4RotationMatrix();
		pRotWedge3->rotateX(90.*deg);
		pRotWedge3->rotateZ(30.*deg);
		new G4PVPlacement(pRotWedge3,G4ThreeVector(),
			SWedgeLogical,"SmallWedgePhys",innerOfChamberLogical,true,0,0);
		G4RotationMatrix* pRotWedge4 = new G4RotationMatrix();
		pRotWedge4->rotateX(90.*deg);
		pRotWedge4->rotateZ(150.*deg);
		new G4PVPlacement(pRotWedge4,G4ThreeVector(),
			SWedgeLogical,"SmallWedgePhys",innerOfChamberLogical,true,1,0);
		G4RotationMatrix* pRotWedge5 = new G4RotationMatrix();
		pRotWedge5->rotateX(90.*deg);
		pRotWedge5->rotateZ(-30.*deg);
		new G4PVPlacement(pRotWedge5,G4ThreeVector(),
			SWedgeLogical,"SmallWedgePhys",innerOfChamberLogical,true,2,0);
		G4RotationMatrix* pRotWedge6 = new G4RotationMatrix();
		pRotWedge6->rotateX(90.*deg);
		pRotWedge6->rotateZ(-150.*deg);
		new G4PVPlacement(pRotWedge6,G4ThreeVector(),
			SWedgeLogical,"SmallWedgePhys",innerOfChamberLogical,true,3,0);

	}

	/////////////////////////
	// LHe Shielding (4k shielding)
	/////////////////////////
	//LHe Shield cylider with r=38.1mm and L=700mm, thickness=0.0015'=0.0381mm
	//double mShieldLHeRin=38.1*mm,mShieldLHeRout=38.1381*mm, mShieldLHeL=700.0*mm;
	startphi=0.0*deg;deltaphi=360.0*deg;
	G4VSolid* shieldLHeSolid = new G4Tubs("shieldLHeTubs",mShieldLHeRin,mShieldLHeRout,
		mShieldLHeL/2.0,startphi,deltaphi);
	G4LogicalVolume* shieldLHeLogical = new G4LogicalVolume(shieldLHeSolid,
		mMaterialManager->aluminum,"shieldLHeLogical",0,0,0);
	shieldLHeLogical->SetVisAttributes(WhiteVisAtt);  //white

	new G4PVPlacement(0,G4ThreeVector(),shieldLHeLogical,"shieldLHePhys",innerOfChamberLogical,0,0);

	return scatChamberContainerPhys;
}


/////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* G2PDetectorConstruction::ConstructG2PTarget(G4LogicalVolume* motherLogical)
{
	G4String SDname;

	G4VSensitiveDetector* targetSD=new HRSDCSD(SDname="targetMaterial");
	G4VSensitiveDetector* targetNoseSD=new HRSDCSD(SDname="targetNose");
	G4VSensitiveDetector* targetEndCapSD=new HRSDCSD(SDname="targetEndCap");
	G4VSensitiveDetector* targetWallSD=new HRSDCSD(SDname="targetWall");

	G4SDManager* SDman = G4SDManager::GetSDMpointer();

	const double inch=2.54*cm;
	double startphi,deltaphi;

	G4RotationMatrix* pRotX90deg = new G4RotationMatrix();
	pRotX90deg->rotateX(90.*deg);
	G4RotationMatrix* pRotX270deg = new G4RotationMatrix();
	pRotX270deg->rotateX(270.*deg);

	//////////////////////////////////
	//target container
	//////////////////////////////////

	//The target containner should stay in the center of Scatter chamber, otherwise
	//there will be some overlapping
	//Note that the target offsets are in the hall coordinate system, I have to shift it 
	//to the target containner system(scat chamber system) 
	//after pRotX90deg, (x,y,z)_ScatChamber==>(x,-z,y)_Hall
	double pTargetOriginX=mTargetXOffset-mScatChamberXOffset;
	double pTargetOriginY=-(mTargetZOffset-mScatChamberZOffset);
	double pTargetOriginZ=mTargetYOffset-mScatChamberYOffset;

	startphi=0.0*deg;deltaphi=360.0*deg;
	G4VSolid* targetContainerSolid = new G4Tubs("targetContainerTubs",0,mShieldLHeRin,
		mShieldLHeL/2.0,startphi,deltaphi);

	G4LogicalVolume* targetContainerLogical = new G4LogicalVolume(targetContainerSolid,
		mMaterialManager->vacuum,"targetContainerLogical",0,0,0);
	targetContainerLogical->SetVisAttributes(HallVisAtt);  //white invisible

	G4VPhysicalVolume* targetContainerPhys=new G4PVPlacement(pRotX90deg,
		G4ThreeVector(mScatChamberXOffset,mScatChamberYOffset,mScatChamberZOffset),
		targetContainerLogical,"targetContainerPhys",motherLogical,0,0);


	//////////////////////////////////
	//tail|nose shiled
	//////////////////////////////////
	//By Jixie: told by Josh that the thickness is 4 mil
	//Toby get the number from Josh: Inner diameter = 1.654", Outter diameter=1.662"

	double pShieldNoseRin=1.654/2*inch,pShieldNoseRout=pShieldNoseRin+0.004*inch; 
	double pShieldNoseL=mShieldLHeL-10.0*mm;
	startphi=0.0*deg;deltaphi=360.0*deg;
	G4VSolid* shieldNoseSolid = new G4Tubs("shieldNoseTubs",pShieldNoseRin,pShieldNoseRout,
		pShieldNoseL/2.0,startphi,deltaphi);

	G4LogicalVolume* shieldNoseLogical = new G4LogicalVolume(shieldNoseSolid,
		mMaterialManager->aluminum,"shieldNoseLogical",0,0,0);
	shieldNoseLogical->SetVisAttributes(WhiteVisAtt);  //white

	//since I place the target inside the nose, I will put the center of the nose at the
	//target position, therefore only the nose should use the target offset
	new G4PVPlacement(0,G4ThreeVector(pTargetOriginX,pTargetOriginY,pTargetOriginZ),
		shieldNoseLogical,"shieldNosePhys",targetContainerLogical,0,0);


	//////////////////////////////////
	//inside the tail|nose 
	//////////////////////////////////

	startphi=0.0*deg;deltaphi=360.0*deg;
	G4VSolid* innerOfNoseSolid = new G4Tubs("innerOfNoseTubs",0,pShieldNoseRin,
		pShieldNoseL/2.0,startphi,deltaphi);
	G4Material* theCoolant=mMaterialManager->liquidHe;
	//for optics, we do not put liquid helium inside the tail|nose as high as normal runs
	//I will use heliumGas4k instead
	if (mSetupG2PTarget>=10) theCoolant=mMaterialManager->heliumGas4k;

	G4LogicalVolume* innerOfNoseLogical = new G4LogicalVolume(innerOfNoseSolid,
		theCoolant,"innerOfNoseLogical",0,0,0);
	innerOfNoseLogical->SetVisAttributes(LightPurpleVisAtt);  //light purple
	
	SDman->AddNewDetector(targetNoseSD);
	innerOfNoseLogical->SetSensitiveDetector(targetNoseSD);

	//By Jixie: Add this step limit can help to calculate the integrated BdL
	double pTargetStepLimit=10;
	gConfig->GetArgument("TargetStepLimit",pTargetStepLimit);
	pTargetStepLimit*=mm;
	G4UserLimits* uTNStepLimits = new G4UserLimits(pTargetStepLimit);
	innerOfNoseLogical->SetUserLimits(uTNStepLimits);

	//since I place the target inside the nose, I will put the center of the nose at the
	//target position, therefore only the nose should use the target offset
	new G4PVPlacement(0,G4ThreeVector(pTargetOriginX,pTargetOriginY,pTargetOriginZ),
		innerOfNoseLogical,"innerOfNosePhys",targetContainerLogical,0,0);

	//////////////////////////////////
	//target insert, target wall and target cap
	//////////////////////////////////

	//#SetupG2PTarget is used to determine whether to setup target wall and also determine the
	//#the radius of the target cell and the target position w.r.t the target nose. 
	//#It can be 0, 1, 2, 3, 6, 10, 20, 30, 60 or 99
	//#0 means do not set up the target, 
	//#1 means it is the 1st target cell(R_out=1.142"/2 and R_in=TargetR=1.072"/2, centered at zero), 
	//#  for NH3He target. Can have target wall. The 4th and 5th target cell are identical to the 1st cell. 
	//#2 is the middle size hole(D=0.397", centered at -10.95mm) for C12. No target wall 
	//#3 is the smallest hole(D=0.312", centered at -10.95mm) for CH2. No target wall 
	//#6 is the 6th target cell(R_out=1.142"/2 and R_in=TargetR=0.960"/2, with the front face aligned to
	//#  the upstream end, z_face=-1.113"/2=-14.135mm) for C12. With target wall, No end caps. 
	//#10, 20, 30, and 60 are the the same as 1,2,3 and 6,respectively, but NO LHe inside the target nose.
	//#99 means use a pvc pipe and a 10-mil C12 target foil located at z=-896.0 mm, which is the situation 
	//#   of commissioning. In this setup, no target wall, no target cell, TargetType, TargetR, TargetL above 
	//#   will be ignored 

	//The target wall and end caps need to be rotated in its mother volumn 
	//in target nose coordinate system, 
	//x'=x_lab, y'=-z_lab, z'=y_lab

	//The 2nd and 3rd hole are smaller, their radii will be set later according to mSetupG2PTarget
	//R_2nd=0.397"/2   R_3rd=0.312"/2
	//Based on the drawing, the hole diameter of the insert is 1.098" 
	//But the outer diameter of the target cell is 1.142"
	//In order to simplify the geometry and avoid overlappings, I use 1.142" here
	double pTargetInsertHoleR=1.142*inch/2; 
	if(mSetupG2PTarget==2 || mSetupG2PTarget==20)  pTargetInsertHoleR=0.397*inch; 
	else if(mSetupG2PTarget==3 || mSetupG2PTarget==30) pTargetInsertHoleR=0.312*inch; 

	//////////////////////////////////////////////////////////////////////////
	//Note that if mSetupG2PTarget==2|3 then no need to put target wall and caps
	
	//#Setup the target wall, will setup target wall only if SetupG2PTarget==1,6 or 10,60
	int pSetupTargetWall=1; 
	if(mSetupG2PTarget!=1 && mSetupG2PTarget!=10 && mSetupG2PTarget!=6 && mSetupG2PTarget!=60) 
	{
		pSetupTargetWall=0;
	}
	int pSetupEndCap=1;
	if(mSetupG2PTarget!=1 && mSetupG2PTarget!=10) pSetupEndCap=0;

	
	//#The target wall is PCTFE, 1st, 3rd and 5th cells are 35 mil or 0.899 mm thick
	//#But 6th cell (for C12) is 90 mil
	//TargetWallThick=0.889;  #2.286 for cell 6 or 60 for C12
	double pTargetWallL=1.113*inch;
	double pTargetWallThick=0.035*inch;
	if(mSetupG2PTarget==6 || mSetupG2PTarget==60) pTargetWallThick=0.090*inch;

	double pTargetR=pTargetInsertHoleR;	
	if(pSetupTargetWall)  pTargetR=pTargetInsertHoleR-pTargetWallThick;
	//////////////////////////////////////////////////////////////////////////


	//////////////////////////////////
	//the insert, just build one block, not 6 blocks
	/////////////////////////////////
	G4VSolid* targetInsertHoleSolid = new G4Tubs("targetInsertHoleTubs",
		0,pTargetInsertHoleR, 0.13*inch/2,0.,360.*deg);

	double pTargetInsertX_Lab=1.156*inch; //in lab system
	double pTargetInsertY_Lab=1.230*inch; //in lab system
	double pTargetInsertZ_Lab=0.125*inch; //in lab system 
	G4VSolid* targetInsertWholeSolid = new G4Box("targetInsertWholeBox",
		pTargetInsertX_Lab/2,pTargetInsertY_Lab/2,pTargetInsertZ_Lab/2);

	//now dig the hole
	G4SubtractionSolid*	targetInsertSolid=new G4SubtractionSolid("targetInsertSolid",
		targetInsertWholeSolid,targetInsertHoleSolid);

	G4LogicalVolume* targetInsertLogical = new G4LogicalVolume(targetInsertSolid,
		mMaterialManager->aluminum, "targetInsertLogical",0,0,0);
	targetInsertLogical->SetVisAttributes(GrayVisAtt);

	//Insert position w.r.t nose, should be positive before 270 deg rotation
	//based on the drawing, it is 1.113/2-0.188+0.125/2 = 0.4310"=10.95mm
	double pTargetInsertInNoseCS=0.4310*inch;  //in nose system, 
											   
	new G4PVPlacement(pRotX270deg,
		G4ThreeVector(0.0,pTargetInsertInNoseCS,0.0),
		targetInsertLogical,"targetInsertPhys",innerOfNoseLogical,0,0);


	//////////////////////////////////////////////////////////////////////////

	if(pSetupTargetWall)
	{
		////////////////////////////////////////
		//target cell holder
		////////////////////////////////////////
		double pTargetCellHolderZ_Lab=0.10*inch; //in lab system 
		G4VSolid* targetCellHolderWholeSolid = new G4Box("targetCellHolderWholeBox",
			pTargetInsertX_Lab/2,pTargetInsertY_Lab/2,pTargetCellHolderZ_Lab/2);

		G4SubtractionSolid*	targetCellHolderSolid=new G4SubtractionSolid("targetCellHolderSolid",
			targetCellHolderWholeSolid,targetInsertHoleSolid);

		G4LogicalVolume* targetCellHolderLogical = new G4LogicalVolume(targetCellHolderSolid,
			mMaterialManager->PCTFE,"targetCellHolderLogical",0,0,0);
		targetCellHolderLogical->SetVisAttributes(LightYellowVisAtt);  

		double pTargetCellHolderInNoseCS=pTargetInsertInNoseCS-
			pTargetInsertZ_Lab/2-pTargetCellHolderZ_Lab/2;
		new G4PVPlacement(pRotX270deg,
			G4ThreeVector(0.0,pTargetCellHolderInNoseCS,0.0),
			targetCellHolderLogical,"targetCellHolderPhys",innerOfNoseLogical,0,0);

		////////////////////////////////////////
		//target wall
		////////////////////////////////////////
		G4VSolid* targetWallSolid = new G4Tubs("targetWallTubs",pTargetR,
			pTargetInsertHoleR,pTargetWallL/2.0,0.,360.*deg);

		G4LogicalVolume* targetWallLogical = new G4LogicalVolume(targetWallSolid,
			mMaterialManager->PCTFE,"targetWallLogical",0,0,0);
		targetWallLogical->SetVisAttributes(LightYellowVisAtt);
	
		SDman->AddNewDetector(targetWallSD);
		targetWallLogical->SetSensitiveDetector(targetWallSD);


		new G4PVPlacement(pRotX270deg,G4ThreeVector(),
			targetWallLogical,"targetWallPhys",innerOfNoseLogical,0,0);
	}

	////////////////////////////////////////
	//target entrance cap and exit cap
	////////////////////////////////////////
	//The end cap is 0.7 mil aluminum foil, told by Josh 
	double pTargetCapThick=0.007*inch;
	if(pSetupEndCap) 
	{
		//double mTargetCapThick=0.0007*inch;
		double pTargetCapZ=pTargetWallL/2.0+pTargetCapThick/2.0;
		G4VSolid* targetCapSolid = new G4Tubs("targetCapTubs",0,pTargetR,
			pTargetCapThick/2.0,0.,360.*deg);

		G4LogicalVolume* targetCapLogical = new G4LogicalVolume(targetCapSolid,
			mMaterialManager->aluminum,"targetCapLogical",0,0,0);
		targetCapLogical->SetVisAttributes(GrayVisAtt);  //Grey invisable

		SDman->AddNewDetector(targetEndCapSD);
		targetCapLogical->SetSensitiveDetector(targetEndCapSD);

		//Entrance cap, 
		new G4PVPlacement(pRotX270deg,
			G4ThreeVector(0.0,pTargetCapZ,0.0),
			targetCapLogical,"targetEntranceCapPhys",innerOfNoseLogical,true,0);
		//Exit cap
		new G4PVPlacement(pRotX270deg,
			G4ThreeVector(0.0,-pTargetCapZ,0.0),
			targetCapLogical,"targetExitCapPhys",innerOfNoseLogical,true,1);
	}
	

	//////////////////////////////////
	//the target,
	//////////////////////////////////

	//major target is 55% solid NH3 mixed with LHe cylinder with r=13.6144mm and z=28.27mm 
	//double pTargetR=13.6144*mm, mTargetL=28.27*mm;
	G4VSolid* targetSolid = new G4Tubs("targetTubs",0.,pTargetR,mTargetL/2.0,0.,360.*deg);

	//for optics, we do not put liquid helium inside the tail|nose as high as normal run
	//I will use vacuum instead 
	G4Material* theTargetMaterial=G2P_NH3He;
	if (mTargetType==0) theTargetMaterial=mMaterialManager->vacuum;
	else if (mTargetType==1) theTargetMaterial=this->G2P_NH3He;
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
	else if (mTargetType==18) theTargetMaterial=mMaterialManager->NH3He;

	G4LogicalVolume* targetLogical = new G4LogicalVolume(targetSolid,theTargetMaterial,
		"targetLogical",0,0,0);
	targetLogical->SetVisAttributes(VioletVisAtt); 
	targetLogical->SetUserLimits(uTNStepLimits);

	//By Jixie: to study the radiator, I make the target sensitive
	SDman->AddNewDetector(targetSD);
	targetLogical->SetSensitiveDetector(targetSD);

	//since I place the target inside the nose, I will put the center of the nose at the
	//target position, therefore only the nose should use the target offset

	//this number is in target nose coordinate system, 
	//x'=x, y'=-z, z'=y
	double pTargetInNoseCS=0;
	if(mSetupG2PTarget==2 || mSetupG2PTarget==20 || mSetupG2PTarget==3 || mSetupG2PTarget==30) 
	{
		pTargetInNoseCS=pTargetInsertInNoseCS; //0.431 *inch;
	}
	else if(mSetupG2PTarget==6 || mSetupG2PTarget==60) 
	{
		pTargetInNoseCS=pTargetWallL/2-mTargetL/2; 
	}
	new G4PVPlacement(pRotX270deg,G4ThreeVector(0,pTargetInNoseCS,0),
		targetLogical,"targetPhys",innerOfNoseLogical,0,0);

	/////////////////////////////////////////////////////////////////////

	return targetContainerPhys;
}


//set SetupG2PTarget to 99 to invoke this PVC target
//G2p commissioning target, a PVC pipe with a 10 mil C12 foil
//all parameters are hard coded, the config file can not change any of them
G4VPhysicalVolume* G2PDetectorConstruction::ConstructG2PPVCTarget(G4LogicalVolume* motherLogical)
{
	const double inch=2.54*cm;
	double startphi,deltaphi;

	//////////////////////////////////
	//target container
	//////////////////////////////////

	double pTargetContainerRout=40.0*mm, pTargetContainerL=1201.0*mm;
	startphi=0.0*deg;deltaphi=360.0*deg;
	G4VSolid* targetContainerSolid = new G4Tubs("pTargetContainerTubs",0,
		pTargetContainerRout,pTargetContainerL/2.0,startphi,deltaphi);

	G4LogicalVolume* targetContainerLogical = new G4LogicalVolume(targetContainerSolid
		,mMaterialManager->heliumGas,"targetContainerLogical",0,0,0);
	targetContainerLogical->SetVisAttributes(HallVisAtt);  

	G4VPhysicalVolume* targetContainerPhys=new G4PVPlacement(0,
		G4ThreeVector(0,0,-896.0*mm),
		targetContainerLogical,"targetContainerPhys",motherLogical,0,0);


	//////////////////////////////////
	//PVC pipe
	//////////////////////////////////
	//PVC pipe, 2.377" outer diameter and 2.027" inner diameter

	double PVCPipeRout=2.377*inch/2, PVCPipeRin=2.027*inch/2, PVCPipeL=1200.0*mm;
	startphi=0.0*deg;deltaphi=360.0*deg;
	G4VSolid* PVCPipeSolid = new G4Tubs("PVCPipeTubs",PVCPipeRin,PVCPipeRout,
		PVCPipeL/2.0,startphi,deltaphi);

	G4LogicalVolume* PVCPipeLogical = new G4LogicalVolume(PVCPipeSolid,mMaterialManager->epoxy,
		"PVCPipeLogical",0,0,0);
	PVCPipeLogical->SetVisAttributes(LightPurpleVisAtt);  

	new G4PVPlacement(0,G4ThreeVector(0,0,0),
		PVCPipeLogical,"PVCPipePhys",targetContainerLogical,0,0);


	//////////////////////////////////
	//PVC pipe holder, downstream one only
	//////////////////////////////////
	//It is a ring union with 2 ears, 
	//The ear is 0.722" in x (from inner of the ring), 2 X 0.174" in y, 0.581" in z
	//can not figure out the position of it, place it 50 cm downstream

	double PVCHolderRin=PVCPipeRout+0.1*mm, PVCHolderRout=PVCHolderRin+0.174*inch, PVCHolderL=0.581*inch;
	startphi=0.0*deg;deltaphi=360.0*deg;
	G4VSolid* PVCHolderRingSolid = new G4Tubs("PVCHolderRingTubs",PVCHolderRin,PVCHolderRout,
		PVCHolderL/2.0,startphi,deltaphi);
	double PVCHolderEarX=0.722*inch;
	G4VSolid* PVCHolderEarSolid = new G4Box("PVCHolderEarBox",PVCHolderEarX/2,0.174*inch,PVCHolderL/2);
	double PVCEarPos_X=PVCHolderRin+PVCHolderEarX/2+0.3*mm;
	G4UnionSolid* PVCHolderRingUEarLSolid = new G4UnionSolid("PVCHolderRingUEarLSolid",
		PVCHolderRingSolid,PVCHolderEarSolid,0,G4ThreeVector(PVCEarPos_X,0,0));
	G4UnionSolid* PVCHolderRingUEarLRSolid = new G4UnionSolid("PVCHolderRingUEarLRSolid",
		PVCHolderRingUEarLSolid,PVCHolderEarSolid,0,G4ThreeVector(-PVCEarPos_X,0,0));

	G4LogicalVolume* PVCHolderLogical = new G4LogicalVolume(PVCHolderRingUEarLRSolid,
		mMaterialManager->stainlesssteel,"PVCHolderLogical",0,0,0);
	PVCHolderLogical->SetVisAttributes(SteelVisAtt);  

	new G4PVPlacement(0,G4ThreeVector(0,0,50*cm),
		PVCHolderLogical,"PVCHolderPhys",targetContainerLogical,0,0);


	//////////////////////////////////
	//the target,
	//////////////////////////////////

	//commissioning target is 10 mil C12 target 
	G4VSolid* targetSolid = new G4Box("targetFoil",1.0*cm,1.0*cm,0.254*mm/2);

	G4LogicalVolume* targetLogical = new G4LogicalVolume(targetSolid,mMaterialManager->carbon,
		"targetLogical",0,0,0);
	targetLogical->SetVisAttributes(VioletVisAtt); 

	//since I place the target inside the nose, I will put the center of the nose at the
	//target position, therefore only the nose should use the target offset
	new G4PVPlacement(0,G4ThreeVector(0,0,0),
		targetLogical,"targetPhys",targetContainerLogical,0,0);


	/////////////////////////////////////////////////////////////////////

	return targetContainerPhys;
}

/////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* G2PDetectorConstruction::ConstructG2PSeptumNSieve(G4LogicalVolume* motherLogical)
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
	//double mPivot2LSieveFace=799.6*cm;
	//double mPivot2RSieveFace=799.6*cm;
	//Below are the results from the G2p sieve slit, collimator, and dump survey carried out on Feb 23rd, 2012. 
	//	Values are relative to the ideal target center with a +X to the beam left, a +Y up, and a +Z downstream. 
	//	Values are in millimeters and degrees. All coordinates are to the upstream face of the components. 
	//	The sieve coordinates are to the center hole.
	//	LOCATION		Z		X		Y		YAW		PITCH		ROLL
	//	Dump			639.21	-0.98	4.90	0.103	0.255		0.313
	//	L.Collimator	639.11	66.33	1.50	0.103	0.255		0.073
	//	R.Collimator	639.35	-66.66	-0.04	0.103	0.255		0.180
	//	L.Sieve CL		795.55	80.35	0.31	5.682	0.141		0.257
	//	R.Sieve CL		795.41	-80.39	0.21	-5.647	-0.183		0.305


	double pSieveSlitX=2.205*inch; //33.13*mm
	double pSieveSlitY=5.134*inch; //130.40*mm
	double pSieveSlitZ=0.2*inch;

	double pSieveSlitHoleR=0.6985*mm;           //radius of small hole 0.055/2 inch
	double pSieveSlitLargeHoleR=1.3462*mm;      //radius of large hole 0.106/2 inch
	double pSieveSlitHoldL=pSieveSlitZ+0.1*mm;  //need to make it longer to avoid round off in the subtraction

	//the whole position relative to the slit center 
	double pSieveSlitHolePosH[7], pSieveSlitHolePosV[7];

	//the big center hole horizontal and vertical offset, From Alan
	double pSieveSlitLargeHoleH=2.967*mm;    //positive means shift away from the beam line
	double pSieveSlitLargeHoleV=2.606*mm;

	//Updated By Jixie@20121129, change these number with survey results
	//left and right sieve are almost equal, so I just use the left sieve survey results
	//These numbers are from the survey,  
	pSieveSlitLargeHoleH=80.35*mm-799.6*mm*sin(mLSeptumAngle);   //1.63 mm for 5.65 degree, or 1.08 mm for
	pSieveSlitLargeHoleV=0.31*mm;

	//please note that the following constants are for the sieve in the right arm only
	//need to mirror(flip) it in order to put in the left arm  
	double pSieveSlitDeltaH[8]={0.537*inch, 0.188*inch, 0.188*inch, 0.188*inch,
		0.241*inch, 0.241*inch, 0.241*inch, 0.381*inch}; //in inch, from left to right
	double pSieveSlitDeltaV[8]={0.496*inch, 0.524*inch, 0.524*inch, 0.524*inch,
		0.524*inch, 0.524*inch, 0.524*inch, 1.494*inch}; //in inch, from top to buttom

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
		//located at z=686.46 mm, no rotation
		double pSeptumX=140.0*cm;
		double pSeptumY=84.4*cm;
		double pSeptumZ=74.0*cm;
		double pSeptumTunnelX=30.4*cm;
		double pSeptumTunnelY=24.4*cm-2.0*inch;  //By Jixie @20120205: Add 2 inches of iron
		double pSeptumBeamHoleX=7.8*cm;
		double pSeptumBeamHoleY=8.0*cm;

		double pSeptumTunnelPos_X=8.4*cm+pSeptumTunnelX/2.0;
		double pSeptumPos_Z=68.646*cm;		
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
		G4RotationMatrix* pSeptumCoilRotBackDwon = new G4RotationMatrix();
		pSeptumCoilRotBackDwon->rotateY(270*deg);
		G4RotationMatrix* pSeptumCoilRotBackUp = new G4RotationMatrix();
		pSeptumCoilRotBackUp->rotateY(270*deg);
		pSeptumCoilRotBackUp->rotateX(180*deg);

		new G4PVPlacement(pSeptumCoilRotBackDwon,
			G4ThreeVector(pSeptumCoilPos_X_in,pSeptumCoilPos_Y,pSeptumCoilPos_Z_down),
			septumCoilLogical,"septumCoilPhys",motherLogical,true,8,0);
		new G4PVPlacement(pSeptumCoilRotBackUp,
			G4ThreeVector(pSeptumCoilPos_X_in,-pSeptumCoilPos_Y,pSeptumCoilPos_Z_down),
			septumCoilLogical,"septumCoilPhys",motherLogical,true,9,0);
		new G4PVPlacement(pSeptumCoilRotBackUp,
			G4ThreeVector(pSeptumCoilPos_X_out,-pSeptumCoilPos_Y,pSeptumCoilPos_Z_down),
			septumCoilLogical,"septumCoilPhys",motherLogical,true,10,0);
		new G4PVPlacement(pSeptumCoilRotBackDwon,
			G4ThreeVector(pSeptumCoilPos_X_out,pSeptumCoilPos_Y,pSeptumCoilPos_Z_down),
			septumCoilLogical,"septumCoilPhys",motherLogical,true,11,0);

		new G4PVPlacement(pSeptumCoilRotBackDwon,
			G4ThreeVector(-pSeptumCoilPos_X_out,pSeptumCoilPos_Y,pSeptumCoilPos_Z_down),
			septumCoilLogical,"septumCoilPhys",motherLogical,true,12,0);
		new G4PVPlacement(pSeptumCoilRotBackUp,
			G4ThreeVector(-pSeptumCoilPos_X_out,-pSeptumCoilPos_Y,pSeptumCoilPos_Z_down),
			septumCoilLogical,"septumCoilPhys",motherLogical,true,13,0);
		new G4PVPlacement(pSeptumCoilRotBackUp,
			G4ThreeVector(-pSeptumCoilPos_X_in,-pSeptumCoilPos_Y,pSeptumCoilPos_Z_down),
			septumCoilLogical,"septumCoilPhys",motherLogical,true,14,0);
		new G4PVPlacement(pSeptumCoilRotBackDwon,
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
			double pLHRSVBPos_X=(pHallCenter2VB-mHRSVBThick/2)*
				sin(mLHRSAngle)+mPivotXOffset;
			double pLHRSVBPos_Y=mPivotYOffset;
			//no need to correct for pivot since the distance is from the hall center
			double pLHRSVBPos_Z=(pHallCenter2VB-mHRSVBThick/2.0)*cos(mLHRSAngle);
			new G4PVPlacement(pRotLHRS,
				G4ThreeVector(pLHRSVBPos_X,pLHRSVBPos_Y,pLHRSVBPos_Z),
				HRSVBLogical,"virtualBoundaryPhys_LHRS",motherLogical,0,0,0);
		}
		if(mSetupRHRS)
		{
			double pRHRSVBPos_X=(pHallCenter2VB-mHRSVBThick/2)*
				sin(mRHRSAngle)+mPivotXOffset;
			double pRHRSVBPos_Y=mPivotYOffset;
			//no need to correct for pivot since the distance is from the hall center
			double pRHRSVBPos_Z=(pHallCenter2VB-mHRSVBThick/2)*cos(mRHRSAngle); 
			new G4PVPlacement(pRotRHRS,
				G4ThreeVector(pRHRSVBPos_X,pRHRSVBPos_Y,pRHRSVBPos_Z),
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


/////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* G2PDetectorConstruction::ConstructHRS(G4LogicalVolume* motherLogical)
{
	//const double inch=2.54*cm;
	G4VPhysicalVolume* theG2PHRSPhys=0;

	G4RotationMatrix *pRotX90deg=new G4RotationMatrix();
	pRotX90deg->rotateX(-90*deg); 

	//Build the HRS, positions were taken fron the NIM paper
	//The Apertures are defined as:
	//Q1 exit is a circle of radius 0.15m
	//Q3 entrance and exit are circles of radius 0.30m
	//The dipole entrance and exit are trapezoids:
	//-0.40 m<x<0.40 m (+x is down) 
	//y=+-(0.125*(1-(1.25*x/8.40)) (smallest gap at positive x, x in m)

	/////////////////////////
	// HRS QQDQ Containner
	/////////////////////////
	//Build a container using polycone, covering 270+/-12 degrees, 1.46m to 20.76m
	//3.05m below the beam line is the ground, need to subtract everything beneath that
	//looks like this:
/*
    HRS container:covering 24 degrees in X-Z plane, 10.05 m height.
    Stuff inside this containner will also use the pivot as the origin, but do not 
    need to worry about the rotation of the HRS.
                        7.0m above beam line 
                       -----------------------------| 20.76m from pivot,  
                      /                             |
                     /                        Q3    |
                    /                               |
                   /                       E        |
                  /                     L           |
          --------                   O              |
         /                        P                 |
    -----                     I                     |
----|----- Q1 --- Q2 ---- D     --------------------|------ beam line -----
    |                                               |   
    ------------------------------------------------|
    1.46m from povot, 3.05 m below beam line 
*/

	double pHRSContainerRin=1.46*m,pHRSContainerRout=20.76*m;
	double pBeamLine2Ground;
	pBeamLine2Ground=-3.05*m;
	//build the container with polycone

	const int kNPlane_HRSContainer=7;
	double rInner_HRSContainer[] = {pHRSContainerRin,pHRSContainerRin,2.5*m,
		3.7*m,9.0*m,pHRSContainerRout-3.0*m,pHRSContainerRout};
	double rOuter_HRSContainer[] = {pHRSContainerRout,pHRSContainerRout,pHRSContainerRout,
		pHRSContainerRout,pHRSContainerRout,pHRSContainerRout,pHRSContainerRout};
	double zPlane_HRSContainer[] = {-2.0*m,1.0*m,1.0*m,
		2.0*m,2.0*m,7.0*m,7.0*m};
	G4Polycone* HRSContainerSolid = new G4Polycone("HRSContainer",258.0*deg,24.0*deg,
		kNPlane_HRSContainer,zPlane_HRSContainer,rInner_HRSContainer,rOuter_HRSContainer);

	////build the container using tube
	//double pHRSContainerHeight=14.0*m;
	//G4VSolid* HRSContainerTub = new G4Tubs("HRSContainerTub",pHRSContainerRin,
	//	pHRSContainerRout,pHRSContainerHeight/2.0,258.0*deg,24.0*deg);
	////ground
	//double pGoundHeight=10*m;
	//G4VSolid* GroundTub = new G4Tubs("GroundTub",0,30*m,pGoundHeight/2,0*deg,360.0*deg);
	////tube subtract the ground
	//G4SubtractionSolid* HRSContainerSolid = new G4SubtractionSolid("HRSContainer",
	//		HRSContainerTub,GroundTub,
	//		0,G4ThreeVector(0,0,pBeamLine2Ground-pGoundHeight/2));

	G4LogicalVolume* LHRSContainerLogical = new G4LogicalVolume(HRSContainerSolid,
		mMaterialManager->vacuum,"LHRSContainerLogical",0,0,0);
	G4LogicalVolume* RHRSContainerLogical = new G4LogicalVolume(HRSContainerSolid,
		mMaterialManager->vacuum,"RHRSContainerLogical",0,0,0);

	LHRSContainerLogical->SetVisAttributes(HallVisAtt); 
	RHRSContainerLogical->SetVisAttributes(HallVisAtt); 

	G4RotationMatrix *pRotLHRSContainer=new G4RotationMatrix();
	pRotLHRSContainer->rotateX(-270*deg);
	pRotLHRSContainer->rotateZ(-mLHRSAngle);  

	G4RotationMatrix *pRotRHRSContainer=new G4RotationMatrix();
	pRotRHRSContainer->rotateX(-270*deg);
	pRotRHRSContainer->rotateZ(-mRHRSAngle);  

	if(mSetupLHRS>=2)
	{
		new G4PVPlacement(pRotLHRSContainer,G4ThreeVector(0,0,0),
			LHRSContainerLogical,"LHRSContainerPhys",motherLogical,0,0,0);
	}
	if(mSetupRHRS>=2)
	{
		new G4PVPlacement(pRotRHRSContainer,G4ThreeVector(0,0,0),
			RHRSContainerLogical,"RHRSContainerPhys",motherLogical,0,0,0);
	}


	/////////////////////////
	// HRS Q1 
	/////////////////////////
	double pHallCenter2LQ1Face=1.69*m;
	double pHallCenter2RQ1Face=1.69*m;
	double pQ1Rin=15.0*cm;
	double pQ1Rout=35.0*cm;
	double pQ1Length=80*cm;
	//double pQ1Length=(1.698M-1.36*m)*2+0.8*m;

	//By Jixie: HepRApp.jar has a bug to view this tubs. I am sure this tube is oriented correctly
	G4VSolid* Q1Solid = new G4Tubs("Q1Tub",pQ1Rin,pQ1Rout,pQ1Length/2.0,0.0,360.0*deg);

	//build 2 copy since there are differennt fields involved in it
	G4LogicalVolume* LQ1Logical = new G4LogicalVolume(Q1Solid,
		mMaterialManager->siliconsteel,"LQ1Logical",0,0,0);
	G4LogicalVolume* RQ1Logical = new G4LogicalVolume(Q1Solid,
		mMaterialManager->siliconsteel,"RQ1Logical",0,0,0);

	LQ1Logical->SetVisAttributes(IronVisAtt); 
	RQ1Logical->SetVisAttributes(IronVisAtt); 

	if(mSetupLHRS>=2)
	{
		//put it in the container, which also center at the hall center
		//therefore only the z_at_lab position need to be considered
		double pLQ1Pos_Z=(pHallCenter2LQ1Face+pQ1Length/2.0);
		new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pLQ1Pos_Z,0),
			LQ1Logical,"LQ1Phys",LHRSContainerLogical,0,0,0);
	}
	if(mSetupRHRS>=2)
	{
		double pRQ1Pos_Z=(pHallCenter2RQ1Face+pQ1Length/2.0);
		new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pRQ1Pos_Z,0),
			RQ1Logical,"RQ1Phys",RHRSContainerLogical,0,0,0);
	}


	/////////////////////////
	// HRS Q2 
	/////////////////////////
	double pHallCenter2LQ2Face=3.74*m;
	double pHallCenter2RQ2Face=3.74*m;
	double pQ2Rin=30.0*cm;
	double pQ2Rout=75.0*cm;
	double pQ2Length=180*cm;

	//By Jixie: HepRApp.jar has a bug to view this tubs. I am sure this tube is oriented correctly
	G4VSolid* Q2Solid = new G4Tubs("Q2Tub",pQ2Rin,pQ2Rout,pQ2Length/2.0,0.0,360.0*deg);

	G4LogicalVolume* LQ2Logical = new G4LogicalVolume(Q2Solid,
		mMaterialManager->siliconsteel,"LQ2Logical",0,0,0);
	G4LogicalVolume* RQ2Logical = new G4LogicalVolume(Q2Solid,
		mMaterialManager->siliconsteel,"RQ2Logical",0,0,0);

	LQ2Logical->SetVisAttributes(IronVisAtt); 
	RQ2Logical->SetVisAttributes(IronVisAtt); 

	if(mSetupLHRS>=3)
	{
		//put it in the container, which also center at the hall center
		double pLQ2Pos_Z=(pHallCenter2LQ2Face+pQ2Length/2.0);
		new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pLQ2Pos_Z,0),
			LQ2Logical,"LQ2Phys",LHRSContainerLogical,0,0,0);
	}
	if(mSetupRHRS>=3)
	{
		double pRQ2Pos_Z=(pHallCenter2RQ2Face+pQ2Length/2.0);
		new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pRQ2Pos_Z,0),
			RQ2Logical,"RQ2Phys",RHRSContainerLogical,0,0,0);
	}

	/////////////////////////
	// HRS Dipole 
	/////////////////////////
	//The dipole is built as a disc subtraced subtract another disc to get the side 
	//then subtract by the tunnel disc
	double pDipoleBendAngle=45*deg, pDipoleFaceAngle=30*deg;
	double pDipoleR=8.4*m;
	double pDipoleRprime=pDipoleR*sin(pDipoleBendAngle/2)/
		sin(180*deg-pDipoleBendAngle/2-pDipoleFaceAngle);
	//double pDipoleRprime=4.0518*m;

	double pDipoleRCenterY=pDipoleR;
	double pDipoleRCenterZ=9.96*m;

	//double pDipoleRprimeCenterY=pDipoleRprime*cos(pDipoleFaceAngle);			// =2.865*m;
	//double pDipoleRprimeCenterZ=9.96*m+pDipoleRprime*sin(pDipoleFaceAngle);	//12.825*m;

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

	G4LogicalVolume* LDipoleLogical = new G4LogicalVolume(DipoleSolid,
		mMaterialManager->siliconsteel,"LDipoleLogical",0,0,0);
	G4LogicalVolume* RDipoleLogical = new G4LogicalVolume(DipoleSolid,
		mMaterialManager->siliconsteel,"RDipoleLogical",0,0,0);

	LDipoleLogical->SetVisAttributes(OrangeVisAtt); 
	RDipoleLogical->SetVisAttributes(OrangeVisAtt); 

	G4RotationMatrix *pRotDipoleInContainer=new G4RotationMatrix();
	pRotDipoleInContainer->rotateY(90*deg); 
	if(mSetupLHRS>=4)
	{
		//put it in the container, which also center at the hall center
		new G4PVPlacement(pRotDipoleInContainer,
			G4ThreeVector(0,-pDipoleRCenterZ,pDipoleRCenterY),
			LDipoleLogical,"LDipolePhys",LHRSContainerLogical,0,0,0);
	}
	if(mSetupRHRS>=4)
	{
		new G4PVPlacement(pRotDipoleInContainer,
			G4ThreeVector(0,-pDipoleRCenterZ,pDipoleRCenterY),
			RDipoleLogical,"RDipolePhys",RHRSContainerLogical,0,0,0);
	}


	/////////////////////////
	// HRS Q3 
	/////////////////////////

	double pQ3CenterY=pDipoleR*(1-cos(pDipoleBendAngle))+2.4*m*sin(pDipoleBendAngle);
	double pQ3CenterZ=9.96*m+pDipoleR*sin(pDipoleBendAngle)+2.4*m*cos(pDipoleBendAngle);
	double pQ3Rin=30.0*cm;
	double pQ3Rout=75.0*cm;
	double pQ3Length=180*cm;

	//By Jixie: HepRApp.jar has a bug to view this tubs. I am sure this tube is oriented correctly
	G4VSolid* Q3Solid = new G4Tubs("Q3Tub",pQ3Rin,pQ3Rout,pQ3Length/2.0,0.0,360.0*deg);

	G4LogicalVolume* LQ3Logical = new G4LogicalVolume(Q3Solid,
		mMaterialManager->siliconsteel,"LQ3Logical",0,0,0);
	G4LogicalVolume* RQ3Logical = new G4LogicalVolume(Q3Solid,
		mMaterialManager->siliconsteel,"RQ3Logical",0,0,0);

	LQ3Logical->SetVisAttributes(IronVisAtt); 
	RQ3Logical->SetVisAttributes(IronVisAtt); 

	G4RotationMatrix *pRotQ3InContainer=new G4RotationMatrix();
	pRotQ3InContainer->rotateX(-45*deg); 
	if(mSetupLHRS>=4)
	{
		//put it in the container, which also center at the hall center
		new G4PVPlacement(pRotQ3InContainer,G4ThreeVector(0,-pQ3CenterZ,pQ3CenterY),
			LQ3Logical,"LQ3Phys",LHRSContainerLogical,0,0,0);
	}
	if(mSetupRHRS>=4)
	{
		new G4PVPlacement(pRotQ3InContainer,G4ThreeVector(0,-pQ3CenterZ,pQ3CenterY),
			RQ3Logical,"RQ3Phys",RHRSContainerLogical,0,0,0);
	}

	//////////////////////////////////////////////////////////
	//#ifdef G4DEBUG_GEOMETRY
	//plot the tunnel to view
	//////////////////////////
	//Q1Q2 tunnel
	//////////////////////////

	const int kNPlane_Q1Q2Tunnel=4;
	double rInner_Q1Q2Tunnel[] = {0,0,0,0};
	double rOuter_Q1Q2Tunnel[] = {pQ1Rin-0.1*mm,pQ1Rin-0.1*mm,pQ2Rin-0.1*mm,pQ2Rin-0.1*mm};
	double zPlane_Q1Q2Tunnel[] = {pHRSContainerRin,pHallCenter2RQ1Face+pQ1Length+10.0*cm,
		pHallCenter2RQ1Face+pQ1Length+30*cm,9.7*m};
	G4Polycone* Q1Q2TunnelSolid = new G4Polycone("Q1Q2TunnelPolycone",0.0,360.0*deg,
		kNPlane_Q1Q2Tunnel,zPlane_Q1Q2Tunnel,rInner_Q1Q2Tunnel,rOuter_Q1Q2Tunnel);

	G4LogicalVolume* LQ1Q2TunnelLogical = new G4LogicalVolume(Q1Q2TunnelSolid,
		mMaterialManager->vacuum,"LQ1Q2TunnelLogical",0,0,0);
	G4LogicalVolume* RQ1Q2TunnelLogical = new G4LogicalVolume(Q1Q2TunnelSolid,
		mMaterialManager->vacuum,"RQ1Q2TunnelLogical",0,0,0);

	LQ1Q2TunnelLogical->SetVisAttributes(YellowVisAtt); 
	RQ1Q2TunnelLogical->SetVisAttributes(YellowVisAtt); 

	if(mSetupLHRS>=3)
	{
		//put it in the container, which also center at the hall center
		new G4PVPlacement(pRotX90deg,G4ThreeVector(0,0,0),
			LQ1Q2TunnelLogical,"LQ1Q2TunnelPhys",LHRSContainerLogical,0,0,0);
	}
	if(mSetupRHRS>=3)
	{
		new G4PVPlacement(pRotX90deg,G4ThreeVector(0,0,0),
			RQ1Q2TunnelLogical,"RQ1Q2TunnelPhys",RHRSContainerLogical,0,0,0);
	}


	//////////////////////////
	//Dipole tunnel
	//////////////////////////
	//G4VSolid* DipoleTunnelSolid = new G4Tubs("DipoleTunnelSolid",pDipoleR-0.4*m+0.1*mm,
	//	pDipoleR+0.4*m-0.1*mm,0.125*m-0.1*mm,180*deg,pDipoleBendAngle);

	G4Polycone* DipoleTunnelSolid = new G4Polycone("DipoleTunnelSolid",
		180.0*deg,pDipoleBendAngle,
		kNPlane_DipoleTunnel,zPlane_DipoleTunnel,rInner_DipoleTunnel,rOuter_DipoleTunnel);

	G4LogicalVolume* LDipoleTunnelLogical = new G4LogicalVolume(DipoleTunnelSolid,
		mMaterialManager->vacuum,"LDipoleTunnelLogical",0,0,0);
	G4LogicalVolume* RDipoleTunnelLogical = new G4LogicalVolume(DipoleTunnelSolid,
		mMaterialManager->vacuum,"RDipoleTunnelLogical",0,0,0);

	LDipoleTunnelLogical->SetVisAttributes(YellowVisAtt); 
	RDipoleTunnelLogical->SetVisAttributes(YellowVisAtt); 

	if(mSetupLHRS>=4)
	{
		//put it in the container, which also center at the hall center
		new G4PVPlacement(pRotDipoleInContainer,
			G4ThreeVector(0,-pDipoleRCenterZ,pDipoleRCenterY),
			LDipoleTunnelLogical,"LDipoleTunnelPhys",LHRSContainerLogical,0,0,0);
	}
	if(mSetupRHRS>=4)
	{
		new G4PVPlacement(pRotDipoleInContainer,
			G4ThreeVector(0,-pDipoleRCenterZ,pDipoleRCenterY),
			RDipoleTunnelLogical,"RDipoleTunnelPhys",RHRSContainerLogical,0,0,0);
	}


	//////////////////////////
	//Q3 tunnel
	//////////////////////////

	G4VSolid* Q3TunnelSolid = new G4Tubs("Q3TunnelTub",0,pQ3Rin-0.1*mm,
		1.8*m,0.0,360.0*deg);

	G4LogicalVolume* LQ3TunnelLogical = new G4LogicalVolume(Q3TunnelSolid,
		mMaterialManager->vacuum,"LQ3TunnelLogical",0,0,0);
	G4LogicalVolume* RQ3TunnelLogical = new G4LogicalVolume(Q3TunnelSolid,
		mMaterialManager->vacuum,"RQ3TunnelLogical",0,0,0);

	LQ3TunnelLogical->SetVisAttributes(YellowVisAtt); 
	RQ3TunnelLogical->SetVisAttributes(YellowVisAtt); 

	double pQ3TunnelPos_Y=pDipoleR*(1-cos(pDipoleBendAngle))+2.0*m*sin(pDipoleBendAngle);
	double pQ3TunnelPos_Z=9.96*m+pDipoleR*sin(pDipoleBendAngle)+2.0*m*cos(pDipoleBendAngle);
	if(mSetupLHRS>=4)
	{
		//put it in the container, which also center at the hall center
		new G4PVPlacement(pRotQ3InContainer,G4ThreeVector(0,-pQ3TunnelPos_Z,pQ3TunnelPos_Y),
			LQ3TunnelLogical,"LQ1Q2TunnelPhys",LHRSContainerLogical,0,0,0);
	}
	if(mSetupRHRS>=4)
	{
		new G4PVPlacement(pRotQ3InContainer,G4ThreeVector(0,-pQ3TunnelPos_Z,pQ3TunnelPos_Y),
			RQ3TunnelLogical,"RQ1Q2TunnelPhys",RHRSContainerLogical,0,0,0);
	}

	//#endif
	//////////////////////////////////////////////////////////

	return theG2PHRSPhys;
}

/////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* G2PDetectorConstruction::ConstructG2PHRS(G4LogicalVolume* motherLogical)
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
	//set up the septum only if there is a angle difference	
	bool mSetupSeptumBlock=((mLHRSAngle-mLSeptumAngle)/deg>0.5)?true:false;

	/////////////////////////////////////////////////
	//  from Hall A NIM Paper, the standard sieve slit
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
	G4VPhysicalVolume* theG2PHRSPhys=0;


	///////////////////////////////////////
	//Sieve Slit for HRS-Angle=5.69
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

	//the big center hole horizontal and vertical offset, From Alan
	double pSieveSlitLargeHoleH=2.967*mm;    //positive means shift away from the beam line
	double pSieveSlitLargeHoleV=2.606*mm;

	//Updated By Jixie@20121129, change these number with survey results
	//left and right sieve are almost equal, so I just use the left sieve survey results
	//These numbers are from the survey,  
	pSieveSlitLargeHoleH=80.35*mm-799.6*mm*sin(mLSeptumAngle);   //1.63 mm for 5.65 degree, or 1.08 mm for
	pSieveSlitLargeHoleV=0.31*mm;

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
		//located at z=686.46 mm, no rotation
		double pSeptumX=140.0*cm;
		double pSeptumY=84.4*cm;
		double pSeptumZ=74.0*cm;
		double pSeptumTunnelX=30.4*cm;
		double pSeptumTunnelY=24.4*cm-2.0*inch;  //By Jixie @20120205: Add 2 inches of iron
		double pSeptumBeamHoleX=7.8*cm;
		double pSeptumBeamHoleY=8.0*cm;

		double pSeptumTunnelPos_X=8.4*cm+pSeptumTunnelX/2.0;
		double pSeptumPos_Z=68.646*cm;		
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
		G4RotationMatrix* pSeptumCoilRotBackDwon = new G4RotationMatrix();
		pSeptumCoilRotBackDwon->rotateY(270*deg);
		G4RotationMatrix* pSeptumCoilRotBackUp = new G4RotationMatrix();
		pSeptumCoilRotBackUp->rotateY(270*deg);
		pSeptumCoilRotBackUp->rotateX(180*deg);

		new G4PVPlacement(pSeptumCoilRotBackDwon,
			G4ThreeVector(pSeptumCoilPos_X_in,pSeptumCoilPos_Y,pSeptumCoilPos_Z_down),
			septumCoilLogical,"septumCoilPhys",motherLogical,true,8,0);
		new G4PVPlacement(pSeptumCoilRotBackUp,
			G4ThreeVector(pSeptumCoilPos_X_in,-pSeptumCoilPos_Y,pSeptumCoilPos_Z_down),
			septumCoilLogical,"septumCoilPhys",motherLogical,true,9,0);
		new G4PVPlacement(pSeptumCoilRotBackUp,
			G4ThreeVector(pSeptumCoilPos_X_out,-pSeptumCoilPos_Y,pSeptumCoilPos_Z_down),
			septumCoilLogical,"septumCoilPhys",motherLogical,true,10,0);
		new G4PVPlacement(pSeptumCoilRotBackDwon,
			G4ThreeVector(pSeptumCoilPos_X_out,pSeptumCoilPos_Y,pSeptumCoilPos_Z_down),
			septumCoilLogical,"septumCoilPhys",motherLogical,true,11,0);

		new G4PVPlacement(pSeptumCoilRotBackDwon,
			G4ThreeVector(-pSeptumCoilPos_X_out,pSeptumCoilPos_Y,pSeptumCoilPos_Z_down),
			septumCoilLogical,"septumCoilPhys",motherLogical,true,12,0);
		new G4PVPlacement(pSeptumCoilRotBackUp,
			G4ThreeVector(-pSeptumCoilPos_X_out,-pSeptumCoilPos_Y,pSeptumCoilPos_Z_down),
			septumCoilLogical,"septumCoilPhys",motherLogical,true,13,0);
		new G4PVPlacement(pSeptumCoilRotBackUp,
			G4ThreeVector(-pSeptumCoilPos_X_in,-pSeptumCoilPos_Y,pSeptumCoilPos_Z_down),
			septumCoilLogical,"septumCoilPhys",motherLogical,true,14,0);
		new G4PVPlacement(pSeptumCoilRotBackDwon,
			G4ThreeVector(-pSeptumCoilPos_X_in,pSeptumCoilPos_Y,pSeptumCoilPos_Z_down),
			septumCoilLogical,"septumCoilPhys",motherLogical,true,15,0);
	}

	//Build the HRS, positions were taken fron the NIM paper
	//The Apertures are defined as:
	//Q1 exit is a circle of radius 0.15m
	//Q3 entrance and exit are circles of radius 0.30m
	//The dipole entrance and exit are trapezoids:
	//-0.40 m<x<0.40 m (+x is down) 
	//y=+-(0.125*(1-(1.25*x/8.40)) (smallest gap at positive x, x in m)

	/////////////////////////
	// HRS QQDQ Containner
	/////////////////////////
	//Build a container using polycone, covering 270+/-12 degrees, 1.46m to 20.76m
	//3.05m below the beam line is the ground, need to subtract everything beneath that
	//looks like this:
/*
    HRS container:covering 24 degrees in X-Z plane, 10.05 m height.
    Stuff inside this containner will also use the pivot as the origin, but do not 
    need to worry about the rotation of the HRS.
                        7.0m above beam line 
                       -----------------------------| 20.76m from pivot,  
                      /                             |
                     /                        Q3    |
                    /                               |
                   /                       E        |
                  /                     L           |
          --------                   O              |
         /                        P                 |
    -----                     I                     |
----|----- Q1 --- Q2 ---- D     --------------------|------ beam line -----
    |                                               |   
    ------------------------------------------------|
    1.46m from povot, 3.05 m below beam line 
*/

	double pHRSContainerRin=1.46*m,pHRSContainerRout=20.76*m;
	double pBeamLine2Ground;
	pBeamLine2Ground=-3.05*m;
	//build the container with polycone

	const int kNPlane_HRSContainer=7;
	double rInner_HRSContainer[] = {pHRSContainerRin,pHRSContainerRin,2.5*m,
		3.7*m,9.0*m,pHRSContainerRout-3.0*m,pHRSContainerRout};
	double rOuter_HRSContainer[] = {pHRSContainerRout,pHRSContainerRout,pHRSContainerRout,
		pHRSContainerRout,pHRSContainerRout,pHRSContainerRout,pHRSContainerRout};
	double zPlane_HRSContainer[] = {-2.0*m,1.0*m,1.0*m,
		2.0*m,2.0*m,7.0*m,7.0*m};
	G4Polycone* HRSContainerSolid = new G4Polycone("HRSContainer",258.0*deg,24.0*deg,
		kNPlane_HRSContainer,zPlane_HRSContainer,rInner_HRSContainer,rOuter_HRSContainer);

	////build the container using tube
	//double pHRSContainerHeight=14.0*m;
	//G4VSolid* HRSContainerTub = new G4Tubs("HRSContainerTub",pHRSContainerRin,
	//	pHRSContainerRout,pHRSContainerHeight/2.0,258.0*deg,24.0*deg);
	////ground
	//double pGoundHeight=10*m;
	//G4VSolid* GroundTub = new G4Tubs("GroundTub",0,30*m,pGoundHeight/2,0*deg,360.0*deg);
	////tube subtract the ground
	//G4SubtractionSolid* HRSContainerSolid = new G4SubtractionSolid("HRSContainer",
	//		HRSContainerTub,GroundTub,
	//		0,G4ThreeVector(0,0,pBeamLine2Ground-pGoundHeight/2));

	G4LogicalVolume* LHRSContainerLogical = new G4LogicalVolume(HRSContainerSolid,
		mMaterialManager->vacuum,"LHRSContainerLogical",0,0,0);
	G4LogicalVolume* RHRSContainerLogical = new G4LogicalVolume(HRSContainerSolid,
		mMaterialManager->vacuum,"RHRSContainerLogical",0,0,0);

	LHRSContainerLogical->SetVisAttributes(HallVisAtt); 
	RHRSContainerLogical->SetVisAttributes(HallVisAtt); 

	G4RotationMatrix *pRotLHRSContainer=new G4RotationMatrix();
	pRotLHRSContainer->rotateX(-270*deg);
	pRotLHRSContainer->rotateZ(-mLHRSAngle);  

	G4RotationMatrix *pRotRHRSContainer=new G4RotationMatrix();
	pRotRHRSContainer->rotateX(-270*deg);
	pRotRHRSContainer->rotateZ(-mRHRSAngle);  

	if(mSetupLHRS>=2)
	{
		new G4PVPlacement(pRotLHRSContainer,G4ThreeVector(0,0,0),
			LHRSContainerLogical,"LHRSContainerPhys",motherLogical,0,0,0);
	}
	if(mSetupRHRS>=2)
	{
		new G4PVPlacement(pRotRHRSContainer,G4ThreeVector(0,0,0),
			RHRSContainerLogical,"RHRSContainerPhys",motherLogical,0,0,0);
	}


	/////////////////////////
	// HRS Q1 
	/////////////////////////
	double pHallCenter2LQ1Face=1.69*m;
	double pHallCenter2RQ1Face=1.69*m;
	double pQ1Rin=15.0*cm;
	double pQ1Rout=35.0*cm;
	double pQ1Length=80*cm;
	//double pQ1Length=(1.698M-1.36*m)*2+0.8*m;

	//By Jixie: HepRApp.jar has a bug to view this tubs. I am sure this tube is oriented correctly
	G4VSolid* Q1Solid = new G4Tubs("Q1Tub",pQ1Rin,pQ1Rout,pQ1Length/2.0,0.0,360.0*deg);

	//build 2 copy since there are differennt fields involved in it
	G4LogicalVolume* LQ1Logical = new G4LogicalVolume(Q1Solid,
		mMaterialManager->siliconsteel,"LQ1Logical",0,0,0);
	G4LogicalVolume* RQ1Logical = new G4LogicalVolume(Q1Solid,
		mMaterialManager->siliconsteel,"RQ1Logical",0,0,0);

	LQ1Logical->SetVisAttributes(IronVisAtt); 
	RQ1Logical->SetVisAttributes(IronVisAtt); 

	if(mSetupLHRS>=2)
	{
		//put it in the container, which also center at the hall center
		//therefore only the z_at_lab position need to be considered
		double pLQ1Pos_Z=(pHallCenter2LQ1Face+pQ1Length/2.0);
		new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pLQ1Pos_Z,0),
			LQ1Logical,"LQ1Phys",LHRSContainerLogical,0,0,0);
	}
	if(mSetupRHRS>=2)
	{
		double pRQ1Pos_Z=(pHallCenter2RQ1Face+pQ1Length/2.0);
		new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pRQ1Pos_Z,0),
			RQ1Logical,"RQ1Phys",RHRSContainerLogical,0,0,0);
	}


	/////////////////////////
	// HRS Q2 
	/////////////////////////
	double pHallCenter2LQ2Face=3.74*m;
	double pHallCenter2RQ2Face=3.74*m;
	double pQ2Rin=30.0*cm;
	double pQ2Rout=75.0*cm;
	double pQ2Length=180*cm;

	//By Jixie: HepRApp.jar has a bug to view this tubs. I am sure this tube is oriented correctly
	G4VSolid* Q2Solid = new G4Tubs("Q2Tub",pQ2Rin,pQ2Rout,pQ2Length/2.0,0.0,360.0*deg);

	G4LogicalVolume* LQ2Logical = new G4LogicalVolume(Q2Solid,
		mMaterialManager->siliconsteel,"LQ2Logical",0,0,0);
	G4LogicalVolume* RQ2Logical = new G4LogicalVolume(Q2Solid,
		mMaterialManager->siliconsteel,"RQ2Logical",0,0,0);

	LQ2Logical->SetVisAttributes(IronVisAtt); 
	RQ2Logical->SetVisAttributes(IronVisAtt); 

	if(mSetupLHRS>=3)
	{
		//put it in the container, which also center at the hall center
		double pLQ2Pos_Z=(pHallCenter2LQ2Face+pQ2Length/2.0);
		new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pLQ2Pos_Z,0),
			LQ2Logical,"LQ2Phys",LHRSContainerLogical,0,0,0);
	}
	if(mSetupRHRS>=3)
	{
		double pRQ2Pos_Z=(pHallCenter2RQ2Face+pQ2Length/2.0);
		new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pRQ2Pos_Z,0),
			RQ2Logical,"RQ2Phys",RHRSContainerLogical,0,0,0);
	}

	/////////////////////////
	// HRS Dipole 
	/////////////////////////
	//The dipole is built as a disc subtraced subtract another disc to get the side 
	//then subtract by the tunnel disc
	double pDipoleBendAngle=45*deg, pDipoleFaceAngle=30*deg;
	double pDipoleR=8.4*m;
	double pDipoleRprime=pDipoleR*sin(pDipoleBendAngle/2)/
		sin(180*deg-pDipoleBendAngle/2-pDipoleFaceAngle);
	//double pDipoleRprime=4.0518*m;

	double pDipoleRCenterY=pDipoleR;
	double pDipoleRCenterZ=9.96*m;

	//double pDipoleRprimeCenterY=pDipoleRprime*cos(pDipoleFaceAngle);			// =2.865*m;
	//double pDipoleRprimeCenterZ=9.96*m+pDipoleRprime*sin(pDipoleFaceAngle);	//12.825*m;

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

	G4LogicalVolume* LDipoleLogical = new G4LogicalVolume(DipoleSolid,
		mMaterialManager->siliconsteel,"LDipoleLogical",0,0,0);
	G4LogicalVolume* RDipoleLogical = new G4LogicalVolume(DipoleSolid,
		mMaterialManager->siliconsteel,"RDipoleLogical",0,0,0);

	LDipoleLogical->SetVisAttributes(OrangeVisAtt); 
	RDipoleLogical->SetVisAttributes(OrangeVisAtt); 

	G4RotationMatrix *pRotDipoleInContainer=new G4RotationMatrix();
	pRotDipoleInContainer->rotateY(90*deg); 
	if(mSetupLHRS>=4)
	{
		//put it in the container, which also center at the hall center
		new G4PVPlacement(pRotDipoleInContainer,
			G4ThreeVector(0,-pDipoleRCenterZ,pDipoleRCenterY),
			LDipoleLogical,"LDipolePhys",LHRSContainerLogical,0,0,0);
	}
	if(mSetupRHRS>=4)
	{
		new G4PVPlacement(pRotDipoleInContainer,
			G4ThreeVector(0,-pDipoleRCenterZ,pDipoleRCenterY),
			RDipoleLogical,"RDipolePhys",RHRSContainerLogical,0,0,0);
	}


	/////////////////////////
	// HRS Q3 
	/////////////////////////

	double pQ3CenterY=pDipoleR*(1-cos(pDipoleBendAngle))+2.4*m*sin(pDipoleBendAngle);
	double pQ3CenterZ=9.96*m+pDipoleR*sin(pDipoleBendAngle)+2.4*m*cos(pDipoleBendAngle);
	double pQ3Rin=30.0*cm;
	double pQ3Rout=75.0*cm;
	double pQ3Length=180*cm;

	//By Jixie: HepRApp.jar has a bug to view this tubs. I am sure this tube is oriented correctly
	G4VSolid* Q3Solid = new G4Tubs("Q3Tub",pQ3Rin,pQ3Rout,pQ3Length/2.0,0.0,360.0*deg);

	G4LogicalVolume* LQ3Logical = new G4LogicalVolume(Q3Solid,
		mMaterialManager->siliconsteel,"LQ3Logical",0,0,0);
	G4LogicalVolume* RQ3Logical = new G4LogicalVolume(Q3Solid,
		mMaterialManager->siliconsteel,"RQ3Logical",0,0,0);

	LQ3Logical->SetVisAttributes(IronVisAtt); 
	RQ3Logical->SetVisAttributes(IronVisAtt); 

	G4RotationMatrix *pRotQ3InContainer=new G4RotationMatrix();
	pRotQ3InContainer->rotateX(-45*deg); 
	if(mSetupLHRS>=4)
	{
		//put it in the container, which also center at the hall center
		new G4PVPlacement(pRotQ3InContainer,G4ThreeVector(0,-pQ3CenterZ,pQ3CenterY),
			LQ3Logical,"LQ3Phys",LHRSContainerLogical,0,0,0);
	}
	if(mSetupRHRS>=4)
	{
		new G4PVPlacement(pRotQ3InContainer,G4ThreeVector(0,-pQ3CenterZ,pQ3CenterY),
			RQ3Logical,"RQ3Phys",RHRSContainerLogical,0,0,0);
	}

	//////////////////////////////////////////////////////////
	//#ifdef G4DEBUG_GEOMETRY
	//plot the tunnel to view
	//////////////////////////
	//Q1Q2 tunnel
	//////////////////////////

	const int kNPlane_Q1Q2Tunnel=4;
	double rInner_Q1Q2Tunnel[] = {0,0,0,0};
	double rOuter_Q1Q2Tunnel[] = {pQ1Rin-0.1*mm,pQ1Rin-0.1*mm,pQ2Rin-0.1*mm,pQ2Rin-0.1*mm};
	double zPlane_Q1Q2Tunnel[] = {pHRSContainerRin,pHallCenter2RQ1Face+pQ1Length+10.0*cm,
		pHallCenter2RQ1Face+pQ1Length+30*cm,9.7*m};
	G4Polycone* Q1Q2TunnelSolid = new G4Polycone("Q1Q2TunnelPolycone",0.0,360.0*deg,
		kNPlane_Q1Q2Tunnel,zPlane_Q1Q2Tunnel,rInner_Q1Q2Tunnel,rOuter_Q1Q2Tunnel);

	G4LogicalVolume* LQ1Q2TunnelLogical = new G4LogicalVolume(Q1Q2TunnelSolid,
		mMaterialManager->vacuum,"LQ1Q2TunnelLogical",0,0,0);
	G4LogicalVolume* RQ1Q2TunnelLogical = new G4LogicalVolume(Q1Q2TunnelSolid,
		mMaterialManager->vacuum,"RQ1Q2TunnelLogical",0,0,0);

	LQ1Q2TunnelLogical->SetVisAttributes(YellowVisAtt); 
	RQ1Q2TunnelLogical->SetVisAttributes(YellowVisAtt); 

	if(mSetupLHRS>=3)
	{
		//put it in the container, which also center at the hall center
		new G4PVPlacement(pRotX90deg,G4ThreeVector(0,0,0),
			LQ1Q2TunnelLogical,"LQ1Q2TunnelPhys",LHRSContainerLogical,0,0,0);
	}
	if(mSetupRHRS>=3)
	{
		new G4PVPlacement(pRotX90deg,G4ThreeVector(0,0,0),
			RQ1Q2TunnelLogical,"RQ1Q2TunnelPhys",RHRSContainerLogical,0,0,0);
	}


	//////////////////////////
	//Dipole tunnel
	//////////////////////////
	//G4VSolid* DipoleTunnelSolid = new G4Tubs("DipoleTunnelSolid",pDipoleR-0.4*m+0.1*mm,
	//	pDipoleR+0.4*m-0.1*mm,0.125*m-0.1*mm,180*deg,pDipoleBendAngle);

	G4Polycone* DipoleTunnelSolid = new G4Polycone("DipoleTunnelSolid",
		180.0*deg,pDipoleBendAngle,
		kNPlane_DipoleTunnel,zPlane_DipoleTunnel,rInner_DipoleTunnel,rOuter_DipoleTunnel);

	G4LogicalVolume* LDipoleTunnelLogical = new G4LogicalVolume(DipoleTunnelSolid,
		mMaterialManager->vacuum,"LDipoleTunnelLogical",0,0,0);
	G4LogicalVolume* RDipoleTunnelLogical = new G4LogicalVolume(DipoleTunnelSolid,
		mMaterialManager->vacuum,"RDipoleTunnelLogical",0,0,0);

	LDipoleTunnelLogical->SetVisAttributes(YellowVisAtt); 
	RDipoleTunnelLogical->SetVisAttributes(YellowVisAtt); 

	if(mSetupLHRS>=4)
	{
		//put it in the container, which also center at the hall center
		new G4PVPlacement(pRotDipoleInContainer,
			G4ThreeVector(0,-pDipoleRCenterZ,pDipoleRCenterY),
			LDipoleTunnelLogical,"LDipoleTunnelPhys",LHRSContainerLogical,0,0,0);
	}
	if(mSetupRHRS>=4)
	{
		new G4PVPlacement(pRotDipoleInContainer,
			G4ThreeVector(0,-pDipoleRCenterZ,pDipoleRCenterY),
			RDipoleTunnelLogical,"RDipoleTunnelPhys",RHRSContainerLogical,0,0,0);
	}


	//////////////////////////
	//Q3 tunnel
	//////////////////////////

	G4VSolid* Q3TunnelSolid = new G4Tubs("Q3TunnelTub",0,pQ3Rin-0.1*mm,
		1.8*m,0.0,360.0*deg);

	G4LogicalVolume* LQ3TunnelLogical = new G4LogicalVolume(Q3TunnelSolid,
		mMaterialManager->vacuum,"LQ3TunnelLogical",0,0,0);
	G4LogicalVolume* RQ3TunnelLogical = new G4LogicalVolume(Q3TunnelSolid,
		mMaterialManager->vacuum,"RQ3TunnelLogical",0,0,0);

	LQ3TunnelLogical->SetVisAttributes(YellowVisAtt); 
	RQ3TunnelLogical->SetVisAttributes(YellowVisAtt); 

	double pQ3TunnelPos_Y=pDipoleR*(1-cos(pDipoleBendAngle))+2.0*m*sin(pDipoleBendAngle);
	double pQ3TunnelPos_Z=9.96*m+pDipoleR*sin(pDipoleBendAngle)+2.0*m*cos(pDipoleBendAngle);
	if(mSetupLHRS>=4)
	{
		//put it in the container, which also center at the hall center
		new G4PVPlacement(pRotQ3InContainer,G4ThreeVector(0,-pQ3TunnelPos_Z,pQ3TunnelPos_Y),
			LQ3TunnelLogical,"LQ1Q2TunnelPhys",LHRSContainerLogical,0,0,0);
	}
	if(mSetupRHRS>=4)
	{
		new G4PVPlacement(pRotQ3InContainer,G4ThreeVector(0,-pQ3TunnelPos_Z,pQ3TunnelPos_Y),
			RQ3TunnelLogical,"RQ1Q2TunnelPhys",RHRSContainerLogical,0,0,0);
	}


	//#endif
	//////////////////////////////////////////////////////////



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

		double pHallCenter2VB=1.40*m; //pHRSContainerRin-6*cm;
		gConfig->SetParameter("Pivot2LHRSVBFace",pHallCenter2VB-mPivotZOffset*cos(mLHRSAngle));
		gConfig->SetParameter("Pivot2RHRSVBFace",pHallCenter2VB-mPivotZOffset*cos(mRHRSAngle)); 

		if(mSetupLHRS)
		{
			double pLHRSVBPos_X=(pHallCenter2VB-mHRSVBThick/2)*
				sin(mLHRSAngle)+mPivotXOffset;
			double pLHRSVBPos_Y=mPivotYOffset;
			//no need to correct for pivot since the distance is from the hall center
			double pLHRSVBPos_Z=(pHallCenter2VB-mHRSVBThick/2.0)*cos(mLHRSAngle);
			new G4PVPlacement(pRotLHRS,
				G4ThreeVector(pLHRSVBPos_X,pLHRSVBPos_Y,pLHRSVBPos_Z),
				HRSVBLogical,"virtualBoundaryPhys_LHRS",motherLogical,0,0,0);
		}
		if(mSetupRHRS)
		{
			double pRHRSVBPos_X=(pHallCenter2VB-mHRSVBThick/2)*
				sin(mRHRSAngle)+mPivotXOffset;
			double pRHRSVBPos_Y=mPivotYOffset;
			//no need to correct for pivot since the distance is from the hall center
			double pRHRSVBPos_Z=(pHallCenter2VB-mHRSVBThick/2)*cos(mRHRSAngle); 
			new G4PVPlacement(pRotRHRS,
				G4ThreeVector(pRHRSVBPos_X,pRHRSVBPos_Y,pRHRSVBPos_Z),
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


	return theG2PHRSPhys;
}



/////////////////////////////////////////////////////////////////////
//This is th eold routine which built sieve, septum, HRSS and HRS VB at a time
 G4VPhysicalVolume* G2PDetectorConstruction::ConstructG2PBeamDump(G4LogicalVolume* motherLogical)
{
	const double inch=2.54*cm;
	//The dump is 23.5cm (width) X 38cm(height) X 15 cm (thick), face2pivot=64cm

	//the top of local dump is 4.358" from the center of the dump, add extra height to make it 
	//symmetric, will chop the extra portion out when building the geometry
	double pBDExtraHeight=mBeamDumpHeight-2*4.358*inch;

	/////////////////////////////////////
	//plastic shielding
	/////////////////////////////////////
	bool pSetupBDPolyShield=true;
	if(pSetupBDPolyShield)
	{
		G4RotationMatrix *pRotBDPolyShield=new G4RotationMatrix();
		//pRotBDPolyShield->rotateY(-90*deg-pBDSeeThirdArmAngle); 
		pRotBDPolyShield->rotateY(-90*deg); 

		//block 1: Al design it as one small box sitting on the big one
		double pBDPolyShield1X=24*inch,pBDPolyShield1Z=12*inch;
		double pBDPolyShield1BottomY=24*inch,pBDPolyShield1TopY=10*inch;

		G4VSolid* BDPolyShield1TopSolid = new G4Box("BDPolyShield1TopBox",12*inch/2.0,
			pBDPolyShield1TopY/2.0,pBDPolyShield1Z/2.0);
		G4VSolid* BDPolyShield1BottomSolid = new G4Box("BDPolyShield1BottomBox",pBDPolyShield1X/2.0,
			pBDPolyShield1BottomY/2.0,pBDPolyShield1Z/2.0);
		G4UnionSolid* BDPolyShield1Solid = new G4UnionSolid("BDPolyShield1Box",
			BDPolyShield1BottomSolid,BDPolyShield1TopSolid,
			0,G4ThreeVector(0,pBDPolyShield1BottomY/2+pBDPolyShield1TopY/2,0));	

		//block 2: one big box subtract a triangle
		double pBDPolyShield2X=13.07*inch,pBDPolyShield2Y=20*inch,pBDPolyShield2Z=8*inch;
		G4VSolid* BDPolyShield2WholeSolid = new G4Box("BDPolyShield2WholeBox",pBDPolyShield2X/2.0,
			pBDPolyShield2Y/2.0,pBDPolyShield2Z/2.0);

		double pBDPolyShield2Y_short = 12*inch;
		double pBDPolyShield2SubX = 20*inch;
		double pBDPolyShield2SubY = 8*inch;
		G4VSolid* BDPolyShield2SubSolid = new G4Box("BDPolyShield2SubBox",pBDPolyShield2SubX/2.0,
			pBDPolyShield2SubY/2.0,pBDPolyShield2Z/2.0+1*mm);
		double pBDPolyShield2SubAngle = atan((pBDPolyShield2Y-pBDPolyShield2Y_short)/pBDPolyShield2X);
		G4RotationMatrix *pRotBDPolyShield2Sub = new G4RotationMatrix();
		pRotBDPolyShield2Sub->rotateZ(pBDPolyShield2SubAngle); 
		//center Y of the subtracted part
		double pBDPolyShield2SubPos_Y = pBDPolyShield2SubY/2 / cos(pBDPolyShield2SubAngle) + 
			(pBDPolyShield2Y/2-(pBDPolyShield2Y-pBDPolyShield2Y_short)/2);
		G4SubtractionSolid* BDPolyShield2Solid = new G4SubtractionSolid("BDPolyShield2Box",
			BDPolyShield2WholeSolid,BDPolyShield2SubSolid,
			pRotBDPolyShield2Sub,G4ThreeVector(0,pBDPolyShield2SubPos_Y,0));

		//block3, inside the aluminum box which holds helium 
		double pBDPolyShield3X=7*inch,pBDPolyShield3Y=14*inch,pBDPolyShield3Z=4*inch;
		G4VSolid* BDPolyShield3Solid = new G4Box("BDPolyShield3Box",pBDPolyShield3X/2.0,
			pBDPolyShield3Y/2.0,pBDPolyShield3Z/2.0);

		G4LogicalVolume* BDPolyShield1Logical = new G4LogicalVolume(BDPolyShield1Solid,
			mMaterialManager->boratedpoly05,"BDPolyShieldLogical",0,0,0);
		BDPolyShield1Logical->SetVisAttributes(LightBlueVisAtt);

		G4LogicalVolume* BDPolyShield2Logical = new G4LogicalVolume(BDPolyShield2Solid,
			mMaterialManager->boratedpoly05,"BDPolyShieldLogical",0,0,0);
		BDPolyShield2Logical->SetVisAttributes(LightBlueVisAtt);

		G4LogicalVolume* BDPolyShield3Logical = new G4LogicalVolume(BDPolyShield3Solid,
			mMaterialManager->boratedpoly05,"BDPolyShieldLogical",0,0,0);
		BDPolyShield3Logical->SetVisAttributes(LightBlueVisAtt);

		//By Jixie @20111208: Updated the position according to drawing A08027-0300-0400
		double pBDPolyShield1Pos_X =  37.741*inch + pBDPolyShield1Z/2.0 + mPivotXOffset;
		double pBDPolyShield1Pos_Y =  -12.0*inch + pBDPolyShield1BottomY/2.0 + mPivotYOffset;
		double pBDPolyShield1Pos_Z =  25.176*inch + pBDPolyShield1X/2.0 + mPivotZOffset;
		if(fabs(mLHRSAngle-mLSeptumAngle)<1.0E-3)  
		{
			//for 12.5 degrees
			pBDPolyShield1Pos_X =  25.228*inch + pBDPolyShield1Z/2.0 + mPivotXOffset;
			pBDPolyShield1Pos_Z =  25.076*inch + pBDPolyShield1X/2.0 + mPivotZOffset;
		}
		new G4PVPlacement(pRotBDPolyShield,
			G4ThreeVector(pBDPolyShield1Pos_X,pBDPolyShield1Pos_Y,pBDPolyShield1Pos_Z),
			BDPolyShield1Logical,"BDPolyShield1Phys",motherLogical,0,0);

		double pBDPolyShield2Pos_X =  12.39*inch + pBDPolyShield2Z/2.0 + mPivotXOffset;
		double pBDPolyShield2Pos_Y = -13.52*inch + pBDPolyShield2Y/2.0 + mPivotYOffset;
		double pBDPolyShield2Pos_Z = -17.42*inch+876.93*mm + pBDPolyShield2X/2.0 + mPivotZOffset;
		new G4PVPlacement(pRotBDPolyShield,
			G4ThreeVector(pBDPolyShield2Pos_X,pBDPolyShield2Pos_Y,pBDPolyShield2Pos_Z),
			BDPolyShield2Logical,"BDPolyShield2Phys",motherLogical,0,0);

		double pBDPolyShield3Pos_X =  7.93*inch + pBDPolyShield3Z/2.0 + mPivotXOffset;
		double pBDPolyShield3Pos_Y = -9.75*inch + pBDPolyShield3Y/2.0 + mPivotYOffset;
		double pBDPolyShield3Pos_Z = -13.94*inch+876.93*mm + pBDPolyShield3X/2.0 + mPivotZOffset;
		new G4PVPlacement(pRotBDPolyShield,
			G4ThreeVector(pBDPolyShield3Pos_X,pBDPolyShield3Pos_Y,pBDPolyShield3Pos_Z),
			BDPolyShield3Logical,"BDPolyShield3Phys",motherLogical,0,0);
	}

	/////////////////////////
	// beam dump Container
	/////////////////////////

	double pBDContainerPos_X=mPivotXOffset;
	double pBDContainerPos_Y=mPivotYOffset;
	double pBDContainerPos_Z=(mPivot2BeamDumpFace+mBeamDumpThick/2.0)+mPivotZOffset;

	G4VSolid* beamDumpContainerWholeSolid = new G4Box("beamDumpContainerWholeBox",
		mBeamDumpWidth/2.0+0.1*mm,(mBeamDumpHeight+pBDExtraHeight)/2.0+0.1*mm,
		mBeamDumpThick/2.0+0.1*mm);
	G4VSolid* beamDumpContainerExtraSolid = new G4Box("beamDumpContainerBox",
		mBeamDumpWidth/2.0+0.2*mm,pBDExtraHeight/2.0,mBeamDumpThick/2.0+0.2*mm);
	//double pBDContainerExtraPos_Y=(mBeamDumpHeight+pBDExtraHeight)/2+0.1*mm-pBDExtraHeight/2.0;
	double pBDContainerExtraPos_Y=(mBeamDumpHeight+pBDExtraHeight)/2+0.1*mm-pBDExtraHeight/2.0;
	G4SubtractionSolid* beamDumpContainerSolid = new G4SubtractionSolid(
		"beamDumpContainerSolid",beamDumpContainerWholeSolid,beamDumpContainerExtraSolid,0,
		G4ThreeVector(0,pBDContainerExtraPos_Y,0));


	G4LogicalVolume* beamDumpContainerLogical = new G4LogicalVolume(beamDumpContainerSolid,
		mMaterialManager->heliumGas,"beamDumpContainerLogical",0,0,0);

	G4VPhysicalVolume* beamDumpPhysical=new G4PVPlacement(0,
		G4ThreeVector(pBDContainerPos_X,pBDContainerPos_Y,pBDContainerPos_Z),
		beamDumpContainerLogical,"beamDumpContainerPhys",motherLogical,0,0);
	beamDumpContainerLogical->SetVisAttributes(HallVisAtt);  

	/////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////
	// beam dump aluminum shell 
	/////////////////////////
	bool pSetupBDAlBox=true;
	if(pSetupBDAlBox)
	{
		//this box was supposed to hold the helium gas
		// left and right panels are placed here
		double pBDAlBoxX=24.00*inch;   //in x
		double pBDAlBoxZ=15.25*inch;   //in z
		double pBDAlBoxY=16.812*inch;   //in y 
		double pBDAlBoxThick=0.0625*inch;    //for side and bottom
		G4VSolid* BDAlBoxSideSolid = new G4Box("BDAlBoxSide",pBDAlBoxThick/2.0,
			pBDAlBoxY/2.0,pBDAlBoxZ/2.0);

		G4LogicalVolume* BDAlBoxSideLogical = new G4LogicalVolume(BDAlBoxSideSolid,
			mMaterialManager->aluminum,"BDAlBoxSideLogical",0,0,0);
		BDAlBoxSideLogical->SetVisAttributes(GrayVisAtt);  

		double pBDAlBoxLPos_X= pBDAlBoxX/2 + mPivotXOffset;
		double pBDAlBoxRPos_X=-pBDAlBoxX/2 + mPivotXOffset;
		//	double pBDAlBoxPos_Y=-11.88*inch + pBDAlBoxY/2 + mPivotYOffset;
		double pBDAlBoxPos_Y=pBDContainerPos_Y+4.358*inch-pBDAlBoxY/2.0;
		double pBDAlBoxPos_Z=pBDContainerPos_Z;
		new G4PVPlacement(0,
			G4ThreeVector(pBDAlBoxLPos_X,pBDAlBoxPos_Y,pBDAlBoxPos_Z),
			BDAlBoxSideLogical,"BDAlBoxLeftSidePhys",motherLogical,0,0);
		new G4PVPlacement(0,
			G4ThreeVector(pBDAlBoxRPos_X,pBDAlBoxPos_Y,pBDAlBoxPos_Z),
			BDAlBoxSideLogical,"BDAlBoxRightSidePhys",motherLogical,0,0);

		/////////////////////////////////////
		//By Jixie @ 20120125, place top, bottom and back panels

		//top plate
		double pBDAlBoxTopThick=0.625*inch;  //for top pannel
		G4VSolid* BDAlBoxTopSolid = new G4Box("BDAlBoxTop",
			pBDAlBoxX/2.0,pBDAlBoxTopThick/2.0,pBDAlBoxZ/2.0);

		G4LogicalVolume* BDAlBoxTopLogical = new G4LogicalVolume(BDAlBoxTopSolid,
			mMaterialManager->aluminum,"BDAlBoxTopLogical",0,0,0);
		BDAlBoxTopLogical->SetVisAttributes(GrayVisAtt);  

		//the top plate is on top of the local dump, which is 4.358" + 0.1 mm from the 
		//center of the dump, not from the beam line
		double pBDAlBoxTopPos_X= pBDContainerPos_X;
		double pBDAlBoxTopPos_Y= pBDContainerPos_Y+4.358*inch+0.1*mm+pBDAlBoxTopThick/2.0;
		double pBDAlBoxTopPos_Z= pBDContainerPos_Z;
		new G4PVPlacement(0,
			G4ThreeVector(pBDAlBoxTopPos_X,pBDAlBoxTopPos_Y,pBDAlBoxTopPos_Z),
			BDAlBoxTopLogical,"BDAlBoxTopPhys",motherLogical,0,0);

		/////////////////////////////////////

		//bottom plate, need to dig rectangle for the top of the stand
		G4VSolid* BDAlBoxBottomWholeSolid = new G4Box("BDAlBoxBottomWhole",
			pBDAlBoxX/2.0-pBDAlBoxThick,pBDAlBoxThick/2.0,pBDAlBoxZ/2.0);
		G4VSolid* BDAlBoxBottomSubRecSolid = new G4Box("BDAlBoxBottomSubRec",
			8.13*inch/2.0,pBDAlBoxThick/2.0+0.1*mm,4.04*inch/2.0);
		G4SubtractionSolid* BDAlBoxBottomSolid = new G4SubtractionSolid("BDAlBoxBottom",
			BDAlBoxBottomWholeSolid,BDAlBoxBottomSubRecSolid);

		G4LogicalVolume* BDAlBoxBottomLogical = new G4LogicalVolume(BDAlBoxBottomSolid,
			mMaterialManager->aluminum,"BDAlBoxBottomLogical",0,0,0);

		BDAlBoxBottomLogical->SetVisAttributes(GrayVisAtt);  

		//the top plate is on top of the local dump, which is 4.358" + 0.1 mm from the 
		//center of the dump, not from the beam line
		double pBDAlBoxBottomPos_X= pBDContainerPos_X;
		double pBDAlBoxBottomPos_Y= pBDAlBoxTopPos_Y-pBDAlBoxY-pBDAlBoxThick/2.0;
		double pBDAlBoxBottomPos_Z= pBDContainerPos_Z;
		new G4PVPlacement(0,
			G4ThreeVector(pBDAlBoxBottomPos_X,pBDAlBoxBottomPos_Y,pBDAlBoxBottomPos_Z),
			BDAlBoxBottomLogical,"BDAlBoxBottomPhys",motherLogical,0,0);

		/////////////////////////////////////

		//back plate, need to dig rectangle and circles for scattered electrons
		double pBDAlBoxBackThick=0.562*inch;  //for back pannel
		G4VSolid* BDAlBoxBackWholeSolid = new G4Box("BDAlBoxBackWhole",
			pBDAlBoxX/2.0-pBDAlBoxThick,pBDAlBoxY/2.0,pBDAlBoxBackThick/2.0);
		G4VSolid* BDAlBoxBackSubRec1Solid = new G4Box("BDAlBoxBackSubRec1",
			2.776*inch/2.0,2.891*inch,pBDAlBoxBackThick/2.0+0.1*mm);
		G4VSolid* BDAlBoxBackSubRec2Solid = new G4Box("BDAlBoxBackSubRec2",
			3.912*inch/2.0,3.693*inch,pBDAlBoxBackThick/2.0+0.1*mm);
		G4VSolid* BDAlBoxBackSubCircleSolid = new G4Tubs("BDAlBoxBackSubCircle",
			0,1.25*inch,pBDAlBoxBackThick/2.0+0.1*mm,0,360*deg);

		double pBDAlBoxBackSubRec1Pos_X = (2.403 + 2.776/2)*inch;
		double pBDAlBoxBackSubRec1Pos_Y = pBDAlBoxY/2-4.929*inch;
		double pBDAlBoxBackSubRec2Pos_X = (6.116 + 3.912/2)*inch;
		double pBDAlBoxBackSubRec2Pos_Y = pBDAlBoxY/2-4.929*inch;
		G4SubtractionSolid* BDAlBoxBackWholeSub1LSolid = new G4SubtractionSolid(
			"BDAlBoxBackWholeSub1LSolid",BDAlBoxBackWholeSolid,BDAlBoxBackSubRec1Solid,0,
			G4ThreeVector(pBDAlBoxBackSubRec1Pos_X,pBDAlBoxBackSubRec1Pos_Y,0));
		G4SubtractionSolid* BDAlBoxBackWholeSub1LRSolid = new G4SubtractionSolid(
			"BDAlBoxBackWholeSub1LRSolid",BDAlBoxBackWholeSub1LSolid,BDAlBoxBackSubRec1Solid,0,
			G4ThreeVector(-pBDAlBoxBackSubRec1Pos_X,pBDAlBoxBackSubRec1Pos_Y,0));
		G4SubtractionSolid* BDAlBoxBackWholeSub1LR2LSolid = new G4SubtractionSolid(
			"BDAlBoxBackWholeSub1LR2LSolid",BDAlBoxBackWholeSub1LRSolid,BDAlBoxBackSubRec2Solid,0,
			G4ThreeVector(pBDAlBoxBackSubRec2Pos_X,pBDAlBoxBackSubRec2Pos_Y,0));
		G4SubtractionSolid* BDAlBoxBackWholeSub1LR2LRSolid = new G4SubtractionSolid(
			"BDAlBoxBackWholeSub1LR2LRSolid",BDAlBoxBackWholeSub1LR2LSolid,BDAlBoxBackSubRec2Solid,0,
			G4ThreeVector(-pBDAlBoxBackSubRec2Pos_X,pBDAlBoxBackSubRec2Pos_Y,0));
		G4SubtractionSolid* BDAlBoxBackSolid = new G4SubtractionSolid(
			"BDAlBoxBackWholeSub1Solid",BDAlBoxBackWholeSub1LR2LRSolid,BDAlBoxBackSubCircleSolid,0,
			G4ThreeVector(0,pBDAlBoxBackSubRec2Pos_Y,0));

		G4LogicalVolume* BDAlBoxBackLogical = new G4LogicalVolume(BDAlBoxBackSolid,
			mMaterialManager->aluminum,"BDAlBoxBackLogical",0,0,0);

		BDAlBoxBackLogical->SetVisAttributes(GrayVisAtt);  

		double pBDAlBoxBackPos_X= pBDContainerPos_X;
		double pBDAlBoxBackPos_Y= pBDAlBoxPos_Y;
		double pBDAlBoxBackPos_Z= pBDContainerPos_Z+pBDAlBoxZ/2.0;
		new G4PVPlacement(0,
			G4ThreeVector(pBDAlBoxBackPos_X,pBDAlBoxBackPos_Y,pBDAlBoxBackPos_Z),
			BDAlBoxBackLogical,"BDAlBoxBackPhys",motherLogical,0,0);

	}

	/////////////////////////
	// beam dump lead block
	/////////////////////////
	//need to subtract 2 trapzoids and 4 rectangles, also both the front side corners for 
	//12.5 degree scattering electrons

	//The local beam dump is not symmetric in height(y),  in order to make things easy I build a
	//symmetric box then chop a box from the top to make it match the real dimention
	//Then I subtract 4 rectangles in order to to dig the opening for the center tungsten core
	//To understand these 4 rectangles you need to check the position of them, from the top to 
	//the bottom they are indexed with 1,2,3 and 4
	//Confirmed by Alan that the beam tunnel opening is 1.18(width) X 1.5(height) inchch, with 
	//the bottom 0.735*inch below the x-z plane 

	//a symmetric bolck (with extra head)
	G4VSolid* BDWithExtraBlockSolid = new G4Box("BDWithExtraBlockBox",mBeamDumpWidth/2.0,
		(mBeamDumpHeight+pBDExtraHeight)/2.0,mBeamDumpThick/2.0);

	//the extra box to be subtracted
	G4VSolid* BDExtraBlockSolid = new G4Box("BDExtraBlockBox",mBeamDumpWidth/2.0+0.1*mm,
		pBDExtraHeight/2.0,mBeamDumpThick/2.0+0.1*mm);

	//the TRUE beam dump whole block without openings
	G4SubtractionSolid* BDBlockSolid=new G4SubtractionSolid("BDBlockBox",
		BDWithExtraBlockSolid,BDExtraBlockSolid,0,
		G4ThreeVector(0,4.358*inch+pBDExtraHeight/2,0));

	//define 2 trapzoids to subtract for openings
	//Opening for 5.69 degrees setup:  
	//  The entrance:  46<X<87  -43<Y<50
	//  The exit:      58<X<106 -53<Y<58
	double pDz=mBeamDumpThick/2.0+0.1*mm;
	double pDy1=(50.0-(-43.0))*mm/2.0, pDx1=(87.0-46.0)*mm/2.0, pDx2=pDx1;
	double pDy2=(58.0-(-53.0))*mm/2.0, pDx3=(106.0-58.0)*mm/2.0, pDx4=pDx3;
	double pAlp1=0.0*rad,pAlp2=0.0*rad;
	double pTheta=atan((58+106-46-87)*mm/2.0/mBeamDumpThick)*rad,pPhi=0.0*rad;

	G4VSolid* BDOpeningLSolid = new G4Trap("BDOpeningLTrap",pDz,pTheta,pPhi,
		pDy1,pDx1,pDx2,pAlp1,pDy2,pDx3,pDx4,pAlp2);
	G4VSolid* BDOpeningRSolid = new G4Trap("BDOpeningRTrap",pDz,-pTheta,pPhi,
		pDy1,pDx1,pDx2,pAlp1,pDy2,pDx3,pDx4,pAlp2);

	double pBDTrapPos_X=((58+106)/2+(46+87)/2)/2.0*mm;
	G4SubtractionSolid* BDSubTrapLSolid=new G4SubtractionSolid("BDSubTrapLSolid",
		BDBlockSolid,BDOpeningLSolid,0,G4ThreeVector(pBDTrapPos_X,0,0));
	G4SubtractionSolid* BDSubTrapLRSolid=new G4SubtractionSolid("BDSubTrapLRSolid",
		BDSubTrapLSolid,BDOpeningRSolid,0,G4ThreeVector(-pBDTrapPos_X,0,0));


	//Need to dig 4 rectangles and then place the Tu core, rad con plates
	//the 1st box is the beam tunnel(1.18 X 1.5), whose bottom is 0.735 inch to the x-z plane
	//or its center is 1.5/2-0.735 = 0.015 inch above the x-z plane

	double pBDBeamTunnelLeftX=1.18*inch/2;
	double pBDBeamTunnelBottomY=-0.735*inch;

	G4VSolid* BDRec1Solid = new G4Box("BDRec1Box",pBDBeamTunnelLeftX,
		1.5*inch/2.0,mBeamDumpThick/2.0+0.1*mm);
	G4ThreeVector pRec1Center(0,0.75*inch+pBDBeamTunnelBottomY,0);

	G4VSolid* BDRec2Solid = new G4Box("BDRec2Box",3.675*inch/2.0,
		2.516*inch/2.0,mBeamDumpThick/2.0+0.1*mm);
	G4ThreeVector pRec2Center(0,-2.516/2*inch,0);

	G4VSolid* BDRec3Solid = new G4Box("BDRec3Box",4.065*inch/2.0,
		2.016*inch/2.0,mBeamDumpThick/2.0+0.1*mm);
	G4ThreeVector pRec3Center(0,-(2.516+2.016/2)*inch,0);

	G4VSolid* BDRec4Solid = new G4Box("BDRec4Box",1.62*inch/2.0,
		1.836*inch/2.0,mBeamDumpThick/2.0+0.1*mm);
	G4ThreeVector pRec4Center(0,-(4.532+1.836/2.0)*inch,0);

	G4SubtractionSolid* BDSubTrapLRNRec1Solid=new G4SubtractionSolid("BDSubTrapLRNRec1Solid",
		BDSubTrapLRSolid,BDRec1Solid,0,pRec1Center);

	G4SubtractionSolid* BDSubTrapLRNRec12Solid=new G4SubtractionSolid("BDSubTrapLRNRec12Solid",
		BDSubTrapLRNRec1Solid,BDRec2Solid,0,pRec2Center);

	G4SubtractionSolid* BDSubTrapLRNRec123Solid=new G4SubtractionSolid("BDSubTrapLRNRec123Solid",
		BDSubTrapLRNRec12Solid,BDRec3Solid,0,pRec3Center);

	G4SubtractionSolid* BDSubTrapLRNRec1234Solid=new G4SubtractionSolid("BDSubTrapLRNRec1234Solid",
		BDSubTrapLRNRec123Solid,BDRec4Solid,0,pRec4Center);

	/////////////////////////////////////////////
	//Now need to subtract a right-angled triangle the front side corner, 
	//the hypotenuse has a angle 10.6 degree w.r.t z axis and its height is 0.214 inch 
	//
	//Need to define a box and calculate its relative position to the local dump block carefully
	//then subtract this box with the right orientation at the right position
	double pSubTriangleHeigh=0.214*inch, pSubTriangleAngle=10.6*deg;
	double pRecCornerX=pSubTriangleHeigh*2, pRecCornerY=5.375*inch, pRecCornerZ=mBeamDumpThick*2;	
	G4VSolid* BDRecCornerSolid = new G4Box("BDRecCornerBox",pRecCornerX/2.0,
		pRecCornerY/2.0,pRecCornerZ/2.0);
	double pRecCornerPos_X=mBeamDumpWidth/2 - pSubTriangleHeigh/cos(pSubTriangleAngle) + 
		mBeamDumpThick/2 * tan(pSubTriangleAngle) + pRecCornerX/2/cos(pSubTriangleAngle);	

	G4RotationMatrix *pRotRecCornerL=new G4RotationMatrix();
	pRotRecCornerL->rotateY(-pSubTriangleAngle); 
	G4RotationMatrix *pRotRecCornerR=new G4RotationMatrix();
	pRotRecCornerR->rotateY( pSubTriangleAngle); 
	G4ThreeVector pRecCornerLCenter( pRecCornerPos_X,0,0);		
	G4ThreeVector pRecCornerRCenter(-pRecCornerPos_X,0,0);

	G4SubtractionSolid* BDSubTrapLRNRec1234NCornerLSolid=new G4SubtractionSolid(
		"BDSubTrapLRNRec1234NCornerLSolid", BDSubTrapLRNRec1234Solid,BDRecCornerSolid,
		pRotRecCornerL,pRecCornerLCenter);

	G4SubtractionSolid* BDSubTrapLRNRec1234NCornerLRSolid=new G4SubtractionSolid(
		"BDSubTrapLRNRec1234NCornerLRSolid", BDSubTrapLRNRec1234NCornerLSolid,BDRecCornerSolid,
		pRotRecCornerR,pRecCornerRCenter);

	////place the subtracted corner rectangle for debug	
	////By Jixie: the rotation and geometry have been debuged and they are all correct
	//G4LogicalVolume* BDRecCornerLogical = new G4LogicalVolume(BDRecCornerSolid,
	//	vacuum,"beamDumpLogical",0,0,0);
	//BDRecCornerLogical->SetVisAttributes(SkyBlueVisAtt);  
	//
	//new G4PVPlacement(pRotRecCornerL,pRecCornerLCenter,
	//	BDRecCornerLogical,"BDRecCornerPhys",beamDumpContainerLogical,true,0);
	//new G4PVPlacement(pRotRecCornerR,pRecCornerRCenter,
	//	BDRecCornerLogical,"BDRecCornerPhys",beamDumpContainerLogical,true,1);
	/////////////////////////////////////////////

	//place the lead block
	G4LogicalVolume* beamDumpLogical = new G4LogicalVolume(BDSubTrapLRNRec1234NCornerLRSolid,
		mMaterialManager->lead,"beamDumpLogical",0,0,0);
	beamDumpLogical->SetVisAttributes(LeadVisAtt);  

	new G4PVPlacement(0,G4ThreeVector(0,0,0),
		beamDumpLogical,"beamDumpPhys",beamDumpContainerLogical,0,0);


	/////////////////////////////////////
	//mMaterialManager->tungsten, cores + sides
	/////////////////////////////////////

	double pCenterBlockX=2.665*inch;
	double pCenterBlockY=2.642*inch;
	double pCenterBlockZ=10.0*cm;

	G4VSolid* BDRec5Solid = new G4Box("BDRec4Box",pCenterBlockX/2.0,
		pCenterBlockY/2.0,mBeamDumpThick/2.0+0.2*mm);

	//up side, rec2-rec5
	G4SubtractionSolid* BDTuRec2SubRec5Solid=new G4SubtractionSolid("BDTuRec2SubRec5Solid",
		BDRec2Solid,BDRec5Solid,0,G4ThreeVector(0,0,0));
	G4LogicalVolume* BDRec2SubRec5Logical = new G4LogicalVolume(BDTuRec2SubRec5Solid,
		mMaterialManager->tungsten,"BDRec2SubRec5Logical",0,0,0);
	new G4PVPlacement(0,pRec2Center,
		BDRec2SubRec5Logical,"BDRec2SubRec5Phys",beamDumpContainerLogical,0,0);
	BDRec2SubRec5Logical->SetVisAttributes(LightGreenVisAtt);

	//bottom side, rec3-rec5
	G4SubtractionSolid* BDTuRec3SubRec5Solid=new G4SubtractionSolid("BDTuRec3SubRec5Solid",
		BDRec3Solid,BDRec5Solid,0,G4ThreeVector(0,(0.8+(2.642-2.016)/2)*inch,0));
	G4LogicalVolume* BDRec3SubRec5Logical = new G4LogicalVolume(BDTuRec3SubRec5Solid,
		mMaterialManager->tungsten,"BDRec3SubRec5Logical",0,0,0);
	new G4PVPlacement(0,pRec3Center,
		BDRec3SubRec5Logical,"BDRec3SubRec5Phys",beamDumpContainerLogical,0,0);
	BDRec3SubRec5Logical->SetVisAttributes(LightGreenVisAtt);  


	//central block, placed 1.09 inch away from the x-z plane
	double pCenterBlockPos_Y=-(1.09*inch+pCenterBlockY/2.0);
	double pCenterBlockPos_Z=mBeamDumpThick/2-pCenterBlockZ/2;
	G4VSolid* BDTuCenterBlockSolid = new G4Box("BDTuCenterBlockBox",pCenterBlockX/2.0,
		pCenterBlockY/2.0,pCenterBlockZ/2.0);

	G4LogicalVolume* BDTuCenterBlockLogical = new G4LogicalVolume(BDTuCenterBlockSolid,
		mMaterialManager->tungsten,"BDTuCenterBlockLogical",0,0,0);
	new G4PVPlacement(0,G4ThreeVector(0,pCenterBlockPos_Y,pCenterBlockPos_Z),
		BDTuCenterBlockLogical,"BDTuCenterBlockPhys",beamDumpContainerLogical,0,0);
	BDTuCenterBlockLogical->SetVisAttributes(LightGreenVisAtt);  

	//place the bottom block, which is for 1.2 GeV beam @ HRS 12.5 degrees 	
	//the real size is 1.62 x 1.961 x 5.9 inches.  I use 1.836 inch in height here
	//since the protruded part is not dig out in the above central block 
	double pBDTuBottomBlockZ=mBeamDumpThick-6.0*mm;  //6mm behind the front face
	double pBDTuBottomBlockPos_Y=-(4.532+1.836/2.0)*inch;
	double pBDTuBottomBlockPos_Z=mBeamDumpThick/2-pBDTuBottomBlockZ/2;
	G4VSolid* BDTuBottomBlockSolid = new G4Box("BDTuBottomBlockBox",1.62*inch/2.0,
		1.836*inch/2.0,pBDTuBottomBlockZ/2.0);

	G4LogicalVolume* BDTuBottomBlockLogical = new G4LogicalVolume(BDTuBottomBlockSolid,
		mMaterialManager->tungsten,"BDTuBottomBlockLogical",0,0,0);
	new G4PVPlacement(0,G4ThreeVector(0,pBDTuBottomBlockPos_Y,pBDTuBottomBlockPos_Z),
		BDTuBottomBlockLogical,"BDTuBottomBlockPhys",beamDumpContainerLogical,0,0);
	BDTuBottomBlockLogical->SetVisAttributes(LightGreenVisAtt);  


	/////////////////////////////////////
	//The mMaterialManager->tungsten, plates in the cartridge 
	/////////////////////////////////////

	//Above the central block (2.665 inch wide), there is another
	//Tusgten plate of 0.47 inch high (see drawing A080270201-0303). 
	//If I built it as 0.47 inch, the top of it to the x-z plane is not 0.735 inch, 
	//I have to make it thinner here:   1.09 - 0.735 = 0.355 inch
	double pBDTuTopBlockX=pCenterBlockX;
	double pBDTuTopBlockY=1.09*inch+pBDBeamTunnelBottomY;
	double pBDTuTopBlockZ=4.0*inch; //pCenterBlockZ; 
	G4VSolid* BDTuTopBlockSolid = new G4Box("BDTuTopBlockBox",pBDTuTopBlockX/2.0,
		pBDTuTopBlockY/2.0,pBDTuTopBlockZ/2.0);
	G4LogicalVolume* BDTuTopBlockLogical = new G4LogicalVolume(BDTuTopBlockSolid,
		mMaterialManager->tungsten,"BDTuTopBlockLogical",0,0,0);
	BDTuTopBlockLogical->SetVisAttributes(LightGreenVisAtt);  

	double pBDTuTopBlockPos_Y=pBDBeamTunnelBottomY-pBDTuTopBlockY/2;
	double pBDTuTopBlockPos_Z=mBeamDumpThick/2-pBDTuTopBlockZ/2;
	new G4PVPlacement(0,G4ThreeVector(0,pBDTuTopBlockPos_Y,pBDTuTopBlockPos_Z),
		BDTuTopBlockLogical,"BDTuTopBlockPhys",beamDumpContainerLogical,0,0);

	//above the TopBlock there are 2 Tu side plates, which build the beam tunnel
	//The top of these 2 blocks is at x-z plane
	double pBDTuTopSideBlockX=0.722*inch;
	double pBDTuTopSideBlockY=fabs(pBDBeamTunnelBottomY);
	double pBDTuTopSideBlockZ=pBDTuTopBlockZ; 
	G4VSolid* BDTuTopSideBlockSolid = new G4Box("BDTuTopSideBlockBox",pBDTuTopSideBlockX/2.0,
		pBDTuTopSideBlockY/2.0,pBDTuTopSideBlockZ/2.0);
	G4LogicalVolume* BDTuTopSideBlockLogical = new G4LogicalVolume(BDTuTopSideBlockSolid,
		mMaterialManager->tungsten,"BDTuTopSideBlockLogical",0,0,0);
	BDTuTopSideBlockLogical->SetVisAttributes(LightGreenVisAtt);  

	double pBDTuTopSideBlockPos_X=pBDBeamTunnelLeftX+pBDTuTopSideBlockX/2;
	double pBDTuTopSideBlockPos_Y=-pBDTuTopSideBlockY/2;
	double pBDTuTopSideBlockPos_Z=mBeamDumpThick/2-pBDTuTopSideBlockZ/2;
	new G4PVPlacement(0,
		G4ThreeVector(pBDTuTopSideBlockPos_X,pBDTuTopSideBlockPos_Y,pBDTuTopSideBlockPos_Z),
		BDTuTopSideBlockLogical,"BDTuTopSideBlockPhys",beamDumpContainerLogical,true,0);
	new G4PVPlacement(0,
		G4ThreeVector(-pBDTuTopSideBlockPos_X,pBDTuTopSideBlockPos_Y,pBDTuTopSideBlockPos_Z),
		BDTuTopSideBlockLogical,"BDTuTopSideBlockPhys",beamDumpContainerLogical,true,1);

	//The lead block of the cartrige have been subtracted more than it need, not I have to
	//add 2 side block backbecause the above 2 side blocks is shorter than the local dump thickness
	double pBDLeadTopSideBlockX=0.722*inch;
	double pBDLeadTopSideBlockY=fabs(pBDBeamTunnelBottomY);
	double pBDLeadTopSideBlockZ=mBeamDumpThick-pBDTuTopSideBlockZ; 
	G4VSolid* BDLeadTopSideBlockSolid = new G4Box("BDLeadTopSideBlockBox",pBDLeadTopSideBlockX/2.0,
		pBDLeadTopSideBlockY/2.0,pBDLeadTopSideBlockZ/2.0);
	G4LogicalVolume* BDLeadTopSideBlockLogical = new G4LogicalVolume(BDLeadTopSideBlockSolid,
		mMaterialManager->lead,"BDLeadTopSideBlockLogical",0,0,0);
	BDLeadTopSideBlockLogical->SetVisAttributes(LeadVisAtt);  

	double pBDLeadTopSideBlockPos_X=pBDBeamTunnelLeftX+pBDLeadTopSideBlockX/2;
	double pBDLeadTopSideBlockPos_Y=-pBDLeadTopSideBlockY/2;
	double pBDLeadTopSideBlockPos_Z=mBeamDumpThick/2-pBDTuTopSideBlockZ-pBDLeadTopSideBlockZ/2;
	new G4PVPlacement(0,
		G4ThreeVector(pBDLeadTopSideBlockPos_X,pBDLeadTopSideBlockPos_Y,pBDLeadTopSideBlockPos_Z),
		BDLeadTopSideBlockLogical,"BDLeadTopSideBlockPhys",beamDumpContainerLogical,true,0);
	new G4PVPlacement(0,
		G4ThreeVector(-pBDLeadTopSideBlockPos_X,pBDLeadTopSideBlockPos_Y,pBDLeadTopSideBlockPos_Z),
		BDLeadTopSideBlockLogical,"BDLeadTopSideBlockPhys",beamDumpContainerLogical,true,1);

	/////////////////////////////////////
	//removalbe mMaterialManager->tungsten, plates 
	/////////////////////////////////////

	//TODO: For GEP, there will be an extra Tungsten block which is placed in the beam hole
	//This block is still not available yet, need to add it here
	//some times this block is not use, I need to check with Al in what situation  it will be used
	bool pSetupRemovalbleTuPlate=false;
	if(pSetupRemovalbleTuPlate)
	{
		double pRemovableTuPlateX=1.165*inch;
		double pRemovableTuPlateY=0.156*inch;
		double pRemovableTuPlateZ=4.000*inch;
		G4VSolid* removableTuPlateSolid = new G4Box("removableTuPlateBox",pRemovableTuPlateX/2.0,
			pRemovableTuPlateY/2.0,pRemovableTuPlateZ/2.0);
		G4LogicalVolume* removableTuPlateLogical = new G4LogicalVolume(removableTuPlateSolid,
			mMaterialManager->tungsten,"removableTuPlateLogical",0,0,0);
		removableTuPlateLogical->SetVisAttributes(LightGreenVisAtt);  

		double pRemovableTuPlatePos_Y=pBDBeamTunnelBottomY+pRemovableTuPlateY/2;
		double pRemovableTuPlatePos_Z=mBeamDumpThick/2-pRemovableTuPlateZ/2;
		new G4PVPlacement(0, G4ThreeVector(0,pRemovableTuPlatePos_Y,pRemovableTuPlatePos_Z),
			removableTuPlateLogical,"removableTuPlatePhys",beamDumpContainerLogical,0,0);
	}

	/////////////////////////////////////
	//rad con test plates 
	/////////////////////////////////////
	//14mm Cu + 0.492inch(12.5mm) Tu + 14mm Al
	//place the rad con Tu plates in front of the center mMaterialManager->tungsten, block
	//the real zide is 2.5(width) X 3.061 (height) X thickness
	//Here I use 2.625 to replace 2.5 since I did not build the aluminum holder box
	//Tu is the middle plate, place it first
	//double pBDRadConAlBoxThick=0.0625*inch;
	//double pRadConTuPlateX=2.5*inch+2*pBDRadConAlBoxThick;
	//double pRadConTuPlateY=3.061*inch+pBDRadConAlBoxThick;
	//By Jixie: I want to make the rad con play just 0.735" away from the x-z plane
	//ignore the dimention given above

	double pRadConTuPlateX=pCenterBlockX;
	double pRadConTuPlateY=pBDTuTopBlockY+pCenterBlockY;
	double pRadConTuPlateZ=12.5*mm;
	G4VSolid* radConTuPlateSolid = new G4Box("radConTuPlateBox",pRadConTuPlateX/2.0,
		pRadConTuPlateY/2.0,pRadConTuPlateZ/2.0);
	G4LogicalVolume* radConTuPlateLogical = new G4LogicalVolume(radConTuPlateSolid,
		mMaterialManager->tungsten,"radConTuPlateLogical",0,0,0);

	//the rad con plate box is 1.91*inch thick, the Tu plate is located at the center of it
	//the beam tunnel opening is 1.18(width) X 1.5 (height) inch. The top of these rad con plates
	//is the bottom of this tunnel
	double pRadConPlatePos_Y=pBDBeamTunnelBottomY-pRadConTuPlateY/2;
	double pRadConTuPlatePos_Z=mBeamDumpThick/2-pCenterBlockZ-1.91*inch/2;
	new G4PVPlacement(0, G4ThreeVector(0,pRadConPlatePos_Y,pRadConTuPlatePos_Z),
		radConTuPlateLogical,"radConTuPlatePhys",beamDumpContainerLogical,0,0);
	radConTuPlateLogical->SetVisAttributes(LightGreenVisAtt);  


	//place the rad con Cu plates in front of the rad con mMaterialManager->tungsten, plate
	double pRadConCuPlateX=pRadConTuPlateX;
	double pRadConCuPlateY=pRadConTuPlateY;
	double pRadConCuPlateZ=1.4*cm;
	G4VSolid* radConCuPlateSolid = new G4Box("radConCuPlateBox",pRadConCuPlateX/2.0,
		pRadConCuPlateY/2.0,pRadConCuPlateZ/2.0);
	G4LogicalVolume* radConCuPlateLogical = new G4LogicalVolume(radConCuPlateSolid,
		mMaterialManager->copper,"radConCuPlateLogical",0,0,0);

	double pRadConCuPlatePos_Z=pRadConTuPlatePos_Z-pRadConTuPlateZ/2-pRadConCuPlateZ/2;
	new G4PVPlacement(0, G4ThreeVector(0,pRadConPlatePos_Y,pRadConCuPlatePos_Z),
		radConCuPlateLogical,"radConCuPlatePhys",beamDumpContainerLogical,0,0);
	radConCuPlateLogical->SetVisAttributes(CuBrownVisAtt);  

	//place the rad con Al plates behind the rad con mMaterialManager->tungsten, plate
	double pRadConAlPlateX=pRadConTuPlateX;
	double pRadConAlPlateY=pRadConTuPlateY;
	double pRadConAlPlateZ=1.4*cm;
	G4VSolid* radConAlPlateSolid = new G4Box("radConAlPlateBox",pRadConAlPlateX/2.0,
		pRadConAlPlateY/2.0,pRadConAlPlateZ/2.0);
	G4LogicalVolume* radConAlPlateLogical = new G4LogicalVolume(radConAlPlateSolid,
		mMaterialManager->aluminum,"radConAlPlateLogical",0,0,0);

	double pRadConAlPlatePos_Z=pRadConTuPlatePos_Z+pRadConTuPlateZ/2+pRadConAlPlateZ/2;
	new G4PVPlacement(0, G4ThreeVector(0,pRadConPlatePos_Y,pRadConAlPlatePos_Z),
		radConAlPlateLogical,"radConAlPlatePhys",beamDumpContainerLogical,0,0);
	radConAlPlateLogical->SetVisAttributes(GrayVisAtt);  



	/////////////////////////////////////
	//local dump stand 
	/////////////////////////////////////
	//This stand is make of aluminum. It have and -| shape foot (made of rectangle pipe)
	//and a rectangle pine body and the head (top) is a plate
	//if mSetupBeamDump<2, this stand will not be setup, which can make the program faster

	if(mSetupBeamDump>=2)
	{
		G4RotationMatrix *pRotX90deg=new G4RotationMatrix();
		pRotX90deg->rotateX(-90*deg); 
		G4RotationMatrix *pRotY90deg=new G4RotationMatrix();
		pRotY90deg->rotateY(-90*deg); 

		G4RotationMatrix *pRotX45deg=new G4RotationMatrix();
		pRotX45deg->rotateX(-45*deg); 
		G4RotationMatrix *pRotX135deg=new G4RotationMatrix();
		pRotX135deg->rotateX(-135*deg); 

		//To bulid the rec pipe, I have to subtract a samll box from the large box
		//This is the small box to be subtracted,
		G4VSolid* BDStandHoleSolid = new G4Box("BDStandHoleBox",3.5*inch/2.0,
			3.5*inch/2.0,20.0*inch/2.0);
		//In roder to chop the end with 45 degree, I have to define this box
		G4VSolid* BDStandEndChopSolid = new G4Box("BDStandEndChopBox",4.01*inch/2.0,
			4.01*inch/2.0,20.0*inch/2.0);

		/////////////////////////////////////
		//local dump stand Head plate
		/////////////////////////////////////
		double pBDStandHeadX=11.0*inch;
		double pBDStandHeadY=0.5*inch;
		double pBDStandHeadZ=7.5*inch;
		G4VSolid* BDStandHeadPlateSolid = new G4Box("BDStandHeadPlateBox",pBDStandHeadX/2.0,
			pBDStandHeadY/2.0,pBDStandHeadZ/2.0);
		//The head plate in the real design have dig a thin box of 9.312(width) x 0.125(y) x 5.966 (z) inches.
		G4VSolid* BDStandHeadSubSolid = new G4Box("BDStandHeadSubBox",9.312*inch/2.0,
			0.125*inch/2.0,5.966*inch/2.0);

		G4SubtractionSolid* BDStandHeadSolid=new G4SubtractionSolid("BDStandHeadSolid",
			BDStandHeadPlateSolid,BDStandHeadSubSolid,0,
			G4ThreeVector(0,pBDStandHeadY/2-0.125*inch/2,0));

		G4LogicalVolume* BDStandHeadLogical = new G4LogicalVolume(BDStandHeadSolid,
			mMaterialManager->aluminum,"BDStandHeadLogical",0,0,0);
		BDStandHeadLogical->SetVisAttributes(GrayVisAtt); 

		//subtrct 0.2 mm for the thickness of the BD containner
		double pBDStandHeadPos_X=pBDContainerPos_X;
		double pBDStandHeadPos_Y=-(mBeamDumpHeight+pBDExtraHeight)/2.0-pBDStandHeadY/2.0+0.125*inch-0.2*mm;
		double pBDStandHeadPos_Z=pBDContainerPos_Z;
		new G4PVPlacement(0, 
			G4ThreeVector(pBDStandHeadPos_X,pBDStandHeadPos_Y,pBDStandHeadPos_Z),
			BDStandHeadLogical,"BDStandHeadPhys",motherLogical,0,0);

		/////////////////////////////////////
		//local dump stand body pipe
		/////////////////////////////////////
		double pBDStandBodyHeight=8.165*inch;
		G4VSolid* BDStandBodySolid = new G4Box("BDStandBodyBox",4.0*inch/2.0,
			4.0*inch/2.0,pBDStandBodyHeight/2.0);
		//subtract the hole from the body
		G4SubtractionSolid* BDStandBodyPipeSolid=new G4SubtractionSolid("BDStandBodyPipeSolid",
			BDStandBodySolid,BDStandHoleSolid);

		G4LogicalVolume* BDStandBodyPipeLogical = new G4LogicalVolume(BDStandBodyPipeSolid,
			mMaterialManager->aluminum,"BDStandBodyPipeLogical",0,0,0);
		BDStandBodyPipeLogical->SetVisAttributes(GrayVisAtt);  
		double pBDStandBodyPos_X=pBDContainerPos_X;
		double pBDStandBodyPos_Y=pBDStandHeadPos_Y-pBDStandHeadY/2-pBDStandBodyHeight/2.0;
		double pBDStandBodyPos_Z=pBDContainerPos_Z;
		new G4PVPlacement(pRotX90deg, 
			G4ThreeVector(pBDStandBodyPos_X,pBDStandBodyPos_Y,pBDStandBodyPos_Z),
			BDStandBodyPipeLogical,"BDStandBodyPipePhys",motherLogical,0,0);


		/////////////////////////////////////
		//local dump stand -| shape foot
		/////////////////////////////////////
		//foot X pipe
		G4VSolid* BDStandFootXBarSolid = new G4Box("BDStandFootXBarBox",4.0*inch/2.0,
			4.0*inch/2.0,19.0*inch/2.0);
		G4SubtractionSolid* BDStandFootXPipeSolid=new G4SubtractionSolid("BDStandFootXPipeSolid",
			BDStandFootXBarSolid,BDStandHoleSolid);
		//chop the end with 45 degree,the end position is z = b + (1+sqrt(2))/2 *a
		//where b is half length of none-subtract part and a is the square side
		double pBDStandFootXSubPos_Z=((19.0-4.0*2)/2 + (1.0+sqrt(2.0))/2*4.0 ) * inch;
		G4SubtractionSolid* BDStandFootXPipeSubEndUSolid=new G4SubtractionSolid(
			"BDStandFootXPipeSubEndUSolid",BDStandFootXPipeSolid,BDStandEndChopSolid,
			pRotX135deg,G4ThreeVector(0,0,-pBDStandFootXSubPos_Z));
		G4SubtractionSolid* BDStandFootXPipeSubEndUDSolid=new G4SubtractionSolid(
			"BDStandFootXPipeSubEndUDSolid",BDStandFootXPipeSubEndUSolid,BDStandEndChopSolid,
			pRotX45deg,G4ThreeVector(0,0,pBDStandFootXSubPos_Z));

		G4LogicalVolume* BDStandFootXLogical = new G4LogicalVolume(BDStandFootXPipeSubEndUDSolid,
			mMaterialManager->aluminum,"BDStandFootXLogical",0,0,0);
		BDStandFootXLogical->SetVisAttributes(GrayVisAtt);  

		double pBDStandFootPos_X=pBDContainerPos_X;
		double pBDStandFootPos_Y=pBDStandBodyPos_Y-pBDStandBodyHeight/2-4*inch/2.0;
		double pBDStandFootXPos_Z=pBDContainerPos_Z+3.906*inch;
		new G4PVPlacement(pRotY90deg, 
			G4ThreeVector(pBDStandFootPos_X,pBDStandFootPos_Y,pBDStandFootXPos_Z),
			BDStandFootXLogical,"BDStandFootXPhys",motherLogical,0,0);


		//foot Z pipe
		G4VSolid* BDStandFootZBarSolid = new G4Box("BDStandFootZBarBox",4.0*inch/2.0,
			4.0*inch/2.0,12.501*inch/2.0);
		G4SubtractionSolid* BDStandFootZPipeSolid=new G4SubtractionSolid("BDStandFootZPipeSolid",
			BDStandFootZBarSolid,BDStandHoleSolid);
		//chop upstream end by 45 degrees
		double pBDStandFootZSubPos_Z=((12.501-4.0*2)/2 + (1.0+sqrt(2.0))/2*4.0) * inch;
		G4SubtractionSolid* BDStandFootZPipeSubEndUSolid=new G4SubtractionSolid(
			"BDStandFootZPipeSubEndUSolid",BDStandFootZPipeSolid,BDStandEndChopSolid,
			pRotX135deg,G4ThreeVector(0,0,-pBDStandFootZSubPos_Z));

		G4LogicalVolume* BDStandFootZLogical = new G4LogicalVolume(BDStandFootZPipeSubEndUSolid,
			mMaterialManager->aluminum,"BDStandFootZLogical",0,0,0);
		BDStandFootZLogical->SetVisAttributes(GrayVisAtt);  

		double pBDStandFootZPos_Z=pBDContainerPos_Z-(12.501/2+(4.0-3.906)-4.0/2)*inch;
		new G4PVPlacement(0, 
			G4ThreeVector(pBDStandFootPos_X,pBDStandFootPos_Y,pBDStandFootZPos_Z),
			BDStandFootZLogical,"BDStandFootZPhys",motherLogical,0,0);
	}
	return beamDumpPhysical;
}


/////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* G2PDetectorConstruction::ConstructG2PThirdArm(G4LogicalVolume* motherLogical)
{
	const double inch=2.54*cm;
	G4String SDname;

	// All G4VSensitiveDetector will be managed (deleted) by SDManager
	//therefore no need to delete it at the end of this subroutine
	G4VSensitiveDetector* thirdArmSC1SD=new HRSStdSD(SDname="thirdArmSC1");
	G4VSensitiveDetector* thirdArmSC2SD=new HRSStdSD(SDname="thirdArmSC2");
	G4VSensitiveDetector* thirdArmShieldingSD=new HRSStdSD(SDname="thirdArmShielding");

	G4SDManager* SDman = G4SDManager::GetSDMpointer();

	//The thirdArm's front face is 2.1 m away from the target center, located at 
	//left side 70 degrees and tilted 0 degrees, but tilted with various angles
	//mThirdArmAngle=70.0*deg,mThirdArmRotXAngle=0.0*deg,mThirdArmRotZAngle=0.0*deg;
	//mPivot2ThirdArmFace=210.0*cm;

	double pThirdArmXoffset,pThirdArmYoffset,pThirdArmZoffset;

	double pBeamEnergy=2.257;
	gConfig->GetArgument("BeamEnergy",pBeamEnergy);

	double pHelm_CurrentRatio=1.0;
	gConfig->GetParameter("Helm_CurrentRatio",pHelm_CurrentRatio);

	//for 1.2 GeV, the 3rd arm will detected 270 MeV<Pvb<410
	//for others,  the 3rd arm will detected 300 MeV<Pvb<450
	//therefore the rotation angle will be different
	if(pBeamEnergy<1.3) mThirdArmRotZAngle=-36.0*deg*pHelm_CurrentRatio;
	else  mThirdArmRotZAngle=-32.0*deg*pHelm_CurrentRatio;

	double pTAHoffset=0.0,pTAVoffset=0.0,pTAShieldThick=0.0;
	gConfig->GetArgument("TAHoffset",pTAHoffset); 
	pTAHoffset*=mm;
	gConfig->GetArgument("TAVoffset",pTAVoffset); 
	pTAVoffset*=mm;
	gConfig->GetArgument("TAShieldThick",pTAShieldThick); 
	pTAShieldThick*=mm;


	if((int)pTAHoffset==-9999 && (int)pTAVoffset==-9999)
	{
		//use default position offset	
		if(pBeamEnergy<1.3)
		{
			pTAHoffset=60*mm;
			pTAVoffset=0*mm;
			pTAShieldThick=0*mm;
		}
		else if(pBeamEnergy<1.8)
		{
			pTAHoffset=0*mm;
			pTAVoffset=70*mm;
			pTAShieldThick=2*mm;
		}
		else if(pBeamEnergy<2.3)
		{
			pTAHoffset=0*mm;
			pTAVoffset=-40*mm;
			pTAShieldThick=2*mm;
		}
	}

	pThirdArmXoffset=pTAHoffset*cos(mThirdArmAngle);
	pThirdArmZoffset=pTAHoffset*sin(mThirdArmAngle);
	pThirdArmYoffset=pTAVoffset;


	//By Jixie and Chao:
	//The rotation matrix will rotate the whole coordinate system
	//For instance: G4RotationMatrix::rotateZ(alpha) then G4RotationMatrix::rotateX(beta)
	//will rotate about z by alpha then rotate about NEW X aixs by beta angle
	G4RotationMatrix *pRotThirdArm=new G4RotationMatrix();
	pRotThirdArm->rotateY(-mThirdArmAngle); 
	if(fabs(mThirdArmRotZAngle)>1.0E-05) pRotThirdArm->rotateZ(-mThirdArmRotZAngle);

	G4RotationMatrix *pRotYThirdArmAngle=new G4RotationMatrix();
	pRotYThirdArmAngle->rotateY(-mThirdArmAngle); 

	/////////////////////////////////////
	//shielding, in front of the detector 
	/////////////////////////////////////

	//A aluminum foil 70 cm away to the target center, covering +/- 20 degrees in theta_tr and
	//also covering +/- 30 degrees in phi_tr 
	bool mSetupThirdArmShield=true;
	if(mSetupThirdArmShield && pTAShieldThick>0.001*mm)
	{
		double pPivot2TAShieldFace=60*cm;
		double pTAShieldX=2*3.1416*pPivot2TAShieldFace*40./360;
		double pTAShieldY=2*3.1416*pPivot2TAShieldFace*60./360;
		double pTAShieldPos_X=(pPivot2TAShieldFace+pTAShieldThick/2.0)*sin(mThirdArmAngle)+mPivotXOffset;
		double pTAShieldPos_Y=mPivotYOffset;
		double pTAShieldPos_Z=(pPivot2TAShieldFace+pTAShieldThick/2.0)*cos(mThirdArmAngle)+mPivotZOffset;

		G4RotationMatrix *pRotTAShield=new G4RotationMatrix();
		pRotTAShield->rotateY(-mThirdArmAngle); 
		G4VSolid* TAShieldSolid = new G4Box("TAShieldBox",pTAShieldX/2.0, pTAShieldY/2.0, pTAShieldThick/2.0);
		G4LogicalVolume* TAShieldLogical = new G4LogicalVolume(TAShieldSolid,
			mMaterialManager->aluminum,"ThirdArmShieldLogical",0,0,0);
		//MuMetal,"ThirdArmShieldLogical",0,0,0);  //test the mu-metal blocking ability
		new G4PVPlacement(pRotTAShield,
			G4ThreeVector(pTAShieldPos_X,pTAShieldPos_Y,pTAShieldPos_Z),
			TAShieldLogical,"TAShieldPhys",motherLogical,0,0,0);
		TAShieldLogical->SetVisAttributes(GrayVisAtt); 
		SDman->AddNewDetector(thirdArmShieldingSD);
		TAShieldLogical->SetSensitiveDetector(thirdArmShieldingSD);
	}

	/////////////////////////
	// ThirdArm Stand
	/////////////////////////
	//By Jixie @20120208: I just added a few bars here to simulate the acceptance for SBS.
	//Need to add the full geometry later

	//horizontal bars
	double pTAStandHX=27.0*inch;
	double pTAStandHY=3.0*inch;
	double pTAStandHZ=1.5*inch;
	G4VSolid* TAStandHSolid = new G4Box("TAStandHBox",
		pTAStandHX/2.0, pTAStandHY/2.0, pTAStandHZ/2.0);

	G4LogicalVolume* TAStandHLogical = new G4LogicalVolume(TAStandHSolid,
		mMaterialManager->aluminum,"TAStandHLogical",0,0,0);
	TAStandHLogical->SetVisAttributes(GrayVisAtt); 

	double pPivot2Face=mPivot2ThirdArmFace+11.6*inch-pTAStandHZ/2.0;
	double pTAStandHPos_Y=-19.52*inch-pTAStandHY/2+mPivotYOffset;

	double pTAStandH1Pos_X=pPivot2Face*sin(mThirdArmAngle) + mPivotXOffset;
	double pTAStandH1Pos_Z=pPivot2Face*cos(mThirdArmAngle) + mPivotZOffset;
	new G4PVPlacement(pRotYThirdArmAngle,
		G4ThreeVector(pTAStandH1Pos_X,pTAStandHPos_Y,pTAStandH1Pos_Z),
		TAStandHLogical,"TAStandH1Phys",motherLogical,0,0,0);

	pPivot2Face=mPivot2ThirdArmFace+11.6*inch+3.0*inch+pTAStandHZ/2.0;
	double pTAStandH2Pos_X=pPivot2Face*sin(mThirdArmAngle) + mPivotXOffset;
	double pTAStandH2Pos_Z=pPivot2Face*cos(mThirdArmAngle) + mPivotZOffset;
	new G4PVPlacement(pRotYThirdArmAngle,
		G4ThreeVector(pTAStandH2Pos_X,pTAStandHPos_Y,pTAStandH2Pos_Z),
		TAStandHLogical,"TAStandH2Phys",motherLogical,0,0,0);

	//vertical poles
	double pTAStandVX=1.5*inch;
	double pTAStandVY=50.0*inch;
	double pTAStandVZ=3.0*inch;
	G4VSolid* TAStandVSolid = new G4Box("TAStandVBox",
		pTAStandVX/2.0, pTAStandVY/2.0, pTAStandVZ/2.0);

	G4LogicalVolume* TAStandVLogical = new G4LogicalVolume(TAStandVSolid,
		mMaterialManager->aluminum,"TAStandVLogical",0,0,0);
	TAStandVLogical->SetVisAttributes(GrayVisAtt); 

	pPivot2Face=mPivot2ThirdArmFace+11.6*inch+pTAStandVZ/2.0;
	double pTAStandVPos_Y=-33.14*inch+pTAStandVY/2 + mPivotYOffset;

	double pTAStandV1Pos_X= pPivot2Face*sin(mThirdArmAngle) + 
		(pTAStandHX/2+pTAStandVX/2)*cos(mThirdArmAngle) + mPivotXOffset;
	double pTAStandV1Pos_Z=pPivot2Face*cos(mThirdArmAngle) - 
		(pTAStandHX/2+pTAStandVX/2)*sin(mThirdArmAngle) + mPivotZOffset;
	new G4PVPlacement(pRotYThirdArmAngle,
		G4ThreeVector(pTAStandV1Pos_X,pTAStandVPos_Y,pTAStandV1Pos_Z),
		TAStandVLogical,"TAStandV1Phys",motherLogical,0,0,0);

	double pTAStandV2Pos_X= pPivot2Face*sin(mThirdArmAngle) - 
		(pTAStandHX/2+pTAStandVX/2)*cos(mThirdArmAngle) + mPivotXOffset;
	double pTAStandV2Pos_Z=pPivot2Face*cos(mThirdArmAngle) + 
		(pTAStandHX/2+pTAStandVX/2)*sin(mThirdArmAngle) + mPivotZOffset;
	new G4PVPlacement(pRotYThirdArmAngle,
		G4ThreeVector(pTAStandV2Pos_X,pTAStandVPos_Y,pTAStandV2Pos_Z),
		TAStandVLogical,"TAStandV2Phys",motherLogical,0,0,0);


	/////////////////////////
	// ThirdArm Container
	/////////////////////////

	//By Jixie: after open the SC, we found that it was only 6 inches
	double pThirdArmSC1X=6.0*inch;
	double pThirdArmSC1Y=24.0*inch;
	double pThirdArmSC1Z=0.27*inch;  //measured by Jixie

	double pThirdArmSC2X=17.0*inch;
	double pThirdArmSC2Y=22.0*inch;
	double pThirdArmSC2Z=2.0*inch;
	double pThirdArmSCGap=6.0*inch;

	double pThirdArmX=max(pThirdArmSC1X,pThirdArmSC2X);
	double pThirdArmY=max(pThirdArmSC1Y,pThirdArmSC2Y);
	double pThirdArmZ=pThirdArmSC1Z+pThirdArmSC2Z+pThirdArmSCGap;

	double pThirdArmPos_X=(mPivot2ThirdArmFace+pThirdArmZ/2.0)*sin(mThirdArmAngle)+
		mPivotXOffset+pThirdArmXoffset;
	double pThirdArmPos_Y=pThirdArmYoffset+mPivotYOffset;
	double pThirdArmPos_Z=(mPivot2ThirdArmFace+pThirdArmZ/2.0)*cos(mThirdArmAngle)+
		mPivotZOffset+pThirdArmZoffset;

	G4VSolid* thirdArmContainerSolid = new G4Box("thirdArmContainerBox",
		pThirdArmX/2.0+1.0*cm,pThirdArmY/2.0+1.0*cm,pThirdArmZ/2.0+5.0*cm);
	G4LogicalVolume* thirdArmContainerLogical = new G4LogicalVolume(thirdArmContainerSolid,
		mMaterialManager->air,"thirdArmContainerLogical",0,0,0);
	G4VPhysicalVolume* thirdArmContainerPhys=new G4PVPlacement(pRotThirdArm,
		G4ThreeVector(pThirdArmPos_X,pThirdArmPos_Y,pThirdArmPos_Z),
		thirdArmContainerLogical,"ThirdArmContainerPhys",motherLogical,0,0,0);
	thirdArmContainerLogical->SetVisAttributes(HallVisAtt);  


	/////////////////////////////////////
	//scintillator planes
	/////////////////////////////////////

	//dE plane
	G4VSolid* thirdArmSC1Solid = new G4Box("ThirdArmSC1Box",pThirdArmSC1X/2.0,
		pThirdArmSC1Y/2.0,pThirdArmSC1Z/2.0);
	G4LogicalVolume* thirdArmSC1Logical = new G4LogicalVolume(thirdArmSC1Solid,
		mMaterialManager->scintillator,"ThirdArmSC1Logical",0,0,0);
	new G4PVPlacement(0,G4ThreeVector(0,0,-pThirdArmSC1Z/2.-pThirdArmSCGap/2.),
		thirdArmSC1Logical,"ThirdArmSC1Phys",thirdArmContainerLogical,0,0,0);
	SDman->AddNewDetector(thirdArmSC1SD);
	thirdArmSC1Logical->SetSensitiveDetector(thirdArmSC1SD);
	thirdArmSC1Logical->SetVisAttributes(YellowGreenVisAtt);


	//E plane
	G4VSolid* thirdArmSC2Solid = new G4Box("ThirdArmSC2Box",pThirdArmSC2X/2.0,
		pThirdArmSC2Y/2.0,pThirdArmSC2Z/2.0);
	G4LogicalVolume* thirdArmSC2Logical = new G4LogicalVolume(thirdArmSC2Solid,
		mMaterialManager->scintillator,"ThirdArmSC2Logical",0,0,0);
	new G4PVPlacement(0,G4ThreeVector(0,0,pThirdArmSC2Z/2.+pThirdArmSCGap/2.),
		thirdArmSC2Logical,"ThirdArmSC2Phys",thirdArmContainerLogical,0,0,0);
	SDman->AddNewDetector(thirdArmSC2SD);
	thirdArmSC2Logical->SetSensitiveDetector(thirdArmSC2SD);
	thirdArmSC2Logical->SetVisAttributes(YellowGreenVisAtt);


	/////////////////////////////////////
	//virtual detector
	/////////////////////////////////////

	if(mSetupThirdArmVD)
	{
		double pThirdArmVDX=min(pThirdArmSC1X,pThirdArmSC2X);
		double pThirdArmVDY=min(pThirdArmSC1Y,pThirdArmSC2Y);
		double ThirdArmVDThick=1.0*mm;
		G4VSolid* thirdArmVDSolid = new G4Box("ThirdArmVDBox",
			pThirdArmVDX/2.0, pThirdArmVDY/2.0, ThirdArmVDThick/2.0);
		G4LogicalVolume* thirdArmVDLogical = new G4LogicalVolume(thirdArmVDSolid, 
			mMaterialManager->air,"ThirdArmVDLogical",0,0,0);

		G4ThreeVector pVDCenter;
		if(mSetupThirdArmVD==1)
		{
			//0.5 cm in front of the 1st SC plane
			pVDCenter.set(0,0,-pThirdArmSC1Z-pThirdArmSCGap/2.-ThirdArmVDThick/2.-0.5*cm);
		}
		else if(mSetupThirdArmVD==2)
		{
			//0.5 cm behind the 2nd SC plane
			pVDCenter.set(0,0,pThirdArmSC2Z+pThirdArmSCGap/2.+ThirdArmVDThick/2.+0.5*cm);
		}
		else
		{
			//0.5 cm behind the 1st SC plane
			pVDCenter.set(0,0,-pThirdArmSCGap/2.+ThirdArmVDThick/2.+0.5*cm);
		}

		new G4PVPlacement(0,pVDCenter,
			thirdArmVDLogical,"virtualBoundaryPhys",thirdArmContainerLogical,0,0,0);
		thirdArmVDLogical->SetVisAttributes(LightYellowVisAtt); 
	}

	return thirdArmContainerPhys;
}


/////////////////////////////////////////////////////////////////////
//By Jixie:
//There are 2 ways to place the FZ fields, A) using global field manager + field map
//B) placing the uniform field into a field logical volumn
//If there is overlaps, the field might now work.
//Note that one logical can take only onen field, therefore we have to build 2 identical logical
//if we use method B, while only one logical + 2 placements is needed in mentod A.
//One can update the field during run time in method A, while method B can only change it before run
//In this routine I choose method A
G4VPhysicalVolume* G2PDetectorConstruction::ConstructG2PChicane(G4LogicalVolume* motherLogical)
{
	const double inch=2.54*cm;
	// All G4VSensitiveDetector will be managed (deleted) by SDManager
	//therefore no need to delete it at the end of this subroutine
	G4String SDname;
	G4VSensitiveDetector* FZB1SD=new HRSStdSD(SDname="FZB1VD");
	G4VSensitiveDetector* FZB2SD=new HRSStdSD(SDname="FZB2VD");
	G4SDManager* SDman = G4SDManager::GetSDMpointer();


	int pUseDefaultFZB1=1,pUseDefaultFZB2=1;
	gConfig->GetArgument("UseDefaultFZB1",pUseDefaultFZB1);
	gConfig->GetArgument("UseDefaultFZB2",pUseDefaultFZB2); 


	///////////////////////////////////////////////////////
	//these parameters is from survey group
	//From the survey
	//From the survey result on 04/19/2012 (A1446.pdf)
	//This is for 1.159 GeV g2p
	//Coordinates WRT G2P target (mm)
	//Deltas in beam following system (mm)
	//Angular components from Optim Data (degrees)
	//run   name  x     y     z      dx   dy   dz    dyaw  dpitch  droll ideal_yaw ideal_pitch
	//1 MFZ1H05A  0.2 -72.6  6158.9 -0.2 -0.1   0.4 0.0041 -0.0017 -0.0032 142.500 -1.14599
	//1 ITV1H05  -0.2 -278.5 4409.5  0.2 -0.3 -11.5 0.0795  0.7345  0.2848 142.500 -8.97736
	//1 MFZ1H05B -0.1 -436.3 2660.0  0.1  0.0  -1.1 0.0287  0.0197 -0.0123 142.500  1.60003

	double pFZB1TiltedAngle=-1.146*deg,pFZB2TiltedAngle=1.6*deg;

	double pFZB1PosX=    0.2*mm+mPivotXOffset;
	double pFZB1PosY=  -72.6*mm+mPivotYOffset;
	double pFZB1PosZ=-6158.9*mm+mPivotZOffset;
	double pFZB2PosX=   -0.1*mm+mPivotXOffset;
	double pFZB2PosY= -436.3*mm+mPivotYOffset;
	double pFZB2PosZ=-2664.8*mm+mPivotZOffset;

	/*
	//from survey of commissioning: A1423
	double pFZB1PosX=0.11*mm;
	double pFZB1PosY=0.0*mm;
	double pFZB1PosZ=-7036*mm;
	double pFZB2PosX=0.04*mm;
	double pFZB2PosY=0.0*mm;
	double pFZB2PosZ=-3783*mm;
	*/

	double pFZB1Bx=0,pFZB1By=0,pFZB1Bz=0;
	double pFZB2Bx=0,pFZB2By=0,pFZB2Bz=0;
	///////////////////////////////////////////////////////
	//Here is all the default settings. will be used based on 
	//the beam energy, HRS angle and target field setting


	//judge which setting it is
	bool bIsSeptaIn=(fabs(mLHRSAngle-mLSeptumAngle)>0.1*deg)?true:false;

	double pBeamEnergy=2.257;
	gConfig->GetArgument("BeamEnergy",pBeamEnergy); 

	double pHelm_CurrentRatio=1.0;
	gConfig->GetParameter("Helm_CurrentRatio",pHelm_CurrentRatio);
	
	int pIsG2P1OrGEP2=0;
	double pHeam_Rot1Angle=270*deg;
	gConfig->GetParameter("Helm_RotAngle1",pHeam_Rot1Angle);  pHeam_Rot1Angle*=deg;
	if(fabs(pHeam_Rot1Angle-270*deg)<5) pIsG2P1OrGEP2=1;
	else if(fabs(pHeam_Rot1Angle-354*deg)<5) pIsG2P1OrGEP2=2;

	/*B field in unit of KG
	B1     B2
	1.159 GeV/c:
	2.5T  90deg  HAdump     -3.02  7.02
	2.5T  7deg   HAdump     -0.36  0.85

	1.706 GeV/c
	7deg   2.5T  HAdump  -0.36  0.85
	90deg  2.5T  HAdump  -2.99  6.99

	2.257 GeV/c:
	2.5T  90deg  HAdump     -2.98   6.98
	2.5T  7deg   HAdump     -0.36   0.85
	5.01T 90deg   Local     -2.94   6.93
	5.01T 7deg   HAdump     -0.73   1.70

	3.355 GeV/c
	5.01T 90deg  Local     -2.94   6.93
	5.01T 7deg   HAdump    -0.73   1.70

	WITHOUT SEPTA (TARGET AT PIVOT)

	1.159 GeV/c
	2.5T  90deg  Local     -1.97   3.95
	2.257 GeV/c
	5.01T 90deg  Local     -3.94   7.92
	3.355 GeV/c
	5.01T 90deg  Local     -3.94   7.92
	*/
	////////////////////////////////////////////////////////
	if(pUseDefaultFZB1)
	{
		pFZB1TiltedAngle=-1.146*deg;
		pFZB1PosX=    0.2*mm+mPivotXOffset;
		pFZB1PosY=  -72.6*mm+mPivotYOffset;
		pFZB1PosZ=-6159.3*mm+mPivotZOffset;
		//for all setups, the only change is the y position
		// Default setting is only good for settings that septum is in. 
		if(pIsG2P1OrGEP2==1)
		{
			pFZB1PosY        = -7.26*cm*1.159/pBeamEnergy*pHelm_CurrentRatio/0.5+mPivotYOffset;
			pFZB1Bx          = 0.302*tesla*pHelm_CurrentRatio/0.5;
		}
		else if(pIsG2P1OrGEP2==2)
		{
			pFZB1PosY        = -7.26*cm*1.159/pBeamEnergy*pHelm_CurrentRatio/1.0+mPivotYOffset;
			pFZB1Bx          = 0.036*tesla*pHelm_CurrentRatio/1.0;
		}

		//I have to updated these parameters in the map so they can be written into the config tree
		gConfig->SetArgument("FZB1TiltedAngle",pFZB1TiltedAngle/rad); 
		gConfig->SetArgument("FZB1PosX",pFZB1PosX/mm); 
		gConfig->SetArgument("FZB1PosY",pFZB1PosY/mm); 
		gConfig->SetArgument("FZB1PosZ",pFZB1PosZ/mm); 
		gConfig->SetArgument("FZB1Bx",pFZB1Bx/tesla); 
		gConfig->SetArgument("FZB1By",pFZB1By/tesla); 
		gConfig->SetArgument("FZB1Bz",pFZB1Bz/tesla); 
	}
	else
	{
		//reading the value from input arguments of option -FZB1 or -ChicaneMagnet1
		gConfig->GetArgument("FZB1TiltedAngle",pFZB1TiltedAngle); pFZB1TiltedAngle*=deg;
		gConfig->GetArgument("FZB1PosX",pFZB1PosX); pFZB1PosX*=mm;
		gConfig->GetArgument("FZB1PosY",pFZB1PosY); pFZB1PosY*=mm;
		gConfig->GetArgument("FZB1PosZ",pFZB1PosZ); pFZB1PosZ*=mm;
		gConfig->GetArgument("FZB1Bx",pFZB1Bx); pFZB1Bx*=tesla;
		gConfig->GetArgument("FZB1By",pFZB1By); pFZB1By*=tesla;
		gConfig->GetArgument("FZB1Bz",pFZB1Bz); pFZB1Bz*=tesla;
	}

	////////////////////////////////////////////////////////
	if(pUseDefaultFZB2)
	{
		if(pUseDefaultFZB2==99)
		{
			//HRS=12.5 degree,
			if(!bIsSeptaIn)
			{
				// TargetFiled=2.5T, Beam=1.159GeV
				if(fabs(pHelm_CurrentRatio)<0.7 && fabs(pBeamEnergy-1.159)<0.3)
				{
					if(pIsG2P1OrGEP2==1)
					{
						pFZB2TiltedAngle = 11.94*deg;
						pFZB2PosY        = -43.63*cm+mPivotYOffset;
						pFZB2Bx          = -0.702*tesla;
					}
					else if(pIsG2P1OrGEP2==2)
					{
						pFZB2TiltedAngle = 1.9*deg;
						pFZB2PosY        = -7.36*cm+mPivotYOffset;
						pFZB2Bx          = -0.085*tesla;
					}
				}
				// TargetFiled=5T, Beam=2.257GeV
				else if(fabs(pHelm_CurrentRatio)>0.9 && fabs(pBeamEnergy-2.257)<0.3) 
				{
					if(pIsG2P1OrGEP2==1)
					{
						pFZB2TiltedAngle = 5.1*deg;
						pFZB2PosY        = -20*cm+mPivotYOffset;
						pFZB2Bx          = -0.693*tesla;
					}
					else if(pIsG2P1OrGEP2==2)
					{
						pFZB2TiltedAngle = 1.48*deg;
						pFZB2PosY        = -10*cm+mPivotYOffset;
						pFZB2Bx          = -0.17*tesla;
					}
				}
				// TargetFiled=5T, Beam=3.359GeV
				else if(fabs(pHelm_CurrentRatio)>0.9 && fabs(pBeamEnergy-3.359)<0.3)
				{
					if(pIsG2P1OrGEP2==1)
					{
						pFZB2TiltedAngle = 3.5*deg;
						pFZB2PosY       = -13.42*cm+mPivotYOffset;
						pFZB2Bx          = -0.693*tesla;
					}
					else if(pIsG2P1OrGEP2==2)
					{
						pFZB2TiltedAngle = 0.99*deg;
						pFZB2PosY       = -12*cm+mPivotYOffset;
						pFZB2Bx          = -0.17*tesla;
					}
				}
			}
			else
			{
				////////////for 6 deg HRS
				// TargetFiled=2.5T, Beam=1.159GeV
				if(fabs(pHelm_CurrentRatio)<0.7 && fabs(pBeamEnergy-1.159)<0.3)
				{
					if(pIsG2P1OrGEP2==1)
					{
						//pFZB2TiltedAngle = 11.94*deg;
						pFZB2PosY        = -44*cm+mPivotYOffset;
						pFZB2Bx          = -0.702*tesla;
					}
					else if(pIsG2P1OrGEP2==2)
					{
						pFZB2TiltedAngle = 1.44*deg;
						pFZB2PosY        = -4*cm+mPivotYOffset;
						pFZB2Bx          = -0.085*tesla;
					}
				}
				// TargetFiled=2.5T, Beam=1.706GeV
				if(fabs(pHelm_CurrentRatio)<0.7 && fabs(pBeamEnergy-1.706)<0.3)
				{
					if(pIsG2P1OrGEP2==1)
					{
						pFZB2TiltedAngle = 8.06*deg;
						pFZB2PosY        = -28*cm+mPivotYOffset;
						pFZB2Bx          = -0.699*tesla;
					}
					else if(pIsG2P1OrGEP2==2)
					{
						pFZB2TiltedAngle = 0.98*deg;
						pFZB2PosY        = -4*cm+mPivotYOffset;
						pFZB2Bx          = -0.085*tesla;
					}
				}
				// TargetFiled=2.5T, Beam=2.257GeV
				if(fabs(pHelm_CurrentRatio)<0.7 && fabs(pBeamEnergy-2.257)<0.3)
				{
					if(pIsG2P1OrGEP2==1)
					{
						pFZB2TiltedAngle = 6.08*deg;
						pFZB2PosY        = -20*cm+mPivotYOffset;
						pFZB2Bx          = -0.699*tesla;
					}
					else if(pIsG2P1OrGEP2==2)
					{
						pFZB2TiltedAngle = 0.74*deg;
						pFZB2PosY        = -4*cm+mPivotYOffset;
						pFZB2Bx          = -0.085*tesla;
					}
				}
				// TargetFiled=5T, Beam=2.257GeV
				else if(fabs(pHelm_CurrentRatio)>0.9 && fabs(pBeamEnergy-2.257)<0.3) 
				{
					if(pIsG2P1OrGEP2==1)
					{
						pFZB2TiltedAngle = 6.08*deg;
						pFZB2PosY        = -18*cm+mPivotYOffset;
						pFZB2Bx          = -0.693*tesla;
					}
					else if(pIsG2P1OrGEP2==2)
					{
						pFZB2TiltedAngle = 1.48*deg;
						pFZB2PosY        = -4*cm+mPivotYOffset;
						pFZB2Bx          = -0.17*tesla;
					}
				}
				// TargetFiled=5T, Beam=3.359GeV
				else if(fabs(pHelm_CurrentRatio)>0.9 && fabs(pBeamEnergy-3.359)<0.3)
				{
					if(pIsG2P1OrGEP2==1)
					{
						pFZB2TiltedAngle = 4.09*deg;
						pFZB2PosY       = -10*cm+mPivotYOffset;
						pFZB2Bx          = -0.693*tesla;
					}
					else if(pIsG2P1OrGEP2==2)
					{
						pFZB2TiltedAngle = 0.99*deg;
						pFZB2PosY       = -4*cm+mPivotYOffset;
						pFZB2Bx          = -0.17*tesla;
					}
				}
			}
		}
		else
		{
			//for all setups, the only change is the y position
			// Default setting is only good for settings that that septum is in. 
			if(pIsG2P1OrGEP2==1)
			{
				pFZB2PosY = -43.63*cm*1.159/pBeamEnergy*pHelm_CurrentRatio/0.5+5*cm+mPivotYOffset;
				pFZB2Bx   = 0.702*tesla*pHelm_CurrentRatio/0.5;
			}
			else if(pIsG2P1OrGEP2==2)
			{
				pFZB2PosY = -7.36*cm*1.159/pBeamEnergy*pHelm_CurrentRatio/1.0+mPivotYOffset;
				pFZB2Bx   = 0.085*tesla*pHelm_CurrentRatio/1.0;
			}
		}

		//I have to updated these parameters in the map so they can be written into the config tree
		gConfig->SetArgument("FZB2TiltedAngle",pFZB2TiltedAngle/rad); 
		gConfig->SetArgument("FZB2PosX",pFZB2PosX/mm); 
		gConfig->SetArgument("FZB2PosY",pFZB2PosY/mm); 
		gConfig->SetArgument("FZB2PosZ",pFZB2PosZ/mm); 
		gConfig->SetArgument("FZB2Bx",pFZB2Bx/tesla); 
		gConfig->SetArgument("FZB2By",pFZB2By/tesla); 
		gConfig->SetArgument("FZB2Bz",pFZB2Bz/tesla); 
	}
	else
	{
		//reading the value from input arguments of option -FZB2 or -ChicaneMagnet2
		gConfig->GetArgument("FZB2TiltedAngle",pFZB2TiltedAngle); pFZB2TiltedAngle*=deg;
		gConfig->GetArgument("FZB2PosX",pFZB2PosX); pFZB2PosX*=mm;
		gConfig->GetArgument("FZB2PosY",pFZB2PosY); pFZB2PosY*=mm;
		gConfig->GetArgument("FZB2PosZ",pFZB2PosZ); pFZB2PosZ*=mm;
		gConfig->GetArgument("FZB2Bx",pFZB2Bx); pFZB2Bx*=tesla;
		gConfig->GetArgument("FZB2By",pFZB2By); pFZB2By*=tesla;
		gConfig->GetArgument("FZB2Bz",pFZB2Bz); pFZB2Bz*=tesla;
	}
	////////////////////////////////////////////////////////

	//now build the vectors
	G4ThreeVector pFZB1Pos3V(pFZB1PosX,pFZB1PosY,pFZB1PosZ);
	G4ThreeVector pFZB1Field3V(pFZB1Bx,pFZB1By,pFZB1Bz);

	G4ThreeVector pFZB2Pos3V(pFZB2PosX,pFZB2PosY,pFZB2PosZ);
	G4ThreeVector pFZB2Field3V(pFZB2Bx,pFZB2By,pFZB2Bz);


	G4RotationMatrix *pRotFZB1=new G4RotationMatrix();
	pRotFZB1->rotateX(pFZB1TiltedAngle); 
	G4RotationMatrix *pRotFZB2=new G4RotationMatrix();
	pRotFZB2->rotateX(pFZB2TiltedAngle); 

	//built the field for the chicane 
	HRSEMFieldSetup* mEMFieldSetup=HRSEMFieldSetup::GetHRSEMFieldSetup();
	mEMFieldSetup->SetBField3VFZB1(pFZB1Field3V); 
	G4FieldManager* FZB1FieldManager=mEMFieldSetup->GetFieldManagerFZB1();

	mEMFieldSetup->SetBField3VFZB2(pFZB2Field3V); 
	G4FieldManager* FZB2FieldManager=mEMFieldSetup->GetFieldManagerFZB2();

	//set the step limit
	double pFZBStepLimit=10;
	gConfig->GetArgument("FZBStepLimit",pFZBStepLimit); 
	pFZBStepLimit*=mm;
	G4UserLimits* uFZBStepLimits = new G4UserLimits(pFZBStepLimit);


	////////////////////////////
	//The field containner
	////////////////////////////
	double pFZBX=17.66*inch, pFZBY=15.37*2*inch, pFZBZ=76.34*inch;

	//stainless steel rectangle pipe, thickness 0.188", there is one flange attached to each end 
	double pFZBVacuumXout=1.65*inch;
	double pFZBVacuumYout=10.0*inch;
	double pFZBVacuumZ=89.51*inch;
	double pFZBVacuumXin=pFZBVacuumXout-2*0.188*inch;
	double pFZBVacuumYin=pFZBVacuumYout-2*0.188*inch;

	//the field container is only the area in between the mid side plates, 1.66" x 11.5" x 76.34"

	double pFZBFieldContainerX=pFZBVacuumXout;
	double pFZBFieldContainerY=11.5*inch;
	double pFZBFieldContainerZ=pFZBZ;
	G4VSolid* FZBFieldContainerSolid = new G4Box("FZBFieldContainerBox",
		pFZBFieldContainerX/2.0,pFZBFieldContainerY/2.0,pFZBFieldContainerZ/2.0);

	//build 2 logical volumn because each one has aits own field manager 
	G4LogicalVolume* FZB1FieldContainerLogical = new G4LogicalVolume(FZBFieldContainerSolid,
		mMaterialManager->vacuum,"FZB1FieldContainerLogical",FZB1FieldManager,0,uFZBStepLimits);
	FZB1FieldContainerLogical->SetVisAttributes(HallVisAtt);  

	G4LogicalVolume* FZB2FieldContainerLogical = new G4LogicalVolume(FZBFieldContainerSolid,
		mMaterialManager->vacuum,"FZB2FieldContainerLogical",FZB2FieldManager,0,uFZBStepLimits);
	FZB2FieldContainerLogical->SetVisAttributes(HallVisAtt);  

	G4VPhysicalVolume* FZB1FieldContainerPhys=new G4PVPlacement(pRotFZB1,pFZB1Pos3V,
		FZB1FieldContainerLogical,"FZB1FieldContainerPhys",motherLogical,0,0,0);

	new G4PVPlacement(pRotFZB2,pFZB2Pos3V,
		FZB2FieldContainerLogical,"FZB2FieldContainerPhys",motherLogical,0,0,0);

	//return FZB1FieldContainerPhys;   //for debug

	/////////////////////////
	// FZB magnet Container
	/////////////////////////

	//in order to avoid overlapping, I have to dig the field containner out from this container
	double pFZBContainerX=pFZBX+2*cm, pFZBContainerY=pFZBY+2*cm, pFZBContainerZ=pFZBZ+15*inch;

	G4VSolid* FZBContainerWholeSolid = new G4Box("FZBContainerWholeBox",
		pFZBContainerX/2.0,pFZBContainerY/2.0,pFZBContainerZ/2.0);

	G4SubtractionSolid *FZBContainerSolid = new G4SubtractionSolid(
		"FZBContainerSolid",FZBContainerWholeSolid,FZBFieldContainerSolid);

	G4LogicalVolume* FZBContainerLogical = new G4LogicalVolume(FZBContainerSolid,
		mMaterialManager->vacuum,"FZBContainerLogical",0,0,0);
	FZBContainerLogical->SetVisAttributes(HallVisAtt);  

	//build one logical volumn but place 2 copies of it
	new G4PVPlacement(pRotFZB1,pFZB1Pos3V,
		FZBContainerLogical,"FZB1ContainerPhys",motherLogical,0,0,0);

	new G4PVPlacement(pRotFZB2,pFZB2Pos3V,
		FZBContainerLogical,"FZB2ContainerPhys",motherLogical,0,0,0);


	/////////////////////////////////////
	//FZB vacuum pipe and flanges
	/////////////////////////////////////

	////stainless steel rectangle pipe, thickness 0.188", there is one flange attached to each end	
	//I have to chop the pipe into 3 pieces, the middle part will be placed in the field container
	//two ends will be placed in the FZB containner

	//has been declared above
	//double pFZBVacuumXout=1.65*inch;
	//double pFZBVacuumYout=10.0*inch;
	//double pFZBVacuumZ=89.51*inch;
	//double pFZBVacuumXin=pFZBVacuumXout-2*0.188*inch;
	//double pFZBVacuumYin=pFZBVacuumYout-2*0.188*inch;

	//the outer box
	G4VSolid* FZBVacuumOutSolid = new G4Box("FZBVacuumOutBox",
		pFZBVacuumXout/2.0,pFZBVacuumYout/2.0,pFZBVacuumZ/2.0);

	//the out box, only the mid part that inside the field containner
	G4VSolid* FZBVacuumOutMidSolid = new G4Box("FZBVacuumOutMidBox",
		pFZBVacuumXout/2.0,pFZBVacuumYout/2.0,pFZBZ/2.0);

	//the inner box
	G4VSolid* FZBVacuumInSolid = new G4Box("FZBVacuumInBox",
		pFZBVacuumXin/2.0,pFZBVacuumYin/2.0,pFZBVacuumZ/2.0+1.0*mm);


	//the whole vacuum pipe = outer - inner
	G4SubtractionSolid *FZBVacuumPipeSolid = new G4SubtractionSolid("FZBVacuumPipeSolid",
		FZBVacuumOutSolid,FZBVacuumInSolid);

	//the mid part of the vacuum pipe = outer - inner
	G4SubtractionSolid *FZBVacuumPipeMidSolid = new G4SubtractionSolid("FZBVacuumPipeMidSolid",
		FZBVacuumOutMidSolid,FZBVacuumInSolid);

	//the 2 end of the vacuum pipe = whole pipe - mid
	G4SubtractionSolid *FZBVacuumPipe2EndsSolid = new G4SubtractionSolid("FZBVacuumPipe2EndsSolid",
		FZBVacuumPipeSolid,FZBVacuumOutMidSolid);


	////////////////////////////////////
	//the end disk (flange)
	////////////////////////////////////
	double pFZBVacuumFlangeR=6.0*inch;
	double pFZBVacuumFlangeZ=1.0*inch;  //TODO: verify this thickness, it is very important
	G4VSolid* FZBVacuumFlangeSolid = new G4Tubs("FZBVacuumFlangeTubs",0,
		pFZBVacuumFlangeR,pFZBVacuumFlangeZ/2,0*deg,360*deg);

	//this flange will be subtracted by a rectange as large as the the vacuum pipe, 
	//then the vacuum pipe will be unioned with 2 flanges to make the total length to 90 inch 
	G4SubtractionSolid *FZBVacuumFlangeSubRecSolid = new G4SubtractionSolid(
		"FZBVacuumFlangeSubRecSolid",FZBVacuumFlangeSolid,FZBVacuumOutSolid);


	//-------------------------------------------------------------------------------------
	////Note: Polyhedron not available for FZBVacuumPipeUnionFlangeUpNDownSolid.
	////This means it cannot be visualized on most systems.
	////Therefore I can not union the 2 ends with their flanges
	////In order to see the flange and vacuum pipe, I have to place them piece by piece
	//
	////pipe union 2 flanges
	//double pFZBVacuumFlangePosZ=90*inch/2-pFZBVacuumFlangeZ/2;
	//G4UnionSolid *FZBVacuumPipeUnionFlangeUpSolid = new G4UnionSolid(
	//	"FZBVacuumPipeUnionFlangeUpSolid",FZBVacuumPipe2EndsSolid,FZBVacuumFlangeSubRecSolid,
	//	0,G4ThreeVector(0,0,-pFZBVacuumFlangePosZ));
	//G4UnionSolid *FZBVacuumPipeUnionFlangeUpNDownSolid = new G4UnionSolid(
	//	"FZBVacuumPipeUnionFlangeUpNDownSolid",FZBVacuumPipeUnionFlangeUpSolid,FZBVacuumFlangeSubRecSolid,
	//	0,G4ThreeVector(0,0,pFZBVacuumFlangePosZ));
	//
	//G4LogicalVolume* FZBVacuumLogical = new G4LogicalVolume(FZBVacuumPipeUnionFlangeUpNDownSolid,
	//	mMaterialManager->stainlesssteel,"FZBVacuumLogical",0,0,0);
	//FZBVacuumLogical->SetVisAttributes(PurpleVisAtt);  
	//
	////place this vacuum assembly into the magnet
	//new G4PVPlacement(0,G4ThreeVector(),
	//	FZBVacuumLogical,"FZBVacuumPhys",FZBContainerLogical,false,0,0);
	//Note: Polyhedron not available for FZBVacuumPipeUnionFlangeUpNDownSolid.
	//This means it cannot be visualized on most systems.
	//In order to see the flange and vacuum pipe, I have to place them piece by piece
	//-------------------------------------------------------------------------------------

	//place the 2 ends part of the vacuum pipe into the magnet containner
	G4LogicalVolume* FZBVacuumPipe2EndsLogical = new G4LogicalVolume(FZBVacuumPipe2EndsSolid,
		mMaterialManager->stainlesssteel,"FZBVacuumPipe2EndsLogical",0,0,0);
	FZBVacuumPipe2EndsLogical->SetVisAttributes(SilverVisAtt);  

	new G4PVPlacement(0,G4ThreeVector(),FZBVacuumPipe2EndsLogical,
		"FZBVacuumPipe2EndsPhys",FZBContainerLogical,false,0,0);


	//place 2 flanges, one at each end into the magnet container
	G4LogicalVolume* FZBVacuumFlangeLogical = new G4LogicalVolume(FZBVacuumFlangeSubRecSolid,
		mMaterialManager->stainlesssteel,"FZBVacuumFlangeLogical",0,0,0);
	FZBVacuumFlangeLogical->SetVisAttributes(GrayVisAtt);  

	double pFZBVacuumFlangePosZ=90*inch/2-pFZBVacuumFlangeZ/2;
	new G4PVPlacement(0,G4ThreeVector(0,0,pFZBVacuumFlangePosZ),
		FZBVacuumFlangeLogical,"FZBVacuumFlangeDownPhys",FZBContainerLogical,true,0,0);
	new G4PVPlacement(0,G4ThreeVector(0,0,-pFZBVacuumFlangePosZ),
		FZBVacuumFlangeLogical,"FZBVacuumFlangeUpPhys",FZBContainerLogical,true,1,0);

	///////////////////////////////////////////////////////////////
	//place the mid part of the vacuum pipe into the field containner
	G4LogicalVolume* FZBVacuumPipeMidLogical = new G4LogicalVolume(FZBVacuumPipeMidSolid,
		mMaterialManager->stainlesssteel,"FZBVacuumPipeLogical",0,0,0);
	FZBVacuumPipeMidLogical->SetVisAttributes(SilverVisAtt);  

	new G4PVPlacement(0,G4ThreeVector(),FZBVacuumPipeMidLogical,
		"FZB1VacuumPipeMidPhys",FZB1FieldContainerLogical,false,0,0);
	new G4PVPlacement(0,G4ThreeVector(),FZBVacuumPipeMidLogical,
		"FZB2VacuumPipeMidPhys",FZB2FieldContainerLogical,false,0,0);



	/////////////////////////////////////
	//FZB block, the silicon steel part
	/////////////////////////////////////
	//The FZ magnet = whole - center plate + 2 middle side plates

	double pFZBSidePlateX=5.85*inch;
	G4VSolid* FZBWholeSolid = new G4Box("FZBWholeBox",pFZBX/2.0,pFZBY/2.0,pFZBZ/2.0);

	//the top and the bottom box, make it longer for subtraction
	double pFZBSubPlateX=pFZBX-2*pFZBSidePlateX;
	double pFZBSubPlateY=pFZBY-2*5.98*inch;
	G4VSolid* FZBSubPlateSolid = new G4Box("FZBSubPlateBox",
		pFZBSubPlateX/2.0,pFZBSubPlateY/2.0,pFZBZ/2.0+1*mm);

	//the middle side plate
	double pFZBPlateMidSideX=pFZBX/2-pFZBVacuumXout/2-pFZBSidePlateX;
	double pFZBPlateMidSideY=11.5*inch;
	G4VSolid* FZBPlateMidSideSolid = new G4Box("FZBPlateMidSideBox",
		pFZBPlateMidSideX/2.0,pFZBPlateMidSideY/2.0,pFZBZ/2.0);

	//subtract the middle rectangle
	G4SubtractionSolid* FZBSubRecSolid=new G4SubtractionSolid("FZBSubRecSolid",
		FZBWholeSolid,FZBSubPlateSolid);

	//Union 2 middle side rectangles
	double pFZBPlateMidSidePosX=pFZBVacuumXout/2+pFZBPlateMidSideX/2;
	G4ThreeVector pFZBMidSideLPos3V=G4ThreeVector( pFZBPlateMidSidePosX,0,0);
	G4ThreeVector pFZBMidSideRPos3V=G4ThreeVector(-pFZBPlateMidSidePosX,0,0);
	G4UnionSolid* FZBSubRecUnionMidLSolid=new G4UnionSolid("FZBSubRecUnionMidLSolid",
		FZBSubRecSolid,FZBPlateMidSideSolid,0,pFZBMidSideLPos3V);
	G4UnionSolid* FZBSubRecUnionMidLRSolid=new G4UnionSolid("FZBSubRecUnionMidLRSolid",
		FZBSubRecUnionMidLSolid,FZBPlateMidSideSolid,0,pFZBMidSideRPos3V);

	G4LogicalVolume* FZBSteelLogical = new G4LogicalVolume(FZBSubRecUnionMidLRSolid,
		mMaterialManager->siliconsteel,"FZBSteelLogical",0,0,0);
	FZBSteelLogical->SetVisAttributes(DarkBlueVisAtt);  

	//place one copy of this silicon steel into each magnet
	new G4PVPlacement(0,G4ThreeVector(),
		FZBSteelLogical,"FZBSteelPhys",FZBContainerLogical,false,0,0);

	/////////////////////////////////////
	//FZB copper coils
	/////////////////////////////////////
	//this part is hard to build,  I simplified it as a big box subtract the siliconsteel block
	//made of copperx

	//the whole box
	//each group of coil is 0.925" wide, place 2 of them with 1mm gap in between
	double pFZBCoilX=2*0.925*inch+1*mm; 
	double pFZBCoilY=pFZBY-2*5.98*inch; 
	double pFZBCoilZ=(78.96+3.35*2)*inch;  
	G4VSolid* FZBCoilWholeSolid = new G4Box("FZBCoilWholeBox",
		pFZBCoilX/2.0,pFZBCoilY/2.0,pFZBCoilZ/2.0);

	//the small box need to be subtracted
	double pFZBCoilSmallRecX=pFZBCoilX+1*mm;
	double pFZBCoilSmallRecY=pFZBPlateMidSideY+1*mm;
	double pFZBCoilSmallRecZ=78.96*inch;  
	G4VSolid* FZBCoilSmallRecSolid = new G4Box("FZBCoilSmallRecBox",
		pFZBCoilSmallRecX/2.0,pFZBCoilSmallRecY/2.0,pFZBCoilSmallRecZ/2.0);

	//coil = whole - small rectangle
	G4SubtractionSolid* FZBCoilSolid=new G4SubtractionSolid("FZBCoilSolid",
		FZBCoilWholeSolid,FZBCoilSmallRecSolid);

	//TODO: cut off the corner ......

	//build the logical volumn
	G4LogicalVolume* FZBCoilLogical = new G4LogicalVolume(FZBCoilSolid,
		mMaterialManager->copper,"FZBCoilLogical",0,0,0);
	FZBCoilLogical->SetVisAttributes(CuBrownVisAtt);  

	//place 2 coils in each FZ magnet
	double pFZBCoilPosX=pFZBX/2-pFZBSidePlateX-pFZBCoilX/2;  //attached to the sife plate
	new G4PVPlacement(0,G4ThreeVector(pFZBCoilPosX,0,0),
		FZBCoilLogical,"FZBCoilLeftPhys",FZBContainerLogical,true,0,0);
	new G4PVPlacement(0,G4ThreeVector(-pFZBCoilPosX,0,0),
		FZBCoilLogical,"FZBCoilRightPhys",FZBContainerLogical,true,1,0);


	/////////////////////////////////////
	//FZB copper coil spacers
	/////////////////////////////////////
	//there are 4 spacers are used to support|fixed the vacuum
	//They thickness and density are of important since they are in the way that bremstrulung
	//photons come out
	//In the engineering drawing, the material is aluminum alloy bond with contact cemont on the surface
	// I am not going to build the cement here

	double pFZBCoilSpacerX=1.75*inch;
	//pFZBCoilSpacerY is 3.25" in the drawing, it should be about 3.64" in order to fill up the gap, 
	//I think the rest is cement, Here I fill it up
	//double pFZBCoilSpacerY=3.25*inch;   
	double pFZBCoilSpacerY=pFZBCoilY/2-pFZBPlateMidSideY/2;
	double pFZBCoilSpacerZ=8.00*inch;
	G4VSolid* FZBCoilSpacerSolid = new G4Box("FZBCoilSpacerBox",
		pFZBCoilSpacerX/2.0,pFZBCoilSpacerY/2.0,pFZBCoilSpacerZ/2.0);

	G4LogicalVolume* FZBCoilSpacerLogical = new G4LogicalVolume(FZBCoilSpacerSolid,
		mMaterialManager->aluminum,"FZBCoilSpacerLogical",0,0,0);
	FZBCoilSpacerLogical->SetVisAttributes(LightBlueVisAtt);  

	//place 4 copies of this spacer at the corner in each FZ magnet
	double pFZBCoilPosY=pFZBPlateMidSideY/2+pFZBCoilSpacerY/2;
	double pFZBCoilPosZ=pFZBZ/2-pFZBCoilSpacerZ/2;
	G4ThreeVector pSpacer1Pos3V=G4ThreeVector(0, pFZBCoilPosY, pFZBCoilPosZ);
	G4ThreeVector pSpacer2Pos3V=G4ThreeVector(0,-pFZBCoilPosY, pFZBCoilPosZ);
	G4ThreeVector pSpacer3Pos3V=G4ThreeVector(0,-pFZBCoilPosY,-pFZBCoilPosZ);
	G4ThreeVector pSpacer4Pos3V=G4ThreeVector(0, pFZBCoilPosY,-pFZBCoilPosZ);
	new G4PVPlacement(0,pSpacer1Pos3V,
		FZBCoilSpacerLogical,"FZBCoilSpacerPhys",FZBContainerLogical,true,0,0);
	new G4PVPlacement(0,pSpacer2Pos3V,
		FZBCoilSpacerLogical,"FZBCoilSpacerPhys",FZBContainerLogical,true,1,0);
	new G4PVPlacement(0,pSpacer3Pos3V,
		FZBCoilSpacerLogical,"FZBCoilSpacerPhys",FZBContainerLogical,true,2,0);
	new G4PVPlacement(0,pSpacer4Pos3V,
		FZBCoilSpacerLogical,"FZBCoilSpacerPhys",FZBContainerLogical,true,3,0);


	//virtual detector
	/////////////////////////////////////

	//mSetupChicaneVD = 0 menas none, 1 is VB1VD, 2 is VB2VD, 3 is both
	if(mSetupChicaneVD)
	{
		double pFZBVDWidth=60.0*cm;
		double pFZBVDHeight=100.0*cm;
		double pFZBVDThick=1.0*cm;
		G4VSolid* FZB1VDSolid = new G4Box("FZB1VDBox",
			pFZBVDWidth/2.0, pFZBVDHeight/2.0, pFZBVDThick/2.0);
		G4LogicalVolume* FZB1VDLogical = new G4LogicalVolume(FZB1VDSolid, 
			mMaterialManager->air,"FZB1VDLogical",0,0,0);
		SDman->AddNewDetector(FZB1SD);
		FZB1VDLogical->SetSensitiveDetector(FZB1SD);
		FZB1VDLogical->SetVisAttributes(LightYellowVisAtt); 
		G4LogicalVolume* FZB2VDLogical = new G4LogicalVolume(FZB1VDSolid, 
			mMaterialManager->air,"FZB2VDLogical",0,0,0); 
		SDman->AddNewDetector(FZB2SD);
		FZB2VDLogical->SetSensitiveDetector(FZB2SD);
		FZB2VDLogical->SetVisAttributes(LightYellowVisAtt);

		G4ThreeVector pVD1Center;
		if(mSetupChicaneVD==1 || mSetupChicaneVD==3)
		{
			new G4PVPlacement(0,G4ThreeVector(0,0,150*cm+pFZB1PosZ),
				FZB1VDLogical,"FZB1VDPhys",motherLogical,0,0,0);
		}
		if(mSetupChicaneVD==2 || mSetupChicaneVD==3)
		{
			new G4PVPlacement(0,G4ThreeVector(0,0,150*cm+pFZB2PosZ),
				FZB2VDLogical,"FZB2VDPhys",motherLogical,0,0,0);
		}
	}

	return FZB1FieldContainerPhys;
}


G4VPhysicalVolume* G2PDetectorConstruction::ConstructG2PPlatform(G4LogicalVolume* motherLogical)
{
	const double inch=2.54*cm;
	G4RotationMatrix *pRotX90deg=new G4RotationMatrix();
	pRotX90deg->rotateX(-90*deg); 
	G4RotationMatrix *pRotY90deg=new G4RotationMatrix();
	pRotY90deg->rotateY(-90*deg); 

	/////////////////////////////////////
	//poles to support 3rd floor
	/////////////////////////////////////
	double PFPoleWidth=4.0*inch;	//x
	double PFPoleHeight=4.0*inch;   //y
	double PFPoleLength=48.0*inch;  //z
	G4VSolid* PFPoleSolid = new G4Box("PFPoleBox",PFPoleWidth/2.0,
		PFPoleHeight/2.0,PFPoleLength/2.0);

	G4LogicalVolume* PFPoleLogical = new G4LogicalVolume(PFPoleSolid,
		mMaterialManager->aluminum,"PFPoleLogical",0,0,0);
	PFPoleLogical->SetVisAttributes(LightGreenVisAtt); 

	double pPFPolePos_X=42.0*inch+PFPoleWidth/2.0;
	double pPFPolePos_Y=-1.0*inch;
	double pPFPoleUpPos_Z=-70.53*inch+PFPoleHeight/2.0;
	double pPFPoleDownPos_Z=32.0*inch+PFPoleHeight/2.0;
	G4VPhysicalVolume *PFPolePhys = new G4PVPlacement(pRotX90deg,
		G4ThreeVector(pPFPolePos_X,pPFPolePos_Y,pPFPoleUpPos_Z),
		PFPoleLogical,"PFPolePhys",motherLogical,true,0,0);	//up left
	new G4PVPlacement(pRotX90deg,
		G4ThreeVector(-pPFPolePos_X,pPFPolePos_Y,pPFPoleUpPos_Z),
		PFPoleLogical,"PFPolePhys",motherLogical,true,1,0);	//up right
	new G4PVPlacement(pRotX90deg,
		G4ThreeVector(-pPFPolePos_X,pPFPolePos_Y,pPFPoleDownPos_Z),
		PFPoleLogical,"PFPolePhys",motherLogical,true,0,0);	//down right
	new G4PVPlacement(pRotX90deg,
		G4ThreeVector(pPFPolePos_X,pPFPolePos_Y,pPFPoleDownPos_Z),
		PFPoleLogical,"PFPolePhys",motherLogical,true,0,0); //down left


	double pSideBoardWidth=0.75*inch;
	double pSideBoardHeight=4.0*inch;


	int pSetupPlatform=0;
	gConfig->GetParameter("SetupPlatform",pSetupPlatform);
	if(pSetupPlatform==1 || pSetupPlatform==3)
	{
		/////////////////////////////////////
		//bars to support 3rd floor
		/////////////////////////////////////
		double PF3FBarWidth=4.0*inch;	//x
		double PF3FBarHeight=4.0*inch;   //y
		double PF3FBarLRLength=98.53*inch;  //z
		double PF3FBarUDLength=170.0*inch;  //z
		G4VSolid* PF3FBarLRSolid = new G4Box("PF3FBarLRBox",PF3FBarWidth/2.0,
			PF3FBarHeight/2.0,PF3FBarLRLength/2.0);
		G4VSolid* PF3FBarUDSolid = new G4Box("PF3FBarUDBox",PF3FBarWidth/2.0,
			PF3FBarHeight/2.0,PF3FBarUDLength/2.0);

		G4LogicalVolume* PF3FBarLRLogical = new G4LogicalVolume(PF3FBarLRSolid,
			mMaterialManager->aluminum,"PF3FBarLRLogical",0,0,0);
		PF3FBarLRLogical->SetVisAttributes(LightGreenVisAtt); 
		G4LogicalVolume* PF3FBarUDLogical = new G4LogicalVolume(PF3FBarUDSolid,
			mMaterialManager->aluminum,"PF3FBarUDLogical",0,0,0);
		PF3FBarUDLogical->SetVisAttributes(LightGreenVisAtt); 

		double pPF3FBarUDPos_X=0.0;
		double pPF3FBarPos_Y=25.5*inch;
		double pPF3FBarUpPos_Z=-70.53*inch+PF3FBarHeight/2.0;
		double pPF3FBarDownPos_Z=32.0*inch+PF3FBarHeight/2.0;
		double pPF3FBarLRPos_X=42.0*inch+PF3FBarWidth/2.0;
		double pPF3FBarLRPos_Z=(pPF3FBarUpPos_Z+pPF3FBarDownPos_Z)/2;

		new G4PVPlacement(0,
			G4ThreeVector(pPF3FBarLRPos_X,pPF3FBarPos_Y,pPF3FBarLRPos_Z),
			PF3FBarLRLogical,"PF3FBarLPhys",motherLogical,true,0,0);	// left
		new G4PVPlacement(0,
			G4ThreeVector(-pPF3FBarLRPos_X,pPF3FBarPos_Y,pPF3FBarLRPos_Z),
			PF3FBarLRLogical,"PF3FBarRPhys",motherLogical,true,1,0);	// right
		new G4PVPlacement(pRotY90deg,
			G4ThreeVector(pPF3FBarUDPos_X,pPF3FBarPos_Y,pPF3FBarUpPos_Z),
			PF3FBarUDLogical,"PF3FBarUpPhys",motherLogical,true,0,0);	//up 
		new G4PVPlacement(pRotY90deg,
			G4ThreeVector(pPF3FBarUDPos_X,pPF3FBarPos_Y,pPF3FBarDownPos_Z),
			PF3FBarUDLogical,"PF3FBarDownPhys",motherLogical,true,0,0); //down 


		/////////////////////////////////////
		// 3rd floor
		/////////////////////////////////////

		double pPF3FX=170.0*inch;
		double pPF3FY=0.5*inch;
		double pPF3FZ=106.53*inch;
		//3rd floor, need to dig hole (R=31") for target chamber
		G4VSolid* PF3rdFloorWholeSolid = new G4Box("PF3rdFloorWholeBox",pPF3FX/2.0,
			pPF3FY/2.0,pPF3FZ/2.0);
		G4VSolid* PF3rdFloorSubCircleSolid = new G4Tubs("PF3rdFloorSubCircleSolid",0,
			39.0*inch/2,0.55*inch/2,0*deg,360*deg);
		G4SubtractionSolid *PF3rdFloorSolid = new G4SubtractionSolid(
			"PF3rdFloorSolid",PF3rdFloorWholeSolid,PF3rdFloorSubCircleSolid,pRotX90deg,
			G4ThreeVector(0,0,mScatChamberZOffset-pPF3FBarLRPos_Z));
		G4LogicalVolume* PF3rdFloorLogical = new G4LogicalVolume(PF3rdFloorSolid,
			mMaterialManager->aluminum,"PF3rdFloorLogical",0,0,0);
		PF3rdFloorLogical->SetVisAttributes(LightGreenVisAtt); 


		new G4PVPlacement(0,G4ThreeVector(0,27.5*inch,pPF3FBarLRPos_Z),
			PF3rdFloorLogical,"PF3rdFloorPhys",motherLogical,0,0,0);	


		/////////////////////////////////////
		// 3rd floor side boards
		/////////////////////////////////////

		G4VSolid* PF3FSideLRSolid = new G4Box("PF3FSideLRSolid",pSideBoardWidth/2.0,
			pSideBoardHeight/2.0,pPF3FZ/2.0);
		G4VSolid* PF3FSideUDSolid = new G4Box("PF3FBarUDBox",pSideBoardWidth/2.0,
			pSideBoardHeight/2.0,pPF3FX/2.0-pSideBoardWidth);

		G4LogicalVolume* PF3FSideLRLogical = new G4LogicalVolume(PF3FSideLRSolid,
			mMaterialManager->aluminum,"PF3FSideLRLogical",0,0,0);
		PF3FSideLRLogical->SetVisAttributes(LightGreenVisAtt); 
		G4LogicalVolume* PF3FSideUDLogical = new G4LogicalVolume(PF3FSideUDSolid,
			mMaterialManager->aluminum,"PF3FSideUDLogical",0,0,0);
		PF3FSideUDLogical->SetVisAttributes(LightGreenVisAtt); 


		double pPF3FSideUDPos_X=0.0;
		double pPF3FSidePos_Y=27.5*inch+pSideBoardHeight/2.0;
		double pPF3FSideUpPos_Z=-70.53*inch+pSideBoardWidth/2.0;
		double pPF3FSideDownPos_Z=36.0*inch-pSideBoardWidth/2.0;
		double pPF3FSideLRPos_X=pPF3FX/2-pSideBoardWidth/2.0;
		double pPF3FSideLRPos_Z=(pPF3FSideUpPos_Z+pPF3FSideDownPos_Z)/2;

		new G4PVPlacement(0,
			G4ThreeVector(pPF3FSideLRPos_X,pPF3FSidePos_Y,pPF3FSideLRPos_Z),
			PF3FSideLRLogical,"PF3FSideLPhys",motherLogical,true,0,0);	// left
		new G4PVPlacement(0,
			G4ThreeVector(-pPF3FSideLRPos_X,pPF3FSidePos_Y,pPF3FSideLRPos_Z),
			PF3FSideLRLogical,"PF3FSideRPhys",motherLogical,true,1,0);	// right
		new G4PVPlacement(pRotY90deg,
			G4ThreeVector(pPF3FSideUDPos_X,pPF3FSidePos_Y,pPF3FSideUpPos_Z),
			PF3FSideUDLogical,"PF3FSideUpPhys",motherLogical,true,0,0);	//up 
		new G4PVPlacement(pRotY90deg,
			G4ThreeVector(pPF3FSideUDPos_X,pPF3FSidePos_Y,pPF3FSideDownPos_Z),
			PF3FSideUDLogical,"PF3FSideDownPhys",motherLogical,true,0,0); //down 
	}

	/////////////////////////////////////
	//rails to support septum
	/////////////////////////////////////

	double PF2FRailWidth=4.0*inch;	//x
	double PF2FRailHeight=4.0*inch;   //y
	double PF2FRailLRLength=98.53*inch;  //z
	double PF2FRailUDLength=92.0*inch;  //z
	G4VSolid* PF2FRailLRSolid = new G4Box("PF2FRailLRBox",PF2FRailWidth/2.0,
		PF2FRailHeight/2.0,PF2FRailLRLength/2.0);
	G4VSolid* PF2FRailUDSolid = new G4Box("PF2FRailUDBox",PF2FRailWidth/2.0,
		PF2FRailHeight/2.0,PF2FRailUDLength/2.0);

	G4LogicalVolume* PF2FRailLRLogical = new G4LogicalVolume(PF2FRailLRSolid,
		mMaterialManager->aluminum,"PF2FRailLRLogical",0,0,0);
	PF2FRailLRLogical->SetVisAttributes(SteelVisAtt); 
	G4LogicalVolume* PF2FRailUDLogical = new G4LogicalVolume(PF2FRailUDSolid,
		mMaterialManager->stainlesssteel,"PF2FRailUDLogical",0,0,0);
	PF2FRailUDLogical->SetVisAttributes(SteelVisAtt); 

	double pPF2FRailUDPos_X=0.0;
	double pPF2FRailPos_Y=-27.5*inch;
	double pPF2FRailUpPos_Z=-70.53*inch+PF2FRailHeight/2.0;
	double pPF2FRailDownPos_Z=32.0*inch+PF2FRailHeight/2.0;
	double pPF2FRailLRPos_X=42.0*inch+PF2FRailWidth/2.0;
	double pPF2FRailLRPos_Z=(pPF2FRailUpPos_Z+pPF2FRailDownPos_Z)/2;

	new G4PVPlacement(0,
		G4ThreeVector(pPF2FRailLRPos_X,pPF2FRailPos_Y,pPF2FRailLRPos_Z),
		PF2FRailLRLogical,"PF2FRailLPhys",motherLogical,true,0,0);	// left
	new G4PVPlacement(0,
		G4ThreeVector(-pPF2FRailLRPos_X,pPF2FRailPos_Y,pPF2FRailLRPos_Z),
		PF2FRailLRLogical,"PF2FRailRPhys",motherLogical,true,1,0);	// right
	new G4PVPlacement(pRotY90deg,
		G4ThreeVector(pPF2FRailUDPos_X,pPF2FRailPos_Y,pPF2FRailUpPos_Z),
		PF2FRailUDLogical,"PF2FRailUpPhys",motherLogical,true,0,0);	//up 
	new G4PVPlacement(pRotY90deg,
		G4ThreeVector(pPF2FRailUDPos_X,pPF2FRailPos_Y,pPF2FRailDownPos_Z),
		PF2FRailUDLogical,"PF2FRailDownPhys",motherLogical,true,0,0); //down 


	if(pSetupPlatform==1 || pSetupPlatform==2)
	{
		/////////////////////////////////////
		// 2nd floor
		/////////////////////////////////////

		double pPF2FX=242.5*inch;
		double pPF2FY=0.5*inch;
		double pPF2FZ=106.53*inch;
		//3rd floor, need to dig hole (R=31") for target chamber
		G4VSolid* PF2ndFloorSolid = new G4Box("PF2ndFloorBox",pPF2FX/2.0,
			pPF2FY/2.0,pPF2FZ/2.0);

		G4LogicalVolume* PF2ndFloorLogical = new G4LogicalVolume(PF2ndFloorSolid,
			mMaterialManager->aluminum,"PF2ndFloorLogical",0,0,0);
		PF2ndFloorLogical->SetVisAttributes(LightGreenVisAtt); 


		new G4PVPlacement(0,G4ThreeVector(0,-38.38*inch,pPF2FRailLRPos_Z),
			PF2ndFloorLogical,"PF2ndFloorPhys",motherLogical,0,0,0);	


		/////////////////////////////////////
		// 2nd floor side boards
		/////////////////////////////////////

		G4VSolid* PF2FSideLRSolid = new G4Box("PF2FSideLRSolid",pSideBoardWidth/2.0,
			pSideBoardHeight/2.0,pPF2FZ/2.0);
		G4VSolid* PF2FSideUDSolid = new G4Box("PF2FBarUDBox",pSideBoardWidth/2.0,
			pSideBoardHeight/2.0,pPF2FX/2.0-pSideBoardWidth);

		G4LogicalVolume* PF2FSideLRLogical = new G4LogicalVolume(PF2FSideLRSolid,
			mMaterialManager->aluminum,"PF2FSideLRLogical",0,0,0);
		PF2FSideLRLogical->SetVisAttributes(LightGreenVisAtt); 
		G4LogicalVolume* PF2FSideUDLogical = new G4LogicalVolume(PF2FSideUDSolid,
			mMaterialManager->aluminum,"PF2FSideUDLogical",0,0,0);
		PF2FSideUDLogical->SetVisAttributes(LightGreenVisAtt); 


		double pPF2FSideUDPos_X=0.0;
		double pPF2FSidePos_Y=-38.13*inch+pSideBoardHeight/2.0;
		double pPF2FSideUpPos_Z=-70.53*inch+pSideBoardWidth/2.0;
		double pPF2FSideDownPos_Z=36.0*inch-pSideBoardWidth/2.0;
		double pPF2FSideLRPos_X=pPF2FX/2-pSideBoardWidth/2.0;
		double pPF2FSideLRPos_Z=(pPF2FSideUpPos_Z+pPF2FSideDownPos_Z)/2;

		new G4PVPlacement(0,
			G4ThreeVector(pPF2FSideLRPos_X,pPF2FSidePos_Y,pPF2FSideLRPos_Z),
			PF2FSideLRLogical,"PF2FSideLPhys",motherLogical,true,0,0);	// left
		new G4PVPlacement(0,
			G4ThreeVector(-pPF2FSideLRPos_X,pPF2FSidePos_Y,pPF2FSideLRPos_Z),
			PF2FSideLRLogical,"PF2FSideRPhys",motherLogical,true,1,0);	// right
		new G4PVPlacement(pRotY90deg,
			G4ThreeVector(pPF2FSideUDPos_X,pPF2FSidePos_Y,pPF2FSideUpPos_Z),
			PF2FSideUDLogical,"PF2FSideUpPhys",motherLogical,true,0,0);	//up 
		new G4PVPlacement(pRotY90deg,
			G4ThreeVector(pPF2FSideUDPos_X,pPF2FSidePos_Y,pPF2FSideDownPos_Z),
			PF2FSideUDLogical,"PF2FSideDownPhys",motherLogical,true,0,0); //down 
	}

	return PFPolePhys;
}

