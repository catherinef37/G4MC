// ********************************************************************
//
// $Id: HRSDetectorConstruction.cc,v 1.02, 2010/12/26 HRS Exp $
// --------------------------------------------------------------
//
// ********************************************************************
//
#include <stdio.h>
#include <math.h>
#include <fstream>

#include "HRSMaterial.hh"
#include "HRSDetectorConstruction.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

#include "G4UniformMagField.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "BField_Dipole.hh"
#include "BField_Dipole_Fringe.hh"

#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4Trap.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4UserLimits.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4Polycone.hh"
#include "G4AssemblyVolume.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"

#include "HRSEMFieldSetup.hh"
#include "HRSStdSD.hh"
#include "UsageManager.hh"

#include "CREXDetectorConstruction.hh"

#include "G4EqMagElectricField.hh"
#include "G4UniformElectricField.hh"
#include "G4ClassicalRK4.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ChordFinder.hh"
//To verify some geometry, I use this flag to place some unit
//just for debugging
//#define G4DEBUG_GEOMETRY 0

//The G4RotationMatrix rotates the whole coordinate system, 
//Looking from top, it always rotate clockwise  

extern UsageManager* gConfig;	


/////////////////////////////////////////////////////////////////////
HRSDetectorConstruction::HRSDetectorConstruction()
{
	this->GetConfig();
	mEMFieldSetup = 0;

	//construct the material manager, this line should behind GetConfig()
	//since it need to access the buffer of gConfig
	mMaterialManager=HRSMaterial::GetHRSMaterialManager();
}

HRSDetectorConstruction::~HRSDetectorConstruction()
{
	if(mEMFieldSetup) delete mEMFieldSetup;
}

/////////////////////////////////////////////////////////////////////
void HRSDetectorConstruction::DumpGeometricalTree(G4VPhysicalVolume* aPhysVol,G4int depth,
												  ostream &fout)
{
	//define the aLogVol to speed up
	G4LogicalVolume* aLogVol=aPhysVol->GetLogicalVolume();
	fout << "#";
	for(int isp=0;isp<depth;isp++) { fout << "--"; }
	fout << aPhysVol->GetName() << "[" << aPhysVol->GetCopyNo() << "] "
		<< aLogVol->GetName() << " "
		<< aLogVol->GetNoDaughters() << " "
		<< aLogVol->GetMaterial()->GetName();
	if(aLogVol->GetSensitiveDetector())
	{
		fout << " " << aLogVol->GetSensitiveDetector()->GetFullPathName();
	}
	fout << endl;
	for(int i=0;i<aLogVol->GetNoDaughters();i++)
	{ 
		DumpGeometricalTree(aLogVol->GetDaughter(i),depth+1,fout); 
	}
}

/////////////////////////////////////////////////////////////////////
void HRSDetectorConstruction::GenerateGeometryMacro(G4VPhysicalVolume* aPhysVol,G4int depth)
{
	ofstream fout;
	if(depth==0) 
	{
		fout.open("geometry.mac",ios_base::out);
		//DumpGeometricalTree(aPhysVol,depth,fout);
	}
	else fout.open("geometry.mac",ios_base::app);


	G4LogicalVolume* aLogVol = aPhysVol->GetLogicalVolume();
	string logName = aLogVol->GetName().c_str();

	//some phys volumn share the same logical volumn
	//each log vol needs to print cmds for just one time
	if(mIsLogVolPrinted.size() > 0)
	{
		if(mIsLogVolPrinted.find(logName) != mIsLogVolPrinted.end()) return;
	}
	//add this log vol into the map
	mIsLogVolPrinted[logName] = 1;

	int visibility = (aLogVol->GetVisAttributes())?1:0;
	if(visibility) visibility=aLogVol->GetVisAttributes()->IsVisible()?1:0;
	G4Colour aColor;
	if(visibility) aColor = aLogVol->GetVisAttributes()->GetColor();
	int nDaughter = aLogVol->GetNoDaughters();

	fout << "#";
	for(int isp=0;isp<depth;isp++) { fout << "--"; }
	fout << aPhysVol->GetName() << "[" << aPhysVol->GetCopyNo() << "] "
		<< " " << logName 
		<< " " << aLogVol->GetMaterial()->GetName()
		<< " Ndaughter=" << nDaughter;
	if(visibility) fout<< " G4Color"<< aColor;
	if(aLogVol->GetSensitiveDetector())
		fout << " SD=" << aLogVol->GetSensitiveDetector()->GetFullPathName();
	fout << endl;

	char space[100];
	sprintf(space," ");
	for(int isp=0;isp<depth;isp++) { sprintf(space,"%s  ",space); }

	//now print the commands
	fout<<space<<"/vis/geometry/set/visibility "<<logName<<" 0 "
		<<(visibility?"true":"false")<<" \n";
	fout<<space<<"/vis/geometry/set/forceSolid "<<logName<<" 0 "<<" true \n";
	fout<<space<<"#/vis/geometry/set/forceWireframe "<<logName<<" 0 "<<" false \n";
	if(nDaughter)
		fout<<space<<"#/vis/geometry/set/daughtersInvisible "<<logName<<" 0 "<<" false \n";
	if(visibility)
		fout<<space<<"#/vis/geometry/set/colour "<<logName<<" 0 "
		<<aColor.GetRed()<<" "
		<<aColor.GetGreen()<<" "
		<<aColor.GetBlue()<<" "
		<<aColor.GetAlpha()<<"\n";

	fout.close();

	for(int i=0;i<nDaughter;i++)
	{ 
		GenerateGeometryMacro(aLogVol->GetDaughter(i),depth+1); 
	}

}
/////////////////////////////////////////////////////////////////////
void HRSDetectorConstruction::GetConfig()
{
	gConfig->GetParameter("HallX",mHallX);
	mHallX*=mm;
	gConfig->GetParameter("HallY",mHallY);
	mHallY*=mm;
	gConfig->GetParameter("HallZ",mHallZ);
	mHallZ*=mm;

	gConfig->GetParameter("FieldX",mFieldX);
	mFieldX*=mm;
	gConfig->GetParameter("FieldY",mFieldY);
	mFieldY*=mm;
	gConfig->GetParameter("FieldZ",mFieldZ);
	mFieldZ*=mm;

	mBMaterialType=2;  
	gConfig->GetParameter("BMaterialType",mBMaterialType);

	mSetupAirAxis=0;
	gConfig->GetParameter("SetupAirAxis",mSetupAirAxis);

	gConfig->GetParameter("SetupVirtualDetector",mSetupVirtualDetector);
	gConfig->GetParameter("VirtualDetectorWidth",mVirtualDetectorWidth);
	mVirtualDetectorWidth*=mm;
	gConfig->GetParameter("VirtualDetectorHeight",mVirtualDetectorHeight);
	mVirtualDetectorHeight*=mm;
	gConfig->GetParameter("VirtualDetectorThick",mVirtualDetectorThick);
	mVirtualDetectorThick*=mm;
	gConfig->GetParameter("VDRotYAngle",mVDRotYAngle);
	mVDRotYAngle*=deg;
	gConfig->GetParameter("VDRotXAngle",mVDRotXAngle);
	mVDRotXAngle*=deg;
	gConfig->GetParameter("Pivot2VDFace",mPivot2VDFace);
	mPivot2VDFace*=mm;

	mVDPhysVolName=gConfig->GetParameter("VDPhysVolName");

	/////////////////////////////////////////////////////////////////////////
	
	gConfig->GetParameter("SetupRadiator",mSetupRadiator);
	gConfig->GetParameter("SetupRadiatorVD",mSetupRadiatorVD);

	gConfig->GetParameter("RadiatorType",mRadiatorType);
	gConfig->GetParameter("RaditorThickness",mRaditorThickness);
	mRaditorThickness*=mm;
	gConfig->GetParameter("Radiator2Pivot",mRadiator2Pivot);
	mRadiator2Pivot*=mm;

	/////////////////////////////////////////////////////////////////////////

	mSetupCREXGeometry=0;
	gConfig->GetParameter("SetupCREXGeometry",mSetupCREXGeometry);

	mSetupRTPCGeometry=0;
	gConfig->GetParameter("SetupRTPCGeometry",mSetupRTPCGeometry);

	//G4cout<<"\n****Load detector config parameters done!***"<<G4endl;

	///////////////////////////////////////////////////////////////////////////
	//global variables

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

	gConfig->GetParameter("ScatChamberXOffset",mScatChamberXOffset);
	mScatChamberXOffset*=mm;
	gConfig->GetParameter("ScatChamberYOffset",mScatChamberYOffset);
	mScatChamberYOffset*=mm;
	gConfig->GetParameter("ScatChamberZOffset",mScatChamberZOffset);
	mScatChamberZOffset*=mm;

	gConfig->GetParameter("PivotXOffset",mPivotXOffset);
	mPivotXOffset*=mm;
	gConfig->GetParameter("PivotYOffset",mPivotYOffset);
	mPivotYOffset*=mm;
	gConfig->GetParameter("PivotZOffset",mPivotZOffset);
	mPivotZOffset*=mm;

	gConfig->GetParameter("TargetXOffset",mTargetXOffset);
	mTargetXOffset*=mm;
	gConfig->GetParameter("TargetYOffset",mTargetYOffset);
	mTargetYOffset*=mm;
	gConfig->GetParameter("TargetZOffset",mTargetZOffset);
	mTargetZOffset*=mm;


	return ;
}


/////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* HRSDetectorConstruction::Construct()
{
	G4String SDname;

	// All G4VSensitiveDetector will be managed (deleted) by SDManager
	//therefore no need to delete it at the end of this subroutine
	//G4VSensitiveDetector* virtualDetectorSD=new HRSStdSD(SDname="virtualDetector");
	G4VSensitiveDetector* virtualBoundarySD=new HRSStdSD(SDname="virtualBoundary");

	// sensitive detectors   
	G4SDManager* SDman = G4SDManager::GetSDMpointer();

	/////////////////////////////////////////////////////////////////////////////

	//set "user limits" for drawing smooth curve, only supported by Jixie's model
	//G4UserLimits* uHallStepLimits = new G4UserLimits(1000.*mm);//Jixie default
	G4UserLimits* uHallStepLimits = new G4UserLimits(.001 * mm);//Nickie default
	//double pBStepLimit=100;//Jixie default
	double pBStepLimit=0.001 * mm;//Nickie's edit
	gConfig->GetArgument("BStepLimit",pBStepLimit);
	pBStepLimit*=mm;
	G4cout << "B limits: " << pBStepLimit << G4endl;
	G4UserLimits* uBStepLimits = new G4UserLimits(pBStepLimit);

	///////////////////////////////////////////////////////////////////////////////
	// Magnetic field ----------------------------------------------------------
	mEMFieldSetup = new HRSEMFieldSetup();  //setup the field, 
	//G4FieldManager* dipoleFieldManager = mEMFieldSetup->GetFieldManagerFZB1();
	///////////////////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////////////
	//Shrink the hall size if bigbite is not built, this will mke it run fast
	//note that this is optimized for g2p only
	//These 2 from Detector_G2P.ini
	//int pSetupChicane=1,pSetupThirdArm=1;
	//if (!mSetupG2PGeometry) 
	//{
	//pSetupChicane=pSetupThirdArm=0;
	//}

	///////////////////////////////////////////////////////////////////////////

	// start to construct the geometries -----------------------------------------

	/////////////////////////
	// experimental hall (hall volume)
	/////////////////////////
	G4VSolid* hallSolid = new G4Box("hallBox",mHallX/2.0,mHallY/2.0,mHallZ/2.0);
	G4LogicalVolume* hallLogical = new G4LogicalVolume(hallSolid,
		mMaterialManager->vacuum,"hallLogical",0,0,uHallStepLimits);
	hallLogical->SetVisAttributes(HallVisAtt);
	mHallPhysVol = new G4PVPlacement(0,G4ThreeVector(),
		hallLogical,"hallPhys",0,0,0);


	/////////////////////////
	// air axes
	/////////////////////////
	//show the axis to help to understand the geometry in in the visualization
	//This part is also used to verify that the G4PVPlacement() will do the rotation at the origin 
	//then put the LogicalVolume at the given position 
	if(mSetupAirAxis)
	{
		G4RotationMatrix* pRotX270deg = new G4RotationMatrix();
		pRotX270deg->rotateX(270.*deg);
		G4RotationMatrix* pRotY270deg = new G4RotationMatrix();
		pRotY270deg->rotateY(270.*deg);

		G4VSolid* xAxisSolid = new G4Tubs("xAxisTubs",0.,5.0*mm,500.0*mm,0.,360.*deg);
		G4LogicalVolume* xAxisLogical = new G4LogicalVolume(xAxisSolid,
			mMaterialManager->vacuum,"xAxisLogical",0,0,0);
		new G4PVPlacement(pRotY270deg,G4ThreeVector(650*mm,0,0),
			xAxisLogical,"xAxisPhys",hallLogical,0,0);
		xAxisLogical->SetVisAttributes(G4VisAttributes(G4Colour::Red()));

		G4VSolid* yAxisSolid = new G4Tubs("yAxisTubs",0.,5.0*mm,500.0*mm,0.,360.*deg);
		G4LogicalVolume* yAxisLogical = new G4LogicalVolume(yAxisSolid,
			mMaterialManager->vacuum,"yAxisLogical",0,0,0);
		new G4PVPlacement(pRotX270deg,G4ThreeVector(0,650*mm,0),
			yAxisLogical,"yAxisPhys",hallLogical,0,0);
		yAxisLogical->SetVisAttributes(G4VisAttributes(G4Colour::Green()));

		G4VSolid* zAxisSolid = new G4Tubs("zAxisTubs",0.,5.0*mm,500.0*mm,0.,360.*deg);
		G4LogicalVolume* zAxisLogical = new G4LogicalVolume(zAxisSolid,
			mMaterialManager->vacuum,"zAxisLogical",0,0,0);
		new G4PVPlacement(0,G4ThreeVector(0,0,650*mm),
			zAxisLogical,"zAxisPhys",hallLogical,0,0);
		zAxisLogical->SetVisAttributes(G4VisAttributes(G4Colour::Blue()));
	}


	/////////////////////////
	// magnetic field region
	/////////////////////////
	G4Material* pBMaterial=0;
	if(mBMaterialType==0) pBMaterial=mMaterialManager->vacuum;
	else if(mBMaterialType==2) pBMaterial=mMaterialManager->heliumGas;
	else pBMaterial=mMaterialManager->air;

	G4VSolid* magneticSolid = new G4Box("magneticBox",mFieldX/2.0,mFieldY/2.0,mFieldZ/2.0);
	pBMaterial=mMaterialManager->vacuum;
	G4LogicalVolume* magneticLogical = new G4LogicalVolume(magneticSolid,
		pBMaterial,"magneticLogical",0,0,0);
	G4VisAttributes* MagneticVisAttTest = new G4VisAttributes(G4Colour(1., 0., 1.));
	magneticLogical->SetUserLimits(uBStepLimits);
	magneticLogical->SetVisAttributes(MagneticVisAtt);
	//magneticLogical->SetVisAttributes(MagneticVisAttTest);

	//the field region is the mother volumn of a lot of dauthers, it should not rotate
	new G4PVPlacement(0,G4ThreeVector(),magneticLogical,"magneticPhys",hallLogical,0,0);



	/////////////////////////
	//the radiator
	/////////////////////////
	if(mSetupRadiator) ConstructRadiator(magneticLogical);



	/////////////////////////
	//CREX geometry
	/////////////////////////	
	//scattering chamber, target, 
	//need to add sieve, septum, HRSVB, 
	if(mSetupCREXGeometry) 
	{
		CREXDetectorConstruction* theCREX = new CREXDetectorConstruction(magneticLogical); 
		theCREX->Construct();
		//update the parameters
		this->GetConfig(); 
	}

	/////////////////////////
	// Sieve slit, septum window and HRS
	/////////////////////////
	if(mSetupLHRS || mSetupRHRS)  this->ConstructHRS(magneticLogical);



	/////////////////////////
	// Cylinder Virtual Boundary
	/////////////////////////
	//-VB or -VirtualBoundary <SetupVirtualBoundary(0)> [VBRin(537)] [VBRout(540)] [VBHeight(2000)] [VBRotAxis(1)]
	//   [VBRotAngle(90)] [VBPosX(0)] [VBPosY(0)] [VBPosZ(0)]:
	//description: A switch to tell the program to place the cylinder virtual boundary. This 
	//solid is a cylinder made of air. One can specify the solid dimention, rotation and center position. 
	//Note that VBRotAxis can be 0(NO rotation) or 1 (X), 2(Y) and 3(Z) axis. The VBRotAngle is in deg. Only
	//one rotation is supported. The center position VBPosX,VBPosY and VBPosZ are all in mm and relative to 
	//the pivot, where is the center of the target in most settings.

	//one can use it for fast reason ...... 

	int pSetupVirtualBoundary=0;
	gConfig->GetArgument("SetupVirtualBoundary",pSetupVirtualBoundary); 
	if(pSetupVirtualBoundary)
	{	
		double pVBRin=537.0,pVBRout=538.0,pVBHeight=2000.0;
		gConfig->GetArgument("VBRin",pVBRin); 
		gConfig->GetArgument("VBRout",pVBRout); 
		gConfig->GetArgument("VBHeight",pVBHeight); 
		pVBRin*=mm,pVBRout*=mm,pVBHeight*=mm;

		G4VSolid* virtualBoundary0Solid = new G4Tubs("virtualBoundary0Tubs",
			pVBRin,(pVBRin+pVBRout)/2,pVBHeight/2.0,-240*deg,300*deg);
		G4VSolid* virtualBoundary1Solid = new G4Tubs("virtualBoundary1Tubs",
			(pVBRin+pVBRout)/2,pVBRout,pVBHeight/2.0,-240*deg,300*deg);
		G4LogicalVolume* virtualBoundary0Logical = new G4LogicalVolume(virtualBoundary0Solid,
			mMaterialManager->air,"virtualBoundary0Logical",0,0,0);
		G4LogicalVolume* virtualBoundary1Logical = new G4LogicalVolume(virtualBoundary1Solid,
			mMaterialManager->air,"virtualBoundary1Logical",0,0,0);

		virtualBoundary0Logical->SetVisAttributes(LightYellowVisAtt);  
		virtualBoundary1Logical->SetVisAttributes(LightYellowVisAtt); 
		SDman->AddNewDetector(virtualBoundarySD);
		virtualBoundary0Logical->SetSensitiveDetector(virtualBoundarySD);

		double pVBPosX=0,pVBPosY=0,pVBPosZ=0;
		gConfig->GetArgument("VBPosX",pVBPosX); 
		gConfig->GetArgument("VBPosY",pVBPosY); 
		gConfig->GetArgument("VBPosZ",pVBPosZ); 
		pVBPosX=pVBPosX*mm+mPivotXOffset;
		pVBPosY=pVBPosX*mm+mPivotYOffset;
		pVBPosZ=pVBPosZ*mm+mPivotZOffset;

		G4RotationMatrix *pRotVB=new G4RotationMatrix();
		int pVBRotAxis=0;
		double pVBRotAngle=0; 
		gConfig->GetArgument("VBRotAxis",pVBRotAxis); 
		gConfig->GetArgument("VBRotAngle",pVBRotAngle); pVBRotAngle*=deg;
		if(pVBRotAxis && fabs(pVBRotAngle)>1.0E-05)
		{
			if(pVBRotAxis==1) pRotVB->rotateX(pVBRotAngle); 
			else if(pVBRotAxis==2) pRotVB->rotateY(pVBRotAngle); 
			else if(pVBRotAxis==3) pRotVB->rotateZ(pVBRotAngle); 
		}

		new G4PVPlacement(pRotVB,G4ThreeVector(pVBPosX,pVBPosY,pVBPosZ),			
			virtualBoundary0Logical,"virtualDetectorPhys",magneticLogical,0,0);
		new G4PVPlacement(pRotVB,G4ThreeVector(pVBPosX,pVBPosY,pVBPosZ),			
			virtualBoundary1Logical,"virtualBoundaryPhys",magneticLogical,0,0);
	}


#ifdef G4DEBUG_GEOMETRY
	if(G4DEBUG_GEOMETRY>=0)
	{
		DumpGeometricalTree(mHallPhysVol);
	}
#endif
	//create the geometry macro to manually modify the geometry in visulization
	GenerateGeometryMacro(mHallPhysVol,0);
	return mHallPhysVol;
}

/////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* HRSDetectorConstruction::ConstructHRS(G4LogicalVolume* motherLogical)
{
	int mSnakeModel;
	gConfig->GetArgument("SnakeModel",mSnakeModel);

	int IS_NIM = 0;//if it's 0, then we go by SNAKE, if it is 1, we go by NIM.
	//As it turns out, I think we want to go by SNAKE, so that is why it is zero right now.
	G4VPhysicalVolume* theHRSPhys=0;

	G4RotationMatrix *pRotX90deg=new G4RotationMatrix();
	pRotX90deg->rotateX(-90*deg); 
	G4RotationMatrix *pRotX45deg=new G4RotationMatrix();
	pRotX45deg->rotateX(-45*deg); 
	G4RotationMatrix *pRotX30deg=new G4RotationMatrix();
	pRotX30deg->rotateX(60*deg); 
	G4RotationMatrix *pRotX105deg=new G4RotationMatrix();
	pRotX105deg->rotateX(165*deg); 
	G4RotationMatrix *pRotXLHRSdeg=new G4RotationMatrix();
	pRotXLHRSdeg->rotateY(mLHRSAngle); 
	G4RotationMatrix *pRotXRHRSdeg=new G4RotationMatrix();
	pRotXRHRSdeg->rotateY(mRHRSAngle); 

	G4SDManager* SDman = G4SDManager::GetSDMpointer();
	G4String SDname;
	G4VSensitiveDetector* Q1WindowSD=new HRSStdSD(SDname="Q1Window");
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
	//                         7.0m above beam line 
	//                       -----------------------------| 20.76m from pivot,
	//                      /                             |
	//                     /                        Q3    |
	//                    /                               |
	//                   /                       E        |
	//                  /                     L           |
	//          --------                   O              |
	//         /                        P                 |
	//    -----                     I                     |
	//----|----- Q1 --- Q2 ---- D     --------------------|------ beam line -----
	//    |                                               |
	//    ------------------------------------------------|
	//    1.46m from povot, 3.05 m below beam line

	*/

	//double pHRSContainerRin=1.46*m,pHRSContainerRout=20.76*m;//JixieMode
	double pHRSContainerRin=1.46*m,pHRSContainerRout=25*m;//NickieMode it just has to be big enough to fit my VB at FP
	//double pBeamLine2Ground;
	//pBeamLine2Ground=-3.05*m;
	//build the container with polycone

	const int kNPlane_HRSContainer=7;
	double rInner_HRSContainer[] = {pHRSContainerRin,pHRSContainerRin,2.5*m,
		3.7*m,9.0*m,pHRSContainerRout-3.0*m,pHRSContainerRout};
	double rOuter_HRSContainer[] = {pHRSContainerRout,pHRSContainerRout,pHRSContainerRout,
		pHRSContainerRout,pHRSContainerRout,pHRSContainerRout,pHRSContainerRout};
	//double zPlane_HRSContainer[] = {-2.0*m,1.0*m,1.0*m,
	//2.0*m,2.0*m,7.0*m,7.0*m};
	double zPlane_HRSContainer[] = {-2.0*m,1.0*m,1.0*m,
					2.0*m,2.0*m,11.0*m,11.0*m};
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

	G4VisAttributes* MagFieldVisAtt = new G4VisAttributes(G4Colour(1., 0., 1.));	

	LHRSContainerLogical->SetVisAttributes(HallVisAtt); 
	RHRSContainerLogical->SetVisAttributes(HallVisAtt);
	//RHRSContainerLogical->SetVisAttributes(MagFieldVisAtt); 

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


	//////////////////////////
	//Q1 Entrance Window //
	//////////////////////////
	bool pSetupQ1EntranceWindow=true;
	//G2P and CREX have their own VB defined at a different location
	if( mSetupCREXGeometry)  pSetupQ1EntranceWindow=false;
	if(pSetupQ1EntranceWindow) 
	{

		//this part is trying to place a virtual boundary at the Q1 entrance

		//Place both left and right VB for HRS, which is pHRSContainerRin+4*mm away from the 
		//hall center(1.462m). This aperture is a round disk of 29.8 cm diameter and 2 mm thick
		//The real Q1 vacumn entrance to hall center is 1.312m, 

		double pHRSQ1WinThick = 2*mm;
		G4VSolid* HRSQ1WinSolid = new G4Tubs("HRSQ1WinTub",0.0,14.9*cm,
			pHRSQ1WinThick/2,0.0,360.0*deg);
		G4LogicalVolume* HRSQ1WinLogical = new G4LogicalVolume(HRSQ1WinSolid,
			mMaterialManager->mylar,"HRSQ1WinLogical",0,0,0);
		SDman->AddNewDetector(Q1WindowSD);
		HRSQ1WinLogical->SetSensitiveDetector(Q1WindowSD);
		HRSQ1WinLogical->SetVisAttributes(LightYellowVisAtt); 

		//since the container has been rotated by 90 deg about x axis,
		//y'=z  z'=-y ==> I have to put this offset as -y
		double pHallCenter2Q1Win=pHRSContainerRin+4*mm;  //place it at the first 1.464 m
		gConfig->SetParameter("Pivot2LHRSVBFace",pHallCenter2Q1Win-mPivotZOffset*cos(mLHRSAngle));
		gConfig->SetParameter("Pivot2RHRSVBFace",pHallCenter2Q1Win-mPivotZOffset*cos(mRHRSAngle)); 

		if(mSetupLHRS)
		{
			new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pHallCenter2Q1Win,0),
				HRSQ1WinLogical,"virtualBoundaryPhys_LHRSQ1Win",LHRSContainerLogical,0,0,0);
		}
		if(mSetupRHRS)
		{
			new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pHallCenter2Q1Win,0),
				HRSQ1WinLogical,"virtualBoundaryPhys_RHRSQ1Win",RHRSContainerLogical,0,0,0);
		}
	}
	//BELOW ARE SNAKE VALUES, NOT NIM VALUES
	double pTarget =             0.0  * cm;
	double pQ1en   = pTarget + 160.0  * cm;
	double pQ1ex   = pQ1en   +  94.13 * cm;
	double pQ2en   = pQ1ex   + 115.58 * cm;
	double pQ2ex   = pQ2en   + 182.66 * cm;
	//ABOVE ARE SNAKE VALUES, NOT NIM VALUES
	
	/////////////////////////
	// HRS Q1              // Nickie is also adding the Q1 collimator here
	/////////////////////////
	double pHallCenter2LQ1Face;//=1.69*m;//NIM
	double pHallCenter2RQ1Face;//=1.69*m;//NIM
	//pHallCenter2LQ1Face = ( IS_NIM == 1 ) ? 169.0 * cm : pQ1en;
	//pHallCenter2RQ1Face = ( IS_NIM == 1 ) ? 169.0 * cm : pQ1en;
	pHallCenter2LQ1Face = pQ1en;
	pHallCenter2RQ1Face = pQ1en;
	//double pHallCenter2LQ1FaceMag=1.600*m;//SNAKE
	//double pHallCenter2RQ1FaceMag=1.600*m;//SNAKE
	//double fringe_extension = 1.25;
	//double fringe_extension = 1.2;
	double fringe_extension = 1.00;

	//flag to trigger sos quad instead of old quad for the proper runs
	//K. Allada and B. Schmookler:
	//Magnetic Length = 70 cm
	//Radius to Pole Tip = 12.827  cm
	
	int    sos   = 1;
	if( mSnakeModel >= 52){
	  sos = 1;
	}
	double q1shift = sos ? 0.0 * m : 0.0 * m;
	double pQ1Rin  = sos ? 12.827 * cm : 15.0*cm;
	double pQ1Rout = sos ? 35.0   * cm : 35.0*cm;//for now, keep the outer same for either
	double pQ1Length;//=80*cm;
	//pQ1Length = ( IS_NIM == 1 ) ? 80.0 * cm : pQ1ex - pQ1en;
	pQ1Length = sos ?  70. * cm : pQ1ex - pQ1en;
	double pQ1PosAct = pQ1ex - pQ1en;
	double pQ1LengthMag=pQ1Length * fringe_extension;
	//double pQ1LengthMag=94.1*cm;
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

	if(mSetupLHRS>=2){
	  //put it in the container, which also center at the hall center
	  //therefore only the z_at_lab position need to be considered
	  double pLQ1Pos_Z=(pHallCenter2LQ1Face+pQ1PosAct/2.0) + q1shift;
	  //double pLQ1Pos_Z=pQ1PosAct + q1shift;
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pLQ1Pos_Z,0),
			    LQ1Logical,"LQ1Phys",LHRSContainerLogical,0,0,0);
	}
	if(mSetupRHRS>=2){
	  double pRQ1Pos_Z=(pHallCenter2RQ1Face+pQ1PosAct/2.0) + q1shift;
	  //double pRQ1Pos_Z=pQ1PosAct + q1shift;
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pRQ1Pos_Z,0),
			    RQ1Logical,"RQ1Phys",RHRSContainerLogical,0,0,0);
	}

	//vac
	double pQ1vacRin  = pQ1Rin;
	double pQ1vacRout = pQ1Rin + 5. * cm;
	double pQ1vacLength = 610.4 * mm; 
	double pQ1vacCenterZ = pHallCenter2LQ1Face + pQ1PosAct + pQ1vacLength / 2.;
	double pQ1vacCenterY = 0;
	G4VSolid* Q1vacSolid = new G4Tubs("Q1vacTub",pQ1vacRin,pQ1vacRout,pQ1vacLength/2.0,0.0,360.0*deg);

	G4LogicalVolume* LQ1vacLogical = new G4LogicalVolume(Q1vacSolid,
		mMaterialManager->siliconsteel,"LQ1vacLogical",0,0,0);
	G4LogicalVolume* RQ1vacLogical = new G4LogicalVolume(Q1vacSolid,
		mMaterialManager->siliconsteel,"RQ1vacLogical",0,0,0);

	LQ1vacLogical->SetVisAttributes(IronVisAtt); 
	RQ1vacLogical->SetVisAttributes(IronVisAtt); 

	G4RotationMatrix *pRotQ1vacInContainer=new G4RotationMatrix();
	pRotQ1vacInContainer->rotateX(90*deg); 
	
	if(mSetupLHRS>=4)
	{
		//put it in the container, which also center at the hall center
		new G4PVPlacement(pRotQ1vacInContainer,G4ThreeVector(0, -pQ1vacCenterZ,pQ1vacCenterY),
			LQ1vacLogical,"LQ1vacPhys",LHRSContainerLogical,0,0,0);
	}
	if(mSetupRHRS>=4)
	{
		new G4PVPlacement(pRotQ1vacInContainer,G4ThreeVector(0, -pQ1vacCenterZ,pQ1vacCenterY),
			RQ1vacLogical,"RQ1vacPhys",RHRSContainerLogical,0,0,0);
	}
	

	/////////////////////////
	// HRS Q2 
	/////////////////////////
	double pHallCenter2LQ2Face;//=3.74*m;//NIM
	double pHallCenter2RQ2Face;//=3.74*m;//NIM
	//pHallCenter2LQ2Face = ( IS_NIM == 1 ) ? 374.0 * cm : pQ2en;
	//pHallCenter2RQ2Face = ( IS_NIM == 1 ) ? 374.0 * cm : pQ2en;
	pHallCenter2LQ2Face = pQ2en;
	pHallCenter2RQ2Face = pQ2en;
	//double pHallCenter2LQ2FaceMag=3.696*m;//SNAKE
	//double pHallCenter2RQ2FaceMag=3.696*m;//SNAKE
	double pQ2Rin=30.0*cm;
	double pQ2Rout=75.0*cm;
	double pQ2Length;//=180*cm;//NIM
	//pQ2Length = ( IS_NIM == 1 ) ? 180.0 * cm : pQ2ex - pQ2en;
	pQ2Length = pQ2ex - pQ2en;
	double pQ2LengthMag=pQ2Length * fringe_extension;
	//double pQ2LengthMag=182.6*cm;//SNAKE

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

	//vac
	
	double pQ2vacRin  = pQ2Rin;
	double pQ2vacRout = pQ2Rin + 5. * cm;
	double pQ2vacLength = ( 1166.06 - 675. ) * mm;
	double pQ2vacCenterZ = pHallCenter2LQ2Face - pQ2vacLength / 2.;
	double pQ2vacCenterY = 0;
	G4VSolid* Q2vacSolid = new G4Tubs("Q2vacTub",pQ2vacRin,pQ2vacRout,pQ2vacLength/2.0,0.0,360.0*deg);

	G4LogicalVolume* LQ2vacLogical = new G4LogicalVolume(Q2vacSolid,
		mMaterialManager->siliconsteel,"LQ2vacLogical",0,0,0);
	G4LogicalVolume* RQ2vacLogical = new G4LogicalVolume(Q2vacSolid,
		mMaterialManager->siliconsteel,"RQ2vacLogical",0,0,0);

	LQ2vacLogical->SetVisAttributes(IronVisAtt); 
	RQ2vacLogical->SetVisAttributes(IronVisAtt); 

	G4RotationMatrix *pRotQ2vacInContainer=new G4RotationMatrix();
	pRotQ2vacInContainer->rotateX(90*deg); 
	if(mSetupLHRS>=4)
	{
		//put it in the container, which also center at the hall center
		new G4PVPlacement(pRotQ2vacInContainer,G4ThreeVector(0, -pQ2vacCenterZ,pQ2vacCenterY),
			LQ2vacLogical,"LQ2vacPhys",LHRSContainerLogical,0,0,0);
	}
	if(mSetupRHRS>=4)
	{
		new G4PVPlacement(pRotQ2vacInContainer,G4ThreeVector(0, -pQ2vacCenterZ,pQ2vacCenterY),
			RQ2vacLogical,"RQ2vacPhys",RHRSContainerLogical,0,0,0);
	}
	

	//vac2
	double pQ2vac2Rin     = pQ2Rin;
	double pQ2vac2Rout    = pQ2Rin + 1. * mm;
	double pQ2vac2Length  = 561.7 * mm;
	double pQ2vac2CenterZ = pHallCenter2RQ2Face + pQ2Length + pQ2vac2Length / 2.;
	double pQ2vac2CenterY = 0;
	G4VSolid* Q2vac2Solid = new G4Tubs("Q2vac2Tub",pQ2vac2Rin,pQ2vac2Rout,pQ2vac2Length/2.0,0.0,360.0*deg);

	G4LogicalVolume* LQ2vac2Logical = new G4LogicalVolume(Q2vac2Solid,
		mMaterialManager->siliconsteel,"LQ2vac2Logical",0,0,0);
	G4LogicalVolume* RQ2vac2Logical = new G4LogicalVolume(Q2vac2Solid,
		mMaterialManager->siliconsteel,"RQ2vac2Logical",0,0,0);

	LQ2vac2Logical->SetVisAttributes(IronVisAtt); 
	RQ2vac2Logical->SetVisAttributes(IronVisAtt); 

	G4RotationMatrix *pRotQ2vac2InContainer=new G4RotationMatrix();
	pRotQ2vac2InContainer->rotateX(90*deg); 
	if(mSetupLHRS>=4)
	{
		//put it in the container, which also center at the hall center
		new G4PVPlacement(pRotQ2vac2InContainer,G4ThreeVector(0, -pQ2vac2CenterZ,pQ2vac2CenterY),
			LQ2vac2Logical,"LQ2vac2Phys",LHRSContainerLogical,0,0,0);
	}
	if(mSetupRHRS>=4)
	{
		new G4PVPlacement(pRotQ2vac2InContainer,G4ThreeVector(0, -pQ2vac2CenterZ,pQ2vac2CenterY),
			RQ2vac2Logical,"RQ2vac2Phys",RHRSContainerLogical,0,0,0);
	}

	/////////////////////////
	// HRS Dipole 
	/////////////////////////
	//The dipole is built as a disc subtraced subtract another disc to get the side 
	//then subtract by the tunnel disc
	double pDipoleBendAngle=45.*deg, pDipoleFaceAngle=30.*deg;
	double pDipoleR=8.4*m;
	double pDipoleRprime=pDipoleR*sin(pDipoleBendAngle/2.)/
		sin(180.*deg-pDipoleBendAngle/2.-pDipoleFaceAngle);
	//double pDipoleRprime=4.0518*m;

	double pDipoleRCenterY=pDipoleR;
	double pDipoleRCenterZ;//=9.96*m;
	//pDipoleRCenterZ = ( IS_NIM == 1 ) ? 9.96 * m : 9.961 * m;//snake
	pDipoleRCenterZ = 9.961 * m;//snake
	//pDipoleRCenterZ = ( IS_NIM == 1 ) ? 9.96 * m : 9.9547 * m;//seamus
	G4cout << pDipoleRCenterZ << G4endl;
	//double pDipoleRprimeCenterY=pDipoleRprime*cos(pDipoleFaceAngle);			// =2.865*m;
	//double pDipoleRprimeCenterZ=9.96*m+pDipoleRprime*sin(pDipoleFaceAngle);	//12.825*m;

	//the center of Rprime relative to R 
	double pRprime2R_Y=pDipoleRprime*cos(pDipoleFaceAngle)-pDipoleR;	// =-5.535*m;
	double pRprime2R_Z=pDipoleRprime*sin(pDipoleFaceAngle);				// =1.865*m;

	//the original disc
	G4VSolid* DipoleWholeTub = new G4Tubs("DipoleWholeTub",
		pDipoleR-0.8*m,pDipoleR+0.8*m,0.4*m,
		172.*deg,pDipoleBendAngle+16.*deg);
	//the disc to be subtracted, musu be thicker 
	G4VSolid* DipolePrimeTub = new G4Tubs("DipolePrimeTub",
		0,pDipoleR,0.5*m,
		180.*deg+pDipoleFaceAngle+pDipoleBendAngle,360.*deg-pDipoleFaceAngle*2.-pDipoleBendAngle);
	//subtract the small tube to form the shape of the sides
	G4SubtractionSolid* DipoleWithSides = new G4SubtractionSolid("DipoleWithSides",
		DipoleWholeTub,DipolePrimeTub,
		0,G4ThreeVector(pRprime2R_Y,-pRprime2R_Z,0));

	//the tunnel disc, I use a rectangle shape here
	//G4VSolid* DipoleTunnelTub = new G4Tubs("DipoleTunnelTub",
	//	pDipoleR-0.4*m,pDipoleR+0.4*m,0.125*m,
	//	170*deg,pDipoleBendAngle+20*deg);
	//The shape of the dipole tunnel is a trapzoid, x=+/-0.4; y=+/-(0.125*(1-(1.25*x/8.40)) Jixie
	//To build this shape, I have to use polycone Jixie
	// -5.22008 < x < -4.98099 J.L.R.
	//  -(-0.192436784*x -0.192436784) < y < -0.192436784*x -0.192436784 J.L.R.
	double dy=0.125*(1.25*0.40/8.40)*m;
	const int kNPlane_DipoleTunnel=4;
	double rInner_DipoleTunnel[] = {pDipoleR-0.4*m,pDipoleR-0.4*m,pDipoleR-0.4*m,pDipoleR-0.4*m};
	double rOuter_DipoleTunnel[] = {pDipoleR-0.4*m,pDipoleR+0.4*m,pDipoleR+0.4*m,pDipoleR-0.4*m};
	double zPlane_DipoleTunnel[] = {-0.125*m-dy,-0.125*m+dy,0.125*m-dy,0.125*m+dy};
	G4Polycone* DipoleTunnelCone = new G4Polycone("DipoleTunnelCone",
		    170.0*deg,pDipoleBendAngle+20.*deg,
		    kNPlane_DipoleTunnel,zPlane_DipoleTunnel,rInner_DipoleTunnel,rOuter_DipoleTunnel);//Jixie's original
	
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
	G4RotationMatrix *pRotDipoleFringe1InContainer=new G4RotationMatrix();
	G4RotationMatrix *pRotDipoleFringe2InContainer=new G4RotationMatrix();
	
	pRotDipoleInContainer->rotateY(90*deg);
	//pRotDipoleFringeInContainer->rotateY(-90*deg);
	pRotDipoleFringe1InContainer->rotateY(180*deg);
	pRotDipoleFringe1InContainer->rotateZ(90*deg);
	pRotDipoleFringe1InContainer->rotateX(90*deg); 
	pRotDipoleFringe2InContainer->rotateX(90*deg);
	pRotDipoleFringe2InContainer->rotateY(90*deg + 45.*deg); 

	G4RotationMatrix *pRot_den=new G4RotationMatrix();
	G4RotationMatrix *pRot_dex=new G4RotationMatrix();	
	pRot_den->rotateY(90*deg);
	pRot_den->rotateX(-60*deg);
	pRot_dex->rotateY(90*deg);
	pRot_dex->rotateX(-90*deg + 105 * deg);

	//pRotDipoleFringe2InContainer->rotateX(90*deg); 
	//if(0)
	if(mSetupLHRS>=4)
	{
	  //put it in the container, which also center at the hall center
	  new G4PVPlacement(pRotDipoleInContainer,
			    G4ThreeVector(0,-pDipoleRCenterZ,pDipoleRCenterY),
			    LDipoleLogical,"LDipolePhys",LHRSContainerLogical,0,0,0);
	}
	//if(0)
	if(mSetupRHRS>=4)
	{
	  new G4PVPlacement(pRotDipoleInContainer,
			    G4ThreeVector(0,-pDipoleRCenterZ,pDipoleRCenterY),
			    RDipoleLogical,"RDipolePhys",RHRSContainerLogical,0,0,0);
	}

	


	/////////////////////////
	// HRS Q3 
	/////////////////////////
	double pQ3Rin=30.0*cm;
	double pQ3Rout=75.0*cm;
	double pQ3Length;//=180*cm;//NIM
	//pQ3Length = ( IS_NIM == 1 ) ? 180.0 * cm : 182.68 * cm; 
	pQ3Length = 182.68 * cm; 
	double pQ3CenterY;//=pDipoleR*(1-cos(pDipoleBendAngle))+2.4*m*sin(pDipoleBendAngle);//NIM
	double pQ3CenterZ;//=9.96*m+pDipoleR*sin(pDipoleBendAngle)+2.4*m*cos(pDipoleBendAngle);//NIM
	//pQ3CenterY = ( IS_NIM == 1 ) ? pDipoleR * ( 1 - cos( pDipoleBendAngle ) ) + 2.4 * m * sin( pDipoleBendAngle ) : 3.5853101 * m  + pQ3Length / sqrt(2.) / 2.;
	//pQ3CenterZ = ( IS_NIM == 1 ) ? 9.96 * m + pDipoleR * sin( pDipoleBendAngle ) + 2.4 * m * cos( pDipoleBendAngle ) : 17.0257042 * m  + pQ3Length / sqrt(2) / 2;
	pQ3CenterY =  3.5853101 * m  + pQ3Length / sqrt(2.) / 2.;
	pQ3CenterZ = 17.0257042 * m  + pQ3Length / sqrt(2.) / 2;

	double pQ3LengthMag=pQ3Length * fringe_extension;
	//double pQ3LengthMag=182.68*cm;//SNAKE

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


	//vacuum pipe
	double pQ3vacRin  = 30.0*cm;
	double pQ3vacRout = 35.0*cm;
	double pQ3vacLength = 575. * mm; 
	double pQ3vacCenterY = pQ3CenterY + pQ3Length / 2. / sqrt( 2 ) + pQ3vacLength / 2. / sqrt( 2 );
	double pQ3vacCenterZ = pQ3CenterZ + pQ3Length / 2. / sqrt( 2 ) + pQ3vacLength / 2. / sqrt( 2 );
	G4VSolid* Q3vacSolid = new G4Tubs("Q3vacTub",pQ3vacRin,pQ3vacRout,pQ3vacLength/2.0,0.0,360.0*deg);

	G4LogicalVolume* LQ3vacLogical = new G4LogicalVolume(Q3vacSolid,
		mMaterialManager->siliconsteel,"LQ3vacLogical",0,0,0);
	G4LogicalVolume* RQ3vacLogical = new G4LogicalVolume(Q3vacSolid,
		mMaterialManager->siliconsteel,"RQ3vacLogical",0,0,0);

	LQ3vacLogical->SetVisAttributes(IronVisAtt); 
	RQ3vacLogical->SetVisAttributes(IronVisAtt); 

	G4RotationMatrix *pRotQ3vacInContainer=new G4RotationMatrix();
	pRotQ3vacInContainer->rotateX(-45*deg); 
	if(mSetupLHRS>=4)
	{
		//put it in the container, which also center at the hall center
		new G4PVPlacement(pRotQ3vacInContainer,G4ThreeVector(0,-pQ3vacCenterZ,pQ3vacCenterY),
			LQ3vacLogical,"LQ3vacPhys",LHRSContainerLogical,0,0,0);
	}
	if(mSetupRHRS>=4)
	{
		new G4PVPlacement(pRotQ3vacInContainer,G4ThreeVector(0,-pQ3vacCenterZ,pQ3vacCenterY),
			RQ3vacLogical,"RQ3vacPhys",RHRSContainerLogical,0,0,0);
	}
	//vacuum pipe
	double pQ3vac2CenterY = pQ3CenterY - pQ3Length / 2. / sqrt( 2 ) - pQ3vacLength / 2. / sqrt( 2 );
	double pQ3vac2CenterZ = pQ3CenterZ - pQ3Length / 2. / sqrt( 2 ) - pQ3vacLength / 2. / sqrt( 2 );
	G4VSolid* Q3vac2Solid = new G4Tubs("Q3vac2Tub",pQ3vacRin,pQ3vacRout,pQ3vacLength/2.0,0.0,360.0*deg);

	G4LogicalVolume* LQ3vac2Logical = new G4LogicalVolume(Q3vac2Solid,
		mMaterialManager->siliconsteel,"LQ3vac2Logical",0,0,0);
	G4LogicalVolume* RQ3vac2Logical = new G4LogicalVolume(Q3vac2Solid,
		mMaterialManager->siliconsteel,"RQ3vac2Logical",0,0,0);

	LQ3vac2Logical->SetVisAttributes(IronVisAtt); 
	RQ3vac2Logical->SetVisAttributes(IronVisAtt); 

	G4RotationMatrix *pRotQ3vac2InContainer=new G4RotationMatrix();
	pRotQ3vac2InContainer->rotateX(-45*deg); 
	if(mSetupLHRS>=4)
	{
		//put it in the container, which also center at the hall center
		new G4PVPlacement(pRotQ3vac2InContainer,G4ThreeVector(0,-pQ3vac2CenterZ,pQ3vac2CenterY),
			LQ3vac2Logical,"LQ3vac2Phys",LHRSContainerLogical,0,0,0);
	}
	if(mSetupRHRS>=4)
	{
		new G4PVPlacement(pRotQ3vac2InContainer,G4ThreeVector(0,-pQ3vac2CenterZ,pQ3vac2CenterY),
			RQ3vac2Logical,"RQ3vac2Phys",RHRSContainerLogical,0,0,0);
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
	double zPlane_Q1Q2Tunnel[] = {pHRSContainerRin+6*mm,pHallCenter2RQ1Face+pQ1PosAct+10.0*cm,
		pHallCenter2RQ1Face+pQ1PosAct+30*cm,9.7*m};
	G4Polycone* Q1Q2TunnelSolid = new G4Polycone("Q1Q2TunnelPolycone",0.0,360.0*deg,
		kNPlane_Q1Q2Tunnel,zPlane_Q1Q2Tunnel,rInner_Q1Q2Tunnel,rOuter_Q1Q2Tunnel);

	G4LogicalVolume* LQ1Q2TunnelLogical = new G4LogicalVolume(Q1Q2TunnelSolid,
		mMaterialManager->vacuum,"LQ1Q2TunnelLogical",0,0,0);
	G4LogicalVolume* RQ1Q2TunnelLogical = new G4LogicalVolume(Q1Q2TunnelSolid,
		mMaterialManager->vacuum,"RQ1Q2TunnelLogical",0,0,0);


	LQ1Q2TunnelLogical->SetVisAttributes(LightYellowVisAtt); 
	RQ1Q2TunnelLogical->SetVisAttributes(LightYellowVisAtt); 
	/*
	if(mSetupLHRS>=3){//This is the old Jixie way
	  //put it in the container, which also center at the hall center
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,0,0),
			    LQ1Q2TunnelLogical,"LQ1Q2TunnelPhys",LHRSContainerLogical,0,0,0);
	}
	if(mSetupRHRS>=3){//This is the old Jixie way
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,0,0),
			    RQ1Q2TunnelLogical,"RQ1Q2TunnelPhys",RHRSContainerLogical,0,0,0);
	}
	*/
 	//This is a copy of Q1, but of smaller radius.//This is the new Nickie way
	//G4VSolid* Q1SolidMag = new G4Tubs("Q1TubMag",0,pQ1Rin,pQ1Length/2.0,0.0,360.0*deg);//NIM
	G4VSolid* Q1SolidMag = new G4Tubs("Q1TubMag",0,pQ1Rin,pQ1LengthMag/2.0,0.0,360.0*deg);//for fringe fields, note longer length

	//build 2 copy since there are different fields involved in it
	G4LogicalVolume* LQ1MagLogical = new G4LogicalVolume(Q1SolidMag,
		mMaterialManager->vacuum,"LQ1MagLogical",0,0,0);
	G4LogicalVolume* RQ1MagLogical = new G4LogicalVolume(Q1SolidMag,
		mMaterialManager->vacuum,"RQ1MagLogical",0,0,0);

	//Nickie's Q1 field
	G4FieldManager* LQ1FieldManager = mEMFieldSetup->GetFieldManagerFZBL1();
	G4FieldManager* RQ1FieldManager = mEMFieldSetup->GetFieldManagerFZBR1();
	G4bool allLocal = true;
	//G4bool allLocal = false;
	LQ1MagLogical->SetFieldManager(LQ1FieldManager,allLocal);
	RQ1MagLogical->SetFieldManager(RQ1FieldManager,allLocal);

	LQ1MagLogical->SetVisAttributes(LightYellowVisAtt); 
	RQ1MagLogical->SetVisAttributes(LightYellowVisAtt); 


	//Nickie's old way
	if(mSetupLHRS>=2){
	  //put it in the container, which also center at the hall center
	  //therefore only the z_at_lab position need to be considered
	  double pLQ1Pos_Z=(pHallCenter2LQ1Face+pQ1PosAct/2.0 + q1shift);//NIM
	  //double pLQ1Pos_Z=(pQ1PosAct + q1shift);//NIM
	  //double pLQ1Pos_Z=(pHallCenter2LQ1FaceMag+pQ1LengthMag/2.0);//SNAKE
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pLQ1Pos_Z,0),
			    LQ1MagLogical,"LQ1MagPhys",LHRSContainerLogical,0,0,0);
	}if(mSetupRHRS>=2){
	  double pRQ1Pos_Z=(pHallCenter2RQ1Face+pQ1PosAct/2.0 + q1shift);//NIM
	  //double pRQ1Pos_Z=(pQ1PosAct + q1shift);//NIM
	  //double pRQ1Pos_Z=(pHallCenter2RQ1FaceMag+pQ1LengthMag/2.0);//SNAKE
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pRQ1Pos_Z,0),
			    RQ1MagLogical,"RQ1MagPhys",RHRSContainerLogical,0,0,0);
	}
	
	//G4VSolid* Q2SolidMag = new G4Tubs("Q2TubMag",0,pQ2Rin,pQ2Length/2.0,0.0,360.0*deg);//NIM
	G4VSolid* Q2SolidMag = new G4Tubs("Q2TubMag",0,pQ2Rin,pQ2LengthMag/2.0,0.0,360.0*deg);//SNAKE

	G4LogicalVolume* LQ2MagLogical = new G4LogicalVolume(Q2SolidMag,
		mMaterialManager->vacuum,"LQ2MagLogical",0,0,0);
	G4LogicalVolume* RQ2MagLogical = new G4LogicalVolume(Q2SolidMag,
		mMaterialManager->vacuum,"RQ2MagLogical",0,0,0);

	//Nickie's Q2 field
	G4FieldManager* LQ2FieldManager = mEMFieldSetup->GetFieldManagerFZBL2();
	G4FieldManager* RQ2FieldManager = mEMFieldSetup->GetFieldManagerFZBR2();
	
	LQ2MagLogical->SetFieldManager(LQ2FieldManager,allLocal);
	RQ2MagLogical->SetFieldManager(RQ2FieldManager,allLocal);

	LQ2MagLogical->SetVisAttributes(LightYellowVisAtt); 
	RQ2MagLogical->SetVisAttributes(LightYellowVisAtt); 

	 //Nickie's old way
	if(mSetupLHRS>=3){
	  //put it in the container, which also center at the hall center
	  double pLQ2Pos_Z=(pHallCenter2LQ2Face+pQ2Length/2.0);//NIM
	  //double pLQ2Pos_Z=(pHallCenter2LQ2FaceMag+pQ2LengthMag/2.0);//SNAKE
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pLQ2Pos_Z,0),
			    LQ2MagLogical,"LQ2MagPhys",LHRSContainerLogical,0,0,0);
	}if(mSetupRHRS>=3){
	  double pRQ2Pos_Z=(pHallCenter2RQ2Face+pQ2Length/2.0);//NIM
	  //double pRQ2Pos_Z=(pHallCenter2RQ2FaceMag+pQ2LengthMag/2.0);//SNAKE
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pRQ2Pos_Z,0),
			    RQ2MagLogical,"RQ2MagPhys",RHRSContainerLogical,0,0,0);
	}
	

	//////////////////////////
	//Dipole tunnel
	//////////////////////////
	//G4VSolid* DipoleTunnelSolid = new G4Tubs("DipoleTunnelSolid",pDipoleR-0.4*m+0.1*mm,
	//	pDipoleR+0.4*m-0.1*mm,0.125*m-0.1*mm,180*deg,pDipoleBendAngle);
	double pLQ1Pos_Z_en2=(pHallCenter2LQ1Face);//NIM
	double pLQ1Pos_Z_ex2=(pHallCenter2LQ1Face + pQ1PosAct);//NIM
	double pLQ2Pos_Z_en2=(pHallCenter2LQ2Face);//NIM
	double pLQ2Pos_Z_ex2=(pHallCenter2LQ2Face + pQ2Length);//NIM
	double pRQ1Pos_Z_en2=(pHallCenter2RQ1Face);//NIM
	double pRQ1Pos_Z_ex2=(pHallCenter2RQ1Face + pQ1PosAct);//NIM
	double pRQ2Pos_Z_en2=(pHallCenter2RQ2Face);//NIM
	double pRQ2Pos_Z_ex2=(pHallCenter2RQ2Face + pQ2Length);//NIM
	double pLDPos_Z_en2 = pLQ2Pos_Z_ex2 + 4.42 * m;
	double pRDPos_Z_en2 = pRQ2Pos_Z_ex2 + 4.42 * m;

	double rInner_DipoleTunnel2[] = {0.,0.,0.,0.};
	double rOuter_DipoleTunnel2[] = {pDipoleR-0.4*m,pDipoleR+0.4*m,pDipoleR+0.4*m,pDipoleR-0.4*m};

	G4Polycone* DipoleTunnelCone1 = new G4Polycone("DipoleTunnelCone1",
	            150.0*deg,105.0*deg,
		    kNPlane_DipoleTunnel,zPlane_DipoleTunnel,rInner_DipoleTunnel2,rOuter_DipoleTunnel2);
	G4Polycone* DipoleTunnelCone2 = new G4Polycone("DipoleTunnelCone2",
	            175.0*deg,pDipoleBendAngle + 10 * deg,
		    kNPlane_DipoleTunnel,zPlane_DipoleTunnel,rInner_DipoleTunnel,rOuter_DipoleTunnel);

	G4Polycone* DipoleTunnelSolid = new G4Polycone("DipoleTunnelSolid",
		180.0*deg,pDipoleBendAngle,
		kNPlane_DipoleTunnel,zPlane_DipoleTunnel,rInner_DipoleTunnel,rOuter_DipoleTunnel);

	G4double psiadj = 180. * deg - 22.5 * deg - 30. * deg;
	G4double badj = 8.4 * m * sin( 30.0 * deg ) / sin( psiadj );
	G4double xadj = badj * sin( 22.5 * deg );
	G4double yadj = badj * cos( 22.5 * deg );
	G4IntersectionSolid* DipoleTunnelCone3 = new G4IntersectionSolid("DipoleTunnelCone3",
			 DipoleTunnelCone2, DipoleTunnelCone1, 0, G4ThreeVector(-yadj, -xadj, 0));
	

	//DipoleTunnelSolid is Jixie with flat entrance
	//DipoleTunnelCone3 is Nickie with 30 deg entrance

	//But that is not good enough. We also have to add extensions to have fringe fields.
	//This is just going to be a trapezoidal extensionpRotDipoleInContaine
	//G4Trd* DipoleFringeSolid = new G4Trd("DipoleFringeSolid", 0.125*m+dy, 0.125*m-dy, 1.4 * m + 0.4 * m * 2. * tan( pi / 6 ), 1.4 * m, 0.4 * m);
	//G4cout << 0.125*m+dy << " and " << 0.125*m-dy << " better be greater than zero!" << G4endl;
	G4Trap* DipoleFringeSolid1 = new G4Trap("DipoleFringeSolid1",
					       0.4 * m , -pi / 12. ,
					       0.0     , 0.125*m-dy,
					       2.05 * m + 0.4 * m * 2. * tan( pi / 12. ), 2.05 * m + 0.4 * m * 2. * tan( pi / 12. ),
					       0.0     , 0.125*m+dy,
					       2.05 * m, 2.05 * m  ,
					       0.0   );
	G4Trap* DipoleFringeSolid2 = new G4Trap("DipoleFringeSolid2",
					       0.4 * m , -pi / 12. ,
					       0.0     , 0.125*m-dy,
					       0.65 * m + 0.4 * m * 2. * tan( pi / 12. ), 0.65 * m + 0.4 * m * 2. * tan( pi / 12. ),
					       0.0     , 0.125*m+dy,
					       0.65 * m, 0.65 * m  ,
					       0.0    );
	
	//pDz  , pTheta,
	//pPhi , pDy1,
	//pDx1 , pDx2,
	//pAlp1, pDy2,
	//pDx3 , pDx4,
	//pAlp2,

	//G4UnionSolid* DipoleTunnelCone4 = new G4UnionSolid("DipoleTunnelCone4",
	//DipoleTunnelCone3, DipoleFringeSolid, pRotDipoleInContainer, G4ThreeVector(-8.4 * m, 0., 0.));
	G4UnionSolid* DipoleTunnelCone4 = new G4UnionSolid("DipoleTunnelCone4",
		   DipoleTunnelCone3, DipoleFringeSolid1, pRotDipoleFringe1InContainer, G4ThreeVector(-8.4 * m, +2.05 * m, 0.));
	G4UnionSolid* DipoleTunnelCone5 = new G4UnionSolid("DipoleTunnelCone4",
		   DipoleTunnelCone4, DipoleFringeSolid2, pRotDipoleFringe2InContainer,
							   G4ThreeVector(-8.4 * m * cos( pi / 4. ) + 0.3 * m / sqrt( 2. ),
									 -8.4 * m * sin( pi / 4. ) - 0.3 * m / sqrt( 2. ), 0.));

	G4LogicalVolume* LDipoleFringe1Logical = new G4LogicalVolume(DipoleFringeSolid1,
	        mMaterialManager->vacuum,"LDipoleFringe1Logical",0,0,0);
	G4LogicalVolume* RDipoleFringe1Logical = new G4LogicalVolume(DipoleFringeSolid1,
	        mMaterialManager->vacuum,"RDipoleFringe1Logical",0,0,0);
	G4LogicalVolume* LDipoleFringe2Logical = new G4LogicalVolume(DipoleFringeSolid2,
	        mMaterialManager->vacuum,"LDipoleFringe2Logical",0,0,0);
	G4LogicalVolume* RDipoleFringe2Logical = new G4LogicalVolume(DipoleFringeSolid2,
	        mMaterialManager->vacuum,"RDipoleFringe2Logical",0,0,0);

	G4LogicalVolume* LDipoleTunnelLogical = new G4LogicalVolume(DipoleTunnelCone5,
		mMaterialManager->vacuum,"LDipoleTunnelLogical",0,0,0);
	G4LogicalVolume* RDipoleTunnelLogical = new G4LogicalVolume(DipoleTunnelCone5,
		mMaterialManager->vacuum,"RDipoleTunnelLogical",0,0,0);
	LDipoleTunnelLogical->SetVisAttributes(YellowVisAtt); 
	RDipoleTunnelLogical->SetVisAttributes(YellowVisAtt); 
	LDipoleFringe1Logical->SetVisAttributes(YellowVisAtt);
	RDipoleFringe1Logical->SetVisAttributes(YellowVisAtt); 
	LDipoleFringe2Logical->SetVisAttributes(YellowVisAtt);
	RDipoleFringe2Logical->SetVisAttributes(YellowVisAtt); 
	//Nickie's dipole field
	G4FieldManager* LdipoleFieldManager = mEMFieldSetup->GetFieldManagerFZBL3();
	G4FieldManager* RdipoleFieldManager = mEMFieldSetup->GetFieldManagerFZBR3();
	//G4double minEps= 1.0e-9;  //   Minimum & value for smallest steps
	//G4double maxEps= 1.0e-8;  //   Maximum & value for largest steps

	//dipoleFieldManager->SetMinimumEpsilonStep( minEps );
	//dipoleFieldManager->SetMaximumEpsilonStep( maxEps );
	//dipoleFieldManager->SetDeltaOneStep( 0.5e-9 * mm );
	//dipoleFieldManager->GetChordFinder()->SetDeltaChord( 0.001 * mm );
	LDipoleTunnelLogical->SetFieldManager(LdipoleFieldManager,allLocal);
	RDipoleTunnelLogical->SetFieldManager(RdipoleFieldManager,allLocal);

        //double pFringeX_ex = ( IS_NIM == 1 ) ?  pQ3CenterY - pQ3Length / sqrt(2.) / 2. - 1.5 * m / sqrt(2) + 0.75 * m / sqrt( 2 ) :   2.4603032 * m + 0.75 * m / sqrt( 2 );
        //double pFringeZ_ex = ( IS_NIM == 1 ) ? -pQ3CenterZ + pQ3Length / sqrt(2.) / 2. + 1.5 * m / sqrt(2) - 0.75 * m / sqrt( 2 ) : -15.9006973 * m - 0.75 * m / sqrt( 2 );

	if(mSetupLHRS>=4)
	{
	  //put it in the container, which also center at the hall center
	  new G4PVPlacement(pRotDipoleInContainer,
			    G4ThreeVector(0,-pDipoleRCenterZ,pDipoleRCenterY),
			    LDipoleTunnelLogical,"LDipoleTunnelPhys",LHRSContainerLogical,0,0,0);
	  /*
	  new G4PVPlacement(pRotDipoleFringe1InContainer,
			    G4ThreeVector(0., -7.75 * m, 0.),
			    LDipoleFringe1Logical,"LDipoleFringe1Phys",LHRSContainerLogical,0,0,0);
	  new G4PVPlacement(pRotDipoleFringe2InContainer,
			    G4ThreeVector(0., pFringeZ_ex, pFringeX_ex),
			    LDipoleFringe2Logical,"LDipoleFringe2Phys",LHRSContainerLogical,0,0,0);
	  */
	}
	if(mSetupRHRS>=4)
	{
	  new G4PVPlacement(pRotDipoleInContainer,
			    G4ThreeVector(0,-pDipoleRCenterZ,pDipoleRCenterY),
			    RDipoleTunnelLogical,"RDipoleTunnelPhys",RHRSContainerLogical,0,0,0);
	  /*
	  new G4PVPlacement(pRotDipoleFringe1InContainer,
			    G4ThreeVector(0., -7.75 * m, 0.),
			    RDipoleFringe1Logical,"RDipoleFringe1Phys",RHRSContainerLogical,0,0,0);
	  new G4PVPlacement(pRotDipoleFringe2InContainer,
			    G4ThreeVector(0., pFringeZ_ex, pFringeX_ex),
			    RDipoleFringe2Logical,"RDipoleFringe2Phys",RHRSContainerLogical,0,0,0);

	  */
	}


	//////////////////////////
	//Q3 tunnel
	//////////////////////////

	G4VSolid* Q3TunnelSolid = new G4Tubs("Q3TunnelTub",0,pQ3Rin-0.1*mm,
		pQ3Length,0.0,360.0*deg);
	G4LogicalVolume* LQ3TunnelLogical = new G4LogicalVolume(Q3TunnelSolid,
		mMaterialManager->vacuum,"LQ3TunnelLogical",0,0,0);
	G4LogicalVolume* RQ3TunnelLogical = new G4LogicalVolume(Q3TunnelSolid,
		mMaterialManager->vacuum,"RQ3TunnelLogical",0,0,0);

	LQ3TunnelLogical->SetVisAttributes(YellowVisAtt); 
	RQ3TunnelLogical->SetVisAttributes(YellowVisAtt); 

	//double pQ3TunnelPos_Y=pDipoleR*(1-cos(pDipoleBendAngle))+2.0*m*sin(pDipoleBendAngle);
	//double pQ3TunnelPos_Z=9.96*m+pDipoleR*sin(pDipoleBendAngle)+2.0*m*cos(pDipoleBendAngle);
	//double pQ3TunnelPos_Z=pDipoleR*sin(pDipoleBendAngle)+2.0*m*cos(pDipoleBendAngle);
	//pQ3TunnelPos_Z += ( IS_NIM == 1 ) ? 9.96*m : pDipoleRCenterZ; 
	//if(mSetupLHRS>=4)
	//{
		//put it in the container, which also center at the hall center
		//new G4PVPlacement(pRotQ3InContainer,G4ThreeVector(0,-pQ3TunnelPos_Z,pQ3TunnelPos_Y),
				  //LQ3TunnelLogical,"LQ1Q2TunnelPhys",LHRSContainerLogical,0,0,0);
	//}
	//if(mSetupRHRS>=4)
	//{
	  //new G4PVPlacement(pRotQ3InContainer,G4ThreeVector(0,-pQ3TunnelPos_Z,pQ3TunnelPos_Y),
	  //RQ3TunnelLogical,"RQ1Q2TunnelPhys",RHRSContainerLogical,0,0,0);
	//}

	//By Jixie: HepRApp.jar has a bug to view this tubs. I am sure this tube is oriented correctly
	//G4VSolid* Q3MagSolid = new G4Tubs("Q3MagTub",0,pQ3Rin,pQ3Length/2.0,0.0,360.0*deg);//NIM
	G4VSolid* Q3MagSolid = new G4Tubs("Q3MagTub",0,pQ3Rin,pQ3LengthMag/2.0,0.0,360.0*deg);//SNAKE

	G4LogicalVolume* LQ3MagLogical = new G4LogicalVolume(Q3MagSolid,
		mMaterialManager->vacuum,"LQ3MagLogical",0,0,0);
	G4LogicalVolume* RQ3MagLogical = new G4LogicalVolume(Q3MagSolid,
		mMaterialManager->vacuum,"RQ3MagLogical",0,0,0);

	//Nickie's Q3 field
	G4FieldManager* LQ3FieldManager = mEMFieldSetup->GetFieldManagerFZBL4();
	G4FieldManager* RQ3FieldManager = mEMFieldSetup->GetFieldManagerFZBR4();
	LQ3MagLogical->SetFieldManager(LQ3FieldManager,allLocal);
	RQ3MagLogical->SetFieldManager(RQ3FieldManager,allLocal);

	LQ3MagLogical->SetVisAttributes(MagFieldVisAtt); 
	RQ3MagLogical->SetVisAttributes(MagFieldVisAtt); 

	G4RotationMatrix *pRotQ3MagInContainer=new G4RotationMatrix();
	pRotQ3MagInContainer->rotateX(-45*deg); 
	if(mSetupLHRS>=4)
	{
	  //put it in the container, which also center at the hall center
	  //new G4PVPlacement(pRotQ3MagInContainer,G4ThreeVector(0,-pQ3CenterZMag,pQ3CenterYMag),//SNAKE
	  //LQ3MagLogical,"LQ3MagPhys",LHRSContainerLogical,0,0,0);//SNAKE
	  new G4PVPlacement(pRotQ3MagInContainer,G4ThreeVector(0,-pQ3CenterZ,pQ3CenterY),//NIM
			    LQ3MagLogical,"LQ3MagPhys",LHRSContainerLogical,0,0,0);//NIM
	}
	if(mSetupRHRS>=4)
	{
	  //new G4PVPlacement(pRotQ3MagInContainer,G4ThreeVector(0,-pQ3CenterZMag,pQ3CenterYMag),//SNAKE
	  //RQ3MagLogical,"RQ3MagPhys",RHRSContainerLogical,0,0,0);//SNAKE
	  new G4PVPlacement(pRotQ3MagInContainer,G4ThreeVector(0,-pQ3CenterZ,pQ3CenterY),//NIM
	    		    RQ3MagLogical,"RQ3MagPhys",RHRSContainerLogical,0,0,0);//NIM
	}


	//Nickie's sensitive detector, at the focal plane, virtual boundary, vb, fp
	//double pFPR      = pDipoleR; //radius of curvature of dipole
	//double pFPA      = 9.96 * m; //distance from pivot to dipole entrance
	//double pFPB      = 2.4  * m; //distance from center of dipole to center of Q3
	//
	//double pFPB      = 2.4  * m + pQ3Length / 2.0 + 3.57 * m + 1.43 * m;
	//double pFPB      = 2.4  * m + pQ3Length / 2.0 + 3.57 * m + 0.7 * m;
	//double pFPCenterY=       pFPR * ( 1 - cos( pDipoleBendAngle ) ) + pFPB * sin( pDipoleBendAngle );
	//double pFPCenterZ=pFPA + pFPR *       sin( pDipoleBendAngle )   + pFPB * cos( pDipoleBendAngle );
        //double pFPA      = 9.96 * m;//distance from pivot to entrance of dipole
        //double pFPCenterX=( pFPA + pFPH + ( pFPH + 1.5 * m + 1.8 * m + 3.57 * m + 1.43 * m ) / sqrt(2) ) *-sin( mRHRSAngle );
        //double pFPCenterZ=( pFPA + pFPH + ( pFPH + 1.5 * m + 1.8 * m + 3.57 * m + 1.43 * m ) / sqrt(2) ) * cos( mRHRSAngle );

	double pFPCenterZ , pFPCenterY ;
	double pVDCCenterZ, pVDCCenterY;
	double pQZ1CenterZ, pQZ1CenterY;
	double pQZ2CenterZ, pQZ2CenterY;
	//if( IS_NIM == 1 ){
	//double pFPR      = 8.4  * m;//radius of curvature of dipole
	//double pFPA      = ( IS_NIM == 1 ) ? 9.96 * m : pDipoleRCenterZ;//distance from pivot to entrance of dipole
	//double pFPH      = pFPR * tan ( 45. / 2. * deg );
	//pFPCenterZ=( pFPA + pFPH + ( pFPH + 1.5 * m + 1.8 * m + 3.57 * m + 1.43 * m ) / sqrt(2) );
	//pFPCenterY=(                 pFPH + 1.5 * m + 1.8 * m + 3.57 * m + 1.43 * m ) / sqrt(2);
	//}else{
	pVDCCenterZ = pQ3CenterZ + ( pQ3Length / 2. + 3.4538 * m ) / sqrt(2) ;
	pVDCCenterY = pQ3CenterY + ( pQ3Length / 2. + 3.4538 * m ) / sqrt(2) ;
	pQZ1CenterZ = pQ3CenterZ + ( pQ3Length / 2. + 3.4538 * m ) / sqrt(2) + 2.6 * cm + 33.5 * cm + 2.6 * cm + 40.0 * cm;
	pQZ1CenterY = pQ3CenterY + ( pQ3Length / 2. + 3.4538 * m ) / sqrt(2) + 2.6 * cm + 33.5 * cm + 2.6 * cm + 40.0 * cm;
	pQZ2CenterZ = pQ3CenterZ + ( pQ3Length / 2. + 3.4538 * m ) / sqrt(2) + 2.6 * cm + 33.5 * cm + 2.6 * cm + 52.0 * cm;
	pQZ2CenterY = pQ3CenterY + ( pQ3Length / 2. + 3.4538 * m ) / sqrt(2) + 2.6 * cm + 33.5 * cm + 2.6 * cm + 52.0 * cm;
	pFPCenterZ  = pQ3CenterZ + ( pQ3Length / 2. + 3.4538 * m + 1.43 * m ) / sqrt(2) ;
	pFPCenterY  = pQ3CenterY + ( pQ3Length / 2. + 3.4538 * m + 1.43 * m ) / sqrt(2) ;
	//}
	
	G4double vb_thickness = .00001 * mm;
	G4VSolid* FPSolid = new G4Tubs("FPTub",0,pQ3Rout * 2,vb_thickness,0.0,360.0*deg);
	G4VSolid* PlaneSolid1 = new G4Tubs("PlaneTub",0,pQ1Rout,vb_thickness,0.0,360.0*deg); //circles
	G4VSolid* PlaneSolid2 = new G4Tubs("PlaneTub",0,pQ2Rout,vb_thickness,0.0,360.0*deg); //circles

	double X_S0 = 29.5 * cm;
	double Y_S0 = 35.5 * cm;
	double Z_S0 = vb_thickness;

	double X_S1 = 29.5 * cm;
	double Y_S1 = 35.5 * cm;
	double Z_S1 = vb_thickness;

	double X_Q1 = 29.5 * cm;
	double Y_Q1 = 35.5 * cm;
	double Z_Q1 = vb_thickness;

	double X_Q2 = 29.5 * cm;
	double Y_Q2 = 35.5 * cm;
	double Z_Q2 = vb_thickness;

	G4VSolid* Scint0  = new G4Box("Scint0" , X_S0 / 2.0, Y_S0 / 2.0, Z_S0 / 2.0);
	G4VSolid* Scint1  = new G4Box("Scint1" , X_S1 / 2.0, Y_S1 / 2.0, Z_S1 / 2.0);
	G4VSolid* Quartz1 = new G4Box("Quartz1", X_Q1 / 2.0, Y_Q1 / 2.0, Z_Q1 / 2.0);
	G4VSolid* Quartz2 = new G4Box("Quartz2", X_Q2 / 2.0, Y_Q2 / 2.0, Z_Q2 / 2.0);

	G4VSolid* magneticSolid = new G4Box("magneticBox",mFieldX/2.0,mFieldY/2.0,mFieldZ/2.0);

	G4LogicalVolume* LFPLogical = new G4LogicalVolume(FPSolid,
		mMaterialManager->vacuum,"LFPLogical",0,0,0);
	G4LogicalVolume* RFPLogical = new G4LogicalVolume(FPSolid,
		mMaterialManager->vacuum,"RFPLogical",0,0,0);

	G4LogicalVolume* LPlaneLogical1 = new G4LogicalVolume(PlaneSolid1,
		mMaterialManager->vacuum,"LPlaneLogical1",0,0,0);
	G4LogicalVolume* RPlaneLogical1 = new G4LogicalVolume(PlaneSolid1,
		mMaterialManager->vacuum,"RPlaneLogical1",0,0,0);
	G4LogicalVolume* LPlaneLogical2 = new G4LogicalVolume(PlaneSolid2,
		mMaterialManager->vacuum,"LPlaneLogical1",0,0,0);
	G4LogicalVolume* RPlaneLogical2 = new G4LogicalVolume(PlaneSolid2,
		mMaterialManager->vacuum,"RPlaneLogical1",0,0,0);

	LFPLogical->SetVisAttributes(MagFieldVisAtt); 
	RFPLogical->SetVisAttributes(MagFieldVisAtt); 
	LPlaneLogical1->SetVisAttributes(MagFieldVisAtt); 
	RPlaneLogical1->SetVisAttributes(MagFieldVisAtt); 
	LPlaneLogical2->SetVisAttributes(MagFieldVisAtt); 
	RPlaneLogical2->SetVisAttributes(MagFieldVisAtt); 

	G4RotationMatrix *pRotVDCInContainer=new G4RotationMatrix();
	pRotVDCInContainer->rotateX(0.*deg); 
	G4RotationMatrix *pRotFPInContainer=new G4RotationMatrix();
	pRotFPInContainer->rotateX(-45*deg); 
	//if(mSetupLHRS>=4){
	if(mSnakeModel == 49 || mSnakeModel == 48 || mSnakeModel > 51 ){
	  double pSeptumX      = 140.0  * cm;
	  double pSeptumY      = 84.4   * cm;
	  double pSeptumZ      = 74.0   * cm;
	  //double pSeptumPlaceZ = 70.414 * cm;
	  double pSeptumPlaceZ = 69.99937 * cm;
	  new G4PVPlacement(0,G4ThreeVector(0,0,pSeptumPlaceZ - 0.5 * pSeptumZ + 2 * vb_thickness),
			    LPlaneLogical2,"virtualBoundaryPhys_sen",motherLogical,0,0,0);//sen
	  new G4PVPlacement(0,G4ThreeVector(0,0,pSeptumPlaceZ),
			    LPlaneLogical2,"virtualBoundaryPhys_sm",motherLogical,0,0,0);//sm
	  new G4PVPlacement(0,G4ThreeVector(0,0,pSeptumPlaceZ + 0.5 * pSeptumZ),
			    LPlaneLogical2,"virtualBoundaryPhys_sex",motherLogical,0,0,0);//sex
	  new G4PVPlacement(0,G4ThreeVector(0,0,36. * cm),
			    LPlaneLogical2,"virtualBoundaryPhys_coil",motherLogical,0,0,0);//coil
	  new G4PVPlacement(0,G4ThreeVector(0,0,-50. * cm),
			    LPlaneLogical2,"virtualBoundaryPhys_mid",motherLogical,0,0,0);//mid
	  
	}
	//double pLQ1Pos_Z_en=(pHallCenter2LQ1Face);//NIM
	//double pLQ1Pos_Z_ex=(pHallCenter2LQ1Face + pQ1Length);//NIM
	//double pLQ2Pos_Z_en=(pHallCenter2LQ2Face);//NIM
	//double pLQ2Pos_Z_ex=(pHallCenter2LQ2Face + pQ2Length);//NIM
	//double pRQ1Pos_Z_en=(pHallCenter2RQ1Face);//NIM
	//double pRQ1Pos_Z_ex=(pHallCenter2RQ1Face + pQ1Length);//NIM
	//double pRQ2Pos_Z_en=(pHallCenter2RQ2Face);//NIM
	//double pRQ2Pos_Z_ex=(pHallCenter2RQ2Face + pQ2Length);//NIM
	//double pLDPos_Z_en = pLQ2Pos_Z_ex + 4.42 * m;
	//double pRDPos_Z_en = pRQ2Pos_Z_ex + 4.42 * m;
	//double pLDPos_X_ex = ( IS_NIM == 1 ) ?  pQ3CenterY - pQ3Length / sqrt(2.) / 2. - 1.5 * m / sqrt(2) : 2.4603032 * m;
	//double pRDPos_X_ex = ( IS_NIM == 1 ) ?  pQ3CenterY - pQ3Length / sqrt(2.) / 2. - 1.5 * m / sqrt(2) : 2.4603032 * m;
	//double pLDPos_Z_ex = ( IS_NIM == 1 ) ? -pQ3CenterZ + pQ3Length / sqrt(2.) / 2. + 1.5 * m / sqrt(2) : -15.9006973 * m;
	//double pRDPos_Z_ex = ( IS_NIM == 1 ) ? -pQ3CenterZ + pQ3Length / sqrt(2.) / 2. + 1.5 * m / sqrt(2) : -15.9006973 * m;

	double pHallCenter2Col = 1.38 * m;
	double pPaulColT        = 0.01 * m;
	double pPaulX = ( - pHallCenter2Col - pPaulColT * 2. ) * cos(mLHRSAngle) ;
	double pPaulY = ( - pHallCenter2Col - pPaulColT * 2. ) * sin(mLHRSAngle);
	if(mSnakeModel == 49 || mSnakeModel == 48 || mSnakeModel > 50 ){
	  //double pLQ1Pos_Z=(pHallCenter2LQ1FaceMag+pQ1LengthMag/1.0);//SNAKE
	  new G4PVPlacement(pRotXLHRSdeg,G4ThreeVector(pPaulY,0,-pPaulX),
			    LPlaneLogical1,"virtualBoundaryPhys_col_LHRS",motherLogical,0,0,0);//q1en
	  /*
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pLQ1Pos_Z_en,0),
			    LPlaneLogical1,"virtualBoundaryPhys_q1en_LHRS",LHRSContainerLogical,0,0,0);//q1en
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pLQ1Pos_Z_ex,0),
			    LPlaneLogical1,"virtualBoundaryPhys_q1ex_LHRS",LHRSContainerLogical,0,0,0);//q1ex
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pLQ2Pos_Z_en,0),
			    LPlaneLogical2,"virtualBoundaryPhys_q2en_LHRS",LHRSContainerLogical,0,0,0);//q2en
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pLQ2Pos_Z_ex,0),
			    LPlaneLogical2,"virtualBoundaryPhys_q2ex_LHRS",LHRSContainerLogical,0,0,0);//q2ex
	  */
	  new G4PVPlacement(0,G4ThreeVector(0, 0,-pQ1Length/2.),
			    LPlaneLogical1,"virtualBoundaryPhys_q1en_LHRS",LQ1MagLogical,0,0,0);//q1en
	  new G4PVPlacement(0,G4ThreeVector(0, 0, pQ1Length/2.),
			    LPlaneLogical1,"virtualBoundaryPhys_q1ex_LHRS",LQ1MagLogical,0,0,0);//q1ex
	  new G4PVPlacement(0,G4ThreeVector(0, 0,-pQ2Length/2.),
			    LPlaneLogical2,"virtualBoundaryPhys_q2en_LHRS",LQ2MagLogical,0,0,0);//q2en
	  new G4PVPlacement(0,G4ThreeVector(0, 0, pQ2Length/2.),
			    LPlaneLogical2,"virtualBoundaryPhys_q2ex_LHRS",LQ2MagLogical,0,0,0);//q2ex
	  //new G4PVPlacement(pRotX30deg,G4ThreeVector(0,-pLDPos_Z_en,0),
	  //LPlaneLogical2,"virtualBoundaryPhys_den_LHRS",LHRSContainerLogical,0,0,0);//den
	  //new G4PVPlacement(pRotX105deg,G4ThreeVector(0, pLDPos_Z_ex,pLDPos_X_ex),
	  //LPlaneLogical2,"virtualBoundaryPhys_dex_LHRS",LHRSContainerLogical,0,0,0);//dex
	  
	  new G4PVPlacement(pRot_den,G4ThreeVector(-8.4*m,0,0),
			    LPlaneLogical2,"virtualBoundaryPhys_den_LHRS",LDipoleTunnelLogical,0,0,0);//den
	  new G4PVPlacement(pRot_dex,G4ThreeVector(-8.4 * m * cos( pi / 4. ), -8.4 * m * sin( pi / 4. ),0),
			    LPlaneLogical2,"virtualBoundaryPhys_dex_LHRS",LDipoleTunnelLogical,0,0,0);//dex
	  
	  //new G4PVPlacement(pRotX45deg,G4ThreeVector(0,-pQ3CenterZ + pQ3Length / sqrt(2) / 2,pQ3CenterY - pQ3Length / sqrt(2) / 2),
	  //LPlaneLogical2,"virtualBoundaryPhys_q3en_LHRS",LHRSContainerLogical,0,0,0);//q3en
	  //new G4PVPlacement(pRotX45deg,G4ThreeVector(0,-pQ3CenterZ - pQ3Length / sqrt(2) / 2,pQ3CenterY + pQ3Length / sqrt(2) / 2),
	  //LPlaneLogical2,"virtualBoundaryPhys_q3ex_LHRS",LHRSContainerLogical,0,0,0);//q3ex
	  new G4PVPlacement(0,G4ThreeVector(0,0,-pQ3Length/2.),
			    LPlaneLogical2,"virtualBoundaryPhys_q3en_LHRS",LQ3MagLogical,0,0,0);//q3en
	  new G4PVPlacement(0,G4ThreeVector(0,0, pQ3Length/2.),
			    LPlaneLogical2,"virtualBoundaryPhys_q3ex_LHRS",LQ3MagLogical,0,0,0);//q3ex
	  new G4PVPlacement(pRotVDCInContainer,G4ThreeVector(0,-pVDCCenterZ,pVDCCenterY),
			    LFPLogical,"virtualBoundaryPhys_vdc_LHRS",LHRSContainerLogical,0,0,0);
	  new G4PVPlacement(pRotFPInContainer,G4ThreeVector(0,-pFPCenterZ,pFPCenterY),
			    LFPLogical,"virtualBoundaryPhys_fp_LHRS",LHRSContainerLogical,0,0,0);
	  new G4PVPlacement(pRotFPInContainer,G4ThreeVector(0,-pQZ1CenterZ,pQZ1CenterY),
			    LFPLogical,"virtualBoundaryPhys_qz1_LHRS",LHRSContainerLogical,0,0,0);
	  new G4PVPlacement(pRotVDCInContainer,G4ThreeVector(0,-pQZ2CenterZ,pQZ2CenterY),
			    LFPLogical,"virtualBoundaryPhys_qz2_LHRS",LHRSContainerLogical,0,0,0);

	}
	//if(mSetupRHRS>=4){
	if(mSnakeModel == 49 || mSnakeModel == 48 || mSnakeModel > 50 ){
	  new G4PVPlacement(pRotXRHRSdeg,G4ThreeVector(-pPaulY, 0, -pPaulX),
			    LPlaneLogical1,"virtualBoundaryPhys_col_RHRS",motherLogical,0,0,0);//q1en

	  //double pRQ1Pos_Z=(pHallCenter2LQ1FaceMag+pQ1LengthMag/1.0);//SNAKE
	  /*
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pRQ1Pos_Z_en,0),
			    RPlaneLogical1,"virtualBoundaryPhys_q1en_RHRS",RHRSContainerLogical,0,0,0);//q1en
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pRQ1Pos_Z_ex,0),
			    RPlaneLogical1,"virtualBoundaryPhys_q1ex_RHRS",RHRSContainerLogical,0,0,0);//q1ex
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pRQ2Pos_Z_en,0),
			    RPlaneLogical2,"virtualBoundaryPhys_q2en_RHRS",RHRSContainerLogical,0,0,0);//q2en
	  new G4PVPlacement(pRotX90deg,G4ThreeVector(0,-pRQ2Pos_Z_ex,0),
			    RPlaneLogical2,"virtualBoundaryPhys_q2ex_RHRS",RHRSContainerLogical,0,0,0);//q2ex
	  */
	  new G4PVPlacement(0,G4ThreeVector(0,0,-pQ1Length/2.),
			    RPlaneLogical1,"virtualBoundaryPhys_q1en_RHRS",RQ1MagLogical,0,0,0);//q1en
	  new G4PVPlacement(0,G4ThreeVector(0,0, pQ1Length/2.),
			    RPlaneLogical1,"virtualBoundaryPhys_q1ex_RHRS",RQ1MagLogical,0,0,0);//q1ex
	  new G4PVPlacement(0,G4ThreeVector(0,0,-pQ2Length/2.),
			    RPlaneLogical2,"virtualBoundaryPhys_q2en_RHRS",RQ2MagLogical,0,0,0);//q2en
	  new G4PVPlacement(0,G4ThreeVector(0,0, pQ2Length/2.),
			    RPlaneLogical2,"virtualBoundaryPhys_q2ex_RHRS",RQ2MagLogical,0,0,0);//q2ex
	  //new G4PVPlacement(pRotX30deg,G4ThreeVector(0,-pRDPos_Z_en,0),
	  //RPlaneLogical2,"virtualBoundaryPhys_den_RHRS",RHRSContainerLogical,0,0,0);//den
	  //new G4PVPlacement(pRotX105deg,G4ThreeVector(0, pRDPos_Z_ex,pRDPos_X_ex),
	  //RPlaneLogical2,"virtualBoundaryPhys_dex_RHRS",RHRSContainerLogical,0,0,0);//dex
	  
	  new G4PVPlacement(pRot_den,G4ThreeVector(-8.4*m,0,0),
			    RPlaneLogical2,"virtualBoundaryPhys_den_RHRS",RDipoleTunnelLogical,0,0,0);//den
	  new G4PVPlacement(pRot_dex,G4ThreeVector(-8.4 * m * cos( pi / 4. ), -8.4 * m * sin( pi / 4. ),0),
			    RPlaneLogical2,"virtualBoundaryPhys_dex_RHRS",RDipoleTunnelLogical,0,0,0);//dex
	  
	  //new G4PVPlacement(pRotX45deg,G4ThreeVector(0,-pQ3CenterZ + pQ3Length / sqrt(2) / 2,pQ3CenterY - pQ3Length / sqrt(2) / 2),
	  //RPlaneLogical2,"virtualBoundaryPhys_q3en_RHRS",RHRSContainerLogical,0,0,0);//q3en
	  //new G4PVPlacement(pRotX45deg,G4ThreeVector(0,-pQ3CenterZ - pQ3Length / sqrt(2) / 2,pQ3CenterY + pQ3Length / sqrt(2) / 2),
	  //RPlaneLogical2,"virtualBoundaryPhys_q3ex_RHRS",RHRSContainerLogical,0,0,0);//q3ex
	  new G4PVPlacement(0,G4ThreeVector(0,0,-pQ3Length/2.),
			    RPlaneLogical2,"virtualBoundaryPhys_q3en_RHRS",RQ3MagLogical,0,0,0);//q3en
	  new G4PVPlacement(0,G4ThreeVector(0,0, pQ3Length/2.),
			    RPlaneLogical2,"virtualBoundaryPhys_q3ex_RHRS",RQ3MagLogical,0,0,0);//q3ex
	  new G4PVPlacement(pRotVDCInContainer,G4ThreeVector(0,-pVDCCenterZ,pVDCCenterY),
			    RFPLogical,"virtualBoundaryPhys_vdc_RHRS",RHRSContainerLogical,0,0,0);
	  new G4PVPlacement(pRotFPInContainer,G4ThreeVector(0,-pFPCenterZ,pFPCenterY),
			    RFPLogical,"virtualBoundaryPhys_fp_RHRS",RHRSContainerLogical,0,0,0);
	  new G4PVPlacement(pRotFPInContainer,G4ThreeVector(0,-pQZ1CenterZ,pQZ1CenterY),
			    RFPLogical,"virtualBoundaryPhys_qz1_RHRS",RHRSContainerLogical,0,0,0);
	  new G4PVPlacement(pRotVDCInContainer,G4ThreeVector(0,-pQZ2CenterZ,pQZ2CenterY),
			    RFPLogical,"virtualBoundaryPhys_qz2_RHRS",RHRSContainerLogical,0,0,0);

	}

	//#endif
	//////////////////////////////////////////////////////////

	return theHRSPhys;
}


/////////////////////////////////////////////////////////////////////


G4VPhysicalVolume* HRSDetectorConstruction::ConstructRadiator(G4LogicalVolume* motherLogical)
{
	G4String SDname;

	// All G4VSensitiveDetector will be managed (deleted) by SDManager
	//therefore no need to delete it at the end of this subroutine
	G4VSensitiveDetector* virtualDetectorSD=new HRSStdSD(SDname="virtualDetector_Rad");

	// sensitive detectors   
	G4SDManager* SDman = G4SDManager::GetSDMpointer();
	G4VPhysicalVolume* theRadiatorPhys = 0;

	////////////////////////
	G4Material *theRadMaterial = mMaterialManager->copper;
	if(mRadiatorType == 1) theRadMaterial = mMaterialManager->iron;
	else if(mRadiatorType == 2) theRadMaterial = mMaterialManager->copper;
	else if(mRadiatorType == 3) theRadMaterial = mMaterialManager->tantalum;
	else if(mRadiatorType == 4) theRadMaterial = mMaterialManager->tungsten;


	/////////////////////////
	//the radiator
	/////////////////////////
	double pRadX=4*cm, pRadY=4*cm;
	double pRadZ=mRaditorThickness;

	G4VSolid* radiatorSolid = new G4Box("radiatorBox",pRadX/2.0,pRadY/2.0,pRadZ/2.0);
	G4LogicalVolume* radiatorLogical = new G4LogicalVolume(radiatorSolid,
		theRadMaterial,"radiatorLogical",0,0,0);
	radiatorLogical->SetVisAttributes(OrangeVisAtt);


	double pRadPos_X=mPivotXOffset, pRadPos_Y=mPivotYOffset;
	double pRadPos_Z=mPivotYOffset + mRadiator2Pivot;
	theRadiatorPhys = new G4PVPlacement(0,
		G4ThreeVector(pRadPos_X,pRadPos_Y,pRadPos_Z),
		radiatorLogical,"radiatorPhys",motherLogical,0,0);


	/////////////////////////
	//the radiator virtual ditector
	/////////////////////////
	if(mSetupRadiatorVD)
	{
		double pRadVDX=20*cm, pRadVDY=20*cm, pRadVDZ=2*mm;

		G4VSolid* radiatorVDSolid = new G4Box("radiatorVDBox",
			pRadVDX/2.0,pRadVDY/2.0,pRadVDZ/2.0);
		G4LogicalVolume* radiatorVDLogical = new G4LogicalVolume(radiatorVDSolid,
			mMaterialManager->air,"radiatorVDLogical",0,0,0);

		SDman->AddNewDetector(virtualDetectorSD);
		radiatorVDLogical->SetSensitiveDetector(virtualDetectorSD);
		radiatorVDLogical->SetVisAttributes(LightYellowVisAtt);


		//place the VD 2 cm downstream from the downstream end plane of the radiator
		double pRadVDPos_X=pRadPos_X, pRadVDPos_Y=pRadPos_Y;
		double pRadVDPos_Z=pRadPos_Z + mRaditorThickness/2 + 2*cm;
		G4String VDPhysName = (mSetupRadiatorVD==2) ? "virtualBoundaryPhys_Rad" : "radiatorVDPhys";
		new G4PVPlacement(0,G4ThreeVector(pRadVDPos_X,pRadVDPos_Y,pRadVDPos_Z),
			radiatorVDLogical,VDPhysName,motherLogical,0,0);
	}

	return theRadiatorPhys;
}

/////////////////////////////////////////////////////////////////////

