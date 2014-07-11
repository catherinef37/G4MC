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

#include "G4VSolid.hh"
#include "G4Box.hh"
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

#include "BigBiteDetectorConstruction.hh"
#include "SBSDetectorConstruction.hh"
#include "G2PDetectorConstruction.hh"
#include "RTPCDetectorConstruction.hh"
#include "CREXDetectorConstruction.hh"
#include "HMSDetectorConstruction.hh"
#include "LACDetectorConstruction.hh"

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
	mSetupG2PGeometry=0;
	gConfig->GetParameter("SetupG2PGeometry",mSetupG2PGeometry);

	mSetupCREXGeometry=0;
	gConfig->GetParameter("SetupCREXGeometry",mSetupCREXGeometry);

	mSetupRTPCGeometry=0;
	gConfig->GetParameter("SetupRTPCGeometry",mSetupRTPCGeometry);

	mSetupBigBite=0;
	gConfig->GetParameter("SetupBigBite",mSetupBigBite);

	gConfig->GetParameter("SetupSuperBigBite",mSetupSuperBigBite);

	gConfig->GetParameter("SetupHMS",mSetupHMS);

	gConfig->GetParameter("SetupLAC",mSetupLAC);

	G4cout<<"\n****Load detector config parameters done!***"<<G4endl;

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
	G4VSensitiveDetector* virtualDetectorSD=new HRSStdSD(SDname="virtualDetector");
	G4VSensitiveDetector* virtualBoundarySD=new HRSStdSD(SDname="virtualBoundary");

	// sensitive detectors   
	G4SDManager* SDman = G4SDManager::GetSDMpointer();

	/////////////////////////////////////////////////////////////////////////////

	//set "user limits" for drawing smooth curve, only supported by Jixie's model
	G4UserLimits* uHallStepLimits = new G4UserLimits(1000.*mm);
	double pBStepLimit=100;
	gConfig->GetArgument("BStepLimit",pBStepLimit);
	pBStepLimit*=mm;
	G4UserLimits* uBStepLimits = new G4UserLimits(pBStepLimit);

	///////////////////////////////////////////////////////////////////////////////
	// Magnetic field ----------------------------------------------------------
	mEMFieldSetup = new HRSEMFieldSetup();  //setup the field, 
	///////////////////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////////////
	//Shrink the hall size if bigbite is not built, this will mke it run fast
	//note that this is optimized for g2p only
	//These 2 from Detector_G2P.ini
	int pSetupChicane=1,pSetupThirdArm=1;
	if (!mSetupG2PGeometry) 
	{
		pSetupChicane=pSetupThirdArm=0;
	}

	if( mSetupBigBite==0 && mSetupSuperBigBite==0 && pSetupChicane==0 && 
		(mSetupVirtualDetector==0 || mPivot2VDFace<8000*mm)  && mSetupHMS<=1 
		&& mSetupLAC<1 )
	{		
		if(mSetupLHRS<2 && mSetupRHRS<2)
		{
			mFieldX=640.0*cm;
			mFieldY=600.0*cm;
			mFieldZ=800.0*cm;
			if(pSetupThirdArm>0)
			{
				mFieldX=800.0*cm;
				mFieldY=800.0*cm;
			}
			if(mSetupVirtualDetector>0 && mPivot2VDFace>4000*mm)
			{
				mFieldX=2*(mPivot2VDFace+2000*mm)*fabs(sin(mVDRotYAngle))+2000*mm;
				mFieldY=mFieldX;
				mFieldZ=2*(mPivot2VDFace+2000*mm)*fabs(cos(mVDRotYAngle))+2000*mm;
				if(mFieldZ<800*cm) mFieldZ=800.0*cm;
			}
		}
		else
		{
			mFieldX=mHallZ*max(fabs(sin(mLHRSAngle)),fabs(sin(mRHRSAngle)));
			mFieldY=mFieldX;
		}

		mHallX=mFieldX+10*cm;
		mHallY=mFieldY+10*cm;
		mHallZ=mFieldZ+10*cm;		
		G4cout<<"By Jixie: Since BigBite, SBS and Chicane are not used, shrink the hall to a smaller\n"
			<<"size of "<< mHallX/m<<" X "<<mHallY/m<<" X "<<mHallZ/m
			<<" (m), which might speed up a lot." <<G4endl;
	}
	///////////////////////////////////////////////////////////////////////////

	// start to construct the geometries -----------------------------------------

	/////////////////////////
	// experimental hall (hall volume)
	/////////////////////////
	G4VSolid* hallSolid = new G4Box("hallBox",mHallX/2.0,mHallY/2.0,mHallZ/2.0);
	G4LogicalVolume* hallLogical = new G4LogicalVolume(hallSolid,
		mMaterialManager->air,"hallLogical",0,0,uHallStepLimits);
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

		G4VSolid* xAxisSolid = new G4Tubs("xAxisTubs",0.,0.5*mm,50.0*mm,0.,360.*deg);
		G4LogicalVolume* xAxisLogical = new G4LogicalVolume(xAxisSolid,
			mMaterialManager->air,"xAxisLogical",0,0,0);
		new G4PVPlacement(pRotY270deg,G4ThreeVector(650*mm,0,0),
			xAxisLogical,"xAxisPhys",hallLogical,0,0);
		xAxisLogical->SetVisAttributes(G4VisAttributes(G4Colour::Red()));

		G4VSolid* yAxisSolid = new G4Tubs("yAxisTubs",0.,0.5*mm,50.0*mm,0.,360.*deg);
		G4LogicalVolume* yAxisLogical = new G4LogicalVolume(yAxisSolid,
			mMaterialManager->air,"yAxisLogical",0,0,0);
		new G4PVPlacement(pRotX270deg,G4ThreeVector(0,650*mm,0),
			yAxisLogical,"yAxisPhys",hallLogical,0,0);
		yAxisLogical->SetVisAttributes(G4VisAttributes(G4Colour::Green()));

		G4VSolid* zAxisSolid = new G4Tubs("zAxisTubs",0.,0.5*mm,50.0*mm,0.,360.*deg);
		G4LogicalVolume* zAxisLogical = new G4LogicalVolume(zAxisSolid,
			mMaterialManager->air,"zAxisLogical",0,0,0);
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
	G4LogicalVolume* magneticLogical = new G4LogicalVolume(magneticSolid,
		pBMaterial,"magneticLogical",0,0,0);
	magneticLogical->SetUserLimits(uBStepLimits);
	magneticLogical->SetVisAttributes(MagneticVisAtt);

	//the field region is the mother volumn of a lot of dauthers, it should not rotate
	new G4PVPlacement(0,G4ThreeVector(),magneticLogical,"magneticPhys",hallLogical,0,0);



	/////////////////////////
	//the radiator
	/////////////////////////
	if(mSetupRadiator) ConstructRadiator(magneticLogical);


	/////////////////////////
	//g2p geometries
	/////////////////////////
	//g2pscattering chamber, target, target coil, local dump, 
	//g2p sieve, g2p septum, thirdarm, HRSVB, chicane, platform
	if(mSetupG2PGeometry)
	{
		G2PDetectorConstruction* theG2P = new G2PDetectorConstruction(magneticLogical); 
		theG2P->Construct();
		//update the parameters
		this->GetConfig(); 
	}

	/////////////////////////
	//RTPC geometries
	/////////////////////////
	if (mSetupRTPCGeometry) 
	{
		RTPCDetectorConstruction* theRTPC = new RTPCDetectorConstruction(magneticLogical); 
		theRTPC->Construct();
		//update the parameters
		this->GetConfig(); 
	}

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
	// BigBite
	/////////////////////////
	if(mSetupBigBite) 
	{
#if defined G4DEBUG_GEOMETRY && (G4DEBUG_GEOMETRY>=2)

		//after 1 day of debugging, I found that G4PVPlacement() only take a 
		//G4RotationMatrix pointer as the first argument, not an object|instance
		//in Geant4_9.X.X. That is why the bigbite detector never shows up in the 
		//visualization in Geant4_9.4.3(but Geant4_8.x do accept an object .....)
		//The following code is for this test
		double dist2targ,detangle;
		dist2targ = 200*cm;//110.*cm; //from front face to the target center 

		detangle=270*deg;
		G4RotationMatrix rotbigbite;
		rotbigbite.rotateY(-detangle+90*deg); 
		G4RotationMatrix* pRotBigBite=new G4RotationMatrix();
		pRotBigBite->rotateY(-90*deg); 

		G4double xmb = 10.*m; //
		G4double ymb = 2.*m;  // Size of this box, large enough but no too big
		G4double zmb = 2.*m;  //  
		G4Box* mboxlong = new G4Box("motherboxlong",xmb/2,ymb/2,zmb/2);
		G4Box* mboxsub = new G4Box("motherboxsub",xmb/2,ymb/2+1.0*cm,zmb/2+1.0*cm);
		//this is the old method
		//G4Box* mbox = new G4Box("motherbox",xmb/2,ymb/2,zmb/2);
		G4SubtractionSolid* mbox = new G4SubtractionSolid("motherbox",
			mboxlong,mboxsub,
			0,G4ThreeVector(-(xmb/2-dist2targ),0,0));

		G4LogicalVolume* logmb = new G4LogicalVolume(mbox, mMaterialManager->air,
			"logmb", 0, 0, 0);
		logmb->SetVisAttributes(OrangeVisAtt);

		new G4PVPlacement(&rotbigbite,G4ThreeVector(0,0,0),
			logmb,"MotherBox",magneticLogical,0,0);
		new G4PVPlacement(pRotBigBite,G4ThreeVector(0,0,-3),
			logmb,"MotherBox",magneticLogical,0,0);			
#endif

		BigBiteDetectorConstruction* theBB = new BigBiteDetectorConstruction(magneticLogical); 
		theBB->Construct();
		//update the parameters
		this->GetConfig(); 
	}

	/////////////////////////

	/////////////////////////
	// SBS 
	/////////////////////////
	if(mSetupSuperBigBite)  
	{
		SBSDetectorConstruction* theSBS = new SBSDetectorConstruction(magneticLogical); 
		theSBS->Construct();
		//update the parameters
		this->GetConfig(); 
	}

	/////////////////////////
	// HMS 
	/////////////////////////
	if(mSetupHMS)  
	{
		HMSDetectorConstruction* theHMS = new HMSDetectorConstruction(magneticLogical); 
		theHMS->Construct();
		//update the parameters
		this->GetConfig(); 
	}

	/////////////////////////
	// LAC
	/////////////////////////
	if(mSetupLAC)  
	{
		LACDetectorConstruction* theLAC = new LACDetectorConstruction(magneticLogical); 
		theLAC->Construct();
		//update the parameters
		this->GetConfig(); 
	}

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

	/////////////////////////
	// G2P Virtual Detector, the 4th arm?
	/////////////////////////
	//bool mSetupVirtualDetector=true;  //will be read from Detector.ini
	if(mSetupVirtualDetector)
	{	
		//double mVirtualDetectorWidth=260.0*mm;
		//double mVirtualDetectorHeight=400.0*mm;
		//double mVirtualDetectorThick=5.0*mm;
		//double mVDRotYAngle=-15.0*deg;
		//double mVDRotXAngle=-8.0*deg;
		//double mPivot2VDFace=635.0*mm;
		//string mVDPhysVolName="virtualBoundaryPhys";

		G4VSolid* virtualDetectorSolid = new G4Box("virtualDetectorBox",mVirtualDetectorWidth/2.0,
			mVirtualDetectorHeight/2.0,mVirtualDetectorThick/2.0);
		G4LogicalVolume* virtualDetectorLogical = new G4LogicalVolume(virtualDetectorSolid,
			mMaterialManager->heliumGas,"virtualDetectorLogical",0,0,0);
		virtualDetectorLogical->SetVisAttributes(LightYellowVisAtt);  
		SDman->AddNewDetector(virtualDetectorSD);
		virtualDetectorLogical->SetSensitiveDetector(virtualDetectorSD);

		G4RotationMatrix *pRotVD=new G4RotationMatrix();
		pRotVD->rotateY(-mVDRotYAngle); 
		pRotVD->rotateX(-mVDRotXAngle); 
		G4ThreeVector pV3VDPos(0.0,0.0,mPivot2VDFace+mVirtualDetectorThick/2.);
		pV3VDPos.transform(pRotVD->inverse());

		//place a VD then the VB
		new G4PVPlacement(pRotVD,G4ThreeVector(pV3VDPos.x()+mPivotXOffset,
			pV3VDPos.y()+mPivotYOffset,pV3VDPos.z()+mPivotZOffset),
			virtualDetectorLogical,"virtualDetectorPhys",magneticLogical,0,0);

		pV3VDPos.set(0.0,0.0,mPivot2VDFace+1.5*mVirtualDetectorThick);
		pV3VDPos.transform(pRotVD->inverse());
		new G4PVPlacement(pRotVD,G4ThreeVector(pV3VDPos.x()+mPivotXOffset,
			pV3VDPos.y()+mPivotYOffset,pV3VDPos.z()+mPivotZOffset),
			virtualDetectorLogical,mVDPhysVolName,magneticLogical,0,0);
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
	G4VPhysicalVolume* theHRSPhys=0;

	G4RotationMatrix *pRotX90deg=new G4RotationMatrix();
	pRotX90deg->rotateX(-90*deg); 

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


	//////////////////////////
	//Q1 Entrance Window
	//////////////////////////
	bool pSetupQ1EntranceWindow=true;
	//G2P and CREX have their own VB defined at a different location
	if(mSetupG2PGeometry || mSetupCREXGeometry)  pSetupQ1EntranceWindow=false;
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
	double zPlane_Q1Q2Tunnel[] = {pHRSContainerRin+6*mm,pHallCenter2RQ1Face+pQ1Length+10.0*cm,
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

