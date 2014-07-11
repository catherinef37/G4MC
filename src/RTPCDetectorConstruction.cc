// ********************************************************************
//
// $Id: RTPCDetectorConstruction.cc,v 1.00, 2013/10/06 HRS Exp $
// --------------------------------------------------------------
//
// ********************************************************************
//
#include "RTPCDetectorConstruction.hh"

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
#include "HRSDCSD.hh"
#include "UsageManager.hh"
#include "HRSMaterial.hh"


extern UsageManager* gConfig;

////////////////////////////////////////////////////////////////////////////////////////////////
RTPCDetectorConstruction::RTPCDetectorConstruction(G4LogicalVolume *mother): 
mMotherLogVol(mother) 
{
	GetConfig();
	ConstructMaterials();
	
	G4cout<<"Contrstruct RTPC geometry ... done! "<<G4endl;
}

RTPCDetectorConstruction::~RTPCDetectorConstruction()
{
	//I might need to delete the materials
	//But it does not hurt if I don't, since this class will have only one instance
	//in the program
	G4cout<<"Delete RTPC geometry ... done! "<<G4endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////
void RTPCDetectorConstruction::GetConfig()
{
	gConfig->ReadFile("Detector_RTPC.ini");
	//////////////////////////////////////////////////////////////////

	const double mH2GasD_STP=0.08988*mg/cm3;
	const double mD2GasD_STP=0.180*mg/cm3;
	const double mHe3GasD_STP=0.1777*mg/cm3;
	const double mHeGasD_STP=0.1786*mg/cm3;
	const double mDMEGasD_STP=1.97*mg/cm3;

	gConfig->GetParameter("TargetXOffset",mTargetXOffset);
	mTargetXOffset*=mm;
	gConfig->GetParameter("TargetYOffset",mTargetYOffset);
	mTargetYOffset*=mm;
	gConfig->GetParameter("TargetZOffset",mTargetZOffset);
	mTargetZOffset*=mm;
	
	gConfig->GetParameter("RTPCLength",mRTPCLength);mRTPCLength*=mm;

	gConfig->GetParameter("D2GasL",mD2GasL);mD2GasL*=mm;
	gConfig->GetParameter("D2GasR",mD2GasR);mD2GasR*=mm;

	gConfig->GetParameter("TargetType",mTargetType);
	gConfig->GetParameter("D2GasT",mD2GasT);mD2GasT*=kelvin;
	gConfig->GetParameter("D2GasP",mD2GasP);mD2GasP*=atmosphere;
	double pTgGasD_STP = mD2GasD_STP;
	if(mTargetType==1) pTgGasD_STP = mH2GasD_STP;
	else if(mTargetType==3) pTgGasD_STP = mHe3GasD_STP;
	else if(mTargetType==4) pTgGasD_STP = mHeGasD_STP;
	mD2GasD=pTgGasD_STP*mD2GasP/atmosphere*273.15/mD2GasT;

	gConfig->GetParameter("HeGasT",mHeGasT);mHeGasT*=kelvin;
	gConfig->GetParameter("HeGasP",mHeGasP);mHeGasP*=atmosphere;
	mHeGasD=mHeGasD_STP*mHeGasP/atmosphere*273.15/mHeGasT;
	
	gConfig->GetParameter("MixGasT",mMixGasT);mMixGasT*=kelvin;
	gConfig->GetParameter("MixGasP",mMixGasP);mMixGasP*=atmosphere;
	mMixHeGasD=mHeGasD_STP*mMixGasP/atmosphere*273.15/mMixGasT;
	mMixDMEGasD=mDMEGasD_STP*mMixGasP/atmosphere*273.15/mMixGasT;

	mRatioHe2DME=82./18.;
	gConfig->GetParameter("RatioHe2DME",mRatioHe2DME);
	mMixDMEHeD=(mMixHeGasD*mRatioHe2DME+mMixDMEGasD)/(mRatioHe2DME+1.0);

	gConfig->GetParameter("TgWallMaterialType",mTgWallMaterialType);
	gConfig->GetParameter("TgWallThick",mTgWallThick);mTgWallThick*=mm;
	
	gConfig->GetParameter("1stMylarR",m1stMylarR);m1stMylarR*=mm;
	gConfig->GetParameter("1stAlThick",m1stAlThick);m1stAlThick*=mm;
	gConfig->GetParameter("1stMylarThick",m1stMylarThick);m1stMylarThick*=mm;

	gConfig->GetParameter("2ndMylarR",m2ndMylarR);m2ndMylarR*=mm;
	gConfig->GetParameter("2ndAlThick",m2ndAlThick);m2ndAlThick*=mm;
	gConfig->GetParameter("2ndMylarThick",m2ndMylarThick);m2ndMylarThick*=mm;

	
	gConfig->GetParameter("GEM1R",mGEM1R);mGEM1R*=mm;
	gConfig->GetParameter("GEM2R",mGEM2R);mGEM2R*=mm;
	gConfig->GetParameter("GEM3R",mGEM3R);mGEM3R*=mm;
	gConfig->GetParameter("PCBReadOutR",mPCBReadOutR);mPCBReadOutR*=mm;

	gConfig->GetParameter("SetupSolenoid",mSetupSolenoid);
	gConfig->GetParameter("SolenoidPosX",mSolenoidPosX);mSolenoidPosX*=mm;
	gConfig->GetParameter("SolenoidPosY",mSolenoidPosY);mSolenoidPosY*=mm;
	gConfig->GetParameter("SolenoidPosZ",mSolenoidPosZ);mSolenoidPosZ*=mm;
	
	gConfig->GetParameter("BStepLimit",mBStepLimit);mBStepLimit*=mm;
	gConfig->GetParameter("DCStepLimit",mDCStepLimit);mDCStepLimit*=mm;

	gConfig->GetParameter("BedPlateThick",mBedPlateThick);mBedPlateThick*=mm;
	gConfig->GetParameter("BedPlateLowEdge",mBedPlateLowEdge);mBedPlateLowEdge*=mm;
	gConfig->GetParameter("BedPlateHighEdge",mBedPlateHighEdge);mBedPlateHighEdge*=mm;

	gConfig->GetParameter("InnerGapSpThick",mInnerGapSpThick);mInnerGapSpThick*=mm;
	gConfig->GetParameter("G10FR4Thick",mG10FR4Thick);mG10FR4Thick*=mm;
	gConfig->GetParameter("GEM1SpThick",mGEM1SpThick);mGEM1SpThick*=mm;
	gConfig->GetParameter("GEM2SpThick",mGEM2SpThick);mGEM2SpThick*=mm;
	gConfig->GetParameter("GEM3SpThick",mGEM3SpThick);mGEM3SpThick*=mm;
	gConfig->GetParameter("ReadOutSpThick",mReadOutSpThick);mReadOutSpThick*=mm;

	gConfig->GetParameter("SetupEntranceNExitCap",mSetupEntranceNExitCap);
	gConfig->GetParameter("SetupEndPlateNCover",mSetupEndPlateNCover);
	gConfig->GetParameter("SetupCableNChip",mSetupCableNChip);

	//sanity check for the mBedPlateHighEdge, which determine the 
	//maximum radius of the whole RTPC
	if(mBedPlateHighEdge<=mPCBReadOutR+15*mm) 
	{
		cout<<"Warning: RTPC BedPlateHighEdge is lower than RTPC read out pad\n"
			<<"         Please check your ini file\n";
		mBedPlateHighEdge=mPCBReadOutR+15*mm;
	}
	if(mD2GasL<=mRTPCLength+75*mm) 
	{
		cout<<"Warning: RTPC target length is shorter than RTPC pad length\n"
			<<"         Please check your ini file\n";
		mD2GasL=mRTPCLength+75*mm;
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////
void RTPCDetectorConstruction::ConstructMaterials()
{
	//This is a stand alone version, most of this materials exist in HRSMaterial 
	G4double a;
	G4double z;
	G4double density;
	G4double weightRatio;
	G4String name;
	G4String symbol;
	G4int nElem,nComponent,nAtoms;
	G4double pressure;
	G4double temperature;

	// elements for mixtures and compounds
	a = 1.01*g/mole;
	G4Element* elH = new G4Element(name="Hydrogen", symbol="H", z=1., a);
	a = 12.01*g/mole;
	G4Element* elC = new G4Element(name="Carbon", symbol="C", z=6., a);
	a = 14.01*g/mole;
	G4Element* elN = new G4Element(name="Nitrogen", symbol="N", z=7., a);
	a = 16.00*g/mole;
	G4Element* elO = new G4Element(name="Oxigen", symbol="O", z=8., a);
	a = 28.09*g/mole;
	G4Element* elSi= new G4Element(name="Silicon", symbol="Si", z=14., a);
	a = 35.45*g/mole;
	G4Element* elCl= new G4Element(name="Chlorine", symbol="Cl", z=17., a);
	a = 51.9961*g/mole;
	G4Element* elCr= new G4Element(name="Chromium", symbol="Cr", z=24., a);
	a = 55.845*g/mole;
	G4Element* elFe= new G4Element(name="Iron", symbol="Fe", z=26., a);
	a = 58.6934*g/mole;
	G4Element* elNi= new G4Element(name="Nickel", symbol="Ni", z=28., a);
	a = 63.55*g/mole;
	G4Element* elCu= new G4Element(name="Copper", symbol="Cu", z=29., a);
	a = 95.94*g/mole;
	G4Element* elMo = new G4Element(name="Molydb", symbol="Mo", z=42., a);

	// Vaccum
	density = universe_mean_density;
	pressure=1.0e-19 *pascal;
	temperature=0.1 *kelvin;
	a=1.01 *g/mole;
	vaccum = new G4Material(name="Galactic",z=1,a,density,
		kStateGas, temperature, pressure);

	// Air
	density = 1.29*mg/cm3;
	air = new G4Material(name="Air", density, nElem=2);
	air->AddElement(elN, weightRatio=.7);
	air->AddElement(elO, weightRatio=.3);

	//deuteriumGas at 7.5 atm & 300K (room temperature);density = 0.12273*mg/cm3;
	density = mD2GasD;
	a=2.014*g/mole;
	pressure=mD2GasP;
	temperature=mD2GasT;
	deuteriumGas = new G4Material(name="DeuteriumGas", z=1., a, density,
		kStateGas, temperature, pressure);

	//kapton
	density = 1.42*g/cm3;
	kapton = new G4Material(name="Kapton", density, nElem=4);
	kapton->AddElement(elH, nAtoms=10);
	kapton->AddElement(elC, nAtoms=22);
	kapton->AddElement(elO, nAtoms=5);
	kapton->AddElement(elN, nAtoms=2);

	/////////////////////////////////////////////////////////
	//ultem C37H24O6N2
	density = 1.27*g/cm3;
	ultem = new G4Material(name="Ultem", density, nElem=4);
	ultem->AddElement(elH, nAtoms=37);
	ultem->AddElement(elC, nAtoms=24);
	ultem->AddElement(elO, nAtoms=6);
	ultem->AddElement(elN, nAtoms=2);

	//fused quartz SiO2
	density = 2.20*g/cm3;
	SiO2 = new G4Material(name="Fused-Quartz", density, nElem=2);
	SiO2->AddElement(elSi, nAtoms=1);
	SiO2->AddElement(elO, nAtoms=2);

	//epoxy-resin , C11H12O3
	//density = 0.95*g/cm3;     //volumn ratio 60%:40%
	density = 1.268*g/cm3;  //weight ratio 60%:40%
	epoxy = new G4Material(name="Epoxy-Resin", density, nElem=3);
	epoxy->AddElement(elH, nAtoms=12);
	epoxy->AddElement(elC, nAtoms=11);
	epoxy->AddElement(elO, nAtoms=3);

	//G10FR4 is 60% of SiO2 and 40% epoxy
	density = 1.7*g/cm3;
	G10FR4 = new G4Material(name="G10FR4", density, nComponent=2);
	G10FR4->AddMaterial(SiO2, weightRatio=0.6);
	G10FR4->AddMaterial(epoxy, weightRatio=0.4);

	//Rohacel71 0.75g/cm^3 ,H-[CH2]n-H
	density = 0.75*g/cm3;
	rohacel71 = new G4Material(name="Rohacel71", density, nElem=2);
	rohacel71->AddElement(elH, nAtoms=2);
	rohacel71->AddElement(elC, nAtoms=1);

	//pcbNchip Cu9Si5C460H506O138 epoxy + silicon + copper
	//weight ratio 57.7:4:1
	density = 2.31*g/cm3;
	pcbNchip = new G4Material(name="Pcb&Chip", density, nElem=5);
	pcbNchip->AddElement(elH, nAtoms=506);
	pcbNchip->AddElement(elC, nAtoms=460);
	pcbNchip->AddElement(elO, nAtoms=138);
	pcbNchip->AddElement(elSi, nAtoms=5);
	pcbNchip->AddElement(elCu, nAtoms=9);

	//cable Cu4C6H9Cl3
	density = 30.27/(2.15*0.04*100.) *g/cm3;
	cable = new G4Material(name="Cable", density, nElem=4);
	cable->AddElement(elH, nAtoms=9);
	cable->AddElement(elC, nAtoms=6);
	cable->AddElement(elCl, nAtoms=3);
	cable->AddElement(elCu, nAtoms=4);

	/////////////////////////////////////////////////////////

	//Target Gas 1=H2, 2=D2,3=He3,4=He4

	a = 1.0081*g/mole;
	H2TgGas = new G4Material(name="H2TgGas", z=1., a, density=mD2GasD,
		kStateGas, temperature=mD2GasT, pressure=mD2GasP);
	a = 2.01410178*g/mole;
	D2TgGas = new G4Material(name="D2TgGas", z=1., a, density=mD2GasD,
		kStateGas, temperature=mD2GasT, pressure=mD2GasP);
	a = 3.0160293*g/mole;
	He3TgGas = new G4Material(name="He3TgGas", z=2., a, density=mD2GasD,
		kStateGas, temperature=mD2GasT, pressure=mD2GasP);
	a = 4.0026*g/mole;
	He4TgGas = new G4Material(name="He4TgGas", z=2., a, density=mD2GasD,
		kStateGas, temperature=mD2GasT, pressure=mD2GasP);
	
	//Helium Gas at 1 atm & room temprature, density=0.163 *mg/cm3;
	a = 4.0026*g/mole;
	density = mHeGasD ;
	pressure= mHeGasP ;
	temperature = mHeGasT ;
	heliumGas = new G4Material(name="HeliumGas", z=2., a, density,
		kStateGas, temperature, pressure);
	//mylar
	density = 1.39*g/cm3;
	mylar = new G4Material(name="Mylar", density, nElem=2);
	mylar->AddElement(elH, nAtoms=8);
	mylar->AddElement(elC, nAtoms=10);

	//Aluminum
	a = 26.982*g/mole;
	density=2.70*g/cm3;
	aluminum = new G4Material(name="Aluminum", z=13., a, density);

	//StainlessSteel, Material Names : stainless steel
	//Material : Fe-Cr-Ni-Mo,  93xx: 3.25% Ni,1.2%Cr,0.12%Mo,
	density = 7.85 *g/cm3;
	stainlesssteel = new G4Material(name="StainlessSteel", density, nElem=4);
	stainlesssteel->AddElement(elFe, weightRatio=0.9543);
	stainlesssteel->AddElement(elCr, weightRatio=0.0325);
	stainlesssteel->AddElement(elNi, weightRatio=0.0120);
	stainlesssteel->AddElement(elMo, weightRatio=0.0012);

	//Mix Helium Gas with DME at 1 atm & room temprature, density=0.163 *mg/cm3;
	a = 4.0026*g/mole;
	density = mMixHeGasD ;
	pressure= mMixGasP ;
	temperature = mMixGasT ;
	mixHeGas = new G4Material(name="MixHeliumGas", z=2., a, density,
		kStateGas, temperature, pressure);

	//DMEGas((CH3-O-CH3) at 1 atm & room temprature, density=1.871 *mg/cm3;
	density = mMixDMEGasD ;
	DMEGas = new G4Material(name="DMEGas", density, nElem=3,
		kStateGas, temperature, pressure);
	//DMEGas->AddElement(elC, nAtoms=2);
	//DMEGas->AddElement(elH, nAtoms=6);
	//DMEGas->AddElement(elO, nAtoms=1);
	DMEGas->AddElement(elH, weightRatio=0.13);
	DMEGas->AddElement(elC, weightRatio=0.522);
	DMEGas->AddElement(elO, weightRatio=0.348);

	//BONUS Gas at 1 atm & room temprature, 82% Helium gas and 18% DME gas
	double WeightRatioDME=46.07/(4.0026*mRatioHe2DME+46.07);
	density = mMixDMEHeD;
	bonusGas = new G4Material(name="BonusGas", density, nComponent=2,
		kStateGas, temperature, pressure);
	bonusGas->AddMaterial(DMEGas, weightRatio=WeightRatioDME);
	bonusGas->AddMaterial(mixHeGas, weightRatio=1.0-WeightRatioDME);

	//copper
	a = 63.546*g/mole;
	density = 8.96*g/cm3;
	copper = new G4Material(name="Copper", z=29., a, density);

}


////////////////////////////////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* RTPCDetectorConstruction::Construct()
{
	G4SDManager* SDMan = G4SDManager::GetSDMpointer();
	G4String SDname;
	
	G4VSensitiveDetector* SDTargetWall = new HRSStdSD(SDname="RTPCTargetWall");
	SDMan->AddNewDetector(SDTargetWall);

	G4VSensitiveDetector* SDMylar1 = new HRSStdSD(SDname="RTPCMylar1");
	SDMan->AddNewDetector(SDMylar1);
	G4VSensitiveDetector* SDMylar2 = new HRSStdSD(SDname="RTPCMylar2");
	SDMan->AddNewDetector(SDMylar2);

	G4VSensitiveDetector* SDDC = new HRSDCSD(SDname="RTPCDC");
	SDMan->AddNewDetector(SDDC);

	G4VSensitiveDetector* SDGEM1 = new HRSStdSD(SDname="RTPCGEM1");
	SDMan->AddNewDetector(SDGEM1);
	G4VSensitiveDetector* SDGEM2 = new HRSStdSD(SDname="RTPCGEM2");
	SDMan->AddNewDetector(SDGEM2);
	G4VSensitiveDetector* SDGEM3 = new HRSStdSD(SDname="RTPCGEM3");
	SDMan->AddNewDetector(SDGEM3);

	G4VSensitiveDetector* SDReadOut = new HRSStdSD(SDname="RTPCReadOut");
	SDMan->AddNewDetector(SDReadOut);
	G4VSensitiveDetector* SDCableNChip = new HRSStdSD(SDname="CableNChip");
	SDMan->AddNewDetector(SDCableNChip);
	
	G4VSensitiveDetector* virtualBoundarySD=new HRSStdSD("virtualBoundary");
	//############################################

	//////////////////////////
	//The mother box
	//////////////////////////

	//In this way we can place everything into the box as if placed them in the hall 
	// Size of this tub, large enough but not too big
	double RTPCContainerL = mD2GasL+120*mm;       
	double RTPCContainerRin = 0.0*m; 
	double RTPCContainerRout = max(mPCBReadOutR+50*mm,mBedPlateHighEdge+12*mm); 

	if(mSetupSolenoid==1)  {RTPCContainerRout=0.50*m; RTPCContainerL=1.2*m;}
	else if(mSetupSolenoid>=2)  {RTPCContainerRout=1.00*m; RTPCContainerL=1.6*m;} 

	G4Tubs* RTPCContainer = new G4Tubs("RTPCContainer",RTPCContainerRin,
		RTPCContainerRout,RTPCContainerL/2,0.0*deg,360.*deg);

	G4LogicalVolume* RTPCContainerLogical = new G4LogicalVolume(RTPCContainer, 
		air, "logRTPCContainer", 0, 0, 0);
	RTPCContainerLogical->SetVisAttributes(HallVisAtt);

	//the position at the hall
	double RTPCContainerPosX=mTargetXOffset;
	double RTPCContainerPosY=mTargetYOffset;
	double RTPCContainerPosZ=mTargetZOffset;
	G4PVPlacement* phyRTPCContainer= new G4PVPlacement(0,
		G4ThreeVector(RTPCContainerPosX,RTPCContainerPosY,RTPCContainerPosZ),
		RTPCContainerLogical,"RTPCContainner",mMotherLogVol,0,0);

	
	//////////////////////////
	//The RTPC ......//prepare some variables
	//////////////////////////
	G4double Z_Half=mRTPCLength/2;		//8.768-0.25-0.25 inch
	G4double R20=m1stMylarR-m1stMylarThick-2.*m1stAlThick;	//19.98373 mm
	G4double R30=m2ndMylarR-m2ndMylarThick-2.*m2ndAlThick; //29.9974 *mm; //1.181*inch
	

	G4double phi_min,phi_max; //start and end phi angle
	G4double r_min,r_max;
	G4double x_max,x_min,y_max,y_min; //start and end radius
	G4double xx_h,yy_h,zz_h;
	G4double z_down,z_up;

	G4RotationMatrix* RotZ90deg = new G4RotationMatrix();
	RotZ90deg->rotateZ(90.*deg);
	G4RotationMatrix* RotZ270deg = new G4RotationMatrix();
	RotZ270deg->rotateZ(270.*deg);	
	G4RotationMatrix* RotZero=new G4RotationMatrix();


	// set "user limits" for drawing smooth curve
	G4UserLimits* uDCStepLimits = new G4UserLimits(mDCStepLimit);
	G4UserLimits* uBStepLimits = new G4UserLimits(mBStepLimit);

	
	//////////////////////////
	//The solenoid
	//////////////////////////
	if(mSetupSolenoid)
	{
		ConstructSolenoid(RTPCContainerLogical);
	}


	/////////////////////
	//Target 
	/////////////////////
	if(mTargetType==1) targetMaterial=H2TgGas;
	else if(mTargetType==3) targetMaterial=He3TgGas;
	else if(mTargetType==4) targetMaterial=He4TgGas;
	else targetMaterial=D2TgGas;
	G4VSolid* targetVesselSolid = new G4Tubs("targetVesselTubs",
		0.,mD2GasR,mD2GasL/2.0,0.,360.*deg);
	G4LogicalVolume* targetVesselLogical = new G4LogicalVolume(targetVesselSolid,
		targetMaterial,"targetVesselLogical",0,0,0);
	new G4PVPlacement(0,G4ThreeVector(),targetVesselLogical,
		"targetPhys",RTPCContainerLogical,0,0);

	/////////////////////
	//Target wall cylinder with 50 um Kapton wall
	/////////////////////
	//mTgWallMaterialType;  //1 is kapton, 2 is aluminum
	if(mTgWallMaterialType==2) targetWallMaterial=aluminum;
	else targetWallMaterial=kapton;
	G4VSolid* targetWallSolid = new G4Tubs("targetWallTubs",mD2GasR,
		mTgWallThick+mD2GasR,mD2GasL/2.0,0.,360.*deg);

	G4LogicalVolume* targetWallLogical = new G4LogicalVolume(targetWallSolid,
		targetWallMaterial,"targetWallLogical",0,0,0);
	
	new G4PVPlacement(0,G4ThreeVector(),targetWallLogical,
		"targetWallPhys",RTPCContainerLogical,0,0);

	targetWallLogical->SetSensitiveDetector(SDTargetWall);
	if(mTgWallMaterialType==2)
		targetWallLogical->SetVisAttributes(GrayVisAtt);  //grey 
	else 
		targetWallLogical->SetVisAttributes(DarkRedVisAtt);  //dark red

	//setup entrance cap and exit cap
	G4VSolid* entranceCoverSolid=0;
	if(mSetupEntranceNExitCap)
	{
		//////////////////////////
		//Downtream exit window aluminum cap 
		//////////////////////////
		//a cup of 3.1 mm height, including 0.1 mm thick bottom, 0.5
		phi_min=0.*deg;
		phi_max=360.*deg;
		const G4double  zPlane_cap[] ={-3*mm+mD2GasL/2,mD2GasL/2,mD2GasL/2,mD2GasL/2+0.1*mm};
		const G4double  rInner_cap[] ={mD2GasR+mTgWallThick,mD2GasR+mTgWallThick,0.0*mm,0.0*mm};
		const G4double rOutner_cap[] ={mD2GasR+0.6*mm,mD2GasR+0.6*mm,mD2GasR+0.6*mm,mD2GasR+0.6*mm};
		G4VSolid* exitCapSolid = new G4Polycone("exitCapPcon",
			phi_min,phi_max,4,zPlane_cap,rInner_cap,rOutner_cap);
		G4LogicalVolume* exitCapLogical = new G4LogicalVolume(exitCapSolid,
			aluminum,"exitCapLogical",0,0,0);
		new G4PVPlacement(0,G4ThreeVector(),exitCapLogical,
			"exitCapPhys",RTPCContainerLogical,0,0);
		exitCapLogical->SetVisAttributes(WhiteVisAtt);

		//////////////////////////
		//Uptream entrance aluminum cover  
		//////////////////////////
		//const G4double  zPlane_u[] ={-58.70*mm,-70.70*mm,-70.70*mm,-72.70*mm,-72.70*mm,
		//	-74.70*mm,-79.45*mm,-79.45*mm,-113.8936*mm,-113.8936*mm,
		//	-123.8936*mm,-123.8936*mm,-145.0*mm,-145.0*mm,-145.1*mm};
		//const G4double  rInner_u[] ={3.125*mm,3.125*mm,3.050*mm,3.050*mm,3.500*mm,
		//	3.500*mm,3.500*mm,3.500*mm,3.500*mm,3.500*mm,
		//	3.500*mm,3.500*mm,3.500*mm,0.000*mm.0.000*mm};
		//const G4double rOutner_u[] ={3.150*mm,7.726*mm,7.726*mm,8.493*mm,8.493*mm,
		//	9.260*mm,9.260*mm,12.00*mm,12.00*mm,15.00*mm,
		//	15.00*mm,12.00*mm,12.00*mm,12.00*mm,12.00*mm};

		phi_min=0.*deg;
		phi_max=360.*deg;
		z_down=-mRTPCLength/2+51.3*mm;  //Its down stream edge is 51.3*mm from RTPC up edge
		r_min=mD2GasR+0.5*mm;  //3.5*mm
		const G4double  zPlane_u[] ={
			z_down,z_down-12*mm,z_down-12*mm,z_down-14*mm,z_down-14*mm,
			z_down-16*mm,z_down-21.75*mm,z_down-21.75*mm,z_down-55.1936*mm,z_down-55.1936*mm,
			z_down-65.1936*mm,z_down-65.1936*mm,z_down-86.3*mm,z_down-86.3*mm,z_down-86.4*mm};
		const G4double  rInner_u[] ={
			r_min-0.375*mm,r_min-0.375*mm,r_min-0.450*mm,r_min-0.450*mm,r_min,
			r_min,r_min,r_min,r_min,r_min,
			r_min,r_min,r_min,0,0};
		const G4double rOutner_u[] ={
			r_min-0.375*mm,r_min+3.226*mm,r_min+3.226*mm,r_min+4.993*mm,r_min+4.993*mm,
			r_min+5.760*mm,r_min+5.760*mm,r_min+8.50*mm,r_min+8.50*mm,r_min+11.50*mm,
			r_min+11.50*mm,r_min+8.50*mm,r_min+8.50*mm,r_min+8.50*mm,r_min+8.50*mm};
		entranceCoverSolid = new G4Polycone("entranceCoverPcon",
			phi_min,phi_max,15,zPlane_u,rInner_u,rOutner_u);
		G4LogicalVolume* entranceCoverLogical = new G4LogicalVolume(entranceCoverSolid,
			aluminum,"entranceCoverLogical",0,0,0);
		new G4PVPlacement(0,G4ThreeVector(),entranceCoverLogical,
			"entranceCoverPhys",RTPCContainerLogical,0,0);

		entranceCoverLogical->SetVisAttributes(WhiteVisAtt); //white
	}


	/////////////////////
	//*****first gap //1 atm helium gas at 2.14mm<r<20mm
	/////////////////////
	//Need to subtract the entrance cover to avoid overlapping
	G4VSolid* firstGapSolid = 0;
	G4VSolid* firstGapSolid_whole = new G4Tubs("firstGapTubs",
		mTgWallThick+mD2GasR,R20,Z_Half,0.,360.*deg);
	if(entranceCoverSolid)
	{
		firstGapSolid = new G4SubtractionSolid("firstGapSolid",
			firstGapSolid_whole,entranceCoverSolid);
	}
	else
	{
		firstGapSolid = firstGapSolid_whole;
	}
	G4LogicalVolume* firstGapLogical = new G4LogicalVolume(firstGapSolid,
		heliumGas,"firstGapLogical",0,0,0);
	new G4PVPlacement(0,G4ThreeVector(),firstGapLogical,
		"firstGapPhys",RTPCContainerLogical,0,0);
	firstGapLogical->SetVisAttributes(HallVisAtt);



	/////////////////////
	//fisrt foils
	/////////////////////
	//*****first cylider foil//0.006mm Mylar and 0.000035mm Al on both sides at r=20mm
	//Thia want to test how low the meomentum threshold can go with this 1at mylar layer
	bool bSetup1stMylar=true;
	if(bSetup1stMylar)
	{
		r_min=R20;	//19.98373
		r_max=R20+m1stAlThick;	//19.98373+0.000035

		phi_min=asin(0.5*mBedPlateThick/r_min);  //10.315*deg;
		phi_max=180.0*deg-2.0*phi_min;
		G4VSolid* firstAlFoil1Solid
			= new G4Tubs("firstAlFoil1Tubs",r_min,r_max,Z_Half,phi_min,phi_max);
		G4LogicalVolume* firstAlFoil1Logical
			= new G4LogicalVolume(firstAlFoil1Solid,aluminum,"firstAlFoil1Logical",0,0,0);
		//left half
		new G4PVPlacement(RotZ90deg,G4ThreeVector(),firstAlFoil1Logical,
			"firstAlFoil1Phys",RTPCContainerLogical,0,0);
		//right half
		new G4PVPlacement(RotZ270deg,G4ThreeVector(),firstAlFoil1Logical,
			"firstAlFoil1Phys",RTPCContainerLogical,0,0);

		r_min=R20+m1stAlThick;	//19.98373+0.000035
		r_max=R20+m1stMylarThick+m1stAlThick; //19.98373+0.000035+0.006
		phi_min=asin(0.5*mBedPlateThick/r_min);    //10.315 *deg;
		phi_max=180.0*deg-2.0*phi_min;
		G4VSolid* firstMylarSolid
			= new G4Tubs("firstMylarTubs",r_min,r_max,Z_Half,phi_min,phi_max);
		G4LogicalVolume* firstMylarLogical
			= new G4LogicalVolume(firstMylarSolid,mylar,"firstMylarLogical",0,0,0);
		//left half
		new G4PVPlacement(RotZ90deg,G4ThreeVector(),firstMylarLogical,
			"firstMylarPhys",RTPCContainerLogical,0,0);
		//right half
		new G4PVPlacement(RotZ270deg,G4ThreeVector(),firstMylarLogical,
			"firstMylarPhys",RTPCContainerLogical,0,0);

		r_min=R20+m1stMylarThick+m1stAlThick;
		r_max=R20+m1stMylarThick+2.0*m1stAlThick;
		phi_min=asin(0.5*mBedPlateThick/r_min); //10.315 *deg;
		phi_max=180.0*deg-2.0*phi_min;
		G4VSolid* firstAlFoil2Solid
			= new G4Tubs("firstAlFoil2Tubs",r_min,r_max,Z_Half,phi_min,phi_max);
		G4LogicalVolume* firstAlFoil2Logical
			= new G4LogicalVolume(firstAlFoil2Solid,aluminum,"firstAlFoil2Logical",0,0,0);
		//left half
		new G4PVPlacement(RotZ90deg,G4ThreeVector(),firstAlFoil2Logical,
			"firstAlFoil2Phys",RTPCContainerLogical,0,0);
		//right half
		new G4PVPlacement(RotZ270deg,G4ThreeVector(),firstAlFoil2Logical,
			"firstAlFoil2Phys",RTPCContainerLogical,0,0);

		firstMylarLogical->SetSensitiveDetector(SDMylar2);
		firstAlFoil1Logical->SetVisAttributes(WhiteVisAtt);   //white
		firstMylarLogical->SetVisAttributes(GrayVisAtt); //grey
		firstAlFoil2Logical->SetVisAttributes(WhiteVisAtt);   //white
	}

	
	/////////////////////
	//Bed plate  	
	/////////////////////
	//*****Bed plate, from r=20 to r=108mm, BedPlateThick=7.1628mm
	//the Bed is 8 mm protured
	double protruedHeight_h=4.0*mm;
	x_min=mBedPlateLowEdge;
	x_max=mBedPlateHighEdge;
	xx_h=(x_max-x_min)/2.0;
	yy_h=mBedPlateThick/2.0;
	zz_h=Z_Half;
	G4ThreeVector Tran_bed(-1.*(xx_h+protruedHeight_h), 0., 0.);
	//G4VSolid* bedPlateSolid = new G4Box("bedPlateBox",xx_h,yy_h,zz_h);
	//end plate thickness 6.35mm
	G4VSolid* Box_bed_edge = new G4Box("Box_bed_edge",protruedHeight_h,yy_h,Z_Half+6.35*mm);
	G4VSolid* Box_bed = new G4Box("Box_bed",xx_h,yy_h,zz_h);
	G4UnionSolid* bedPlateSolid = new G4UnionSolid("bedPlateSolid",
		Box_bed ,Box_bed_edge,RotZero,Tran_bed);
	G4LogicalVolume* bedPlateLogical = new G4LogicalVolume(bedPlateSolid,
		ultem,"bedPlateLogical",0,0,0);
	//left half
	new G4PVPlacement(RotZ90deg,G4ThreeVector(0.,x_min+xx_h,0.),
		bedPlateLogical,"bedPlatePhysL",RTPCContainerLogical,0,0);
	//right half
	new G4PVPlacement(RotZ270deg,G4ThreeVector(0.,-x_min-xx_h,0.),
		bedPlateLogical,"bedPlatePhysR",RTPCContainerLogical,0,0);

	bedPlateLogical->SetVisAttributes(YellowVisAtt);


	///////////////////// 
	//PCBBedPlate
	/////////////////////
	//PCB board attach to the bed plate
	x_min=R30+2.*m2ndAlThick+m2ndMylarThick;
	//in case the bed plate does not start from the 2nd marler foil
	if(x_min<mBedPlateLowEdge) x_min=mBedPlateLowEdge;

	x_max=mBedPlateHighEdge;
	y_min=mBedPlateThick/2.0;
	y_max=mG10FR4Thick+mBedPlateThick/2.0;
	xx_h=(x_max-x_min)/2.0;
	yy_h=mG10FR4Thick/2.0;
	zz_h=Z_Half;
	G4ThreeVector Tran_pcbBedSp(-1.*(xx_h+protruedHeight_h), 0., 0.);
	G4VSolid* Box_pcbBedSp = new G4Box("Box_pcbBedSp",xx_h,yy_h,zz_h);
	G4VSolid* Box_pcbBedSp_edge = new G4Box("Box_pcbBedSp_edge",
		protruedHeight_h,yy_h,Z_Half+6.35*mm);
	G4UnionSolid* pcbBedSpSolid = new G4UnionSolid("pcbBedSpBox",
		Box_pcbBedSp,Box_pcbBedSp_edge,RotZero,Tran_pcbBedSp);
	G4LogicalVolume* pcbBedSpLogical
		= new G4LogicalVolume(pcbBedSpSolid,G10FR4,"pcbBedSpLogical",0,0,0);
	//left half
	new G4PVPlacement(RotZ270deg,G4ThreeVector(y_min+yy_h,-x_min-xx_h,0.),
		pcbBedSpLogical,"pcbBedSpPhysL1",RTPCContainerLogical,0,0);
	new G4PVPlacement(RotZ90deg,G4ThreeVector(y_min+yy_h,x_min+xx_h,0.),
		pcbBedSpLogical,"pcbBedSpPhysR1",RTPCContainerLogical,0,0);
	//right half
	new G4PVPlacement(RotZ270deg,G4ThreeVector(-y_min-yy_h,-x_min-xx_h,0.),
		pcbBedSpLogical,"pcbBedSpPhysL2",RTPCContainerLogical,0,0);
	new G4PVPlacement(RotZ90deg,G4ThreeVector(-y_min-yy_h,x_min+xx_h,0.),
		pcbBedSpLogical,"pcbBedSpPhysR2",RTPCContainerLogical,0,0);

	pcbBedSpLogical->SetVisAttributes(PcbGreenVisAtt);  


	///////////////////// 
	//InnerGap support
	///////////////////// 
	x_max=R30;
	x_min=R20;
	y_min=mBedPlateThick/2.0;
	y_max=(mBedPlateThick/2.0)+mInnerGapSpThick;
	xx_h=(x_max-x_min)/2.0;
	yy_h=mInnerGapSpThick/2.0;
	zz_h=Z_Half;
	G4VSolid* innerGapSpSolid = new G4Box("innerGapSpBox",xx_h,yy_h,zz_h);
	G4LogicalVolume* innerGapSpLogical
		= new G4LogicalVolume(innerGapSpSolid,ultem,"innerGapSpLogical",0,0,0);
	//left half
	new G4PVPlacement(RotZ90deg,G4ThreeVector(y_min+yy_h,-x_min-xx_h,0.),
		innerGapSpLogical, "innerGapSpPhysL1",RTPCContainerLogical,0,0);
	new G4PVPlacement(RotZ270deg,G4ThreeVector(-y_min-yy_h,-x_min-xx_h,0.),
		innerGapSpLogical, "innerGapSpPhysL2",RTPCContainerLogical,0,0);
	//right half
	new G4PVPlacement(RotZ90deg,G4ThreeVector(y_min+yy_h,x_min+xx_h,0.),
		innerGapSpLogical, "innerGapSpPhysR1",RTPCContainerLogical,0,0);
	new G4PVPlacement(RotZ270deg,G4ThreeVector(-y_min-yy_h,x_min+xx_h,0.),
		innerGapSpLogical, "innerGapSpPhysR2",RTPCContainerLogical,0,0);

	innerGapSpLogical->SetVisAttributes(YellowVisAtt); //ultem yellow

	///////////////////// 
	//InnerGap drift region
	///////////////////// 
	//second drift region //1 atm 300k mixture gases (bonusGas) at 20.08907mm<r<30mm
	r_min=R20+m1stMylarThick+2.0*m1stAlThick;
	r_max=R30;
	phi_min=asin((0.5*mBedPlateThick+mInnerGapSpThick)/r_max); //9.896 *deg;
	phi_max=180.0*deg-2.0*phi_min;
	G4VSolid* innerGapSolid
		= new G4Tubs("innerGapTubs",r_min,r_max,Z_Half,phi_min,phi_max);
	G4LogicalVolume* innerGapLogical
		= new G4LogicalVolume(innerGapSolid,bonusGas,"innerGapLogical",0,0,0);
	new G4PVPlacement(RotZ90deg,G4ThreeVector(),innerGapLogical,
		"innerGapPhys",RTPCContainerLogical,0,0);
	new G4PVPlacement(RotZ270deg,G4ThreeVector(),innerGapLogical,
		"innerGapPhys",RTPCContainerLogical,0,0);
	innerGapLogical->SetUserLimits(uBStepLimits);
	innerGapLogical->SetVisAttributes(HallVisAtt);

	///////////////////// 
	//2nd mylar foil
	///////////////////// 
	//second cylider foil//0.0254mm Mylar and 0.000035 Al on both sides at r=30mm
	r_min=R30;	//29.9974 *mm; //2.362*inch
	r_max=R30+m2ndAlThick;	//29.9974+0.006 mm;
	phi_min=asin((0.5*mBedPlateThick)/r_min);//6.856 *deg;
	phi_max=180.0*deg-2.0*phi_min;
	G4VSolid* secondAlFoil1Solid
		= new G4Tubs("secondAlFoil1Tubs",r_min,r_max,Z_Half,phi_min,phi_max);
	G4LogicalVolume* secondAlFoil1Logical
		= new G4LogicalVolume(secondAlFoil1Solid,aluminum,"secondAlFoil1Logical",0,0,0);
	new G4PVPlacement(RotZ90deg,G4ThreeVector(),secondAlFoil1Logical,
		"secondAlFoil1Phys",RTPCContainerLogical,0,0);
	new G4PVPlacement(RotZ270deg,G4ThreeVector(),secondAlFoil1Logical,
		"secondAlFoil1Phys",RTPCContainerLogical,0,0);

	r_min=R30+m2ndAlThick;	//29.9974+0.006
	r_max=R30+m2ndAlThick+m2ndMylarThick; //29.9974+0.006+0.000035
	phi_min=asin((0.5*mBedPlateThick)/r_min);//6.856 *deg;
	phi_max=180.0*deg-2.0*phi_min;
	G4VSolid* secondMylarSolid
		= new G4Tubs("secondMylarTubs",r_min,r_max,Z_Half,phi_min,phi_max);
	G4LogicalVolume* secondMylarLogical
		= new G4LogicalVolume(secondMylarSolid,mylar,"secondMylarLogical",0,0,0);
	new G4PVPlacement(RotZ90deg,G4ThreeVector(),secondMylarLogical,
		"secondMylarPhys",RTPCContainerLogical,0,0);
	new G4PVPlacement(RotZ270deg,G4ThreeVector(),secondMylarLogical,
		"secondMylarPhys",RTPCContainerLogical,0,0);

	r_min=R30+m2ndAlThick+m2ndMylarThick;
	r_max=R30+2.*m2ndAlThick+m2ndMylarThick;
	phi_min=asin((0.5*mBedPlateThick)/r_min);//6.856 *deg;
	phi_max=180.0*deg-2.0*phi_min;
	G4VSolid* secondAlFoil2Solid
		= new G4Tubs("secondAlFoil2Tubs",r_min,r_max,Z_Half,phi_min,phi_max);
	G4LogicalVolume* secondAlFoil2Logical
		= new G4LogicalVolume(secondAlFoil2Solid,aluminum,"secondAlFoil2Logical",0,0,0);
	new G4PVPlacement(RotZ90deg,G4ThreeVector(),secondAlFoil2Logical,
		"secondAlFoil2Phys",RTPCContainerLogical,0,0);
	new G4PVPlacement(RotZ270deg,G4ThreeVector(),secondAlFoil2Logical,
		"secondAlFoil2Phys",RTPCContainerLogical,0,0);

	secondMylarLogical->SetSensitiveDetector(SDMylar2);
	secondAlFoil1Logical->SetVisAttributes(WhiteVisAtt);  //white
	secondMylarLogical->SetVisAttributes(GrayVisAtt);//grey
	secondAlFoil2Logical->SetVisAttributes(WhiteVisAtt);  //white


	//////////////////////////
	//main Drift Region
	//////////////////////////
	//1 atm 300k mixture gases (bonusGas) at 30.02547mm<r<60mm
	r_min=R30+2.*m2ndAlThick+m2ndMylarThick;
	r_max=mGEM1R;
	phi_min=asin((0.5*mBedPlateThick+mG10FR4Thick)/r_min);//9.895 *deg;
	phi_max=180.0*deg-2.0*phi_min;
	G4VSolid* driftRegionSolid
		= new G4Tubs("driftRegionTubs",r_min,r_max,Z_Half,phi_min,phi_max);
	G4LogicalVolume* driftRegionLogical = new G4LogicalVolume(driftRegionSolid,
		bonusGas,"driftRegionLogical",0,0,uDCStepLimits);
	new G4PVPlacement(RotZ90deg,G4ThreeVector(),driftRegionLogical,
		"driftRegionPhys",RTPCContainerLogical,0,0);
	new G4PVPlacement(RotZ270deg,G4ThreeVector(),driftRegionLogical,
		"driftRegionPhys",RTPCContainerLogical,0,0);
	// Set additional contraints on the track, with G4UserSpecialCuts
	driftRegionLogical->SetUserLimits(uDCStepLimits);
	driftRegionLogical->SetSensitiveDetector(SDDC);
	driftRegionLogical->SetVisAttributes(HallVisAtt);
	//main Drift Region end  //main Drift Region end  //main Drift Region end

	//////////////////////////
	//GEM layer 1
	//////////////////////////
	//first GEM foil//0.05mm Kapton and 0.005 Cu on both sides at r=60mm
	r_min=mGEM1R;
	r_max=mGEM1R+0.005*mm;
	phi_min=asin((0.5*mBedPlateThick+mG10FR4Thick)/r_min);//4.930 *deg;
	phi_max=180.0*deg-2.0*phi_min;
	G4VSolid* GEM1CuFoil1Solid
		= new G4Tubs("GEM1CuFoil1Tubs",r_min,r_max,Z_Half,phi_min,phi_max);
	G4LogicalVolume* GEM1CuFoil1Logical
		= new G4LogicalVolume(GEM1CuFoil1Solid,copper,"GEM1CuFoil1Logical",0,0,0);
	new G4PVPlacement(RotZ90deg,G4ThreeVector(),GEM1CuFoil1Logical,
		"GEM1CuFoil1Phys",RTPCContainerLogical,0,0);
	new G4PVPlacement(RotZ270deg,G4ThreeVector(),GEM1CuFoil1Logical,
		"GEM1CuFoil1Phys",RTPCContainerLogical,0,0);

	r_min=mGEM1R+0.005*mm;
	r_max=mGEM1R+0.005*mm+0.05*mm;
	phi_min=asin((0.5*mBedPlateThick+mG10FR4Thick)/r_min);//4.930 *deg;
	phi_max=180.0*deg-2.0*phi_min;
	G4VSolid* GEM1KaptonSolid
		= new G4Tubs("GEM1KaptonTubs",r_min,r_max,Z_Half,phi_min,phi_max);
	G4LogicalVolume* GEM1KaptonLogical
		= new G4LogicalVolume(GEM1KaptonSolid,kapton,"GEM1KaptonLogical",0,0,0);
	new G4PVPlacement(RotZ90deg,G4ThreeVector(),GEM1KaptonLogical,
		"GEM1KaptonPhys",RTPCContainerLogical,0,0);
	new G4PVPlacement(RotZ270deg,G4ThreeVector(),GEM1KaptonLogical,
		"GEM1KaptonPhys",RTPCContainerLogical,0,0);

	r_min=mGEM1R+0.005*mm+0.05*mm;
	r_max=mGEM1R+0.005*mm+0.05*mm+0.005*mm;
	phi_min=asin((0.5*mBedPlateThick+mG10FR4Thick)/r_min);//4.930 *deg;
	phi_max=180.0*deg-2.0*phi_min;
	G4VSolid* GEM1CuFoil2Solid
		= new G4Tubs("GEM1CuFoil2Tubs",r_min,r_max,Z_Half,phi_min,phi_max);
	G4LogicalVolume* GEM1CuFoil2Logical
		= new G4LogicalVolume(GEM1CuFoil2Solid,copper,"GEM1CuFoil2Logical",0,0,0);
	new G4PVPlacement(RotZ90deg,G4ThreeVector(),GEM1CuFoil2Logical,
		"GEM1CuFoil2Phys",RTPCContainerLogical,0,0);
	new G4PVPlacement(RotZ270deg,G4ThreeVector(),GEM1CuFoil2Logical,
		"GEM1CuFoil2Phys",RTPCContainerLogical,0,0);

	GEM1CuFoil1Logical->SetSensitiveDetector(SDGEM1);
	GEM1KaptonLogical->SetSensitiveDetector(SDGEM1);
	GEM1CuFoil2Logical->SetSensitiveDetector(SDGEM1);
	GEM1CuFoil1Logical->SetVisAttributes(CuBrownVisAtt);
	GEM1KaptonLogical->SetVisAttributes(CuBrownVisAtt);
	GEM1CuFoil2Logical->SetVisAttributes(CuBrownVisAtt);

	//////////////////////////
	//GEM1 support
	//////////////////////////
	x_min=mGEM1R+0.005*mm+0.05*mm+0.005*mm;
	x_max=mBedPlateHighEdge;
	y_min=mBedPlateThick/2.0+mG10FR4Thick;
	xx_h=(x_max-x_min)/2.0;
	yy_h=mGEM1SpThick/2.0;
	zz_h=Z_Half;
	G4VSolid* GEM1SpSolid = new G4Box("GEM1SpBox",xx_h,yy_h,zz_h);
	G4LogicalVolume* GEM1SpLogical
		= new G4LogicalVolume(GEM1SpSolid,ultem,"GEM1SpLogical",0,0,0);
	//left half
	new G4PVPlacement(RotZ90deg,G4ThreeVector(y_min+yy_h,-x_min-xx_h,0.),
		GEM1SpLogical,"GEM1SpPhysL1",RTPCContainerLogical,0,0);
	new G4PVPlacement(RotZ270deg,G4ThreeVector(-y_min-yy_h,-x_min-xx_h,0.),
		GEM1SpLogical,"GEM1SpPhysL2",RTPCContainerLogical,0,0);
	//right half
	new G4PVPlacement(RotZ90deg,G4ThreeVector(y_min+yy_h,x_min+xx_h,0.),
		GEM1SpLogical,"GEM1SpPhysR1",RTPCContainerLogical,0,0);
	new G4PVPlacement(RotZ270deg,G4ThreeVector(-y_min-yy_h,x_min+xx_h,0.),
		GEM1SpLogical,"GEM1SpPhysR2",RTPCContainerLogical,0,0);

	GEM1SpLogical->SetVisAttributes(YellowVisAtt); //ultem yellow //color closed to Cu

	//////////////////////////
	//Drift Region between gem1 and gem2 
	//////////////////////////
	//1 atm 300k mixture gases (bonusGas) at 60.06mm<r<63mm
	r_min=mGEM1R+0.005*mm+0.05*mm+0.005*mm;
	r_max=mGEM2R;
	phi_min=asin((0.5*mBedPlateThick+mG10FR4Thick+mGEM1SpThick)/r_min);
	phi_max=180.0*deg-2.0*phi_min;
	G4VSolid* driftRegionG1G2Solid
		= new G4Tubs("driftRegionG1G2Tubs",r_min,r_max,Z_Half,phi_min,phi_max);
	G4LogicalVolume* driftRegionG1G2Logical = new G4LogicalVolume(driftRegionG1G2Solid,
		bonusGas,"driftRegionG1G2Logical",0,0,0);
	new G4PVPlacement(RotZ90deg,G4ThreeVector(),driftRegionG1G2Logical,
		"driftG1G2Phys",RTPCContainerLogical,0,0);
	new G4PVPlacement(RotZ270deg,G4ThreeVector(),driftRegionG1G2Logical,
		"driftG1G2Phys",RTPCContainerLogical,0,0);
	driftRegionG1G2Logical->SetVisAttributes(HallVisAtt);
	//Drift Region between gem1 and gem2 end  


	//////////////////////////
	//GEM layer 2
	//////////////////////////
	//second GEM foil//0.05mm Kapton and 0.005 Cu on both sides at r=63mm
	r_min=mGEM2R;
	r_max=mGEM2R+0.005*mm;
	phi_min=asin((0.5*mBedPlateThick+mG10FR4Thick+mGEM1SpThick)/r_min);//6.621 *deg;
	phi_max=180.0*deg-2.0*phi_min;
	G4VSolid* GEM2CuFoil1Solid
		= new G4Tubs("GEM2CuFoil1Tubs",r_min,r_max,Z_Half,phi_min,phi_max);
	G4LogicalVolume* GEM2CuFoil1Logical
		= new G4LogicalVolume(GEM2CuFoil1Solid,copper,"GEM2CuFoil1Logical",0,0,0);
	new G4PVPlacement(RotZ90deg,G4ThreeVector(),GEM2CuFoil1Logical,
		"GEM2CuFoil1Phys",RTPCContainerLogical,0,0);
	new G4PVPlacement(RotZ270deg,G4ThreeVector(),GEM2CuFoil1Logical,
		"GEM2CuFoil1Phys",RTPCContainerLogical,0,0);

	r_min=mGEM2R+0.005*mm;
	r_max=mGEM2R+0.005*mm+0.050*mm;
	phi_min=asin((0.5*mBedPlateThick+mG10FR4Thick+mGEM1SpThick)/r_min);//6.621 *deg;
	phi_max=180.0*deg-2.0*phi_min;
	G4VSolid* GEM2KaptonSolid
		= new G4Tubs("GEM2KaptonTubs",r_min,r_max,Z_Half,phi_min,phi_max);
	G4LogicalVolume* GEM2KaptonLogical
		= new G4LogicalVolume(GEM2KaptonSolid,kapton,"GEM2KaptonLogical",0,0,0);
	new G4PVPlacement(RotZ90deg,G4ThreeVector(),GEM2KaptonLogical,
		"GEM2KaptonPhys",RTPCContainerLogical,0,0);
	new G4PVPlacement(RotZ270deg,G4ThreeVector(),GEM2KaptonLogical,
		"GEM2KaptonPhys",RTPCContainerLogical,0,0);

	r_min=mGEM2R+0.005*mm+0.050*mm;
	r_max=mGEM2R+0.005*mm+0.050*mm+0.005*mm;
	phi_min=asin((0.5*mBedPlateThick+mG10FR4Thick+mGEM1SpThick)/r_min);//6.621 *deg;
	phi_max=180.0*deg-2.0*phi_min;
	G4VSolid* GEM2CuFoil2Solid
		= new G4Tubs("GEM2CuFoil2Tubs",r_min,r_max,Z_Half,phi_min,phi_max);
	G4LogicalVolume* GEM2CuFoil2Logical
		= new G4LogicalVolume(GEM2CuFoil2Solid,copper,"GEM2CuFoil2Logical",0,0,0);
	new G4PVPlacement(RotZ90deg,G4ThreeVector(),GEM2CuFoil2Logical,
		"GEM2CuFoil2Phys",RTPCContainerLogical,0,0);
	new G4PVPlacement(RotZ270deg,G4ThreeVector(),GEM2CuFoil2Logical,
		"GEM2CuFoil2Phys",RTPCContainerLogical,0,0);

	GEM2CuFoil1Logical->SetSensitiveDetector(SDGEM2);
	GEM2KaptonLogical->SetSensitiveDetector(SDGEM2);
	GEM2CuFoil2Logical->SetSensitiveDetector(SDGEM2);
	GEM2CuFoil1Logical->SetVisAttributes(CuBrownVisAtt);
	GEM2KaptonLogical->SetVisAttributes(CuBrownVisAtt);
	GEM2CuFoil2Logical->SetVisAttributes(CuBrownVisAtt);

	//////////////////////////
	//GEM2 support
	//////////////////////////
	x_min=mGEM2R+0.005*mm+0.050*mm+0.005*mm;
	x_max=mBedPlateHighEdge;
	y_min=mBedPlateThick/2.0+mG10FR4Thick+mGEM1SpThick;
	xx_h=(x_max-x_min)/2.0;
	yy_h=(mGEM2SpThick/2.0) *mm;
	zz_h=Z_Half;
	G4VSolid* GEM2SpSolid = new G4Box("GEM2SpBox",xx_h,yy_h,zz_h);
	G4LogicalVolume* GEM2SpLogical
		= new G4LogicalVolume(GEM2SpSolid,ultem,"GEM2SpLogical",0,0,0);
	//left half
	new G4PVPlacement(RotZ90deg,G4ThreeVector(y_min+yy_h,-x_min-xx_h,0.),GEM2SpLogical,
		"GEM2SpPhysL1",RTPCContainerLogical,0,0);
	new G4PVPlacement(RotZ270deg,G4ThreeVector(-y_min-yy_h,-x_min-xx_h,0.),GEM2SpLogical,
		"GEM2SpPhysL2",RTPCContainerLogical,0,0);
	//right half
	new G4PVPlacement(RotZ90deg,G4ThreeVector(y_min+yy_h,x_min+xx_h,0.),GEM2SpLogical,
		"GEM2SpPhysR1",RTPCContainerLogical,0,0);
	new G4PVPlacement(RotZ270deg,G4ThreeVector(-y_min-yy_h,x_min+xx_h,0.),GEM2SpLogical,
		"GEM2SpPhysR2",RTPCContainerLogical,0,0);

	GEM2SpLogical->SetVisAttributes(YellowVisAtt); //ultem yellow

	//////////////////////////
	//Drift Region between gem2 and gem3 
	//////////////////////////
	//1 atm 300k mixture gases (bonusGas) at 63.06mm<r<66mm
	r_min=mGEM2R+0.005*mm+0.050*mm+0.005*mm;
	r_max=mGEM3R;
	phi_min=asin((0.5*mBedPlateThick+mG10FR4Thick+mGEM1SpThick+mGEM1SpThick)/r_min);
	phi_max=180.0*deg-2.0*phi_min;
	G4VSolid* driftRegionG2G3Solid
		= new G4Tubs("driftRegionG2G3Tubs",r_min,r_max,Z_Half,phi_min,phi_max);
	G4LogicalVolume* driftRegionG2G3Logical = new G4LogicalVolume(driftRegionG2G3Solid,
		bonusGas,"driftRegionG2G3Logical",0,0,0);
	new G4PVPlacement(RotZ90deg,G4ThreeVector(),driftRegionG2G3Logical,
		"driftG2G3Phys",RTPCContainerLogical,0,0);
	new G4PVPlacement(RotZ270deg,G4ThreeVector(),driftRegionG2G3Logical,
		"driftG2G3Phys",RTPCContainerLogical,0,0);
	driftRegionG2G3Logical->SetVisAttributes(HallVisAtt);
	//Drift Region between gem2 and gem3 end  


	//////////////////////////
	//GEM layer 3
	//////////////////////////
	//third GEM foil//0.05mm Kapton and 0.005 Cu on both sides at r=66mm
	r_min=mGEM3R;
	r_max=mGEM3R+0.005*mm;
	phi_min=asin((0.5*mBedPlateThick+mG10FR4Thick+
		mGEM1SpThick+mGEM2SpThick)/r_min);//10.376 *deg;
	phi_max=180.0*deg-2.0*phi_min;
	G4VSolid* GEM3CuFoil1Solid
		= new G4Tubs("GEM3CuFoil1Tubs",r_min,r_max,Z_Half,phi_min,phi_max);
	G4LogicalVolume* GEM3CuFoil1Logical
		= new G4LogicalVolume(GEM3CuFoil1Solid,copper,"GEM3CuFoil1Logical",0,0,0);
	new G4PVPlacement(RotZ90deg,G4ThreeVector(),GEM3CuFoil1Logical,
		"GEM3CuFoil1Phys",RTPCContainerLogical,0,0);
	new G4PVPlacement(RotZ270deg,G4ThreeVector(),GEM3CuFoil1Logical,
		"GEM3CuFoil1Phys",RTPCContainerLogical,0,0);

	r_min=mGEM3R+0.005*mm;
	r_max=mGEM3R+0.005*mm+0.050*mm;
	phi_min=asin((0.5*mBedPlateThick+mG10FR4Thick+
		mGEM1SpThick+mGEM2SpThick)/r_min);//10.376 *deg;
	phi_max=180.0*deg-2.0*phi_min;
	G4VSolid* GEM3KaptonSolid
		= new G4Tubs("GEM3KaptonTubs",r_min,r_max,Z_Half,phi_min,phi_max);
	G4LogicalVolume* GEM3KaptonLogical
		= new G4LogicalVolume(GEM3KaptonSolid,kapton,"GEM3KaptonLogical",0,0,0);
	new G4PVPlacement(RotZ90deg,G4ThreeVector(),GEM3KaptonLogical,
		"GEM3KaptonPhys",RTPCContainerLogical,0,0);
	new G4PVPlacement(RotZ270deg,G4ThreeVector(),GEM3KaptonLogical,
		"GEM3KaptonPhys",RTPCContainerLogical,0,0);

	r_min=mGEM3R+0.005*mm+0.050*mm;
	r_max=mGEM3R+0.005*mm+0.050*mm+0.005*mm;
	phi_min=asin((0.5*mBedPlateThick+mG10FR4Thick
		+mGEM1SpThick+mGEM2SpThick)/r_min);//10.376 *deg;
	phi_max=180.0*deg-2.0*phi_min;
	G4VSolid* GEM3CuFoil2Solid
		= new G4Tubs("GEM3CuFoil2Tubs",r_min,r_max,Z_Half,phi_min,phi_max);
	G4LogicalVolume* GEM3CuFoil2Logical
		= new G4LogicalVolume(GEM3CuFoil2Solid,copper,"GEM3CuFoil2Logical",0,0,0);
	new G4PVPlacement(RotZ90deg,G4ThreeVector(),GEM3CuFoil2Logical,
		"GEM3CuFoil2Phys",RTPCContainerLogical,0,0);
	new G4PVPlacement(RotZ270deg,G4ThreeVector(),GEM3CuFoil2Logical,
		"GEM3CuFoil2Phys",RTPCContainerLogical,0,0);

	GEM3CuFoil1Logical->SetSensitiveDetector(SDGEM3);
	GEM3KaptonLogical->SetSensitiveDetector(SDGEM3);
	GEM3CuFoil2Logical->SetSensitiveDetector(SDGEM3);
	GEM3CuFoil1Logical->SetVisAttributes(CuBrownVisAtt);
	GEM3KaptonLogical->SetVisAttributes(CuBrownVisAtt);
	GEM3CuFoil2Logical->SetVisAttributes(CuBrownVisAtt);

	
	//////////////////////////
	//GEM3 support
	//////////////////////////
	x_min=mGEM3R+0.005*mm+0.050*mm+0.005*mm;
	x_max=mBedPlateHighEdge;
	y_min=mBedPlateThick/2.0+mG10FR4Thick+mGEM1SpThick+mGEM2SpThick;
	xx_h=(x_max-x_min)/2.0;
	yy_h=mGEM3SpThick/2.0;
	zz_h=Z_Half;
	G4VSolid* GEM3SpSolid = new G4Box("GEM3SpBox",xx_h,yy_h,zz_h);
	G4LogicalVolume* GEM3SpLogical
		= new G4LogicalVolume(GEM3SpSolid,ultem,"GEM3SpLogical",0,0,0);
	//left half
	new G4PVPlacement(RotZ90deg,G4ThreeVector(y_min+yy_h,-x_min-xx_h,0.),GEM3SpLogical,
		"GEM3SpPhysL1",RTPCContainerLogical,0,0);
	new G4PVPlacement(RotZ270deg,G4ThreeVector(-y_min-yy_h,-x_min-xx_h,0.),GEM3SpLogical,
		"GEM3SpPhysL2",RTPCContainerLogical,0,0);
	//right half
	new G4PVPlacement(RotZ90deg,G4ThreeVector(y_min+yy_h,x_min+xx_h,0.),GEM3SpLogical,
		"GEM3SpPhysR1",RTPCContainerLogical,0,0);
	new G4PVPlacement(RotZ270deg,G4ThreeVector(-y_min-yy_h,x_min+xx_h,0.),GEM3SpLogical,
		"GEM3SpPhysR2",RTPCContainerLogical,0,0);

	GEM3SpLogical->SetVisAttributes(YellowVisAtt); //ultem yellow


	//////////////////////////
	//Drift Region between gem3 and PCB 
	//////////////////////////
	//1 atm 300k mixture gases (bonusGas) at 66.06mm<r<69mm
	r_min=mGEM3R+0.005*mm+0.050*mm+0.005*mm;
	r_max=mPCBReadOutR;
	phi_min=asin((0.5*mBedPlateThick+mG10FR4Thick
		+mGEM1SpThick+mGEM2SpThick+mGEM3SpThick)/r_min);
	phi_max=180.0*deg-2.0*phi_min;
	G4VSolid* driftRegionG3PCBSolid
		= new G4Tubs("driftRegionG3PCBTubs",r_min,r_max,Z_Half,phi_min,phi_max);
	G4LogicalVolume* driftRegionG3PCBLogical = new G4LogicalVolume(driftRegionG3PCBSolid,
		bonusGas,"driftRegionG3PCBLogical",0,0,0);
	new G4PVPlacement(RotZ90deg,G4ThreeVector(),driftRegionG3PCBLogical,
		"driftG3PCBPhys",RTPCContainerLogical,0,0);
	new G4PVPlacement(RotZ270deg,G4ThreeVector(),driftRegionG3PCBLogical,
		"driftG3PCBPhys",RTPCContainerLogical,0,0);
	driftRegionG3PCBLogical->SetVisAttributes(HallVisAtt);
	//Drift Region between gem3 and PCB end  


	//////////////////////////
	//pcbReadOut
	//////////////////////////
	r_min=mPCBReadOutR;
	r_max=mPCBReadOutR+mG10FR4Thick;	//2.731'=6.93674cm
	phi_min=asin((0.5*mBedPlateThick+mG10FR4Thick+mGEM1SpThick+
		mGEM2SpThick+mGEM3SpThick)/r_min);
	//phi_min=13.642 *deg;
	phi_max=180.0*deg-2.0*phi_min;

	G4VSolid* pcbReadOutSolid
		= new G4Tubs("pcbReadOutTubs",r_min,r_max,Z_Half,phi_min,phi_max);
	G4LogicalVolume* pcbReadOutLogical
		= new G4LogicalVolume(pcbReadOutSolid,G10FR4,"pcbReadOutLogical",0,0,0);
	new G4PVPlacement(RotZ90deg,G4ThreeVector(),pcbReadOutLogical,
		"pcbReadOutPhysL",RTPCContainerLogical,0,0);
	new G4PVPlacement(RotZ270deg,G4ThreeVector(),pcbReadOutLogical,
		"pcbReadOutPhysR",RTPCContainerLogical,0,0);

	pcbReadOutLogical->SetSensitiveDetector(SDReadOut);
	pcbReadOutLogical->SetVisAttributes(PcbGreenVisAtt); //green, color for pcb
	//pcbReadOut end  //pcbReadOut end  //pcbReadOut end  //pcbReadOut end

	//////////////////////////
	//ReadOut Support
	//////////////////////////
	//use boolean method to construct the readOutSp
	//readOutSp= Box -tub (sector start from 16 deg) - tub_small
	//build a box, then cut it at phi=16 deg, then subtract a small cylinder to 
	//form the shape of the bottom
	//readOutSp //readOutSp  //readOutSp //readOutSp
	r_min=mPCBReadOutR+mG10FR4Thick; //2.731'=6.93674cm
	r_max=r_min+150*mm; //make it 5mm larger ennough to cut at phi=phi_min
	G4double mReadOutSpThick=6.35*mm;

	y_min=0.5*mBedPlateThick+mG10FR4Thick+mGEM1SpThick+mGEM2SpThick+mGEM3SpThick;
	y_max=y_min+mReadOutSpThick;
	phi_min=asin(y_max/r_min);
	//16.*deg is BoNuS RTPC coverage start, if radius change, this angle will change
	//I want to keep this minimum angle
	if(phi_min>16*deg) phi_min=16*deg;
	phi_max=phi_min+15.*deg;
	x_min=r_min*cos(phi_min);
	x_max=mBedPlateHighEdge;
	xx_h=(x_max-x_min)/2.;
	yy_h=mReadOutSpThick/2.0;
	G4ThreeVector yTran(-x_min-xx_h, -y_min-yy_h, 0.);
	G4RotationMatrix* RotX180Z90deg = new G4RotationMatrix();
	RotX180Z90deg->rotateX(180.*deg); RotX180Z90deg->rotateZ(90.*deg);
	G4RotationMatrix* RotX180Z270deg = new G4RotationMatrix();
	RotX180Z270deg->rotateX(180.*deg); RotX180Z270deg->rotateZ(270.*deg);

	G4VSolid* Tub_readOutSp = new G4Tubs("Tub_readOutSp",
		r_min,r_max,Z_Half+5*mm,phi_min,phi_max);
	G4VSolid* Tub_readOutSp_small = new G4Tubs("Tub_readOutSp_small",
		0,r_min,Z_Half+5*mm,0.*deg,phi_min);
	G4VSolid* Box_readOutSp = new G4Box("Box_readOutSp",xx_h,yy_h,Z_Half);

	G4SubtractionSolid* BoxSubTub_readOutsp = new G4SubtractionSolid("BoxSubTub_readOutsp",
		Box_readOutSp,Tub_readOutSp,RotZero,yTran);
	G4SubtractionSolid* readOutSpSolid = new G4SubtractionSolid("readOutSp",
		BoxSubTub_readOutsp, Tub_readOutSp_small,RotZero,yTran);

	//rotate first then transform
	G4LogicalVolume* readOutSpLogical = new G4LogicalVolume(readOutSpSolid,
		ultem,"readOutSpLogical",0,0,0);

	//left half
	new G4PVPlacement(RotX180Z90deg,G4ThreeVector(y_min+yy_h,x_min+xx_h,0.),readOutSpLogical,
		"readOutSpPhysL2",RTPCContainerLogical,0,0);
	new G4PVPlacement(RotZ90deg,G4ThreeVector(y_min+yy_h,-x_min-xx_h,0.),readOutSpLogical,
		"readOutSpPhysL1",RTPCContainerLogical,0,0);
	//right half
	new G4PVPlacement(RotZ270deg,G4ThreeVector(-y_min-yy_h,x_min+xx_h,0.),readOutSpLogical,
		"readOutSpPhysR2",RTPCContainerLogical,0,0);
	new G4PVPlacement(RotX180Z270deg,G4ThreeVector(-y_min-yy_h,-x_min-xx_h,0.),readOutSpLogical,
		"readOutSpPhysR1",RTPCContainerLogical,0,0);

	readOutSpLogical->SetVisAttributes(YellowVisAtt);
	//readOutSp end  //readOutSp end //readOutSp end //readOutSp end

	//############################################
	//////////////////////////
	// virtual boundary foil
	//////////////////////////
	//To speed up,
	//place virtual boundary foil here to kill all particles which penetrate the readout pcb
	//I put 7 z planes here such that I can choose the first 4 or the whole 7 planes
	phi_min=0.*deg;
	phi_max=360.*deg;
	//good for BoNuS
	double pDZ = 50.0*mm+mRTPCLength/2; 
	double pRmax=mPCBReadOutR+45.8*mm;  //114.8 mm
	G4double zPlane_vb[]={-1.0*mm-pDZ,-pDZ,-pDZ,pDZ-32.0*mm,pDZ,pDZ,pDZ+1.0*mm};
	G4double rInner_vb[]={15.0*mm,15.0*mm,pRmax-0.1*mm,pRmax-0.1*mm,pRmax-0.1*mm,15.0*mm,15.0*mm};
	G4double rOuter_vb[]={pRmax,pRmax,pRmax,pRmax,pRmax,pRmax,pRmax};	//

	int nZplane=7;
	bool DoNotCoverDownEndPlate=false;
	int pSetupSuperBigBite=0;
	gConfig->GetParameter("SetupSuperBigBite",pSetupSuperBigBite);
	int pSetupBigBite=0;
	gConfig->GetParameter("SetupBigBite",pSetupBigBite); 
	if(pSetupBigBite || pSetupSuperBigBite)  DoNotCoverDownEndPlate=true;
	if(DoNotCoverDownEndPlate)  nZplane=4;

	G4VSolid* virtualBoundarySolid = new G4Polycone("virtualBoundaryPcon",
		phi_min,phi_max,nZplane,zPlane_vb,rInner_vb,rOuter_vb);
	G4LogicalVolume* virtualBoundaryLogical
		= new G4LogicalVolume(virtualBoundarySolid,air,"virtualBoundaryLogical",0,0,0);
	new G4PVPlacement(0,G4ThreeVector(),virtualBoundaryLogical,
		"virtualBoundaryPhys_RTPC",RTPCContainerLogical,0,0);

	virtualBoundaryLogical->SetSensitiveDetector(virtualBoundarySD);
	if(DoNotCoverDownEndPlate)
		virtualBoundaryLogical->SetVisAttributes(LightYellowVisAtt);
	else 
		virtualBoundaryLogical->SetVisAttributes(HallVisAtt);
	//virtual boundary foil end  //virtual boundary foil end  //virtual boundary foil end

	//############################################


	/////////////#ifdef ENDPLATENCOVER
	if (mSetupEndPlateNCover)
	{
		//There are many pieces of "ring" stuff on both ends of the detector. Since all
		//these rings construct a perfect plate with the same material, I just simply
		//construct a plate on both ends. (It is boring to add from piece to piece)
		
		//////////////////////////
		// End plate  
		//////////////////////////
		r_min=R20;
		r_max=mPCBReadOutR+mG10FR4Thick+9.63*mm;	//3.110'
		phi_min=0.*deg;
		phi_max=360.*deg;
		xx_h=mBedPlateHighEdge;
		yy_h=0.5*mBedPlateThick+mG10FR4Thick+mGEM1SpThick+mGEM2SpThick
			+mGEM3SpThick+mReadOutSpThick;
		double pEndPlateThick=6.35*mm;
		z_down=Z_Half+pEndPlateThick/2;
		z_up=-1.*(Z_Half+pEndPlateThick/2);
		// G4VSolid* endPlateSolid
		//  = new G4Tubs("endPlateTubs",r_min,r_max,pEndPlateThick/2,phi_min,phi_max);
		G4VSolid* Tub_endPlate_max = new G4Tubs("Tub_endPlate_large",
			0.,r_max,pEndPlateThick/2,phi_min,phi_max);
		G4VSolid* Tub_endPlate_min = new G4Tubs("Tub_endPlate_small",
			0.,r_min,pEndPlateThick/2,phi_min,phi_max);
		G4VSolid* Box_endPlate = new G4Box("Box_endPlate",xx_h,yy_h,pEndPlateThick/2);
		G4UnionSolid* UBoxNTub_max= new G4UnionSolid("Box+Tub_max", 
			Box_endPlate, Tub_endPlate_max);
		G4SubtractionSolid* endPlateSolid= new G4SubtractionSolid("endPlateSolid",
			UBoxNTub_max,Tub_endPlate_min);
		G4LogicalVolume* endPlateLogical = new G4LogicalVolume(endPlateSolid,
			ultem,"endPlateLogical",0,0,0);
		new G4PVPlacement(RotZ90deg,G4ThreeVector(0.,0.,z_down),endPlateLogical,
			"downEndPlatePhys",RTPCContainerLogical,0,0);
		new G4PVPlacement(RotZ90deg,G4ThreeVector(0.,0.,z_up),endPlateLogical,
			"upEndPlatePhys",RTPCContainerLogical,0,0);

		endPlateLogical->SetVisAttributes(YellowVisAtt); //endPlate ultem yellow
		// End plate end // End plate end // End plate end // End plate end

		//////////////////////////
		//Downstream End plate pcb cover 
		//////////////////////////
		r_min=R20;
		r_max=mPCBReadOutR+mG10FR4Thick+9.63*mm;
		phi_min=0.*deg;
		phi_max=360.*deg;
		xx_h=mBedPlateHighEdge;
		yy_h=0.5*mBedPlateThick+mG10FR4Thick+mGEM1SpThick+mGEM2SpThick
			+mGEM3SpThick+mReadOutSpThick;
		double pDownEndCoverThick=0.063*25.4*mm;   // 0.063'
		z_down=Z_Half+pEndPlateThick+pDownEndCoverThick/2;	//105.0036+0.25'+ 0.063'/2
		//G4VSolid* downEndCoverSolid
		//  = new G4Tubs("downEndCoverTubs",r_min,r_max,0.8001*mm,phi_min,phi_max);
		G4VSolid* Tub_endCover_max = new G4Tubs("Tub_endCover_large",
			0.,r_max,pDownEndCoverThick/2,phi_min,phi_max);
		G4VSolid* Tub_endCover_min = new G4Tubs("Tub_endCover_small",
			0.,r_min,pDownEndCoverThick/2,phi_min,phi_max);
		G4VSolid* Box_endCover = new G4Box("Box_endCover",xx_h,yy_h,pDownEndCoverThick/2);
		G4UnionSolid* UBoxNTub_max_cover= new G4UnionSolid("Box+Tub_max_cover", 
			Box_endCover, Tub_endCover_max);
		G4SubtractionSolid* downEndCoverSolid = new G4SubtractionSolid("downEndCoverSolid",
			UBoxNTub_max_cover,Tub_endCover_min);

		G4LogicalVolume* downEndCoverLogical = new G4LogicalVolume(downEndCoverSolid,
			G10FR4,"downEndCoverLogical",0,0,0);
		new G4PVPlacement(RotZ90deg,G4ThreeVector(0.,0.,z_down),downEndCoverLogical,
			"downEndCoverPhys",RTPCContainerLogical,0,0);

		downEndCoverLogical->SetVisAttributes(CuBrownVisAtt); //color of Cu
		//Downstream End plate pcb cover end

		//////////////////////////
		//Uptream End plate stainlesssteel cover
		//////////////////////////
		//const G4double zPlane[] ={-111.3536*mm,-113.8936*mm,-113.89536*mm,-117.7036*mm};
		//const G4double rInner[] ={12.0523*mm,12.0523*mm,33.401*mm,33.401*mm};
		//const G4double rOutner[]={79.000*mm,79.000*mm,66.675*mm,66.675*mm};

		phi_min=0.*deg;
		phi_max=360.*deg;
		z_down=-Z_Half-6.3536*mm;
		r_max=mPCBReadOutR+10*mm;
		const G4double zPlane[] ={z_down,z_down-2.5*mm,z_down-2.5*mm,z_down-6.35*mm};
		const G4double rInner[] ={12.0523*mm,12.0523*mm,33.401*mm,33.401*mm};
		const G4double rOutner[]={r_max,r_max,r_max-12.325*mm,r_max-12.325*mm};
		G4VSolid* upEndCoverSolid
			= new G4Polycone("upEndCoverPcon",phi_min,phi_max,4,zPlane,rInner,rOutner);
		G4LogicalVolume* upEndCoverLogical = new G4LogicalVolume(upEndCoverSolid,
			stainlesssteel,"upEndCoverLogical",0,0,0);
		new G4PVPlacement(0,G4ThreeVector(),upEndCoverLogical,
			"upEndCoverPhys",RTPCContainerLogical,0,0);

		upEndCoverLogical->SetVisAttributes(SteelVisAtt); //color of steel, white
		//Uptream End plate stainlesssteel cover end //Uptream End plate stainlesssteel cover end

		//////////////////////////
		//Downtream exit window rohacel71 cover
		//////////////////////////
		
		//const G4double  zPlane_d[] ={110.9538*mm,112.9538*mm,112.9538*mm,120.9538*mm,120.9538*mm,132.9538*mm};
		//const G4double  rInner_d[] ={3.10*mm,3.10*mm,3.10*mm,3.10*mm,3.10*mm,3.10*mm};
		//const G4double rOutner_d[] ={19.90*mm,19.90*mm,30.00*mm,30.00*mm,12.70*mm,12.70*mm};

		phi_min=0.*deg;
		phi_max=360.*deg;
		z_up=Z_Half+5.9538*mm; //110.9538*mm
		r_min=mD2GasR+0.1*mm;  //3.1*mm
		const G4double  zPlane_d[] ={z_up,z_up+2*mm,z_up+2*mm,z_up+10*mm,z_up+10*mm,z_up+22*mm};
		const G4double  rInner_d[] ={r_min,r_min,r_min,r_min,r_min,r_min};
		const G4double rOutner_d[] ={m1stMylarR,m1stMylarR,m2ndMylarR,m2ndMylarR,12.70*mm,12.70*mm};
		G4VSolid* exitCoverSolid = new G4Polycone("exitCoverPcon",
			phi_min,phi_max,6,zPlane_d,rInner_d,rOutner_d);
		G4LogicalVolume* exitCoverLogical = new G4LogicalVolume(exitCoverSolid,
			rohacel71,"exitCoverLogical",0,0,0);
		new G4PVPlacement(0,G4ThreeVector(),exitCoverLogical,
			"exitCoverPhys",RTPCContainerLogical,0,0);

		//the following are aliminum of phoa stuff, so let their color to be white
		exitCoverLogical->SetVisAttributes(WhiteVisAtt);
		//Downtream exit window rohacel71 cover end //Downtream exit window rohacel71 cover end

		//////////////////////////
		//Downtream helium wall & tape
		//////////////////////////
		r_min=12.70*mm;
		r_max=r_min+0.001*25.4*mm;
		phi_min=0.*deg;
		phi_max=360.*deg;
		z_up=Z_Half+15.9538*mm;
		G4VSolid* heWallSolid = new G4Tubs("heWallTub",r_min,r_max,14.523*mm,phi_min,phi_max);

		G4LogicalVolume* heWallLogical = new G4LogicalVolume(heWallSolid,
			kapton,"heTubWallLogical",0,0,0);
		new G4PVPlacement(0,G4ThreeVector(0,0,z_up+14.523*mm),heWallLogical,
			"heTubWallPhys",RTPCContainerLogical,0,0);

		// 1 mil mylar tape , 1 mil = 0.0001 inch =0.0254 mm, 19.3 mm  wide
		zz_h=6.65*mm;
		G4VSolid* heWallTapeSolid = new G4Tubs("heWallTub",
			r_max,r_max+0.0254*mm,zz_h,phi_min,phi_max);

		G4LogicalVolume* heWallTapeLogical = new G4LogicalVolume(heWallTapeSolid,
			mylar,"heWallTapeLogical",0,0,0);
		new G4PVPlacement(0,G4ThreeVector(0,0,z_up+zz_h),heWallTapeLogical,
			"heWallTapePhys",RTPCContainerLogical,0,0);
		heWallLogical->SetVisAttributes(LightYellowVisAtt);
		heWallTapeLogical->SetVisAttributes(YellowVisAtt);
		///Downtream helium wall & tape end  ///Downtream helium wall end
	}
	/////////////#endif

	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	////////////////#ifdef SetupCableNChip
	//version 2.0 2007/3/12
	// along Z axis, in BoNuS, place 5 chips in each row, 40 rows in total
	//if Half_Z is long enough, can place more
	if(mSetupCableNChip)
	{
		//////////////////////////
		//pcbNchip
		//////////////////////////
		//build the pcbNchip box xyz=35 x 2.5 x 31 mm
		//locate at (r_min=69.3674, phi=0, z_min=-80.0mm)

		double pPCBRout=mPCBReadOutR+mG10FR4Thick;

		//atan(1.25/69.3674)=1.023*deg;
		G4double phi_off_L=(-70.8-1.2)*deg;
		G4double phi_off_R=(109.2-1.2)*deg;
		G4VSolid* pcbNchipSolid = new G4Box("pcbNchipBox",17.5*mm,1.25*mm,15.5*mm);
		G4LogicalVolume* pcbNchipLogical
			= new G4LogicalVolume(pcbNchipSolid,pcbNchip,"pcbNchipLogical",0,0,0);

		const int nrow=40;
		const double colZSpan=40*mm;
		const int ncol=int((mRTPCLength-8*mm)/colZSpan);

		//z for the 1st col(most upstream)
		double z_col0 = -mRTPCLength/2 + (mRTPCLength-ncol*colZSpan)/2 + colZSpan/2;
		int i,j;
		G4double r_chip=pPCBRout+17.5*mm;
		G4double r_conn=pPCBRout+2.75*mm;
		G4double r_cable=pPCBRout+13.0*mm+10.75*mm;
		G4double phi[nrow];
		G4double d_phi=atan((2.5+1.25+1.25)/pPCBRout);
		G4double dphi_conn=atan((1.25+1.25)/pPCBRout);
		G4RotationMatrix *Rm_chip[nrow], *Rm_cable[nrow];
		for(int i=0;i<nrow;i++)
		{
			//rotate the assembly inside the mother volume
			if(i<nrow/2) phi[i]=phi_off_L+(7.4*i)*deg;
			else phi[i]=phi_off_R+(7.4*(i-nrow/2))*deg;
			Rm_chip[i]=new G4RotationMatrix();
			Rm_chip[i]->rotateZ(phi[i]*-1);
			//G4PVPlacement rotate clockwise but MakeImprint  rotate anti-clockwise //
			//the phi position of the cable is off 2.5 degree (center to center)
			Rm_cable[i]=new G4RotationMatrix();
			//Rm_cable[i]->rotateZ(phi[i]+d_phi);
			Rm_cable[i]->rotateZ(phi[i]);
		}

		//now place the chips
		//this is a stupid way, I should use G4PVReplica
		for(i=0;i<nrow;i++)
		{	// along phi
			for (j=0;j<ncol;j++)
			{	
				G4ThreeVector Tm(r_chip*cos(phi[i]),r_chip*sin(phi[i]),z_col0+j*colZSpan);
				new G4PVPlacement(Rm_chip[i],Tm,pcbNchipLogical,
					"pcbNchipPhys",RTPCContainerLogical,true,i*nrow+j);
			}
		}

		//////////////////////////
		//connectors
		//////////////////////////
		//
		G4ThreeVector tmpTran(21.0*mm, 0.*mm, -12.75*mm);
		G4VSolid* Box_cnn_cable = new G4Box("Box_cnn_cable",10.75*mm,1.25*mm,2.75*mm);
		G4VSolid* Box_cnn_chip = new G4Box("Box_cnn_chip",2.75*mm,1.25*mm,7.0*mm);
		G4UnionSolid* connectorSolid = new G4UnionSolid("connectorSolid",
			Box_cnn_chip,Box_cnn_cable,RotZero,tmpTran);
		G4LogicalVolume* connectorLogical = new G4LogicalVolume(connectorSolid,
			cable,"connectorLogical",0,0,0);

		for(i=0;i<nrow;i++)
		{	// along phi
			for (j=0;j<ncol;j++)
			{	// along Z axis
				G4ThreeVector Tm(r_conn*cos(phi[i]+dphi_conn),
					r_conn*sin(phi[i]+dphi_conn),z_col0+j*colZSpan);
				new G4PVPlacement(Rm_chip[i],Tm,connectorLogical,
					"connectorPhys",RTPCContainerLogical,true,i*nrow+j);
			}
		}
		//connectors end //connectors end	//connectors end	//connectors end

		//////////////////////////
		//cable
		//////////////////////////
		//just keep this for information
		////cable 1 xyz=21.5mm * 0.4mm * 200.0mm locate at (r_min=80.mm, phi=0, z_min=-125.0mm)
		//G4VSolid* cable1Solid = new G4Box("cable1BoxSolid",1.075*cm,0.02*cm,10.0 *cm);
		//G4LogicalVolume* cable1LV
		//	= new G4LogicalVolume(cable1Solid,cable,"cb1LV",0,0,0);

		////cable 2 xyz=21.5mm * 0.4mm * 160.0mm locate at (r=80, phi=0, z_min=-125.0mm)
		//G4VSolid* cable2Solid = new G4Box("cable2BoxSolid",1.075*cm,0.02*cm,8.0 *cm);
		//G4LogicalVolume* cable2LV
		//	= new G4LogicalVolume(cable2Solid,cable,"cb2LV",0,0,0);

		////cable 3 xyz=21.5mm * 0.4mm * 120.0mm locate at (r=80, phi=0, z_min=-125.0mm)
		//G4VSolid* cable3Solid = new G4Box("cable3BoxSolid",1.075*cm,0.02*cm,6.0 *cm);
		//G4LogicalVolume* cable3LV
		//	= new G4LogicalVolume(cable3Solid,cable,"cb3LV",0,0,0);

		////cable 4 xyz=21.5mm * 0.4mm * 80.0mm locate at (r=80, phi=0, z_min=-125.0mm)
		//G4VSolid* cable4Solid = new G4Box("cable4BoxSolid",1.075*cm,0.02*cm,4.0 *cm);
		//G4LogicalVolume* cable4LV
		//	= new G4LogicalVolume(cable4Solid,cable,"cb4LV",0,0,0);

		////cable 5 xyz=21.5mm * 0.4mm * 40.0mm locate at (r=80, phi=0, z_min=-125.0mm)
		//G4VSolid* cable5Solid = new G4Box("cable5BoxSolid",1.075*cm,0.02*cm,2.0 *cm);
		//G4LogicalVolume* cable5LV
		//	= new G4LogicalVolume(cable5Solid,cable,"cb5LV",0,0,0);

		//G4AssemblyVolume* cableAssembly=new G4AssemblyVolume();
		//G4RotationMatrix* Ra = new G4RotationMatrix();
		//G4ThreeVector Ta(0.,0., 0.);
		//cableAssembly->AddPlacedVolume(cable1LV,Ta,Ra);
		//Ta.set(0.,-0.5*mm, -20.0*mm);
		//cableAssembly->AddPlacedVolume(cable2LV,Ta,Ra);
		//Ta.set(0.,-1.0*mm, -40.0*mm);
		//cableAssembly->AddPlacedVolume(cable3LV,Ta,Ra);
		//Ta.set(0.,-1.5*mm, -60.0*mm);
		//cableAssembly->AddPlacedVolume(cable4LV,Ta,Ra);
		//Ta.set(0.,-2.0*mm, -80.0*mm);
		//cableAssembly->AddPlacedVolume(cable5LV,Ta,Ra);

		////instantiate the cable assembly
		//for(i=0;i<40;i++)
		//{
		//	G4ThreeVector Tm(r_cable*cos(phi[i]+d_phi),r_cable*sin(phi[i]+d_phi),-34.0*mm);
		//	cableAssembly->MakeImprint(RTPCContainerLogical,Tm,Rm_cable[i],i);
		//}

		////sensortive the detector
		//pcbNchipLogical->SetSensitiveDetector(SDCableNChip);
		//connectorLogical->SetSensitiveDetector(SDCableNChip);
		//cable1LV->SetSensitiveDetector(SDCableNChip);
		//cable2LV->SetSensitiveDetector(SDCableNChip);
		//cable3LV->SetSensitiveDetector(SDCableNChip);
		//cable4LV->SetSensitiveDetector(SDCableNChip);
		//cable5LV->SetSensitiveDetector(SDCableNChip);

		//// visualization attributes ------------------------------------
		//connectorLogical->SetVisAttributes(BlackVisAtt);
		//pcbNchipLogical->SetVisAttributes(PcbGreenVisAtt);
		//cable1LV->SetVisAttributes(CableVisAtt);
		//cable2LV->SetVisAttributes(CableVisAtt);
		//cable3LV->SetVisAttributes(CableVisAtt);
		//cable4LV->SetVisAttributes(CableVisAtt);
		//cable5LV->SetVisAttributes(CableVisAtt);

		//a new way to build the cable
		//////////////////////////
		//cable
		//////////////////////////

		//cable Assembly, xyz=21.5mm * 0.4mm * cablelength_# 
		//located at (r=80, phi=0, z_min=-125.0mm)
		
		double cableLength_0=45*mm;  //length of the first cable
		G4AssemblyVolume* cableAssembly=new G4AssemblyVolume();
		G4RotationMatrix* Ra = new G4RotationMatrix();
		G4ThreeVector Ta(0.,0., 0.);
		G4VSolid** cableSolid=new G4VSolid* [ncol]; 
		G4LogicalVolume** cableLogical=new G4LogicalVolume* [ncol];
		double cableYSpan=0.5*mm;
		double cableThick=0.3*mm;

		char tmpStr[255];
		for (j=0;j<ncol;j++)
		{
			//j=0 is the most upstream cable, most short one
			//each cable is attached to the upstream end of a chip
			//all cables start from one cablelength_0(45mm here) upstream of the most upstream chip
			//therefore the z locations of cables are:  -cablelength_0/2+z_col0+j*colZSpan;
			double cableLength=cableLength_0+j*colZSpan;
			double chip_z_up=z_col0+(j-0.5)*colZSpan; //the upstream edge of the chip 
			cableSolid[j] = new G4Box("cable5BoxSolid",1.075*cm,cableThick,cableLength/2);
			sprintf(tmpStr,"cb%dLV",j+1);
			cableLogical[j] = new G4LogicalVolume(cableSolid[j],cable,tmpStr,0,0,0);
			Ta.set(0.,-2.5*mm+cableYSpan*j,chip_z_up-cableLength/2);
			cableAssembly->AddPlacedVolume(cableLogical[j],Ta,Ra);
		}

		//instantiate the cable assembly
		for(i=0;i<nrow;i++)
		{
			//shift the cable by 4.5 mm since the chip is only 31 mm and th ecolspan is 40 mm
			G4ThreeVector Tm(r_cable*cos(phi[i]+d_phi),r_cable*sin(phi[i]+d_phi),4.5*mm);
			cableAssembly->MakeImprint(RTPCContainerLogical,Tm,Rm_cable[i],i);
		}
		//sensortive the detector
		pcbNchipLogical->SetSensitiveDetector(SDCableNChip);
		connectorLogical->SetSensitiveDetector(SDCableNChip);
		for (j=0;j<ncol;j++) cableLogical[j]->SetSensitiveDetector(SDCableNChip);

		// visualization attributes ------------------------------------
		connectorLogical->SetVisAttributes(BlackVisAtt);
		pcbNchipLogical->SetVisAttributes(PcbGreenVisAtt);
		for (j=0;j<ncol;j++) cableLogical[j]->SetVisAttributes(CableVisAtt);

	}
	///////////#endif


	return phyRTPCContainer;
}




//setup solenoid geometry
//if mSetupSolenoid==1, build DVCS solenoid, which is used in CLAS DVCS and BoNuS
//if mSetupSolenoid==2, build UVA solenoid, provided by Gorden Gates
G4VPhysicalVolume* RTPCDetectorConstruction::ConstructSolenoid(G4LogicalVolume *pMotherLogVol)
{
	double startphi=0.*deg, deltaphi=360.*deg;

	G4VPhysicalVolume* thePhysVol=0;

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

	G4VSolid* solenoidSolid = 0;
	G4LogicalVolume* solenoidLogical = 0;

	if(mSetupSolenoid==1)
	{
		//DVCS solenoid
		startphi=0.*deg; deltaphi=360.*deg;
		const int kNPlane_DVCSCoil=3;
		double rInner_DVCSCoil[] = {115*mm,115*mm,313.3*mm};
		double rOuter_DVCSCoil[] = {456*mm,456*mm,456*mm};
		double zPlane_DVCSCoil[] = {-510*mm,106.2*mm,228.8*mm};

		solenoidSolid = new G4Polycone("DVCSCoilPolycone",startphi,deltaphi,
			kNPlane_DVCSCoil,zPlane_DVCSCoil,rInner_DVCSCoil,rOuter_DVCSCoil);

	}
	else if(mSetupSolenoid==2)
	{
		//UVA Gorden Cates's solenoid
		startphi=0.*deg; deltaphi=360.*deg;
		//assuming it is a cylinder of Rin=200mm and Rout=851mm, L=1527mm
		solenoidSolid = new G4Tubs("UVASolenoidTubs",200.*mm,851.*mm,
			1527.*mm/2.0,startphi,deltaphi);
	}

	solenoidLogical = new G4LogicalVolume(solenoidSolid,
		stainlesssteel,"RTPCCoilLogical",0,0,0);
	//solenoidLogical->SetVisAttributes(SteelVisAtt); 
	solenoidLogical->SetVisAttributes(DarkBlueVisAtt);   //dark blue looks better in 3D view

	//note that mSolenoidPosX,mSolenoidPosY,mSolenoidPosZ are in Hall coordinate,
	//if the pMotherLogVol is not the worldLogVol, one need to do the transfromation
	thePhysVol=new G4PVPlacement(0,
		G4ThreeVector(mSolenoidPosX-mTargetXOffset,mSolenoidPosY-mTargetYOffset,
		mSolenoidPosZ-mTargetZOffset),
		solenoidLogical,"SolenoidPhys",pMotherLogVol,false,0);

	return thePhysVol;
}

