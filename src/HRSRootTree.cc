// ********************************************************************
//
// $Id: HRSRootTree.cc,v 1.3, 2013/02/26 HRS Exp $
// GEANT4 tag $Name: geant4-09-04 $
// 
// ********************************************************************

#include "HRSRootTree.hh"
#include "UsageManager.hh"
#include "GlobalDebuger.hh"

#include "G4RunManager.hh"
#include "HRSPrimaryGeneratorAction.hh"

#include "TObject.h"
#include "math.h"
#include <string>
#include <iostream>
#include "HRSTransform_TCSNHCS.hh"
#include "HRSTranUseSNAKE.hh"
#include "HRSRand.hh"
#include "TMath.h"
#include "HRSRecUseDB.hh"

#include "HRSGlobal.hh"					//incldue some string routines
#include "../XSModel/XSModel.hh"	    //for Elas, PBosted, QFS_N_EPC, Wiser and Compton XS

#include "Drift_Sieve2Tg.hh"

using namespace std;

//use this option to save space, some of these raw steps will not be written 
//into the ntuple: StepTL, StepDsty, StepKE, StepRadlen,
#define HRS_MINI_TREE 1
#define DEBUG_HRS_RECONSTRUCTION 1

//use this option to put the hit at the virtual boundary into the ntuple
//#define HRS_TREE_DEBUG 2

extern UsageManager* gConfig;


double GetInelasXS(int PID,double Eb, double Theta, double Pf, G4Material* theMaterial)
{	
	/*
	//here is the debug detail for G4Material
	//one can use this following to get the element list in the vertex material
	G4Material* theMaterial = theTrack->GetLogicalVolumeAtVertex()->GetMaterial(); 
	const G4ElementVector* theElementV = theMaterial->GetElementVector();
	const G4double* theMassFractionV = theMaterial->GetFractionVector();  //fraction of mass

	G4cout<<"Material Name = "<<theMaterial->GetName()<<G4endl;
	for(int ii=0;ii<theElementV->size();ii++)
	{
	G4Element* element = (*theElementV)[ii];	
	cout<< " Element: " << std::setw(20)<< element->GetName() 
	<< "   Z = " << std::setw(4) << std::setprecision(2) <<  element->GetZ() 
	<< "   N = " << std::setw(5) << std::setprecision(2) <<  element->GetN()
	<< "   A = " << std::setw(6) << std::setprecision(4)
	<< element->GetA()/(g/mole) << " g/mole" ;
	cout<<"  MassFration="<<theMassFractionV[ii]<<endl;

	}
	//The output will look like the following:
	//Material Name = SolidNH3+LiquidHe
	//Element:                    H   Z =    1   N =     1   A =  1.008 g/mole  MassFration=0.09765
	//Element:                    N   Z =    7   N =    14   A =  14.01 g/mole  MassFration=0.4523
	//Element:             LiquidHe   Z =    2   N =     4   A =  4.003 g/mole  MassFration=0.45
	*/	

	if(Eb<Pf) return 0.0;
	const G4ElementVector* theElementV = theMaterial->GetElementVector();
	const G4double* theMassFractionV = theMaterial->GetFractionVector();  //fraction of mass

	int nElement = theElementV->size();
	if(nElement<1) return 0.0;


	double *pXS = new double [nElement];
	double *pNofMol = new double [nElement];
	double pNofMol_total=0;
	double pXS_eff=0;

	G4Element* element = 0;
	for(int ii=0;ii<nElement;ii++)
	{	
		element = (*theElementV)[ii];
		int Z = int(element->GetZ());
		int N = int(element->GetN())-Z;
		//By Jixie: this step need very long time to finish in some cases
		//pXS[ii] = QFS_N_EPC::getQElXS(PID, Z, N, Eb, Theta, Pf);
		pXS[ii] = PBosted::GetXS(PID, Z, N, Eb, Theta, Pf);
		//assuming there is 1 unit mass (g) of this material, get the number of mol of this element
		pNofMol[ii] = theMassFractionV[ii]/element->GetA(); 
		pNofMol_total += pNofMol[ii];

		pXS_eff += pXS[ii]*pNofMol[ii]; 
	}
	if(pNofMol_total>0.0) pXS_eff /= pNofMol_total;

	delete [] pXS;
	delete [] pNofMol;

	return pXS_eff;
}


HRSRootTree::HRSRootTree(int pRunNumber)
{
#ifdef HRS_TREE_DEBUG
	if( HRS_TREE_DEBUG > Global_Debug_Level ) 
		SetGlobalDebugLevel("HRSRootTree",(int) HRS_TREE_DEBUG);		
#endif

	iRunNumber=pRunNumber;

	iNoDetectorResponse=0;
	gConfig->GetArgument("NoDetectorResponse",iNoDetectorResponse); 

	iEvtCounter=0;
	iRootEvtID=0;
	bConfigTreeFilled=false;

	//for cmd lime option -o or -out
	//$t description: Sepecify the output file name and create mode.This option will be invoked only if 
	//$t the length of OutFileName larger than 3. The program will ovverwrite the exist files if 
	//$t Recreate is not zero. If Recreate is zero, will automaticly add a subrun index and increase this subrun 
	//$t index to make a new name for the output file. If OutFileName is not provided or shorter than 4, 
	//$t the program will use 'G2P_G4Sim_nt_$Run_$Subrun.root' with none-overwritten mode.

	char strRawName[100];
	std::string pOutFileName=gConfig->GetArgument("OutFileName");
	int pRecreate=0;
	gConfig->GetArgument("Recreate",pRecreate);
	if(pOutFileName.length()>3) 
	{ 
		if(!pRecreate) 
		{
			sprintf(strRawName,"%s",pOutFileName.c_str());
			CreateFileName(strRawName,rootfile);
		}
		else
		{			
			strcpy(rootfile,pOutFileName.c_str()); 
		}
	}
	else
	{
		sprintf(strRawName,"G4Sim_nt_%02d.root",iRunNumber);
		CreateFileName(strRawName,rootfile);
	}

	for (int i=0;i<MaxPrimaryNum;i++) {tree[i] = 0; track[i] = new MyTrack();}

	//set the initial value to 'counts'
	TA1_N=TA2_N=T_N=SD_N=0;
	Reset();

	file=0;  //TFile set to null use this to identified if the tree has been Initilized

	std::cout<<"HRSRootTree() construction done!"<<std::endl;
}

void HRSRootTree::Initilize()
{
	printf("HRSRootTree::HRSRootTree(): Initializing HRSRootTree......\n");


	/////////////////////////////////////////////////////
	mSetupLHRS=mSetupRHRS=0;
	gConfig->GetParameter("SetupLHRS",mSetupLHRS);
	gConfig->GetParameter("SetupRHRS",mSetupRHRS);

	gConfig->GetParameter("TargetXOffset",mTargetXOffset); //in mm
	gConfig->GetParameter("TargetYOffset",mTargetYOffset); //in mm
	gConfig->GetParameter("TargetZOffset",mTargetZOffset); //in mm
	gConfig->GetParameter("PivotXOffset",mPivotXOffset);   //in mm
	gConfig->GetParameter("PivotYOffset",mPivotYOffset);   //in mm
	gConfig->GetParameter("PivotZOffset",mPivotZOffset);   //in mm

	/////////////////////////////////////////////////////
	//this block takes care of arguments
	gConfig->GetArgument("BeamEnergy",mBeamEnergy); //GeV
	gConfig->GetArgument("LHRSMomentum",mLHRSMomentum); //GeV
	gConfig->GetArgument("RHRSMomentum",mRHRSMomentum); //GeV
	gConfig->GetArgument("BeamTiltedAngle",mBeamTiltedAngle); //deg
	mBeamTiltedAngle*=deg;

	gConfig->GetArgument("TargetMass",mTargetMass); //GeV
	gConfig->GetArgument("TargetAtomicNumber",mTargetAtomicNumber); 
	gConfig->GetArgument("TargetNeutronNumber",mTargetNeutronNumber); 

	mUseSeptumPlusStdHRS=0;
	gConfig->GetArgument("UseSeptumPlusStdHRS",mUseSeptumPlusStdHRS);


	gConfig->GetArgument("FZB1TiltedAngle",mFZB1TiltedAngle); 
	mFZB1TiltedAngle*=deg;
	gConfig->GetArgument("FZB1PosX",mFZB1PosX); 
	gConfig->GetArgument("FZB1PosY",mFZB1PosY); 
	gConfig->GetArgument("FZB1PosZ",mFZB1PosZ);
	gConfig->GetArgument("FZB1Bx",mFZB1Bx); 
	gConfig->GetArgument("FZB1By",mFZB1By); 
	gConfig->GetArgument("FZB1Bz",mFZB1Bz); 


	gConfig->GetArgument("FZB2TiltedAngle",mFZB2TiltedAngle); 
	mFZB2TiltedAngle*=deg;
	gConfig->GetArgument("FZB2PosX",mFZB2PosX); 
	gConfig->GetArgument("FZB2PosY",mFZB2PosY); 
	gConfig->GetArgument("FZB2PosZ",mFZB2PosZ);
	gConfig->GetArgument("FZB2Bx",mFZB2Bx); 
	gConfig->GetArgument("FZB2By",mFZB2By); 
	gConfig->GetArgument("FZB2Bz",mFZB2Bz); 


	mBPMXRes=0.5; //mm;	
	gConfig->GetArgument("BPMXRes",mBPMXRes);	
	mBPMYRes=1.0; //mm;	
	gConfig->GetArgument("BPMYRes",mBPMYRes);
	//specify where to stop reconstruction
	//0 is at target plane, 1 is at vertex plane, 2 is at exact vertex z
	mWhereToStopRecon=0;
	gConfig->GetArgument("WhereToStopRecon",mWhereToStopRecon); 

	mSnakeModel=11;
	gConfig->GetArgument("SnakeModel",mSnakeModel);
	if(mUseSeptumPlusStdHRS)  
	{
		cout<<"You have specified to use option '-UseSeptumPlusStdHRS', so the standard 12.5 SNAKE model will be used\n";
		mSnakeModel=1;
	}

	//Do not use HRS optics reconstruction is HRS is not used
	mUseOpticsDB=0;
	if(mSetupLHRS || mSetupRHRS)
	{
		gConfig->GetArgument("UseOpticsDB", mUseOpticsDB);
		if(mUseOpticsDB==1)
		{
			//string pHRSOpticsDBL,pHRSOpticsDBR;
			string pHRSOpticsDBL=gConfig->GetArgument("HRSOpticsDBL");
			string pHRSOpticsDBR=gConfig->GetArgument("HRSOpticsDBR");
			mRecUseDBL = new HRSRecUseDB("L",pHRSOpticsDBL.c_str());
			mRecUseDBR = new HRSRecUseDB("R",pHRSOpticsDBR.c_str());

#if defined DEBUG_HRS_RECONSTRUCTION  && (DEBUG_HRS_RECONSTRUCTION >= 2)
			mRecUseDBL->PrintDataBase();
			mRecUseDBR->PrintDataBase();
#endif
		}
	}

	/////////////////////////////////////////////////////
	//For HRS parameters
	mPivot2LHRSVBFace=mPivot2RHRSVBFace=1300.0;

	gConfig->GetParameter("LSeptumAngle",mLSeptumAngle);   //deg
	mLSeptumAngle*=deg;
	gConfig->GetParameter("RSeptumAngle",mRSeptumAngle);   //deg
	mRSeptumAngle*=deg;
	gConfig->GetParameter("LHRSAngle",mLHRSAngle);         //deg
	mLHRSAngle*=deg;
	gConfig->GetParameter("RHRSAngle",mRHRSAngle);         //deg
	mRHRSAngle*=deg;

	if(mSetupLHRS || mSetupRHRS)
	{
		gConfig->GetParameter("Pivot2LHRSVBFace",mPivot2LHRSVBFace);
		gConfig->GetParameter("Pivot2RHRSVBFace",mPivot2RHRSVBFace);

		int pSeptum_UseUniformB=0;
		gConfig->GetParameter("Septum_UseUniformB",pSeptum_UseUniformB); 
		mSeptumCurrentRatioL=mSeptumCurrentRatioR=1.0;  
		gConfig->GetParameter("Septum_CurrentRatioL",mSeptumCurrentRatioL);
		gConfig->GetParameter("Septum_CurrentRatioR",mSeptumCurrentRatioR);
		mUseSeptumField=(pSeptum_UseUniformB==0 && ( fabs(mSeptumCurrentRatioL)>1.0E-08 || 
			fabs(mSeptumCurrentRatioR)>1.0E-08) )?1:0;
		gConfig->GetParameter("Septum_OriginX",mSeptumXOffset); //in cm
		mSeptumXOffset*=10.0;
		gConfig->GetParameter("Septum_OriginY",mSeptumYOffset); //in cm
		mSeptumYOffset*=10.0;
		gConfig->GetParameter("Septum_OriginZ",mSeptumZOffset); //in cm
		mSeptumZOffset*=10.0;
		gConfig->GetParameter("Septum_RotAxis1",mSeptumRotAxis1); 
		gConfig->GetParameter("Septum_RotAngle1",mSeptumRotAngle1); //deg
		mSeptumRotAngle1*=deg;
		gConfig->GetParameter("Septum_RotAxis2",mSeptumRotAxis2); 
		gConfig->GetParameter("Septum_RotAngle2",mSeptumRotAngle2); //deg
		mSeptumRotAngle2*=deg;
		gConfig->GetParameter("Septum_RotAxis3",mSeptumRotAxis3); 
		gConfig->GetParameter("Septum_RotAngle3",mSeptumRotAngle3); //deg
		mSeptumRotAngle3*=deg;
	}

	//////////////////////////////////////////////////////
	//target field parameters, note that the unit for length is cm!!!
	int pHelm_UseUniformB=0;
	gConfig->GetParameter("Helm_UseUniformB",pHelm_UseUniformB); 
	mHelmCurrentRatio=1.0;  
	gConfig->GetParameter("Helm_CurrentRatio",mHelmCurrentRatio);
	mUseHelmField=(pHelm_UseUniformB==0 && fabs(mHelmCurrentRatio)>1.0E-08)?1:0;
	gConfig->GetParameter("Helm_OriginX",mHelmXOffset); //in cm
	mHelmXOffset*=10.0;
	gConfig->GetParameter("Helm_OriginY",mHelmYOffset); //in cm
	mHelmYOffset*=10.0;
	gConfig->GetParameter("Helm_OriginZ",mHelmZOffset); //in cm
	mHelmZOffset*=10.0;
	gConfig->GetParameter("Helm_RotAxis1",mHelmRotAxis1); 
	gConfig->GetParameter("Helm_RotAngle1",mHelmRotAngle1); //deg
	mHelmRotAngle1*=deg;
	gConfig->GetParameter("Helm_RotAxis2",mHelmRotAxis2); 
	gConfig->GetParameter("Helm_RotAngle2",mHelmRotAngle2); //deg
	mHelmRotAngle2*=deg;
	gConfig->GetParameter("Helm_RotAxis3",mHelmRotAxis3); 
	gConfig->GetParameter("Helm_RotAngle3",mHelmRotAngle3); //deg
	mHelmRotAngle3*=deg;

	/////////////////////////////////////////////////////
	//VD
	gConfig->GetParameter("SetupVirtualDetector",mSetupVD); 
	if(mSetupVD)
	{
		gConfig->GetParameter("VDRotYAngle",mVDAngle); //deg
		mVDAngle*=deg;
		gConfig->GetParameter("VDRotXAngle",mVDTiltAngle); //deg
		mVDTiltAngle*=deg;
		gConfig->GetParameter("Pivot2VDFace",mPivot2VDFace); //mm
	}

	/////////////////////////////////////////////////////
	//LAC
	gConfig->GetParameter("SetupLAC",mSetupLAC); 
	if(mSetupLAC)
	{
		gConfig->GetParameter("LACAngle",mLACAngle); //deg
		mLACAngle*=deg;
		gConfig->GetParameter("LACTiltAngle",mLACTiltAngle); //deg
		mLACTiltAngle*=deg;
		gConfig->GetParameter("Pivot2LACFace",mPivot2LACFace); //mm
	}

	/////////////////////////////////////////////////////
	//BigBite
	gConfig->GetParameter("SetupBigBite",mSetupBigBite); 
	if(mSetupBigBite)
	{
		gConfig->GetParameter("BigBiteAngle",mBigBiteAngle); //deg
		mBigBiteAngle*=deg;
		gConfig->GetParameter("BigBiteTiltAngle",mBigBiteTiltAngle); //deg
		mBigBiteTiltAngle*=deg;
		gConfig->GetParameter("Pivot2BigBiteFace",mPivot2BigBiteFace); //mm

		gConfig->GetParameter("SetupHAND",mSetupHAND); 
		gConfig->GetParameter("Pivot2HANDLeadWall",mPivot2HANDLeadWall); //mm
	}

	/////////////////////////////////////////////////////
	//SBS
	gConfig->GetParameter("SetupSuperBigBite",mSetupSuperBigBite); 
	if(mSetupSuperBigBite)
	{
		gConfig->GetParameter("SuperBigBiteAngle",mSuperBigBiteAngle); //deg
		mSuperBigBiteAngle*=deg;
		gConfig->GetParameter("Pivot2SuperBigBiteFace",mPivot2SuperBigBiteFace); //mm

		int pTosca_SBS_UseUniformB=0;
		gConfig->GetParameter("Tosca_SBS_UseUniformB",pTosca_SBS_UseUniformB); 
		mSBSFieldCurrentRatio=1.0;  
		gConfig->GetParameter("Tosca_SBS_CurrentRatio",mSBSFieldCurrentRatio);
		mUseSBSField=( pTosca_SBS_UseUniformB==0 && fabs(mSBSFieldCurrentRatio)>1.0E-08 )?1:0;
		gConfig->GetParameter("Tosca_SBS_OriginX",mSBSFieldXOffset); //in cm
		mSBSFieldXOffset*=10.0;
		gConfig->GetParameter("Tosca_SBS_OriginY",mSBSFieldYOffset); //in cm
		mSBSFieldYOffset*=10.0;
		gConfig->GetParameter("Tosca_SBS_OriginZ",mSBSFieldZOffset); //in cm
		mSBSFieldZOffset*=10.0;
		gConfig->GetParameter("Tosca_SBS_RotAxis1",mSBSFieldRotAxis1); 
		gConfig->GetParameter("Tosca_SBS_RotAngle1",mSBSFieldRotAngle1); //deg
		mSBSFieldRotAngle1*=deg;
		gConfig->GetParameter("Tosca_SBS_RotAxis2",mSBSFieldRotAxis2); 
		gConfig->GetParameter("Tosca_SBS_RotAngle2",mSBSFieldRotAngle2); //deg
		mSBSFieldRotAngle2*=deg;
		gConfig->GetParameter("Tosca_SBS_RotAxis3",mSBSFieldRotAxis3); 
		gConfig->GetParameter("Tosca_SBS_RotAngle3",mSBSFieldRotAngle3); //deg
		mSBSFieldRotAngle3*=deg;
	}

	/////////////////////////////////////////////////////
	//HMS
	gConfig->GetParameter("SetupHMS",mSetupHMS); 
	if(mSetupHMS)
	{
		gConfig->GetParameter("HMSAngle",mHMSAngle); //deg
		mHMSAngle*=deg;
		gConfig->GetParameter("Pivot2HMSFace",mPivot2HMSFace); //mm
		gConfig->GetParameter("HMSMomentum",mHMSMomentum); //GeV
	}

	/////////////////////////////////////////////////////
	//G2P stuff
	mSetupG2PGeometry=0;
	gConfig->GetParameter("SetupG2PGeometry",mSetupG2PGeometry);

	mThirdArmAngle=mThirdArmRotZAngle=mPivot2ThirdArmFace=0;
	mSetupG2PTarget=mTargetType=0;
	mSetupChicane=mSetupChicaneVD=0;

	if(mSetupG2PGeometry)
	{
		gConfig->GetParameter("ThirdArmAngle",mThirdArmAngle); //deg
		mThirdArmAngle*=deg;
		gConfig->GetParameter("SetupThirdArmVD",mSetupThirdArmVD); 
		gConfig->GetParameter("ThirdArmRotZAngle",mThirdArmRotZAngle); //deg
		mThirdArmRotZAngle*=deg;
		mThirdArmRotZAngle=-36.5*deg*mHelmCurrentRatio;
		gConfig->GetParameter("Pivot2ThirdArmFace",mPivot2ThirdArmFace); //mm

		gConfig->GetParameter("SetupG2PTarget",mSetupG2PTarget); 
		gConfig->GetParameter("TargetType",mTargetType); 
		gConfig->GetParameter("TargetL",mTargetL); //mm

		gConfig->GetParameter("SetupChicane",mSetupChicane); 
		gConfig->GetParameter("SetupChicaneVD",mSetupChicaneVD);
	}

	/////////////////////////////////////////////////////
	mSetupCREXGeometry=0;
	gConfig->GetParameter("SetupCREXGeometry",mSetupCREXGeometry);
	if(mSetupCREXGeometry)
	{
		gConfig->GetParameter("TargetType",mTargetType); 
		gConfig->GetParameter("TargetL",mTargetL); //mm
	}

	/////////////////////////////////////////////////////
	mSetupRTPC=0;
	gConfig->GetParameter("SetupRTPCGeometry",mSetupRTPC);
	if(mSetupRTPC)
	{
		gConfig->GetParameter("RatioHe2DME",mRatioHe2DME); 
		gConfig->GetParameter("TargetL",mTargetL); //mm
	}


	//need to update these parameters
	G4RunManager *theRunManager=G4RunManager::GetRunManager(); 
	mPrimaryGeneratorAction = 
		(HRSPrimaryGeneratorAction*) theRunManager->GetUserPrimaryGeneratorAction();

	mBeamEnergy=mPrimaryGeneratorAction->GetBeamEnergy()/GeV;
	mLHRSMomentum=mPrimaryGeneratorAction->GetLeftHRSMomentum()/GeV;
	mRHRSMomentum=mPrimaryGeneratorAction->GetRightHRSMomentum()/GeV;
	mTargetMass=mPrimaryGeneratorAction->GetTargetMass()/GeV;
	mTargetAtomicNumber=mPrimaryGeneratorAction->GetTargetAtomicNumber();

#ifdef HRS_TREE_DEBUG
	if(Global_Debug_Level>=2) gConfig->PrintOpt(); 
	if(Global_Debug_Level>=1) 
	{
		cout<<"///////////////////////////////////////////////////////////"<<endl;
		cout<<"In this run, TargetMass="<<mTargetMass
			<<"  TargetAtomicNumber="<<mTargetAtomicNumber<<endl;
		cout<<"///////////////////////////////////////////////////////////"<<endl;
		if(Global_Debug_Level>=2) STOP4DEBUG;
	}
#endif

	/////////////////////////////////////////////////////////////////////////////
	mGenHistoOnly=0;
	gConfig->GetArgument("GenHistoOnly",mGenHistoOnly);
	gConfig->GetArgument("SkimLevel",iSkimLevel);
	gConfig->GetParameter("BookTrees",iBookTrees); 
	if(iBookTrees>MaxPrimaryNum) iBookTrees=MaxPrimaryNum;
	gConfig->GetParameter("BookHistos",iBookHistos);
	int pStoreTrajectory=0;
	gConfig->GetArgument("StoreTrajectory",pStoreTrajectory);
	mCalculateXS=1;
	gConfig->GetArgument("CalculateXS",mCalculateXS);

	bConfigTreeFilled=false;
	file=new TFile(rootfile,"RECREATE");

	// create the trees only if iBookTrees >= 1 and mGenHistoOnly!=0
	if (iBookTrees >= 1)
	{
		config=new TTree("config","run configuration");

		//from argument

		config->Branch("Run",&iRunNumber,"Run/I");
		config->Branch("SkimLevel",&iSkimLevel,"SkimLevel/I");
		config->Branch("BookTrees",&iBookTrees,"BookTrees/I");
		config->Branch("Beam",&mBeamEnergy,"Beam/D");
		config->Branch("BeamTiltedAngle",&mBeamTiltedAngle,"BeamTiltedAngle/D");
		config->Branch("TargetM",&mTargetMass,"TargetM/D");
		config->Branch("TargetL",&mTargetL,"TargetL/D");
		config->Branch("TargetAtomicNumber",&mTargetAtomicNumber,"TargetAtomicNumber/D");
		config->Branch("TargetNeutronNumber",&mTargetNeutronNumber,"TargetNeutronNumber/D");


		//some variable for reconstruction 
		config->Branch("WhereToStopRecon",&mWhereToStopRecon,"WhereToStopRecon/I");
		config->Branch("SnakeModel",&mSnakeModel,"SnakeModel/I");
		config->Branch("BPMYRes",&mBPMYRes,"BPMYRes/D");
		config->Branch("BPMXRes",&mBPMXRes,"BPMXRes/D");

		config->Branch("UseSeptumPlusStdHRS",&mUseSeptumPlusStdHRS,"UseSeptumPlusStdHRS/I");
		config->Branch("UseOpticsDB",&mUseOpticsDB,"UseOpticsDB/I");	

		/////////////////////////////////////////////////////
		//global config variable, located in detector.ini
		config->Branch("TargetXOffset",&mTargetXOffset,"TargetXOffset/D");
		config->Branch("TargetYOffset",&mTargetYOffset,"TargetYOffset/D");
		config->Branch("TargetZOffset",&mTargetZOffset,"TargetZOffset/D");
		config->Branch("PivotXOffset",&mPivotXOffset,"PivotXOffset/D");
		config->Branch("PivotYOffset",&mPivotYOffset,"PivotYOffset/D");
		config->Branch("PivotZOffset",&mPivotZOffset,"PivotZOffset/D");

		config->Branch("SetupLHRS",&mSetupLHRS,"SetupLHRS/I");		
		config->Branch("SetupRHRS",&mSetupRHRS,"SetupRHRS/I");

		config->Branch("LHRSMomentum",&mLHRSMomentum,"LHRSMomentum/D");
		config->Branch("RHRSMomentum",&mRHRSMomentum,"RHRSMomentum/D");
		config->Branch("LHRSAngle",&mLSeptumAngle,"LHRSAngle/D");
		config->Branch("RHRSAngle",&mRSeptumAngle,"RHRSAngle/D");

		config->Branch("Pivot2LHRSVBFace",&mPivot2LHRSVBFace,"Pivot2LHRSVBFace/D");
		config->Branch("Pivot2RHRSVBFace",&mPivot2RHRSVBFace,"Pivot2RHRSVBFace/D");

		/////////////////////////////////////////////////////
		config->Branch("UseHelmField",&mUseHelmField,"UseHelmField/I");
		config->Branch("HelmXOffset",&mHelmXOffset,"HelmXOffset/D");
		config->Branch("HelmYOffset",&mHelmYOffset,"HelmYOffset/D");
		config->Branch("HelmZOffset",&mHelmZOffset,"HelmZOffset/D");
		config->Branch("HelmRotAxis1",&mHelmRotAxis1,"HelmRotAxis1/D");
		config->Branch("HelmRotAxis2",&mHelmRotAxis2,"HelmRotAxis2/D");
		config->Branch("HelmRotAxis3",&mHelmRotAxis3,"HelmRotAxis3/D");
		config->Branch("HelmRotAngle1",&mHelmRotAngle1,"HelmRotAngle1/D");
		config->Branch("HelmRotAngle2",&mHelmRotAngle2,"HelmRotAngle2/D");
		config->Branch("HelmRotAngle3",&mHelmRotAngle3,"HelmRotAngle3/D");
		config->Branch("HelmCurrentRatio",&mHelmCurrentRatio,"HelmCurrentRatio/D");

		config->Branch("UseSeptumField",&mUseSeptumField,"UseSeptumField/I");
		config->Branch("SeptumXOffset",&mSeptumXOffset,"SeptumXOffset/D");
		config->Branch("SeptumYOffset",&mSeptumYOffset,"SeptumYOffset/D");
		config->Branch("SeptumZOffset",&mSeptumZOffset,"SeptumZOffset/D");
		config->Branch("SeptumRotAxis1",&mSeptumRotAxis1,"SeptumRotAxis1/D");
		config->Branch("SeptumRotAxis2",&mSeptumRotAxis2,"SeptumRotAxis2/D");
		config->Branch("SeptumRotAxis3",&mSeptumRotAxis3,"SeptumRotAxis3/D");
		config->Branch("SeptumRotAngle1",&mSeptumRotAngle1,"SeptumRotAngle1/D");
		config->Branch("SeptumRotAngle2",&mSeptumRotAngle2,"SeptumRotAngle2/D");
		config->Branch("SeptumRotAngle3",&mSeptumRotAngle3,"SeptumRotAngle3/D");
		config->Branch("SeptumCurrentRatioL",&mSeptumCurrentRatioL,"SeptumCurrentRatioL/D");
		config->Branch("SeptumCurrentRatioR",&mSeptumCurrentRatioR,"SeptumCurrentRatioR/D");


		config->Branch("SetupG2PTarget",&mSetupG2PTarget,"SetupG2PTarget/I");
		config->Branch("TargetType",&mTargetType,"TargetType/I");

		config->Branch("ThirdArmAngle",&mThirdArmAngle,"ThirdArmAngle/D");
		config->Branch("SetupThirdArmVD",&mSetupThirdArmVD,"SetupThirdArmVD/D");
		config->Branch("ThirdArmRotZAngle",&mThirdArmRotZAngle,"ThirdArmRotZAngle/D");
		config->Branch("Pivot2ThirdArmFace",&mPivot2ThirdArmFace,"Pivot2ThirdArmFace/D");

		config->Branch("SetupChicane",&mSetupChicane,"SetupChicane/I");
		config->Branch("SetupChicaneVD",&mSetupChicaneVD,"SetupChicaneVD/I");

		config->Branch("FZB1TiltedAngle",&mFZB1TiltedAngle,"FZB1TiltedAngle/D");
		config->Branch("FZB1PosX",&mFZB1PosX,"FZB1PosX/D");
		config->Branch("FZB1PosY",&mFZB1PosY,"FZB1PosY/D");
		config->Branch("FZB1PosZ",&mFZB1PosZ,"FZB1PosZ/D");
		config->Branch("FZB1Bx",&mFZB1Bx,"FZB1Bx/D");
		config->Branch("FZB1By",&mFZB1By,"FZB1By/D");
		config->Branch("FZB1Bz",&mFZB1Bz,"FZB1Bz/D");

		config->Branch("FZB2TiltedAngle",&mFZB2TiltedAngle,"FZB2TiltedAngle/D");
		config->Branch("FZB2PosX",&mFZB2PosX,"FZB2PosX/D");
		config->Branch("FZB2PosY",&mFZB2PosY,"FZB2PosY/D");
		config->Branch("FZB2PosZ",&mFZB2PosZ,"FZB2PosZ/D");
		config->Branch("FZB2Bx",&mFZB2Bx,"FZB2Bx/D");
		config->Branch("FZB2By",&mFZB2By,"FZB2By/D");
		config->Branch("FZB2Bz",&mFZB2Bz,"FZB2Bz/D");


		/////////////////////////////////////////////////////
		config->Branch("SetupVD",&mSetupVD,"SetupVD/I");
		config->Branch("VDAngle",&mVDAngle,"VDAngle/D");
		config->Branch("VDTiltAngle",&mVDTiltAngle,"VDTiltAngle/D");
		config->Branch("Pivot2VDFace",&mPivot2VDFace,"Pivot2VDFace/D");

		/////////////////////////////////////////////////////
		config->Branch("SetupLAC",&mSetupLAC,"SetupLAC/I");
		config->Branch("LACAngle",&mLACAngle,"LACAngle/D");
		config->Branch("LACTiltAngle",&mLACTiltAngle,"LACTiltAngle/D");
		config->Branch("Pivot2LACFace",&mPivot2LACFace,"Pivot2LACFace/D");

		/////////////////////////////////////////////////////
		config->Branch("SetupBigBite",&mSetupBigBite,"SetupBigBite/I");
		config->Branch("BigBiteAngle",&mBigBiteAngle,"BigBiteAngle/D");
		config->Branch("BigBiteTiltAngle",&mBigBiteTiltAngle,"BigBiteTiltAngle/D");
		config->Branch("Pivot2BigBiteFace",&mPivot2BigBiteFace,"Pivot2BigBiteFace/D");

		config->Branch("SetupHAND",&mSetupHAND,"SetupHAND/I");
		config->Branch("Pivot2HANDLeadWall",&mPivot2HANDLeadWall,"Pivot2HANDLeadWall/D");

		/////////////////////////////////////////////////////
		config->Branch("SetupSuperBigBite",&mSetupSuperBigBite,"SetupSuperBigBite/I");
		config->Branch("SuperBigBiteAngle",&mSuperBigBiteAngle,"SuperBigBiteAngle/D");
		config->Branch("Pivot2SuperBigBiteFace",&mPivot2SuperBigBiteFace,"Pivot2SuperBigBiteFace/D");


		config->Branch("UseSBSField",&mUseSBSField,"UseSBSField/I");
		config->Branch("SBSFieldXOffset",&mSBSFieldXOffset,"SBSFieldXOffset/D");
		config->Branch("SBSFieldYOffset",&mSBSFieldYOffset,"SBSFieldYOffset/D");
		config->Branch("SBSFieldZOffset",&mSBSFieldZOffset,"SBSFieldZOffset/D");
		config->Branch("SBSFieldRotAxis1",&mSBSFieldRotAxis1,"SBSFieldRotAxis1/D");
		config->Branch("SBSFieldRotAxis2",&mSBSFieldRotAxis2,"SBSFieldRotAxis2/D");
		config->Branch("SBSFieldRotAxis3",&mSBSFieldRotAxis3,"SBSFieldRotAxis3/D");
		config->Branch("SBSFieldRotAngle1",&mSBSFieldRotAngle1,"SBSFieldRotAngle1/D");
		config->Branch("SBSFieldRotAngle2",&mSBSFieldRotAngle2,"SBSFieldRotAngle2/D");
		config->Branch("SBSFieldRotAngle3",&mSBSFieldRotAngle3,"SBSFieldRotAngle3/D");
		config->Branch("SBSFieldCurrentRatio",&mSBSFieldCurrentRatio,"SBSFieldCurrentRatio/D");


		/////////////////////////////////////////////////////
		config->Branch("SetupHMS",&mSetupHMS,"SetupHMS/I");
		config->Branch("HMSAngle",&mHMSAngle,"HMSAngle/D");
		config->Branch("Pivot2HMSFace",&mPivot2HMSFace,"Pivot2HMSFace/D");
		config->Branch("HMSMomentum",&mHMSMomentum,"HMSMomentum/D");


		/////////////////////////////////////////////////////
		config->Branch("SetupRTPC",&mSetupRTPC,"SetupRTPC/I");
		config->Branch("RatioHe2DME",&mRatioHe2DME,"RatioHe2DME/D");

		/////////////////////////////////////////////////////
		if(!iNoDetectorResponse)
		{
			detector=new TTree("D","detector response tree");

			detector->Branch("Index",&iRootEvtID,"Index/I");	

			detector->Branch("TA1_N",&TA1_N,"TA1_N/I");
			detector->Branch("TA1_Pid",TA1_Pid,"TA1_Pid[TA1_N]/I");
			detector->Branch("TA1_Tid",TA1_Tid,"TA1_Tid[TA1_N]/I");
			detector->Branch("TA1_ParentTid",TA1_ParentTid,"TA1_ParentTid[TA1_N]/I");
			detector->Branch("TA1_T",TA1_T,"TA1_T[TA1_N]/D");
			detector->Branch("TA1_X",TA1_X,"TA1_X[TA1_N]/D");
			detector->Branch("TA1_Y",TA1_Y,"TA1_Y[TA1_N]/D");
			detector->Branch("TA1_Z",TA1_Z,"TA1_Z[TA1_N]/D");
			detector->Branch("TA1_Edep",TA1_Edep,"TA1_Edep[TA1_N]/D");;
			detector->Branch("TA1_NonIonEdep",TA1_NonIonEdep,"TA1_NonIonEdep[TA1_N]/D");
			detector->Branch("TA1_P",TA1_P,"TA1_P[TA1_N]/D");
			detector->Branch("TA1_Theta",TA1_Theta,"TA1_Theta[TA1_N]/D");
			detector->Branch("TA1_Phi",TA1_Phi,"TA1_Phi[TA1_N]/D");
			detector->Branch("TA1_Pout",TA1_Pout,"TA1_Pout[TA1_N]/D");

			detector->Branch("TA2_N",&TA2_N,"TA2_N/I");
			detector->Branch("TA2_Pid",TA2_Pid,"TA2_Pid[TA2_N]/I");
			detector->Branch("TA2_Tid",TA2_Tid,"TA2_Tid[TA2_N]/I");
			detector->Branch("TA2_ParentTid",TA2_ParentTid,"TA2_ParentTid[TA2_N]/I");
			detector->Branch("TA2_T",TA2_T,"TA2_T[TA2_N]/D");
			detector->Branch("TA2_X",TA2_X,"TA2_X[TA2_N]/D");
			detector->Branch("TA2_Y",TA2_Y,"TA2_Y[TA2_N]/D");
			detector->Branch("TA2_Z",TA2_Z,"TA2_Z[TA2_N]/D");
			detector->Branch("TA2_Edep",TA2_Edep,"TA2_Edep[TA2_N]/D");
			detector->Branch("TA2_NonIonEdep",TA2_NonIonEdep,"TA2_NonIonEdep[TA2_N]/D");		
			detector->Branch("TA2_P",TA2_P,"TA2_P[TA2_N]/D");
			detector->Branch("TA2_Theta",TA2_Theta,"TA2_Theta[TA2_N]/D");
			detector->Branch("TA2_Phi",TA2_Phi,"TA2_Phi[TA2_N]/D");
			detector->Branch("TA2_Pout",TA2_Pout,"TA2_Pout[TA2_N]/D");

			detector->Branch("SD_N",&SD_N,"SD_N/I");
			detector->Branch("SD_Id",SD_Id,"SD_Id[SD_N]/I");
			detector->Branch("SD_Pid",SD_Pid,"SD_Pid[SD_N]/I");
			detector->Branch("SD_Tid",SD_Tid,"SD_Tid[SD_N]/I");
			detector->Branch("SD_ParentTid",SD_ParentTid,"SD_ParentTid[SD_N]/I");
			detector->Branch("SD_T",SD_T,"SD_T[SD_N]/D");
			detector->Branch("SD_X",SD_X,"SD_X[SD_N]/D");
			detector->Branch("SD_Y",SD_Y,"SD_Y[SD_N]/D");
			detector->Branch("SD_Z",SD_Z,"SD_Z[SD_N]/D");
			detector->Branch("SD_Edep",SD_Edep,"SD_Edep[SD_N]/D");	
			detector->Branch("SD_NonIonEdep",SD_NonIonEdep,"SD_NonIonEdep[SD_N]/D");	
			detector->Branch("SD_P",SD_P,"SD_P[SD_N]/D");
			detector->Branch("SD_Theta",SD_Theta,"SD_Theta[SD_N]/D");
			detector->Branch("SD_Phi",SD_Phi,"SD_Phi[SD_N]/D");
			detector->Branch("SD_Pout",SD_Pout,"SD_Pout[SD_N]/D");

			detector->Branch("T_N",&T_N,"T_N/I");
			detector->Branch("T_Pid",T_Pid,"T_Pid[T_N]/I");
			detector->Branch("T_Tid",T_Tid,"T_Tid[T_N]/I");
			detector->Branch("T_ParentTid",T_ParentTid,"T_ParentTid[T_N]/I");
			detector->Branch("T_T",T_T,"T_T[T_N]/D");
			detector->Branch("T_X",T_X,"T_X[T_N]/D");
			detector->Branch("T_Y",T_Y,"T_Y[T_N]/D");
			detector->Branch("T_Z",T_Z,"T_Z[T_N]/D");
			detector->Branch("T_P",T_P,"T_P[T_N]/D");
			detector->Branch("T_Theta",T_Theta,"T_Theta[T_N]/D");
			detector->Branch("T_Phi",T_Phi,"T_Phi[T_N]/D");

			if(pStoreTrajectory)
			{
				detector->Branch("T_StepN",T_StepN,"T_StepN[T_N]/I");
				detector->Branch("T_StepX",T_StepX,Form("T_StepX[T_N][%d]/D",MaxTrackHit));
				detector->Branch("T_StepY",T_StepY,Form("T_StepY[T_N][%d]/D",MaxTrackHit));
				detector->Branch("T_StepZ",T_StepZ,Form("T_StepZ[T_N][%d]/D",MaxTrackHit));
			}
		}

		//create branch with standard tracks, no reconstruction variables included yet 
		char treename[100],treetitle[100];
		for(int j=0;j<iBookTrees;j++)
		{
			if(!tree[j])
			{
				sprintf(treename,"track%d",j);
				sprintf(treetitle,"track %d",j+1);
				tree[j]=new TTree(treename,treetitle);
			}
			if(!mGenHistoOnly)
			{
				tree[j]->Branch("Index",&iRootEvtID,"Index/I");		

				tree[j]->Branch("PdgId",&(track[j]->PdgId),"PdgId/I");
				tree[j]->Branch("TrackId",&(track[j]->TrackId),"TrackId/I");
				tree[j]->Branch("TrackClass",&(track[j]->TrackClass),"TrackClass/I");

				tree[j]->Branch("X0",&(track[j]->X0),"X0/D");
				tree[j]->Branch("Y0",&(track[j]->Y0),"Y0/D");
				tree[j]->Branch("Z0",&(track[j]->Z0),"Z0/D");
				tree[j]->Branch("P0",&(track[j]->P0),"P0/D");
				tree[j]->Branch("Theta0",&(track[j]->Theta0),"Theta0/D");
				tree[j]->Branch("Phi0",&(track[j]->Phi0),"Phi0/D");
				tree[j]->Branch("X0_tr",&(track[j]->X0_tr),"X0_tr/D");
				tree[j]->Branch("Y0_tr",&(track[j]->Y0_tr),"Y0_tr/D");
				tree[j]->Branch("Z0_tr",&(track[j]->Z0_tr),"Z0_tr/D");
				tree[j]->Branch("Theta0_tr",&(track[j]->Theta0_tr),"Theta0_tr/D");
				tree[j]->Branch("Phi0_tr",&(track[j]->Phi0_tr),"Phi0_tr/D");

				tree[j]->Branch("Xtg_tr",&(track[j]->Xtg_tr),"Xtg_tr/D");
				tree[j]->Branch("Ytg_tr",&(track[j]->Ytg_tr),"Ytg_tr/D");
				tree[j]->Branch("Thetatg_tr",&(track[j]->Thetatg_tr),"Thetatg_tr/D");
				tree[j]->Branch("Phitg_tr",&(track[j]->Phitg_tr),"Phitg_tr/D");

				tree[j]->Branch("Xvb",&(track[j]->Xvb),"Xvb/D");
				tree[j]->Branch("Yvb",&(track[j]->Yvb),"Yvb/D");
				tree[j]->Branch("Zvb",&(track[j]->Zvb),"Zvb/D");
				tree[j]->Branch("Pvb",&(track[j]->Pvb),"Pvb/D");
				tree[j]->Branch("Thetavb",&(track[j]->Thetavb),"Thetavb/D");
				tree[j]->Branch("Phivb",&(track[j]->Phivb),"Phivb/D");
				tree[j]->Branch("Xvb_tr",&(track[j]->Xvb_tr),"Xvb_tr/D");
				tree[j]->Branch("Yvb_tr",&(track[j]->Yvb_tr),"Yvb_tr/D");
				tree[j]->Branch("Zvb_tr",&(track[j]->Zvb_tr),"Zvb_tr/D");
				tree[j]->Branch("Thetavb_tr",&(track[j]->Thetavb_tr),"Thetavb_tr/D");
				tree[j]->Branch("Phivb_tr",&(track[j]->Phivb_tr),"Phivb_tr/D");

				tree[j]->Branch("Xfp_tr",&(track[j]->Xfp_tr),"Xfp_tr/D");
				tree[j]->Branch("Yfp_tr",&(track[j]->Yfp_tr),"Yfp_tr/D");
				tree[j]->Branch("Thetafp_tr",&(track[j]->Thetafp_tr),"Thetafp_tr/D");
				tree[j]->Branch("Phifp_tr",&(track[j]->Phifp_tr),"Phifp_tr/D");

				tree[j]->Branch("Xtg_rec_tr",&(track[j]->Xtg_rec_tr),"Xtg_rec_tr/D");
				tree[j]->Branch("Ytg_rec_tr",&(track[j]->Ytg_rec_tr),"Ytg_rec_tr/D");
				tree[j]->Branch("Thetatg_rec_tr",&(track[j]->Thetatg_rec_tr),"Thetatg_rec_tr/D");
				tree[j]->Branch("Phitg_rec_tr",&(track[j]->Phitg_rec_tr),"Phitg_rec_tr/D");

				tree[j]->Branch("X_rec_tr",&(track[j]->X_rec_tr),"X_rec_tr/D");
				tree[j]->Branch("Y_rec_tr",&(track[j]->Y_rec_tr),"Y_rec_tr/D");
				tree[j]->Branch("Z_rec_tr",&(track[j]->Z_rec_tr),"Z_rec_tr/D");
				tree[j]->Branch("Theta_rec_tr",&(track[j]->Theta_rec_tr),"Theta_rec_tr/D");
				tree[j]->Branch("Phi_rec_tr",&(track[j]->Phi_rec_tr),"Phi_rec_tr/D");
				tree[j]->Branch("X_rec",&(track[j]->X_rec),"X_rec/D");
				tree[j]->Branch("Y_rec",&(track[j]->Y_rec),"Y_rec/D");
				tree[j]->Branch("Z_rec",&(track[j]->Z_rec),"Z_rec/D");
				tree[j]->Branch("P_rec",&(track[j]->P_rec),"P_rec/D");
				tree[j]->Branch("Theta_rec",&(track[j]->Theta_rec),"Theta_rec/D");
				tree[j]->Branch("Phi_rec",&(track[j]->Phi_rec),"Phi_rec/D");

				tree[j]->Branch("Delta",&(track[j]->Delta),"Delta/D");
				tree[j]->Branch("Delta_rec",&(track[j]->Delta_rec),"Delta_rec/D");

				//Reconstructed target variable using analyzer database
				tree[j]->Branch("Xtg_rec_db_tr",&(track[j]->Xtg_rec_db_tr),"Xtg_rec_db_tr/D");
				tree[j]->Branch("Ytg_rec_db_tr",&(track[j]->Ytg_rec_db_tr),"Ytg_rec_db_tr/D");
				tree[j]->Branch("Thetatg_rec_db_tr",&(track[j]->Thetatg_rec_db_tr),"Thetatg_rec_db_tr/D");
				tree[j]->Branch("Phitg_rec_db_tr",&(track[j]->Phitg_rec_db_tr),"Phitg_rec_db_tr/D");
				tree[j]->Branch("Delta_rec_db",&(track[j]->Delta_rec_db),"Delta_rec_db/D");

				//To debug the reconstruction, I added the following leaves
				tree[j]->Branch("X_proj2tg_tr",&(track[j]->X_proj2tg_tr),"X_proj2tg_tr/D");
				tree[j]->Branch("Y_proj2tg_tr",&(track[j]->Y_proj2tg_tr),"Y_proj2tg_tr/D");

				tree[j]->Branch("X_rec2tg_tr",&(track[j]->X_rec2tg_tr),"X_rec2tg_tr/D");
				tree[j]->Branch("Y_rec2tg_tr",&(track[j]->Y_rec2tg_tr),"Y_rec2tg_tr/D");
				tree[j]->Branch("Theta_rec2tg_tr",&(track[j]->Theta_rec2tg_tr),"Theta_rec2tg_tr/D");
				tree[j]->Branch("Phi_rec2tg_tr",&(track[j]->Phi_rec2tg_tr),"Phi_rec2tg_tr/D");
				tree[j]->Branch("P_rec2tg",&(track[j]->P_rec2tg),"P_rec2tg/D");

				tree[j]->Branch("X_proj2sl_tr",&(track[j]->X_proj2sl_tr),"X_proj2sl_tr/D");
				tree[j]->Branch("Y_proj2sl_tr",&(track[j]->Y_proj2sl_tr),"Y_proj2sl_tr/D");

				tree[j]->Branch("TrackRadlen",&(track[j]->TrackRadlen),"TrackRadlen/D");
				tree[j]->Branch("Theta0Eff",&(track[j]->Theta0Eff),"Theta0Eff/D");
				tree[j]->Branch("ElasXS",&(track[j]->ElasXS),"ElasXS/D"); //elastic XS
				tree[j]->Branch("XS",&(track[j]->XS),"XS/D");  //inelastic electro-production XS from epc or qfs

#ifndef USE_MINI_TREE
				//because dynamic array in a branch of a tree[j] can not be read easily,
				//so I put the "step array" into seperate branches
				tree[j]->Branch("StepNum",&(track[j]->StepNum),"StepNum/I");
				tree[j]->Branch("StepX",&(track[j]->StepX[j]),"StepX[StepNum]/D");
				tree[j]->Branch("StepY",&(track[j]->StepY[j]),"StepY[StepNum]/D");
				tree[j]->Branch("StepZ",&(track[j]->StepZ[j]),"StepZ[StepNum]/D");
				tree[j]->Branch("StepdE",&(track[j]->StepdE[j]),"StepdE[StepNum]/D");
				tree[j]->Branch("StepL",&(track[j]->StepL[j]),"StepL[StepNum]/D");

				tree[j]->Branch("StepEkin",&(track[j]->StepEkin[j]),"StepEkin[StepNum]/D");
				tree[j]->Branch("StepTL",&(track[j]->StepTL[j]),"StepTL[StepNum]/D");
				tree[j]->Branch("StepRadlen",&(track[j]->StepRadlen[j]),"StepRadlen[StepNum]/D");
				tree[j]->Branch("StepDsty",&(track[j]->StepDsty[j]),"StepDsty[StepNum]/D");

				tree[j]->Branch("StepBx",&(track[j]->StepBx[j]),"StepBx[StepNum]/D");
				tree[j]->Branch("StepBy",&(track[j]->StepBy[j]),"StepBy[StepNum]/D");
				tree[j]->Branch("StepBz",&(track[j]->StepBz[j]),"StepBz[StepNum]/D");

				tree[j]->Branch("TrackBdLx",&(track[j]->TrackBdLx),"TrackBdLx/D");
				tree[j]->Branch("TrackBdLy",&(track[j]->TrackBdLy),"TrackBdLy/D");
				tree[j]->Branch("TrackBdLz",&(track[j]->TrackBdLz),"TrackBdLz/D");

				tree[j]->Branch("R0",&(track[j]->R0),"R0/D");			
				tree[j]->Branch("A0",&(track[j]->A0),"A0/D");
				tree[j]->Branch("B0",&(track[j]->B0),"B0/D");
#endif			

				tree[j]->Branch("Ei",&Ei,"Ei/D");
				tree[j]->Branch("Helicity",&Helicity,"Helicity/I");
				tree[j]->Branch("Pol",&(track[j]->Pol),"Pol/D");	

				tree[j]->SetMaxTreeSize((Long64_t)(9.0E+09));

			}
		}
	}

	// creat the histograms only if iBookHistos = 1
	//these histos will be used to randomize particles 
	if (iBookHistos)
	{
		cout << endl << " Creating histograms for HRSHisto engine... ";
		for(int i=0;i<iBookTrees;i++)
		{
			//declare histogram
			hP0[i]     =new TH1D(Form("P0_%d",i),"P0",4000,0.,4.);
			hTheta0[i] =new TH1D(Form("Theta0_%d",i),"Theta0",180,0.,TMath::Pi());
			hPhi0[i]   =new TH1D(Form("Phi0_%d",i),"Phi0",720,-TMath::Pi(),TMath::Pi());

			const double kMaxX=15.0, kMaxY=15.0, kMaxZ=15.0;
			hX0[i]     =new TH1D(Form("X0_%d",i),"X0",2*int(kMaxX),-kMaxX,kMaxX);
			hY0[i]     =new TH1D(Form("Y0_%d",i),"Y0",2*int(kMaxY),-kMaxY,kMaxY);
			hZ0[i]     =new TH1D(Form("Z0_%d",i),"Z0",2*int(kMaxZ),mTargetZOffset-kMaxZ,mTargetZOffset-kMaxZ); 

			h2Z0VSY0[i]=new TH2D(Form("Z0VSY0_%d",i),"Z0 VS Y0",2*int(kMaxY),-kMaxY,kMaxY,
				2*int(kMaxZ),mTargetZOffset-kMaxZ,mTargetZOffset-kMaxZ);	
			h2Y0VSX0[i]=new TH2D(Form("Y0VSX0_%d",i),"Y0 VS X0",
				2*int(kMaxX),-kMaxX,kMaxX,2*int(kMaxY),-kMaxY,kMaxY);

			h2P0VSTheta0[i]=new TH2D(Form("P0VSTheta0_%d",i),"P0 VS Theta0",
				180,0.,TMath::Pi(),4000,0.,4.);
			h2Theta0VSPhi0[i]=new TH2D(Form("Theta0VSPhi0_%d",i),"Theta0 VS Phi0",
				180,0.,TMath::Pi(),360,-TMath::Pi(),TMath::Pi());

			h2Theta0_trVSPhi0_tr[i]=new TH2D(Form("Theta0_trVSPhi0_tr_%d",i),
				"Theta0_tr VS Phi0_tr",1000,-0.5,0.5,1000,-0.5,0.5);
			h2Thetavb_trVSPhivb_tr[i]=new TH2D(Form("Thetavb_trVSPhivb_tr_%d",i),
				"Thetavb_tr VS Phivb_tr",200,-0.1,0.1,200,-0.1,0.1);

			h3Z0Y0X0[i]=new TH3D(Form("Z0Y0X0_%d",i),"Z0 VS Y0 VS X0", 2*int(kMaxX),-kMaxX,kMaxX,
				2*int(kMaxY),-kMaxY,kMaxY, 2*int(kMaxZ),mTargetZOffset-kMaxZ,mTargetZOffset-kMaxZ);	

			//h3P0Theta0Phi0[i]=new TH3D(Form("P0Theta0Phi0_%d",i),"P0 VS Theta0 VS Phi0",
			//	360,-TMath::Pi(),TMath::Pi(),180,0.,TMath::Pi(),4000,0.,4.);
		}
		cout << "...done!" << endl;
	}

	printf("HRSRootTree::Initialize(): Initializing HRSRootTree... done!\n");


	mHRSTranModel = new HRSTransport();
	mHRSTranModel->ChangeModel(mSnakeModel);
	printf("HRSRootTree::Initialize(): Initializing SNAKE Models... done!\n");
}

//////////////////////////////////////////////////////////////////////////////////////


HRSRootTree::~HRSRootTree()
{
	if(file)
	{
		if (!bConfigTreeFilled)	config->Fill();
		file->Write("",TObject::kOverwrite);
		file->Close();
		cout << "\nClose root file " << rootfile <<endl;  
		file->Delete();
		for (int i=0;i<MaxPrimaryNum;i++)
		{
			if (track[i]) delete track[i];
		}
	}

	if(mHRSTranModel) delete mHRSTranModel;

	std::cout<<"delete HRSRootTree ... done!"<<std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////
//Fill primary track tree, tree[0-7], only if track[0]->TrackClass >= iSkimLevel
void HRSRootTree::FillTree()
{
	//don't fill the histo if it is a bad event
	if (track[0]->TrackClass<iSkimLevel) return;
	//////////////////////////////////////////
	//Fill only good or golden track into mornitoring histo
	if (iBookHistos==1)
	{
		for(int t=0;t<iBookTrees;t++)
		{
			if( (t==0 && track[t]->TrackClass>=0) || (t>0 && track[t]->TrackClass>=0) )
			{
				hP0[t]->Fill(track[t]->P0);
				hTheta0[t]->Fill(track[t]->Theta0);
				hPhi0[t]->Fill(track[t]->Phi0);

				hX0[t]->Fill(track[t]->X0);
				hY0[t]->Fill(track[t]->Y0);
				hZ0[t]->Fill(track[t]->Z0);

				h2Z0VSY0[t]->Fill(track[t]->Y0,track[t]->Z0);
				h2Y0VSX0[t]->Fill(track[t]->X0,track[t]->Y0);

				h2P0VSTheta0[t]->Fill(track[t]->Theta0,track[t]->P0);
				h2Theta0VSPhi0[t]->Fill(track[t]->Phi0,track[t]->Theta0);

				h2Theta0_trVSPhi0_tr[t]->Fill(track[t]->Phi0_tr,track[t]->Theta0_tr);
				h2Thetavb_trVSPhivb_tr[t]->Fill(track[t]->Phivb_tr,track[t]->Thetavb_tr); 

				h3Z0Y0X0[t]->Fill(track[t]->X0,track[t]->Y0,track[t]->Z0);
				//h3P0Theta0Phi0[t]->Fill(track[t]->Phi0,track[t]->Theta0,track[t]->P0);
			}
		}
	}
	//////////////////////////////////////////
	//fill config tree only once
	if (iBookTrees >= 1)
	{
		if (!bConfigTreeFilled)
		{
			//Fill SD table, the value have been filled in HRSEventAction

			config->Branch("SDNum",&mSDNum,"SDNum/I");
			//I try to build the SD table in this way but it kep only the 1st str
			//config->Branch("SDName",&mSDName,"SDName[SDNum][100]/C");
			//so I build it in this way:
			char tmpName[100];
			char tmpTitle[100];
			const char *strN,*strT;
			strN=tmpName;
			strT=tmpTitle;
			for(int ihc=0;ihc<mSDNum;ihc++)
			{
				sprintf(tmpName,"SDNameForSDID%d",ihc);
				sprintf(tmpTitle,"SDNameForSDID%d[100]/C",ihc);
				//just for information
				//the next 2 lines do not work, due to the fact that 
				//virtual TBranch* TTree::Branch(const char*, void*, const char*, Int_t)
				// only take const char* as input, not an arry
				//config->Branch(tmpName,mSDName[ihc],tmpTitle);
				//config->Branch(Form("SDName_%d",ihc),&mSDName[ihc][0],Form("SDName_%d[100]/C",ihc));

				config->Branch(strN,mSDName[ihc],strT);
			}

			config->Fill();
			bConfigTreeFilled=true;
		}

		if(!mGenHistoOnly) tree[0]->Fill();
		iRootEvtID+=1;
		for (int i=1;i<MaxPrimaryNum;i++)//iSkimLevel only apply to track0
		{
			if (track[i]->TrackClass>=-1 && iBookTrees>i && !mGenHistoOnly)  
				tree[i]->Fill();
		}		

#ifdef HRS_TREE_DEBUG
		if (Global_Debug_Level>=1)
		{
			cout<<"One root event added. iRootEvtID="<<iRootEvtID<<" track 0 ";
			for (int i=1;i<MaxPrimaryNum;i++)
			{
				if (track[i]->TrackClass>=0) cout<<i<<" ";
			}
			cout<<"Filled!\n";
		}

		if (Global_Debug_Level>=1) Stop4Debug(1);
#endif
		if (!(iEvtCounter%1000) && iRootEvtID) file->Write("",TObject::kOverwrite);
	}
}



//////////////////////////////////////////////////////////////////////////////////////
void HRSRootTree::Reset()//reset the variables
{
	for (int i=0;i<MaxPrimaryNum;i++)
	{
		for (int j=0;j<track[i]->StepNum;j++)
		{
			track[i]->StepX[j]=-150.0;
			track[i]->StepY[j]=-150.0;
			track[i]->StepZ[j]=-150.0;
			track[i]->StepdE[j]=-1.;
			track[i]->StepEkin[j]=-1.;
			track[i]->StepL[j]=-1.;
			track[i]->StepTL[j]=-1.;
			track[i]->StepRadlen[j]=-1.;
			track[i]->StepDsty[j]=-1.;

			track[i]->StepBx[j]=0.0;
			track[i]->StepBy[j]=0.0;
			track[i]->StepBz[j]=0.0;
		}
		track[i]->PdgId =-1;
		track[i]->TrackId = i;
		track[i]->TrackClass=-1;

		track[i]->X0=-150.0;
		track[i]->Y0=-150.0;
		track[i]->Z0=-150.0;
		track[i]->P0=-10.0;
		track[i]->Theta0=-10.0;
		track[i]->Phi0=-10.0;
		track[i]->X0_tr=-150.0;
		track[i]->Y0_tr=-150.0;
		track[i]->Z0_tr=-150.0;
		track[i]->Theta0_tr=-10.0;
		track[i]->Phi0_tr=-10.0;

		track[i]->Xtg_tr=-150.0;
		track[i]->Ytg_tr=-150.0;
		track[i]->Thetatg_tr=-10.0;
		track[i]->Phitg_tr=-10.0;
		track[i]->Xtg_rec_tr=-150.0;
		track[i]->Ytg_rec_tr=-150.0;
		track[i]->Thetatg_rec_tr=-10.0;
		track[i]->Phitg_rec_tr=-10.0;

		track[i]->Xtg_rec_db_tr=-150.0;
		track[i]->Ytg_rec_db_tr=-150.0;
		track[i]->Thetatg_rec_db_tr=-10.0;
		track[i]->Phitg_rec_db_tr=-10.0;
		track[i]->Delta_rec_db=-10.0;

		track[i]->Xvb=-150.0;
		track[i]->Yvb=-150.0;
		track[i]->Zvb=-150.0;
		track[i]->Pvb=-10.0;
		track[i]->Thetavb=-10.0;
		track[i]->Phivb=-10.0;
		track[i]->Xvb_tr=-150.0;
		track[i]->Yvb_tr=-150.0;
		track[i]->Zvb_tr=-150.0;
		track[i]->Thetavb_tr=-10.0;
		track[i]->Phivb_tr=-10.0;


		track[i]->Xfp_tr=-150.0;
		track[i]->Yfp_tr=-150.0;
		track[i]->Thetafp_tr=-10.0;
		track[i]->Phifp_tr=-10.0;

		track[i]->X_rec_tr=-150.0;
		track[i]->Y_rec_tr=-150.0;
		track[i]->Z_rec_tr=-150.0;
		track[i]->Theta_rec_tr=-10.0;
		track[i]->Phi_rec_tr=-10.0;
		track[i]->X_rec=-150.0;
		track[i]->Y_rec=-150.0;
		track[i]->Z_rec=-150.0;
		track[i]->P_rec=0.0;
		track[i]->Theta_rec=-10.0;
		track[i]->Phi_rec=-10.0;

		track[i]->Delta=-10.0;
		track[i]->Delta_rec=-10.0;

		track[i]->TrackBdLx=0.0;
		track[i]->TrackBdLy=0.0;
		track[i]->TrackBdLz=0.0;

		track[i]->Theta0Eff=0.0;

		track[i]->R0=0.0;
		track[i]->A0=0.0;
		track[i]->B0=0.0;
		track[i]->TrackRadlen=0.0;

		track[i]->StepNum=0;

		track[i]->ElasXS=0.0;
		track[i]->XS=0.0;

		track[i]->X_proj2tg_tr=-150.0;
		track[i]->Y_proj2tg_tr=-150.0;

		track[i]->X_rec2tg_tr=-150.0;
		track[i]->Y_rec2tg_tr=-150.0;
		track[i]->Theta_rec2tg_tr=-10.0;
		track[i]->Phi_rec2tg_tr=-10.0;
		track[i]->P_rec2tg=0.0;

		track[i]->X_proj2sl_tr=-150.0;
		track[i]->Y_proj2sl_tr=-150.0;

		track[i]->Pol=0.0;
	}
	PartNum=0;	  //primary particle number

	if(!iNoDetectorResponse)
	{
		for(int i=0;i<TA1_N;i++)
		{
			TA1_Pid[i]=0;
			TA1_Tid[i]=0;
			TA1_ParentTid[i]=0;
			TA1_T[i]=0;
			TA1_X[i]=-10;
			TA1_Y[i]=-10;
			TA1_Z[i]=-10;
			TA1_Edep[i]=0;
			TA1_NonIonEdep[i]=0;
			TA1_P[i]=-10;
			TA1_Theta[i]=-10;
			TA1_Phi[i]=-10;
			TA1_Pout[i]=-10;
		}
		TA1_N=0;

		for(int i=0;i<TA2_N;i++)
		{
			TA2_Pid[i]=0;
			TA2_Tid[i]=0;
			TA2_ParentTid[i]=0;
			TA2_T[i]=0;
			TA2_X[i]=-10;
			TA2_Y[i]=-10;
			TA2_Z[i]=-10;
			TA2_Edep[i]=0;
			TA2_NonIonEdep[i]=0;
			TA2_P[i]=-10;
			TA2_Theta[i]=-10;
			TA2_Phi[i]=-10;
			TA2_Pout[i]=-10;
		}
		TA2_N=0;

		for(int i=0;i<SD_N;i++)
		{
			SD_Id[i]=0;
			SD_Pid[i]=0;
			SD_Tid[i]=0;
			SD_ParentTid[i]=0;
			SD_T[i]=0;
			SD_X[i]=-10;
			SD_Y[i]=-10;
			SD_Z[i]=-10;
			SD_Edep[i]=0;
			SD_NonIonEdep[i]=0;
			SD_P[i]=-10;
			SD_Theta[i]=-10;
			SD_Phi[i]=-10;
			SD_Pout[i]=-10;
		}
		SD_N=0;

		for(int i=0;i<T_N;i++)
		{
			T_Pid[i]=0;
			T_Tid[i]=0;
			T_ParentTid[i]=0;
			T_T[i]=0;
			T_X[i]=-10;
			T_Y[i]=-10;
			T_Z[i]=-10;
			T_P[i]=-10;
			T_Theta[i]=-10;
			T_Phi[i]=-10;						
		}

		int pStoreTrajectory=0;
		gConfig->GetArgument("StoreTrajectory",pStoreTrajectory);
		if(pStoreTrajectory)
		{
			for(int i=0;i<T_N;i++)
			{
				T_StepN[i]=0;
				for(int ss=0;ss<MaxTrackHit;ss++)
				{
					T_StepX[i][ss]=0;
					T_StepY[i][ss]=0;
					T_StepZ[i][ss]=0;
				}
			}
		}
		T_N=0;
	}//End of noDetector Response
}

//////////////////////////////////////////////////////////////////////////////////////

bool HRSRootTree::TransportThruVD(int i, double pEndPlaneAngle)
{
	MyTrack *pTrack=track[i];
	if(pTrack->Pvb/keV<1.0) return false;

	//convert the initial point (vertex) from HCS to TCS
	Transform::X_HCS2TCS(pTrack->X0-mPivotXOffset,pTrack->Y0-mPivotYOffset,
		pTrack->Z0-mPivotZOffset, pEndPlaneAngle, pTrack->X0_tr, pTrack->Y0_tr,pTrack->Z0_tr);
	Transform::P_HCS2TCS(pTrack->Theta0, pTrack->Phi0, pEndPlaneAngle, 
		pTrack->Theta0_tr, pTrack->Phi0_tr);

	////////////////////////////////////////////////////////////////////
	//Project|Drift this 0_tr to target plane tg_tr
	TVector3 X,P;  //ThreeVector in Hall system for Drift(), in unit of m, GeV
	double pCharge = mPrimaryGeneratorAction->GetPDGCharge(i);  
	double pMass = mPrimaryGeneratorAction->GetPDGMass(i)/GeV;  
#ifdef HRS_TREE_DEBUG
	if(HRS_TREE_DEBUG >= 3)
		cout<<"TransportThruVD(): Charge="<<pCharge<<"  Mass="<<pMass
			<<"  Name="<<mPrimaryGeneratorAction->GetParticleName()<<endl;
#endif
	if(!mUseHelmField || fabs(pCharge)<1.0E-8)
	{
		double tmpZtg_tr=0;
		Transform::Project(pTrack->X0_tr, pTrack->Y0_tr,pTrack->Z0_tr, -pTrack->Z0_tr, pTrack->Theta0_tr, 
			pTrack->Phi0_tr,pTrack->Xtg_tr,pTrack->Ytg_tr,tmpZtg_tr);

		pTrack->Thetatg_tr=pTrack->Theta0_tr;
		pTrack->Phitg_tr=pTrack->Phi0_tr;
	}
	else
	{
		//Drift from vertex plane to target plane
		//need to drift forward if Z0_tr<0 or backward if Z0_tr>0
		X.SetXYZ(pTrack->X0/1000.,pTrack->Y0/1000.,pTrack->Z0/1000.);
		P.SetMagThetaPhi(pTrack->P0,pTrack->Theta0,pTrack->Phi0); 

		int ret=0;
		double pZtrLimit=0.0, pTLLimit=1.0;
		if(pTrack->Z0_tr-pZtrLimit<-1.0E-6)
		{
			//forward
			DriftSieve2Tg::Drift2Ztr(X,P,pEndPlaneAngle,pMass, pCharge,pZtrLimit,
				pTLLimit,0.00005,mTargetZOffset/1000);				
		}
		else if(pTrack->Z0_tr-pZtrLimit>1.0E-6)
		{
			//backward
			P*=-1;	//flip P for backward, make sure the particle is positron
			DriftSieve2Tg::Drift2Ztr(X,P,pEndPlaneAngle,pMass,-pCharge,pZtrLimit,
				pTLLimit,0.00005,mTargetZOffset/1000);				
			P*=-1;  //flip it back
		}
		if(ret<0) return false; 

		//now convert the result from HCS to TCS, pay attention to the unit 
		double tmpZtg_tr=0;
		Transform::X_HCS2TCS(X.x()*1000-mPivotXOffset,X.y()*1000-mPivotYOffset, X.z()*1000-mPivotZOffset,
			pEndPlaneAngle, pTrack->Xtg_tr,pTrack->Ytg_tr,tmpZtg_tr);
		Transform::P_HCS2TCS(P.Theta(), P.Phi(), pEndPlaneAngle,pTrack->Thetatg_tr,pTrack->Phitg_tr);		
	}

	////////////////////////////////////////////////////////////////////
	//Convert HCS to TCS at VB
	Transform::X_HCS2TCS(pTrack->Xvb-mPivotXOffset,pTrack->Yvb-mPivotYOffset,
		pTrack->Zvb-mPivotZOffset, pEndPlaneAngle, pTrack->Xvb_tr, pTrack->Yvb_tr,pTrack->Zvb_tr);
	Transform::P_HCS2TCS(pTrack->Thetavb, pTrack->Phivb, pEndPlaneAngle, 
		pTrack->Thetavb_tr, pTrack->Phivb_tr);

	//determine the TrackClass
	if( pTrack->Pvb/GeV > 0.000015)  pTrack->TrackClass=1;

	if(pTrack->TrackClass==1)
	{
		if( (pTrack->P0 - pTrack->Pvb) / pTrack->P0 < 0.20 ) pTrack->TrackClass=2;
	}

	/*
	//there is no reconstruction right now, the following should be replaced 
	//with the real reconstruction code
	pTrack->Theta_rec_tr=pTrack->Theta0_tr;
	pTrack->Phi_rec_tr=pTrack->Phi0_tr;
	Transform::P_TCS2HCS(pTrack->Theta_rec_tr, pTrack->Phi_rec_tr, pEndPlaneAngle, 
	pTrack->Theta_rec, pTrack->Phi_rec);
	*/
	
	return true;
}

////////////////////////////////////////////////////////////////////////////////////////////

bool HRSRootTree::TransportThruBigBite(int i)
{
	MyTrack *pTrack=track[i];
	if(pTrack->Pvb/keV<1.0) return false;

	//double pEndPlaneAngle=mPGA->GetDetectorAngle(i);
	double pEndPlaneAngle=mBigBiteAngle;

	Transform::X_HCS2TCS(pTrack->X0-mPivotXOffset,pTrack->Y0-mPivotYOffset,
		pTrack->Z0-mPivotZOffset, pEndPlaneAngle, pTrack->X0_tr, pTrack->Y0_tr,pTrack->Z0_tr);
	Transform::P_HCS2TCS(pTrack->Theta0, pTrack->Phi0, pEndPlaneAngle, 
		pTrack->Theta0_tr, pTrack->Phi0_tr);

	////////////////////////////////////////////////////////////////////
	//Project|Drift this X0_tr,Y0_tr to target plane
	TVector3 X,P;  //ThreeVector in Hall system for Drift(), in unit of m, GeV
	double pCharge = mPrimaryGeneratorAction->GetPDGCharge(i);  
	double pMass = mPrimaryGeneratorAction->GetPDGMass(i)/GeV;  
	if(!mUseHelmField || fabs(pCharge)<1.0E-8)
	{
		double tmpZtg_tr=0;
		Transform::Project(pTrack->X0_tr, pTrack->Y0_tr,pTrack->Z0_tr, -pTrack->Z0_tr, pTrack->Theta0_tr, 
			pTrack->Phi0_tr,pTrack->Xtg_tr,pTrack->Ytg_tr,tmpZtg_tr);

		pTrack->Thetatg_tr=pTrack->Theta0_tr;
		pTrack->Phitg_tr=pTrack->Phi0_tr;
	}
	else
	{
		//Drift from vertex plane to target plane
		//need to drift forward if Z0_tr<0 or backward if Z0_tr>0

		X.SetXYZ(pTrack->X0/1000.,pTrack->Y0/1000.,pTrack->Z0/1000.);
		P.SetMagThetaPhi(pTrack->P0,pTrack->Theta0,pTrack->Phi0); 

		int ret=0;
		double pZtrLimit=0.0, pTLLimit=1.0;
		if(pTrack->Z0_tr-pZtrLimit<-1.0E-6)
		{
			//forward
			DriftSieve2Tg::Drift2Ztr(X,P,pEndPlaneAngle,pMass, pCharge,pZtrLimit,
				pTLLimit,0.00005,mTargetZOffset/1000);				
		}
		else if(pTrack->Z0_tr-pZtrLimit>1.0E-6)
		{
			//backward
			P*=-1;	//flip P for backward, make sure the particle is positron
			DriftSieve2Tg::Drift2Ztr(X,P,pEndPlaneAngle,pMass,-pCharge,pZtrLimit,
				pTLLimit,0.00005,mTargetZOffset/1000);				
			P*=-1;  //flip it back
		}
		if(ret<0) return false;  

		//now convert the result from HCS to TCS 
		double tmpZtg_tr=0;
		Transform::X_HCS2TCS(X.x()*1000-mPivotXOffset,X.y()*1000-mPivotYOffset, X.z()*1000-mPivotZOffset,
			pEndPlaneAngle, pTrack->Xtg_tr,pTrack->Ytg_tr,tmpZtg_tr);
		Transform::P_HCS2TCS(P.Theta(), P.Phi(), pEndPlaneAngle,pTrack->Thetatg_tr,pTrack->Phitg_tr);

	}

	////////////////////////////////////////////////////////////////////
	//Convert HCS to TCS at VB
	Transform::X_HCS2TCS(pTrack->Xvb-mPivotXOffset,pTrack->Yvb-mPivotYOffset,
		pTrack->Zvb-mPivotZOffset, pEndPlaneAngle, pTrack->Xvb_tr, pTrack->Yvb_tr,pTrack->Zvb_tr);
	Transform::P_HCS2TCS(pTrack->Thetavb, pTrack->Phivb, pEndPlaneAngle, 
		pTrack->Thetavb_tr, pTrack->Phivb_tr);

	//determine the TrackClass
	if( pTrack->Pvb/GeV > 0.15)  pTrack->TrackClass=1;

	if(pTrack->TrackClass==1)
	{
		if( (pTrack->P0 - pTrack->Pvb) / pTrack->P0 < 0.20 ) pTrack->TrackClass=2;
	}

	/*
	//there is no reconstruction right now, the following should be replaced 
	//with the real reconstruction code
	pTrack->Theta_rec_tr=pTrack->Theta0_tr;
	pTrack->Phi_rec_tr=pTrack->Phi0_tr;
	Transform::P_TCS2HCS(pTrack->Theta_rec_tr, pTrack->Phi_rec_tr, pEndPlaneAngle, 
	pTrack->Theta_rec, pTrack->Phi_rec);
	*/
	return true;
}

////////////////////////////////////////////////////////////////////////////////////////////
//Get the effective BPM vertical position when target field is on
//input: true BPM vertical position in meter and particle momentum in GeV
//return: effective BPM vertical position in meter
double HRSRootTree::GetEffBPMVertical(double pX_BPM_tr_m,double pMomentum_GeV)
{
	//Apply the correction to get the effective X_BPM_tr at target plane
	//the following is the result of fitting "[0]/x+[1]" to 
	//(X_proj2tg_tr - Xtg_tr) vs Pvb @ 5.0T target field
	//need to fit for 2.5Ttarget field later
	// EXT PARAMETER                                   STEP         FIRST   
	//NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
	// 1  p0          -6.39438e+00   5.37354e-03   1.74320e-04  -2.30577e-05
	// 2  p1           8.06795e-04   4.92714e-03   1.02122e-05   6.92299e-07
	double X_BPM_tr_eff=pX_BPM_tr_m;	
	if(mUseHelmField)
	{
		X_BPM_tr_eff += (-6.39438/pMomentum_GeV + 8.06795e-04)/1000 * mHelmCurrentRatio;
	}
	return X_BPM_tr_eff;
}


//////////////////////////////////////////////////////////////////////////////////////////////
//do the energy loss correction for the electron, only good for G2P
double HRSRootTree::ElossCorr(double pMomentum_GeV)
{
	if(!mSetupG2PGeometry) return pMomentum_GeV;
	//fit pol3 to P0-Pvb vs P0

	double xx = pMomentum_GeV, pMomentum_corr=pMomentum_GeV,pEloss=0;

	//the following is for 1 mm C12 target, With He4
	//	Minimizer is Linear
	//Chi2                      =      146.358
	//NDf                       =           82
	//p0                        =  0.000998404   +/-   0.000189069 
	//p1                        =  0.000852627   +/-   0.000463301 
	//p2                        =  -0.00028837   +/-   0.000316421 
	//p3                        =  6.29983e-05   +/-   6.30386e-05 
	double pPara_TgType3[]={0.000998404, 0.000852627, -0.00028837, 6.29983e-05};

	//the following is the fit for 1 mm C12 and no He4
	//	Chi2                      =      91.5804
	//NDf                       =           82
	//p0                        =  0.000585605   +/-   0.000165723 
	//p1                        =  0.000206183   +/-   0.000410247 
	//p2                        =  0.000117377   +/-   0.000278337 
	//p3                        = -2.56432e-05   +/-   5.44367e-05 
	double pPara_TgType30[]={0.000585605, 0.000206183, 0.000117377, -2.56432e-05};


	//the following is for NH3 target, With He4
	//NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
	//Minimizer is Linear
	//Chi2                      =      66.6492
	//NDf                       =           82
	//p0                        =   0.00258358   +/-   0.000302646 
	//p1                        =  0.000299926   +/-   0.000748339 
	//p2                        =  0.000687104   +/-   0.000513689 
	//p3                        = -0.000137723   +/-   0.000101663 
	double pPara_TgType1[]={0.00258358, 0.000299926, 0.000687104, -0.000137723};

	double *pPara=pPara_TgType30;
	if(mSetupG2PTarget>=10) 
	{
		//for No He4 
		pPara=pPara_TgType30;
	}
	else if(mSetupG2PTarget>=1)
	{	
		//With He4, not separate by target type any more, just use the target thickness
		if(mTargetL<2.0)  pPara=pPara_TgType3;
		else pPara=pPara_TgType1;
	}

	pEloss = pPara[0] + pPara[1]*xx + pPara[2]*xx*xx + pPara[3]*xx*xx*xx;
	pMomentum_corr += pEloss;
	return pMomentum_corr;
}

//////////////////////////////////////////////////////////////////////////////////////////////
//HRSRootTree::TransportThruHMS() has been rewritten. 
//Here is the 'plane' it goes thrught and the key word of the variables names:
//vertex  -> target-> sieve(vb) -> target     -> focal -> target     -> sieve       -> target    -> vertex
//$0,$0_tr-> $tg_tr-> $vb,$vb_tr->$_proj2tg_tr-> $fp_tr-> $_rec2tg_tr-> $_proj2sl_tr-> $tg_rec_tr-> $_rec,$_rec_tr
//	
//This routine will do the following:

//Section 1: 
//	A) Convert vertex, vb variable from HTC to TCS, 
//	B) Project vertex plane to target plane	
//	C) Determine TrackClass, stop if it less than SkimLevel
//Section 2:
//	A) Project Sieve plane (or VB plane) (_vb_tr) to target plane (_proj2tg_tr)
//	B) Use snake model to transport the _proj2tg_tr to focal plane (_fp), and also
//   reconstruct them to target plane (_rec2tg_tr)
//	C) Stop if fail to call SNAKE model
//Section 3:
//  If the reconstrction is good, set the following tree variables:
//  $fp_tr -> $_rec2tg_tr -> $_proj2sl_tr -> $tg_rec_tr -> $_rec,$_rec_tr
//  A) Stored $fp_tr and $_rec2tg_tr
//  B) If target field is off, set $tg_rec_tr=$rec2tg_tr, then calculate vertexZ, then
//     project from target plane to vertex plane
//  C) If target field is on, do the following:
//     I)  Project $rec2tg_tr to sieve plane: $_proj2sl_tr 
//     II) Call Drift2Ztr() to get them back to target plane: $tg_rec_tr 
//  III)Calculate vertexZ, call Drift2Z() move them from target plane to vertex plane: $_rec,$rec_tr 

bool HRSRootTree::TransportThruHMS(int i)
{
	MyTrack *pTrack=track[i];
	if(pTrack->Pvb/keV<1.0) return false;

	//I do not check if this particle hit HMS or not, 
	//caller should check this
	
	double pEndPlaneAngle=mHMSAngle;
	double pDetectorP0=mHMSMomentum;
	double pTrackL_vb2tg=mPivot2HMSFace;

	//in this code, default pDetectorP0=mHMSMomentum=0;
	//This is used to simulate the whole momentum range [0,8.0) Gev
	//In this case, I would like to set pDetectorP0 = Pvb
	//if one want to check the transportation without considering the delta, just set
	//pDetectorP0 less than 10 Kev
	if(pDetectorP0<1.0E-05) pDetectorP0=pTrack->Pvb;

	// *************************************************************************
	//Section 1: 
	//A) Convert vertex, vb variable from HTC to TCS, 
	//B) Project vertex plane to target plane	
	//C) Determine TrackClass, stop if it less than SkimLevel
	// *************************************************************************

	//Convert HCS to TCS at vertex
	Transform::X_HCS2TCS(pTrack->X0-mPivotXOffset,pTrack->Y0-mPivotYOffset,
		pTrack->Z0-mPivotZOffset, pEndPlaneAngle, pTrack->X0_tr, pTrack->Y0_tr,pTrack->Z0_tr);
	Transform::P_HCS2TCS(pTrack->Theta0, pTrack->Phi0, pEndPlaneAngle, 
		pTrack->Theta0_tr, pTrack->Phi0_tr);

	////////////////////////////////////////////////////////////////////
	//Project|Drift this X0_tr,Y0_tr to target plane
	TVector3 X,P;  //ThreeVector in Hall system for Drift(), in unit of m, GeV
	double pCharge = mPrimaryGeneratorAction->GetPDGCharge(i);  
	double pMass = mPrimaryGeneratorAction->GetPDGMass(i);
	if(!mUseHelmField || fabs(pCharge)<1.0E-8)
	{
		double tmpZtg_tr=0;
		Transform::Project(pTrack->X0_tr, pTrack->Y0_tr,pTrack->Z0_tr, -pTrack->Z0_tr, pTrack->Theta0_tr, 
			pTrack->Phi0_tr,pTrack->Xtg_tr,pTrack->Ytg_tr,tmpZtg_tr);

		pTrack->Thetatg_tr=pTrack->Theta0_tr;
		pTrack->Phitg_tr=pTrack->Phi0_tr;
	}
	else
	{
		//Drift from vertex plane to target plane
		//need to drift forward if Z0_tr<0 or backward if Z0_tr>0

		X.SetXYZ(pTrack->X0/1000.,pTrack->Y0/1000.,pTrack->Z0/1000.);
		P.SetMagThetaPhi(pTrack->P0,pTrack->Theta0,pTrack->Phi0); 

		int ret=0;
		double pZtrLimit=0.0, pTLLimit=1.0;
		if(pTrack->Z0_tr-pZtrLimit<-1.0E-6)
		{
			//forward
			DriftSieve2Tg::Drift2Ztr(X,P,pEndPlaneAngle,pMass, pCharge,pZtrLimit,
				pTLLimit,0.00005,mTargetZOffset/1000);				
		}
		else if(pTrack->Z0_tr-pZtrLimit>1.0E-6)
		{
			//backward
			P*=-1;	//flip P for backward, make sure the particle is positron
			DriftSieve2Tg::Drift2Ztr(X,P,pEndPlaneAngle,pMass,-pCharge,pZtrLimit,
				pTLLimit,0.00005,mTargetZOffset/1000);				
			P*=-1;  //flip it back
		}
		if(ret<0) return false; 

		//now convert the result from HCS to TCS 
		double tmpZtg_tr=0;
		Transform::X_HCS2TCS(X.x()*1000-mPivotXOffset,X.y()*1000-mPivotYOffset, X.z()*1000-mPivotZOffset,
			pEndPlaneAngle, pTrack->Xtg_tr,pTrack->Ytg_tr,tmpZtg_tr);
		Transform::P_HCS2TCS(P.Theta(), P.Phi(), pEndPlaneAngle,pTrack->Thetatg_tr,pTrack->Phitg_tr);
	}

	////////////////////////////////////////////////////////////////////
	//Convert HCS to TCS at VB
	//For HMS, no septum will be used
	DriftSieve2Tg::SetUseSeptumInDrift(0);
	Transform::X_HCS2TCS(pTrack->Xvb-mPivotXOffset,pTrack->Yvb-mPivotYOffset,
		pTrack->Zvb-mPivotZOffset, pEndPlaneAngle, pTrack->Xvb_tr, pTrack->Yvb_tr,pTrack->Zvb_tr);

	Transform::P_HCS2TCS(pTrack->Thetavb, pTrack->Phivb, pEndPlaneAngle, 
		pTrack->Thetavb_tr, pTrack->Phivb_tr);

	////////////////////////////////////////////////////////////////////
	//Determine the TrackClass
	//choose eplice shape
	if(pow((pTrack->Phivb_tr-0.0)/0.035,2.0) + pow((pTrack->Thetavb_tr-0.0)/0.070,2.0) < 1.0) 
	{
		pTrack->TrackClass=1;

		if( (pTrack->P0 - pTrack->Pvb) / pTrack->P0 < 0.10 )
		{
			pTrack->TrackClass=2;

			if( pow((pTrack->Phivb_tr-0.0)/0.028,2.0) + pow((pTrack->Thetavb_tr-0.0)/0.060,2.0) < 1.0 &&
				(pTrack->P0 - pTrack->Pvb) / pTrack->P0 < 0.08 ) pTrack->TrackClass=3;
		}
	}
	//Jixie: to make it run fast, stop tracking if less than skim level
	if (track[0]->TrackClass<iSkimLevel) return false;


	// *************************************************************************
	//Section 2:
	//  A) Project Sieve plane (or VB plane) (_vb_tr) to target plane (_proj2tg_tr)
	//  B) Use snake model to transport the _proj2tg_tr to focal plane (_fp), and also
	//     reconstruct them to target plane (_rec2tg_tr)
	//  C) Stop if fail to call SNAKE model
	// *************************************************************************

	pTrack->Delta = (pTrack->Pvb - pDetectorP0) / pDetectorP0;
	//We know that this particle will not go thrught HRS, stop here to save time
	if(fabs(pTrack->Delta)>0.1) return false;


	// To use the transportation functions (from target plane to focus plane), 
	// we have to project back to the target plane (Z_tg=0), all length in unit of mm
	double Z_proj2tg_tr;
	Transform::Project(pTrack->Xvb_tr, pTrack->Yvb_tr, pTrack->Zvb_tr, -pTrack->Zvb_tr, pTrack->Thetavb_tr,
		pTrack->Phivb_tr, pTrack->X_proj2tg_tr, pTrack->Y_proj2tg_tr, Z_proj2tg_tr);


	////////////////////////////////////////////////////////////////////////
	//transport from tg_tr to fp_tr
	//currently there is no HMS transportation package, so I skip this transport 
	//will add it back once it is ready

	bool bGoodParticle=false;
	double pV5_fp[5]={0,0,0,0,0};
	//temporally comment this part out, will add it back later
	//Fill the vector to be used in forward transportation functions 
	//in meter and rad
	//double pV5_proj2tg[5]={pTrack->X_proj2tg_tr/1000.,pTrack->Thetavb_tr,
	//	pTrack->Y_proj2tg_tr/1000.,pTrack->Phivb_tr,pTrack->Delta};
	//bGoodParticle=mHMSTranModel->Forward(pV5_proj2tg, pV5_fp);
	//bool bApplyVDCSmearing=false;
	//if(bApplyVDCSmearing)
	//{
	//	// Resolutions 
	//	double mWireChamberRes_x = 0.0013; //m;
	//	double mWireChamberRes_y = 0.0013; //m;
	//	double mWireChamberRes_theta = 0.0003; //rad;
	//	double mWireChamberRes_phi = 0.0003; //rad;
	//	pV5_fp[0] += mWireChamberRes_x * HRSRand::fGaus();
	//	pV5_fp[2] += mWireChamberRes_y * HRSRand::fGaus();
	//	pV5_fp[1] += mWireChamberRes_theta * HRSRand::fGaus();
	//	pV5_fp[3] += mWireChamberRes_phi * HRSRand::fGaus();
	//}
	bGoodParticle=true;
	if(!bGoodParticle) return false;

	/////////////////////////////////////////////////////////////////////////
	////snake backward
	////Get the smeared BPM vetical position
	//double X_BPM_tr_m=pTrack->Xtg_tr/1000. + mBPMYRes/1000.*HRSRand::fGaus();

	////to get the effective X_BPM_tr, we need to go thrught 2 steps: 
	////1) Use P0 to get X_BPM_tr_AtP0, use it in the snake backward model to get P_rec
	////2) Use P_rec to get X_BPM_tr_eff, then put it into snake backward model for reconstruction

	////Apply the correction to get the effective X_BPM_tr at target plane
	//double X_BPM_tr_eff=GetEffBPMVertical(X_BPM_tr_m, pDetectorP0);
	//pV5_fp[4]=X_BPM_tr_eff;
	//mHMSTranModel->Backward(pV5_fp,pV5_rec);
	//double tmpPrec=pDetectorP0*(1.0+pV5_rec[4]);
	//X_BPM_tr_eff=GetEffBPMVertical(X_BPM_tr_m, tmpPrec);

	//pV5_fp[4]=X_BPM_tr_eff;
	////////////////////////////////////////////////////////////////////////
	//Fill the vector to be used in backward transportation functions 
	double pV5_rec[5]={pTrack->X_proj2tg_tr/1000.,pTrack->Thetavb_tr,
		pTrack->Y_proj2tg_tr/1000.,pTrack->Phivb_tr,pTrack->Delta};
	//mHMSTranModel->Backward(pV5_fp,pV5_rec);



	//////////////////////////////////////////////////////////////////////

	// *************************************************************************
	//Section 3:
	//  If the reconstrction is good, set the following tree variables:
	//  $fp_tr -> $_rec2tg_tr -> $_proj2sl_tr -> $tg_rec_tr -> $_rec,$_rec_tr
	//  A) Stored $fp_tr and $_rec2tg_tr
	//  B) If target field is off, set $tg_rec_tr=$rec2tg_tr, then calculate vertexZ, then
	//     project from target plane to vertex plane
	//  C) If target field is on, do the following:
	//     I)  Project $rec2tg_tr to sieve plane: $_proj2sl_tr 
	//     II) Call Drift2Ztr() to get them back to target plane: $tg_rec_tr 
	//  III)Calculate vertexZ, call Drift2Z() move them from target plane to vertex plane: $_rec,$rec_tr 
	// *************************************************************************


	//Update TrackClass and store focal plane variables
	pTrack->TrackClass+=4;
	// in unit of mm and rad
	pTrack->Xfp_tr=pV5_fp[0]*1000.;
	pTrack->Thetafp_tr=pV5_fp[1];
	pTrack->Yfp_tr=pV5_fp[2]*1000.;
	pTrack->Phifp_tr=pV5_fp[3];

	
	//currently there is no HMS transportation package, so I simply do these smearing  
	double pDetRes_x = 0.001; //m; 
	double pDetRes_y = 0.001; //m; 
	double pDetRes_theta = 0.001; //rad; 
	double pDetRes_phi = 0.001; //rad;
	double pDetRes_delta = 0.002; //rad;
	pV5_rec[0]=pTrack->X_proj2tg_tr/1000.+HRSRand::fGaus(0,pDetRes_x);
	pV5_rec[2]=pTrack->Y_proj2tg_tr/1000.+HRSRand::fGaus(0,pDetRes_y);
	pV5_rec[1]=pTrack->Thetavb_tr+HRSRand::fGaus(0,pDetRes_theta);
	pV5_rec[3]=pTrack->Phivb_tr+HRSRand::fGaus(0,pDetRes_phi);
	pV5_rec[4]=pTrack->Delta+HRSRand::fGaus(0,pDetRes_delta);

	pTrack->X_rec2tg_tr     = pV5_rec[0]*1000.;  //turn into mm
	pTrack->Theta_rec2tg_tr = pV5_rec[1];
	pTrack->Y_rec2tg_tr     = pV5_rec[2]*1000.;
	pTrack->Phi_rec2tg_tr   = pV5_rec[3];
	pTrack->Delta_rec       = pV5_rec[4];

	pTrack->P_rec2tg = pDetectorP0*(1.0+pV5_rec[4]);	

	//By Jixie:
	//People will be very interested to compare the reconstruction results among stop
	//drifting at the following location: target plane, vertex plane and exact z0
	//therefore I introduce a variable mWhereToStopRecon to 
	//specify where to stop reconstruction
	//0 is at target plane, 1 is at vertex plane, 2 is at exact vertex z 

	//ReconstructToExactZ0 is to test how good the uncertainty will be if we know 
	//vertexZ perfectly (HRS always has trouble in vertex z resolution).	


	/////////////////////////////////////////////////
	//If target field off, just calculate the variables
	//When target field is on, linearly project to sieve slit
	//then Drift() or SNAKE_SL2TG() it back to target plane
	/////////////////////////////////////////////////

	if(!mUseHelmField)
	{
		////////////////////////
		//Target field is off ....
		////////////////////////

		//I do not need to deal with $_proj2sl_tr here, just for completeness
		//can comment it out to speed up
		double Z_proj2sl_tr;		
		Transform::Project(pTrack->X_rec2tg_tr,pTrack->Y_rec2tg_tr,0.0,pTrackL_vb2tg,pTrack->Theta_rec2tg_tr,
			pTrack->Phi_rec2tg_tr,pTrack->X_proj2sl_tr,pTrack->Y_proj2sl_tr,Z_proj2sl_tr);

		pTrack->Xtg_rec_tr     = pTrack->X_rec2tg_tr;  
		pTrack->Thetatg_rec_tr = pTrack->Theta_rec2tg_tr; 
		pTrack->Ytg_rec_tr     = pTrack->Y_rec2tg_tr; 
		pTrack->Phitg_rec_tr   = pTrack->Phi_rec2tg_tr; 

		if(mWhereToStopRecon>=3)
		{
			//Reconstruct To BPMXRes/sin(pEndPlaneAngle) smeared Z0 
			const double pZRes_mm = mBPMXRes/sin(pEndPlaneAngle); 
			pTrack->Z_rec = pTrack->Z0 + pZRes_mm*HRSRand::fGaus();
		}
		else if(mWhereToStopRecon==2)
		{
			//ReconstructToExactZ0 
			pTrack->Z_rec = pTrack->Z0;
		}
		else if(mWhereToStopRecon==1)
		{			
			double Y_BPM_tr_mm=pTrack->X0 + mBPMXRes*HRSRand::fGaus();
			pTrack->Z_rec = -pTrack->Ytg_rec_tr*cos(pTrack->Phitg_rec_tr) / sin(pEndPlaneAngle+pTrack->Phitg_rec_tr) + 
				Y_BPM_tr_mm/tan(pEndPlaneAngle+pTrack->Phitg_rec_tr);  //in unit of mm
			pTrack->Z_rec += mTargetZOffset;
		}
		else if(mWhereToStopRecon==0) pTrack->Z_rec = mTargetZOffset;

		//NO field, do linear projection to vertex plane
		double pVZ2TargetPlane=(pTrack->Z_rec-mTargetZOffset)*cos(pEndPlaneAngle);

		//The snake model has the target at the target center, Got confirmed by Min and John 
		double tmpZ_tg_snake=0;  //will be useful in case there is an offset
		Transform::Project(pTrack->Xtg_rec_tr,pTrack->Ytg_rec_tr,tmpZ_tg_snake,pVZ2TargetPlane-tmpZ_tg_snake,
			pTrack->Thetatg_rec_tr,pTrack->Phitg_rec_tr,pTrack->X_rec_tr,pTrack->Y_rec_tr,pTrack->Z_rec_tr);

		Transform::X_TCS2HCS(pTrack->X_rec_tr, pTrack->Y_rec_tr, pTrack->Z_rec_tr, pEndPlaneAngle, 
			pTrack->X_rec, pTrack->Y_rec, pTrack->Z_rec);
		pTrack->X_rec+=this->mTargetXOffset;
		pTrack->Y_rec+=this->mTargetYOffset;
		pTrack->Z_rec+=this->mTargetZOffset; 

		//This is my recostruction for Theta_rec and Phi_rec
		pTrack->Theta_rec_tr = pTrack->Thetatg_rec_tr; 
		pTrack->Phi_rec_tr = pTrack->Phitg_rec_tr; 
		Transform::P_TCS2HCS(pTrack->Theta_rec_tr, pTrack->Phi_rec_tr, pEndPlaneAngle, 
			pTrack->Theta_rec, pTrack->Phi_rec);
	}
	else
	{
		////////////////////////
		//Target field is on ...
		////////////////////////

		//By Jixie: HMS reconstruction only reconstruct to the target plane
		//Here I project it to the sieve slit plane, which is pTrackL_vb2tg, going forward
		//Note that pTrackL_vb2tg has been set above
		double Z_proj2sl_tr=0.0;
		Transform::Project(pTrack->X_rec2tg_tr,pTrack->Y_rec2tg_tr,0.0,pTrackL_vb2tg,pTrack->Theta_rec2tg_tr,
			pTrack->Phi_rec2tg_tr,pTrack->X_proj2sl_tr,pTrack->Y_proj2sl_tr,Z_proj2sl_tr);


		//Now reconstruct from sieve slit plane to target plane using Drift package

#if defined DEBUG_HRS_RECONSTRUCTION && (DEBUG_HRS_RECONSTRUCTION>=6)
		//use the exact simulated values to debug this drift routine
		X.SetXYZ(pTrack->Xvb/1000.,pTrack->Yvb/1000.,pTrack->Zvb/1000.); 
		P.SetMagThetaPhi(pTrack->Pvb,pTrack->Thetavb,pTrack->Phivb);
#else
		//convert proj2sl variables from TCS to HCS, 
		//should compatible for both non-septum and SeptumPlusSTDHRS situation
		double tmpX,tmpY,tmpZ,tmpTheta,tmpPhi;	
		Transform::X_TCS2HCS(pTrack->X_proj2sl_tr,pTrack->Y_proj2sl_tr,Z_proj2sl_tr,
			pEndPlaneAngle,tmpX,tmpY,tmpZ);
		Transform::P_TCS2HCS(pTrack->Theta_rec2tg_tr,pTrack->Phi_rec2tg_tr,
			pEndPlaneAngle,tmpTheta,tmpPhi);
		//convert to meter 
		X.SetXYZ((tmpX+mPivotXOffset)/1000.,(tmpY+mPivotYOffset)/1000.,(tmpZ+mPivotZOffset)/1000.);
		P.SetMagThetaPhi(pTrack->P_rec2tg,tmpTheta,tmpPhi); 
#endif


		//Drift from SL to Target, must be in Hall Coordinate system and in meter, GeV
		//backward
		double tmpZtrLimit=0.0;  
		P*=-1.0;	//flip momentum	for backward drift
		int ret=DriftSieve2Tg::Drift2Ztr(X,P,pEndPlaneAngle,pMass,-pCharge,tmpZtrLimit,
			2.0*pTrackL_vb2tg/1000,0.00001,mTargetZOffset/1000.0);				
		P*=-1.0;	//flip momentum	for forward drift
		//if something wrong with the drift, stop
		if(ret<0)  return false;


		//now convert the result from HCS to TCS and store them into pV5tg_rec_tr, in unit of m
		double tmpZtg_tr_mm;
		Transform::X_HCS2TCS(X.x()*1000-mPivotXOffset,X.y()*1000-mPivotYOffset, X.z()*1000-mPivotZOffset,
			pEndPlaneAngle, pTrack->Xtg_rec_tr, pTrack->Ytg_rec_tr, tmpZtg_tr_mm);
		Transform::P_HCS2TCS(P.Theta(), P.Phi(), pEndPlaneAngle, pTrack->Thetatg_rec_tr, pTrack->Phitg_rec_tr);


		////////////////////////////////////////////////////////
		//Reconstruction to the target plane is completed. Now
		//reconstruct theta, phi, x, y, z at vertex in HCS and TCS
		//
		//The snake always reconstruct to target plane
		//Now I have to calculate vertexZ 
		//then Drift from target plane to vertexZ plane 
		////////////////////////////////////////////////////////

		if(mWhereToStopRecon>=1)
		{
			//mWhereToStopRecon=0 is stop drifting at target plane, 1 is at vertex plane, 
			//2 is at exact vertex z0, 3 is exact z0 smeared with BPMXRes/sin(pEndPlaneAngle)

			//Get Vertex Z
			if(mWhereToStopRecon>=3)
			{
				//erconstruct to smeared thrown Z0
				//Reconstruct To BPMXRes/sin(pEndPlaneAngle) smeared Z0 
				const double pHRSZRes_mm = mBPMXRes/sin(pEndPlaneAngle); 
				pTrack->Z_rec=pTrack->Z0 + pHRSZRes_mm*HRSRand::fGaus();
			}
			else if(mWhereToStopRecon==2)
			{
				//Reconstruct To Exact Z0 
				pTrack->Z_rec = pTrack->Z0;
			}
			else
			{
				//reconstruct to V5tg_rec_tr predicted Z0
				//This line is from Nilanga's optics note, page 9
				//It has been debugged
				//double Y_BPM_tr_mm=pTrack->X0;
				double Y_BPM_tr_mm=pTrack->X0 + mBPMXRes*HRSRand::fGaus();
				pTrack->Z_rec = -pTrack->Ytg_rec_tr*cos(pTrack->Phitg_rec_tr) / sin(pEndPlaneAngle+pTrack->Phitg_rec_tr) + 
					Y_BPM_tr_mm/tan(pEndPlaneAngle+pTrack->Phitg_rec_tr);  //in unit of mm
				pTrack->Z_rec+=mTargetZOffset;
			}

			//////////////////////////////////////////
			//Drift from target plane to ideal vertexZ  
			//this part is complicate: if Z_rec is at downstream of target plane, we have to drift forward
			//otherwise drift backward

			if(pTrack->Z_rec/1000 > X.z())
			{
				//forward
				//set the tracklength limit to 0.2 m to save time
				//P is already ready for forward, make sure the particle is electron
				DriftSieve2Tg::Drift2Z(X,P,pMass,pCharge,pTrack->Z_rec/1000,0.2,0.0001);				
			}
			else if(pTrack->Z_rec/1000 < X.z())
			{
				//backward
				//set the tracklength limit to 0.2 m to save time
				P*=-1;	//flip P for backward, make sure the particle is positron
				DriftSieve2Tg::Drift2Z(X,P,pMass,-pCharge,pTrack->Z_rec/1000,0.2,0.0001);	
				P*=-1;  //flip it back
			}

		}

		//now convert the result from HCS to TCS, converted to mm at the same time 
		Transform::X_HCS2TCS(X.x()*1000-mPivotXOffset,X.y()*1000-mPivotYOffset, X.z()*1000-mPivotZOffset,
			pEndPlaneAngle, pTrack->X_rec_tr,pTrack->Y_rec_tr,pTrack->Z_rec_tr);
		Transform::P_HCS2TCS(P.Theta(), P.Phi(), pEndPlaneAngle,pTrack->Theta_rec_tr,pTrack->Phi_rec_tr);

		//////////////////////////////////////////
		//now store the vertex plane variables in mm

		pTrack->X_rec=X.x()*1000;
		pTrack->Y_rec=X.y()*1000;
		pTrack->Z_rec=X.z()*1000;  
		pTrack->Theta_rec=P.Theta();
		pTrack->Phi_rec=P.Phi();

	}


	////////////////////////////////////////////////////////
	//Energy loss correction
	////////////////////////////////////////////////////////	

	//pTrack->P_rec = ElossCorr(pTrack->P_rec2tg);
	//Delta is defined by (Pvb-pDetectorP0)/pDetectorP0
	//therefore I do not want to correct Delta_rec after Eloss corr.
	//pTrack->Delta_rec = (pTrack->P_rec-pDetectorP0)/pDetectorP0;

	return bGoodParticle;
}



//////////////////////////////////////////////////////////////////////////////////////////////
//HRSRootTree::TransportThruHRS() has been rewritten. 
//Here is the 'plane' it goes thrught and the key word of the variables names:
//vertex  -> target-> sieve(vb) -> target     -> focal -> target     -> sieve       -> target    -> vertex
//$0,$0_tr-> $tg_tr-> $vb,$vb_tr->$_proj2tg_tr-> $fp_tr-> $_rec2tg_tr-> $_proj2sl_tr-> $tg_rec_tr-> $_rec,$_rec_tr
//	
//This routine will do the following:

//Section 1: 
//	A) Convert vertex, vb variable from HTC to TCS, 
//	B) Project vertex plane to target plane	
//	C) Determine TrackClass, stop if it less than SkimLevel
//Section 2:
//	A) Project Sieve plane (or VB plane) (_vb_tr) to target plane (_proj2tg_tr)
//	B) Use snake model to transport the _proj2tg_tr to focal plane (_fp), and also
//   reconstruct them to target plane (_rec2tg_tr)
//	C) Stop if fail to call SNAKE model
//Section 3:
//  If the reconstrction is good, set the following tree variables:
//  $fp_tr -> $_rec2tg_tr -> $_proj2sl_tr -> $tg_rec_tr -> $_rec,$_rec_tr
//  A) Stored $fp_tr and $_rec2tg_tr
//  B) If target field is off, set $tg_rec_tr=$rec2tg_tr, then calculate vertexZ, then
//     project from target plane to vertex plane
//  C) If target field is on, do the following:
//     I)  Project $rec2tg_tr to sieve plane: $_proj2sl_tr 
//     II) Call Drift2Ztr() to get them back to target plane: $tg_rec_tr 
//  III)Calculate vertexZ, call Drift2Z() move them from target plane to vertex plane: $_rec,$rec_tr 

bool HRSRootTree::TransportThruHRS(int i)
{
	MyTrack *pTrack=track[i];
	if(pTrack->Pvb/keV<1.0) return false;
/*
	//I do not check if this particle hit HRS or not, 
	//caller should check this
	//Need to identify if this particle goes into hrs or 3rd arm or bigbite
	//When target field is on, it is hard to use angle at target to tell
	//Use VB position might be better
	//bigbite was placed at 72 degrees

	if(fabs(pTrack->Theta0 - mBigBiteAngle)/deg <15.0)  return false;
	double pVBAngle=atan2(pTrack->Xvb-mPivotXOffset, pTrack->Zvb-mPivotZOffset);
	if(fabs(pVBAngle - mBigBiteAngle)/deg <15.0)  return false;
	if(mThirdArmAngle/deg>15 && fabs(pVBAngle - mThirdArmAngle)/deg <15.0)  return false;
*/

	int    pIsLeftArm=(pTrack->Xvb>0)?1:0; //left arm or right arm
	double pEndPlaneAngle=(pIsLeftArm)?mLSeptumAngle:mRSeptumAngle;
	double pHRSMomentum=(pIsLeftArm)?mLHRSMomentum:mRHRSMomentum;
	double pTrackL_vb2tg=(pIsLeftArm)?mPivot2LHRSVBFace:mPivot2RHRSVBFace;

	//in this code, default mLHRSMomentum=mRHRSMomentum=0;
	//This is used to simulate the whole momentum range [0,4.0) Gev
	//In this case, I would like to set pHRSMomentum = Pvb
	//if one want to check the transportation without considering the delta, just set
	//pHRSMomentum less than 10 Kev
	if(pHRSMomentum<1.0E-05) pHRSMomentum=pTrack->Pvb;

	// *************************************************************************
	//Section 1: 
	//A) Convert vertex, vb variable from HTC to TCS, 
	//B) Project vertex plane to target plane	
	//C) Determine TrackClass, stop if it less than SkimLevel
	// *************************************************************************

	//Convert HCS to TCS at vertex
	Transform::X_HCS2TCS(pTrack->X0-mPivotXOffset,pTrack->Y0-mPivotYOffset,
		pTrack->Z0-mPivotZOffset, pEndPlaneAngle, pTrack->X0_tr, pTrack->Y0_tr,pTrack->Z0_tr);
	Transform::P_HCS2TCS(pTrack->Theta0, pTrack->Phi0, pEndPlaneAngle, 
		pTrack->Theta0_tr, pTrack->Phi0_tr);

	////////////////////////////////////////////////////////////////////
	//Project|Drift this X0_tr,Y0_tr to target plane
	TVector3 X,P;  //ThreeVector in Hall system for Drift(), in unit of m, GeV
	double pCharge = mPrimaryGeneratorAction->GetPDGCharge(i);  
	double pMass = mPrimaryGeneratorAction->GetPDGMass(i);
	if(!mUseHelmField || abs(pCharge)<1.0E-8)
	{
		double tmpZtg_tr=0;
		Transform::Project(pTrack->X0_tr, pTrack->Y0_tr,pTrack->Z0_tr, -pTrack->Z0_tr, pTrack->Theta0_tr, 
			pTrack->Phi0_tr,pTrack->Xtg_tr,pTrack->Ytg_tr,tmpZtg_tr);

		pTrack->Thetatg_tr=pTrack->Theta0_tr;
		pTrack->Phitg_tr=pTrack->Phi0_tr;
	}
	else
	{
		//Drift from vertex plane to target plane
		//need to drift forward if Z0_tr<0 or backward if Z0_tr>0

		X.SetXYZ(pTrack->X0/1000.,pTrack->Y0/1000.,pTrack->Z0/1000.);
		P.SetMagThetaPhi(pTrack->P0,pTrack->Theta0,pTrack->Phi0); 

		int ret=0;
		double pZtrLimit=0.0, pTLLimit=1.0;
		if(pTrack->Z0_tr-pZtrLimit<-1.0E-6)
		{
			//forward
			DriftSieve2Tg::Drift2Ztr(X,P,pEndPlaneAngle,pMass, pCharge,pZtrLimit,
				pTLLimit,0.00005,mTargetZOffset/1000);				
		}
		else if(pTrack->Z0_tr-pZtrLimit>1.0E-6)
		{
			//backward
			P*=-1;	//flip P for backward, make sure the particle is positron
			DriftSieve2Tg::Drift2Ztr(X,P,pEndPlaneAngle,pMass,-pCharge,pZtrLimit,
				pTLLimit,0.00005,mTargetZOffset/1000);				
			P*=-1;  //flip it back
		}
		if(ret<0) return false; 

		//now convert the result from HCS to TCS 
		double tmpZtg_tr=0;
		Transform::X_HCS2TCS(X.x()*1000-mPivotXOffset,X.y()*1000-mPivotYOffset, X.z()*1000-mPivotZOffset,
			pEndPlaneAngle, pTrack->Xtg_tr,pTrack->Ytg_tr,tmpZtg_tr);
		Transform::P_HCS2TCS(P.Theta(), P.Phi(), pEndPlaneAngle,pTrack->Thetatg_tr,pTrack->Phitg_tr);
	}

	////////////////////////////////////////////////////////////////////
	//Convert HCS to TCS at VB
	if(mUseSeptumField && mUseSeptumPlusStdHRS) 
	{
		//In this case I will use the septum field to propagate particles and place the Virtual Boundary at
		//the Q1 entrance, therefore the end plane angle is HRS angle and the origin is hall center
		//I should project the track back to the hall center, not the pivot any more. 
		//To do that I should not apply the mPivotZOffset
		DriftSieve2Tg::SetUseSeptumInDrift(1);
		pEndPlaneAngle=(pIsLeftArm)?mLHRSAngle:mRHRSAngle;
		Transform::X_HCS2TCS(pTrack->Xvb-mPivotXOffset,pTrack->Yvb-mPivotYOffset,
			pTrack->Zvb, pEndPlaneAngle, pTrack->Xvb_tr, pTrack->Yvb_tr,pTrack->Zvb_tr);
	}
	else
	{
		DriftSieve2Tg::SetUseSeptumInDrift(0);
		Transform::X_HCS2TCS(pTrack->Xvb-mPivotXOffset,pTrack->Yvb-mPivotYOffset,
			pTrack->Zvb-mPivotZOffset, pEndPlaneAngle, pTrack->Xvb_tr, pTrack->Yvb_tr,pTrack->Zvb_tr);
	}

	Transform::P_HCS2TCS(pTrack->Thetavb, pTrack->Phivb, pEndPlaneAngle, 
		pTrack->Thetavb_tr, pTrack->Phivb_tr);

	////////////////////////////////////////////////////////////////////
	//Determine the TrackClass
	//choose eplice shape
	if(pow((pTrack->Phivb_tr-0.0)/0.035,2.0) + pow((pTrack->Thetavb_tr-0.0)/0.060,2.0) < 1.0) 
	{
		pTrack->TrackClass=1;

		if( (pTrack->P0 - pTrack->Pvb) / pTrack->P0 < 0.05 ) 
		{
			pTrack->TrackClass=2;

			if( pow((pTrack->Phivb_tr-0.0)/0.028,2.0) + pow((pTrack->Thetavb_tr-0.0)/0.050,2.0) < 1.0 &&
				(pTrack->P0 - pTrack->Pvb) / pTrack->P0 < 0.04 ) pTrack->TrackClass=3;
		}
	}
	//Jixie: to make it run fast, stop tracking if less than skim level
	if (track[0]->TrackClass<iSkimLevel) return false;


	// *************************************************************************
	//Section 2:
	//  A) Project Sieve plane (or VB plane) (_vb_tr) to target plane (_proj2tg_tr)
	//  B) Use snake model to transport the _proj2tg_tr to focal plane (_fp), and also
	//     reconstruct them to target plane (_rec2tg_tr)
	//  C) Stop if fail to call SNAKE model
	// *************************************************************************

	pTrack->Delta = (pTrack->Pvb - pHRSMomentum) / pHRSMomentum;
	//We know that this particle will not go thrught HRS, stop here to save time
	if(fabs(pTrack->Delta)>0.1) return false;


	// To use the transportation functions (from target plane to focus plane), 
	// we have to project back to the target plane (Z_tg=0), all length in unit of mm
	double Z_proj2tg_tr;
	Transform::Project(pTrack->Xvb_tr, pTrack->Yvb_tr, pTrack->Zvb_tr, -pTrack->Zvb_tr, pTrack->Thetavb_tr,
		pTrack->Phivb_tr, pTrack->X_proj2tg_tr, pTrack->Y_proj2tg_tr, Z_proj2tg_tr);


	//Fill the vector to be used in forward transportation functions 
	//in meter and rad
	double pV5_rec[5]={pTrack->X_proj2tg_tr/1000.,pTrack->Thetavb_tr,
		pTrack->Y_proj2tg_tr/1000.,pTrack->Phivb_tr,pTrack->Delta};
	double pV5_fp[5]={0,0,0,0,0};

#if defined DEBUG_HRS_RECONSTRUCTION && (DEBUG_HRS_RECONSTRUCTION>=5)
	//debug the routine, set the initial at (0,0,0,0,0)
	pV5_rec[0]=0;
	pV5_rec[1]=0;
	pV5_rec[2]=0;
	pV5_rec[3]=0;
	pV5_rec[4]=pTrack->Delta;
#endif

	//G2P transverse setting
	const int kG2pHelmAngle=90, kGepHelmAngle=6;
	int iFieldRotation=90;
	if (fabs(sin(mHelmRotAngle1)-sin(270.*deg))<0.001) iFieldRotation=kG2pHelmAngle;
	else if (fabs(sin(mHelmRotAngle1)-sin(6.*deg))<0.001) iFieldRotation=kGepHelmAngle;
	else if (fabs(sin(mHelmRotAngle1)-sin(0.*deg))<0.001) iFieldRotation=0;


	bool bGoodParticle=false;
	//need to set the HRS arm and angle, 
	mHRSTranModel->SetArm((pIsLeftArm>0)?true:false); 
	mHRSTranModel->SetHRSAngle(pEndPlaneAngle); 

	//////////////////////////////////////////////////////////////////////
	//input are in m or rad, will convert angle into tan(angle) inside 
	//and also convert tan(angle) back to angle before return
	//bool SNAKE::SNAKEThruHRS(HRSTransport* pSnakeModel,double pX0_tg_m,
	//			  const double *V5_tg,double* pV5_fp, double* pV5_rec);
	//reconstructed values stored in pV5_rec

	double pV5_proj2tg[5]={pTrack->X_proj2tg_tr/1000.,pTrack->Thetavb_tr,
		pTrack->Y_proj2tg_tr/1000.,pTrack->Phivb_tr,pTrack->Delta};

	bGoodParticle=mHRSTranModel->Forward(pV5_proj2tg, pV5_fp);
	if(!bGoodParticle) return false;

	bool bApplyVDCSmearing=true;
	if(bApplyVDCSmearing)
	{
		// Resolutions 
		double mWireChamberRes_x = 0.0013; //m;
		double mWireChamberRes_y = 0.0013; //m;
		double mWireChamberRes_theta = 0.0003; //rad;
		double mWireChamberRes_phi = 0.0003; //rad;
		pV5_fp[0] += mWireChamberRes_x * HRSRand::fGaus();
		pV5_fp[2] += mWireChamberRes_y * HRSRand::fGaus();
		pV5_fp[1] += mWireChamberRes_theta * HRSRand::fGaus();
		pV5_fp[3] += mWireChamberRes_phi * HRSRand::fGaus();
	}

	///////////////////////////////////////////////////////////////////////
	//Get the smeared BPM vetical position
	double X_BPM_tr_m=pTrack->Xtg_tr/1000. + mBPMYRes/1000.*HRSRand::fGaus();

	//to get the effective X_BPM_tr, we need to go thrught 2 steps: 
	//1) Use P0 to get X_BPM_tr_AtP0, use it in the snake backward model to get P_rec
	//2) Use P_rec to get X_BPM_tr_eff, then put it into snake backward model for reconstruction

	//Apply the correction to get the effective X_BPM_tr at target plane
	double X_BPM_tr_eff=GetEffBPMVertical(X_BPM_tr_m, pHRSMomentum);
	pV5_fp[4]=X_BPM_tr_eff;
	mHRSTranModel->Backward(pV5_fp,pV5_rec);
	double tmpPrec=pHRSMomentum*(1.0+pV5_rec[4]);
	X_BPM_tr_eff=GetEffBPMVertical(X_BPM_tr_m, tmpPrec);

	pV5_fp[4]=X_BPM_tr_eff;
	mHRSTranModel->Backward(pV5_fp,pV5_rec);


	//Reconstruct use analyzer database
	double pV5_rec_db[5];
	if (mUseOpticsDB==1) 
	{
		if (pIsLeftArm) mRecUseDBL->CalcTargetCoords(pV5_fp, pV5_rec_db);
		else mRecUseDBR->CalcTargetCoords(pV5_fp, pV5_rec_db);
	}

	//////////////////////////////////////////////////////////////////////

	// *************************************************************************
	//Section 3:
	//  If the reconstrction is good, set the following tree variables:
	//  $fp_tr -> $_rec2tg_tr -> $_proj2sl_tr -> $tg_rec_tr -> $_rec,$_rec_tr
	//  A) Stored $fp_tr and $_rec2tg_tr
	//  B) If target field is off, set $tg_rec_tr=$rec2tg_tr, then calculate vertexZ, then
	//     project from target plane to vertex plane
	//  C) If target field is on, do the following:
	//     I)  Project $rec2tg_tr to sieve plane: $_proj2sl_tr 
	//     II) Call Drift2Ztr() to get them back to target plane: $tg_rec_tr 
	//  III)Calculate vertexZ, call Drift2Z() move them from target plane to vertex plane: $_rec,$rec_tr 
	// *************************************************************************


	//Update TrackClass and store focal plane variables
	pTrack->TrackClass+=4;
	// in unit of mm and rad
	pTrack->Xfp_tr=pV5_fp[0]*1000.;
	pTrack->Thetafp_tr=pV5_fp[1];
	pTrack->Yfp_tr=pV5_fp[2]*1000.;
	pTrack->Phifp_tr=pV5_fp[3];

	if (mUseOpticsDB==1) {
		pTrack->Xtg_rec_db_tr = pV5_rec_db[0]*1000.;
		pTrack->Thetatg_rec_db_tr = pV5_rec_db[1];
		pTrack->Ytg_rec_db_tr = pV5_rec_db[2]*1000.;
		pTrack->Phitg_rec_db_tr = pV5_rec_db[3];
		pTrack->Delta_rec_db = pV5_rec_db[4];
	}

#if defined DEBUG_HRS_RECONSTRUCTION && (DEBUG_HRS_RECONSTRUCTION>=4)
	//debug the routine, set the snake reconstructed result	as the proj2tg result
	//to remove snake resolution, this will be used to determine the resulution due to
	//my own Drift() routine or vertex Z effect 
	pV5_rec[0]=pTrack->X_proj2tg_tr/1000.;
	pV5_rec[1]=pTrack->Thetavb_tr;
	pV5_rec[2]=pTrack->Y_proj2tg_tr/1000.;
	pV5_rec[2]=pTrack->Phivb_tr;
	pV5_rec[4]=pTrack->Delta;
#endif


	pTrack->X_rec2tg_tr     = pV5_rec[0]*1000.;  //turn into mm
	pTrack->Theta_rec2tg_tr = pV5_rec[1];
	pTrack->Y_rec2tg_tr     = pV5_rec[2]*1000.;
	pTrack->Phi_rec2tg_tr   = pV5_rec[3];
	pTrack->Delta_rec       = pV5_rec[4];

	pTrack->P_rec2tg = pHRSMomentum*(1.0+pV5_rec[4]);	

	//By Jixie:
	//People will be very interested to compare the reconstruction results among stop
	//drifting at the following location: target plane, vertex plane and exact z0
	//therefore I introduce a variable mWhereToStopRecon to 
	//specify where to stop reconstruction
	//0 is at target plane, 1 is at vertex plane, 2 is at exact vertex z 

	//ReconstructToExactZ0 is to test how good the uncertainty will be if we know 
	//vertexZ perfectly (HRS always has trouble in vertex z resolution).	


	/////////////////////////////////////////////////
	//If target field off, just calculate the variables
	//When target field is on, linearly project to sieve slit
	//then Drift() or SNAKE_SL2TG() it back to target plane
	/////////////////////////////////////////////////

	if(!mUseHelmField)
	{
		////////////////////////
		//Target field is off ....
		////////////////////////

		//I do not need to deal with $_proj2sl_tr here, just for completeness
		//can comment it out to speed up
		double Z_proj2sl_tr;		
		Transform::Project(pTrack->X_rec2tg_tr,pTrack->Y_rec2tg_tr,0.0,pTrackL_vb2tg,pTrack->Theta_rec2tg_tr,
			pTrack->Phi_rec2tg_tr,pTrack->X_proj2sl_tr,pTrack->Y_proj2sl_tr,Z_proj2sl_tr);

		pTrack->Xtg_rec_tr     = pTrack->X_rec2tg_tr;  
		pTrack->Thetatg_rec_tr = pTrack->Theta_rec2tg_tr; 
		pTrack->Ytg_rec_tr     = pTrack->Y_rec2tg_tr; 
		pTrack->Phitg_rec_tr   = pTrack->Phi_rec2tg_tr; 

		if(mWhereToStopRecon>=3)
		{
			//Reconstruct To BPMXRes/sin(pEndPlaneAngle) smeared Z0 
			const double pHRSZRes_mm = mBPMXRes/sin(pEndPlaneAngle); 
			pTrack->Z_rec = pTrack->Z0 + pHRSZRes_mm*HRSRand::fGaus();
		}
		else if(mWhereToStopRecon==2)
		{
			//ReconstructToExactZ0 
			pTrack->Z_rec = pTrack->Z0;
		}
		else if(mWhereToStopRecon==1)
		{
			/////////////////////////////////////////////////////////////
			//This part is the debug result of Z_rec calculation
			//this line is from Nilanga's optics note, page 9
			//The tested result shows this method have a resolution of 10.8mm 
			//NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
			//	1  Constant     2.72251e+02   6.09607e+00  -2.74515e-03   4.19840e-07
			//	2  Mean        -5.61764e+00   2.45684e-01  -1.17190e-05   5.89404e-06
			//	3  Sigma        1.07950e+01   3.41003e-01   1.69630e-06  -1.33610e-05

			//pTrack->Z_rec = -pTrack->Ytg_rec_tr*cos(pTrack->Phitg_rec_tr) / sin(pEndPlaneAngle+pTrack->Phitg_rec_tr) + 
			//	pTrack->X0/tan(pEndPlaneAngle+pTrack->Phitg_rec_tr);  //in unit of mm
			//pTrack->X_rec+=this->mTargetXOffset;
			//pTrack->Y_rec+=this->mTargetYOffset;
			//pTrack->Z_rec+=this->mTargetZOffset;  
			///////////////////////////////////////////////////////////////

			//this line is from Nilanga's optics note, page 9
			//By Jixie: This calculation does not work if target field is on,  do not know how bad it is yet
			//It has been debugged
			//double Y_BPM_tr_mm=pTrack->X0;
			double Y_BPM_tr_mm=pTrack->X0 + mBPMXRes*HRSRand::fGaus();
			pTrack->Z_rec = -pTrack->Ytg_rec_tr*cos(pTrack->Phitg_rec_tr) / sin(pEndPlaneAngle+pTrack->Phitg_rec_tr) + 
				Y_BPM_tr_mm/tan(pEndPlaneAngle+pTrack->Phitg_rec_tr);  //in unit of mm
			pTrack->Z_rec += mTargetZOffset;
		}
		else if(mWhereToStopRecon==0) pTrack->Z_rec = mTargetZOffset;

		//NO field, do linear projection to vertex plane
		double pVZ2TargetPlane=(pTrack->Z_rec-mTargetZOffset)*cos(pEndPlaneAngle);

		//The snake model has the target at the target center, Got confirmed by Min and John 
		double tmpZ_tg_snake=0;  //will be useful in case there is an offset
		Transform::Project(pTrack->Xtg_rec_tr,pTrack->Ytg_rec_tr,tmpZ_tg_snake,pVZ2TargetPlane-tmpZ_tg_snake,
			pTrack->Thetatg_rec_tr,pTrack->Phitg_rec_tr,pTrack->X_rec_tr,pTrack->Y_rec_tr,pTrack->Z_rec_tr);

		Transform::X_TCS2HCS(pTrack->X_rec_tr, pTrack->Y_rec_tr, pTrack->Z_rec_tr, pEndPlaneAngle, 
			pTrack->X_rec, pTrack->Y_rec, pTrack->Z_rec);
		pTrack->X_rec+=this->mTargetXOffset;
		pTrack->Y_rec+=this->mTargetYOffset;
		pTrack->Z_rec+=this->mTargetZOffset; 

		//This is my recostruction for Theta_rec and Phi_rec
		pTrack->Theta_rec_tr = pTrack->Thetatg_rec_tr; 
		pTrack->Phi_rec_tr = pTrack->Phitg_rec_tr; 
		Transform::P_TCS2HCS(pTrack->Theta_rec_tr, pTrack->Phi_rec_tr, pEndPlaneAngle, 
			pTrack->Theta_rec, pTrack->Phi_rec);
	}
	else
	{
		////////////////////////
		//Target field is on ...
		////////////////////////

		//By Jixie: HRS reconstruction only reconstruct to the target plane
		//Here I project it to the sieve slit plane, which is pTrackL_vb2tg=806 mm forward
		//Note that pTrackL_vb2tg has been set to 806.0 mm above
		double Z_proj2sl_tr=0.0;
		Transform::Project(pTrack->X_rec2tg_tr,pTrack->Y_rec2tg_tr,0.0,pTrackL_vb2tg,pTrack->Theta_rec2tg_tr,
			pTrack->Phi_rec2tg_tr,pTrack->X_proj2sl_tr,pTrack->Y_proj2sl_tr,Z_proj2sl_tr);


		//Now reconstruct from sieve slit plane to target plane using Drift package

#if defined DEBUG_HRS_RECONSTRUCTION && (DEBUG_HRS_RECONSTRUCTION>=6)
		//use the exact simulated values to debug this drift routine
		X.SetXYZ(pTrack->Xvb/1000.,pTrack->Yvb/1000.,pTrack->Zvb/1000.); 
		P.SetMagThetaPhi(pTrack->Pvb,pTrack->Thetavb,pTrack->Phivb);
#else
		//convert proj2sl variables from TCS to HCS, 
		//should compatible for both non-septum and SeptumPlusSTDHRS situation
		double tmpX,tmpY,tmpZ,tmpTheta,tmpPhi;	
		Transform::X_TCS2HCS(pTrack->X_proj2sl_tr,pTrack->Y_proj2sl_tr,Z_proj2sl_tr,
			pEndPlaneAngle,tmpX,tmpY,tmpZ);
		Transform::P_TCS2HCS(pTrack->Theta_rec2tg_tr,pTrack->Phi_rec2tg_tr,
			pEndPlaneAngle,tmpTheta,tmpPhi);
		//convert to meter 
		X.SetXYZ((tmpX+mPivotXOffset)/1000.,(tmpY+mPivotYOffset)/1000.,(tmpZ+mPivotZOffset)/1000.);
		P.SetMagThetaPhi(pTrack->P_rec2tg,tmpTheta,tmpPhi); 
#endif


		//Drift from SL to Target, must be in Hall Coordinate system and in meter, GeV
		//backward
		double tmpZtrLimit=0.0;  
		P*=-1.0;	//flip momentum	for backward drift
		int ret=DriftSieve2Tg::Drift2Ztr(X,P,pEndPlaneAngle,pMass,-pCharge,tmpZtrLimit,
			2.0*pTrackL_vb2tg/1000,0.00001,mTargetZOffset/1000.0);				
		P*=-1.0;	//flip momentum	for forward drift
		//if something wrong with the drift, stop
		if(ret<0)  return false;


		//now convert the result from HCS to TCS and store them into pV5tg_rec_tr, in unit of m
		double tmpZtg_tr_mm;
		Transform::X_HCS2TCS(X.x()*1000-mPivotXOffset,X.y()*1000-mPivotYOffset, X.z()*1000-mPivotZOffset,
			pEndPlaneAngle, pTrack->Xtg_rec_tr, pTrack->Ytg_rec_tr, tmpZtg_tr_mm);
		Transform::P_HCS2TCS(P.Theta(), P.Phi(), pEndPlaneAngle, pTrack->Thetatg_rec_tr, pTrack->Phitg_rec_tr);


		////////////////////////////////////////////////////////
		//Reconstruction to the target plane is completed. Now
		//reconstruct theta, phi, x, y, z at vertex in HCS and TCS
		//
		//The snake always reconstruct to target plane
		//Now I have to calculate vertexZ 
		//then Drift from target plane to vertexZ plane 
		////////////////////////////////////////////////////////

		if(mWhereToStopRecon>=1)
		{
			//mWhereToStopRecon=0 is stop drifting at target plane, 1 is at vertex plane, 
			//2 is at exact vertex z0, 3 is exact z0 smeared with BPMXRes/sin(pEndPlaneAngle)

			//Get Vertex Z
			if(mWhereToStopRecon>=3)
			{
				//erconstruct to smeared thrown Z0
				//Reconstruct To BPMXRes/sin(pEndPlaneAngle) smeared Z0 
				const double pHRSZRes_mm = mBPMXRes/sin(pEndPlaneAngle); 
				pTrack->Z_rec=pTrack->Z0 + pHRSZRes_mm*HRSRand::fGaus();
			}
			else if(mWhereToStopRecon==2)
			{
				//Reconstruct To Exact Z0 
				pTrack->Z_rec = pTrack->Z0;
			}
			else
			{
				//reconstruct to V5tg_rec_tr predicted Z0
				//This line is from Nilanga's optics note, page 9
				//It has been debugged
				//double Y_BPM_tr_mm=pTrack->X0;
				double Y_BPM_tr_mm=pTrack->X0 + mBPMXRes*HRSRand::fGaus();
				pTrack->Z_rec = -pTrack->Ytg_rec_tr*cos(pTrack->Phitg_rec_tr) / sin(pEndPlaneAngle+pTrack->Phitg_rec_tr) + 
					Y_BPM_tr_mm/tan(pEndPlaneAngle+pTrack->Phitg_rec_tr);  //in unit of mm
				pTrack->Z_rec+=mTargetZOffset;
			}

			//////////////////////////////////////////
			//Drift from target plane to ideal vertexZ  
			//this part is complicate: if Z_rec is at downstream of target plane, we have to drift forward
			//otherwise drift backward

			if(pTrack->Z_rec/1000 > X.z())
			{
				//forward
				//set the tracklength limit to 0.2 m to save time
				//P is already ready for forward, make sure the particle is electron
				DriftSieve2Tg::Drift2Z(X,P,pMass,pCharge,pTrack->Z_rec/1000,0.2,0.0001);				
			}
			else if(pTrack->Z_rec/1000 < X.z())
			{
				//backward
				//set the tracklength limit to 0.2 m to save time
				P*=-1;	//flip P for backward, make sure the particle is positron
				DriftSieve2Tg::Drift2Z(X,P,pMass,-pCharge,pTrack->Z_rec/1000,0.2,0.0001);	
				P*=-1;  //flip it back
			}

		}

		//now convert the result from HCS to TCS, converted to mm at the same time 
		Transform::X_HCS2TCS(X.x()*1000-mPivotXOffset,X.y()*1000-mPivotYOffset, X.z()*1000-mPivotZOffset,
			pEndPlaneAngle, pTrack->X_rec_tr,pTrack->Y_rec_tr,pTrack->Z_rec_tr);
		Transform::P_HCS2TCS(P.Theta(), P.Phi(), pEndPlaneAngle,pTrack->Theta_rec_tr,pTrack->Phi_rec_tr);

		//////////////////////////////////////////
		//now store the vertex plane variables in mm

		pTrack->X_rec=X.x()*1000;
		pTrack->Y_rec=X.y()*1000;
		pTrack->Z_rec=X.z()*1000;  
		pTrack->Theta_rec=P.Theta();
		pTrack->Phi_rec=P.Phi();

	}


	////////////////////////////////////////////////////////
	//Energy loss correction
	////////////////////////////////////////////////////////	

	pTrack->P_rec = ElossCorr(pTrack->P_rec2tg);
	//Delta is defined by (Pvb-pHRSMomentum)/pHRSMomentum
	//therefore I do not want to correct Delta_rec after Eloss corr.
	//pTrack->Delta_rec = (pTrack->P_rec-pHRSMomentum)/pHRSMomentum;

	return bGoodParticle;
}


//////////////////////////////////////////////////////////////////////////////////////////////
//HRSRootTree::TransportThruHRS_NoTgField() has been rewritten. 
//Here is the 'plane' it goes thrught and the key word of the variables names:
//vertex  -> target-> sieve(vb) -> target     -> focal -> target     -> sieve       -> target    -> vertex
//$0,$0_tr-> $tg_tr-> $vb,$vb_tr->$_proj2tg_tr-> $fp_tr-> $_rec2tg_tr-> $_proj2sl_tr-> $tg_rec_tr-> $_rec,$_rec_tr
//	
//This routine will do the following:

//Section 1: 
//A) Convert vertex, vb variable from HTC to TCS, 
//B) Project vertex plane to target plane	
//C) Determine TrackClass, stop if it less than SkimLevel
//Section 2:
//A) Project Sieve plane (or VB plane) (_vb_tr) to target plane (_proj2tg_tr)
//B) Use snake model to transport the _proj2tg_tr to focal plane (_fp), and also
//   reconstruct them to target plane (_rec2tg_tr)
//C) Stop if fail to call SNAKE model
//Section 3:
//If the reconstrction is good, set the following tree variables:
//$fp_tr -> $_rec2tg_tr -> $_proj2sl_tr -> $tg_rec_tr -> $_rec,$_rec_tr
//A) Stored $fp_tr and $_rec2tg_tr
//B) If target field is off, set $tg_rec_tr=$rec2tg_tr, then calculate vertexZ, then
//   project from target plane to vertex plane
//C) If target field is on, do the following:
//   I)  Project $rec2tg_tr to sieve plane: $_proj2sl_tr 
//   II) Call Drift() or SNAKE_Sl2Tg() to get them back to target plane: $tg_rec_tr 
//   III)Calculate vertexZ, call Drift() move them from target plane to vertex plane: $_rec,$rec_tr 

bool HRSRootTree::TransportThruHRS_NoTgField(int i)
{
	MyTrack *pTrack=track[i];
	if(pTrack->Pvb/keV<1.0) return false;

	int    pIsLeftArm=(pTrack->Xvb>0)?1:0; //left arm or right arm
	double pEndPlaneAngle=(pIsLeftArm)?mLSeptumAngle:mRSeptumAngle;
	double pHRSMomentum=(pIsLeftArm)?mLHRSMomentum:mRHRSMomentum;
	double pTrackL_vb2tg=(pIsLeftArm)?mPivot2LHRSVBFace:mPivot2RHRSVBFace;

	//in this code, default mLHRSMomentum=mRHRSMomentum=0;
	//The purpose is to simulate the whole momentum range [0,4.0) Gev ignoring the delta effect of HRS
	//In this case, I have to set pHRSMomentum = Pvb event by event
	if(pHRSMomentum<1.0E-05) pHRSMomentum=pTrack->Pvb;

	// *************************************************************************
	//Section 1: 
	//A) Convert vertex, vb variable from HTC to TCS, 
	//B) Project vertex plane to target plane	
	//C) Determine TrackClass, stop if it less than SkimLevel
	// *************************************************************************


	//Convert HCS to TCS at vertex
	Transform::X_HCS2TCS(pTrack->X0-mPivotXOffset,pTrack->Y0-mPivotYOffset,
		pTrack->Z0-mPivotZOffset, pEndPlaneAngle, pTrack->X0_tr, pTrack->Y0_tr,pTrack->Z0_tr);
	Transform::P_HCS2TCS(pTrack->Theta0, pTrack->Phi0, pEndPlaneAngle, 
		pTrack->Theta0_tr, pTrack->Phi0_tr);

	////////////////////////////////////////////////////////////////////
	//Project|Drift this X0_tr,Y0_tr to target plane
	TVector3 X,P;  //ThreeVector in Hall system for Drift(), in unit of m, GeV

	double tmpZtg_tr=0;
	Transform::Project(pTrack->X0_tr, pTrack->Y0_tr,pTrack->Z0_tr, -pTrack->Z0_tr, pTrack->Theta0_tr, 
		pTrack->Phi0_tr,pTrack->Xtg_tr,pTrack->Ytg_tr,tmpZtg_tr);

	pTrack->Thetatg_tr=pTrack->Theta0_tr;
	pTrack->Phitg_tr=pTrack->Phi0_tr;

	////////////////////////////////////////////////////////////////////
	//Convert HCS to TCS at VB
	if(mUseSeptumField && mUseSeptumPlusStdHRS) 
	{
		//In this case I will use the septum field to propagate particles and place the Virtual Boundary at
		//the Q1 entrance, therefore the end plane angle is HRS angle and I should project the track back
		//to the hall center, not the pivot. To do that I should not apply the mPivotZOffset
		DriftSieve2Tg::SetUseSeptumInDrift(1);
		pEndPlaneAngle=(pIsLeftArm)?mLHRSAngle:mRHRSAngle;
		Transform::X_HCS2TCS(pTrack->Xvb-mPivotXOffset,pTrack->Yvb-mPivotYOffset,
			pTrack->Zvb, pEndPlaneAngle, pTrack->Xvb_tr, pTrack->Yvb_tr,pTrack->Zvb_tr);
	}
	else
	{
		DriftSieve2Tg::SetUseSeptumInDrift(0);
		Transform::X_HCS2TCS(pTrack->Xvb-mPivotXOffset,pTrack->Yvb-mPivotYOffset,
			pTrack->Zvb-mPivotZOffset, pEndPlaneAngle, pTrack->Xvb_tr, pTrack->Yvb_tr,pTrack->Zvb_tr);
	}

	Transform::P_HCS2TCS(pTrack->Thetavb, pTrack->Phivb, pEndPlaneAngle, 
		pTrack->Thetavb_tr, pTrack->Phivb_tr);

	////////////////////////////////////////////////////////////////////
	//Determine the TrackClass
	//choose eplice shape
	if(pow((pTrack->Phivb_tr-0.0)/0.035,2.0) + pow((pTrack->Thetavb_tr-0.0)/0.060,2.0) < 1.0) 
		pTrack->TrackClass=1;

	if(pTrack->TrackClass==1)
	{
		if( (pTrack->P0 - pTrack->Pvb) / pTrack->P0 < 0.05 ) pTrack->TrackClass=2;
	}

	if(pTrack->TrackClass==2)
	{
		if( pow((pTrack->Phivb_tr-0.0)/0.028,2.0) + pow((pTrack->Thetavb_tr-0.0)/0.050,2.0) < 1.0 &&
			(pTrack->P0 - pTrack->Pvb) / pTrack->P0 < 0.04 ) pTrack->TrackClass=3;
	}

	//Jixie: to make it run fast, stop tracking if less than skim level
	if (track[0]->TrackClass<iSkimLevel) return false;


	// *************************************************************************
	//Section 2:
	//  A) Project Sieve plane (or VB plane) (_vb_tr) to target plane (_proj2tg_tr)
	//  B) Use snake model to transport the _proj2tg_tr to focal plane (_fp), and also
	//     reconstruct them to target plane (_rec2tg_tr)
	//  C) Stop if fail to call SNAKE model
	// *************************************************************************

	pTrack->Delta = (pTrack->Pvb - pHRSMomentum) / pHRSMomentum;
	//We know that this particle will not go thrught HRS, stop here to save time
	if(fabs(pTrack->Delta)>0.1) return false;


	// To use the transportation functions (from target plane to focus plane), 
	// we have to project it back to the target plane (Z_tg=0), all length in unit of mm
	double Z_proj2tg_tr;
	Transform::Project(pTrack->Xvb_tr, pTrack->Yvb_tr, pTrack->Zvb_tr, -pTrack->Zvb_tr, pTrack->Thetavb_tr,
		pTrack->Phivb_tr, pTrack->X_proj2tg_tr, pTrack->Y_proj2tg_tr, Z_proj2tg_tr);


	//Fill the vector to be used in forward transportation functions 
	//in meter and rad
	double pV5_rec[5],pV5_fp[5];
	double pV5_proj2tg[5]={pTrack->X_proj2tg_tr/1000.,pTrack->Thetavb_tr,
		pTrack->Y_proj2tg_tr/1000.,pTrack->Phivb_tr,pTrack->Delta};

	bool bGoodParticle=false;
	//need to set the HRS arm and angle, 
	mHRSTranModel->SetArm((pIsLeftArm>0)?true:false); 
	mHRSTranModel->SetHRSAngle(pEndPlaneAngle); 

	//////////////////////////////////////////////////////////////////////
	//input are in m or rad, will convert angle into tan(angle) inside 
	//and also convert tan(angle) back to angle before return
	//bool SNAKE::SNAKEThruHRS(HRSTransport* pSnakeModel,double pX0_tg_m,
	//			  const double *V5_tg,double* pV5_fp, double* pV5_rec);
	//reconstructed values stored in pV5_rec

	bGoodParticle=mHRSTranModel->Forward(pV5_proj2tg, pV5_fp);
	if(!bGoodParticle) return false;

	bool bApplyVDCSmearing=true;
	if(bApplyVDCSmearing)
	{
		// Resolutions 
		double mWireChamberRes_x = 0.0013; //m;
		double mWireChamberRes_y = 0.0013; //m;
		double mWireChamberRes_theta = 0.0003; //rad;
		double mWireChamberRes_phi = 0.0003; //rad;
		pV5_fp[0] += mWireChamberRes_x * HRSRand::fGaus();
		pV5_fp[2] += mWireChamberRes_y * HRSRand::fGaus();
		pV5_fp[1] += mWireChamberRes_theta * HRSRand::fGaus();
		pV5_fp[3] += mWireChamberRes_phi * HRSRand::fGaus();
	}

	///////////////////////////////////////////////////////////////////////
	//Get the smeared BPM vetical position
	double X_BPM_tr_m=pTrack->Xtg_tr/1000. + mBPMYRes/1000.*HRSRand::fGaus();
	pV5_fp[4]=X_BPM_tr_m;
	mHRSTranModel->Backward(pV5_fp,pV5_rec);


	//Reconstruct use analyzer database
	double pV5_rec_db[5];
	if (mUseOpticsDB==1) 
	{
		if (pIsLeftArm) mRecUseDBL->CalcTargetCoords(pV5_fp, pV5_rec_db);
		else mRecUseDBR->CalcTargetCoords(pV5_fp, pV5_rec_db);
	}

	//////////////////////////////////////////////////////////////////////

	// *************************************************************************
	//Section 3:
	//  If the reconstrction is good, set the following tree variables:
	//  $fp_tr -> $_rec2tg_tr -> $_proj2sl_tr -> $tg_rec_tr -> $_rec,$_rec_tr
	//  A) Stored $fp_tr and $_rec2tg_tr
	//  B) If target field is off, set $tg_rec_tr=$rec2tg_tr, then calculate vertexZ, then
	//     project from target plane to vertex plane
	//  C) If target field is on, do the following:
	//     I)  Project $rec2tg_tr to sieve plane: $_proj2sl_tr 
	//     II) Call Drift() or SNAKE_Sl2Tg() to get them back to target plane: $tg_rec_tr 
	//  III)Calculate vertexZ, call Drift() move them from target plane to vertex plane: $_rec,$rec_tr 
	// *************************************************************************


	//Update TrackClass and store focal plane variables
	pTrack->TrackClass+=4;
	// in unit of mm and rad
	pTrack->Xfp_tr=pV5_fp[0]*1000.;
	pTrack->Thetafp_tr=pV5_fp[1];
	pTrack->Yfp_tr=pV5_fp[2]*1000.;
	pTrack->Phifp_tr=pV5_fp[3];

	if (mUseOpticsDB==1) {
		pTrack->Xtg_rec_db_tr = pV5_rec_db[0]*1000.;
		pTrack->Thetatg_rec_db_tr = pV5_rec_db[1];
		pTrack->Ytg_rec_db_tr = pV5_rec_db[2]*1000.;
		pTrack->Phitg_rec_db_tr = pV5_rec_db[3];
		pTrack->Delta_rec_db = pV5_rec_db[4];
	}

	pTrack->X_rec2tg_tr     = pV5_rec[0]*1000.;  //turn into mm
	pTrack->Theta_rec2tg_tr = pV5_rec[1];
	pTrack->Y_rec2tg_tr     = pV5_rec[2]*1000.;
	pTrack->Phi_rec2tg_tr   = pV5_rec[3];
	pTrack->Delta_rec       = pV5_rec[4];

	pTrack->P_rec2tg = pHRSMomentum*(1.0+pV5_rec[4]);	

	//By Jixie:
	//People will be very interested to compare the reconstruction results among stop
	//drifting at the following location: target plane, vertex plane and exact z0
	//therefore I introduce a variable mWhereToStopRecon to 
	//specify where to stop reconstruction
	//0 is at target plane, 1 is at vertex plane, 2 is at exact vertex z 

	//ReconstructToExactZ0 is to test how good the uncertainty will be if we know 
	//vertexZ perfectly (HRS always has trouble in vertex z resolutio, especially for small HRS angle).	

	//since no target field, I do not have to project the particle to sieve plane then Drift it back
	//But for completeness, I still want to set values to thse dumb variables
	//You can comment it out to speed up
	double Z_proj2sl_tr;		
	Transform::Project(pTrack->X_rec2tg_tr,pTrack->Y_rec2tg_tr,0.0,pTrackL_vb2tg,pTrack->Theta_rec2tg_tr,
		pTrack->Phi_rec2tg_tr,pTrack->X_proj2sl_tr,pTrack->Y_proj2sl_tr,Z_proj2sl_tr);

	pTrack->Xtg_rec_tr     = pTrack->X_rec2tg_tr;  
	pTrack->Thetatg_rec_tr = pTrack->Theta_rec2tg_tr; 
	pTrack->Ytg_rec_tr     = pTrack->Y_rec2tg_tr; 
	pTrack->Phitg_rec_tr   = pTrack->Phi_rec2tg_tr; 

	if(mWhereToStopRecon>=3)
	{
		//Reconstruct To BPMXRes/sin(pEndPlaneAngle) smeared Z0 
		const double pHRSZRes_mm = mBPMXRes/sin(pEndPlaneAngle); 
		pTrack->Z_rec = pTrack->Z0 + pHRSZRes_mm*HRSRand::fGaus();
	}
	else if(mWhereToStopRecon==2)
	{
		//ReconstructToExactZ0 
		pTrack->Z_rec = pTrack->Z0;
	}
	else if(mWhereToStopRecon==1)
	{
		//this line is from Nilanga's optics note, page 9
		double Y_BPM_tr_mm=pTrack->X0 + mBPMXRes*HRSRand::fGaus();
		pTrack->Z_rec = -pTrack->Ytg_rec_tr*cos(pTrack->Phitg_rec_tr) / sin(pEndPlaneAngle+pTrack->Phitg_rec_tr) + 
			Y_BPM_tr_mm/tan(pEndPlaneAngle+pTrack->Phitg_rec_tr);  //in unit of mm
		pTrack->Z_rec += mTargetZOffset;
	}
	else if(mWhereToStopRecon==0) pTrack->Z_rec = mTargetZOffset;

	//NO field, do linear projection to vertex plane
	double pVZ2TargetPlane=(pTrack->Z_rec-mTargetZOffset)*cos(pEndPlaneAngle);

	//The snake model has the target at the target center, Got confirmed by Min and John 
	double tmpZ_tg_snake=0;  //will be useful in case there is an offset
	Transform::Project(pTrack->Xtg_rec_tr,pTrack->Ytg_rec_tr,tmpZ_tg_snake,pVZ2TargetPlane-tmpZ_tg_snake,
		pTrack->Thetatg_rec_tr,pTrack->Phitg_rec_tr,pTrack->X_rec_tr,pTrack->Y_rec_tr,pTrack->Z_rec_tr);

	Transform::X_TCS2HCS(pTrack->X_rec_tr, pTrack->Y_rec_tr, pTrack->Z_rec_tr, pEndPlaneAngle, 
		pTrack->X_rec, pTrack->Y_rec, pTrack->Z_rec);
	pTrack->X_rec+=this->mTargetXOffset;
	pTrack->Y_rec+=this->mTargetYOffset;
	pTrack->Z_rec+=this->mTargetZOffset; 

	//This is my recostruction for Theta_rec and Phi_rec
	pTrack->Theta_rec_tr = pTrack->Thetatg_rec_tr; 
	pTrack->Phi_rec_tr = pTrack->Phitg_rec_tr; 
	Transform::P_TCS2HCS(pTrack->Theta_rec_tr, pTrack->Phi_rec_tr, pEndPlaneAngle, 
		pTrack->Theta_rec, pTrack->Phi_rec);



	////////////////////////////////////////////////////////
	//Energy loss correction
	////////////////////////////////////////////////////////	

	pTrack->P_rec = ElossCorr(pTrack->P_rec2tg);
	//Delta is defined by (Pvb-pHRSMomentum)/pHRSMomentum
	//therefore I do not want to correct Delta_rec after Eloss corr.
	//pTrack->Delta_rec = (pTrack->P_rec-pHRSMomentum)/pHRSMomentum;

	return bGoodParticle;
}


void HRSRootTree::DoRootTree()
{
	//process each primary tracks
	if(!file) Initilize();

	PartNum=mPrimaryGeneratorAction->GetParticleNum();
	Ei=mPrimaryGeneratorAction->GetIncidentEnergy()/GeV;
	Helicity=mPrimaryGeneratorAction->GetHelicity();

	//for (int j=0;j<MaxPrimaryNum;j++)   
	for (int j=0;j<PartNum;j++)   //use this line to speed up
	{	
		if (track[j]->Pvb/keV>1.0) track[j]->TrackClass=0;
		if(track[j]->TrackClass>=0)
		{
			bool bGoodParticle=false;

			//check which VB this particle hit
			if(track[j]->VBName.find("HRS") != string::npos)
			{ 
				// Transportation through HRS
				if(mUseHelmField) bGoodParticle=TransportThruHRS(j);
				else bGoodParticle=TransportThruHRS_NoTgField(j);
			}
			else if(track[j]->VBName.find("SBS") != string::npos)  
			{
				// Transportation through SuperBigBite
				bGoodParticle=TransportThruVD(j,mSuperBigBiteAngle);
			}
			else if(track[j]->VBName.find("HMS") != string::npos)  
			{
				// Transportation through HMS
				bGoodParticle=TransportThruHMS(j);
			}
			else if(track[j]->VBName.find("VD") != string::npos)  
			{
				// Transportation through a Virtual Plane
				bGoodParticle=TransportThruVD(j,mVDAngle);
			}
			else if(track[j]->VBName.find("LAC") != string::npos)  
			{
				// Transportation through a Virtual Plane
				bGoodParticle=TransportThruVD(j,mLACAngle);
			}
			else if(track[j]->VBName.find("BB") != string::npos)  
			{
				// Transportation through BigBite
				bGoodParticle=TransportThruBigBite(j);
			}
			else if(track[j]->VBName.find("RTPC") != string::npos)  
			{
				// Transportation through RTPC, 
				// Digitization: for each hit in the drift region, need not drift ionized electron from 
				// ionization location to the PCB pad to get PadID (->ID), DriftTime (->TDC)
				// convert the deposit energy of this hit to ADC 
				// Then using aother reconstruct class to do the reconstruction using chain of (ID, TDC, ADC)  
				// the above feature is not built yet
				// bGoodParticle=TransportThruRTPC(j);
			}
			else
			{
				//TBD, do nothing now
				//double pEndPlaneAngle=mPrimaryGeneratorAction->GetDetectorAngle(i);
				//bGoodParticle=TransportThruVD(j,mVDAngle);
			}

			//in order to speed up this program
			//do not get this value if this is a bad particle
			if(bGoodParticle) 
			{
				track[j]->Theta0Eff = mPrimaryGeneratorAction->GetThetaEff(j); 
				track[j]->Pol = mPrimaryGeneratorAction->GetPolarization(j);
			}

			//calculate elas and inelas XS, to save time, only do it for particle 0
			//if this track is not electron, need to find out the theta of the coupled electron
			if(j==0 || j<3)
			{
				if(mCalculateXS==1 || mCalculateXS==3)
				{

					double pTheta=track[j]->Theta0;
					if(track[j]->PdgId!=11 && track[j]->PdgId!=22) 
					{
						//Ei - P_p * cos(Theta_p) = P_e * cos(Theta_e)   (1)
						//P_p * sin(Theta_p) = P_e * sin(Theta_e)        (2) 
						// (2)/(1) ==> tan(Theta_e) = (P_p * sin(Theta_p)) / (Ei - P_p * cos(Theta_p))  
						pTheta = atan2(track[j]->P0 * sin(track[j]->Theta0), Ei - track[j]->P0 * cos(track[j]->Theta0));
					}
					//for elastic cross section, all in GeV, return in microbarn
					//calculate the elas XS based on the given target mass in cmd line
					
					if(track[j]->PdgId == 22 || mPrimaryGeneratorAction->GetPrimaryEngine(j)=="Compton")
						track[j]->ElasXS=Compton::GetRCSXS_lab(Ei,pTheta)*0.001;  //turn nb into ub
					else
						track[j]->ElasXS=ElasModel::GetXS(mTargetAtomicNumber,
							mTargetNeutronNumber,mBeamEnergy,pTheta,mTargetMass);
					if(track[j]->ElasXS<0.0) track[j]->ElasXS=0.0;
				}

				//for inelastic cross section, all in GeV, return in microbarn/GeV
				//sometimes QElXS are too slow, I do not want to do this calculate by default
				if(mCalculateXS==2 || mCalculateXS==3)
				{
					int pid=track[j]->PdgId;
					if( pid==2212 || pid==2112 || abs(pid)==211 || pid==111 || pid==11 || pid==22 
						/*|| abs(pid)==321*/)
					{
						//calculate the XS based on the given target in cmd line
						//if(mTargetNeutronNumber<0.1 && mTargetAtomicNumber>1.0) 
						//{
						//const double a2mev=931.49404 * MeV;   //atomic mass unit, transfer A to MeV
						//mTargetNeutronNumber= ceil(mTargetMass*GeV/a2mev-mTargetAtomicNumber-0.4);
						//}
						//track[j]->XS=GetInelasXS(track[j]->PdgId,(int)mTargetAtomicNumber,(int)mTargetNeutronNumber,
						//mBeamEnergy,track[j]->Theta0,track[j]->P0);
						//

						//calculate the inelas XS based on the material at vertex
						//double getQElXS(int PID,double Eb, double Theta, double Pf, G4Material* theMaterial)
						track[j]->XS=GetInelasXS(track[j]->PdgId,mBeamEnergy,track[j]->Theta0,track[j]->P0,material[j]);

					}
					if(track[j]->XS<0.0) track[j]->XS=0.0;
				}
			}

		}

		//do not process the other track if this track is garbage
		if (track[0]->TrackClass<iSkimLevel || track[j]->P0<0) break; 
	}

	//do not fill the tree if the first track is garbage 
	if (track[0]->TrackClass>=iSkimLevel && track[0]->P0>0) 
	{
		FillTree();
	} 

	//fill the detector response tree
	if(!iNoDetectorResponse )
	{
		detector->Fill();
	}

	iEvtCounter+=1;
	if (iEvtCounter==1)	printf("\nHRSRootTree():: run %d started......\n",iRunNumber);
	else if (!(iEvtCounter%1000) && iEvtCounter)
	{
		printf("Root#=%6d Total#=%6d events processed......\r",iRootEvtID,iEvtCounter);
	}

	Reset(); 
	return;
}

