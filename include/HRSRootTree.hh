// ********************************************************************
//
// $Id: HRSRootTree.hh,v 1.0, 2010/12/26  12:14:12 HRS Exp $
// GEANT4 tag $Name: geant4-09-04 $ 
//
//..............................................................................

#ifndef HRSRootTree_h
#define HRSRootTree_h 1

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include "G4RotationMatrix.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"
#include "HRSRecUseDB.hh"

class HRSTransport;
class HRSPrimaryGeneratorAction;

//maximum number of steps in a track
#define MaxStepPerTrack  1024

//maximum number of hits in a SD
#define MaxSDHit 1024

//maximum number of hits in a track
#define MaxTrackHit 128
//maximum number of tracks in the whole event
#define MaxParticle 10240


//#define ROTATION_DEBUG 1

//this number must equal to the one in HRSPrimaryGeneratorAction.hh
#ifndef MaxPrimaryNum
#define MaxPrimaryNum 8
#endif

//The normal unit in the root tree is GeV, mm, rad
class MyTrack
{
public:
  int	   PdgId;	//pdg code for this particle
  int	   TrackId;
  int    TrackClass;
  /*
    #    Flag to identify a good track. For all HRS track. the conditions are:
    #    (A) hit virtual boundary;
    #    (B) pow(Thetavb_tr-0.0)/0.060,2.0)+pow(Phivb_tr-0.0)/0.035,2.0)<1.0 
    #    (C) dPOverP0 Cut: (P0-Pvb)/P0<0.05
    #    (D) pow(Thetavb_tr-0.0)/0.050,2.0)+pow(Phivb_tr-0.0)/0.028,2.0)<1.0 && (P0-Pvb)/P0<0.04
    #    -1 garbage: can not hit the virtual boundary
    #    0 barely normal: A satisfied 
    #    1 normal: A, B satisfied
    #    2 good: A, B, C satisfied
    #    3 gold: A, B, C, D satisfied
    #    will plus 4 if it can be transported to the focus plane 
  */
  
  double  XS_208Pb;
  
  double  X0;
  double  Y0;
  double  Z0;
  double  P0;
  double  Theta0;
  double  Phi0;
  double  X0_tr;
  double  Y0_tr;
  double  Z0_tr;
  double  Theta0_tr;
  double  Phi0_tr;
  double  P0_tr;
  
  //the image of X0_tr,Y0_tr on target plane
  double  Xtg_tr;		
  double  Ytg_tr;
  double  Thetatg_tr;		
  double  Phitg_tr;
  double  Ptg_tr;
  
  double  X_vb;
  double  Y_vb;
  double  Z_vb;
  double  P_vb;
  double  Theta_vb;
  double  Phi_vb;
  double  X_vb_tr;
  double  Y_vb_tr;
  double  Z_vb_tr;
  double  Theta_vb_tr;
  double  Phi_vb_tr;
  double  P_vb_tr;
  
  double  X_vdc;
  double  Y_vdc;
  double  Z_vdc;
  double  Theta_vdc;
  double  Phi_vdc;
  double  P_vdc;
  double  X_vdc_tr;
  double  Y_vdc_tr;
  double  Z_vdc_tr;
  double  Theta_vdc_tr;
  double  Phi_vdc_tr;
  double  P_vdc_tr;

  double  X_qz1;
  double  Y_qz1;
  double  Z_qz1;
  double  Theta_qz1;
  double  Phi_qz1;
  double  P_qz1;
  double  X_qz1_tr;
  double  Y_qz1_tr;
  double  Z_qz1_tr;
  double  Theta_qz1_tr;
  double  Phi_qz1_tr;
  double  P_qz1_tr;

  double  X_qz2;
  double  Y_qz2;
  double  Z_qz2;
  double  Theta_qz2;
  double  Phi_qz2;
  double  P_qz2;
  double  X_qz2_tr;
  double  Y_qz2_tr;
  double  Z_qz2_tr;
  double  Theta_qz2_tr;
  double  Phi_qz2_tr;
  double  P_qz2_tr;

  double  X_fp;
  double  Y_fp;
  double  Z_fp;
  double  Theta_fp;
  double  Phi_fp;
  double  P_fp;
  double  X_fp_tr;
  double  Y_fp_tr;
  double  Z_fp_tr;
  double  Theta_fp_tr;
  double  Phi_fp_tr;
  double  P_fp_tr;
  
  double  X_fpr_tr;
  double  Y_fpr_tr;
  double  Z_fpr_tr;
  double  Theta_fpr_tr;
  double  Phi_fpr_tr;
  double  P_fpr_tr;
  
  double  X_pr_tr;
  double  Y_pr_tr;
  double  Z_pr_tr;
  double  Theta_pr_tr;
  double  Phi_pr_tr;
  double  P_pr_tr;
  
  double  X_q1en_tr;
  double  Y_q1en_tr;
  double  Z_q1en_tr;
  double  Theta_q1en_tr;
  double  Phi_q1en_tr;
  double  P_q1en_tr;
  int     a_q1en_tr;

  double  X_q1ex_tr;
  double  Y_q1ex_tr;
  double  Z_q1ex_tr;
  double  Theta_q1ex_tr;
  double  Phi_q1ex_tr;
  double  P_q1ex_tr;
  int     a_q1ex_tr;

  double  X_q2en_tr;
  double  Y_q2en_tr;
  double  Z_q2en_tr;
  double  Theta_q2en_tr;
  double  Phi_q2en_tr;
  double  P_q2en_tr;
  int     a_q2en_tr;

  double  X_q2ex_tr;
  double  Y_q2ex_tr;
  double  Z_q2ex_tr;
  double  Theta_q2ex_tr;
  double  Phi_q2ex_tr;
  double  P_q2ex_tr;
  int     a_q2ex_tr;

  double  X_den_tr;
  double  Y_den_tr;
  double  Z_den_tr;
  double  Theta_den_tr;
  double  Phi_den_tr;
  double  P_den_tr;
  int     a_den_tr;

  double  X_dex_tr;
  double  Y_dex_tr;
  double  Z_dex_tr;
  double  Theta_dex_tr;
  double  Phi_dex_tr;
  double  P_dex_tr;
  int     a_dex_tr;

  double  X_q3en_tr;
  double  Y_q3en_tr;
  double  Z_q3en_tr;
  double  Theta_q3en_tr;
  double  Phi_q3en_tr;
  double  P_q3en_tr;
  int     a_q3en_tr;

  double  X_q3ex_tr;
  double  Y_q3ex_tr;
  double  Z_q3ex_tr;
  double  Theta_q3ex_tr;
  double  Phi_q3ex_tr;
  double  P_q3ex_tr;
  int     a_q3ex_tr;

  double  X_sen_tr;
  double  Y_sen_tr;
  double  Z_sen_tr;
  double  Theta_sen_tr;
  double  Phi_sen_tr;
  double  P_sen_tr;
  int     a_sen_tr;

  double  X_sm_tr;
  double  Y_sm_tr;
  double  Z_sm_tr;
  double  Theta_sm_tr;
  double  Phi_sm_tr;
  double  P_sm_tr;
  int     a_sm_tr;

  double  X_sex_tr;
  double  Y_sex_tr;
  double  Z_sex_tr;
  double  Theta_sex_tr;
  double  Phi_sex_tr;
  double  P_sex_tr;
  int     a_sex_tr;

  double  X_coil_tr;
  double  Y_coil_tr;
  double  Z_coil_tr;
  double  Theta_coil_tr;
  double  Phi_coil_tr;
  double  P_coil_tr;
  int     a_coil_tr;

  double  X_mid_tr;
  double  Y_mid_tr;
  double  Z_mid_tr;
  double  Theta_mid_tr;
  double  Phi_mid_tr;
  double  P_mid_tr;
  int     a_mid_tr;

  double  X_col_tr;
  double  Y_col_tr;
  double  Z_col_tr;
  double  Theta_col_tr;
  double  Phi_col_tr;
  double  P_col_tr;
  int     a_col_tr;

  double  X_q1en;
  double  Y_q1en;
  double  Z_q1en;
  double  Theta_q1en;
  double  Phi_q1en;
  double  X_q1ex;
  double  Y_q1ex;
  double  Z_q1ex;
  double  Theta_q1ex;
  double  Phi_q1ex;
  double  X_q2en;
  double  Y_q2en;
  double  Z_q2en;
  double  Theta_q2en;
  double  Phi_q2en;
  double  X_q2ex;
  double  Y_q2ex;
  double  Z_q2ex;
  double  Theta_q2ex;
  double  Phi_q2ex;
  double  X_den;
  double  Y_den;
  double  Z_den;
  double  Theta_den;
  double  Phi_den;
  double  X_dex;
  double  Y_dex;
  double  Z_dex;
  double  Theta_dex;
  double  Phi_dex;
  double  X_q3en;
  double  Y_q3en;
  double  Z_q3en;
  double  Theta_q3en;
  double  Phi_q3en;
  double  X_q3ex;
  double  Y_q3ex;
  double  Z_q3ex;
  double  Theta_q3ex;
  double  Phi_q3ex;
  double  X_sen;
  double  Y_sen;
  double  Z_sen;
  double  Theta_sen;
  double  Phi_sen;
  double  X_sm;
  double  Y_sm;
  double  Z_sm;
  double  Theta_sm;
  double  Phi_sm;
  double  X_sex;
  double  Y_sex;
  double  Z_sex;
  double  Theta_sex;
  double  Phi_sex;
  double  X_coil;
  double  Y_coil;
  double  Z_coil;
  double  Theta_coil;
  double  Phi_coil;
  double  X_mid;
  double  Y_mid;
  double  Z_mid;
  double  Theta_mid;
  double  Phi_mid;
  double  X_col;
  double  Y_col;
  double  Z_col;
  double  Theta_col;
  double  Phi_col;
  
  
  //reconstructed variables at target plane
  double  Xtg_rec_tr;
  double  Ytg_rec_tr;
  double	Thetatg_rec_tr;
  double  Phitg_rec_tr;
  
  double  Xtg_rec_db_tr;
  double  Ytg_rec_db_tr;
  double	Thetatg_rec_db_tr;
  double  Phitg_rec_db_tr;
  
  double  X_rec_tr;
  double  Y_rec_tr;
  double  Z_rec_tr;
  double	Theta_rec_tr;
  double  Phi_rec_tr;
  double  X_rec;
  double  Y_rec;
  double  Z_rec;
  double  P_rec;
  double  Theta_rec;
  double  Phi_rec;
  
	double  Delta;
	double  Delta_rec;
	double  Delta_rec_db;

	double  TrackRadlen; //this is the sum of StepLen/RadLen for the whole track
	double  Theta0Eff;
	double  ElasXS;
	double  XS;

	double  R0;
	double  A0;
	double  B0;

	///////////////////////////////////////
	//project from sieve plane to target plane for calling snake forward 
	double  X_proj2tg_tr;
	double  Y_proj2tg_tr;

	//to store snake backward result
	double  X_rec2tg_tr;
	double  Y_rec2tg_tr;
	double	Theta_rec2tg_tr;
	double  Phi_rec2tg_tr;
	double  P_rec2tg;

	//project snake result to sieve plane
	double  X_proj2sl_tr;
	double  Y_proj2sl_tr;
	
	double  Pol;
	///////////////////////////////////////

	int     StepNum;
	double  StepX[MaxStepPerTrack];
	double  StepY[MaxStepPerTrack];
	double  StepZ[MaxStepPerTrack];
	double  StepEkin[MaxStepPerTrack];
	double  StepdE[MaxStepPerTrack];
	double  StepL[MaxStepPerTrack];

	//non-mini tree
	double  StepTL[MaxStepPerTrack];
	double  StepRadlen[MaxStepPerTrack];
	double  StepDsty[MaxStepPerTrack];

	//added by Jixie @20110525	
	double  StepBx[MaxStepPerTrack];
	double  StepBy[MaxStepPerTrack];
	double  StepBz[MaxStepPerTrack];
	//TrackBdLx = #inte{v[1]*b[2]-b[1]*v[2]}dt=L[1]*b[2]-b[1]*L[2];
	//TrackBdLy = #inte{v[2]*b[0]-b[2]*v[0]}dt=L[2]*b[0]-b[2]*L[0];
	//TrackBdLz = #inte{v[0]*b[1]-b[0]*v[1]}dt=L[0]*b[1]-b[0]*L[1];
	double  TrackBdLx; //this is the sum of B*dL in x for the whole track
	double  TrackBdLy; //this is the sum of B*dL in y for the whole track
	double  TrackBdLz; //this is the sum of B*dL in z for the whole track	

	string VBName;
};



class HRSRootTree
{
public:
	int     TA1_N;	//must be initialized before using
	int     TA1_Pid[MaxSDHit];
	int     TA1_Tid[MaxSDHit];
	int     TA1_ParentTid[MaxSDHit];
	double  TA1_T[MaxSDHit];				//in ns
	double  TA1_X[MaxSDHit];				//in mm
	double  TA1_Y[MaxSDHit];
	double  TA1_Z[MaxSDHit];
	double  TA1_Edep[MaxSDHit];			//in MeV
	double  TA1_NonIonEdep[MaxSDHit];	//in MeV
	double  TA1_P[MaxSDHit];
	double  TA1_Theta[MaxSDHit];
	double  TA1_Phi[MaxSDHit];
	double  TA1_Pout[MaxSDHit];

	int     TA2_N;	//must be initialized before using
	int     TA2_Pid[MaxSDHit];
	int     TA2_Tid[MaxSDHit];
	int     TA2_ParentTid[MaxSDHit];
	double  TA2_T[MaxSDHit];
	double  TA2_X[MaxSDHit];
	double  TA2_Y[MaxSDHit];
	double  TA2_Z[MaxSDHit];
	double  TA2_Edep[MaxSDHit];
	double  TA2_NonIonEdep[MaxSDHit];	//in MeV
	double  TA2_P[MaxSDHit];
	double  TA2_Theta[MaxSDHit];
	double  TA2_Phi[MaxSDHit];
	double  TA2_Pout[MaxSDHit];

	//for general sensitive detectors
	int     SD_N;	//must be initialized before using
	int     SD_Id[MaxSDHit];
	int     SD_Pid[MaxSDHit];
	int     SD_Tid[MaxSDHit];
	int     SD_ParentTid[MaxSDHit];
	double  SD_T[MaxSDHit];
	double  SD_X[MaxSDHit];
	double  SD_Y[MaxSDHit];
	double  SD_Z[MaxSDHit];
	double  SD_Edep[MaxSDHit];
	double  SD_NonIonEdep[MaxSDHit];	//in MeV
	double  SD_P[MaxSDHit];
	double  SD_Theta[MaxSDHit];
	double  SD_Phi[MaxSDHit];
	double  SD_Pout[MaxSDHit];

	//store the trajectory of every track
	int     T_N;	//must be initialized before using
	int     T_Pid[MaxParticle];
	int     T_Tid[MaxParticle];
	int     T_ParentTid[MaxParticle];
	double  T_T[MaxParticle];
	double  T_X[MaxParticle];
	double  T_Y[MaxParticle];
	double  T_Z[MaxParticle];
	double  T_P[MaxParticle];
	double  T_Theta[MaxParticle];
	double  T_Phi[MaxParticle];

	int     T_StepN[MaxParticle];
	double  T_StepX[MaxParticle][MaxTrackHit];
	double  T_StepY[MaxParticle][MaxTrackHit];
	double  T_StepZ[MaxParticle][MaxTrackHit];

public:
	HRSRootTree(int pRunNumber=1);
	~HRSRootTree();

	void Initilize();

	void FillTree(); // fill tree
	void Reset();
	void DoRootTree();

	void SetRunNumber(int val){iRunNumber=val;};
	int  GetRunNumber(){return iRunNumber;};

	int  GetRootEvtID(){return iRootEvtID;};
	int  GetTotalEvtID(){return iEvtCounter;};

	char* GetRootFileName(){return rootfile;};

	//transport particles through HRS, calling fortan routines
	bool TransportThruHRS(int i);
	bool TransportThruHMS(int i);
	bool TransportThruBigBite(int i);
	bool TransportThruHRS_NoTgField(int i);
	bool TransportThruVD(int i, double EndPlaneAngle);

	//Get the effective BPM vertical position when target field is on
	//input: true BPM vertical position in meter and particle momentum in GeV
	//return: effective BPM vertical position in meter
	double GetEffBPMVertical(double pX_BPM_tr_m,double pMomentum_GeV);

	//do the energy loss correction for the electron
	double ElossCorr(double pMomentum_GeV);


	void SetSDNum(int val){mSDNum=val;};	
	void SetSDName(int i, const char* val){if(i<MaxSDHit) sprintf(mSDName[i],"%s",val);};

	bool GetConfigTreeFilledFlag(){return bConfigTreeFilled;};
  float InterpolateNickie(float Q2_in);

public:
	// tree variables
	int    iRunNumber;
	MyTrack* track[MaxPrimaryNum];
	int    iNoDetectorResponse;
	G4Material* material[MaxPrimaryNum]; //material @ vertex

private:
	TTree *tree[MaxPrimaryNum]; // hits info, event
	TTree *config; // run configuration information
	TTree *detector; //detector tree
	TFile *file;
	char rootfile[200];

	/////////////////////////////////////////////////////
	HRSRecUseDB *mRecUseDBL,*mRecUseDBR;

	//from argument
	int    mUseOpticsDB;

	double mBeamEnergy,mLHRSMomentum,mRHRSMomentum,mBeamTiltedAngle;
	double mTargetMass,mTargetL,mTargetAtomicNumber,mTargetNeutronNumber;
	int    mUseSeptumPlusStdHRS;

	/////////////////////////////////////////////////////
	//global config variable, located in detector.ini
	double mTargetXOffset,mTargetYOffset,mTargetZOffset;
	double mPivotXOffset,mPivotYOffset, mPivotZOffset;


	int    mSetupLHRS,mSetupRHRS;
	double mLSeptumAngle,mRSeptumAngle;
	double mLHRSAngle,mRHRSAngle;

	double mPivot2LHRSVBFace,mPivot2RHRSVBFace;

	int    mSetupG2PGeometry;
	int    mSetupCREXGeometry;

	/////////////////////////////////////////////////////
	int    mSetupG2PTarget,mTargetType;
	double mThirdArmAngle,mSetupThirdArmVD,mThirdArmRotZAngle,mPivot2ThirdArmFace;

	int    mUseHelmField;
	double mHelmXOffset, mHelmYOffset, mHelmZOffset;
	double mHelmRotAxis1,mHelmRotAxis2,mHelmRotAxis3;
	double mHelmRotAngle1,mHelmRotAngle2,mHelmRotAngle3;
	double mHelmCurrentRatio;

	int    mUseSeptumField;
	double mSeptumXOffset, mSeptumYOffset, mSeptumZOffset;
	double mSeptumRotAxis1,mSeptumRotAxis2,mSeptumRotAxis3;
	double mSeptumRotAngle1,mSeptumRotAngle2,mSeptumRotAngle3;
	double mSeptumCurrentRatioL,mSeptumCurrentRatioR;


	int    mSetupChicane,mSetupChicaneVD;
	double mFZB1TiltedAngle,mFZB1PosX,mFZB1PosY,mFZB1PosZ;
	double mFZB1Bx,mFZB1By,mFZB1Bz;
	double mFZB2TiltedAngle,mFZB2PosX,mFZB2PosY,mFZB2PosZ;
	double mFZB2Bx,mFZB2By,mFZB2Bz;

	/////////////////////////////////////////////////////
	int    mSetupVD;
	double mVDAngle,mVDTiltAngle,mPivot2VDFace;

	/////////////////////////////////////////////////////
	int    mSetupLAC;
	double mLACAngle,mLACTiltAngle,mPivot2LACFace;

	/////////////////////////////////////////////////////
	int    mSetupBigBite;
	double mBigBiteAngle,mBigBiteTiltAngle,mPivot2BigBiteFace;

	int    mSetupHAND;
	double mPivot2HANDLeadWall;
	/////////////////////////////////////////////////////
	int    mSetupRTPC;
	double mRatioHe2DME;

	/////////////////////////////////////////////////////
	int    mSetupSuperBigBite;
	double mSuperBigBiteAngle,mPivot2SuperBigBiteFace;

	int    mUseSBSField;
	double mSBSFieldXOffset, mSBSFieldYOffset, mSBSFieldZOffset;
	double mSBSFieldRotAxis1,mSBSFieldRotAxis2,mSBSFieldRotAxis3;
	double mSBSFieldRotAngle1,mSBSFieldRotAngle2,mSBSFieldRotAngle3;
	double mSBSFieldCurrentRatio;

	int    mSetupHMS;
	double mHMSAngle,mPivot2HMSFace,mHMSMomentum;

	/////////////////////////////////////////////////////

	int  mSDNum;
	char mSDName[MaxSDHit][100];
	//some variable for reconstruction 
	//specify where to stop reconstruction, 0 is at target plane, 1 is at vertex plane,
	//2 is at exact vertex z0, 3 is ezact z0 smeared with BPMXRes/sin(pEndPlaneAngle) 
	int    mWhereToStopRecon;	

	//Use SnakeModel to identify which SNAKE packages
	//will be use, SnakeModel==10-19 is for g2p, 20 for E97110 GDH experiment
	//Here is a list of SnakeModel candidates:	
	//SnakeModel=1: Standard HRS setting, no septum field, no target field
	//SnakeModel=11: g2p septum 484816+shim, 5.65 deg, 3cm raster, Created by John, 
	//SnakeModel=12: g2p septum 403216+shim, 5.65 deg, SNAKE Model not ready yet 
	//SnakeModel=13: g2p septum 400016+shim, 5.65 deg, SNAKE Model not ready yet 
	//SnakeModel=19: g2p septum 484816, no shim, 5.65 deg, Min's version for normal septum without shim
	//SnakeModel=20: GDH exp with normal X0 version
	//SnakeModel=21: GDH exp with large X0 version
	int    mSnakeModel;
	HRSTransport* mHRSTranModel;

	double mBPMXRes,mBPMYRes;

	//---------------------------------------
	bool     bConfigTreeFilled;

	int      iSkimLevel,iBookTrees,iBookHistos;
	int      iEvtCounter;
	int      iRootEvtID;     //event index
	//---------------------------------------

	//root variables, do not follow name rules 
	int	 PartNum;	//primary particle number
	//This is the incident energy, for electron beam it is the beam nergy
	//for photon beam it is the radiated beam energy
	double  Ei;   
	int     Helicity;

	int  mGenHistoOnly;
	int  mCalculateXS;   //0 none, 1 Elas, 2 Inelas XS, 3 both

	//histogram
	TH1D *hP0[MaxPrimaryNum], *hTheta0[MaxPrimaryNum], *hPhi0[MaxPrimaryNum];
	TH1D *hX0[MaxPrimaryNum], *hY0[MaxPrimaryNum], *hZ0[MaxPrimaryNum];
	TH2D *h2Z0VSY0[MaxPrimaryNum], *h2Y0VSX0[MaxPrimaryNum]; 
	TH2D *h2Theta0_trVSPhi0_tr[MaxPrimaryNum], *h2Theta_vb_trVSPhi_vb_tr[MaxPrimaryNum];
	TH3D *h3Z0Y0X0[MaxPrimaryNum];
	TH2D *h2P0VSTheta0[MaxPrimaryNum], *h2Theta0VSPhi0[MaxPrimaryNum];
	//TH3D *h3P0Theta0Phi0[MaxPrimaryNum];

	HRSPrimaryGeneratorAction* mPrimaryGeneratorAction;

};

#endif


