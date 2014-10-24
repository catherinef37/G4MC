
#ifndef HRSNtReader_h
#define HRSNtReader_h

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TCanvas.h"
#include "TObject.h"
#include "TMath.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include "math.h"
#include <string>
#include <iostream>
using namespace std;

class HRSNtReader {

public :
	HRSNtReader(const char *filename, const char *treename, int loadconfigtree=0);
	virtual ~HRSNtReader();

	virtual Bool_t   ReadConfig(const char *filename);
	virtual void     Init();
	virtual Int_t    Cut();
	virtual Int_t    GetEntry(Long64_t entry);
	virtual void     Loop();

public :
	TChain          *mChain;   //!pointer to the analyzed TTree or TChain
	Long64_t        mEntries;  //number of events in this chain
	bool            mIsCombinedTree;


	//Declaration of leaves types for config tree
	Int_t           Run;
	Int_t           SkimLevel;
	Int_t           BookTrees;
	Double_t        Beam;
	Double_t        TargetM;
	Double_t        TargetAtomicNumber;
	Double_t        LHRSAngle;
	Double_t        RHRSAngle;
	Double_t        TargetXOffset;
	Double_t        TargetYOffset;
	Double_t        TargetZOffset;
	Double_t        PivotXOffset;
	Double_t        PivotYOffset;
	Double_t        PivotZOffset;
	Double_t        LHRSMomentum;
	Double_t        RHRSMomentum;
	Int_t           UseHelmField;
	Double_t        HelmXOffset;
	Double_t        HelmYOffset;
	Double_t        HelmZOffset;
	Double_t        HelmRotAxis1;
	Double_t        HelmRotAxis2;
	Double_t        HelmRotAxis3;
	Double_t        HelmRotAngle1;
	Double_t        HelmRotAngle2;
	Double_t        HelmRotAngle3;
	Double_t        HelmCurrentRatio;
	Int_t           UseSeptumField;
	Double_t        SeptumXOffset;
	Double_t        SeptumYOffset;
	Double_t        SeptumZOffset;
	Double_t        SeptumRotAxis1;
	Double_t        SeptumRotAxis2;
	Double_t        SeptumRotAxis3;
	Double_t        SeptumRotAngle1;
	Double_t        SeptumRotAngle2;
	Double_t        SeptumRotAngle3;
	Double_t        SeptumCurrentRatioL;
	Double_t        SeptumCurrentRatioR;
	Double_t        BigBiteAngle;
	Double_t        BigBiteTiltAngle;
	Double_t        Pivot2BigBiteFace;

	
	// Declaration of leaf types
	Int_t           Index;
	Int_t           PdgId;
	Int_t           TrackId;
	Int_t           TrackClass;
	Double_t        X0;
	Double_t        Y0;
	Double_t        Z0;
	Double_t        P0;
	Double_t        Theta0;
	Double_t        Phi0;
	Double_t        X0_tr;
	Double_t        Y0_tr;
	Double_t        Z0_tr;
	Double_t        Theta0_tr;
	Double_t        Phi0_tr;
	Double_t        Xvb;
	Double_t        Yvb;
	Double_t        Zvb;
	Double_t        Pvb;
	Double_t        Thetavb;
	Double_t        Phivb;
	Double_t        Xvb_tr;
	Double_t        Yvb_tr;
	Double_t        Zvb_tr;
	Double_t        Thetavb_tr;
	Double_t        Phivb_tr;
	Double_t        Xfp_tr;
	Double_t        Yfp_tr;
	Double_t        Thetafp_tr;
	Double_t        Phifp_tr;
	Double_t        Delta;
	Double_t        Delta_rec;
	Double_t        X_rec_tr;
	Double_t        Y_rec_tr;
	Double_t        Theta_rec_tr;
	Double_t        Phi_rec_tr;
	Double_t        Z_rec;
	Double_t        P_rec;
	Double_t        Theta_rec;
	Double_t        Phi_rec;
	Double_t        TrackRadlen;	
	Double_t        Theta0Eff;	
	Double_t        ElasXS;
	Double_t        XS;

	Int_t           StepNum;
	Double_t        StepX[1024];   //[StepNum]
	Double_t        StepY[1024];   //[StepNum]
	Double_t        StepZ[1024];   //[StepNum]
	Double_t        StepdE[1024];   //[StepNum]
	Double_t        StepL[1024];   //[StepNum]
	Double_t        StepEkin[1024];   //[StepNum]
	Double_t        StepTL[1024];   //[StepNum]
	Double_t        StepRadlen[1024];   //[StepNum]
	Double_t        StepDsty[1024];   //[StepNum]
	Double_t        StepBx[1024];   //[StepNum]
	Double_t        StepBy[1024];   //[StepNum]
	Double_t        StepBz[1024];   //[StepNum]
	Double_t        TrackBdLx;
	Double_t        TrackBdLy;
	Double_t        TrackBdLz;
	Double_t        R0;
	Double_t        A0;
	Double_t        B0;

};

#endif

