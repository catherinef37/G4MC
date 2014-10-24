//this script is used to plot the acceptance for one perticular 
//sensertive detector, one has to provide the sdid, which can be
//found in the config tree. simply type config->Show(0) and then all 
//sdid will be shown
#include "stdlib.h" 
#include <iostream>
#include "math.h"
using namespace std;

#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TQObject.h>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TString.h"
#include "TCut.h"
#include "TCutG.h"
#include "TPaveText.h"
#include "TText.h"
#include "TPad.h"
#include "TF1.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"


//Declaration of leaves types
Int_t           Run=0;
Int_t           SkimLevel;
Int_t           BookTrees;
Double_t        Beam;
Double_t        BeamTiltedAngle;
Double_t        TargetM;
Double_t        TargetL;
Double_t        TargetAtomicNumber;
Double_t        TargetNeutronNumber;
Int_t           WhereToStopRecon;
Int_t           SnakeModel;
Double_t        BPMYRes;
Double_t        BPMXRes;
Int_t           UseSeptumPlusStdHRS;
Int_t           UseOpticsDB;
Double_t        TargetXOffset;
Double_t        TargetYOffset;
Double_t        TargetZOffset;
Double_t        PivotXOffset;
Double_t        PivotYOffset;
Double_t        PivotZOffset;
Int_t           SetupLHRS;
Int_t           SetupRHRS;
Double_t        LHRSMomentum;
Double_t        RHRSMomentum;
Double_t        LHRSAngle;
Double_t        RHRSAngle;
Double_t        Pivot2LHRSVBFace;
Double_t        Pivot2RHRSVBFace;
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
Int_t           SetupG2PTarget;
Int_t           TargetType;
Double_t        ThirdArmAngle;
Double_t        SetupThirdArmVD;
Double_t        ThirdArmRotZAngle;
Double_t        Pivot2ThirdArmFace;
Int_t           SetupChicane;
Int_t           SetupChicaneVD;
Double_t        FZB1TiltedAngle;
Double_t        FZB1PosX;
Double_t        FZB1PosY;
Double_t        FZB1PosZ;
Double_t        FZB1Bx;
Double_t        FZB1By;
Double_t        FZB1Bz;
Double_t        FZB2TiltedAngle;
Double_t        FZB2PosX;
Double_t        FZB2PosY;
Double_t        FZB2PosZ;
Double_t        FZB2Bx;
Double_t        FZB2By;
Double_t        FZB2Bz;
Int_t           SetupVD;
Double_t        VDAngle;
Double_t        VDTiltAngle;
Double_t        Pivot2VDFace;
Int_t           SetupBigBite;
Double_t        BigBiteAngle;
Double_t        BigBiteTiltAngle;
Double_t        Pivot2BigBiteFace;
Int_t           SetupHAND;
Double_t        Pivot2HANDLeadWall;
Int_t           SetupSuperBigBite;
Double_t        SuperBigBiteAngle;
Double_t        Pivot2SuperBigBiteFace;
Int_t           UseSBSField;
Double_t        SBSFieldXOffset;
Double_t        SBSFieldYOffset;
Double_t        SBSFieldZOffset;
Double_t        SBSFieldRotAxis1;
Double_t        SBSFieldRotAxis2;
Double_t        SBSFieldRotAxis3;
Double_t        SBSFieldRotAngle1;
Double_t        SBSFieldRotAngle2;
Double_t        SBSFieldRotAngle3;
Double_t        SBSFieldCurrentRatio;
Int_t           SetupHMS;
Double_t        HMSAngle;
Double_t        Pivot2HMSFace;
Int_t           SetupRTPC;
Double_t        RatioHe2DME;
Int_t           SDNum;
Char_t          SDName[20][100];

bool		bIsCombinedTree;	//to ideneify if this is a combined ntuple



bool ReadConfig()
{
  TTree *config = (TTree*)gDirectory->Get("config");
  // Set branch addresses.
  config->SetBranchAddress("Run",&Run);
  config->SetBranchAddress("SkimLevel",&SkimLevel);
  config->SetBranchAddress("BookTrees",&BookTrees);
  config->SetBranchAddress("Beam",&Beam);
  config->SetBranchAddress("BeamTiltedAngle",&BeamTiltedAngle);
  config->SetBranchAddress("TargetM",&TargetM);
  config->SetBranchAddress("TargetL",&TargetL);
  config->SetBranchAddress("TargetAtomicNumber",&TargetAtomicNumber);
  config->SetBranchAddress("TargetNeutronNumber",&TargetNeutronNumber);
  config->SetBranchAddress("WhereToStopRecon",&WhereToStopRecon);
  config->SetBranchAddress("SnakeModel",&SnakeModel);
  config->SetBranchAddress("BPMYRes",&BPMYRes);
  config->SetBranchAddress("BPMXRes",&BPMXRes);
  config->SetBranchAddress("UseSeptumPlusStdHRS",&UseSeptumPlusStdHRS);
  config->SetBranchAddress("UseOpticsDB",&UseOpticsDB);
  config->SetBranchAddress("TargetXOffset",&TargetXOffset);
  config->SetBranchAddress("TargetYOffset",&TargetYOffset);
  config->SetBranchAddress("TargetZOffset",&TargetZOffset);
  config->SetBranchAddress("PivotXOffset",&PivotXOffset);
  config->SetBranchAddress("PivotYOffset",&PivotYOffset);
  config->SetBranchAddress("PivotZOffset",&PivotZOffset);
  config->SetBranchAddress("SetupLHRS",&SetupLHRS);
  config->SetBranchAddress("SetupRHRS",&SetupRHRS);
  config->SetBranchAddress("LHRSMomentum",&LHRSMomentum);
  config->SetBranchAddress("RHRSMomentum",&RHRSMomentum);
  config->SetBranchAddress("LHRSAngle",&LHRSAngle);
  config->SetBranchAddress("RHRSAngle",&RHRSAngle);
  config->SetBranchAddress("Pivot2LHRSVBFace",&Pivot2LHRSVBFace);
  config->SetBranchAddress("Pivot2RHRSVBFace",&Pivot2RHRSVBFace);
  config->SetBranchAddress("UseHelmField",&UseHelmField);
  config->SetBranchAddress("HelmXOffset",&HelmXOffset);
  config->SetBranchAddress("HelmYOffset",&HelmYOffset);
  config->SetBranchAddress("HelmZOffset",&HelmZOffset);
  config->SetBranchAddress("HelmRotAxis1",&HelmRotAxis1);
  config->SetBranchAddress("HelmRotAxis2",&HelmRotAxis2);
  config->SetBranchAddress("HelmRotAxis3",&HelmRotAxis3);
  config->SetBranchAddress("HelmRotAngle1",&HelmRotAngle1);
  config->SetBranchAddress("HelmRotAngle2",&HelmRotAngle2);
  config->SetBranchAddress("HelmRotAngle3",&HelmRotAngle3);
  config->SetBranchAddress("HelmCurrentRatio",&HelmCurrentRatio);
  if(config->GetBranch("UseSeptumField")) 
    {
      config->SetBranchAddress("UseSeptumField",&UseSeptumField);
      config->SetBranchAddress("SeptumXOffset",&SeptumXOffset);
      config->SetBranchAddress("SeptumYOffset",&SeptumYOffset);
      config->SetBranchAddress("SeptumZOffset",&SeptumZOffset);
      config->SetBranchAddress("SeptumRotAxis1",&SeptumRotAxis1);
      config->SetBranchAddress("SeptumRotAxis2",&SeptumRotAxis2);
      config->SetBranchAddress("SeptumRotAxis3",&SeptumRotAxis3);
      config->SetBranchAddress("SeptumRotAngle1",&SeptumRotAngle1);
      config->SetBranchAddress("SeptumRotAngle2",&SeptumRotAngle2);
      config->SetBranchAddress("SeptumRotAngle3",&SeptumRotAngle3);
      config->SetBranchAddress("SeptumCurrentRatioL",&SeptumCurrentRatioL);
      config->SetBranchAddress("SeptumCurrentRatioR",&SeptumCurrentRatioR);
    }
  config->SetBranchAddress("SetupG2PTarget",&SetupG2PTarget);
  config->SetBranchAddress("TargetType",&TargetType);
  config->SetBranchAddress("ThirdArmAngle",&ThirdArmAngle);
  config->SetBranchAddress("SetupThirdArmVD",&SetupThirdArmVD);
  config->SetBranchAddress("ThirdArmRotZAngle",&ThirdArmRotZAngle);
  config->SetBranchAddress("Pivot2ThirdArmFace",&Pivot2ThirdArmFace);
  config->SetBranchAddress("SetupChicane",&SetupChicane);
  config->SetBranchAddress("SetupChicaneVD",&SetupChicaneVD);
  config->SetBranchAddress("FZB1TiltedAngle",&FZB1TiltedAngle);
  config->SetBranchAddress("FZB1PosX",&FZB1PosX);
  config->SetBranchAddress("FZB1PosY",&FZB1PosY);
  config->SetBranchAddress("FZB1PosZ",&FZB1PosZ);
  config->SetBranchAddress("FZB1Bx",&FZB1Bx);
  config->SetBranchAddress("FZB1By",&FZB1By);
  config->SetBranchAddress("FZB1Bz",&FZB1Bz);
  config->SetBranchAddress("FZB2TiltedAngle",&FZB2TiltedAngle);
  config->SetBranchAddress("FZB2PosX",&FZB2PosX);
  config->SetBranchAddress("FZB2PosY",&FZB2PosY);
  config->SetBranchAddress("FZB2PosZ",&FZB2PosZ);
  config->SetBranchAddress("FZB2Bx",&FZB2Bx);
  config->SetBranchAddress("FZB2By",&FZB2By);
  config->SetBranchAddress("FZB2Bz",&FZB2Bz);
  config->SetBranchAddress("SetupVD",&SetupVD);
  config->SetBranchAddress("VDAngle",&VDAngle);
  config->SetBranchAddress("VDTiltAngle",&VDTiltAngle);
  config->SetBranchAddress("Pivot2VDFace",&Pivot2VDFace);
  config->SetBranchAddress("SetupBigBite",&SetupBigBite);
  config->SetBranchAddress("BigBiteAngle",&BigBiteAngle);
  config->SetBranchAddress("BigBiteTiltAngle",&BigBiteTiltAngle);
  config->SetBranchAddress("Pivot2BigBiteFace",&Pivot2BigBiteFace);
  config->SetBranchAddress("SetupHAND",&SetupHAND);
  config->SetBranchAddress("Pivot2HANDLeadWall",&Pivot2HANDLeadWall);
  config->SetBranchAddress("SetupSuperBigBite",&SetupSuperBigBite);
  config->SetBranchAddress("SuperBigBiteAngle",&SuperBigBiteAngle);
  config->SetBranchAddress("Pivot2SuperBigBiteFace",&Pivot2SuperBigBiteFace);
  config->SetBranchAddress("UseSBSField",&UseSBSField);
  config->SetBranchAddress("SBSFieldXOffset",&SBSFieldXOffset);
  config->SetBranchAddress("SBSFieldYOffset",&SBSFieldYOffset);
  config->SetBranchAddress("SBSFieldZOffset",&SBSFieldZOffset);
  config->SetBranchAddress("SBSFieldRotAxis1",&SBSFieldRotAxis1);
  config->SetBranchAddress("SBSFieldRotAxis2",&SBSFieldRotAxis2);
  config->SetBranchAddress("SBSFieldRotAxis3",&SBSFieldRotAxis3);
  config->SetBranchAddress("SBSFieldRotAngle1",&SBSFieldRotAngle1);
  config->SetBranchAddress("SBSFieldRotAngle2",&SBSFieldRotAngle2);
  config->SetBranchAddress("SBSFieldRotAngle3",&SBSFieldRotAngle3);
  config->SetBranchAddress("SBSFieldCurrentRatio",&SBSFieldCurrentRatio);
  config->SetBranchAddress("SetupHMS",&SetupHMS);
  config->SetBranchAddress("HMSAngle",&HMSAngle);
  config->SetBranchAddress("Pivot2HMSFace",&Pivot2HMSFace);
  config->SetBranchAddress("SetupRTPC",&SetupRTPC);
  config->SetBranchAddress("RatioHe2DME",&RatioHe2DME);
  config->SetBranchAddress("SDNum",&SDNum);
  for(int ss=0;ss<20;ss++)
    {
      if(config->GetBranch(Form("SDNameForSDID%d",ss)))
	config->SetBranchAddress(Form("SDNameForSDID%d",ss),SDName[ss]);
    }


  //     This is the loop skeleton
  //       To read only selected branches, Insert statements like:
  // config->SetBranchStatus("*",0);  // disable all branches
  // TTreePlayer->SetBranchStatus("branchname",1);  // activate branchname
  
  Long64_t nentries = config->GetEntries();
  
  bIsCombinedTree=false;
  if(nentries==0) return false;
  Long64_t nbytes = 0;
  double tmpAtg,tmpE,tmpZOff;
  for (Long64_t i=0; i<nentries;i++) 
    {
      nbytes += config->GetEntry(i);
      cout<<"SetupHMS="<<SetupHMS<<endl;
      //check if this is a combined tree with various beam energies
      if(i>0 && (fabs(tmpAtg-TargetAtomicNumber)>0.01 || fabs(tmpE-Beam)>0.01 ||
		 fabs(tmpZOff-PivotZOffset)>30.0 ) )  
	{
	  bIsCombinedTree=true;
	  break;
	}
      else
	{
	  tmpAtg=TargetAtomicNumber;
	  tmpE=Beam;
	  tmpZOff=PivotZOffset;
	}
    }
  return true;
}


void PlotAcc_PhiTheta(const char* key="", int sdid=9, int trackid=0, const char* detector="SBS")
{
  gStyle->SetOptStat(0);
  //set the draw-option "text" format
  gStyle->SetPaintTextFormat(".0f");

  TTree *D = (TTree*)gDirectory->Get("D");
  TCanvas *c1=new TCanvas("c1","",800,600);
  c1->SetRightMargin(0.15);
  c1->cd(0);
  D->Draw(Form("T_Phi[%d]*57.3:T_Theta[%d]*57.3>>hD(25,5,30,40,-60,60)",trackid,trackid),"","");
  D->Draw(Form("T_Phi[%d]*57.3:T_Theta[%d]*57.3>>hN(25,5,30,40,-60,60)",trackid,trackid),
	  Form("SD_Tid==%d && SD_Id==%d",trackid+1,sdid),"");
  TH2D *hD=(TH2D*) gROOT->FindObject("hD");
  TH2D *hN=(TH2D*) gROOT->FindObject("hN");
  TH2D *h2N = (TH2D*) hN->Clone("h2N");
  h2N->Divide(hD);
  h2N->SetTitle(Form("%s Acc in \% ; #theta (deg) ;#phi (deg) ",detector));
  h2N->Scale(100.0);
  h2N->Draw("colz text");

  char strName[100];
  sprintf(strName,"Acc_%s_PhiTheta%s.png",detector,key);
  c1->SaveAs(strName);
  sprintf(strName,"Acc_%s_PhiTheta%s.C",detector,key);
  c1->SaveAs(strName);
}

void PlotAcc_PhiTheta_HMS(const char* key="", int sdid=9, int trackid=0, const char* detector="SBS")
{
  TTree *D = (TTree*)gDirectory->Get("D");

  gStyle->SetOptStat(0);
  //set the draw-option "text" format
  gStyle->SetPaintTextFormat(".0f");

  TCanvas *c1=new TCanvas("c1","",800,600);
  c1->SetRightMargin(0.15);
  c1->cd(0);
  //D->Draw("T_Phi[1]*57.3:T_Theta[1]*57.3>>hD(25,0,50,90,-180,180)","","");
  //D->Draw("T_Phi[1]*57.3:T_Theta[1]*57.3>>hN(25,0,50,90,-180,180)",Form("SD_Tid==2 && SD_Id==%d",sdid),"");
  D->Draw(Form("T_Phi[%d]*57.3:T_Theta[%d]*57.3>>hD(25,5,55,90,180,180)",trackid,trackid),"","");
  D->Draw(Form("T_Phi[%d]*57.3:T_Theta[%d]*57.3>>hN(25,5,55,90,-180,180)",trackid,trackid),
	  Form("SD_Tid==%d && SD_Id==%d",trackid+1,sdid),"");
  TH2D *hD=(TH2D*) gROOT->FindObject("hD");
  TH2D *hN=(TH2D*) gROOT->FindObject("hN");
  TH2D *h2N = (TH2D*) hN->Clone("h2N");
  h2N->Divide(hD);
  h2N->SetTitle(Form("%s Acc in \% ; #theta (deg) ;#phi (deg) ",detector));
  h2N->Scale(100.0);
  h2N->Draw("colz text");

  char strName[100];
  sprintf(strName,"Acc_%s_PhiTheta%s.png",detector,key);
  c1->SaveAs(strName);
  sprintf(strName,"Acc_%s_PhiTheta%s.C",detector,key);
  c1->SaveAs(strName);
}

void PlotAcc_PTheta(const char* key="", int sdid=9, int trackid=0, const char* detector="SBS")
{
  gStyle->SetOptStat(0);
  //set the draw-option "text" format
  gStyle->SetPaintTextFormat(".0f");
  

  TTree *D = (TTree*)gDirectory->Get("D");
  TCanvas *c2=new TCanvas("c2","",800,600);
  c2->cd(0);
  c2->SetRightMargin(0.15);
 
  D->Draw(Form("T_P[%d]*57.3:T_Theta[%d]*57.3>>hDpt(25,5,30,48,0,12)",trackid,trackid),"","");
  D->Draw(Form("T_P[%d]*57.3:T_Theta[%d]*57.3>>hNpt(25,5,30,48,0,12)",trackid,trackid),
	  Form("SD_Tid==%d && SD_Id==%d",trackid+1,sdid),"");
  TH2D *hD=(TH2D*) gROOT->FindObject("hDpt");
  TH2D *hN=(TH2D*) gROOT->FindObject("hNpt");
  TH2D *h2N = (TH2D*) hN->Clone("h2Npt");
  h2N->Divide(hD);
  h2N->SetTitle("RTPC+ABS Acc in \% ; #theta (deg) ; P (GeV/c) ");
  h2N->Scale(100.0);
  h2N->Draw("colz text");

  char strName[100];
  sprintf(strName,"Acc_%s_PTheta%s.png",detector,key);
  c2->SaveAs(strName);
  sprintf(strName,"Acc_%s_PTheta%s.C",detector,key);
  c2->SaveAs(strName);
}


void PlotAcc_RTPC_N_SBS(const char* key="",int sdid=9, int trackid=2)
{
  const char *detector="RTPC_N_SBS";
  PlotAcc_PhiTheta(key,sdid,trackid,detector);
  //PlotAcc_PTheta(key,sdid,trackid,detector);
}


void PlotAcc_LAC_N_SBS()
{
  //const char *detector="NPS_N_HMS";
  int sdid=0, trackid=1;

  PlotAcc_PhiTheta("_SBS_20Deg",sdid=0,trackid=1,"LAC_N_SBS");
  PlotAcc_PhiTheta("_LAC_15Deg",sdid=1,trackid=0,"LAC_N_SBS");
}

void PlotAcc_NPS_N_HMS()
{
  //const char *detector="NPS_N_HMS";
  int sdid=0, trackid=1;

  PlotAcc_PhiTheta_HMS("_HMS_20Deg",sdid=0,trackid=1,"NPS_N_HMS");
  PlotAcc_PhiTheta_HMS("_NPS_15Deg",sdid=1,trackid=0,"NPS_N_HMS");
}

void PlotThrown(const char* inkey="")
{
  system("mkdir -p Graph");
  TTree *D = (TTree*)gDirectory->Get("D");
  if(Run==0) ReadConfig();
  
  TH1F *h1A,*h1B; h1A=h1B=0;
  TH2F *h2A,*h2B; h2A=h2B=0;
  TGraph *grA,*grB; grA=grB=0;
  
  char key[100];
  if(SetupHMS) 
    {
      sprintf(key,"NPS_%.0fdeg_HMS_%.0f_deg%s",VDAngle*57.3,HMSAngle*57.3,inkey);
    }
  else
    {
      sprintf(key,"LAC_%.0fdeg_SBS_%.0f_deg%s",VDAngle*57.3,SuperBigBiteAngle*57.3,inkey);
    }

  TCanvas *c22 = new TCanvas("c22","",600,800);
  c22->Divide(1,2);
  c22->cd(1);
  D->Draw("T_P[0]:T_Theta[0]*57.3>>hPT_g","","");
  //h2A = (TH2F*) gROOT->FindObject("hPT_g");
  //h2A->SetTitle("Thrown RCS Photon; #theta_{#gamma} (deg);E_{#gamma} (GeV)");
  grA = (TGraph*) gROOT->FindObject("Graph")->Clone("GrA");
  grA->SetTitle("Thrown RCS Photon; #theta_{#gamma} (deg);E_{#gamma} (GeV)");
  grA->Draw("AP");

  c22->cd(2);
  D->Draw("T_P[1]:T_Theta[1]*57.3>>hPT_p","","");
  //h2B = (TH2F*) gROOT->FindObject("hPT_p");
  //h2B->SetTitle("Thrown RCS Proton; #theta_{p} (deg);P_{p} (GeV/c)");
  grB = (TGraph*) gROOT->FindObject("Graph")->Clone("GrB");
  grB->SetTitle("Thrown RCS Proton; #theta_{p} (deg);P_{p} (GeV/c)");
  grB->Draw("AP");


  c22->SaveAs(Form("Graph/ThrownDistr_%s.png",key));

  //////////////////////////////////////////////////////////////////

  TCanvas *c21 = new TCanvas("c21","",600,400);

  c21->cd(0);
  D->Draw("T_Theta[0]*57.3:T_Theta[1]*57.3>>hTT_gvsp","","");
  //h2A = (TH2F*) gROOT->FindObject("hTT_gvsp");
  //h2A->SetTitle("Thrown #theta:  Photon Vs Proton; #theta_{p} (deg); #theta_{#gamma} (deg)");
  grA = (TGraph*) gROOT->FindObject("Graph")->Clone("GrTT");
  grA->SetTitle("Thrown #theta:  Photon Vs Proton; #theta_{p} (deg); #theta_{#gamma} (deg)");
  grA->Draw("AP");

  c21->SaveAs(Form("Graph/Theta_gvsp_%s.png",key));

  //////////////////////////////////////////////////////////////////

}



void PlotPhoton(const char* inkey="")
{
  system("mkdir -p Graph");
  //TTree *D = (TTree*)gDirectory->Get("D");
  TTree *track0 = (TTree*)gDirectory->Get("track0");
  if(Run==0) ReadConfig();
  
  TH1F *h1A,*h1B; h1A=h1B=0;
  TH2F *h2A,*h2B, *h2C; h2A=h2B=h2C=0;
  TGraph *grA,*grB; grA=grB=0;
  
  char key[100];
  if(SetupHMS) 
    {
      sprintf(key,"NPS_%.0fdeg_HMS_%.0f_deg%s",VDAngle*57.3,HMSAngle*57.3,inkey);
    }
  else
    {
      sprintf(key,"LAC_%.0fdeg_SBS_%.0f_deg%s",VDAngle*57.3,SuperBigBiteAngle*57.3,inkey);
    }
  char detector[100];
  if(SetupHMS) 
    {
      sprintf(detector,"NPS(%.0f deg)",VDAngle*57.3);
    }
  else
    {
      sprintf(detector,"LAC (%.0f deg)",VDAngle*57.3);
    }


  TCanvas *c22 = new TCanvas("c22","",600,800);
  c22->Divide(1,2);
  c22->cd(1);
  gPad->SetRightMargin(0.15);
  track0->Draw("Phi0*57.3:Theta0*57.3>>hPTdeg_g(20,10,30,50,-50,50)","","colz");
  h2A = (TH2F*) gROOT->FindObject("hPTdeg_g");
  h2A->SetTitle("Thrown RCS Photon; #theta_{#gamma} (deg);#phi_{#gamma} (deg)");
 

  c22->cd(2);
  gPad->SetRightMargin(0.15);
  track0->Draw("Phi0*57.3:Theta0*57.3>>hPTdeg_g_det(20,10,30,50,-50,50)","Pvb>0.2","colz");
  h2B = (TH2F*) gROOT->FindObject("hPTdeg_g_det");
  h2B->SetTitle("Detected RCS Photon; #theta_{#gamma} (deg);#phi_{#gamma} (deg)");

  //c22->SaveAs(Form("Graph/PhotonAcc_raw_%s.png",key));

  //////////////////////////////////////////////////////////////////

  TCanvas *c21 = new TCanvas("c21","",800,600);

  c21->cd();
  gPad->SetRightMargin(0.15);
  h2C = (TH2F*) h2B->Clone("acc_g");
  h2C->Divide(h2A);
  h2C->Scale(100.0);
  h2C->SetTitle(Form("%s Photon Acc in \% ; #theta_{#gamma} (deg);#phi_{#gamma} (deg)",detector));
  h2C->Draw("colz text");
  c21->SaveAs(Form("Graph/PhotonAcc_%s.png",key));

  //////////////////////////////////////////////////////////////////
}

void PlotElectron(const char* inkey="")
{
  system("mkdir -p Graph");
  //TTree *D = (TTree*)gDirectory->Get("D");
  TTree *track2 = (TTree*)gDirectory->Get("track2");
  if(Run==0) ReadConfig();
  
  TH1F *h1A,*h1B; h1A=h1B=0;
  TH2F *h2A,*h2B, *h2C; h2A=h2B=h2C=0;
  TGraph *grA,*grB; grA=grB=0;
  
  char key[100];
  if(SetupHMS) 
    {
      sprintf(key,"NPS_%.0fdeg_HMS_%.0f_deg%s",VDAngle*57.3,HMSAngle*57.3,inkey);
    }
  else
    {
      sprintf(key,"LAC_%.0fdeg_SBS_%.0f_deg%s",VDAngle*57.3,SuperBigBiteAngle*57.3,inkey);
    }
  char detector[100];
  if(SetupHMS) 
    {
      sprintf(detector,"NPS(%.0f deg)",VDAngle*57.3);
    }
  else
    {
      sprintf(detector,"LAC (%.0f deg)",VDAngle*57.3);
    }


  TCanvas *c22 = new TCanvas("c22","",600,800);
  c22->Divide(1,2);
  c22->cd(1);
  gPad->SetRightMargin(0.15);
  track2->Draw("Phi0*57.3:Theta0*57.3>>hPTdeg_e(20,10,30,50,-50,50)","","colz");
  h2A = (TH2F*) gROOT->FindObject("hPTdeg_e");
  h2A->SetTitle("Thrown RCS Electron; #theta_{e} (deg);#phi_{e} (deg)");
 

  c22->cd(2);
  gPad->SetRightMargin(0.15);
  track2->Draw("Phi0*57.3:Theta0*57.3>>hPTdeg_e_det(20,10,30,50,-50,50)","Pvb>0.2","colz");
  h2B = (TH2F*) gROOT->FindObject("hPTdeg_e_det");
  h2B->SetTitle("Detected RCS Electron; #theta_{e} (deg);#phi_{e} (deg)");

  //c22->SaveAs(Form("Graph/ElectronAcc_raw_%s.png",key));

  //////////////////////////////////////////////////////////////////

  TCanvas *c21 = new TCanvas("c21","",800,600);

  c21->cd();
  gPad->SetRightMargin(0.15);
  h2C = (TH2F*) h2B->Clone("acc_e");
  h2C->Divide(h2A);
  h2C->Scale(100.0);
  h2C->SetTitle(Form("%s Electron Acc in \% ; #theta_{e} (deg);#phi_{e} (deg)",detector));
  h2C->Draw("colz text");
  c21->SaveAs(Form("Graph/ElectronAcc_%s.png",key));

  //////////////////////////////////////////////////////////////////
}

void PlotProton(const char* inkey="")
{
  system("mkdir -p Graph");
  //TTree *D = (TTree*)gDirectory->Get("D");
  TTree *track1 = (TTree*)gDirectory->Get("track1");
  if(Run==0) ReadConfig();
  
  TH1F *h1A,*h1B; h1A=h1B=0;
  TH2F *h2A,*h2B, *h2C; h2A=h2B=h2C=0;
  TGraph *grA,*grB; grA=grB=0;
  
  char key[100];
  if(SetupHMS) 
    {
      sprintf(key,"NPS_%.0fdeg_HMS_%.0f_deg%s",VDAngle*57.3,HMSAngle*57.3,inkey);
    }
  else
    {
      sprintf(key,"LAC_%.0fdeg_SBS_%.0f_deg%s",VDAngle*57.3,SuperBigBiteAngle*57.3,inkey);
    }
  char detector[100];
  if(SetupHMS) 
    {
      sprintf(detector,"HMS(%.0f deg)",HMSAngle*57.3);
    }
  else
    {
      sprintf(detector,"SBS(%.0f deg)",SuperBigBiteAngle*57.3);
    }

  TCanvas *c22 = new TCanvas("c22","",600,800);
  c22->Divide(1,2);
  c22->cd(1);
  gPad->SetRightMargin(0.15);
  track1->Draw("fmod(Phi0*57.3+360,360.):Theta0*57.3>>hPTdeg_p(20,20,40,40,160,200)","","colz");
  h2A = (TH2F*) gROOT->FindObject("hPTdeg_p");
  h2A->SetTitle("Thrown RCS Proton; #theta_{p} (deg);#phi_{p} (deg)");
 

  c22->cd(2);
  gPad->SetRightMargin(0.15);
  track1->Draw("fmod(Phi0*57.3+360,360.):Theta0*57.3>>hPTdeg_p_det(20,20,40,40,160,200)","Pvb>0.2","colz");
  h2B = (TH2F*) gROOT->FindObject("hPTdeg_p_det");
  h2B->SetTitle("Detected RCS Proton; #theta_{p} (deg);#phi_{p} (deg)");

  //c22->SaveAs(Form("Graph/ProtonAcc_raw_%s.png",key));

  //////////////////////////////////////////////////////////////////

  TCanvas *c21 = new TCanvas("c21","",800,600);

  c21->cd(0);
  gPad->SetRightMargin(0.15);
  h2C = (TH2F*) h2B->Clone("acc_p");
  h2C->Divide(h2A);
  h2C->Scale(100.0);
  h2C->SetTitle(Form("%s Proton Acc in \% ; #theta_{p} (deg);#phi_{p} (deg)",detector));
  h2C->Draw("colz text");
  c21->SaveAs(Form("Graph/ProtonAcc_%s.png",key));

  //////////////////////////////////////////////////////////////////

}


void PlotRCS(const char* inkey="")
{
  system("mkdir -p Graph");
  //TTree *D = (TTree*)gDirectory->Get("D");
  TTree *track0 = (TTree*)gDirectory->Get("track0");
  track0->AddFriend("track1");

  if(Run==0) ReadConfig();
  
  TH1F *h1A,*h1B; h1A=h1B=0;
  TH2F *h2A,*h2B, *h2C; h2A=h2B=h2C=0;
  TGraph *grA,*grB,*grC,*grD; grA=grB=grC=grD=0;
  
  char key[100];
  if(SetupHMS) 
    {
      sprintf(key,"NPS_%.0fdeg_HMS_%.0f_deg%s",VDAngle*57.3,HMSAngle*57.3,inkey);
    }
  else
    {
      sprintf(key,"LAC_%.0fdeg_SBS_%.0f_deg%s",VDAngle*57.3,SuperBigBiteAngle*57.3,inkey);
    }
  char detector[100];
  if(SetupHMS) 
    {
      sprintf(detector,"HMS(%.0f deg)+NPS(%.0f deg)",HMSAngle*57.3,VDAngle*57.3);
    }
  else
    {
      sprintf(detector,"SBS(%.0f deg)+LAC(%.0f deg)",SuperBigBiteAngle*57.3,VDAngle*57.3);
    }

  TCanvas *c31 = new TCanvas("c31","",800,600);
  //c31->Divide(1,2);
  //c31->cd(1);
  //gPad->SetRightMargin(0.15);
  char theTg[255],theBin[255];
  char strTg[4][255];
  const char *strCut[4]={"",
			 "track1.Pvb>0",
			 "track0.Pvb>0",
			 "track1.Pvb>0 && track1.Pvb>0"
  };
  const char *strTitle[4]={
    "Thrown RCS photon; #theta_{#gamma} (deg);#phi_{#gamma} (deg)",
    "Thrown RCS photon: proton detected; #theta_{#gamma} (deg);#phi_{#gamma} (deg)",
    "Detected RCS photon; #theta_{#gamma} (deg);#phi_{#gamma} (deg)",
    "Both photon and proton detected; #theta_{#gamma} (deg);#phi_{#gamma} (deg)"
  };
  if(SetupHMS) 
    {
      sprintf(theTg,"track0.Phi0*57.3:track0.Theta0*57.3");
      sprintf(theBin,"(20,6,26,40,-40,40)");
    }
  else
    {
      sprintf(theTg,"fmod(track0.Phi0*57.3+360,360.):track0.Theta0*57.3");
      sprintf(theBin,"(25,5,55,40,100,260)");
    }
  sprintf(strTg[0],"%s >> h2PT_acc_A%s",theTg,theBin);
  sprintf(strTg[1],"%s >> h2PT_acc_B%s",theTg,theBin);
  sprintf(strTg[2],"%s >> h2PT_acc_C%s",theTg,theBin);
  sprintf(strTg[3],"%s >> h2PT_acc_D%s",theTg,theBin);
  
  //track0->Draw("fmod(track0.Phi0*57.3+360,360.):track0.Theta0*57.3","");
  track0->Draw(strTg[0],strCut[0],"");
  h2A = (TH2F*) gROOT->FindObject("h2PT_acc_A")->Clone("h2A");
  h2A->SetTitle("Thrown photon; #theta_{#gamma} (deg);#phi_{#gamma} (deg)");
  grA = (TGraph*) gROOT->FindObject("Graph")->Clone("GrA");
  grA->SetTitle(strTitle[0]);
  grA->SetMarkerStyle(1); 
  grA->SetMarkerColor(12);
  grA->Draw("AP");
      
  //track0->Draw("fmod(track0.Phi0*57.3+360,360.):track0.Theta0*57.3","track1.Pvb>0 ","same");
  track0->Draw(strTg[1],strCut[1],"");
  grB = (TGraph*) gROOT->FindObject("Graph")->Clone("GrB");
  //grB->SetTitle("Thrown Photon: proton detected; #theta_{#gamma} (deg);#phi_{#gamma} (deg)");
  grB->SetTitle(strTitle[1]);
  grB->SetMarkerStyle(1); 
  grB->SetMarkerColor(9);

  //track0->Draw("fmod(track0.Phi0*57.3+360,360.):track0.Theta0*57.3","track0.Pvb>0","same");
  track0->Draw(strTg[2],strCut[2],"");
  grC = (TGraph*) gROOT->FindObject("Graph")->Clone("GrC");
  //grC->SetTitle("Photon detected; #theta_{#gamma} (deg);#phi_{#gamma} (deg)");
  grC->SetTitle(strTitle[2]);
  grC->SetMarkerStyle(1); 
  grC->SetMarkerColor(4);

  //track0->Draw("fmod(track0.Phi0*57.3+360,360.):track0.Theta0*57.3","track1.Pvb>0 && track0.Pvb>0","same");
  track0->Draw(strTg[3],strCut[3],"");
  h2B = (TH2F*) gROOT->FindObject("h2PT_acc_D")->Clone("h2B");
  h2B->SetTitle("Detect both photon and proton; #theta_{#gamma} (deg);#phi_{#gamma} (deg)");
  grD = (TGraph*) gROOT->FindObject("Graph")->Clone("GrD");
  grD->SetTitle(strTitle[3]);
  grD->SetMarkerStyle(1); 
  grD->SetMarkerColor(2);

  grA->Draw("AP");
  grB->Draw("Psame");
  grC->Draw("Psame");
  grD->Draw("Psame");

  c31->SaveAs(Form("Graph/Acc_gp_raw_%s.png",key));

  //////////////////////////////////////////////////////////////////

  TCanvas *c32 = new TCanvas("c32","",800,600);

  c32->cd(0);
  gPad->SetRightMargin(0.15);
  h2C = (TH2F*) h2B->Clone("acc_g");
  h2C->Divide(h2A);
  h2C->Scale(100.0);
  h2C->SetTitle(Form("%s #gamma-p coincident Acc (\%) ; #theta_{#gamma} (deg);#phi_{#gamma} (deg)",detector));
  h2C->Draw("colz text");
  if(SetupHMS)
    {
      h2C->GetXaxis()->SetRangeUser(10,30);
      h2C->GetYaxis()->SetRangeUser(-20,20);
    }
  else
    {
      h2C->GetXaxis()->SetRangeUser(15,35);
      h2C->GetXaxis()->SetRangeUser(140,230);
    }
  h2C->Draw("colz text");
  c32->Modified();

  c32->SaveAs(Form("Graph/Acc_gp_%s.png",key));

 //////////////////////////////////////////////////////////////////
  //plot RCS coincident events
  TGraph *grE;
  TCanvas *c33 = new TCanvas("c33","",600,800);
  c33->Divide(1,2);

  TCut RCSCut = "track1.Pvb/track1.P0>0.9 && track0.Pvb/track0.P0>0.9";
  c33->cd(1);
  track0->Draw("track0.Pvb:track0.Theta0*57.3 >> PTheta_g",RCSCut,"");
  grE = (TGraph*) gROOT->FindObject("Graph")->Clone("GrE_1");
  grE->SetTitle(Form("%s: #gamma-p coincident; #theta_{#gamma} (deg);E_{#gamma} (GeV)",detector));
  grE->Draw("AP");

  c33->cd(2);
  track0->Draw("track1.Pvb:track1.Theta0*57.3 >> PTheta_p",RCSCut,"");
  grE = (TGraph*) gROOT->FindObject("Graph")->Clone("GrE_2");
  grE->SetTitle(Form("%s: #gamma-p coincident; #theta_{p} (deg);P_{p} (GeV)",detector));
  grE->Draw("AP");
  c33->Modified();
  c33->SaveAs(Form("Graph/PTheta_gp_%s.png",key));
}


void PlotEP(const char* inkey="")
{
  system("mkdir -p Graph");
  //TTree *D = (TTree*)gDirectory->Get("D");
  TTree *track2 = (TTree*)gDirectory->Get("track2");
  track2->AddFriend("track1");

  if(Run==0) ReadConfig();
  
  TH1F *h1A,*h1B; h1A=h1B=0;
  TH2F *h2A,*h2B, *h2C; h2A=h2B=h2C=0;
  TGraph *grA,*grB,*grC,*grD; grA=grB=grC=grD=0;
  
  char key[100];
  if(SetupHMS) 
    {
      sprintf(key,"NPS_%.0fdeg_HMS_%.0f_deg%s",VDAngle*57.3,HMSAngle*57.3,inkey);
    }
  else
    {
      sprintf(key,"LAC_%.0fdeg_SBS_%.0f_deg%s",VDAngle*57.3,SuperBigBiteAngle*57.3,inkey);
    }
  char detector[100];
  if(SetupHMS) 
    {
      sprintf(detector,"HMS(%.0f deg)+NPS(%.0f deg)",HMSAngle*57.3,VDAngle*57.3);
    }
  else
    {
      sprintf(detector,"SBS(%.0f deg)+LAC(%.0f deg)",SuperBigBiteAngle*57.3,VDAngle*57.3);
    }

  TCanvas *c31 = new TCanvas("c31","",800,600);
  //c31->Divide(1,2);
  //c31->cd(1);
  //gPad->SetRightMargin(0.15);
  char theTg[255],theBin[255];
  char strTg[4][255];
  const char *strCut[4]={"",
			 "track1.Pvb>0",
			 "track2.Pvb>0",
			 "track1.Pvb>0 && track1.Pvb>0"
  };
  const char *strTitle[4]={
    "Thrown electron; #theta_{e} (deg);#phi_{e} (deg)",
    "Thrown electron: proton detected; #theta_{e} (deg);#phi_{e} (deg)",
    "Detected electron; #theta_{e} (deg);#phi_{e} (deg)",
    "Both electron and proton detected; #theta_{e} (deg);#phi_{e} (deg)"
  };
  if(SetupHMS) 
    {
      sprintf(theTg,"track2.Phi0*57.3:track2.Theta0*57.3");
      sprintf(theBin,"(25,5,30,50,-50,50)");
    }
  else
    {
      sprintf(theTg,"fmod(track2.Phi0*57.3+360,360.):track2.Theta0*57.3");
      sprintf(theBin,"(25,6,56,40,100,260)");
    }
  sprintf(strTg[0],"%s >> h2PT_e_acc_A%s",theTg,theBin);
  sprintf(strTg[1],"%s >> h2PT_e_acc_B%s",theTg,theBin);
  sprintf(strTg[2],"%s >> h2PT_e_acc_C%s",theTg,theBin);
  sprintf(strTg[3],"%s >> h2PT_e_acc_D%s",theTg,theBin);
  
  //track2->Draw("fmod(track2.Phi0*57.3+360,360.):track2.Theta0*57.3","");
  track2->Draw(strTg[0],strCut[0],"");
  h2A = (TH2F*) gROOT->FindObject("h2PT_e_acc_A")->Clone("h2A");
  h2A->SetTitle("Thrown electron; #theta_{e} (deg);#phi_{e} (deg)");
  grA = (TGraph*) gROOT->FindObject("Graph")->Clone("GrA");
  grA->SetTitle(strTitle[0]);
  grA->SetMarkerStyle(1); 
  grA->SetMarkerColor(12);
  grA->Draw("AP");
      
  //track2->Draw("fmod(track2.Phi0*57.3+360,360.):track2.Theta0*57.3","track1.Pvb>0 ","same");
  track2->Draw(strTg[1],strCut[1],"");
  grB = (TGraph*) gROOT->FindObject("Graph")->Clone("GrB");
  //grB->SetTitle("Thrown Electron: proton detected; #theta_{e} (deg);#phi_{e} (deg)");
  grB->SetTitle(strTitle[1]);
  grB->SetMarkerStyle(1); 
  grB->SetMarkerColor(9);

  //track2->Draw("fmod(track2.Phi0*57.3+360,360.):track2.Theta0*57.3","track2.Pvb>0","same");
  track2->Draw(strTg[2],strCut[2],"");
  grC = (TGraph*) gROOT->FindObject("Graph")->Clone("GrC");
  //grC->SetTitle("Electron detected; #theta_{e} (deg);#phi_{e} (deg)");
  grC->SetTitle(strTitle[2]);
  grC->SetMarkerStyle(1); 
  grC->SetMarkerColor(4);

  //track2->Draw("fmod(track2.Phi0*57.3+360,360.):track2.Theta0*57.3","track1.Pvb>0 && track2.Pvb>0","same");
  track2->Draw(strTg[3],strCut[3],"");
  h2B = (TH2F*) gROOT->FindObject("h2PT_e_acc_D")->Clone("h2B");
  h2B->SetTitle("Detect both electron and proton; #theta_{e} (deg);#phi_{e} (deg)");
  grD = (TGraph*) gROOT->FindObject("Graph")->Clone("GrD");
  grD->SetTitle(strTitle[3]);
  grD->SetMarkerStyle(1); 
  grD->SetMarkerColor(2);

  grA->Draw("AP");
  grB->Draw("Psame");
  grC->Draw("Psame");
  grD->Draw("Psame");

  c31->SaveAs(Form("Graph/Acc_ep_raw_%s.png",key));

  //////////////////////////////////////////////////////////////////

  TCanvas *c32 = new TCanvas("c32","",800,600);

  c32->cd(0);
  gPad->SetRightMargin(0.15);
  h2C = (TH2F*) h2B->Clone("acc_g");
  h2C->Divide(h2A);
  h2C->Scale(100.0);
  h2C->SetTitle(Form("%s e-p coincident Acc(\%) ; #theta_{e} (deg);#phi_{e} (deg)",detector));
  h2C->Draw("colz text");
  if(SetupHMS)
    {
      h2C->GetXaxis()->SetRangeUser(10,30);
      h2C->GetYaxis()->SetRangeUser(-20,20);
    }
  else
    {
      h2C->GetXaxis()->SetRangeUser(10,40);
      h2C->GetXaxis()->SetRangeUser(132,232);
    }
  h2C->Draw("colz text");
  c32->Modified();

  c32->SaveAs(Form("Graph/Acc_ep_%s.png",key));

 //////////////////////////////////////////////////////////////////
}

void PlotEGDiff(const char* inkey="")
{
  system("mkdir -p Graph");
  //TTree *D = (TTree*)gDirectory->Get("D");
  TTree *track1 = (TTree*)gDirectory->Get("track1");
  track1->AddFriend("track0");
  track1->AddFriend("track2");

  if(Run==0) ReadConfig();
  
  TH1F *h1A,*h1B; h1A=h1B=0;
  TH2F *h2A,*h2B, *h2C; h2A=h2B=h2C=0;
  TGraph *grA,*grB,*grC; grA=grB=grC=0;
  
  char key[100]; 
  char theCut[255];
  if(SetupHMS) 
    {
      sprintf(key,"NPS_%.0fdeg_HMS_%.0f_deg%s",VDAngle*57.3,HMSAngle*57.3,inkey);
    }
  else
    {
      sprintf(key,"LAC_%.0fdeg_SBS_%.0f_deg%s",VDAngle*57.3,SuperBigBiteAngle*57.3,inkey);
    }
  char detector[100], strPhi[100];
  if(SetupHMS) 
    {
      sprintf(detector,"HMS(%.0f deg)",HMSAngle*57.3);
      sprintf(strPhi,"track0.Phi0*57.3");
      sprintf(theCut,"track0.Pvb/track0.P0>0.9 && track1.Pvb/track1.P0>0.9 && track2.Pvb/track2.P0>0.9 && track0.Yvb<track2.Yvb");
    }
  else
    {
      sprintf(detector,"SBS(%.0f deg)",SuperBigBiteAngle*57.3);
      sprintf(strPhi,"fmod(track0.Phi0*57.3+360,360)");
      sprintf(theCut,"track0.Pvb/track0.P0>0.9 && track1.Pvb/track1.P0>0.9 && track2.Pvb/track2.P0>0.9 && track0.Yvb>track2.Yvb");
    }

  TCanvas *c42 = new TCanvas("c22","",400,900);

  c42->Divide(1,3);
  c42->cd(1);
  track1->Draw(Form("track0.Yvb-track2.Yvb:%s>>h2dYPhi0",strPhi),theCut,"");
  grA = (TGraph*) gROOT->FindObject("Graph")->Clone("h2dYPhi0");
  grA->SetTitle("Vertical offset of electron and photon; #phi_{#gamma} (deg);y_{#gamma}-y_{e} (mm)");
  grA->Draw("AP");

  c42->cd(2);
  track1->Draw("track0.Yvb-track2.Yvb:track0.P0>>h2dYP0",theCut,"");
  grB = (TGraph*) gROOT->FindObject("Graph")->Clone("h2dYP0");
  grB->SetTitle("Vertical offset of electron and photon; P_{#gamma} (GeV/c);y_{#gamma}-y_{e} (mm)");
  grB->Draw("AP");
 
  c42->cd(3);
  track1->Draw("track0.Yvb-track2.Yvb:track0.Theta0*57.3>>h2dYTheta0",theCut,"");
  grC = (TGraph*) gROOT->FindObject("Graph")->Clone("h2dYTheta0");
  grC->SetTitle("Vertical offset of electron and photon; #theta_{#gamma} (deg);y_{#gamma}-y_{e} (mm)");
  grC->Draw("AP");

  c42->Modified();
  c42->SaveAs(Form("Graph/dY_eg_%s.png",key));
}


void PlotWACS()
{
  gStyle->SetOptStat(0);
  //set the draw-option "text" format
  gStyle->SetPaintTextFormat(".0f");
  system("mkdir -p Graph");
  ReadConfig();

  PlotThrown();
  PlotPhoton();
  PlotProton();
  PlotElectron();

  PlotRCS();
  PlotEP();
  //return;
  PlotEGDiff();
  //return;
}
