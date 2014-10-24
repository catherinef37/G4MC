//plot mobitoring histo for G2P G4Sim
#include "TROOT.h"
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

#include "math.h"
#include <string>
#include <iostream>

const double deg = asin(1.0)/90.;

//Declaration of leaves types tree
Int_t           Run;
Int_t           SkimLevel;
Int_t           BookTrees;
Double_t        Beam;
Double_t        BeamTiltedAngle;
Double_t        TargetM;
Double_t        TargetL=2.54;
Double_t        TargetAtomicNumber;
Double_t        TargetNeutronNumber;
Int_t           SetupG2PTarget;
Int_t           TargetType;
Double_t        TargetXOffset;
Double_t        TargetYOffset;
Double_t        TargetZOffset;
Double_t        PivotXOffset;
Double_t        PivotYOffset;
Double_t        PivotZOffset;
Double_t        LHRSMomentum;
Double_t        RHRSMomentum;
Double_t        LHRSAngle;
Double_t        RHRSAngle;
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
Double_t        HelmCurrentRatio=1.0;
Int_t           UseSeptumField;
Int_t           UseSeptumPlusStdHRS;
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
Double_t        ThirdArmAngle;
Double_t        SetupThirdArmVD;
Double_t        ThirdArmRotZAngle;
Double_t        Pivot2ThirdArmFace;
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
  config->SetBranchAddress("SetupG2PTarget",&SetupG2PTarget);
  config->SetBranchAddress("TargetType",&TargetType);
  config->SetBranchAddress("TargetXOffset",&TargetXOffset);
  config->SetBranchAddress("TargetYOffset",&TargetYOffset);
  config->SetBranchAddress("TargetZOffset",&TargetZOffset);
  config->SetBranchAddress("PivotXOffset",&PivotXOffset);
  config->SetBranchAddress("PivotYOffset",&PivotYOffset);
  config->SetBranchAddress("PivotZOffset",&PivotZOffset);
  config->SetBranchAddress("LHRSMomentum",&LHRSMomentum);
  config->SetBranchAddress("RHRSMomentum",&RHRSMomentum);
  config->SetBranchAddress("LHRSAngle",&LHRSAngle);
  config->SetBranchAddress("RHRSAngle",&RHRSAngle);
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
  config->SetBranchAddress("UseSeptumField",&UseSeptumField);
  config->SetBranchAddress("UseSeptumPlusStdHRS",&UseSeptumPlusStdHRS);
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
  config->SetBranchAddress("BigBiteAngle",&BigBiteAngle);
  config->SetBranchAddress("BigBiteTiltAngle",&BigBiteTiltAngle);
  config->SetBranchAddress("Pivot2BigBiteFace",&Pivot2BigBiteFace);
  config->SetBranchAddress("ThirdArmAngle",&ThirdArmAngle);
  config->SetBranchAddress("SetupThirdArmVD",&SetupThirdArmVD);
  config->SetBranchAddress("ThirdArmRotZAngle",&ThirdArmRotZAngle);
  config->SetBranchAddress("Pivot2ThirdArmFace",&Pivot2ThirdArmFace);


  Long64_t nentries = config->GetEntries();

  bIsCombinedTree=false;
  if(nentries==0) return false;
  Long64_t nbytes = 0;
  double tmpAtg,tmpE,tmpZOff;
  for (Long64_t i=0; i<nentries;i++) 
    {
      nbytes += config->GetEntry(i);
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


//cutshape: 1 rectangle; 2 ellipse
void Monitor(char *filename="", int trackindex=0, int cutshape=2)
{
  gStyle->SetStatW(0.28);
  gStyle->SetStatH(0.06);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);

  if(strlen(filename)>5)
    {
      cout<<"trying to open root file "<<filename<<endl;
      TFile* _file0=TFile::Open(filename);
    }
  ReadConfig();
 char strTg[50];

  if(TargetAtomicNumber==1 && TargetNeutronNumber==0) sprintf(strTg,"LH2");
  else if(TargetAtomicNumber==1 && TargetNeutronNumber==1) sprintf(strTg,"LD2"); 
  else if(TargetAtomicNumber==2 && TargetNeutronNumber==2) sprintf(strTg,"LHe"); 
  else if(TargetAtomicNumber==6 && TargetNeutronNumber==6) sprintf(strTg,"C12");   
  else if(TargetAtomicNumber==13 && TargetNeutronNumber==14) sprintf(strTg,"Al"); 
  else if(TargetAtomicNumber==1) sprintf(strTg,"NH3+LHe");
  else  
    sprintf(strTg,"^{%d}A%d",TargetAtomicNumber+TargetNeutronNumber,TargetAtomicNumber);

  TTree *tree = (TTree*)gDirectory->Get(Form("track%d",trackindex));
  if(trackindex>0) tree->AddFriend("track0");

  char Cut[1000],XSCut[1000];
  double kVHigh=0.045,kVLow=-0.045,kHHigh=0.010,kHLow=-0.022, kECut=0.04;
  if(fabs(LHRSAngle)>10.0*deg)
    {
      kVHigh=0.05,kVLow=-0.05;
      kHHigh=0.026;kHLow=-0.026;
    }
  //From Alan
  kVHigh=0.05,kVLow=-0.05;
  kHHigh=0.026;kHLow=-0.026;

  TCanvas *pCan = new TCanvas("pCan","Pass Septum|Q1 ACC Cut",1200,800);
  pCan->Divide(4,3,0.001,0.001);

  //try to identify if it is a good HRS track
  double MaxXS=10000.0;
  int IsHRSTrack=0;		
  sprintf(Cut,"TrackClass>3");
  pCan->cd(1);tree->Draw("ElasXS>>hElasXS_T1",Cut);
  if(hElasXS_T1->GetEntries()>1)
    {
      MaxXS=hElasXS_T1->GetMean()*40.0;
      IsHRSTrack=2;
    }
  else
    {
      //cutshape: 1 rectangle; 2 ellipse
      if (cutshape==1)
	{
	  sprintf(Cut,"TrackClass>1 && fabs(Phivb_tr-(%.3f))<%.3f && fabs(Thetavb_tr-(%.3f))<%.3f && (P0-Pvb)/P0<%.3f",
		  (kHHigh+kHLow)/2.0,(kHHigh-kHLow)/2.0,(kVHigh+kVLow)/2.0,(kVHigh-kVLow)/2.0,kECut);
	}
      else
	{
	  sprintf(Cut,"TrackClass>1 && pow((Phivb_tr-(%.3f))/%.3f,2.0)+pow((Thetavb_tr-(%.3f))/%.3f,2.0)<1.0 && (P0-Pvb)/P0<%.3f",
		  (kHHigh+kHLow)/2.0,(kHHigh-kHLow)/2.0,(kVHigh+kVLow)/2.0,(kVHigh-kVLow)/2.0,kECut);
	}
      tree->Draw("ElasXS>>hElasXS_T2",Cut);
      if(hElasXS_T2->GetEntries()>1)
	{
	  MaxXS=hElasXS_T2->GetMean()*40.0;
	  IsHRSTrack=1;
	}
      else
	{
	  kECut=0.50;
	  sprintf(Cut,Form("TrackClass>=0 && (P0-Pvb)/P0<%.3f",kECut));
	}
    }

  if(MaxXS<100) MaxXS=100.0;
  sprintf(XSCut,"%s && ElasXS<%.2f",Cut,MaxXS);

  char Cut_reject[]="TrackClass<=0"; 

  //Draw Theta0:Phi0
  pCan->cd(1);
  tree->Draw("Y0:X0>>h2Y0X0",Cut,"");
  h2Y0X0->SetTitle("Y0:X0 (cut);X0 (mm);Y0 (mm)");
  h2Y0X0->SetMarkerColor(2);
  h2Y0X0->Draw(); //need to draw again to show the color and title
  gPad->Modified();

  pCan->cd(2);
  tree->Draw("Z0 >> hZ0","","");
  hZ0->SetTitle("Z0: all(black) and cut(red); Z0 (mm)");
  tree->Draw("Z0 >> hZ0Cut",Cut,"same");
  hZ0Cut->SetLineColor(2);
  gPad->Modified();

  pCan->cd(3);
  tree->Draw("Theta0_tr:Phi0_tr>>h2Theta0_trPhi0_tr_all",
	     "abs(Theta0_tr)<0.5 && abs(Phi0_tr)<0.12","");
  h2Theta0_trPhi0_tr_all->SetTitle("Theta0_tr:Phi0_tr all (black) and cut(red);Phi0_tr (rad);Theta0_tr (rad)");
  tree->Draw("Theta0_tr:Phi0_tr>>h2Theta0_trPhi0_tr_cut",Cut,"same");
  //h2Theta0_trPhi0_tr_cut->SetTitle("Theta0_tr:Phi0_tr (cut)");
  h2Theta0_trPhi0_tr_cut->SetMarkerColor(2);
  h2Theta0_trPhi0_tr_all->Draw();
  h2Theta0_trPhi0_tr_cut->Draw("same");
  gPad->Modified();

  pCan->cd(4);
  tree->Draw("Theta0:Phi0>>h2Theta0Phi0","","");
  h2Theta0Phi0->SetTitle("Theta0:Phi0, all(black) and cut(red);Phi0 (rad);Theta0 (rad)");
  tree->Draw("Theta0:Phi0 >> h2Theta0Phi0Cut",Cut,"same");
  //h2Theta0Phi0Cut->SetTitle("Theta0:Phi0 (cut)");
  h2Theta0Phi0Cut->SetMarkerColor(2);
  h2Theta0Phi0->Draw();
  h2Theta0Phi0Cut->Draw("same");
  double pMeanTheta0=h2Theta0Phi0Cut->GetMean(2);
  gPad->Modified();

  pCan->cd(5);
  if(pMeanTheta0/deg<30.)
    {
      tree->Draw("-Xvb_tr:-Yvb_tr>>h2YvbXvb","TrackClass>=0 && abs(Xvb_tr)<150 && abs(Yvb_tr)<150","");
      h2YvbXvb->GetYaxis()->SetRangeUser(-150,150);
    }
  else
    {
      tree->Draw("-Xvb_tr:Yvb_tr>>h2YvbXvb","TrackClass>=0","");
    }
  h2YvbXvb->SetTitle(" VD front view: all(black) and cut(blue);-Yvb_tr (mm);-Xvb_tr (mm)");
  tree->Draw("-Xvb_tr:-Yvb_tr>>h2YvbXvbCut",Cut,"same");
  h2YvbXvbCut->SetMarkerColor(4);
  h2YvbXvb->Draw();
  h2YvbXvbCut->Draw("same");
	
	
  pCan->cd(6);
  
  tree->Draw("Zvb_tr>> hZvbtr", Cut,"");
  double pMeanZvbtr=hZvbtr->GetMean();
		 
  tree->Draw("Thetavb_tr:Phivb_tr >> h2TPvbReject",
	     "abs(Thetavb_tr)<0.12 && abs(Phivb_tr)<0.12","");
  //h2TPvbReject->SetTitle("Thetavb_tr:Phivb_tr all{black};Phivb_tr (rad);Thetavb_tr (rad)");
  h2TPvbReject->GetXaxis()->SetRangeUser(-0.12,0.12);
  h2TPvbReject->GetYaxis()->SetRangeUser(-0.12,0.12);
  h2TPvbReject->SetTitle("Thetavb_tr:Phivb_tr all(black) and cut(blue);Phivb_tr;Thetavb_tr");
  tree->Draw("Thetavb_tr:Phivb_tr >> h2TPvbCut",Cut,"same");
  h2TPvbCut->SetMarkerColor(4);
  h2TPvbReject->Draw();
  h2TPvbCut->Draw("same");
  gPad->Modified();


  pCan->cd(8);
  TPaveText *ptInfo = new TPaveText(0.02,0.02,0.98,0.98,"blNDC");
  ptInfo->SetName("info");
  ptInfo->SetBorderSize(1);
  ptInfo->SetFillColor(0);
  ptInfo->SetTextAlign(12);
  TText *text=0;
  char line[3][100];
  text=ptInfo->AddText("All '_tr' means in TCS, otherwise in HCS");

  sprintf(line[0],"#splitline{ #splitline{All '0' means measured at vertex, }{");
  sprintf(line[1],"'vb' means at virtual boundary, which is} }{");
  //sprintf(line[2]," #color[2]{%.1f} mm from the target}",pMeanZvbtr);
  sprintf(line[2],"#color[2]{%s} collimator.}",
	  (UseSeptumPlusStdHRS || LHRSAngle/deg>12.49)?"Q1 Entrance":((pMeanZvbtr<850.0)?"Sieve Slit":"SEPTUM Entrance"));
  ptInfo->AddText(Form("%s%s%s",line[0],line[1],line[2]));

  text=ptInfo->AddText(Form("Target: %.3f mm %s in cell %d",TargetL,strTg,SetupG2PTarget));
  text->SetTextColor(4);
  text=ptInfo->AddText(Form("#splitline{Helm Field: Angle=%.0f ^{o}, Ratio=%.2f}{Septum Field Ratio: L=%.4f, R=%.4f}",360.0-HelmRotAngle1/deg,HelmCurrentRatio,(UseSeptumField?SeptumCurrentRatioL:0.0),(UseSeptumField?SeptumCurrentRatioR:0.0)));
  text->SetTextColor(4);
  if(IsHRSTrack==2)
    {
      text=ptInfo->AddText(Form("#splitline{Beam=%.4f GeV}{%s transportation is used.}",Beam, ((UseSeptumPlusStdHRS)?"STD HRS":"HRS")));
      text->SetTextColor(4);
      if(LHRSMomentum<0.001)
	{
	  text=ptInfo->AddText("But the delta is set to ZERO");
	} 
      else
	{
	  text=ptInfo->AddText(Form("HRS P0: L=%.3f, R=%.3f, Angle=%.2f^{o}",
				    LHRSMomentum,RHRSMomentum,LHRSAngle/deg));
	}		
      text->SetTextColor(4);
    }
  else if(IsHRSTrack==1)
    {
      text=ptInfo->AddText(Form("The Acceptance Cut is defined as %s:",((cutshape==2)?"ellipse":"rectangle")));
      text->SetTextColor(4);
      text=ptInfo->AddText(Form("1) Vertical cut: %.3f<Thetavb_tr<%.3f",kVLow,kVHigh));
      text->SetTextColor(4);
      text=ptInfo->AddText(Form("2) Honrizontal cut: %.3f<Phivb_tr<%.3f",kHLow,kHHigh));
      text->SetTextColor(4);
      text=ptInfo->AddText(Form("3) Energy Cut: (P0-Pvb)/P0 < %.3f",kECut));
      text->SetTextColor(4);
    }
  else
    {
      text=ptInfo->AddText("The Acceptance Cut is defined as:");
      text->SetTextColor(4);
      text=ptInfo->AddText(Form(" Energy Cut: (P0-Pvb)/P0 < %.3f",kECut));
      text->SetTextColor(4);
    }
  ptInfo->Draw();
  gPad->Modified();

  pCan->cd(9);
  tree->Draw("P0>>hP0","","");
  hP0->SetTitle("P0: all(black) and cut(red) and Pvb, cut(blue); P0 and Pvb (GeV)");
  tree->Draw("P0>>hP0Cut",Cut,"same");
  hP0Cut->SetLineColor(2);
  tree->Draw("Pvb>>hPvbCut",Cut,"same");
  hPvbCut->SetLineColor(4);
  gPad->Modified();

  pCan->cd(10);
  tree->Draw("Theta0>>hTheta0","","");
  hTheta0->SetTitle("Theta0: all(black) and cut(red);Theta0 (rad)");
  tree->Draw("Theta0>>hTheta0Cut",Cut,"same");
  hTheta0Cut->SetLineColor(2);
  gPad->Modified();

  pCan->cd(11);
  tree->Draw("Phi0>>hPhi0","","");
  hPhi0->SetTitle("Phi0: all(black) and cut(red);Phi0 (rad)");
  tree->Draw("Phi0>>hPhi0Cut",Cut,"same");
  hPhi0Cut->SetLineColor(2);
  gPad->Modified();


  pCan->cd(7);
  gPad->SetRightMargin(0.15);
  double pMeanPvb,pMeanP0,pMeanElasXS;
  double pElasE=Beam/(1+2*Beam/TargetM*pow(sin(pMeanTheta0/2),2.0));
  
  tree->Draw("P0>>hP0Cut_w",Form("(%s)*ElasXS",Cut),"");
  pMeanP0=hP0Cut_w->GetMean();

	tree->Draw("Theta0>>hTheta0Cut_w",Form("(%s)*ElasXS",Cut),"");
	pMeanTheta0=hTheta0Cut_w->GetMean();
	pMeanElasXS=hTheta0Cut_w->Integral()/hTheta0Cut->GetEntries();
	tree->Draw("Pvb>>hPvb_w",Form("(%s)*ElasXS",Cut),"");
	hPvb_w->SetTitle("XS weighted Pvb (cut);Pvb (GeV)");
	pMeanPvb=hPvb_w->GetMean();
 

  if(pMeanP0/pElasE > 0.8)
    {
      tree->Draw("Pvb:Theta0>>h2PvbTheta0",Form("(%s)*ElasXS",Cut),"contz");
      h2PvbTheta0->SetTitle("Elastic XS weighted P_{vb} Vs #theta_{0} (cut);#theta_{0} (rad);P_{vb} (GeV) ");
    }
  else
    {
      tree->Draw("Pvb:Theta0>>h2PvbTheta0",Form("(%s)*XS",Cut),"contz");
      h2PvbTheta0->SetTitle("Inelas XS weighted P_{vb} Vs #theta{0} (cut);#theta_{0} (rad);P_{vb} (GeV) ");
    }
  //This method seems to give a wrong average
  //pMeanPvb=h2PvbTheta0->GetMean(2);
  //pMeanTheta0=h2PvbTheta0->GetMean(1);
  //TH1 *h1=(TH1*) (h2PvbTheta0->ProjectionX());
  //pMeanElasXS=h1->Integral()/h1->GetEntries();
  
  if(h2PvbTheta0->GetMaximum()>1000.0) gPad->SetLogz(1);
  h2PvbTheta0->Draw("contz");
  
  //show the stats
  TPaveStats *pPS = new TPaveStats(0.11,0.11,0.52,0.32,"brNDC");
  pPS->SetBorderSize(1);
  pPS->SetFillColor(0);
  pPS->SetTextAlign(12);
  pPS->SetOptStat(110);
  pPS->AddText(Form("<XS> = %.3f #mub",pMeanElasXS));
  pPS->AddText(Form("<P_{vb}> = %6.4f",pMeanPvb));
  pPS->AddText(Form("<#theta_{0}> = %6.2f^{o}",pMeanTheta0*57.3));
  pPS->Draw();
  gPad->Modified();


  pCan->cd(12);
  //tree->Draw("Pvb",Cut,"");
  tree->Draw("ElasXS>>hElasXS",XSCut,"");
  hElasXS->SetTitle("Elas XS (cut); XS (ub)");
  hElasXS->SetLineColor(2);

  TPaveText *ptRate = new TPaveText(0.30,0.4,0.99,0.90,"blNDC");
  ptRate->SetName("rate");
  ptRate->SetBorderSize(1);
  ptRate->SetFillColor(0);
  ptRate->SetTextAlign(12);

  double pHRSAcc=(LHRSAngle/deg<10.0)? 4.5.0:6.0;  //msr

  //Lumi=3.759E+33 * *I_na * density * length / molmass , (in g, cm)
  double pLumi=18.03*TargetL/2.54;  //100 mil C12
  if(TargetAtomicNumber<1.5) pLumi= 6.41*2*TargetL/2.54;  //100 mil H in CH2
  else if(TargetAtomicNumber<2.5) pLumi= 6.00*TargetL/44.0;  //4.4cm LHe
  else if(int(TargetAtomicNumber)==13) pLumi= 0.955*TargetL/0.254;  //10 mil Al


  TText *pTx=0;	
  double pRate=pLumi*pHRSAcc*pMeanElasXS;
  double pSieveRate=pRate/32.0;
  double pMinute=500000.0/pSieveRate/60.0;
	
  pTx=ptRate->AddText(Form("100 nA %.3f GeV, 100 mil %s",Beam,strTg));
  pTx->SetTextColor(2);
  pTx=ptRate->AddText(Form("HRS=%.2f deg, <XS>=%.2f ub",LHRSAngle/deg,pMeanElasXS));
  pTx->SetTextColor(2);
  pTx=ptRate->AddText(Form("Lumi=%.2fE+33/cm2/s, Acc=%.1f msr",pLumi,pHRSAcc));
  pTx->SetTextColor(2);
  pTx=ptRate->AddText(Form("Rate=%.3f kHz",pRate/1000.));
  pTx->SetTextColor(2);
  pTx=ptRate->AddText(Form("Sieve_Rate=%.2f Hz",pSieveRate));
  pTx->SetTextColor(2);
  pTx=ptRate->AddText(Form("500k_Sieve = %.1f minutes",pMinute));
  pTx->SetTextColor(2);
  ptRate->Draw();
  gPad->Modified();

  pCan->Update();
  char grName[255];
  if(trackindex>0)
    {
      if(bIsCombinedTree)
	{
	  sprintf(grName,"MonitorHisto_Track%d_Z%.0f_EAll_R%.1f.png",trackindex,
		  TargetAtomicNumber,HelmCurrentRatio);
	}
      else
	{
	  sprintf(grName,"MonitorHisto_Track%d_Z%.0f_E%.1f_R%.1f.png",trackindex,
		  TargetAtomicNumber,Beam,HelmCurrentRatio);
	}
    }
  else
    {
      if(bIsCombinedTree)
	{
	  sprintf(grName,"MonitorHisto_HRS%.1fdeg_Z%.0f_EAll_R%.1f.png",LHRSAngle/deg,
		  TargetAtomicNumber,HelmCurrentRatio);
	}
      else
	{
	  sprintf(grName,"MonitorHisto_HRS%.1fdeg_Z%.0f_E%.1f_R%.1f.png",LHRSAngle/deg,
		  TargetAtomicNumber,Beam,HelmCurrentRatio);
	}
    }
  pCan->SaveAs(grName);
  //print the summary
  char strLine0[255],strLine1[255];
  sprintf(strLine0,"%6s %6s %9s %8s %9s %7s %5s %6s %8s %11s %7s %8s %13s %12s\n",
	  "Target","Beam","BeamAngle","HRSAngle","HelmRatio","<Theta0>","<P0>","<Pvb>","<XS>(ub)",
	  "Lumi(10^33)","Acc(msr)","Rate(Hz)","SieveRate(Hz)","Minute(500k)");
  sprintf(strLine1,"%6s %6.3f %9.2f %8.2f %9.2f %7.2f %6.4f %6.4f %8.2f %11.2f %7.1f %8.2f %13.2f %12.3f\n",
	  strTg,Beam,BeamTiltedAngle/deg,LHRSAngle/deg,HelmCurrentRatio,
	  pMeanTheta0/deg,pMeanP0,pMeanPvb,pMeanElasXS,pLumi,pHRSAcc,pRate,pSieveRate,pMinute);
  cout<<strLine0<<strLine1;

  ofstream outfile;
  outfile.open (Form("result_track%d.txt",trackindex),ios_base::app);
  outfile<<strLine1;
  outfile.close();
}

