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

//Declaration of leaves types for config tree
Int_t           Run;
Int_t           SkimLevel;
Int_t           BookTrees;
Double_t        HelmCurrentRatio;
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

bool			bIsCombinedTree;
bool ReadConfig()
{
	TTree *config = (TTree*)gDirectory->Get("config");

	// Set branch addresses.// Set branch addresses.
	config->SetBranchAddress("Run",&Run);
	config->SetBranchAddress("SkimLevel",&SkimLevel);
	config->SetBranchAddress("BookTrees",&BookTrees);
	config->SetBranchAddress("HelmCurrentRatio",&HelmCurrentRatio);
	config->SetBranchAddress("Beam",&Beam);
	config->SetBranchAddress("TargetM",&TargetM);
	config->SetBranchAddress("TargetAtomicNumber",&TargetAtomicNumber);
	config->SetBranchAddress("LHRSAngle",&LHRSAngle);
	config->SetBranchAddress("RHRSAngle",&RHRSAngle);
	config->SetBranchAddress("TargetXOffset",&TargetXOffset);
	config->SetBranchAddress("TargetYOffset",&TargetYOffset);
	config->SetBranchAddress("TargetZOffset",&TargetZOffset);
	config->SetBranchAddress("PivotXOffset",&PivotXOffset);
	config->SetBranchAddress("PivotYOffset",&PivotYOffset);
	config->SetBranchAddress("PivotZOffset",&PivotZOffset);
	config->SetBranchAddress("LHRSMomentum",&LHRSMomentum);
	config->SetBranchAddress("RHRSMomentum",&RHRSMomentum);

	// This is the loop skeleton
	// To read only selected branches, Insert statements like:
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

double fcn(double *x, double *par)
{
	return par[0]/x[0]+par[1]+par[2]*x[0];
}

void BFieldRatio_Theta0(double xx=3.0)
{
	gStyle->SetStatW(0.16);
	gStyle->SetStatH(0.06);
	gStyle->SetStatX(0.90);
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);

	//ReadConfig();

	TChain *T0=new TChain("track0");
	T0->Add("scan_nt_6_00.root");

	TChain *T1=new TChain("track0");
	T1->Add("nt_scan6_fieldratio0.5.root");

	TChain *T2=new TChain("track0");
	T2->Add("nt_scan6_fieldratio2.0.root");


	TCanvas *pCan = new TCanvas("pCan","BFieldRatio",900,720);
	pCan->Divide(2,2);
	char Cut[]="TrackClass>3 && abs(P0-Pvb)/P0<0.05";
	TF1 *f0=new TF1("f0","[0]/x+[1]+[2]*x",0.3,3.5);
	TF1 *f1=new TF1("f1","[0]/x+[1]+[2]*x",0.3,3.5);
	TF1 *f2=new TF1("f2","[0]/x+[1]+[2]*x",0.3,3.5);
	f0->SetLineColor(1);
	f1->SetLineColor(2);
	f2->SetLineColor(4);
	
	pCan->cd(1);
	T0->Draw("Theta0:P0>>h2TP_0(70,0.,3.5,80,0,0.8)",Cut,"Prof");
	h2TP_0->SetTitle("Scale = 1.0 ;P_{0} {GeV/c} ; #theta_{0} (rad) ");
	h2TP_0->SetMarkerColor(1);
	h2TP_0->Fit(f0,"R","",0.3,3.4);
	
	pCan->cd(2);
	T1->Draw("Theta0:P0>>h2TP_1(70,0.,3.5,80,0,0.8)",Cut,"Prof");
	h2TP_1->SetTitle("Scale = 0.5 ;P_{0} {GeV/c} ; #theta_{0} (rad) ");
	h2TP_1->SetMarkerColor(2);
	h2TP_1->Fit(f1,"R","",0.3,3.4);
	
	pCan->cd(3);
	T2->Draw("Theta0:P0>>h2TP_2(70,0.,3.5,80,0,0.8)",Cut,"Prof");
	h2TP_2->SetTitle("Scale = 2.0 ;P_{0} {GeV/c} ; #theta_{0} (rad) ");
	h2TP_2->SetMarkerColor(4);
	h2TP_2->Fit(f2,"R","",0.5,3.4);
	pCan->cd(4);
	//gStyle->SetOptFit(0);
	TH2D *H2Frame=new TH2D("Frame","Blue(2.0), Black(1.0), Red(0.5) ;P_{0} {GeV/c} ; #theta_{0} (rad)",
		70,0.,3.5,80,0,0.8);
	H2Frame->Draw();
	h2TP_0->Draw("profsame");
	h2TP_1->Draw("profsame");
	h2TP_2->Draw("profsame");
	
	//double xx=1.6;
	double yy = fcn(&xx,f2->GetParameters());
	TLine *LH=new TLine(0.,yy,3.5,yy);
	TLine *LV1=new TLine(xx/4.,0.6,xx/4.,0.);
	TLine *LV2=new TLine(xx/2.,0.6,xx/2.,0.);
	TLine *LV3=new TLine(xx,0.6,xx,0.);
	LH->Draw("same");
	LV1->Draw("same");
	LV2->Draw("same");
	LV3->Draw("same");
	
	pCan->SaveAs("Theta0VsP0_BFieldRatio_HRS6.png");
}


void BFieldRatio_Theta0_tr(double xx=3.0)
{
	gStyle->SetStatW(0.16);
	gStyle->SetStatH(0.06);
	gStyle->SetStatX(0.40);
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);

	//ReadConfig();

	TChain *T0=new TChain("track0");
	T0->Add("scan_nt_6_00.root");

	TChain *T1=new TChain("track0");
	T1->Add("nt_scan6_fieldratio0.5.root");

	TChain *T2=new TChain("track0");
	T2->Add("nt_scan6_fieldratio2.0.root");


	TCanvas *pCan = new TCanvas("pCan","BFieldRatio",900,720);
	pCan->Divide(2,2);
	char Cut[]="TrackClass>3 && abs(P0-Pvb)/P0<0.05";
	TF1 *f0=new TF1("f0","[0]/x+[1]+[2]*x",0.3,3.5);
	TF1 *f1=new TF1("f1","[0]/x+[1]+[2]*x",0.3,3.5);
	TF1 *f2=new TF1("f2","[0]/x+[1]+[2]*x",0.3,3.5);
	f0->SetLineColor(1);
	f1->SetLineColor(2);
	f2->SetLineColor(4);
	
	pCan->cd(1);
	T0->Draw("Theta0_tr:P0>>h2TtrP_0(70,0.,3.5,50,-0.5,0)",Cut,"Prof");
	h2TtrP_0->SetTitle("Scale = 1.0 ;P_{0} {GeV/c} ; #theta_{0}^{tr} (rad) ");
	h2TtrP_0->SetMarkerColor(1);
	h2TtrP_0->Fit(f0,"R","",0.3,3.4);
	
	pCan->cd(2);
	T1->Draw("Theta0_tr:P0>>h2TtrP_1(70,0.,3.5,50,-0.5,0)",Cut,"Prof");
	h2TtrP_1->SetTitle("Scale = 0.5 ;P_{0} {GeV/c} ; #theta_{0}^{tr} (rad) ");
	h2TtrP_1->SetMarkerColor(2);
	h2TtrP_1->Fit(f1,"R","",0.3,3.4);
	
	pCan->cd(3);
	T2->Draw("Theta0_tr:P0>>h2TtrP_2(70,0.,3.5,50,-0.5,0)",Cut,"Prof");
	h2TtrP_2->SetTitle("Scale = 2.0 ;P_{0} {GeV/c} ; #theta_{0}^{tr} (rad) ");
	h2TtrP_2->SetMarkerColor(4);
	h2TtrP_2->Fit(f2,"R","",0.5,3.4);
	pCan->cd(4);
	//gStyle->SetOptFit(0);
	TH2D *H2Frame_tr=new TH2D("Frame_tr","Blue(2.0), Black(1.0), Red(0.5) ;P_{0} {GeV/c} ; #theta_{0}^{tr} (rad)",
		70,0.,3.5,50,-0.4,0);
	Frame_tr->Draw();
	h2TtrP_0->Draw("profsame");
	h2TtrP_1->Draw("profsame");
	h2TtrP_2->Draw("profsame");
	
	//double xx=1.6;
	double yy = fcn(&xx,f2->GetParameters());
	TLine *LH=new TLine(0.,yy,3.5,yy);
	TLine *LV1=new TLine(xx/4.,-0.4,xx/4.,0.);
	TLine *LV2=new TLine(xx/2.,-0.4,xx/2.,0.);
	TLine *LV3=new TLine(xx,-0.4,xx,0.);
	LH->Draw("same");
	LV1->Draw("same");
	LV2->Draw("same");
	LV3->Draw("same");
	
	pCan->SaveAs("Theta0_trVsP0_BFieldRatio_HRS6.png");
}

void BFieldRatio(double xx=3.0)
{
	BFieldRatio_Theta0_tr(xx);
	BFieldRatio_Theta0(xx);
}
