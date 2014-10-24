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

//cutshape: 1 rectangle; 2 ellipse
void DrawCoil(char *filename="", int trackindex=0)
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

	TTree *tree = (TTree*)gDirectory->Get(Form("track%d",trackindex));
	if(trackindex>0) tree->AddFriend("track0");

	TCanvas *pCan = new TCanvas("pCan","ACC Due to Coils",900,600);
	pCan->Divide(3,2,0.001,0.001);

	char Cut[1000],CoilCut[1000];

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
		sprintf(Cut,"TrackClass>1");

		tree->Draw("ElasXS>>hElasXS_T2",Cut);
		if(hElasXS_T2->GetEntries()>1) 
		{
			MaxXS=hElasXS_T1->GetMean()*40.0;
			IsHRSTrack=1;
		}
		else
		{			
			sprintf(Cut,"TrackClass>=0");
		}
	}

	sprintf(CoilCut,"StepDsty>7 && (%s)",Cut);

	////////////////////////////////////////////////////////////////////////

	pCan->cd(1);
	tree->Draw("StepX:StepZ>>h2StepXZCoil",CoilCut,"");
	h2StepXZCoil->SetTitle("StepX:StepZ @ coils; StepZ (mm); StepX (mm)");
	h2StepXZCoil->SetMarkerStyle(4);
	h2StepXZCoil->SetMarkerColor(4);
	h2StepXZCoil->Draw();

	pCan->cd(2);
	tree->Draw("P0>>hP0",Cut,"");
	hP0->SetTitle("P0: all(black) and hit coils (blue); P0 (GeV/c)");
	tree->Draw("P0>>hP0Coil",CoilCut,"same");
	hP0Coil->SetLineColor(4);

	//show the stats
	TPaveStats *pPS = new TPaveStats(0.65,0.75,0.99,0.90,"brNDC");
	pPS->SetName("mywords");
	pPS->SetBorderSize(2);
	pPS->SetFillColor(0);
	pPS->SetTextAlign(12);
	pPS->SetOptStat(110);

	double pNAll=hP0->GetEntries();
	double pNCoil=hP0Coil->GetEntries();
	pPS->AddText(Form("All = %.0f",pNAll));
	pPS->AddText(Form("Coil = %.0f",pNCoil));
	//hP0->GetListOfFunctions()->Add(pPS);
	//pPS->SetParent(hP0->GetListOfFunctions());
	pPS->Draw();

	pCan->cd(3);
	tree->Draw("Pvb>>hPvbCoil",CoilCut,"");
	hPvbCoil->SetTitle("Pvb, hit coils; Pvb (GeV/c)");
	hPvbCoil->SetLineColor(4);
	//show the stats
	TPaveStats *pInfo = new TPaveStats(0.65,0.75,0.99,0.90,"brNDC");
	pInfo->SetName("mywords");
	pInfo->SetBorderSize(2);
	pInfo->SetFillColor(0);
	pInfo->SetTextAlign(12);
	pInfo->SetOptStat(110);
	pInfo->AddText(Form("HRS = %.1f",LHRSAngle*180./3.14159));
	pInfo->AddText(Form("E = %.3f",Beam));
	//hP0->GetListOfFunctions()->Add(pInfo);
	//pInfo->SetParent(hPvbCoil->GetListOfFunctions());
	pInfo->Draw();

	pCan->cd(4);
	tree->Draw("Theta0*57.3:Phi0*57.3>>h2T0P0",Cut,"");
	h2T0P0->SetTitle("Theta0:Phi0 all(black) Hit coils (blue); Phi0 (deg);Theta0 (deg)");
	tree->Draw("Theta0*57.3:Phi0*57.3>>h2T0P0Coil",CoilCut,"*same");
	h2T0P0Coil->SetMarkerStyle(4);
	h2T0P0Coil->SetMarkerColor(4);
	h2T0P0->Draw();
	h2T0P0Coil->Draw("same");


	pCan->cd(5);
	tree->Draw("X0:Y0>>h2X0Y0",Cut,"");
	h2X0Y0->SetTitle("X0:Y0 all(black) Hit coils (blue); Y0 (mm);X0 (mm)");
	tree->Draw("X0:Y0>>h2X0Y0Coil",CoilCut,"same");
	h2X0Y0Coil->SetMarkerStyle(4);
	h2X0Y0Coil->SetMarkerColor(4);
	h2X0Y0->Draw();
	h2X0Y0Coil->Draw("same");

	pCan->cd(6);
	tree->Draw("Z0>>hZ0",Cut,"");
	hZ0->SetTitle("Z0: all(black) and hit coils (blue); Z0 (mm)");
	tree->Draw("Z0>>hZ0Coil",CoilCut,"same");
	hZ0Coil->SetLineColor(4);


	pCan->cd(0);
	pCan->Update();


	char grName[255];

	if(bIsCombinedTree)
	{
		sprintf(grName,"Coil_Track%d_HRS%.1fdeg_EAll.gif",trackindex,
			LHRSAngle*180./3.14159);
	}
	else
	{
		sprintf(grName,"Coil_Track%d_HRS%.1fdeg_E%.1f.gif",trackindex,
			LHRSAngle*180./3.14159,Beam);
	}


	pCan->SaveAs(grName);
	//pCan->SaveAs("HRS12.5_HitCoils.png");
	//pCan->SaveAs("HRS12.5_HitCoils.C");

}
