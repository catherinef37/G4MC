
//Declaration of leaves types tree
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
Double_t        HelmCurrentRatio=1.0;
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
Double_t        SeptumCurrentRatio;
Double_t        BigBiteAngle;
Double_t        BigBiteTiltAngle;
Double_t        Pivot2BigBiteFace;

bool			bIsCombinedTree;	//to ideneify if this is a combined ntuple



bool ReadConfig()
{
  TTree *config = (TTree*)gDirectory->Get("config");

  // Set branch addresses.
  config->SetBranchAddress("Run",&Run);
  config->SetBranchAddress("SkimLevel",&SkimLevel);
  config->SetBranchAddress("BookTrees",&BookTrees);
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
  if(config->GetBranch("HelmCurrentRatio"))
    {
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
      config->SetBranchAddress("SeptumXOffset",&SeptumXOffset);
      config->SetBranchAddress("SeptumYOffset",&SeptumYOffset);
      config->SetBranchAddress("SeptumZOffset",&SeptumZOffset);
      config->SetBranchAddress("SeptumRotAxis1",&SeptumRotAxis1);
      config->SetBranchAddress("SeptumRotAxis2",&SeptumRotAxis2);
      config->SetBranchAddress("SeptumRotAxis3",&SeptumRotAxis3);
      config->SetBranchAddress("SeptumRotAngle1",&SeptumRotAngle1);
      config->SetBranchAddress("SeptumRotAngle2",&SeptumRotAngle2);
      config->SetBranchAddress("SeptumRotAngle3",&SeptumRotAngle3);
      config->SetBranchAddress("SeptumCurrentRatio",&SeptumCurrentRatio);
      config->SetBranchAddress("BigBiteAngle",&BigBiteAngle);
      config->SetBranchAddress("BigBiteTiltAngle",&BigBiteTiltAngle);
      config->SetBranchAddress("Pivot2BigBiteFace",&Pivot2BigBiteFace);
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

void Phi0VsPvb(char *filename="")
{
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetOptFit(1);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetStatH(0.2);
  gStyle->SetStatW(0.17);
  gStyle->SetStatX(1.0-gStyle->GetPadRightMargin());
	

  if(strlen(filename)>5)
    {
      cout<<"trying to open root file "<<filename<<endl;
      TFile* _file0=TFile::Open(filename);
    }
  ReadConfig();

  TF1 *f1 = new TF1("f1","[0]/x+[1]+[2]*x+[3]*x*x",0.2,3.3);
  f1->SetFillColor(0);
  f1->SetMarkerStyle(20);
  f1->SetMarkerColor(2);
  f1->SetLineColor(2);
  f1->SetLineWidth(3);
	
  TF1 *f2 = new TF1("f2","[0]/x+[1]+[2]*x+[3]*x*x",0.2,3.3);
  f2->SetFillColor(0);
  f2->SetMarkerStyle(20);
  f2->SetMarkerColor(4);
  f2->SetLineColor(4);
  f2->SetLineWidth(3);

  TText *text=0;
  char HRSCut[255];
  sprintf(HRSCut,"abs((P0-Pvb)/P0)<0.05 && TrackClass>3");

  TCanvas *pCan=new TCanvas("pCan","Phi0 in G2P",600,800);
  pCan->Divide(1,2);
  pCan->cd(1);
  gPad->SetRightMargin(0.05);

  track0->Draw("Phi0:Pvb>>h2Phi0Pvb(175,0.,3.50,100,0,1.0)",HRSCut,"prof");
  h2Phi0Pvb->SetMarkerStyle(20);
  h2Phi0Pvb->SetMarkerColor(1);
  h2Phi0Pvb->SetTitle(Form("#phi_{0} Vs P_{vb} for HRS=%.1f^{o}; P_{vb} (GeV); #phi_{0} (rad)",LHRSAngle*180./3.14159));
  h2Phi0Pvb->GetYaxis()->SetRangeUser(0,1.5);
  h2Phi0Pvb->Draw("prof");
  h2Phi0Pvb->Fit(f1,"R","", 0.3,3.5);

  TPaveText *pt1 = new TPaveText(1.5,0.65,3.3,0.80,"br");
  pt1->SetBorderSize(0);
  pt1->SetFillColor(0);
  pt1->AddText("#color[2]{f(x)=p0/x+p1+p2*x+p3*x*x}");
  pt1->Draw();

  pCan->cd(2);
  gPad->SetRightMargin(0.05);
	
  track0->Draw("Phi0:P0>>h2Phi0P0(175,0.,3.50,100,0,1.0)",HRSCut,"prof");
  h2Phi0P0->SetMarkerStyle(20);
  h2Phi0P0->SetMarkerColor(1);
  h2Phi0P0->SetTitle(Form("#phi_{0} Vs P_{0} for HRS=%.1f^{o}; P_{0} (GeV); #phi_{0} (rad)",LHRSAngle*180./3.14159));
  h2Phi0P0->GetYaxis()->SetRangeUser(0,1.5);
  h2Phi0P0->Draw("prof");
  h2Phi0P0->Fit(f2,"R","", 0.3,3.5);
  f1->Draw("same");
	
  TPaveText *pt2 = new TPaveText(1.5,0.65,3.3,0.80,"br");
  pt2->SetBorderSize(0);
  pt2->SetFillColor(0);
  text=pt2->AddText("f(x)=p0/x+p1+p2*x+p3*x*x");
  text->SetTextColor(4);
  pt2->Draw();

  pCan->cd();
  pCan->SaveAs(Form("PhiVsP_HRS%.1f.gif",LHRSAngle*180./3.14159));	
}
 
void Theta0VsPvb(char *filename="")
{
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetOptFit(1);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetStatH(0.2);
  gStyle->SetStatW(0.17);
  gStyle->SetStatX(1.0-gStyle->GetPadRightMargin());
	

  if(strlen(filename)>5)
    {
      cout<<"trying to open root file "<<filename<<endl;
      TFile* _file0=TFile::Open(filename);
    }
  ReadConfig();

  TF1 *f1 = new TF1("f1","exp([0]+[1]*x)+[2]+[3]*x",0.2,3.3);
  f1->SetFillColor(0);
  f1->SetMarkerStyle(20);
  f1->SetMarkerColor(2);
  f1->SetLineColor(2);
  f1->SetLineWidth(3);
	
  TF1 *f2 = new TF1("f2","exp([0]+[1]*x)+[2]+[3]*x",0.2,3.3);
  f2->SetFillColor(0);
  f2->SetMarkerStyle(20);
  f2->SetMarkerColor(4);
  f2->SetLineColor(4);
  f2->SetLineWidth(3);

  TText *text=0;
  char HRSCut[255];
  sprintf(HRSCut,"abs((P0-Pvb)/P0)<0.05 && TrackClass>3");

  TCanvas *pCan=new TCanvas("pCan","Theta0 in G2P",600,800);
  pCan->Divide(1,2);
  pCan->cd(1);
  gPad->SetRightMargin(0.05);

  track0->Draw("Theta0:Pvb>>h2Theta0Pvb(350,0.,3.50,100,0,1.0)",HRSCut,"prof");
  h2Theta0Pvb->SetMarkerStyle(20);
  h2Theta0Pvb->SetMarkerColor(1);
  h2Theta0Pvb->SetTitle(Form("#theta_{0} Vs P_{vb} for HRS=%.1f^{o}; P_{vb} (GeV); #theta_{0} (rad)",LHRSAngle*180./3.14159));
  h2Theta0Pvb->GetYaxis()->SetRangeUser(0,0.8);
  h2Theta0Pvb->Draw("prof");
  h2Theta0Pvb->Fit(f1,"R","", 0.3,3.3);

  TPaveText *pt1 = new TPaveText(1.5,0.35,3.3,0.45,"br");
  pt1->SetBorderSize(0);
  pt1->SetFillColor(0);
  pt1->AddText("#color[2]{f(x)=exp(p0+p1*x)+p2+p3*x}");
  pt1->Draw();

  pCan->cd(2);
  gPad->SetRightMargin(0.05);
	
  track0->Draw("Theta0:P0>>h2Theta0P0(350,0.,3.50,100,0,1.0)",HRSCut,"prof");
  h2Theta0P0->SetMarkerStyle(20);
  h2Theta0P0->SetMarkerColor(1);
  h2Theta0P0->SetTitle(Form("#theta_{0} Vs P_{0} for HRS=%.1f^{o}; P_{0} (GeV); #theta_{0} (rad)",LHRSAngle*180./3.14159));
  h2Theta0P0->GetYaxis()->SetRangeUser(0,0.8);
  h2Theta0P0->Draw("prof");
  h2Theta0P0->Fit(f2,"R","", 0.3,3.3);
  f1->Draw("same");
	
  TPaveText *pt2 = new TPaveText(1.5,0.35,3.3,0.45,"br");
  pt2->SetBorderSize(0);
  pt2->SetFillColor(0);
  text=pt2->AddText("f(x)=exp(p0+p1*x)+p2+p3*x");
  text->SetTextColor(4);
  pt2->Draw();

  pCan->cd();
  pCan->SaveAs(Form("ThetaVsP_HRS%.1f.gif",LHRSAngle*180./3.14159));	
}
 
void Theta0VsPvb_newfcn(char *filename="")
{
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetOptFit(1);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetStatH(0.2);
  gStyle->SetStatW(0.17);
  gStyle->SetStatX(1.0-gStyle->GetPadRightMargin());
	

  if(strlen(filename)>5)
    {
      cout<<"trying to open root file "<<filename<<endl;
      TFile* _file0=TFile::Open(filename);
    }
  ReadConfig();

  TF1 *f1 = new TF1("f1","[0]/x+[1]+[2]*x",0.2,3.3);
  f1->SetFillColor(0);
  f1->SetMarkerStyle(20);
  f1->SetMarkerColor(2);
  f1->SetLineColor(2);
  f1->SetLineWidth(3);
	
  TF1 *f2 = new TF1("f2","[0]/x+[1]+[2]*x",0.2,3.3);
  f2->SetFillColor(0);
  f2->SetMarkerStyle(20);
  f2->SetMarkerColor(4);
  f2->SetLineColor(4);
  f2->SetLineWidth(3);

  TText *text=0;
  char HRSCut[255];
  sprintf(HRSCut,"abs((P0-Pvb)/P0)<0.05 && TrackClass>3");

  TCanvas *pCan=new TCanvas("pCan","Theta0 in G2P",600,800);
  pCan->Divide(1,2);
  pCan->cd(1);
  gPad->SetRightMargin(0.05);

  track0->Draw("Theta0:Pvb>>h2Theta0Pvb(175,0.,3.50,100,0,1.0)",HRSCut,"prof");
  h2Theta0Pvb->SetMarkerStyle(20);
  h2Theta0Pvb->SetMarkerColor(1);
  h2Theta0Pvb->SetTitle(Form("#theta_{0} Vs P_{vb} for HRS=%.1f^{o}; P_{vb} (GeV); #theta_{0} (rad)",LHRSAngle*180./3.14159));
  h2Theta0Pvb->GetYaxis()->SetRangeUser(0,0.8);
  h2Theta0Pvb->Draw("prof");
  h2Theta0Pvb->Fit(f1,"R","", 0.3,3.3);

  TPaveText *pt1 = new TPaveText(1.5,0.35,3.3,0.45,"br");
  pt1->SetBorderSize(0);
  pt1->SetFillColor(0);
  pt1->AddText("#color[2]{f(x)=p0/x+p1+p2*x}");
  pt1->Draw();

  pCan->cd(2);
  gPad->SetRightMargin(0.05);
	
  track0->Draw("Theta0:P0>>h2Theta0P0(175,0.,3.50,100,0,1.0)",HRSCut,"prof");
  h2Theta0P0->SetMarkerStyle(20);
  h2Theta0P0->SetMarkerColor(1);
  h2Theta0P0->SetTitle(Form("#theta_{0} Vs P_{0} for HRS=%.1f^{o}; P_{0} (GeV); #theta_{0} (rad)",LHRSAngle*180./3.14159));
  h2Theta0P0->GetYaxis()->SetRangeUser(0,0.8);
  h2Theta0P0->Draw("prof");
  h2Theta0P0->Fit(f2,"R","", 0.3,3.3);
  f1->Draw("same");
	
  TPaveText *pt2 = new TPaveText(1.5,0.35,3.3,0.45,"br");
  pt2->SetBorderSize(0);
  pt2->SetFillColor(0);
  text=pt2->AddText("f(x)=p0/x+p1+p2*x");
  text->SetTextColor(4);
  pt2->Draw();

  pCan->cd();
  pCan->SaveAs(Form("ThetaVsP_HRS%.1f_newfcn.gif",LHRSAngle*180./3.14159));	
}
////////////////////////////////////////////////////////////////
 
void Theta0_trVsPvb(char *filename="")
{
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetOptFit(1);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetStatH(0.2);
  gStyle->SetStatW(0.17);
  //gStyle->SetStatStyle(4000);
  gStyle->SetStatX(0.95);
  gStyle->SetStatY(0.65);
	

  if(strlen(filename)>5)
    {
      cout<<"trying to open root file "<<filename<<endl;
      TFile* _file0=TFile::Open(filename);
    }
  ReadConfig();

  TF1 *f1 = new TF1("f1","[0]/x+[1]+[2]*x",0.2,3.3);
  f1->SetFillColor(0);
  f1->SetMarkerStyle(20);
  f1->SetMarkerColor(2);
  f1->SetLineColor(2);
  f1->SetLineWidth(3);
	
  TF1 *f2 = new TF1("f2","[0]/x+[1]+[2]*x",0.2,3.3);
  f2->SetFillColor(0);
  f2->SetMarkerStyle(20);
  f2->SetMarkerColor(4);
  f2->SetLineColor(4);
  f2->SetLineWidth(3);

  TText *text=0;
  char HRSCut[255];
  sprintf(HRSCut,"abs((P0-Pvb)/P0)<0.05 && TrackClass>3");

  TCanvas *pCan=new TCanvas("pCan","Theta0_tr in G2P",600,800);
  pCan->Divide(1,2);
  pCan->cd(1);
  gPad->SetRightMargin(0.05);

  track0->Draw("Theta0_tr:Pvb>>h2Theta0_trPvb(175,0.,3.50,100,0,1.0)",HRSCut,"prof");
  h2Theta0_trPvb->SetMarkerStyle(20);
  h2Theta0_trPvb->SetMarkerColor(1);
  h2Theta0_trPvb->SetTitle(Form("#theta_{0}^{tr} Vs P_{vb} for HRS=%.1f^{o}; P_{vb} (GeV); #theta_{0}^{tr} (rad)",LHRSAngle*180./3.14159));
  h2Theta0_trPvb->GetYaxis()->SetRangeUser(-0.5,0.1);
  h2Theta0_trPvb->Draw("prof");
  h2Theta0_trPvb->Fit(f1,"R","", 0.3,3.3);

  TPaveText *pt1 = new TPaveText(1.5,-0.35,3.3,-0.25,"br");
  pt1->SetBorderSize(0);
  pt1->SetFillColor(0);
  pt1->AddText("#color[2]{f(x)=p0/x+p1+p2*x}");
  pt1->Draw();

  pCan->cd(2);
  gPad->SetRightMargin(0.05);
	
  track0->Draw("Theta0_tr:P0>>h2Theta0_trP0(175,0.,3.50,100,0,1.0)",HRSCut,"prof");
  h2Theta0_trP0->SetMarkerStyle(20);
  h2Theta0_trP0->SetMarkerColor(1);
  h2Theta0_trP0->SetTitle(Form("#theta_{0}^{tr} Vs P_{0} for HRS=%.1f^{o}; P_{0} (GeV); #theta_{0}^{tr} (rad)",LHRSAngle*180./3.14159));
  h2Theta0_trP0->GetYaxis()->SetRangeUser(-0.5,0.1);
  h2Theta0_trP0->Draw("prof");
  h2Theta0_trP0->Fit(f2,"R","", 0.3,3.3);
  f1->Draw("same");
	
  TPaveText *pt2 = new TPaveText(1.5,-0.35,3.3,-0.25,"br");
  pt2->SetBorderSize(0);
  pt2->SetFillColor(0);
  text=pt2->AddText("f(x)=p0/x+p1+p2*x");
  text->SetTextColor(4);
  pt2->Draw();

  pCan->cd();
  pCan->SaveAs(Form("Theta0_trVsP_HRS%.1f.gif",LHRSAngle*180./3.14159));	
}

void ThetaVsPvb(const char *filename="")
{
  //plot them all
  Theta0VsPvb(filename);
  Theta0VsPvb_newfcn(filename);
  Phi0VsPvb(filename);
  Theta0_trVsPvb(filename);
}
 

