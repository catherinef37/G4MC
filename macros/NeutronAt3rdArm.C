
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
Double_t        SeptumCurrentRatio;
Double_t        BigBiteAngle;
Double_t        BigBiteTiltAngle;
Double_t        Pivot2BigBiteFace;

bool		bIsCombinedTree;  //to ideneify if this is a combined ntuple


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

void NeutronAt3rdArm(char *filename="")
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

  TCanvas *c1=new TCanvas("c1","",600,800);
  c1->Divide(1,2);
  int nEl=100000,nN=0,nNdet=0;

  //nN = track0->GetEntries();
  c1->cd(1);
  track0->Draw("Z0>>hZ0(1500,-1400,100)","");
  nN=hZ0->GetEntries();
  if(nN>10) gPad->SetLogy(1);
  hZ0->SetTitle("vertex Z, all(black) and detected(red); Z0(mm) ");
  track0->Draw("Z0>>hZ0_det(1500,-1400,100)","abs(Xvb-2000)<1000","same");
  hZ0_det->SetLineColor(2);
  nNdet=hZ0_det->GetEntries();
  pt1= new TPaveText(0.11,0.74,0.45,0.895,"brNDC");
  pt1->SetBorderSize(0);
  pt1->SetFillColor(0);
  pt1->AddText("detector: 20cmx60cm,210cm,70^{o}");
  pt1->AddText(Form("all: %6d",nN));
  pt1->AddText(Form("#color[2]{detected: %6d}",nNdet));
  pt1->Draw("same");
  
  c1->cd(2);
  track0->Draw("P0>>hP0(200,0,0.4)");
  if(nN>10) gPad->SetLogy(1);
  hP0->SetTitle("P0(black) and Pvb(red); P0 or Pvb (GeV) ");
  track0->Draw("Pvb>>hPvb","Pvb>0 && abs(Xvb-2000)<1000","same");
  hPvb->SetLineColor(2);

  c1->SaveAs(Form("n_E%.3f.png",Beam));
}

