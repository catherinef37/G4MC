
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
Double_t        ThirdArmAngle;
Double_t        SetupThirdArmVD;
Double_t        ThirdArmRotZAngle;
Double_t        Pivot2ThirdArmFace;


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
   config->SetBranchAddress("ThirdArmAngle",&ThirdArmAngle);
   config->SetBranchAddress("SetupThirdArmVD",&SetupThirdArmVD);
   config->SetBranchAddress("ThirdArmRotZAngle",&ThirdArmRotZAngle);
   config->SetBranchAddress("Pivot2ThirdArmFace",&Pivot2ThirdArmFace);
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

void GetXS(char *filename="")
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

  track0->AddFriend("D");
  //TCut theCut="Pvb>0";
  TCut theCut="TA1_Edep>7.5 && TA2_Edep>39";
 

  TCanvas *c1=new TCanvas("c1","",1000,800);
  c1->Divide(2,2);

  c1->cd(1);
  track0->Draw("XS:P0>>h2xsp0(100,0,1.0,200,0,0.2)","","prof");
  h2xsp0->SetTitle("Averaged inelastic cross section; P0(GeV) ; XS (#mub) ");
  TH1 *hxsprof=(TH1*)(gROOT->FindObject("h2xsp0"));

  c1->cd(2);
  track0->Draw("P0>>hp0(100,0,1.0)","");
  hp0->SetTitle("Thrown (black) and detected(red); P0 or Pvb (GeV) ");
  TH1 *hp0=(TH1*)(gROOT->FindObject("hp0"));
  track0->Draw("P0>>hp0det(100,0,1.0)",theCut,"same");
  TH1 *hp0det=(TH1*)(gROOT->FindObject("hp0det"));
  hp0det->SetLineColor(2);
  TH1 *heff = hp0det->Clone("heff");

  c1->cd(3);
  heff->Divide(hp0);
  heff->Draw();

  /*
  pt1= new TPaveText(0.11,0.74,0.45,0.895,"brNDC");
  pt1->SetBorderSize(0);
  pt1->SetFillColor(0);
  pt1->AddText("detector: 20cmx60cm,210cm,70^{o}");
  pt1->Draw("same");
  */

  c1->cd(4);
  TH1 *hxs=heff->Clone("hxs");
  //hxs->Sumw2();
  hxs->Multiply(hxsprof);
  hxs->SetTitle("Detected inelastic XS; P0 (GeV) ; inelastic XS (#mub)");
  double pErr[1000]={0,...,0};
  hxs->SetError(pErr);
  //hxs->Sumw2();
  hxs->SetMarkerColor(2);
  hxs->SetMarkerStyle(2);
  hxs->Draw("prof");

  c1->cd(1);
  //hxsprof->Draw();
  hxs->Draw("profsame");

  c1->SaveAs(Form("InelasXS_E%.3f.png",Beam));
}

