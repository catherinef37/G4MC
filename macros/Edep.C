
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

void Edep(char *filename="")
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


  D->AddFriend("track0");
  TCanvas *c1=new TCanvas("c1","",900,700);
  c1->Divide(2,2);
  //c1->cd(1);
  //config->Draw("Beam");
  //double Beam=htemp->GetMean();

  c1->cd(1);
  D->Draw("TA1_Pid%10000 >> hPid1","T_Tid>1 && T_X>1500","");
  hPid1->SetTitle("particles in dE plane; TA1_Pid");

  c1->cd(2);
  D->Draw("TA2_Pid%10000 >> hPid2","T_Tid>1 && T_X>1500","");
  hPid2->SetTitle("particles in E plane; TA2_Pid");

  c1->cd(3);
  D->Draw("TA1_Edep >> hEin","T_Tid>1 && T_X>1500","");
  double pEin=hEin->GetMean()*hEin->GetEntries();
  hEin->SetTitle(Form("dE plane: %.2f MeV/160ns (100k e^{-} @ 100nA); TA1_Edep (MeV)",pEin));

  c1->cd(4);
  D->Draw("TA2_Edep >> hEout","T_Tid>1 && T_X>1500","");
  double pEout=hEout->GetMean()*hEout->GetEntries();
  hEout->SetTitle(Form("E plane: %.2f MeV/160ns (100k e^{-} @ 100nA); TA2_Edep (MeV)",pEout));

  c1->SaveAs(Form("n_Edep_E%.3f.png",Beam));
}

