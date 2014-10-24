
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

void GetThirdArmEff(char *filename="")
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

  double pTA1_PMin=0.27,pTA1_PMax=0.41;
  pTA1_PMin=0.1,pTA1_PMax=0.6;
  double dECut1=8.5,dECut2=20.0;
  if(Beam>1.7)
    {
      pTA1_PMin=0.32,pTA1_PMax=0.46;
      pTA1_PMin=0.1,pTA1_PMax=0.6;
      dECut1=7.5,dECut2=39.0;
    }
  else if(Beam>2.2)
    {
      pTA1_PMin=0.32,pTA1_PMax=0.46;
      pTA1_PMin=0.1,pTA1_PMax=0.6;
      dECut1=7.5,dECut2=39.0;
    }

  TCut TA1PCut=Form("TA1_P>%.3f && TA1_P<%.3f",pTA1_PMin,pTA1_PMax);
  TCut ADC1Cut = Form("TA1_Edep>%.1f && TA1_N>0",dECut1);
  TCut ADC2Cut = Form("TA2_Edep>%.1f && TA2_N>0",dECut2);


  TCanvas *c1=new TCanvas("c1","",1000,800);
  c1->Divide(2,2);

  c1->cd(1);
  D->Draw("P0>>hp0(50,0.2,0.7)","","");
  D->Draw("P0>>hp0TA1(50,0.2,0.7)","TA1_P>0","same");
  TH1 *hp0=(TH1*)(gROOT->FindObject("hp0"));
  TH1 *hp0TA1=(TH1*)(gROOT->FindObject("hp0TA1"));
  hp0->SetTitle("Thrown (black) and Layer1 detected (red); P0 (GeV) ");
  hp0TA1->SetLineColor(2);

  c1->cd(2);
  TH1* hpene=hp0TA1->Clone("hpene");
  hpene->Divide(hp0);
  hpene->SetMarkerStyle(2);
  hpene->SetMarkerColor(2);
  hpene->SetTitle("Target Chamber Penetration Efficiency; P0 (GeV) ");
  hpene->Draw("prof");

  c1->cd(3);
  D->Draw("TA1_P>>hp1(50,0.1,0.6)",TA1PCut,"");
  TH1 *hp1=(TH1*)(gROOT->FindObject("hp1"));
  hp1->SetTitle("Detected (black) accepted(red); TA1_P (GeV)");

  D->Draw("TA1_P>>hp1det(50,0.1,0.6)",TA1PCut && ADC1Cut && ADC2Cut,"same");
  TH1 *hp1det=(TH1*)(gROOT->FindObject("hp1det"));
  hp1det->SetLineColor(2);


  c1->cd(4);
  TH1 *hdeteff = hp1det->Clone("hdeteff");
  hdeteff->Divide(hp1);
  hdeteff->SetMarkerStyle(2);
  hdeteff->SetMarkerColor(2);
  hdeteff->SetTitle(Form("Detector Eff: ADC1Cut=%.1f MeV, ADC2Cut=%.1f MeV; TA1_P (GeV) ",dECut1,dECut2));
  hdeteff->Draw("prof");

  c1->SaveAs(Form("ThirdArmEff_E%.3f.png",Beam));
}

