
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


void ThirdArm(char *filename="")
{
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetOptFit(1);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetStatH(0.2);
  gStyle->SetStatW(0.17);
  //gStyle->SetStatStyle(4000);
  gStyle->SetStatX(1.0-gStyle->GetPadRightMargin());
  //gStyle->SetStatY(0.65);
	

  if(strlen(filename)>5)
    {
      cout<<"trying to open root file "<<filename<<endl;
      TFile* _file0=TFile::Open(filename);
    }
  ReadConfig();

  //========= Macro generated for ThirdArmCut   
  TCutG *cutg = new TCutG("cutg",5);
  cutg->SetVarX("-Yvb_tr");
  cutg->SetVarY("-Xvb_tr");
  cutg->SetTitle("ThirArmCut");
  cutg->SetFillColor(1);
  cutg->SetLineColor(2);
  cutg->SetLineStyle(2);
  cutg->SetLineWidth(3);
  cutg->SetMarkerStyle(20);
   
  //y_up=tan(-54)*x-230
  //y_down=tan(-54)*x-30
  double beta=-54.0/57.3;
  double yoff_u=-10.0,yoff_d=-250.0;

  double bandheight=600.0,bandwidth,xhalf,xshift;
  double x=0;

  xhalf=fabs(bandheight/2.0*cos(beta));
  bandwidth=fabs((yoff_u-yoff_d)*cos(beta));
  xshift=fabs(bandwidth*sin(beta)/2.0);

  cutg->SetPoint(0,x=-xhalf+xshift,tan(beta)*x+yoff_u);
  cutg->SetPoint(1,x=-xhalf-xshift,tan(beta)*x+yoff_d);
  cutg->SetPoint(2,x= xhalf-xshift,tan(beta)*x+yoff_d);
  cutg->SetPoint(3,x= xhalf+xshift,tan(beta)*x+yoff_u);
  cutg->SetPoint(4,x=-xhalf+xshift,tan(beta)*x+yoff_u);


  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);

  //TCut mycut="abs(Pvb-0.35)<0.07 && abs(Yvb_tr)<300 && abs(-Xvb_tr+100)<400";
  TCut pEnergyCut0,pEnergyCut;
  if(Beam<1.2)
    {
      pEnergyCut0="abs(Pvb-0.34)<0.07";
      pEnergyCut="abs(Pvb-0.34)<0.07 && abs(Yvb_tr)<300 && abs(-Xvb_tr+0)<400";
    }  
  else if(Beam<1.8) 
    {   
      pEnergyCut0="abs(Pvb-0.39)<0.07";
      pEnergyCut="abs(Pvb-0.39)<0.07 && abs(Yvb_tr)<300 && abs(-Xvb_tr+100)<400";
      //pEnergyCut0="abs(Theta0*57.3-70.5)<2.5";
      //pEnergyCut="abs(Theta0*57.3-70.5)<2.5 && abs(Yvb_tr)<300 && abs(-Xvb_tr+100)<400";
    }
  else
    {
      //pEnergyCut0="abs(Theta0*57.3-72.5)<2.5";
      //pEnergyCut="abs(Theta0*57.3-72.5)<2.5 && abs(Yvb_tr)<300 && abs(-Xvb_tr+100)<400";
      pEnergyCut0="abs(Pvb-0.39)<0.07";
      pEnergyCut="abs(Pvb-0.39)<0.07 && abs(Yvb_tr)<300 && abs(-Xvb_tr+100)<400";
    }

  TCanvas *c1=new TCanvas("c1","",700,800);
  c1->Divide(2,2,0.01,0.001);


  c1->cd(2);
  TF1 *fcn=new TF1("fcn","[0]+[1]*x",-200,200);
  track0->Draw("-Xvb_tr:-Yvb_tr >> h2(60,-300,300,120,-700,500)",pEnergyCut0,"prof");
  
  double xmin=-50,xmax=100;
  if(Beam>1.2) {xmin=-100;xmax=100;}

  h2->Fit(fcn,"RQ0","",xmin,xmax);
  double yoff=fcn->GetParameter(0);
  beta=atan(fcn->GetParameter(1));
  //y=tan(beta)*x+yoff ==> xoff=yoff/tan(beta);
  double xoff=-yoff/tan(beta);

  h2->Fit(fcn,"R","",xmin+xoff,xmax+xoff);
  yoff=fcn->GetParameter(0);
  beta=atan(fcn->GetParameter(1));
  xoff=-yoff/tan(beta);
  
  cout<<"xoff="<<xoff<<"  yoff="<<yoff<<endl;

  if(Beam<1.2) 
    {
      yoff_u=100,yoff_d=-300;
      yoff_u=500,yoff_d=-150;    //good for 0.27<Pvb<0.41 GeV
    }
  else if(Beam<1.8)
    {
      yoff_u=360,yoff_d=-420;     //good for 0.32<Pvb<0.46 GeV
      //yoff_u=360,yoff_d=-360;    //good for 68<Theta0<73 degrees   
    }
  else if(Beam<2.3)
    {
      yoff_u=70,yoff_d=-210;     //good for 70<Theta0<75 degrees
      yoff_u=130,yoff_d=-230;    //good for 0.30<Pvb<0.48 GeV
    }
  xhalf=fabs(bandheight/2.0*cos(beta));
  bandwidth=fabs((yoff_u-yoff_d)*cos(beta));
  xshift=fabs(bandwidth*sin(beta)/2.0);

  bool pKeepAtHorizontal=true;
  if(pKeepAtHorizontal)
    {
      cutg->SetPoint(4,x=-xhalf+xshift+xoff,tan(beta)*x+yoff_u);
      cutg->SetPoint(0,x=-xhalf+xshift+xoff,tan(beta)*x+yoff_u);
      cutg->SetPoint(1,x=-xhalf-xshift+xoff,tan(beta)*x+yoff_d);
      cutg->SetPoint(2,x= xhalf-xshift+xoff,tan(beta)*x+yoff_d);
      cutg->SetPoint(3,x= xhalf+xshift+xoff,tan(beta)*x+yoff_u);
    }
  else
    {
      cutg->SetPoint(4,x=-xhalf+xshift,tan(beta)*x+yoff_u);
      cutg->SetPoint(0,x=-xhalf+xshift,tan(beta)*x+yoff_u);
      cutg->SetPoint(1,x=-xhalf-xshift,tan(beta)*x+yoff_d);
      cutg->SetPoint(2,x= xhalf-xshift,tan(beta)*x+yoff_d);
      cutg->SetPoint(3,x= xhalf+xshift,tan(beta)*x+yoff_u);
    }
  cutg->Draw();

  c1->cd(1);
  track0->Draw("-Xvb_tr:-Yvb_tr >> hvb(70,-350,350,80,-400,400)",pEnergyCut,"");
  hvb->SetTitle(Form("#splitline{3rd arm: #theta=70^{o}, 2100 mm, %.1f^{o} rotated}{size: %.0fmm X %.0fmm, xoff=%.0f mm, yoff=%.0fmm}; horizontal (mm) ; vertical (mm) ",90.0-fabs(beta*57.3),bandwidth,bandheight,xoff,yoff));
  hvb->Draw("contz");
  cutg->Draw();
  fcn->Draw("same");
  
  c1->cd(3);
  //track0->Draw("-Xvb_tr:-Yvb_tr >> hvb_rej(70,-350,350,80,-400,400)","cutg" && !pEnergyCut0,"contz");
  //hvb_rej->SetTitle("Rejected hits; honrizontal (mm) ; vertical (mm) ");
  track0->Draw("P0 >> hP0(80,0.2,0.6)","cutg" && pEnergyCut0,"");
  hP0->SetTitle("P0(black) and Pvb(red) ; P0 or Pvb (GeV)");
  track0->Draw("Pvb >> hPvb(80,0.2,0.6)","cutg" && pEnergyCut0,"same");
  hPvb->SetLineColor(2);

  c1->cd(4);
  track0->Draw("Theta0*57.3 >> htall(60,60,90)","cutg","");
  track0->Draw("Theta0*57.3 >> ht","cutg" && pEnergyCut,"same");
  ht->SetLineColor(2);
  htall->SetTitle("#theta_{0} all(black) and accepted(red)");

  c1->SaveAs(Form("ThirdArmHits_E%.3f.png",Beam));
}


void ThirdArmCut(char *filename="")
{
  ThirdArm(filename);
}

