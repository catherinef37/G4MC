
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
Double_t        SeptumCurrentRatioL;
Double_t        SeptumCurrentRatioR;
Double_t        BigBiteAngle;
Double_t        BigBiteTiltAngle;
Double_t        Pivot2BigBiteFace;

bool		bIsCombinedTree;	//to ideneify if this is a combined ntuple



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

      if(config->GetBranch("SeptumCurrentRatioL"))
	{
	  config->SetBranchAddress("SeptumCurrentRatioL",&SeptumCurrentRatioL);
	  config->SetBranchAddress("SeptumCurrentRatioR",&SeptumCurrentRatioR);
	}
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

void SieveSlit(char *filename="", int WeightByXS=1)
{
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetOptFit(1);
  gStyle->SetPadRightMargin(0.15);
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
	
  const int kNRaster=12;
  char strDest[kNRaster][255];
  TText *text=0;
  char RasterCut[kNRaster][255];
  TH2* h2=0;

  TCanvas *pCan=new TCanvas("pCan","Sieve Slit",1200,700);
  pCan->Divide(kNRaster/3,3);
  /*
  TF1 *f1 = new TF1("f1","[0]/x+[1]+[2]*x+[3]*x*x",0.2,3.3);
  f1->SetFillColor(0);
  f1->SetMarkerStyle(20);
  f1->SetMarkerColor(2);
  f1->SetLineColor(2);
  f1->SetLineWidth(3);
  */

  //for(int i=0;i<kNRaster;i++)
  for(int i=0;i<11;i++)
  {
    pCan->cd(i+1);
    sprintf(strDest[i],"-Xvb_tr:-Yvb_tr >> hRaster_%d(25,-24.4,14.4,21,-44,50)",i);
    if(WeightByXS)
      {
	//sprintf(RasterCut[i],"(abs(X0)<%.1f && abs(Y0)<%.1f && TrackClass>4 && abs(Z0+887.74)<1.27)*ElasXS/1000.",i+1.0,i+1.0);
	sprintf(RasterCut[i],"(abs(sqrt(X0*X0+Y0*Y0)-%.1f)<=0.5 && TrackClass>4)*ElasXS/1000.",i+0.0);
      }
    else
      {    
	//sprintf(RasterCut[i],"abs(X0)<%.1f && abs(Y0)<%.1f && TrackClass>4 && abs(Z0+887.74)<1.27",i+1.0,i+1.0);
	sprintf(RasterCut[i],"abs(sqrt(X0*X0+Y0*Y0)-%.1f)<=0.5 && TrackClass>4",i+0.0);
      }
    track0->Draw(strDest[i],RasterCut[i],"contz");
    h2 = (TH2*) (gROOT->FindObject(Form("hRaster_%d",i)));
    h2->SetTitle(Form("%.1f<R_{beam}<=%.1f (mm); horizontal (mm) ; vertical (mm) ",(i>0)?i-0.5:0.0,i+0.5));
    /*  
	TPaveText *pt = new TPaveText(1.5,0.65,3.3,0.80,"brNDC");
	pt->SetBorderSize(0);
	pt->SetFillColor(0);
	text=pt->AddText(Form("Raster %d x %d mm",i,i));
	pt->Draw();
    */
  }

  pCan->cd();
 if(WeightByXS)
   {
     pCan->SaveAs(Form("SeiveSlit_%.2fdeg_E%.3f_VariousRaster_XSWeighted.png",
		       LHRSAngle*180./3.14159, Beam));	
   }
 else
   {  
     pCan->SaveAs(Form("SeiveSlit_%.2fdeg_E%.3f_VariousRaster.png",
		       LHRSAngle*180./3.14159, Beam));	
   }
}
 
////////////////////////////////////////////////////////////////
 
