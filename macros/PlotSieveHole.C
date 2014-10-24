//script to plot variants for one single sieve hole by a single click
//left click is for 2x3.5 elipse, middle click is for 3x6 elipse
 
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TQObject.h>
#include <TROOT.h>

char gKey[255];
double gX0,gY0;

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
      config->SetBranchAddress("SeptumCurrentRatioL",&SeptumCurrentRatioL);
      config->SetBranchAddress("SeptumCurrentRatioR",&SeptumCurrentRatioR);
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

//get the index and return the center position
int GetHoleIndex(double x, double y, double &x0, double &y0)
{
  double xbinsL[8]={-24.0,-18.0,-12.0,-6.0,0.0,4.5,9.0,13.5};
  double xbinsR[8]={-13.5,-9.0,-4.5,0.0,6.0,12.0,18.0,24.0};
  double *xbins=xbinsL;
  
  int left=1;
  //determine if it is left or right arm
  TVirtualPad *padsav = gPad;
  pCan->cd(9);
  track0->Draw("Xvb>>hXvb","TrackClass>5");
  padsav->cd();
  double MeanXvb=hXvb->GetMean();  
  if(MeanXvb<-40) left=0;
  if(left==0) xbins=xbinsR;

  //vertical span=13.3mm x 7=93.1mm, from -43.6 to 49.6
  double ybins[8]={-43.6,-30.3,-17.0,-3.7,9.6,22.9,36.2,49.5};
  int h=-1,v=-1;
  for(int i=0;i<8;i++)
    {
      if(x>=xbins[i] && x<xbins[i+1]) {h=i;x0=(xbins[h]+xbins[h+1])/2;break;} 
    }
  for(int j=0;j<8;j++)
    {
      if(y>=ybins[j] && y<ybins[j+1]) {v=j;y0=(ybins[v]+ybins[v+1])/2;break;} 
    }

    int idx=(h<0 || v<0)?-1:10*h+v; 
    cout<<"x="<<x<<" y="<<y<<" ==> holeindex="<<idx<<endl;
    return idx;
}

void VarInHole(TCanvas *pCan,int holeindex=33)
{
  if(!gROOT->FindObject("CUTG")) return;
  TCut all = "TrackClass>5 && Xfp_tr>-20 && abs(Yvb)<100";
  TCut theCut="CUTG && TrackClass>5 && Xfp_tr>-20 && abs(Yvb)<100";
  
  ///////////////////////////////////////////////////////////////////	
 
  TH1* h1=0;
  TH2 *h2=0;
  char dest[255],hname[255];
  /*
  pCan->cd(1);
  track0->Draw("-Xvb_tr:-Yvb_tr >> h2VH",all,"contz");
  h2 = (TH2*) (gROOT->FindObject("h2VH"));
  h2->SetTitle("electrons pass sieve; horizontal (mm) ; vertical (mm) ");
  theHoleCut->Draw("same");
  */

  ////////////////////////////////////
  pCan->cd(2);
  sprintf(hname,"h2VH_vb_%02d",holeindex);
  sprintf(dest,"-Xvb_tr:-Yvb_tr >> %s(12,%.1f,%.1f,12,%.1f,%.1f)",
	  hname,gX0-4.5,gX0+4.5,gY0-13.3,gY0+13.3);
  h2 = (TH2*) (gROOT->FindObject(hname));
  if(h2)  {cout<<"delete "<<h2->GetName();delete h2;}
  track0->Draw(dest,theCut,"contz");
  h2 = (TH2*) (gROOT->FindObject(hname));
  h2->SetTitle(Form("Virtual Boundary: Hole %02d; -Y_{vb}^{tr} (mm) ; -X_{vb}^{tr} (mm) ",holeindex));

  TText *text=0;
  TPaveText *pt = new TPaveText(0.50,0.74,0.85,0.89,"brNDC");
  pt->SetBorderSize(0);
  pt->SetFillColor(0);
  text=pt->AddText(Form("Beam = %.3f",Beam));
  text->SetTextColor(2);
  if(!bIsCombinedTree)
  {  
    text=pt->AddText(Form("Target_A = %.0f",TargetAtomicNumber));
  } 
  else
  {  
    text=pt->AddText("Multi-Targets");	
  }
  text->SetTextColor(2);
  pt->Draw("same");
  
  ////////////////////////////////////
  pCan->cd(3);

  sprintf(hname,"h2TP_vb_%02d",holeindex);
  sprintf(dest,"Thetavb_tr:Phivb_tr >> %s",hname);
  h2 = (TH2*) (gROOT->FindObject(hname));
  if(h2)  {cout<<"delete "<<h2->GetName();delete h2;}
  track0->Draw(dest,theCut,"contz");
  h2 = (TH2*) (gROOT->FindObject(hname));
  h2->SetTitle("Virtual Boundary: #theta_{vb}^{tr} Vs #phi_{vb}^{tr}; #phi_{vb}^{tr} (rad) ; #theta_{vb}^{tr} (rad) ");

  ////////////////////////////////////
  pCan->cd(4);    
  sprintf(hname,"h1Pvb_%02d",holeindex);
  sprintf(dest,"Pvb >> %s",hname);
  h1 = (TH1*) (gROOT->FindObject(hname));
  if(h1)  {cout<<"delete "<<h1->GetName();delete h1;}

  if(gROOT->FindObject("h1Pvb")) delete gROOT->FindObject("h1Pvb");
  track0->Draw("Pvb >> h1Pvb",theCut);
  double MeanPvb=h1Pvb->GetMean();
  if(!bIsCombinedTree)
  {
    sprintf(dest,"Pvb >> %s",hname);
    track0->Draw(Form("Pvb >> %s(100,%.3f,%.3f)",hname,
		      MeanPvb-0.01,MeanPvb+0.01),theCut);
  }
  else
  {
    if(Beam<1.3)
      track0->Draw(Form("Pvb >> %s(60,%.3f,%.3f)",hname,
			MeanPvb-0.02,MeanPvb+0.01),theCut);
    else if (Beam<1.8)
      track0->Draw(Form("Pvb >> %s(100,%.3f,%.3f)",hname,
			MeanPvb-0.04,MeanPvb+0.01),theCut);
    else
      track0->Draw(Form("Pvb >> %s(140,%.3f,%.3f)",hname,
			MeanPvb-0.06,MeanPvb+0.01),theCut);
  }
  track0->Draw(dest,theCut);
  h1 = (TH1*) (gROOT->FindObject(hname));
  h1->SetTitle("P0(black), P_{vb}(red) and P_{rec}(blue); P (GeV)");
  h1->SetLineColor(2);
  
  track0->Draw("P0",theCut,"same");
  
  TH1 *h1Prec=0; 
  h1Prec = (TH1*) (gROOT->FindObject("h1Prec"));
  if(h1Prec) delete h1Prec;
  track0->Draw("P_rec >> h1Prec",theCut,"same");
  h1Prec = (TH1*) (gROOT->FindObject("h1Prec"));
  h1Prec->SetLineColor(4);

  ////////////////////////////////////
  pCan->cd(5);
  sprintf(hname,"h2Theta0Phi0_%02d",holeindex);
  sprintf(dest,"Theta0*57.3:Phi0*57.3 >> %s",hname);
  h2 = (TH2*) (gROOT->FindObject(hname));
  if(h2)  {cout<<"delete "<<h2->GetName();delete h2;}
  track0->Draw(dest,theCut,"contz");
  h2 = (TH2*) (gROOT->FindObject(hname));
  h2->SetTitle("Target: #theta_{0} Vs #phi_{0} ; #phi_{0} (deg) ;#theta_{0} (deg) ");

  ////////////////////////////////////
  pCan->cd(6);
  sprintf(hname,"h2Theta0Phi0_tr_%02d",holeindex);
  sprintf(dest,"Theta0_tr:Phi0_tr >> %s",hname);
  h2 = (TH2*) (gROOT->FindObject(hname));
  if(h2)  {cout<<"delete "<<h2->GetName();delete h2;}
  track0->Draw(dest,theCut,"contz");
  h2 = (TH2*) (gROOT->FindObject(hname));
  h2->SetTitle("Target: #theta_{0}^{tr} Vs #phi_{0}^{tr} ; #phi_{0}^{tr} (rad) ;#theta_{0}^{tr} (rad) ");

  ////////////////////////////////////
  pCan->cd(7);
  sprintf(hname,"h2Y0X0_%02d",holeindex);
  sprintf(dest,"Y0:X0 >> %s",hname);
  h2 = (TH2*) (gROOT->FindObject(hname));
  if(h2)  {cout<<"delete "<<h2->GetName();delete h2;}
  track0->Draw(dest,theCut,"contz");
  h2 = (TH2*) (gROOT->FindObject(hname));
  h2->SetTitle("Target: Y0 Vs X0; X0 (mm) ; Y0 (mm) ");

  ////////////////////////////////////
  pCan->cd(8);
  sprintf(hname,"h2Y0X0_tr_%02d",holeindex);
  sprintf(dest,"-X0_tr:-Y0_tr >> %s",hname);
  h2 = (TH2*) (gROOT->FindObject(hname));
  if(h2)  {cout<<"delete "<<h2->GetName();delete h2;}
  track0->Draw(dest,theCut,"contz");
  h2 = (TH2*) (gROOT->FindObject(hname));
  h2->SetTitle("Target: -X_{0}^{tr} Vs -Y_{0}^{tr}; -Y_{0}^{tr} (mm) ; -X_{0}^{tr} (mm) ");

  ////////////////////////////////////
  pCan->cd(9);
  sprintf(hname,"h2VH_proj2tg_tr_%02d",holeindex);
  sprintf(dest,"-X_proj2tg_tr:-Y_proj2tg_tr >> %s",hname);
  h2 = (TH2*) (gROOT->FindObject(hname));
  if(h2)  {cout<<"delete "<<h2->GetName();delete h2;}
  track0->Draw(dest,theCut,"contz");
  h2 = (TH2*) (gROOT->FindObject(hname));
  h2->SetTitle("Project to Target: -X_{tg}^{tr} Vs -Y_{tg}^{tr}; -Y_{tg}^{tr} (mm) ; -X_{tg}^{tr} (mm) ");

  ////////////////////////////////////
  pCan->cd(10);
  sprintf(hname,"h2VH_proj2sl_%02d",holeindex);
  sprintf(dest,"-X_proj2sl_tr:-Y_proj2sl_tr >> %s",hname);
  h2 = (TH2*) (gROOT->FindObject(hname));
  if(h2)  {cout<<"delete "<<h2->GetName();delete h2;}
  track0->Draw(dest,theCut,"contz");
  h2 = (TH2*) (gROOT->FindObject(hname));
  h2->SetTitle("Project to Sieve: -X_{sl}^{tr} Vs -Y_{sl}^{tr}; -Y_{sl}^{tr} (mm) ; -X_{sl}^{tr} (mm) ");

  ////////////////////////////////////
  pCan->cd(11);
  sprintf(hname,"h2VH_fp_%02d",holeindex);
  sprintf(dest,"-Xfp_tr:-Yfp_tr >> %s",hname);
  h2 = (TH2*) (gROOT->FindObject(hname));
  if(h2)  {cout<<"delete "<<h2->GetName();delete h2;}
  track0->Draw(dest,theCut,"contz");
  h2 = (TH2*) (gROOT->FindObject(hname));
  h2->SetTitle(Form("Focus Plane: Hole %02d; -Yfp (mm) ; -Xfp (mm) ",holeindex));

  ////////////////////////////////////
  pCan->cd(12);
  sprintf(hname,"h2TP_fp_%02d",holeindex);
  sprintf(dest,"Thetafp_tr:Phifp_tr >> %s",hname);
  h2 = (TH2*) (gROOT->FindObject(hname));
  if(h2)  {cout<<"delete "<<h2->GetName();delete h2;}
  track0->Draw(dest,theCut,"contz");
  h2 = (TH2*) (gROOT->FindObject(hname));
  h2->SetTitle("Focus Plane: #theta_{fp}^{tr} Vs #phi_{fp}^{tr}; #phi_{fp}^{tr} (rad) ; #theta_{fp}^{tr} (rad) ");

  ////////////////////////////////////
  pCan->cd(13);
  sprintf(hname,"h2VH_rec_%02d",holeindex);
  sprintf(dest,"-X_rec_tr:-Y_rec_tr >> %s",hname);
  h2 = (TH2*) (gROOT->FindObject(hname));
  if(h2)  {cout<<"delete "<<h2->GetName();delete h2;}
  track0->Draw(dest,theCut,"contz");
  h2 = (TH2*) (gROOT->FindObject(hname));
  h2->SetTitle("Reconstructed: -X_{rec}^{tr} Vs -Y_{rec}^{tr}; -Y_{rec}^{tr} (mm) ; -X_{rec}^{tr} (mm) ");

  ////////////////////////////////////
  pCan->cd(14);
  sprintf(hname,"h2TP_rec_%02d",holeindex);
  sprintf(dest,"Theta_rec_tr:Phi_rec_tr >> %s",hname);
  h2 = (TH2*) (gROOT->FindObject(hname));
  if(h2)  {cout<<"delete "<<h2->GetName();delete h2;}
  track0->Draw(dest,theCut,"contz");
  h2 = (TH2*) (gROOT->FindObject(hname));
  h2->SetTitle("Reconstructed: #theta_{rec}^{tr} Vs #phi_{rec}^{tr}; #phi_{rec}^{tr} (rad) ; #theta_{rec}^{tr} (rad) ");


  ////////////////////////////////////
  pCan->cd(15);
  sprintf(hname,"h1Z_rec_%02d",holeindex);
  sprintf(dest,"Z_rec >> %s",hname);
  h1 = (TH1*) (gROOT->FindObject(hname));
  if(h1)  {cout<<"delete "<<h1->GetName();delete h1;}
  track0->Draw(dest,theCut,"");
  h1 = (TH1*) (gROOT->FindObject(hname));
  h1->SetTitle("Reconstructed Z_{rec}(blue) and Z0(black); Z_{rec} (mm) ");
  h1->SetLineColor(4);
  track0->Draw("Z0",theCut,"same");

  ////////////////////////////////////
  pCan->cd();
  system("mkdir -p graph");
  if(!bIsCombinedTree)
  {  
    pCan->SaveAs(Form("graph/SieveHole%02d_E%.3f_Helm%.0fdeg_R%.1f_A%.0f.png",holeindex,
		      Beam,HelmRotAngle1*180/3.141593,HelmCurrentRatio,TargetAtomicNumber));	
  } 
  else
  { 
    pCan->SaveAs(Form("graph/SieveHole%02d_E%.3f_Helm%.0fdeg_R%.1f_combined.png",holeindex,
		      Beam,HelmRotAngle1*180/3.141593,HelmCurrentRatio));
  }
}
////////////////////////////////////////////////////////////////
 

void PlotSieveHole(char *filename="",char *key="")
{
  // using signal/slot in TCanvas/TPad to get feedback about processed events. 


  if(strlen(filename)>5)
    {
      cout<<"trying to open root file "<<filename<<endl;
      TFile* _file0=TFile::Open(filename);
    }

  if(gROOT->GetListOfFiles()->GetEntries()<1) {
    cout<<"no root file opened yet, quit ... \n"; 
    return;
  }
  ReadConfig();

  TString file=_file0->GetName();
  if(strlen(key)<2)
    {
      int start=0;
      if(file.Last('/')>=0) file.Remove(0,file.Last('/')+1);
      file.Remove(file.Length()-5,5);
      key=file.Data();
    }
  sprintf(gKey,"%s",key);

  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetOptFit(1);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetStatH(0.2);
  gStyle->SetStatW(0.17);
  gStyle->SetStatX(1.0-gStyle->GetPadRightMargin());
	
  
  TCanvas *pCan=new TCanvas("pCan","Sieve Slit",1350,800);
  pCan->Divide(5,3,0.001,0.001);
  
  pCan->cd(1);
  TH2F *h2VH=new TH2F("h2VH",
		      Form("%s; -Yvb_tr (mm); -Xvb_tr (mm)",key),
		      50,-24.0,13.5,42,-43.6,49.5);
  track0->Draw("-Xvb_tr:-Yvb_tr >> h2VH","TrackClass>5","contz"); 

  pCan->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", 0, 0,
		"DoEvent(Int_t,Int_t,Int_t,TObject*)");
}


void DoEvent(Int_t event, Int_t x, Int_t y, TObject *selected)
{
  //do not response to mouse movement
  if(event>20 || event<=10) return;
  
  TCanvas *c = (TCanvas *) gTQSender;
  double px = gPad->AbsPixeltoX(x);
  gX0 = gPad->PadtoX(px);
  double py = gPad->AbsPixeltoY(y);
  gY0 = gPad->PadtoY(py);
  
  //printf("Canvas %s: event=%d, x=%d, y=%d, selected_type=%s, selected_name=%s, p2x=%f, p2y=%f\n", c->GetName(),event, x, y, selected->IsA()->GetName(),selected->GetName(),px,py);
  
  TVirtualPad *padsav = 0;
  
  int holeindex=-1;
  if(strcmp("h2VH",selected->GetName())==0) 
    {
      padsav = gPad;
      if(event==11) DrawElipse(x,y);
      else if(event==12) DrawElipse(x,y,3.0,6.0);
      holeindex=GetHoleIndex(px,py,gX0,gY0);
      if(holeindex>0) VarInHole(c,holeindex);
   }
  else if(strcmp("CUTG",selected->GetName())==0) 
    {
      padsav = gPad;
      holeindex=GetHoleIndex(px,py,gX0,gY0);
      if(holeindex>0) VarInHole(c,holeindex);
    }

  if(padsav) padsav->cd();
}

void DrawElipse(double px,double py,double a=2.0,double b=3.5)
{
  double x = gPad->AbsPixeltoX(px);
  double y = gPad->AbsPixeltoY(py);
 
  //redefine CUTG
  if(gROOT->FindObject("CUTG")) delete gROOT->FindObject("CUTG");
  const int kNpt=9;
  TCutG *cutg = new TCutG("CUTG",kNpt);
  cutg->SetVarX("-Yvb_tr");
  cutg->SetVarY("-Xvb_tr");
  cutg->SetTitle("SelectedHoleCut");
  cutg->SetFillColor(0);
  cutg->SetMarkerStyle(20);
  cutg->SetLineWidth(2);
  
  for(int i=0;i<=kNpt;i++)
    {
      double phi=i*4*acos(0)/kNpt;       
      cutg->SetPoint(i,x+a*cos(phi),y+b*sin(phi));
    }
  
  cutg->Draw("same");

  gPad->Update();
}

