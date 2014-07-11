//this script is used to find angle acceptance such that one can 
//determeine the spectrometer's angle
#include "stdlib.h"
#include <iostream>
#include "math.h"
using namespace std;

#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TQObject.h>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TString.h"
#include "TCut.h"
#include "TCutG.h"
#include "TPaveText.h"
#include "TText.h"
#include "TPad.h"
#include "TF1.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TGraphErrors.h"

char* GetKeyWords(int &SetupHMS, double &Beam, double &ProtonAngle, double &GammaAngle)
{
 TTree *gp = (TTree*) gROOT->FindObject("gp");
  
  const double deg = atan(1.0)/45.;
  
  //int SetupHMS=0;
  //double Beam, ProtonAngle, GammaAngle;
  gp->SetBranchAddress("SetupHMS",&SetupHMS);
  gp->SetBranchAddress("Beam",&Beam);
  gp->SetBranchAddress("GammaAngle",&GammaAngle);
  gp->SetBranchAddress("ProtonAngle",&ProtonAngle);
  
  gp->GetEntry(0);
  
  static char *key = new char [255];  
  sprintf(key,"g%.0f_pr%.0f_E%.1f",GammaAngle/deg,ProtonAngle/deg,Beam);
  cout<<"key="<<key<<endl;
  return key;
}

double GetRCSTheta_cm(double Ei, double Theta_lab_rad)
{
  const double M = 0.9383;
  double cosTheta = cos(Theta_lab_rad);
  double Ef = Ei/(1+Ei/M*(1.0-cosTheta));
  double t  = -2.0*Ei*Ef*(1.0-cosTheta);
  double s  = M*M + 2.0*M*Ei;
  double sinHalfTheta_cm = sqrt(-t*s)/(s-M*M);
  double Theta_cm = 2.0*asin(sinHalfTheta_cm);
  return Theta_cm;
}

//this routine only work for gp tree
//void gpAngleAcc(int nbin=15, double th_g_deg_min=10, double th_g_deg_max=40, double p0_gev=0.0)
void gpAngleAcc(double p0_gev=0.0)
{
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadGridX(true);  
  gStyle->SetPadGridY(true);
  
  TTree *gp = (TTree*) gROOT->FindObject("gp");
  
  const double deg = atan(1.0)/45.;
  
  int SetupHMS=0;
  double Beam, ProtonAngle, GammaAngle;
  char *key = GetKeyWords(SetupHMS, Beam, ProtonAngle, GammaAngle); 
  cout<<"key="<<key<<endl;
  ProtonAngle = (ProtonAngle<3.14) ? ProtonAngle/deg : (360-ProtonAngle/deg); 
  GammaAngle = (GammaAngle<3.14) ? GammaAngle/deg : (360-GammaAngle/deg); 


  //due to small acceptance of HMS, 
  //there is no electron-photon coincident events
  if(!SetupHMS)
    {
      TCanvas *c1 = new TCanvas("c1","",800,600);
      gp->Draw("sqrt(pow(Yvb_e-Yvb_g,2)+pow(Xvb_e-Xvb_g,2)):Ei>>hdist","Pvb_e>0 && Pvb_g>0");
      TGraph *gr = (TGraph *)((TGraph *) gROOT->FindObject("Graph"))->Clone("gedist");
      //hdist->SetTitle("Distance of e and #gamma; Ei (GeV); Distance (mm)");
      gr->SetTitle("Distance of e and #gamma; Ei (GeV); Distance (mm)");
      gr->Draw("AP");
      c1->SaveAs(Form("Graph/ge_dist_%s.png",key));
    }

  char xbin[200];
  //sprintf(xbin,"15,10,40");  //for photon arm at 25 deg
  //sprintf(xbin,"15,20,50");  //for photon arm at 35 deg
  int nbin=15; double th_g_deg_min=10, th_g_deg_max=40;
  th_g_deg_min = GammaAngle - 15;
  th_g_deg_max = GammaAngle + 15;

  sprintf(xbin,"%d,%.0f,%.0f",nbin,th_g_deg_min,th_g_deg_max);    

  char tpbin[200];
  double th_p_deg_min=10, th_p_deg_max=40;
  th_p_deg_min = ProtonAngle - 15;
  th_p_deg_max = ProtonAngle + 15;
  sprintf(tpbin,"%d,%.0f,%.0f",nbin,th_p_deg_min,th_p_deg_max);   

  char tgcmbin[200];
  double th_g_cm_deg_min=50, th_g_cm_deg_max=110;
  double th_g_cm_deg = 50;
  
  if(fabs(GammaAngle-20.)<3) th_g_cm_deg = 65;
  else if(fabs(GammaAngle-20.)<3) th_g_cm_deg = 65;
  else if(fabs(GammaAngle-25.)<3) th_g_cm_deg = 75;
  else if(fabs(GammaAngle-30.)<3) th_g_cm_deg = 85;
  else if(fabs(GammaAngle-35.)<3) th_g_cm_deg = 95;
  else if(fabs(GammaAngle-40.)<3) th_g_cm_deg = 105;
  else if(GammaAngle>43.) th_g_cm_deg = 120;

  th_g_cm_deg_min = th_g_cm_deg - 20;
  th_g_cm_deg_max = th_g_cm_deg + 20;
  sprintf(tgcmbin,"%d,%.0f,%.0f",10,th_g_cm_deg_min,th_g_cm_deg_max);  
  

  char weight[255];
  if(p0_gev<0.1) sprintf(weight,"XS*P0_g*P0_g/3.142");
  else sprintf(weight,"XS*P0_g*P0_g/3.142*(abs(P0_p-%.3f)/%.3f<0.1)",p0_gev,p0_gev);

  TCanvas *c2 = new TCanvas("c2","",800,600);
  
  c2->cd();
  gPad->SetRightMargin(0.12);
  gPad->SetGridx();  gPad->SetGridy();
  gp->Draw("Theta0_cm_g*57.3>>htcm",weight,"");
  TH1F *htcm = (TH1F*) gROOT->FindObject("htcm");
  double new_th_g_cm_deg = htcm->GetMean();
  if(fabs(new_th_g_cm_deg-th_g_cm_deg)>5.0)
    {
      th_g_cm_deg += int((new_th_g_cm_deg-th_g_cm_deg)/5.0)*5.0;
      th_g_cm_deg_min = th_g_cm_deg - 20;
      th_g_cm_deg_max = th_g_cm_deg + 20;
      sprintf(tgcmbin,"%d,%.0f,%.0f",10,th_g_cm_deg_min,th_g_cm_deg_max);  
    }
  gp->Draw(Form("Theta0_cm_g*57.3:Theta0_g*57.3>>htcmlab(%s,%s)",xbin,tgcmbin),weight,"colz text");
  TH2F *htcmlab = (TH2F*) gROOT->FindObject("htcmlab");
  htcmlab->SetMinimum(0.5);
  htcmlab->SetTitle(";#theta_{#gamma}^{lab} (deg);#theta_{#gamma}^{cm} (deg)");
  c2->SaveAs("Graph/RCS_Theta_cmVslab.png");

  gp->Draw(Form("Theta0_p*57.3:Theta0_g*57.3>>htptg(%s,%s)",xbin,tpbin),weight,"colz text");
  TH2F *htptg = (TH2F*) gROOT->FindObject("htptg");
  htptg->SetMinimum(0.5);
  htptg->SetTitle(";#theta_{#gamma}^{lab} (deg);#theta_{p}^{lab} (deg)");
  c2->SaveAs("Graph/RCS_Theta_lab.png");

  gp->Draw(Form("P0_g:Theta0_g*57.3>>hpgtg(%s,30,1,4)",xbin),weight,"colz text"); 
  TH2F *hpgtg = (TH2F*) gROOT->FindObject("hpgtg");
  hpgtg->SetMinimum(0.5);
  hpgtg->SetTitle(";#theta_{#gamma}^{lab} (deg);P_{#gamma} (GeV/c)");
  c2->SaveAs("Graph/RCS_Pg_Thetag.png");


  gp->Draw("P0_p>>hpp","","");
  TH1F *hpp = (TH1F*) gROOT->FindObject("hpp");
  double pp_mean = hpp->GetMean();
  char ppbin[255];
  sprintf(ppbin,"%d,%.1f,%.1f",20,pp_mean-2,pp_mean+2);  
  gp->Draw(Form("P0_p:Theta0_g*57.3>>hpptg(%s,%s)",xbin,ppbin),weight,"colz text");
  TH2F *hpptg = (TH2F*) gROOT->FindObject("hpptg");
  hpptg->SetMinimum(0.5);
  hpptg->SetTitle(";#theta_{#gamma}^{lab} (deg);P_{p} (GeV/c)");
  c2->SaveAs("Graph/RCS_Pp_Thetag.png");
 

  gp->Draw(Form("Ei:Theta0_g*57.3>>heitg(%s,40,1,5)",xbin),weight,"colz text"); 
  TH2F *heitg = (TH2F*) gROOT->FindObject("heitg");
  heitg->SetMinimum(0.5);
  heitg->SetTitle(";#theta_{#gamma}^{lab} (deg);Incident Energy (GeV)");
  c2->SaveAs("Graph/RCS_Ei_Thetag.png");

  gp->Draw(Form("t:Theta0_g*57.3>>httg(%s,10,-5.0,0)",xbin),weight,"colz text"); 
  TH2F *httg = (TH2F*) gROOT->FindObject("httg");
  httg->SetMinimum(0.5);
  httg->SetTitle(";#theta_{#gamma}^{lab} (deg);t (GeV^{2})");
  c2->SaveAs("Graph/RCS_t_Thetag.png");

  gp->Draw(Form("u:Theta0_g*57.3>>hutg(%s,10,-10.0,0.0)",xbin),weight,"colz text"); 
  TH2F *hutg = (TH2F*) gROOT->FindObject("hutg");
  hutg->SetMinimum(0.5);
  hutg->SetTitle(";#theta_{#gamma}^{lab} (deg);u (GeV^{2})");
  c2->SaveAs("Graph/RCS_u_Thetag.png");

  gp->Draw(Form("s:Theta0_g*57.3>>hstg(%s,10,0,20)",xbin),weight,"colz text"); 
  TH2F *hstg = (TH2F*) gROOT->FindObject("hstg");
  hstg->SetMinimum(0.5);
  hstg->SetTitle(";#theta_{#gamma}^{lab} (deg);s (GeV^{2})");
  c2->SaveAs("Graph/RCS_s_Thetag.png");
  ///////////////////////////////////////

  TCanvas *c3 = new TCanvas("c3","",1100,800);
  c3->Divide(3,2);
  c3->cd(1);  
  gPad->SetRightMargin(0.12);htcmlab->Draw("colz text");
  c3->cd(2);  
  gPad->SetRightMargin(0.12);htptg->Draw("colz text");
  c3->cd(3);  
  gPad->SetRightMargin(0.12);hpgtg->Draw("colz text");
  c3->cd(4);  
  gPad->SetRightMargin(0.12);hpptg->Draw("colz text");
  c3->cd(5);  
  gPad->SetRightMargin(0.12);httg->Draw("colz text");
  c3->cd(6);  
  gPad->SetRightMargin(0.12);hutg->Draw("colz text");

  if(p0_gev>0.1) c3->SaveAs(Form("Graph/RCS_Kine_%s_P0%.3f.png",key,p0_gev));
  else c3->SaveAs(Form("Graph/RCS_Kine_%s.png",key));
}


void pr1_g2(double beam=6.6){
  TCanvas *c1 = new TCanvas("c1","",800,600);
  TTree *track1 = (TTree*) gROOT->FindObject("track1");
  track1->AddFriend("track0");

  track1->Draw("Theta0*57.3:Phi0*57.3>>hTP0","","");
  TH2F* hTP0 = (TH2F*) gROOT->FindObject("hTP0");
  hTP0->SetMarkerColor(1);
  track1->Draw("Theta0*57.3:Phi0*57.3>>hTP1","track1.Pvb>0","same");
  TH2F* hTP1 = (TH2F*) gROOT->FindObject("hTP1");
  hTP1->SetMarkerColor(2);
  track1->Draw("Theta0*57.3:Phi0*57.3>>hTP2","track1.Pvb>0 && track0.Pvb>0","same");
  TH2F* hTP2 = (TH2F*) gROOT->FindObject("hTP2");
  hTP2->SetMarkerColor(4);
  hTP2->SetMarkerStyle(2);
  hTP0->SetTitle("Thrown(black), #color[2]{#gamma detected (red)}, #color[4]{#gamma-p detected (blue)};#phi_{#gamma} (deg);#theta_{#gamma} (deg)");

  hTP0->Draw("");
  hTP1->Draw("same");
  hTP2->Draw("same");

  c1->SaveAs(Form("Graph/gp_acc_E%.1f.png",beam));
 }

double GetBeamEnergy()
{
  double pBeam=0.0;
  TTree *config = (TTree*) gROOT->FindObject("config");
  config->SetBranchAddress("Beam",&pBeam);
  config->GetEntry(0);
  cout<< "The beam energy for this run is "<< pBeam <<" GeV"<<endl;
  return pBeam;
}

double GetConfigLeaf(const char *leaf="VDAngle")
{
  double val=0.0;
  TTree *config = (TTree*) gROOT->FindObject("config");
  config->SetBranchAddress(leaf,&val);
  config->GetEntry(0);
  //cout<< "In this run,  "<< leaf <<" = "<<val<<endl;
  return val;
}

int GetConfigLeaf_Int(const char *leaf="VDAngle")
{
  int val=0;
  TTree *config =  (TTree*) gROOT->FindObject("config");
  config->SetBranchAddress(leaf,&val);
  config->GetEntry(0);
  //cout<< "In this run,  "<< leaf <<" = "<<val<<endl;
  return val;
}


void CheckAngleAcc(double Eimin=4.0, double Eimax=6.5){
  gROOT->Reset();
  double beam =  GetBeamEnergy();
  if(beam<1) 
    {
      cout<<"\nWrong ntuple file! config tree not exist...I quit...\n";
      return;
    }

  /////////////////////////////////////////////////////////////
  const double deg = atan(1.0)/45.;
  double PhotonAngle = GetConfigLeaf("VDAngle")/deg;

  //By Jixie: These 2 line will change the value of pSBS
  //There must be some bugs in root cint, I have to set pSBS value after these  2 lines
  double ProtonAngle = 0;
  if(PhotonAngle>180.) ProtonAngle = GetConfigLeaf("SuperBigBiteAngle")/deg;
  else ProtonAngle = GetConfigLeaf("HMSAngle")/deg;

  int pSBS = (PhotonAngle>180.) ? 1 : 0;
  //cout<<" pSBS="<<pSBS<<endl;

  const char* ProtonDetector = (pSBS==1)?"SBS":"HMS";
  const char* GammaDetector = (pSBS==1)?"LAC":"NPS";
  //cout<<" pSBS="<<pSBS<<endl;

  //By Jixie: These 2 line will change the value of pSBS
  //There must be some bugs in root cint, I have to set pSBS value after these  2 lines
  //if(pSBS==1) ProtonAngle = GetConfigLeaf("SuperBigBiteAngle")/deg;
  //else ProtonAngle = GetConfigLeaf("HMSAngle")/deg;
  //cout<<" pSBS="<<pSBS<<endl;

  char target[255];
  if(pSBS==1) sprintf(target,"Theta0*57.3:fmod(Phi0*57.3+360.,360)");
  else sprintf(target,"Theta0*57.3:Phi0*57.3");
  
  cout<<"GammaAngle="<<PhotonAngle<<"  "<<ProtonDetector<<"_Angle="<<ProtonAngle
      <<" tgstr="<<target<<endl;

  /////////////////////////////////////////////////////////////
  TTree *track0 = (TTree*) gROOT->FindObject("track0");
  track0->AddFriend("track1");
  track0->AddFriend("track2");


  TCut EiCut = Form("(Ei>%.3f && Ei<%.3f)",Eimin,Eimax);

  TCanvas *c1 = new TCanvas("c1","",0,20,800,600);

  track0->Draw(Form("%s>>hTPFrame",target),"","");
  TH2F *hTPFrame=(TH2F*) gROOT->FindObject("hTPFrame");
  double x1,x2,y1,y2;
  x1=hTPFrame->GetXaxis()->GetXmin();
  x2=hTPFrame->GetXaxis()->GetXmax();
  y1=hTPFrame->GetYaxis()->GetXmin();
  y2=hTPFrame->GetYaxis()->GetXmax()*1.05;
  char bin2D[255];
  sprintf(bin2D,"40,%.2f,%.2f,40,%.2f,%.2f",x1,x2,y1,y2);

  c1->Clear();
  track0->Draw(Form("%s>>hTP0(%s)",target,bin2D),EiCut && "","");
  TGraph* gr0 = (TGraph*) gROOT->FindObject("Graph")->Clone("gr0");
  TH2F *hTP0=(TH2F*) gROOT->FindObject("hTP0");
  hTP0->SetMarkerColor(1);
  gr0->SetMarkerColor(1);

  track0->Draw(Form("%s>>hTP1(%s)",target,bin2D),EiCut && "track0.Pvb>0","");  
  TGraph* gr1 = (TGraph*) gROOT->FindObject("Graph")->Clone("gr1");
  TH2F *hTP1=(TH2F*) gROOT->FindObject("hTP1");
  hTP1->SetMarkerColor(2);  
  gr1->SetMarkerColor(2);

  track0->Draw(Form("%s>>hTP3(%s)",target,bin2D),EiCut && "track1.Pvb>0","");  
  TGraph* gr3 = (TGraph*) ((TGraph*) gROOT->FindObject("Graph"))->Clone("gr3");
  TH2F *hTP3=(TH2F*) gROOT->FindObject("hTP3");
  hTP3->SetMarkerColor(3);  hTP3->SetMarkerStyle(2);
  gr3->SetMarkerColor(3);   gr3->SetMarkerStyle(2);

  track0->Draw(Form("%s>>hTP2(%s)",target,bin2D),EiCut && "track1.Pvb>0 && track0.Pvb>0","");
  TGraph* gr2 = (TGraph*) gROOT->FindObject("Graph")->Clone("gr2");
  TH2F *hTP2=(TH2F*) gROOT->FindObject("hTP2");
  hTP2->SetMarkerColor(4);  hTP2->SetMarkerStyle(2);
  gr2->SetMarkerColor(4);   gr2->SetMarkerStyle(2);
  
  
  TPaveText *pt2 = new TPaveText(0.4,0.83,0.85,0.898,"brNDC");
  pt2->SetFillStyle(4000);
  pt2->SetBorderSize(0);
  pt2->SetFillColor(0);
  if(pSBS)
    pt2->AddText(Form("#color[2]{%s = %.0f^{o}}  #color[3]{%s = %.0f^{o}}",
		      GammaDetector,PhotonAngle-360.,ProtonDetector,ProtonAngle));
  else
    pt2->AddText(Form("#color[2]{%s = %.0f^{o}}  #color[3]{%s = %.0f^{o}}",
		      GammaDetector,PhotonAngle,ProtonDetector,ProtonAngle-360.));

  /*
  hTP0->SetTitle("Thrown(black), #color[2]{#gamma detected (red)}, #color[3]{p detected (green)}, #color[4]{#gamma-p detected (blue)};#phi_{#gamma} (deg);#theta_{#gamma} (deg)");
  hTP0->Draw("");
  hTP1->Draw("same");
  hTP3->Draw("same");
  hTP2->Draw("same");
  pt2->Draw("same");
  c1->SaveAs(Form("Graph/gp_acc_g%.0f_pr%.0f_E%.1f_histo.png",PhotonAngle,ProtonAngle,beam));
  */

  gr0->SetTitle("Thrown(black), #color[2]{#gamma detected (red)}, #color[3]{p detected (green)}, #color[4]{#gamma-p detected (blue)};#phi_{#gamma} (deg);#theta_{#gamma} (deg)");
  gr0->Draw("AP");
  gr1->Draw("psame");
  gr3->Draw("psame");
  gr2->Draw("psame");
  pt2->Draw("same");


  c1->SaveAs(Form("Graph/gp_acc_g%.0f_pr%.0f_E%.1f.png",PhotonAngle,ProtonAngle,beam));

  
  TCanvas *c11 = new TCanvas("c11","",700,40,600,400);

  c11->cd();  
  track0->Draw("Theta0*57.3>>htg00",EiCut && "","");
  track0->Draw("Theta0*57.3>>htg10",EiCut && "track1.Pvb>0","same");
  track0->Draw("Theta0*57.3>>htg01",EiCut && "track0.Pvb>0","same");
  track0->Draw("Theta0*57.3>>htg11",EiCut && "track0.Pvb>0 && track1.Pvb>0","same");
  TH1F *htg00=(TH1F*) gROOT->FindObject("htg00");
  TH1F *htg01=(TH1F*) gROOT->FindObject("htg01");
  TH1F *htg10=(TH1F*) gROOT->FindObject("htg10");
  TH1F *htg11=(TH1F*) gROOT->FindObject("htg11");
  htg00->SetTitle("Thrown(black), #color[2]{#gamma detected (red)}, #color[3]{p detected (green)}, #color[4]{#gamma-p detected (blue)};#theta_{#gamma} (deg)");
  htg00->SetLineColor(1); htg00->SetLineWidth(2);
  htg01->SetLineColor(2); htg01->SetLineWidth(2);
  htg10->SetLineColor(3); htg10->SetLineWidth(2);
  htg11->SetLineColor(4); htg11->SetLineWidth(2);
  
  TPaveText *pt1 = new TPaveText(0.4,0.83,0.85,0.898,"brNDC");
  pt1->SetFillStyle(4000);
  pt1->SetBorderSize(0);
  pt1->SetFillColor(0);
  if(pSBS)
    pt1->AddText(Form("#color[2]{%s = %.0f^{o}}  #color[3]{%s = %.0f^{o}}",
		      GammaDetector,PhotonAngle-360.,ProtonDetector,ProtonAngle));
  else
    pt1->AddText(Form("#color[2]{%s = %.0f^{o}}  #color[3]{%s = %.0f^{o}}",
		      GammaDetector,PhotonAngle,ProtonDetector,ProtonAngle-360.));
  pt1->Draw("same");
  gPad->SetLogy(1); 
  c11->SaveAs(Form("Graph/thetag_acc_g%.0f_pr%.0f_E%.1f.png",PhotonAngle,ProtonAngle,beam));

  TCanvas *c12 = new TCanvas("c12","",700,500,600,400);
  c12->cd();
  track0->Draw("track1.Theta0*57.3>>htp00",EiCut && "","");
  track0->Draw("track1.Theta0*57.3>>htp10",EiCut && "track1.Pvb>0","same");
  track0->Draw("track1.Theta0*57.3>>htp01",EiCut && "track0.Pvb>0","same");
  track0->Draw("track1.Theta0*57.3>>htp11",EiCut && "track0.Pvb>0 && track1.Pvb>0","same");
  TH1F *htp00=(TH1F*) gROOT->FindObject("htp00");
  TH1F *htp01=(TH1F*) gROOT->FindObject("htp01");
  TH1F *htp10=(TH1F*) gROOT->FindObject("htp10");
  TH1F *htp11=(TH1F*) gROOT->FindObject("htp11");
  htp00->SetTitle("Thrown(black), #color[2]{#gamma detected (red)}, #color[3]{p detected (green)}, #color[4]{#gamma-p detected (blue)};#theta_{p} (deg)");
  htp00->SetLineColor(1); htp00->SetLineWidth(2);
  htp01->SetLineColor(2); htp01->SetLineWidth(2);
  htp10->SetLineColor(3); htp10->SetLineWidth(2);
  htp11->SetLineColor(4); htp11->SetLineWidth(2);
  
  pt1->Draw("same");
  gPad->SetLogy(1);

  c12->SaveAs(Form("Graph/thetap_acc_g%.0f_pr%.0f_E%.1f.png",PhotonAngle,ProtonAngle,beam));

  /*
  TCanvas *c2 = new TCanvas("c2","",800,200,800,600);

  c2->cd(); gPad->SetRightMargin(0.12);
  track0->Draw("sqrt(pow(track2.Yvb-track0.Yvb,2.0)+pow(track2.Xvb-track0.Xvb,2.0)):Ei>>hdis","track0.Pvb>0 && track1.Pvb>0 && track2.Pvb>0","*");
  TH2F *hdis=(TH2F*) gROOT->FindObject("hdis");
  hdis->SetTitle("distance between photon and electron ; Ei (GeV); Distance (mm)");
  hdis->Draw("colz");
  c2->SaveAs(Form("Graph/ge_dist_E%.1f.png",beam));
 
  c2->Clear();
  c2->Divide(2,1);
  c2->cd(1);
  track0->Draw("Theta0*57.3>>hht",EiCut && "track1.Pvb>0","");
  TH1F *hht=(TH1F*) gROOT->FindObject("hht");
  hht->SetTitle("p detected; #theta_{#gamma} (deg)");
  c2->cd(2);
  if(pSBS==1) track0->Draw("fmod(Phi0*57.3+360.,360.)>>hhp",EiCut && "track1.Pvb>0","");
  else   track0->Draw("Phi0*57.3>>hhp","track1.Pvb>0","");
  TH1F *hhp=(TH1F*) gROOT->FindObject("hhp");
  hhp->SetTitle("p detected; #phi_{#gamma} (deg)");

  c2->Modified();
  c2->Update();
  c2->SaveAs(Form("Graph/Expected_TP_gamma_E%.1f.png",beam));
  
  TCanvas *c3 = new TCanvas("c3","",800,900);
  c3->Divide(1,2);
  c3->cd(1);
  gPad->SetLogz();
  gPad->SetRightMargin(0.1);
  track0->Draw("track1.Theta0*57.3:track0.Theta0*57.3>>ht0t1(25,10,60,30,10,70)","(EiCut && Theta0*57.3>5 && Pvb>1.0 && track1.Pvb>0.20)","colztext");
  TH2F *ht0t1=(TH2F*) gROOT->FindObject("ht0t1");
  ht0t1->SetTitle("Not weighted;#theta_{#gamma} (deg); #theta_{p} (deg)");
  c3->cd(2);
  gPad->SetLogz();
  gPad->SetRightMargin(0.1);
  track0->Draw("track1.Theta0*57.3:track0.Theta0*57.3>>ht0t1w(25,10,60,30,10,70)","ElasXS*P0*P0*(EiCut && Theta0*57.3>5 && Pvb>1.0 && track1.Pvb>0.20)","colztext");
  TH2F *ht0t1w=(TH2F*) gROOT->FindObject("ht0t1w");
  ht0t1w->SetTitle("Weighted by RCS XS;#theta_{#gamma} (deg); #theta_{p} (deg)");
  c3->SaveAs(Form("Graph/Thrown_T0T1_cmp_E%.1f.png",beam));
  */
 
}

void FitGauss(TH1F *h1)
{
  if(!h1) return;
  h1->Fit("gaus","","Q");
  TF1 *f = (TF1 *)h1->GetListOfFunctions()->FindObject("gaus");
  double mean=f->GetParameter(1);
  double sigma=f->GetParameter(2);
  
  h1->Fit(f,"","R+",mean-1.5*sigma,mean+1.5*sigma);
  h1->Draw();
}


void PhotonRes(){
  gROOT->Reset();
  gStyle->SetStatX(1.0-gStyle->GetPadRightMargin());
  gStyle->SetStatW(0.12);
  double beam =  GetBeamEnergy();
  if(beam<1) 
    {
      cout<<"\nWrong ntuple file! config tree not exist...I quit...\n";
      return;
    }

  /////////////////////////////////////////////////////////////
  const double deg = atan(1.0)/45.;
  double PhotonAngle = GetConfigLeaf("VDAngle")/deg;

  //By Jixie: These 2 line will change the value of pSBS
  //There must be some bugs in root cint, I have to set pSBS value after these  2 lines
  double ProtonAngle = 0;
  if(PhotonAngle>180.) ProtonAngle = GetConfigLeaf("SuperBigBiteAngle")/deg;
  else ProtonAngle = GetConfigLeaf("HMSAngle")/deg;

  int pSBS = (PhotonAngle>180.) ? 1 : 0;
  //cout<<" pSBS="<<pSBS<<endl;

  //const char* ProtonDetector = (pSBS==1)?"SBS":"HMS";
  const char* PhotonDetector = (pSBS==1)?"LAC":"NPS";
  //cout<<" pSBS="<<pSBS<<endl;

  /////////////////////////////////////////////////////////////
  TTree *track0 = (TTree*) gROOT->FindObject("track0");
  track0->AddFriend("track1");


  double PRes = 0.002, ThetaRes = 0.001;
  if(pSBS==1) PRes = 0.005;

  TCanvas *c3 = new TCanvas("c3","",600,900);
  c3->Divide(1,2);
  c3->cd(1);
  TH1F *htempT=(TH1F*) gROOT->FindObject("htempT");
  if(htempT) delete htempT;
  track0->Draw("track1.Theta0>>htempT","track1.Pvb>1.0");
  htempT=(TH1F*) gROOT->FindObject("htempT");
  double pMeanAngle=htempT->GetMean();
  TCut dThetaCut = Form("abs(track1.Theta0-%.3f) < %.4f && track1.Pvb>1", pMeanAngle, ThetaRes);
  c3->cd(2);
  TH1F *htempP=(TH1F*) gROOT->FindObject("htempP");
  if(htempP) delete htempP;
  track0->Draw("track1.P0>>htempP",dThetaCut);
  htempP=(TH1F*) gROOT->FindObject("htempP");
  double pMeanP=htempP->GetMean();
  TCut dPCut = Form("abs(track1.P0-%.3f)/%.3f < %.3f",pMeanP,pMeanP,PRes);
		       
  char strTitle[255];
  sprintf(strTitle,"P_{p}=%.3f +/- %.1f%%, #theta_{p}=%.4f +/- %.4f",
	  pMeanP,PRes*100,pMeanAngle,ThetaRes);
  TCut resCut =  dThetaCut + dPCut ;
  resCut += "track0.Pvb>0" ;

  c3->cd(1);
  track0->Draw("Theta0_tr",resCut,"");
  c3->cd(2);
  track0->Draw("Phi0_tr",resCut,"");

  TCanvas *c4 = new TCanvas("c4","",600,900);
  c4->Divide(1,3);
  c4->cd(1); 

  c4->cd(1); 
  track0->Draw("Ei>>hEiRange",resCut,"");
  TH1F *hEiRange=(TH1F*) gROOT->FindObject("hEiRange");
  hEiRange->SetTitle(Form("%s; E_{i} (GeV);",strTitle));
  FitGauss(hEiRange);
		     
  c4->cd(2); 
  track0->Draw("P0>>hPRange",resCut,"");
  TH1F *hPRange=(TH1F*) gROOT->FindObject("hPRange");
  hPRange->SetTitle(Form("%s; P_{#gamma} (GeV);",strTitle));
  FitGauss(hPRange);
		    
  c4->cd(3); 
  track0->Draw("Theta0>>hTRange",resCut,"");
  TH1F *hTRange=(TH1F*) gROOT->FindObject("hTRange");
  hTRange->SetTitle(Form("%s; #theta_{#gamma} (rad);",strTitle));
  FitGauss(hTRange);
		    
  c4->SaveAs(Form("Graph/%s_Res_%.0f.png",PhotonDetector,PhotonAngle));
  
}

