//this script is used to find proton arm acceptance (HMS or SBS)
// 1) find proton acceptance only, cosTh vs phi
// 2) find gamma-proton coincident acceptance, cosTh vs phi
// 3) find proton acceptance only, theta_tr vs phi_tr

//determeine the spectrometer's angle
#include "stdlib.h"
#include <iostream>
#include <fstream>
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



void Write(TCutG *cutg, double GammaAngle=330, double ProtonAngle=30)
{
  ofstream fout;
  char key[255], file[255];
  const char *detector[]={"HMS","SBS"};
  int pSBS = (GammaAngle>180.) ? 1 : 0;
  sprintf(key,"g%.0f_p%.0f",GammaAngle,ProtonAngle);
  sprintf(file,"%sSolidAngle_%s.cc",detector[pSBS],key);

  fout.open(file);
  char var[100],sub1[255], sub2[255];
  sprintf(var,"g%sTCutG_%s",detector[pSBS],key);      //gHMSTCutG_g20_p320
  sprintf(sub1,"Init_%sTCutG_%s",detector[pSBS],key);
  sprintf(sub2,"%sSolidAngle_%s",detector[pSBS],key);

  int N = cutg->GetN();

  fout<<"#include <TROOT.h>"<<endl;
  fout<<"#include <TCutG.h>"<<endl;

  fout<<"\n"<<"static TCutG* "<<var<<" = 0;\n"<<endl;

  //the 1st routine: TCutG* Init_TCutG_g20_p320()
  fout<<"TCutG* "<<sub1<<"() {\n" <<endl;
  fout<<"\t"<<"if("<<var<<") return "<<var<<";\n" <<endl;
  fout<<"\t"<<"TCutG *cutg = new TCutG(\""<<var<<"\","<<N<<");"<<endl;
  fout<<"\t"<<"cutg->SetVarX(\"Phi0\");"<<endl;
  fout<<"\t"<<"cutg->SetVarY(\"CosTheta0\");"<<endl;
  fout<<"\t"<<"cutg->SetTitle(\"TCutG_"<<key<<"\");"<<endl;
  fout<<"\t"<<"cutg->SetFillColor(1);"<<endl;
  fout<<"\t"<<"cutg->SetMarkerStyle(20);"<<endl;

  double x,y;
  for(int i=0;i<N;i++){
    cutg->GetPoint(i,x,y);
    fout<<"\t"<<"cutg->SetPoint("<<i<<","<<x<<","<<y<<");"<<endl;
    
  }

  fout<<"\t"<<"//cutg->Draw(\"AP\");"<<endl;
  fout<<"\t"<<var<<" = cutg;"<<endl;
  fout<<"\t"<<"return cutg;"<<endl;
  
  fout<<"}\n"<<endl;

  //the 2nd routine:  int HMSSolidAngle_g20_p320(double x,double y)
  fout<<"int "<<sub2<<"(double x,double y) {\n" <<endl;
  fout<<"\t"<<"if ( ! "<<var<<" )  "<<sub1<<"();"<<endl;
  fout<<"\t"<<"return "<<var<<"->IsInside(x,y);"<<endl;
  fout<<"}\n"<<endl;

  fout.close();
}


void SolidAcc()
{
   const double deg = atan(1.0)/45.;
   double GammaAngle = GetConfigLeaf("VDAngle")/deg; 
   double ProtonAngle = 0.0;
   if(GammaAngle>180.) ProtonAngle = GetConfigLeaf("SuperBigBiteAngle")/deg;
   else ProtonAngle = GetConfigLeaf("HMSAngle")/deg;

   const char *detector[]={"HMS","SBS"};
   int pSBS = (GammaAngle>180.) ? 1 : 0;
   //double Beam = GetConfigLeaf("Beam"); 


   gStyle->SetPalette(1);  
   TTree *track1 = (TTree*)gROOT->FindObject("track1");
   track1->AddFriend("track0");

   TCanvas *c1 = new TCanvas("c1","Contours",10,10,800,600);
   c1->cd();
   int  min1=1, min2=5;
   if(pSBS==0)  {min1=1;min2=3;}

   gPad->SetRightMargin(0.12);

   char target[255];
   if(pSBS==1)  sprintf(target,"cos(Theta0):Phi0");
   else sprintf(target,"cos(Theta0):fmod(Phi0+2*3.14159,2*3.14159)");

   track1->Draw(Form("%s>>hTPFrame",target),"track1.Pvb>1.0");
   TH2F *hTPFrame=(TH2F*) gROOT->FindObject("hTPFrame");
   if(hTPFrame->GetEntries()<100) return;

   double x1,x2,y1,y2;
   x1=hTPFrame->GetXaxis()->GetXmin();
   x2=hTPFrame->GetXaxis()->GetXmax();
   y1=hTPFrame->GetYaxis()->GetXmin();
   y2=hTPFrame->GetYaxis()->GetXmax();
  
   if(pSBS==1){ 
     if(x1<-0.6) x1=-0.6;
     if(x2> 0.6) x2= 0.6;
   }
   else{
     if(x1<2.8) x1=2.8;
     if(x2>3.5) x2=3.5;
   }
   if(y1<0.0) y1=0.0;
   if(pSBS==1 && y1<cos((ProtonAngle+12)*deg)) y1=cos((ProtonAngle+12)*deg);
   if(y2>1.0) y2=1.0;
   //x1=2.8;x2=3.5;y1=0.6;y2=1.0;
   cout<<"  x1="<<x1<<"  x2="<<x2<<"  y1="<<y1<<"  y2="<<y2<<endl;

   int nbinx=40,nbiny=40;
   if(pSBS==0) nbinx=nbiny=20;
   TH2F *h2tp0, *h2tp1, *h2tp2;
   h2tp0 = new TH2F("h2tp0",
		    Form("%s (Angle=%.0f^{o}); #phi (rad); cos#theta",
			 detector[pSBS],ProtonAngle),
		    nbinx,x1,x2,nbiny,y1,y2);
   h2tp1 = new TH2F("h2tp1",
		    Form("%s acceptance (Angle=%.0f^{o}); #phi (rad); cos#theta",
			 detector[pSBS],ProtonAngle),
		    nbinx,x1,x2,nbiny,y1,y2);
   h2tp2 = new TH2F("h2tp2",
		    Form("%s #gamma-p acceptance (Pr=%.0f^{o} #gamma=%.0f^{o}); #phi (rad); cos#theta",
			 detector[pSBS],ProtonAngle, GammaAngle),
		    nbinx,x1,x2,nbiny,y1,y2);
   //if(pSBS==1) track1->Project("h2tp0",target,"");
   //else track1->Project("h2tp0",target,"track1.Pvb>0");

   track1->Project("h2tp0",target,"");

   track1->Project("h2tp1",target,"track1.Pvb>0");

   //Form("N1_g%.0f_p%.0f",GammaAngle,ProtonAngle)

   cout<<"pSBS="<<pSBS<<" nevent="<<h2tp1->GetEntries()<<endl;
   h2tp1->Divide(h2tp0);h2tp1->Scale(100.);

   //this line will create object contours
   h2tp1->Draw("contztext,list");
   h2tp1->SetMinimum(10);
   
   //we must call Update to force the canvas to be painted.  When 
   //painting the contour plot, the list of contours is generated
   //and a reference to it added to the Root list of special objects
   c1->Update();
   c1->SaveAs(Form("Graph/%sAcc_%.0f.png",detector[pSBS],ProtonAngle));
   h2tp1->SaveAs(Form("Graph/%sAcc_%.0f.root",detector[pSBS],ProtonAngle));

   TCanvas *c2 = new TCanvas("c2","First contour",100,100,800,600);      
   c2->cd();
   TObjArray *contours = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
   if (!contours) return;

   double pArea1,pArea2;
   TCutG *cutg1=0, *cutg2=0;
   TGraph *gc1,*gc2;
   bool valid = false;
   while(!valid)
     {
       cout<<"min1="<<min1<<"  min2="<<min2<<endl;
       TList *lcontour1 = (TList*)contours->At(min1);
       TList *lcontour2 = (TList*)contours->At(min2);
       if (!lcontour1 || ! lcontour2) return;
       gc1 = (TGraph*)lcontour1->First();
       gc2 = (TGraph*)lcontour2->First();
       if (!gc1 || !gc2) return;
       if (gc2->GetN() < 10) return;
       
       gc1->SetTitle(Form("%s acceptance (Angle=%.0f^{o}) ; #phi (rad); cos#theta",detector[pSBS],ProtonAngle));
       gc1->SetMarkerStyle(21);
       gc1->SetMarkerColor(1);
       gc1->Draw("alp");
       gc2->SetMarkerStyle(20);
       gc2->SetMarkerColor(4);
       gc2->Draw("lp same");


       //To get the area, 
       pArea1 = gc1->Integral();
       pArea2 = gc2->Integral();
       cout<<"TGraph::Integral():  Area1="<<pArea1<<"  Area2="<<pArea2<<endl;
       
       //We make a TCutG object with the array obtained from this graph
       //It turns out these 2 results are identical, but TCUTG return negative values 
       if(cutg1) delete cutg1;
       cutg1 = new TCutG("cutg1",gc1->GetN(),gc1->GetX(),gc1->GetY());
       if(cutg2) delete cutg2;
       cutg2 = new TCutG("cutg2",gc2->GetN(),gc2->GetX(),gc2->GetY());
       pArea1=pArea2=0.0;
       pArea1 = -cutg1->Area(); pArea2 = -cutg2->Area();
       cout<<"TCUTG::Area():       Area1="<<pArea1<<"  Area2="<<pArea2<<endl;
       
       if(pArea1<0.001) min1 +=1; 
       if(pArea2<0.001) min2 += 1;
       else valid = true;

       if(min2>20) return;
     }

   cutg2->SaveAs(Form("%sSolidAngleCutG.C",detector[pSBS]));

   Write(cutg2,GammaAngle,ProtonAngle);

   TPaveText *pt1 = new TPaveText(0.40,0.84,0.85,0.89,"brNDC");
   pt1->SetBorderSize(0);
   pt1->SetFillColor(0);
   pt1->AddText(Form("#color[1]{Area1 = %.4f}, #color[4]{Area2 = %.4f}",pArea1,pArea2));
   pt1->Draw("same");
 
   c2->Update();
   c2->SaveAs(Form("Graph/%sSolidAngle_%.0f.png",detector[pSBS],ProtonAngle));
   //c2->SaveAs(Form("Graph/%sSolidAngle_%.0f.C",detector[pSBS],ProtonAngle));

   //Now redraw histo  then create figure
   c1->cd();
   TCutG *myCut = (TCutG*) cutg2->Clone("myCut");
   myCut->SetName("myCut");
   if(pSBS==1)   myCut->SetVarX("Phi0");
   else myCut->SetVarX("fmod(Phi0+2*3.14159,2*3.14159)");
   myCut->SetVarY("cos(Theta0)");
   
   myCut->SetLineWidth(3);myCut->SetLineColor(2);myCut->Draw("same");
   c1->Update();
   c1->SaveAs(Form("Graph/%sAcc_%.0f.png",detector[pSBS],ProtonAngle));
   h2tp1->SaveAs(Form("Graph/%sAcc_%.0f.root",detector[pSBS],ProtonAngle));


   c2->Clear();
   c2->SetRightMargin(0.12);
   track1->Draw(Form("%s>>h2tp0",target),"myCut");
   track1->Draw(Form("%s>>h2tp2",target),"myCut && track1.Pvb>1.0 && track0.Pvb>1.0");
  
   cout<<"pSBS="<<pSBS<<" nevent="<<h2tp2->GetEntries()<<endl;
   h2tp2->Divide(h2tp0);h2tp2->Scale(100.);

   //this line will create object contours
   h2tp2->Draw("contztext");
   myCut->Draw("same"); 
   TPaveText *pt2 = new TPaveText(0.60,0.84,0.85,0.89,"brNDC");
   pt2->SetBorderSize(0);
   pt2->SetFillColor(0);
   pt2->AddText(Form("SolidAngle = %.1f msr",pArea2*1000));
   pt2->Draw("same");
 
   h2tp2->SetMinimum(10);

   c2->Update();
   c2->SaveAs(Form("Graph/%sAcc_p%.0f_g%.0f.png",detector[pSBS],ProtonAngle,GammaAngle));
   h2tp2->SaveAs(Form("Graph/%sAcc_p%.0f_g%.0f.root",detector[pSBS],ProtonAngle,GammaAngle));

   TCanvas *c3 = new TCanvas("c3","Acceptance",200,200,800,600);  
   c3->SetRightMargin(0.12);
   c3->cd();
   TH2F *h2tp0_tr, *h2tp1_tr;
   if(pSBS==1) {
     nbinx=30;x1=-0.15;x2=0.15;
     nbiny=30;y1=-0.40;y2=0.20;  
   }
   else{
     nbinx=14;x1=-0.07;x2=0.07;
     nbiny=30;y1=-0.10;y2=0.20;  
   }
   h2tp0_tr = new TH2F("h2tp0_tr",
		       Form("%s (Angle=%.0f^{o}); #phi_{tr} (rad); #theta_{tr} (rad)",
			    detector[pSBS],ProtonAngle),
		       nbinx,x1,x2,nbiny,y1,y2);
   h2tp1_tr = new TH2F("h2tp1_tr",
		       Form("%s #gamma-p acceptance (Pr=%.0f^{o}, #gamma=%.0f); #phi_{tr} (rad); #theta_{tr} (rad)",
			    detector[pSBS],ProtonAngle,GammaAngle),
		       nbinx,x1,x2,nbiny,y1,y2);

   track1->Draw("Theta0_tr:Phi0_tr>>h2tp0_tr","myCut");
   track1->Draw("Theta0_tr:Phi0_tr>>h2tp1_tr","myCut && track1.Pvb>1.0 && track0.Pvb>1.0");
  
   cout<<"pSBS="<<pSBS<<" nevent="<<h2tp1_tr->GetEntries()<<endl;
   h2tp1_tr->Divide(h2tp0_tr);h2tp1_tr->Scale(100.);

   //this line will create object contours
   h2tp1_tr->Draw("contztext");
   pt2->Draw("same");
 
   h2tp1_tr->SetMinimum(10);

   c3->Update();
   c3->SaveAs(Form("Graph/%sAcc_tr_p%.0f_g%.0f.png",detector[pSBS],ProtonAngle,GammaAngle));
   h2tp1_tr->SaveAs(Form("Graph/%sAcc_tr_p%.0f_g%.0f.root",detector[pSBS],ProtonAngle,GammaAngle));

}
