//20080324  By Jixie Zhang
//Fit a TH2 to get the upper and lower boundaries then integral it to
//get the area

#define FITTH2_DEBUG 1

#include <math.h>
#include <string>
#include <time.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>

#include <TROOT.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TObject.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TProfile.h>
#include <TLine.h>
#include <TStyle.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TCutG.h>
#include <TAxis.h>
#include <TDirectory.h>
#include <TPaveText.h>
using namespace std;


void GetTH1BoundaryValues(TH1 *h1,int &pXBinStart,int &pXBinEnd, double &pYStart, 
			  double &pYEnd, double entriescut=-1)
{
  // Get the start bin and end bin 
  int pNXBin=h1->GetNbinsX();
  pXBinStart=0, pXBinEnd=-1;

  if(entriescut<=0) entriescut=0.1*h1->GetMaximum(); 
  for(int i=1;i<=pNXBin;i++)
    {
      //continuous 4 bins none zero
      if(pXBinStart<=0 && i+3<=pNXBin)
	{
	  if (h1->GetBinContent(i)>entriescut   && h1->GetBinContent(i+1)>0 &&
	      h1->GetBinContent(i+2)>0 && h1->GetBinContent(i+3)>0 )
	    pXBinStart=i;
	}  

      if(pXBinEnd<=0 && pNXBin-i-2>=1)
	{
	  if (h1->GetBinContent(pNXBin-i+1)>entriescut && h1->GetBinContent(pNXBin-i)>0 &&
	      h1->GetBinContent(pNXBin-i-1)>0 && h1->GetBinContent(pNXBin-i-2)>0)
	    pXBinEnd=pNXBin-i+1;
	}
      if(pXBinStart>0 && pXBinEnd>0) break;
    }
  if(pXBinEnd<=0) pXBinEnd=pNXBin;

  pYStart=h1->GetBinLowEdge(pXBinStart);
  pYEnd=h1->GetBinLowEdge(pXBinEnd)+h1->GetBinWidth(pXBinEnd);

  double pYmax=h1->GetMaximum();	
  TLine *L1=new TLine(pYStart,0,pYStart,pYmax);
  L1->SetLineColor(2);
  h1->GetListOfFunctions()->Add(L1); 
  TLine *L2=new TLine(pYEnd,0,pYEnd,pYmax);
  h1->GetListOfFunctions()->Add(L2); 
  L2->SetLineColor(2);

}
// arguments: 
//    TH2 *h2:              2-D histo;
//    int entriescut=10:    the minimum entries of the 1-D y histo, if less than this 
//                          number,just ignore this group
//                          if entriescut<=0, ignore this cut;
//    int firstbin=0:       the 1st bin of the fitted region
//    int lastbin=-1:       the last bin of the fitted region
//    bool  save1dhisto=false:      whether to write the 1-D temparaly histos 
//    int nbpg=1; //number of bins per group, just like rebin
//
//return value: ngps, number of good groups which has a peak value
//              -2  , too few entries (<100) in this 2-D h2
//              -1  , too few groups (<3)
double GetTH2Area(TH2 *h2,double entriescut=0, int firstbin=0, int lastbin=-1, 
		  bool save1dhisto=false, int nbpg=1)
{

  gStyle->SetOptFit(0);
  gStyle->SetPadRightMargin(0.15);
  if (!h2) return -3;
  if(h2->GetEntries()<100) return -2;

  char name[100],title[200];
  // Get the start bin and end bin 
  int pNXBin=h2->GetNbinsX();
  //int pNYBin=h2->GetNbinsY();  //not used
  int pXBinStart=0, pXBinEnd=-1;
  double pXStart,pXEnd;

  system("mkdir tmp_FitTH2");
  sprintf(name,"%s_prjx",h2->GetName());
  TH1D *h1prjx=h2->ProjectionX(name);
  GetTH1BoundaryValues((TH1*)h1prjx,pXBinStart,pXBinEnd,pXStart,pXEnd,entriescut);
	

  //reset the fitted interval:   firstbin and   lastbin  
#ifdef FITTH2_DEBUG
  cout<<"pXBinStart="<<pXBinStart<<"  pXBinEnd="<<pXBinEnd<<endl;
#endif
  if( firstbin<pXBinStart || firstbin>pXBinEnd) firstbin=pXBinStart;
  if( lastbin<pXBinStart || lastbin>pXBinEnd) lastbin=pXBinEnd;

  int pNGroup=(int)ceil(double(lastbin-firstbin)/double(nbpg));
#ifdef FITTH2_DEBUG
  cout<<"firstbin="<<firstbin<<"  lastbin="<<lastbin<<" pNGroup="<<pNGroup<<endl;
#endif
  if(firstbin<0 || lastbin>pNXBin || pNGroup<3) return -1;

  //now we can start to get the min and max 
  TH1D **h1prjy_g;
  h1prjy_g=new TH1D* [pNGroup];
  double *pXcenter,*pYmin,*pYmax,*pDeltaY;
  pXcenter=new double[pNGroup];
  pYmin=new double[pNGroup];
  pYmax=new double[pNGroup];
  pDeltaY=new double[pNGroup];

  //create the canvas
  sprintf(title,"%s",h2->GetYaxis()->GetTitle());
  int ncol=4,nrow=ceil(double(pNGroup)/double(ncol));
  int canw=256*ncol,canh=200*nrow;
  TCanvas *c1=new TCanvas("c1",title,canw,canh);
  c1->Divide(ncol,nrow,0.001,0.001);

  //Get the Y min and Y max
  double XBinW=h1prjx->GetBinWidth(1);
  TH1D *h1=0;
  int pNpt=0;
  double IntegralArea = 0.0;
  for(int i=0;i<pNGroup;i++)
    {
      c1->cd(i+1);
      //get the 1-D y distribution for this group
      sprintf(name,"%s_prjy_group%02d",h2->GetName(),i);
      double Xmin=h1prjx->GetBinLowEdge(firstbin+i*nbpg);
      double Xmax=Xmin+XBinW*nbpg;
      sprintf(title,"%f<=%s<%f",Xmin,h2->GetXaxis()->GetTitle(),Xmax);
      h1=h2->ProjectionY(name,firstbin+i*nbpg,firstbin+(i+1)*nbpg);
      h1->SetTitle(title);      
      h1->SetXTitle(h2->GetYaxis()->GetTitle());

      if(h1->GetEntries()>2.0*entriescut)
	{
	  GetTH1BoundaryValues((TH1*)h1,pXBinStart,pXBinEnd,pYmin[pNpt],pYmax[pNpt],entriescut);
	  pXcenter[pNpt]=(Xmax+Xmin)/2.0;
	  pDeltaY[pNpt]=pYmax[pNpt]-pYmin[pNpt];
	  pNpt++;
	  //calculate area
	  IntegralArea += pDeltaY[pNpt] * XBinW*nbpg;
	}
      h1->Draw();
      h1prjy_g[i]=h1;
    }  
  c1->cd();
  sprintf(name,"tmp_FitTH2/%s_1Ddistr.png",h2->GetName());
  c1->SaveAs(name);
  sprintf(name,"tmp_FitTH2/%s_1Ddistr.root",h2->GetName());
  c1->SaveAs(name);

  cout<<"\nIntegralArea = "<<IntegralArea<<endl;
  ///////////////////////////////////
  if(save1dhisto)
    {
      sprintf(name,"%s_1Ddistr",h2->GetName());
      c1->Write(name);
    }
  else
    { //remove these temperate 1-D histos
      // for(int i=0;i<pNGroup;i++)
      // {
      // if(h1prjy_g[i]) h1prjy_g[i]->Delete();
      // }
      ;
    }
  ///////////////////////////////////


  //create the graph
  TCanvas *c2gr=new TCanvas("c2gr","boundary on top",800,600);
  //c2gr->Divide(1,2,0.001,0.001);

  TGraph *grlow=new TGraph(pNpt,pXcenter,pYmin);    
  grlow->SetMarkerStyle(22);
  grlow->SetMarkerColor(1); 
  TGraph *grhigh=new TGraph(pNpt,pXcenter,pYmax);    
  grhigh->SetMarkerStyle(23);
  grhigh->SetMarkerColor(1); 

  TF1 *f1=new TF1("f1","[0]+[1]*x+[2]*x*x+[3]*pow(x,3.0)+[4]*pow(x,4.0)+[5]*pow(x,5.0)",
		  pXcenter[0],pXcenter[pNpt-1]);
  TF1 *f2=new TF1("f2","[0]+[1]*x+[2]*x*x+[3]*pow(x,3.0)+[4]*pow(x,4.0)+[5]*pow(x,5.0)",
		  pXcenter[0],pXcenter[pNpt-1]);
  f1->SetLineColor(1);
  f2->SetLineColor(1);

  c2gr->cd();
  grlow->Draw("AP");
  grlow->Fit(f1);
  c2gr->SaveAs(Form("tmp_FitTH2/%s_LowBoundaryGraph.png",h2->GetName()));

  c2gr->Clear(); 
  grhigh->Draw("AP");
  grhigh->Fit(f2);
  c2gr->SaveAs(Form("tmp_FitTH2/%s_HighBoundaryGraph.png",h2->GetName()));

  double pArea=f2->Integral(pXcenter[0],pXcenter[pNpt-1]) - f1->Integral(pXcenter[0],pXcenter[pNpt-1]);
  //double pArea=f2->Integral(pXStart,pXEnd) - f1->Integral(pXStart,pXEnd);

  c2gr->Clear(); 
  h2->Draw("contz");
  grlow->Draw("Psame");
  grhigh->Draw("Psame");
	
  TLine *L1=new TLine(pXcenter[0],pYmin[0],pXcenter[0],pYmax[0]);
  TLine *L2=new TLine(pXcenter[pNpt-1],pYmin[pNpt-1],pXcenter[pNpt-1],pYmax[pNpt-1]);
  L1->SetLineWidth(3);
  L2->SetLineWidth(3);
  L1->Draw("same");
  L2->Draw("same");
  TPaveText *pt2 = new TPaveText(0.7,0.82,0.89,0.89,"brNDC");
  pt2->SetBorderSize(0);
  pt2->SetFillColor(0);
  pt2->AddText(Form("#color[4]{Area = %g}",pArea));
  pt2->AddText(Form("#color[4]{IntegralArea = %g}",IntegralArea));
  pt2->Draw();

  c2gr->cd(); 
  c2gr->SaveAs(Form("tmp_FitTH2/%s_BoundaryGraph.png",h2->GetName()));

  //free the memory
  delete [] pXcenter;
  delete [] pYmin;
  delete [] pYmax;
  delete [] pDeltaY;

  return pArea;
}


////////////////////////////////////////////////////////////////////////////////
//  this method requires the h2 have very fine bins, otherwise it will give very 
//  large uncertainty as large as = (Nx+NY)/ Nx*Ny = 1/Nx + 1/Ny
// arguments: 
//    TH2 *h2:              2-D histo;
//    TF1 *fcn:             the fitting function
//    int firstbin=0:       the 1st bin of the fitted region
//    int lastbin=-1:       the last bin of the fitted region
//    int entriescut:       number of entries required in each slice 
//    Option_t *option;     Fitting option, "QNRG2|3|4|5", by default is "QNRG4"
//
////////////////////////////////////////////////////////////////////////////////
double GetTH2AreaW(TH2 *h2, double entriescut=-1, int WeightByZ=0)
{
  system("mkdir -p tmp_FitTH2");
  gStyle->SetOptFit(0);
  gStyle->SetPadRightMargin(0.15);
  if(entriescut<0) entriescut = 0.1 * h2->GetMaximum();
  double pGridArea = h2->GetXaxis()->GetBinWidth(1) * h2->GetYaxis()->GetBinWidth(1); 
  int Nx=h2->GetNbinsX(),Ny=h2->GetNbinsY();
  double pArea=0.0;
  for(int i=1;i<=Nx;i++)
    {
      for(int j=1;j<=Ny;j++)
	{
	  double pGridEntries = h2->GetBinContent(i*Nx+j);
	  if(pGridEntries>entriescut) 
	    {
	      if(WeightByZ) pArea += pGridArea*pGridEntries;
	      else pArea += pGridArea;
	    }
	}
    }

  //create the graph
  TCanvas *c1=new TCanvas("c1","",800,600);
  c1->cd();
  h2->Draw("contz");
  TPaveText *pt2 = new TPaveText(0.50,0.83,0.85,0.89,"brNDC");
  pt2->SetBorderSize(0);
  pt2->SetFillColor(0);
  if(WeightByZ) pt2->AddText(Form("#color[4]{Weighted Area = %f}",pArea));
  else pt2->AddText(Form("#color[4]{Area = %f}",pArea));
  pt2->Draw();


  c1->cd(); 
  if(WeightByZ) c1->SaveAs(Form("tmp_FitTH2/%s_Area_WeightByZ.png",h2->GetName()));
  else c1->SaveAs(Form("tmp_FitTH2/%s_Area.png",h2->GetName()));
  return pArea;
}

//to integrate the XS, the z must be averaged XS, not accumulated XS
void Inte_dSigmadOmega_SBS()
{
  double Lumi_NH3 = 2.8 / 100. ;  //L_g2p_N14 = 2.8E34,  in unit of 10^36
  double Lumi_He  = 3.6 / 100. ;  //L_g2p_He4 = 3.6E34,  in unit of 10^36

  TTree* gp = (TTree*)gROOT->FindObject("gp");
  gStyle->SetOptFit(0);
  gStyle->SetPadRightMargin(0.15);

  gp->Draw("P0_g*P0_g/3.142*XS*1000>>hW","XS>0","");
  TH1F* hW = (TH1F*)gROOT->FindObject("hW");
  double MeanW = hW->GetMean();   //the weighting factor

  gp->Draw("Theta0_g:fmod(Phi0_g+2*3.1416,2*3.1416)>>h2PhiTheta(30,2.4,3.9,35,0.2,0.9)","XS>0","colz");
  TH2F* h2D = (TH2F*)gROOT->FindObject("h2PhiTheta");
  gp->Draw("Theta0_g:fmod(Phi0_g+2*3.1416,2*3.1416)>>h2PhiThetaXS(30,2.4,3.9,35,0.2,0.9)","P0_g*P0_g/3.142*XS*1000*(XS>0)","colz");
  TH2F* h2N = (TH2F*)gROOT->FindObject("h2PhiThetaXS");
  TH2F* h2MeanXS = (TH2F*)h2N->Clone("h2MeanXS");
  h2MeanXS->Divide(h2D);
  cout<<"Number of Entries = "<<h2D->GetEntries()<<endl;
  h2MeanXS->SetTitle("SBS+LAC RCS Acceptance: XS Weighted (pb); #phi_{#gamma} (rad) ; #theta_{#gamma} (rad) ");
  h2MeanXS->Draw("colz");
  //.L GetTH2Area.cc
  //GetTH2Area(h2MeanXS,0.1,0,-1,true);
  double SolidAngle = 0;
  SolidAngle = GetTH2AreaW(h2MeanXS,0.2*2,0);  
  cout<<"MeanW = "<< MeanW <<" SolidAngle = "<< SolidAngle <<endl;
  double rate_H_mean = SolidAngle * MeanW * 3 * Lumi_NH3;
  double rate_All_mean = SolidAngle *MeanW * ( 17 * Lumi_NH3 + 4 * Lumi_He);
  cout<<"Rate_h_mean = "<< rate_H_mean <<" Hz  Rate_all_mean = "<< rate_All_mean<<" Hz"<<endl;

  double Inte_dSigmadOmega = 0 ;
  Inte_dSigmadOmega = GetTH2AreaW(h2MeanXS,0.2*2,1);   //in unit of pb, or 10^-36 cm^2

  double rate_H = Inte_dSigmadOmega * 3 * Lumi_NH3;
  double rate_All = Inte_dSigmadOmega * ( 17 * Lumi_NH3 + 4 * Lumi_He);
  cout<<"Rate_h = "<< rate_H <<" Hz  Rate_all = "<< rate_All<<" Hz"<<endl;
}

void Inte_dSigmadOmega_HMS()
{
  double Lumi_NH3 = 2.8 / 100. ;  //L_g2p_N14 = 2.8E34,  in unit of 10^36
  double Lumi_He  = 3.6 / 100. ;  //L_g2p_He4 = 3.6E34,  in unit of 10^36

  TTree* gp = (TTree*)gROOT->FindObject("gp");
  gStyle->SetOptFit(0);
  gStyle->SetPadRightMargin(0.15);

  gp->Draw("P0_g*P0_g/3.142*XS*1000>>hW","XS>0","");
  TH1F* hW = (TH1F*)gROOT->FindObject("hW");
  double MeanW = hW->GetMean();   //the weighting factor

  gp->Draw("Theta0_g:Phi0_g>>h2PhiTheta(40,-0.16,0.34,40,0.2,0.6)","XS>0","colz");
  TH2F* h2D = (TH2F*)gROOT->FindObject("h2PhiTheta");
  gp->Draw("Theta0_g:Phi0_g>>h2PhiThetaXS(40,-0.16,0.34,40,0.2,0.6)","P0_g*P0_g/3.142*XS*1000*(XS>0)","colz");
  TH2F* h2N = (TH2F*)gROOT->FindObject("h2PhiThetaXS");
  TH2F* h2MeanXS = (TH2F*)h2N->Clone("h2MeanXS");
  h2MeanXS->Divide(h2D);
  cout<<"Number of Entries = "<<h2D->GetEntries()<<endl;
  h2MeanXS->SetTitle("SBS+LAC RCS Acceptance: XS Weighted (pb); #phi_{#gamma} (rad) ; #theta_{#gamma} (rad) ");
  h2MeanXS->Draw("colz");
  //.L GetTH2Area.cc
  //GetTH2Area(h2MeanXS,0.1,0,-1,true);
  double SolidAngle = 0;
  SolidAngle = GetTH2AreaW(h2MeanXS,13.0*2,0);
  cout<<"MeanW = "<< MeanW <<" SolidAngle = "<< SolidAngle <<endl;
  double rate_H_mean = SolidAngle * MeanW * 3 * Lumi_NH3;
  double rate_All_mean = SolidAngle *MeanW * ( 17 * Lumi_NH3 + 4 * Lumi_He);
  cout<<"Rate_h_mean = "<< rate_H_mean <<" Hz  Rate_all_mean = "<< rate_All_mean<<" Hz"<<endl;

  double  Inte_dSigmadOmega = 0;
  Inte_dSigmadOmega = GetTH2AreaW(h2MeanXS,13.0*2,1); //in unit of pb, or 10^-36 cm^2

  double rate_H = Inte_dSigmadOmega * 3 * Lumi_NH3;
  double rate_All = Inte_dSigmadOmega * ( 17 * Lumi_NH3 + 4 * Lumi_He);
  cout<<"Rate_h = "<< rate_H <<" Hz  Rate_all = "<< rate_All<<" Hz"<<endl;
}

//return integral BremXS :  #inte{dSigmadK * dK}
double BremXS(double Beam, double E_g_min, double E_g_max)
{
  //from PDG (2010)  EQ 27.27
  double A = 64; //copper atomic number
  double X0= 1.43993; //radlength of copper, in cm
  double Na = 6.02E23; //avogadro number
  double k = 0.0, y = 0.0;   // y=E_g/Beam
  double dSigmadK = 0;
  double Inte_dSigmadK = 0;
  double barn = 1.0E-028;  // in unit of cm^2

  int Nbin = 1000;
  double dK = (E_g_max-E_g_min)/Nbin;
  for(int i=0;i<Nbin;i++)
    {
      k = E_g_min + i*dK;
      y = k / Beam;
      dSigmadK = A/(X0*Na*k)*(4.0/3.0-4.0*y/3.0+y*y);
      Inte_dSigmadK += dSigmadK * dK;
    }
  cout<<"BremXS(Beam="<<Beam<<", Kmin="<<E_g_min<<", Kmax="<< E_g_max<<") = "
      <<Inte_dSigmadK<<"\n";

  double Normalization = 0;
  Nbin = int((Beam-0.00001)/dK);
  for(int i=0;i<Nbin;i++)
    {
      k = 0.00001+double(i)*dK;
      y = k / Beam;
      dSigmadK = A/(X0*Na*k)*(4.0/3.0-4.0*y/3.0+y*y);
      Normalization += dSigmadK * dK;
    }
  cout<<"Normailized BremXS(Beam="<<Beam<<", Kmin="<<E_g_min<<", Kmax="<< E_g_max<<") = "
      <<Inte_dSigmadK/Normalization<<"\n";

  return Inte_dSigmadK/barn;
}


double GammaFlux(double Thickness_in_radlen, double Beam, double E_g_min, double E_g_max)
{
  double Term1 = 4.*log(E_g_max/E_g_min)/3.;
  double Term2 = 4.*(E_g_max-E_g_min)/3./Beam;
  double Term3 = (E_g_max*E_g_max-E_g_min*E_g_min)/2./Beam/Beam;
  double NGamma = Thickness_in_radlen * (Term1-Term2+Term3);
  return NGamma;
}


//20140128
//Jixie: I found a very good way to get the area, just use class TCUTG::Area
//or just use TGraph::Integral()
//see root tutorial: $ROOTSYS/tutorial/hist/FirstContour.C  for details


void Write(TCutG *cutg, double GammaAngle=330, double ProtonAngle=30)
{
  ofstream fout;
  char key[255], file[255];
  sprintf(key,"SBSSolidAngle_g%.0f_p%.0f",GammaAngle,ProtonAngle);
  sprintf(file,"%s.Ccc,key);
  
  fout.open(file);

  int N = cutg->GetN();

  fout<<"#include <TROOT.h>"<<endl;
  fout<<"#include <TCutG.h>"<<endl;

  fout<<"int "<<key<<"(double x,double y) {\n" <<endl;
  fout<<"\t"<<"TCutG *cutg = new TCutG(\"cutg2\","<<N<<");"<<endl;
  fout<<"\t"<<"cutg->SetVarX(\"Phi0\");"<<endl;
  fout<<"\t"<<"cutg->SetVarY(\"CosTheta0\");"<<endl;
  fout<<"\t"<<"cutg->SetTitle(\""<<key<<"\");"<<endl;
  fout<<"\t"<<"cutg->SetFillColor(1);"<<endl;
  fout<<"\t"<<"cutg->SetMarkerStyle(20);"<<endl;

  double x,y;
  for(int i=0;i<N;i++){
    cutg->GetPoint(i,x,y);
    fout<<"\t"<<"cutg->SetPoint("<<i<<","<<x<<","<<y<<");"<<endl;
    
  }

  fout<<"\t"<<"//cutg->Draw(\"AP\");"<<endl;
  fout<<"\t"<<"delete cutg;"<<endl;
  fout<<"\t"<<"return cutg->IsInside(x,y);"<<endl;
  
  fout<<"}\n"<<endl;
  
  fout.close();
}


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
   
  static char key[255];  
  sprintf(key,"g%.0f_pr%.0f_E%.1f",GammaAngle/deg,ProtonAngle/deg,Beam);
  cout<<"key="<<key<<endl;
  return key;
}

void Area(int weight=0)
{
  int SetupHMS=0;
  double Beam, ProtonAngle, GammaAngle;
  const double deg = atan(1.0)/45.;
  GetKeyWords(SetupHMS, Beam, ProtonAngle, GammaAngle);  


  gStyle->SetPalette(1);  
  TTree *gp = (TTree*)gROOT->FindObject("gp");
  TCanvas *c1 = new TCanvas("c1","Contours",10,10,800,600);
  TCanvas *c2 = new TCanvas("c2","First contour",100,100,800,600);

  c1->cd();
  int  min1=0, min2=2;

  gPad->SetRightMargin(0.12);
  //this line will create object contours
  if(!weight) gp->Draw("cos(Theta0_p):Phi0_p","", "contz,list");
  else gp->Draw("cos(Theta0_p):Phi0_p","XS*P0_g*P0_g*1000", "contz,list");
 
  //gp->Draw("Theta0_tr_p:Phi0_tr_p","", "contz,list");

  cout<<"min1="<<min1<<"  min2="<<min2<<endl;
   
  //we must call Update to force the canvas to be painted.  When 
  //painting the contour plot, the list of contours is generated
  //and a reference to it added to the Root list of special objects
  c1->Update();
      
  c2->cd();
  TObjArray *contours = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
  if (!contours) return;
  TList *lcontour1 = (TList*)contours->At(min1);
  TList *lcontour2 = (TList*)contours->At(min2);
  if (!lcontour1 || ! lcontour2) return;
  TGraph *gc1 = (TGraph*)lcontour1->First();
  TGraph *gc2 = (TGraph*)lcontour2->First();
  if (!gc1 || !gc2) return;
  if (gc2->GetN() < 10) return;

  gc1->SetTitle("SBS acceptance; #phi (rad); cos#theta");
  //gc1->SetTitle("SBS acceptance; #phi_{0}^{tr} (rad); #theta_{0}^{tr} (rad)");
  gc1->SetMarkerStyle(21);
  gc1->SetMarkerColor(1);
  gc1->Draw("alp");
  gc2->SetMarkerStyle(20);
  gc2->SetMarkerColor(4);
  gc2->Draw("lp same");

  //To get the area, 
  double pArea1 = gc1->Integral();
  double pArea2 = gc2->Integral();
  cout<<"TGraph::Integral():  Area1="<<pArea1<<"  Area2="<<pArea2<<endl;


  //We make a TCutG object with the array obtained from this graph
  //It turns out these 2 results are identical, but TCUTG return negative values 
  TCutG *cutg1 = new TCutG("cutg1",gc1->GetN(),gc1->GetX(),gc1->GetY());
  TCutG *cutg2 = new TCutG("cutg2",gc2->GetN(),gc2->GetX(),gc2->GetY());
  pArea1=pArea2=0.0;
  //pArea1 = cutg1->Integral(); pArea2 = cutg2->Integral();
  pArea1 = -cutg1->Area(); pArea2 = -cutg2->Area();
  cout<<"TCUTG::Area():       Area1="<<pArea1<<"  Area2="<<pArea2<<endl;
  //cutg2->SaveAs("SBSSolidAngleTCutG.C");  

  Write(cutg2,GammaAngle/deg,ProtonAngle/deg);

  TPaveText *pt2 = new TPaveText(0.50,0.83,0.89,0.89,"brNDC");
  pt2->SetBorderSize(0);
  pt2->SetFillColor(0);
  pt2->AddText(Form("#color[1]{Area1 = %.4f}, #color[4]{Area2 = %.4f}",pArea1,pArea2));
  pt2->Draw("same");
 
  /*
   
  //We create a polymarker object with npmax points.
  const Int_t npmax = 5000;
  TPolyMarker *pm1 = new TPolyMarker(npmax);
  Int_t np = 0;
  while(1) {
  Double_t x = -4 +8*gRandom->Rndm();
  Double_t y = -4 +8*gRandom->Rndm();
  if (cutg1->IsInside(x,y)) {
  pm1->SetPoint(np,x,y);
  np++;
  if (np >= npmax) break;
  }
  }
  pm1->Draw("same");    

  
   
  //We create a polymarker object with npmax points.
  TPolyMarker *pm2 = new TPolyMarker(npmax);
  np = 0;
  while(1) {
  Double_t x = -4 +8*gRandom->Rndm();
  Double_t y = -4 +8*gRandom->Rndm();
  if (cutg2->IsInside(x,y)) {
  pm2->SetPoint(np,x,y);
  np++;
  if (np >= npmax) break;
  }
  }
  pm2->Draw("same");    
  */
  c2->Update();

  c2->SaveAs("Graph/SBSSolidAngle.png");
  //c2->SaveAs("Graph/SBSSolidAngle.cc");  

  //c1->cd();
  //cutg1->Draw("p same");
  //c1->Update();
}
