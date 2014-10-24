//this code is used to plot RCS physics
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

double GetRCSXS(double Ei, double Theta_cm_g);

#ifndef WIN32
#include "RCSXS.cc"
#endif

//input: E0,Theta_lab_g in GeV and rad
//output: E_g,P_p,Theat_p,Theta_cm_g, in GeV and Rad, and dSigmaOverdt in nb/GeV2/Sr
void RCSPhys(double E0, double Theta_lab_g, double& E_g, double& P_p, 
	     double& Theat_lab_p, double& Theta_cm_g, double& dSigmaOverdt)
{
	const double kMassPr = 0.9383;

	double cosTh_g = cos(Theta_lab_g);
	E_g = E0 / (1.0 + E0/kMassPr * (1.0-cosTh_g) );

	TVector3 V3G;
	V3G.SetMagThetaPhi(E_g,Theta_lab_g,0.0); 
	TVector3 V3P(-V3G.X(),-V3G.Y(),E0-V3G.Z());
	P_p = V3P.Mag();
	Theat_lab_p = V3P.Theta();

	TLorentzVector V4Beam(0,0,E0,E0);
	TLorentzVector V4Tg(0,0,0,kMassPr);

	TLorentzVector V4Beta = V4Beam + V4Tg;
	TVector3 V3Beta(0,0,-V4Beta.P()/V4Beta.E());
	TLorentzVector V4G_cm;
	V4G_cm.SetVectM(V3G,0.0);
	V4G_cm.Boost(V3Beta);

	Theta_cm_g = V4G_cm.Theta();
	dSigmaOverdt = GetRCSXS(E0,Theta_cm_g);
}



void testRCSPhys()
{
	//input: E0,Theta_lab_g in GeV and rad
	//output: E_g,P_p,Theat_p,Theta_cm_g, in GeV and Rad, and dSigmaOverdt in nb/GeV2/Sr
	//extern void RCSPhys(double E0, double Theta_lab_g, double& E_g, double& P_p, 
	//	double& Theat_lab_p, double& Theta_cm_g, double& dSigmaOverdt);

	const int N = 50;
	double E0=11.0,Theta_lab_g[N]; //in GeV and rad
	double E_g[N],P_p[N],Theta_lab_p[N],Theta_cm_g[N],dSigmaOverdt[N];

	for(int i=0;i<N;i++)
	{
		Theta_lab_g[i] = (3.0 + 3.0 * i)/57.3;
		RCSPhys(E0,Theta_lab_g[i],
			E_g[i],P_p[i],Theta_lab_p[i],Theta_cm_g[i],dSigmaOverdt[i]);
		Theta_lab_g[i]*=57.3;
		Theta_lab_p[i]*=57.3;
		Theta_cm_g[i] *=57.3;
		dSigmaOverdt[i]*=1000.0;
	}

	double *xx[]={Theta_lab_g, Theta_lab_g, Theta_lab_g,Theta_lab_p};	
	double *yy[]={E_g,Theta_lab_p,dSigmaOverdt,P_p};

	const char *Title[]={
		"E = 11.0 GeV; #theta_{#gamma} (deg); E_{#gamma}",
		"E = 11.0 GeV; #theta_{#gamma} (deg); #theta_{p}",
		"E = 11.0 GeV; #theta_{#gamma} (deg); d#sigma/dt (pb/GeV^{2})",
		"E = 11.0 GeV; #theta_{p} (deg); P_{p} (GeV)"
	};
	const char *Name[] ={"EgTheatg_E11.png","TpThetag_E11.png","XSThetag_E11.png","PpThetap_E11.png"}; 
	TGraph *Gr=0; 

	TCanvas *cg= new TCanvas("cg","",800,600);

	for(int ig=0;ig<4;ig++)
	{
		cg->Clear();
		cg->cd();
		Gr = new TGraph(N,xx[ig],yy[ig]);
		Gr->SetTitle(Title[ig]);
		Gr->Draw("AL");
		cg->SaveAs(Name[ig]);
	}

}

double GammaFlux(double Thickness_in_radlen, double Beam, double E_g_min, double E_g_max)
{
  double Term1 = 4.*log(E_g_max/E_g_min)/3.;
  double Term2 = 4.*(E_g_max-E_g_min)/3./Beam;
  double Term3 = (E_g_max*E_g_max-E_g_min*E_g_min)/2./Beam/Beam;
  double NGamma = Thickness_in_radlen * (Term1-Term2+Term3);
  return NGamma;
}

//this routine calculate rates for a given incident energy
//Assuming the given theta_g is the average value for the whole photon arm, 
//calculate RCS rates for gamma-pr elas events   
void testRCSPhys(double Beam)
{
	//input: E0,Theta_lab_g in GeV and rad
	//output: E_g,P_p,Theat_p,Theta_cm_g, in GeV and Rad, and dSigmaOverdt in nb/GeV2/Sr
	//extern void RCSPhys(double E0, double Theta_lab_g, double& E_g, double& P_p, 
	//	double& Theat_lab_p, double& Theta_cm_g, double& dSigmaOverdt);

  double Lumi_NH3 = 2.8 / 100. ;  //L_g2p_N14 = 2.8E34,  in unit of 10^36
  double Lumi_He  = 3.6 / 100. ;  //L_g2p_He4 = 3.6E34,  in unit of 10^36
  double SolidAngle_gp_HMS = 0.01875;  //msr
  double SolidAngle_gp_SBS = 0.35;     //msr, change from 0.25 to 0.40
  double GF = GammaFlux(0.06,Beam,0.3*Beam,0.999*Beam);
  //for 4.3 GeV, I want to caculate the rate for E99114
  if(fabs(Beam-4.3)<0.001)   
    {
      SolidAngle_gp_SBS = SolidAngle_gp_HMS;
      Lumi_NH3 = 120./3.;  //luminosity of E99114 = 1.2E38
      GF = GammaFlux(0.06,Beam,0.3*Beam,0.999*Beam);
    }
        const int N = 58;
	double E0=Beam,Theta_lab_g[N]; //in GeV and rad
	double E_g[N],P_p[N],Theta_lab_p[N],Theta_cm_g[N],dSigmaOverdt[N];
	double t[N], rate_h[N], rate_all[N];

	ofstream fout;
	fout.open(Form("rcs_phys_E%.3f.txt",E0));
	char str[255];
	sprintf(str,"  Ei_g   Ef_g  Theta_g    P_p  Theta_p Theta_cm_g    t  dSigma(pb)  rate_h rate_all(hz)\n");
	cout<<str;fout<<str;
	for(int i=0;i<N;i++)
	{
		Theta_lab_g[i] = (6.0 + 2.0 * i)/57.3;
		RCSPhys(E0,Theta_lab_g[i],
			E_g[i],P_p[i],Theta_lab_p[i],Theta_cm_g[i],dSigmaOverdt[i]);
		t[i] = -2 * E0 * E_g[i] * ( 1.0 - cos(Theta_lab_g[i]) );
		Theta_lab_g[i]*=57.3;
		Theta_lab_p[i]*=57.3;
		Theta_cm_g[i] *=57.3;
		dSigmaOverdt[i] *= 1000.0;  //turn into pb
		rate_h[i] = dSigmaOverdt[i] * (E_g[i]*E_g[i]/3.1416) * SolidAngle_gp_SBS * GF * (3*Lumi_NH3);
		rate_all[i] = dSigmaOverdt[i] * (E_g[i]*E_g[i]/3.1416) * SolidAngle_gp_SBS *
		  GF *(17*Lumi_NH3+4*Lumi_He);
		if( (Theta_lab_g[i]<54) && 
		    (Theta_lab_p[i]<54 || fabs(Theta_lab_p[i]-90)<18) )
		  {
		    sprintf(str,"%6.3f %6.3f %8.1f %6.3f %8.1f %8.1f %8.3f %8.3f %8.3f %8.3f\n",
			    E0,E_g[i],Theta_lab_g[i],P_p[i],Theta_lab_p[i],Theta_cm_g[i],t[i],
			    dSigmaOverdt[i],rate_h[i],rate_all[i]);
		    cout<<str;fout<<str;
		  }
	}
	fout.close();

	double *xx[]={Theta_lab_g, Theta_lab_g, Theta_lab_g, Theta_lab_p, Theta_lab_g};	
	double *yy[]={E_g, Theta_lab_p, Theta_cm_g, P_p, dSigmaOverdt};

	char Title[][255]={
		"E = BEAM GeV; #theta_{#gamma} (deg); E_{#gamma}",
		"E = BEAM GeV; #theta_{#gamma} (deg); #theta_{p}",
		"E = BEAM GeV; #theta_{#gamma} (deg); #theta_{#gamma}^{cm}",
		"E = BEAM GeV; #theta_{p} (deg); P_{p} (GeV)",
		"E = BEAM GeV; #theta_{#gamma} (deg); d#sigma/dt (pb/GeV^{2})"
	};
	//char Name[][255] ={"EgTheatg_EBEAM.png","TpThetag_EBEAM.png",
	//"TgcmThetag_EBEAM.png","PpThetap_EBEAM.png","XSThetag_EBEAM.png"}; 
	char Name[][255] ={"EgTheatg_EBEAM.png","TpThetag_EBEAM.png",
		"TgcmThetag_EBEAM.png","PpThetap_EBEAM.png","XSThetag_EBEAM.png"}; 

	char strBeam[100];
	sprintf(strBeam,"%.1f",Beam);

	TString strTitle,strName;

	TGraph *Gr=0; 
	TCanvas *cg= new TCanvas("cg","",800,600);

	for(int ig=0;ig<5;ig++)
	{
		cg->Clear();
		cg->cd();
		Gr = new TGraph(N,xx[ig],yy[ig]);
		strTitle = Title[ig];
		strTitle.ReplaceAll("BEAM",strBeam); 
		Gr->SetTitle(strTitle.Data());
		Gr->Draw("AL");
		if(ig==4) cg->SetLogy(1);
		strName = Name[ig];
		strName.ReplaceAll("BEAM",strBeam); 
		strName.Insert(0,"Graph/");
		cg->SaveAs(strName.Data());
	}
}
