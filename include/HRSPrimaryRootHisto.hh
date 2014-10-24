// ********************************************************************
// $Id: HRSPrimaryRootHisto.hh,v 1.0, 2010/12/26 HRS Exp $
// --------------------------------------------------------------
// This class is written as a root event generator. Jixie created the root histograms
// "primaryEvent.root" which is based on the BoNuS real data. This class will generater events
// using these histogram. So far only inclusive events e(D,e'p_rtpc) and exclusive pion 
// minus events e(D,e'pi-p_rtpc)p_clas are available. Other channels may be added later.
// 2008/9/18: Jixie removed all G4 stuff in order to test this class in root. Turn off
// the smearing feature for clas proton and NO cuts on RTPC particle 
// Convertion: All energy in GeV and All length in cm unit

#ifndef HRSPrimaryRootHisto_h
#define HRSPrimaryRootHisto_h 1

#include <TROOT.h>
#include <TF1.h>
#include <TRandom.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TObject.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TProfile.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TAxis.h>
#include <TDirectory.h>

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iostream>
using namespace std;

class HRSPrimaryRootHisto 
{
public:
	HRSPrimaryRootHisto();
	virtual ~HRSPrimaryRootHisto();
	void Initialize();

	//e(D,e'p_rtpc)X
	void GetInclusiveEvent(TVector3 &vertex, TVector3 &V3P_rtpc, TVector3 &V3P_el);


	//e(D,e'pi-p_rtpc)p_clas
	void GetExclusiveEvent(TVector3 &vertex, TVector3 &V3P_rtpc, TVector3 &V3P_el,
						   TVector3 &V3P_pi, TVector3 &V3P_clas);


private:		
	double RTPC_P_Threshold(double theta_deg);//return mementum in GeV unit
	void   GetXbyY(TH2 *h2, double y, double &x);
	void   GetYbyX(TH2 *h2, double x, double &y);

public:	
	bool   bHistoOpened;

private:
	TFile *mFile;
	TH1F *hE,*hT;
	TH2F *h2ThetaPhi_el,*h2WprimeQ2_Inc,*h2WprimeQ2_Exc,*h2ThetaZ_rtpc;
	TH3F *h3P_rtpc,*h3P_el,*h3P_pi,*h3P_pr; //Phi(x) Theta(y) P(z) 
	double mBeamE, mTargetMass;  //beam energy abd target mass

};

#endif


