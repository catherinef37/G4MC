// ********************************************************************
// $Id: HRSPrimaryRootEvent.hh,v 1.0, 2010/12/26 HRS Exp $
// --------------------------------------------------------------
// This class is written as a root event generator. Jixie created the root histograms
// "primaryEvent.root" which is based on the BoNuS real data. This class will generater events
// using these histogram. So far only inclusive events e(D,e'p_rtpc) and exclusive pion 
// minus events e(D,e'pi-p_rtpc)p_clas are available. Other channels may be added later.
// 2008/9/18: Jixie removed all G4 stuff in order to test this class in root. Turn off
// the smearing feature for clas proton and NO cuts on RTPC particle 
// Convertion: All energy in GeV and All length in cm unit

#ifndef HRSPrimaryRootEvent_h
#define HRSPrimaryRootEvent_h 1

#include <TROOT.h>
#include <TF1.h>
#include <TRandom.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TChain.h>
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

#include "HRSNtReader.hh"
#define MaxTrackNum 8

class HRSPrimaryRootEvent 
{
public:
	HRSPrimaryRootEvent();
	virtual ~HRSPrimaryRootEvent();
	
	void LoadNtuple(int trackid, const char* filename, const char* treename, 
		int trignum, int skipnum);

	//read vertex and initial momentum from the tree
	//return   -1: reach end of file
	//       <=-2: error
	//          0: get one event correctly
	int  GetParticle(int trackid,int &pdgid,TVector3& v3x,TVector3& v3p);


private:
	bool EventCut(int trackid);

public:
	HRSNtReader* mNtuple[MaxTrackNum];
	int  mCurrentEntry[MaxTrackNum];

};

#endif


