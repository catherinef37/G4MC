//this script is used to plot the acceptance for one perticular 
//sensertive detector, one has to provide the sdid, which can be
//found in the config tree. simply type config->Show(0) and then all 
//sdid will be shown
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
#include "TGraphErrors.h"

#include "track0.cc"
#include "track1.cc"
#include "track2.cc"
#include "config.cc"


track0 *Gamma=0;
track1 *Proton=0;
track2 *Electron=0;
config *Config=0;

#include "RCSXS.cc"

void DoPlot(const char* inkey="", TTree *gp=0)
{
	system("mkdir -p Graph");          
	if(!gp) gp = (TTree*)gDirectory->Get("gp");

	TH1F *h1A,*h1B; h1A=h1B=0;
	TH2F *h2A,*h2B, *h2C; h2A=h2B=h2C=0;
	TGraph *grA,*grB,*grC; grA=grB=grC=0;

	TCanvas *c42 = new TCanvas("c42","",500,800);
	gp->Draw("Beam>>hb");
	h1A =  (TH1F*) gROOT->FindObject("hb");
	double Beam = h1A->GetMean();

	gp->Draw("SetupHMS>>hh");
	h1A =  (TH1F*) gROOT->FindObject("hh");
	int SetupHMS = h1A->GetMean();

	gp->Draw("GammaAngle>>hga");
	h1B =  (TH1F*) gROOT->FindObject("hga");
	double GammaAngle = h1B->GetMean(); 

	gp->Draw("ProtonAngle>>hpa");
	h1B =  (TH1F*) gROOT->FindObject("hpa");
	double ProtonAngle = h1B->GetMean();
	c42->Clear();

	char key[100];                       
	char theCut[255];                    
	if(SetupHMS)                         
	{                                  
		sprintf(key,"NPS_%.0fdeg_HMS_%.0fdeg%s",GammaAngle*57.3,ProtonAngle*57.3,inkey);
	}                                                                            
	else                                                                           
	{                                                                            
		sprintf(key,"LAC_%.0fdeg_SBS_%.0fdeg%s",GammaAngle*57.3,ProtonAngle*57.3,inkey);
	}                             

	char detector[100];
	if(SetupHMS)                                                                            
	{                                                                                     
		sprintf(detector,"HMS(%.0f deg)",ProtonAngle*57.3);                                
		sprintf(theCut,"Pvb_g/P0_g>0.9 && Pvb_p/P0_p>0.9");
	}   
	else    
	{
		sprintf(detector,"SBS(%.0f deg)",ProtonAngle*57.3); 
		sprintf(theCut,"Pvb_g/P0_g>0.9 && Pvb_p/P0_p>0.9"); 
	}

	c42->Divide(1,2);
	c42->cd(1);      
	gp->Draw("Theta0_cm_g*57.3:Theta0_g*57.3>>h2TcmT_g",theCut,"");
	grA = (TGraph*) gROOT->FindObject("Graph")->Clone("grA");            
	grA->SetTitle(Form("E=%.1f GeV; #theta_{#gamma}^{lab} (deg);#theta_{#gamma}^{cm} (deg)",Beam));
	grA->Draw("AP");                                                                                   

	c42->cd(2);
	gp->Draw("XS*1000.:Theta0_g*57.3>>h2XST_g",theCut,"");
	grB = (TGraph*) gROOT->FindObject("Graph")->Clone("grB");            
	grB->SetTitle(Form("E=%.1f GeV; #theta_{#gamma}^{lab} (deg);d#sigma/dt (pb/GeV^{2})",Beam));
	grB->Draw("AP");                                                       


	c42->Modified();
	c42->SaveAs(Form("Graph/XS_E%.1f_%s.png",Beam,key));
}

void AnaWACS()
{
	gStyle->SetOptStat(0);
	//set the draw-option "text" format
	gStyle->SetPaintTextFormat(".0f");
	system("mkdir -p Graph");

	const double deg = atan(1.0)/45.;

	Gamma  = new track0();
	Proton = new track1();
	Electron = new track2();
	Config = new config();

	const double kMassPr=0.9383;

	double Beam, Ei;
	int    SetupHMS=0;
	double GammaAngle,ProtonAngle;

	double X0,Y0,Z0;

	//GAMMA
	double  Pol;

	double P0_g,Theta0_g,Phi0_g;
	double Xvb_g,Yvb_g,Zvb_g;
	double Pvb_g,Thetavb_g,Phivb_g;
	double P0_cm_g,Theta0_cm_g,Phi0_cm_g;

	double X0_tr_g,Y0_tr_g,Z0_tr_g;
	double Theta0_tr_g,Phi0_tr_g;
	double Xtg_tr_g,Ytg_tr_g;
	double Thetatg_tr_g,Phitg_tr_g;
	double Xvb_tr_g,Yvb_tr_g,Zvb_tr_g;
	double Thetavb_tr_g,Phivb_tr_g;

	//electron
	double P0_e,Theta0_e,Phi0_e;
	double Xvb_e,Yvb_e,Zvb_e;
	double Pvb_e,Thetavb_e,Phivb_e;
	double P0_cm_e,Theta0_cm_e,Phi0_cm_e;

	double X0_tr_e,Y0_tr_e,Z0_tr_e;
	double Theta0_tr_e,Phi0_tr_e;
	double Xtg_tr_e,Ytg_tr_e;
	double Thetatg_tr_e,Phitg_tr_e;
	double Xvb_tr_e,Yvb_tr_e,Zvb_tr_e;
	double Thetavb_tr_e,Phivb_tr_e;

	//PROTON
	double P0_p,Theta0_p,Phi0_p;
	double Xvb_p,Yvb_p,Zvb_p;
	double Pvb_p,Thetavb_p,Phivb_p;
	double P0_cm_p,Theta0_cm_p,Phi0_cm_p;

	double X0_tr_p,Y0_tr_p,Z0_tr_p;
	double Theta0_tr_p,Phi0_tr_p;
	double Xtg_tr_p,Ytg_tr_p;
	double Thetatg_tr_p,Phitg_tr_p;
	double Xvb_tr_p,Yvb_tr_p,Zvb_tr_p;
	double Thetavb_tr_p,Phivb_tr_p;

	double s,t,u;
	double XS;

	///////////////////////////////////////////////////////////////////////
	//read config
	Config->GetEntry(0);
	Beam=Config->Beam;
	SetupHMS=Config->SetupHMS;
	GammaAngle=Config->VDAngle;
	ProtonAngle=(SetupHMS?Config->HMSAngle:Config->SuperBigBiteAngle);

	char pOutFileName[100];
	sprintf(pOutFileName, "nt_gp_g%.0f_p%.0f_E%.1f.root",
		GammaAngle/deg,ProtonAngle/deg,Beam);
	TFile *pFile=new TFile(pOutFileName,"RECREATE");
	TTree *pTree=new TTree("gp","RCS gamma-p coincident events");     

	pTree->Branch("Beam",&Beam,"Beam/D"); 
	pTree->Branch("Ei",&Ei,"Ei/D"); 
	pTree->Branch("SetupHMS",&SetupHMS,"SetupHMS/I");
	pTree->Branch("GammaAngle",&GammaAngle,"GammaAngle/D");
	pTree->Branch("ProtonAngle",&ProtonAngle,"ProtonAngle/D");  

	pTree->Branch("X0",&X0,"X0/D");        
	pTree->Branch("Y0",&Y0,"Y0/D");        
	pTree->Branch("Z0",&Z0,"Z0/D");  

	//gamma
	pTree->Branch("Pol_g",&Pol_g,"Pol_g/D"); 
	pTree->Branch("P0_g",&P0_g,"P0_g/D");  
	pTree->Branch("Theta0_g",&Theta0_g,"Theta0_g/D");
	pTree->Branch("Phi0_g",&Phi0_g,"Phi0_g/D");
	pTree->Branch("Xvb_g",&Xvb_g,"Xvb_g/D");        
	pTree->Branch("Yvb_g",&Yvb_g,"Yvb_g/D");        
	pTree->Branch("Zvb_g",&Zvb_g,"Zvb_g/D");  
	pTree->Branch("Pvb_g",&Pvb_g,"Pvb_g/D");  
	pTree->Branch("Thetavb_g",&Thetavb_g,"Thetavb_g/D");
	pTree->Branch("Phivb_g",&Phivb_g,"Phivb_g/D");
	pTree->Branch("P0_cm_g",&P0_cm_g,"P0_cm_g/D");  
	pTree->Branch("Theta0_cm_g",&Theta0_cm_g,"Theta0_cm_g/D");
	pTree->Branch("Phi0_cm_g",&Phi0_cm_g,"Phi0_cm_g/D");

	pTree->Branch("X0_tr_g",&X0_tr_g,"X0_tr_g/D");        
	pTree->Branch("Y0_tr_g",&Y0_tr_g,"Y0_tr_g/D");        
	pTree->Branch("Z0_tr_g",&Z0_tr_g,"Z0_tr_g/D");  
	pTree->Branch("Theta0_tr_g",&Theta0_tr_g,"Theta0_tr_g/D");
	pTree->Branch("Phi0_tr_g",&Phi0_tr_g,"Phi0_tr_g/D");

	pTree->Branch("Xtg_tr_g",&Xtg_tr_g,"Xtg_tr_g/D");        
	pTree->Branch("Ytg_tr_g",&Ytg_tr_g,"Ytg_tr_g/D");        
	pTree->Branch("Thetatg_tr_g",&Thetatg_tr_g,"Thetatg_tr_g/D");
	pTree->Branch("Phitg_tr_g",&Phitg_tr_g,"Phitg_tr_g/D");

	pTree->Branch("Xvb_tr_g",&Xvb_tr_g,"Xvb_tr_g/D");        
	pTree->Branch("Yvb_tr_g",&Yvb_tr_g,"Yvb_tr_g/D");        
	pTree->Branch("Zvb_tr_g",&Zvb_tr_g,"Zvb_tr_g/D");  
	pTree->Branch("Thetavb_tr_g",&Thetavb_tr_g,"Thetavb_tr_g/D");
	pTree->Branch("Phivb_tr_g",&Phivb_tr_g,"Phivb_tr_g/D");

	//proton
	pTree->Branch("P0_p",&P0_p,"P0_p/D");  
	pTree->Branch("Theta0_p",&Theta0_p,"Theta0_p/D");
	pTree->Branch("Phi0_p",&Phi0_p,"Phi0_p/D");
	pTree->Branch("Xvb_p",&Xvb_p,"Xvb_p/D");        
	pTree->Branch("Yvb_p",&Yvb_p,"Yvb_p/D");        
	pTree->Branch("Zvb_p",&Zvb_p,"Zvb_p/D");  
	pTree->Branch("Pvb_p",&Pvb_p,"Pvb_p/D");  
	pTree->Branch("Thetavb_p",&Thetavb_p,"Thetavb_p/D");
	pTree->Branch("Phivb_p",&Phivb_p,"Phivb_p/D");
	pTree->Branch("P0_cm_p",&P0_cm_p,"P0_cm_p/D");  
	pTree->Branch("Theta0_cm_p",&Theta0_cm_p,"Theta0_cm_p/D");
	pTree->Branch("Phi0_cm_p",&Phi0_cm_p,"Phi0_cm_p/D");

	pTree->Branch("X0_tr_p",&X0_tr_p,"X0_tr_p/D");        
	pTree->Branch("Y0_tr_p",&Y0_tr_p,"Y0_tr_p/D");        
	pTree->Branch("Z0_tr_p",&Z0_tr_p,"Z0_tr_p/D");  
	pTree->Branch("Theta0_tr_p",&Theta0_tr_p,"Theta0_tr_p/D");
	pTree->Branch("Phi0_tr_p",&Phi0_tr_p,"Phi0_tr_p/D");

	pTree->Branch("Xtg_tr_p",&Xtg_tr_p,"Xtg_tr_p/D");        
	pTree->Branch("Ytg_tr_p",&Ytg_tr_p,"Ytg_tr_p/D");        
	pTree->Branch("Thetatg_tr_p",&Thetatg_tr_p,"Thetatg_tr_p/D");
	pTree->Branch("Phitg_tr_p",&Phitg_tr_p,"Phitg_tr_p/D");

	pTree->Branch("Xvb_tr_p",&Xvb_tr_p,"Xvb_tr_p/D");        
	pTree->Branch("Yvb_tr_p",&Yvb_tr_p,"Yvb_tr_p/D");        
	pTree->Branch("Zvb_tr_p",&Zvb_tr_p,"Zvb_tr_p/D");  
	pTree->Branch("Thetavb_tr_p",&Thetavb_tr_p,"Thetavb_tr_p/D");
	pTree->Branch("Phivb_tr_p",&Phivb_tr_p,"Phivb_tr_p/D");

	//electron which match to the proton,
	//this is used to study the electron background
	pTree->Branch("P0_e",&P0_e,"P0_e/D");  
	pTree->Branch("Theta0_e",&Theta0_e,"Theta0_e/D");
	pTree->Branch("Phi0_e",&Phi0_e,"Phi0_e/D");
	pTree->Branch("Xvb_e",&Xvb_e,"Xvb_e/D");        
	pTree->Branch("Yvb_e",&Yvb_e,"Yvb_e/D");        
	pTree->Branch("Zvb_e",&Zvb_e,"Zvb_e/D");  
	pTree->Branch("Pvb_e",&Pvb_e,"Pvb_e/D");  
	pTree->Branch("Thetavb_e",&Thetavb_e,"Thetavb_e/D");
	pTree->Branch("Phivb_e",&Phivb_e,"Phivb_e/D");
	pTree->Branch("P0_cm_e",&P0_cm_e,"P0_cm_e/D");  
	pTree->Branch("Theta0_cm_e",&Theta0_cm_e,"Theta0_cm_e/D");
	pTree->Branch("Phi0_cm_e",&Phi0_cm_e,"Phi0_cm_e/D");

	pTree->Branch("X0_tr_e",&X0_tr_e,"X0_tr_e/D");        
	pTree->Branch("Y0_tr_e",&Y0_tr_e,"Y0_tr_e/D");        
	pTree->Branch("Z0_tr_e",&Z0_tr_e,"Z0_tr_e/D");  
	pTree->Branch("Theta0_tr_e",&Theta0_tr_e,"Theta0_tr_e/D");
	pTree->Branch("Phi0_tr_e",&Phi0_tr_e,"Phi0_tr_e/D");

	pTree->Branch("Xtg_tr_e",&Xtg_tr_e,"Xtg_tr_e/D");        
	pTree->Branch("Ytg_tr_e",&Ytg_tr_e,"Ytg_tr_e/D");        
	pTree->Branch("Thetatg_tr_e",&Thetatg_tr_e,"Thetatg_tr_e/D");
	pTree->Branch("Phitg_tr_e",&Phitg_tr_e,"Phitg_tr_e/D");

	pTree->Branch("Xvb_tr_e",&Xvb_tr_e,"Xvb_tr_e/D");        
	pTree->Branch("Yvb_tr_e",&Yvb_tr_e,"Yvb_tr_e/D");        
	pTree->Branch("Zvb_tr_e",&Zvb_tr_e,"Zvb_tr_e/D");  
	pTree->Branch("Thetavb_tr_e",&Thetavb_tr_e,"Thetavb_tr_e/D");
	pTree->Branch("Phivb_tr_e",&Phivb_tr_e,"Phivb_tr_e/D");

	pTree->Branch("s",&s,"s/D");      
	pTree->Branch("t",&t,"t/D");      
	pTree->Branch("u",&u,"u/D");      
	pTree->Branch("XS",&XS,"XS/D");

	///////////////////////////////////////////////////////////////////////

	//start processing event by event
	Long64_t nentries = Gamma->fChain->GetEntriesFast();
	TLorentzVector V4Gamma, V4Proton, V4Ei, V4Tg;
	TLorentzVector V4Beta, V4Gamma_cm, V4Proton_cm;
	TVector3 V3Gamma, V3Proton, V3Beta;

	V4Tg.SetXYZT(0,0,0,kMassPr);

	Long64_t nb0 = 0, nb1 = 0, nb2 = 0;
	for (Long64_t i=0; i<nentries;i++) 
	{
		if(!((i+1)%1000))
			printf("processing event %6d / %6d \r",int(i+1),int(nentries));

		//do proton first since its acceptance is smaller
		//read proton
		nb1 = Proton->fChain->GetEntry(i);        
		if(nb1<=0) break;
		//apply cuts
		if ( Proton->Cut(i) < 0) continue;
		//apply other cuts
		if( Proton->Pvb/Proton->P0 < 0.8 )
		{
			continue;
		}

		//read gamma
		nb0 = Gamma->fChain->GetEntry(i); 
		if(nb0<=0) break;
		//apply cuts
		if (Gamma->Cut(i) < 0) continue;
		//apply other cuts
		if(Gamma->Pvb/Gamma->P0 < 0.8 )
		{
			continue;
		}

		//reach here, it should now be an  e-p elastic events
		//read electron tree
		nb2 = Electron->fChain->GetEntry(i); 
		//if(nb2<=0) break;
		//apply cuts
		if (Electron->Cut(i) < 0) continue;
	

		//load variables

		X0=Gamma->X0;
		Y0=Gamma->Y0;
		Z0=Gamma->Z0;
		Ei=Gamma->Ei;

		P0_g=Gamma->P0;
		Theta0_g=Gamma->Theta0;
		Phi0_g=Gamma->Phi0;
		Xvb_g=Gamma->Xvb;
		Yvb_g=Gamma->Yvb;
		Zvb_g=Gamma->Zvb;
		Pvb_g=Gamma->Pvb;
		Thetavb_g=Gamma->Thetavb;
		Phivb_g=Gamma->Phivb;

		X0_tr_g=Gamma->X0_tr;
		Y0_tr_g=Gamma->Y0_tr;
		Z0_tr_g=Gamma->Z0_tr;
		Theta0_tr_g=Gamma->Theta0_tr;
		Phi0_tr_g=Gamma->Phi0_tr;

		Xtg_tr_g=Gamma->Xtg_tr;
		Ytg_tr_g=Gamma->Ytg_tr;
		Thetatg_tr_g=Gamma->Thetatg_tr;
		Phitg_tr_g=Gamma->Phitg_tr;

		Xvb_tr_g=Gamma->Xvb_tr;
		Yvb_tr_g=Gamma->Yvb_tr;
		Zvb_tr_g=Gamma->Zvb_tr;
		Thetavb_tr_g=Gamma->Thetavb_tr;
		Phivb_tr_g=Gamma->Phivb_tr;

		double y=Ei/Beam;
		double Pol_e = 0.9;
		Pol_g =	 Pol_e*(4*y-y*y)/(4-4*y+3*y*y)

		//electron
		P0_e=Electron->P0;
		Theta0_e=Electron->Theta0;
		Phi0_e=Electron->Phi0;
		Xvb_e=Electron->Xvb;
		Yvb_e=Electron->Yvb;
		Zvb_e=Electron->Zvb;
		Pvb_e=Electron->Pvb;
		Thetavb_e=Electron->Thetavb;
		Phivb_e=Electron->Phivb;

		X0_tr_e=Electron->X0_tr;
		Y0_tr_e=Electron->Y0_tr;
		Z0_tr_e=Electron->Z0_tr;
		Theta0_tr_e=Electron->Theta0_tr;
		Phi0_tr_e=Electron->Phi0_tr;

		Xtg_tr_e=Electron->Xtg_tr;
		Ytg_tr_e=Electron->Ytg_tr;
		Thetatg_tr_e=Electron->Thetatg_tr;
		Phitg_tr_e=Electron->Phitg_tr;

		Xvb_tr_e=Electron->Xvb_tr;
		Yvb_tr_e=Electron->Yvb_tr;
		Zvb_tr_e=Electron->Zvb_tr;
		Thetavb_tr_e=Electron->Thetavb_tr;
		Phivb_tr_e=Electron->Phivb_tr;

		//proton
		P0_p=Proton->P0;
		Theta0_p=Proton->Theta0;
		Phi0_p=Proton->Phi0;
		Xvb_p=Proton->Xvb;
		Yvb_p=Proton->Yvb;
		Zvb_p=Proton->Zvb;
		Pvb_p=Proton->Pvb;
		Thetavb_p=Proton->Thetavb;
		Phivb_p=Proton->Phivb;

		X0_tr_p=Proton->X0_tr;
		Y0_tr_p=Proton->Y0_tr;
		Z0_tr_p=Proton->Z0_tr;
		Theta0_tr_p=Proton->Theta0_tr;
		Phi0_tr_p=Proton->Phi0_tr;

		Xtg_tr_p=Proton->Xtg_tr;
		Ytg_tr_p=Proton->Ytg_tr;
		Thetatg_tr_p=Proton->Thetatg_tr;
		Phitg_tr_p=Proton->Phitg_tr;

		Xvb_tr_p=Proton->Xvb_tr;
		Yvb_tr_p=Proton->Yvb_tr;
		Zvb_tr_p=Proton->Zvb_tr;
		Thetavb_tr_p=Proton->Thetavb_tr;
		Phivb_tr_p=Proton->Phivb_tr;

		//cout<<"P0_g="<<P0_g<<"  P0_p="<<P0_p<<"  Pvb_g="<<Pvb_g<<"  Pvb_p="<<Pvb_p<<endl;

		//Reach here: both gamma and proton detected
		//fill 4-vectors then boost them to cm frame
		V4Ei.SetXYZT(0,0,Ei,Ei);
		V4Beta = V4Ei + V4Tg;
		V3Beta.SetXYZ(0,0,-V4Beta.P()/V4Beta.T());

		V3Gamma.SetMagThetaPhi(P0_g,Theta0_g,Phi0_g);
		V4Gamma.SetVectM(V3Gamma,0);

		V3Proton.SetMagThetaPhi(P0_p,Theta0_p,Phi0_p);
		V4Proton.SetVectM(V3Proton,kMassPr);

		s = (V4Ei+V4Tg).M2();
		t = (V4Ei-V4Gamma).M2();
		u = (V4Ei-V4Proton).M2();

		//boost
		V4Gamma_cm = V4Gamma;
		V4Gamma_cm.Boost(V3Beta);
		P0_cm_g = V4Gamma_cm.P();
		//P0_cm_g = V4Gamma_cm.Vect().Mag();  //work but slow
		Theta0_cm_g = V4Gamma_cm.Theta();
		Phi0_cm_g = V4Gamma_cm.Phi();

		V4Proton_cm = V4Proton;
		V4Proton_cm.Boost(V3Beta); 
		P0_cm_p = V4Proton_cm.P();
		//P0_cm_p = V4Proton_cm.Vect().Mag(); //work but slow
		Theta0_cm_p = V4Proton_cm.Theta();
		Phi0_cm_p = V4Proton_cm.Phi();

		//cout<<"Theta0_g="<<Theta0_cm_g*57.3<<"  P0_g="<<P0_g
		//    <<"  P0_cm_g="<<Pvb_g<<"  Ei="<<Ei<<"  t="<<t<<endl;

		//calculate the Xs
		XS = GetRCSXS(Ei,Theta0_cm_g);

		pTree->Fill();

		//reset
		P0_g=Theta0_g=Phi0_g=-99.0;
		Xvb_g=Yvb_g=Zvb_g=-99.0;
		Pvb_g=Thetavb_g=Phivb_g=-99.0;
		P0_cm_g=Theta0_cm_g=Phi0_cm_g=-99.0;

		P0_p=Theta0_p=Phi0_p=-99.0;
		Xvb_p=Yvb_p=Zvb_p=-99.0;
		Pvb_p=Thetavb_p=Phivb_p=-99.0;
		P0_cm_p=Theta0_cm_p=Phi0_cm_p=-99.0;

	}

	pFile->Write("", TObject::kOverwrite);

	//DoPlot("",pTree);

	pFile->Close();

	printf("\n\n");

}
