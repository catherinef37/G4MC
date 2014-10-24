/*
//macro to skim HRSMC ntuple by various levels
//input: histfile.root output: histfile_skimmed.root 
//usage: root -b -q Skim.C\(1,"track0"\)
//__________________________________________________________________________
*/
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <TROOT.h>
#include <TDirectory.h>
#include <TMath.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
using namespace std;

//##########################config block start#############################
//#define SKIM_DEBUG 0

bool IsLinux=true;


//jlab config
char InPath[]=".";
char OutPath[]=".";

//##########################config block end################################
//-------------------------global variable block start----------------------

bool DoCut(int level=1);
void Skim(int level=1,char *treename="track0");
void SetBranchAddress(TTree *T);


//Declaration of leaves types
Int_t           Index;
Int_t           PdgId;
Int_t           TrackId;
Int_t           TrackClass;
Double_t        X0;
Double_t        Y0;
Double_t        Z0;
Double_t        P0;
Double_t        Theta0;
Double_t        Phi0;
Double_t        X0_tr;
Double_t        Y0_tr;
Double_t        Z0_tr;
Double_t        Theta0_tr;
Double_t        Phi0_tr;
Double_t        Xtg_tr;
Double_t        Ytg_tr;
Double_t        Thetatg_tr;
Double_t        Phitg_tr;
Double_t        Xvb;
Double_t        Yvb;
Double_t        Zvb;
Double_t        Pvb;
Double_t        Thetavb;
Double_t        Phivb;
Double_t        Xvb_tr;
Double_t        Yvb_tr;
Double_t        Zvb_tr;
Double_t        Thetavb_tr;
Double_t        Phivb_tr;
Double_t        Xfp_tr;
Double_t        Yfp_tr;
Double_t        Thetafp_tr;
Double_t        Phifp_tr;
Double_t        Xtg_rec_tr;
Double_t        Ytg_rec_tr;
Double_t        Thetatg_rec_tr;
Double_t        Phitg_rec_tr;
Double_t        X_rec_tr;
Double_t        Y_rec_tr;
Double_t        Z_rec_tr;
Double_t        Theta_rec_tr;
Double_t        Phi_rec_tr;
Double_t        X_rec;
Double_t        Y_rec;
Double_t        Z_rec;
Double_t        P_rec;
Double_t        Theta_rec;
Double_t        Phi_rec;
Double_t        Delta;
Double_t        Delta_rec;
Double_t        X_proj2tg_tr;
Double_t        Y_proj2tg_tr;
Double_t        X_rec2tg_tr;
Double_t        Y_rec2tg_tr;
Double_t        Theta_rec2tg_tr;
Double_t        Phi_rec2tg_tr;
Double_t        P_rec2tg;
Double_t        X_proj2sl_tr;
Double_t        Y_proj2sl_tr;
Double_t        TrackRadlen;
Double_t        Theta0Eff;
Double_t        ElasXS;
Double_t        XS;
Int_t           StepNum;
Double_t        StepX[256];
Double_t        StepY[256];
Double_t        StepZ[256];
Double_t        StepdE[256];
Double_t        StepL[256];
Double_t        StepEkin[256];
Double_t        StepTL[256];
Double_t        StepRadlen[256];
Double_t        StepDsty[256];
Double_t        StepBx[256];
Double_t        StepBy[256];
Double_t        StepBz[256];
Double_t        TrackBdLx;
Double_t        TrackBdLy;
Double_t        TrackBdLz;
Double_t        R0;
Double_t        A0;
Double_t        B0;

void SetBranchAddress(TTree *track0)
{
  if(!track0) return;
   
  // Set branch addresses.
  track0->SetBranchAddress("Index",&Index);
  track0->SetBranchAddress("PdgId",&PdgId);
  track0->SetBranchAddress("TrackId",&TrackId);
  track0->SetBranchAddress("TrackClass",&TrackClass);
  track0->SetBranchAddress("X0",&X0);
  track0->SetBranchAddress("Y0",&Y0);
  track0->SetBranchAddress("Z0",&Z0);
  track0->SetBranchAddress("P0",&P0);
  track0->SetBranchAddress("Theta0",&Theta0);
  track0->SetBranchAddress("Phi0",&Phi0);
  track0->SetBranchAddress("X0_tr",&X0_tr);
  track0->SetBranchAddress("Y0_tr",&Y0_tr);
  track0->SetBranchAddress("Z0_tr",&Z0_tr);
  track0->SetBranchAddress("Theta0_tr",&Theta0_tr);
  track0->SetBranchAddress("Phi0_tr",&Phi0_tr);

  track0->SetBranchAddress("Xvb",&Xvb);
  track0->SetBranchAddress("Yvb",&Yvb);
  track0->SetBranchAddress("Zvb",&Zvb);
  track0->SetBranchAddress("Pvb",&Pvb);
  track0->SetBranchAddress("Thetavb",&Thetavb);
  track0->SetBranchAddress("Phivb",&Phivb);
  track0->SetBranchAddress("Xvb_tr",&Xvb_tr);
  track0->SetBranchAddress("Yvb_tr",&Yvb_tr);
  track0->SetBranchAddress("Zvb_tr",&Zvb_tr);
  track0->SetBranchAddress("Thetavb_tr",&Thetavb_tr);
  track0->SetBranchAddress("Phivb_tr",&Phivb_tr);
  track0->SetBranchAddress("Xfp_tr",&Xfp_tr);
  track0->SetBranchAddress("Yfp_tr",&Yfp_tr);
  track0->SetBranchAddress("Thetafp_tr",&Thetafp_tr);
  track0->SetBranchAddress("Phifp_tr",&Phifp_tr);

  track0->SetBranchAddress("TrackRadlen",&TrackRadlen);
  track0->SetBranchAddress("Theta0Eff",&Theta0Eff);
  track0->SetBranchAddress("ElasXS",&ElasXS);
  track0->SetBranchAddress("XS",&XS);

  //the following leavies might not available in old version
  if(track0->GetBranch("Xtg_tr"))
    {
      track0->SetBranchAddress("Xtg_tr",&Xtg_tr);
      track0->SetBranchAddress("Ytg_tr",&Ytg_tr);
      track0->SetBranchAddress("Thetatg_tr",&Thetatg_tr);
      track0->SetBranchAddress("Phitg_tr",&Phitg_tr);
      track0->SetBranchAddress("Xtg_rec_tr",&Xtg_rec_tr);
      track0->SetBranchAddress("Ytg_rec_tr",&Ytg_rec_tr);
      track0->SetBranchAddress("Thetatg_rec_tr",&Thetatg_rec_tr);
      track0->SetBranchAddress("Phitg_rec_tr",&Phitg_rec_tr);
      track0->SetBranchAddress("X_rec_tr",&X_rec_tr);
      track0->SetBranchAddress("Y_rec_tr",&Y_rec_tr);
      track0->SetBranchAddress("Z_rec_tr",&Z_rec_tr);
      track0->SetBranchAddress("Theta_rec_tr",&Theta_rec_tr);
      track0->SetBranchAddress("Phi_rec_tr",&Phi_rec_tr);
      track0->SetBranchAddress("X_rec",&X_rec);
      track0->SetBranchAddress("Y_rec",&Y_rec);
      track0->SetBranchAddress("Z_rec",&Z_rec);
      track0->SetBranchAddress("P_rec",&P_rec);
      track0->SetBranchAddress("Theta_rec",&Theta_rec);
      track0->SetBranchAddress("Phi_rec",&Phi_rec);
      track0->SetBranchAddress("Delta",&Delta);
      track0->SetBranchAddress("Delta_rec",&Delta_rec);
      track0->SetBranchAddress("X_proj2tg_tr",&X_proj2tg_tr);
      track0->SetBranchAddress("Y_proj2tg_tr",&Y_proj2tg_tr);
      track0->SetBranchAddress("X_rec2tg_tr",&X_rec2tg_tr);
      track0->SetBranchAddress("Y_rec2tg_tr",&Y_rec2tg_tr);
      track0->SetBranchAddress("Theta_rec2tg_tr",&Theta_rec2tg_tr);
      track0->SetBranchAddress("Phi_rec2tg_tr",&Phi_rec2tg_tr);
      track0->SetBranchAddress("P_rec2tg",&P_rec2tg);
      track0->SetBranchAddress("X_proj2sl_tr",&X_proj2sl_tr);
      track0->SetBranchAddress("Y_proj2sl_tr",&Y_proj2sl_tr);
    }


  //the following will be removed by skim
  if(track0->GetBranch("StepNum"))
    {
      track0->SetBranchAddress("StepNum",&StepNum);
      track0->SetBranchAddress("StepX",StepX);
      track0->SetBranchAddress("StepY",StepY);
      track0->SetBranchAddress("StepZ",StepZ);
      track0->SetBranchAddress("StepdE",StepdE);
      track0->SetBranchAddress("StepL",StepL);
      track0->SetBranchAddress("StepEkin",StepEkin);
      track0->SetBranchAddress("StepTL",StepTL);
      track0->SetBranchAddress("StepRadlen",StepRadlen);
      track0->SetBranchAddress("StepDsty",StepDsty);
      track0->SetBranchAddress("StepBx",StepBx);
      track0->SetBranchAddress("StepBy",StepBy);
      track0->SetBranchAddress("StepBz",StepBz);
      track0->SetBranchAddress("TrackBdLx",&TrackBdLx);
      track0->SetBranchAddress("TrackBdLy",&TrackBdLy);
      track0->SetBranchAddress("TrackBdLz",&TrackBdLz);
      track0->SetBranchAddress("R0",&R0);
      track0->SetBranchAddress("A0",&A0);
      track0->SetBranchAddress("B0",&B0);
    }
}


//////////////////////////////////////////////////////////////////

//counters
unsigned int iCounter =0;
unsigned int iCounter0=0;
unsigned int iCounter1=0;
unsigned int iCounter2=0;
unsigned int iCounter3=0;
unsigned int iCounter4=0;
unsigned int iCounter5=0;
unsigned int iCounter6=0;
unsigned int iCounter7=0;
unsigned int iCounter8=0;



//-------------------------global variable block end------------------------

void Skim(int level, char *treename)
{

  if(level<0 ) level=0;

  char cmd[255];
  if(!IsLinux) sprintf(cmd,"mkdir %s",OutPath);
  else sprintf(cmd,"mkdir -p %s",OutPath);
  system(cmd);

  //////////////////////////////////////////////////
  char OldFileName[200],NewFileName[200];

  sprintf(OldFileName,"histfile.root");
  //check this file exist or not
  Long_t id,size,flags,mt;
  Int_t iFound = gSystem->GetPathInfo(OldFileName,&id,&size,&flags,&mt);
  if(iFound!=0) 
    {
      cout<<"histfile.root not found. exit ...\n";
      continue; //File not exist
    }
  sprintf(cmd,"ls -lh %s",OldFileName);
  system(cmd);
  cout <<"\n*****************************New File Start*******************************\n";
  cout << "Skim(): try to skim file "<<OldFileName<<" ..."<<endl;

      
  TFile *oldfile = new TFile(OldFileName);
  TTree *oldtree = (TTree*)oldfile->Get(treename); //get the tree address
 
  SetBranchAddress(oldtree);
  oldtree->SetBranchStatus("*",1);
  if(level>=2)
    {
      oldtree->SetBranchStatus("Step*",0);
      oldtree->SetBranchStatus("R0",0);
      oldtree->SetBranchStatus("A0",0);
      oldtree->SetBranchStatus("B0",0);
      oldtree->SetBranchStatus("TrackBdL*",0);
    }

  //Create a new file + a clone of old tree header. Do not copy events	  
  sprintf(NewFileName,"histfile_skimmed.root");      
  TFile *newfile = new TFile(NewFileName,"recreate");
  TTree *newtree = oldtree->CloneTree(0);

  UInt_t nentries = (UInt_t)oldtree->GetEntries();
  cout<<"Number of entries ="<<nentries<<endl;
#ifdef BONUSSKIM_DEBUG
  if(SKIM_DEBUG>=1)
    {		
      nentries=(nentries<10000) ? nentries:10000; 		
      cout << "Skim(): try to skim "<<nentries<<" for debug..."<<endl;
    }
#endif

  //reset counters
  iCounter=0;
  iCounter0=0;
  iCounter1=0;
  iCounter2=0;
  iCounter3=0;
  iCounter4=0;
  iCounter5=0;
  iCounter6=0;
  iCounter7=0;
  iCounter8=0;

  for (UInt_t i=0;i<nentries; i++) 
    {
      oldtree->GetEntry(i);

      if (DoCut(level)) {newtree->Fill(); iCounter++;}

      if( !((i+1)%1000) || i+1==nentries )
	printf("%6d/%6d events have been processed, %6d pass level %d skim......\r",
	       i+1,nentries,iCounter,level);

    }
  cout<<endl;

  //print the result
#ifdef SKIM_DEBUG
  if(SKIM_DEBUG>=4)	newtree->Print();	
#endif	
  cout<<setw(15)<<"  iCounter="<<setw(6)<<iCounter
      <<setw(15)<<" iCounter0="<<setw(6)<<iCounter0
      <<setw(15)<<" iCounter1="<<setw(6)<<iCounter1
      <<setw(15)<<" iCounter2="<<setw(6)<<iCounter2<<endl
      <<setw(15)<<" iCounter3="<<setw(6)<<iCounter3
      <<setw(15)<<" iCounter4="<<setw(6)<<iCounter4
      <<setw(15)<<" iCounter5="<<setw(6)<<iCounter5
      <<setw(15)<<" iCounter6="<<setw(6)<<iCounter6<<endl
      <<setw(15)<<" iCounter7="<<setw(6)<<iCounter7
      <<setw(15)<<" iCounter8="<<setw(6)<<iCounter8<<endl
      <<endl;

  cout<<"level 0: TrackClass>=0, keep only track0 and config tree"<<endl;
  cout<<"level 1: TrackClass>=1"<<endl;
  cout<<"level 2: TrackClass>=2, remove Step*, TrackBdL*, R, A and B from track0"<<endl;
  cout<<"level 3: TrackClass>=3"<<endl;
  cout<<"level 4: TrackClass>=4"<<endl;
  cout<<"level 5: TrackClass>=5"<<endl;
  cout<<"level 6: TrackClass>=6"<<endl;
  cout<<"level 7: TrackClass>=7"<<endl;
     
  //copy the config tree
  TTree *configtree = (TTree*)oldfile->Get("config"); //get the config tree address
  if(configtree) 
    {
      cout<<"config tree is found, copy config tree ...\n";
      TTree *newconfigtree=configtree->CloneTree(0);

      configtree->SetBranchStatus("*",1);
      configtree->GetEntry(0);

      newconfigtree->Fill();
    }
  newfile->Write();
  delete newfile;
  delete oldfile;

  //sprintf(cmd,"ls -lh %s",NewFileName);
  //system(cmd);
}

	
bool DoCut(int level)
{
  //input:
  /*
  cout<<"level 0: TrackClass>=0"<<endl;
  cout<<"level 1: TrackClass>=1"<<endl;
  cout<<"level 2: TrackClass>=2"<<endl;
  cout<<"level 3: TrackClass>=3"<<endl;
  cout<<"level 4: TrackClass>=4"<<endl;
  cout<<"level 5: TrackClass>=5"<<endl;
  cout<<"level 6: TrackClass>=6"<<endl;
  cout<<"level 7: TrackClass>=7"<<endl;
  */
  
  if (TrackClass<0) return false;

  if(level>=1) 
    {
      if (TrackClass<1) return false;
      iCounter1++;
    }

  if(level>=2)
    {
      if (TrackClass<2) return false;
      iCounter2++;
    }
  if(level>=3)
    {
      if (TrackClass<3) return false;
      iCounter3++;
      if(level>=4)
	{ 
	  if (TrackClass<4) return false;
	  iCounter4++;
	}
      if(level>=5)
	{ 
	  if (TrackClass<5) return false;
	  iCounter5++;
	}
      if(level>=6)
	{ 
	  if (TrackClass<6) return false;
	  iCounter6++;
	}
      if(level>=7)
	{ 
	  if (TrackClass<7) return false;
	  iCounter7++;
	}
    }
  return true;
}

