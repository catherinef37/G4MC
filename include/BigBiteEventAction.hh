/***********************************************
    
    BigBiteEventAction.hh
    
    Create the user class BigBiteEventAction for
    the (mainly end-of) event actions management.
    
***********************************************/

#ifndef ndEA_VAR
#define ndEA_VAR 1


#include "G4SDManager.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"

#include "TCint.h"
#include "TTree.h"

class BigBiteEventAction : public G4UserEventAction
{
public:
  BigBiteEventAction(TTree*,char);


  ~BigBiteEventAction();
  
  char a; // 'a' if we have scintillators, 'c' if we have wire chambers
  
  // For Tree
  TTree* tnd;
  
  
  // These are the variable used for tree filling, 
  // for signification of each one, see ndEA.cc

  G4int eventn, istree;
  G4double theta,phi,ppart;
  G4double prec,abser,reler,recthe,recphi,prec0,outthe,hbna,hbnE;
  G4double ptof,ptofr,zekina;
  G4double rthe, rphi;
  G4double  mesrpp,  mestpp,mesrpp0,  thetmp, rtofpp;

//   // No multihit ADC assumed
  G4double ekinapta[56]; 
  G4double ekinaptde[24];
  G4double ekinapte[24];
  G4double ekinaptp0[30];  
  G4double ekinaptp1[30];
  G4double ekinaptp2[24];
  G4double ekinaptp3[6];
  G4double ekinaptp4[4];
  G4double ekinaptp5[2];
  G4double ekinaptp6[4];  
  G4double ekinaptp7[6];
  G4double ekinaptp8[12];

  //Multihit TDC assumed. 
  G4double auxhtimeta[56][6];
  G4double lmestimede[24][6];
  G4double lmestimee[24][6];
  G4double lmestimep0[30][6]; 
  G4double lmestimep1[30][6];
  G4double lmestimep2[24][6];
  G4double lmestimep3[6][6];
  G4double lmestimep4[4][6];

  G4double lmestimep5[2][6];

  G4double lmestimep6[4][6]; 
  G4double lmestimep7[6][6];
  G4double lmestimep8[12][6];
  G4double rmestimede[24][6];
  G4double rmestimee[24][6];
  G4double rmestimep0[30][6]; 
  G4double rmestimep1[30][6];
  G4double rmestimep2[24][6];
  G4double rmestimep3[6][6];
  G4double rmestimep4[4][6];
  G4double rmestimep5[2][6];
  G4double rmestimep6[4][6]; 
  G4double rmestimep7[6][6];
  G4double rmestimep8[12][6];
   
public:
  // what has to be done at the beginning of each event
  void BeginOfEventAction(const G4Event*);

  // what has to be done at the end of each event
  void EndOfEventAction(const G4Event*);    


private:
  G4SDManager* sdman;
  G4int collidwc1;   //
  G4int collidwc2;   //
  G4int collidaux;   // pointers to get collectionsID of the events
  G4int colliddlt;
  G4int collidmain;
  G4int collidp0;
  G4int collidp1;
  G4int collidp2;
  G4int collidp3;
  G4int collidp4; 
  G4int collidp5;
  G4int collidp6;
  G4int collidp7;
  G4int collidp8;
  
  
  G4ThreeVector hrealpwc1;  // real hit position in wire chamber 1
  G4ThreeVector hrealpwc2;  // real hit position in wire chamber 2
  G4ThreeVector hrealpaux;  // real hit position in auxilliary plane
  G4ThreeVector hrealpdlt;  // real hit position in delta-E plane
  G4ThreeVector hrealpmain; // real hit position in main-E plane
  
};

#endif

