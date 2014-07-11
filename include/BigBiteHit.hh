/***********************************************
    
    With the help of the work of
    J.Annand and D.Hammilton.
    
    Here we are dealing with hits and sensitive
    detectors ; that's why we create the 
    'sensitive detector class' here.

***********************************************/

#ifndef ndHSD_VAR
#define ndHSD_VAR 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4Track.hh"
#include "G4SDManager.hh"


class G4Step;
class G4HCofThisEvent;

//# Part 1 : Hits class
class BigBiteHit : public G4VHit
{
public:

BigBiteHit();
~BigBiteHit();

  BigBiteHit(const BigBiteHit &BigBiteHit);

  G4int barnumber;  // bar number in scintillator planes
  G4double hitime;
  G4double kinenergy;
  G4double kineloss;
  G4double ltime;  // time of left PM
  G4double rtime;  // time of right PM
  G4double lmestime;  // 'measured' time of left PM
  G4double rmestime;  // 'measured' time of right PM
  G4ThreeVector posofhit;  // position of the hit

  static int nbhtotal;
  static G4double tof;  //time of flight
  static G4double realp;  // real impulsion of the particule in the event
			  // (before energy losses in scintillators)
  static G4double measp;  // 'measured' impulsion with time-of-flight
//  void draw(G4double,G4double,G4double);
//  void Print();
};

typedef G4THitsCollection<BigBiteHit> ndhitscollection;
// Don't know how it works, but if needed ...

extern G4Allocator<BigBiteHit> ndhitsallocator;
// Seems to be needed, but don't know why.


//# Part 2 : detector sensitivity
class BigBiteSD : public G4VSensitiveDetector
{
 //As G4VSD is an abstract class, we have to create this
public:
  BigBiteSD(G4String sdname);
  ~BigBiteSD();

  G4bool ProcessHits(G4Step* step,G4TouchableHistory* touchist);
  void Initialize(G4HCofThisEvent* hcote);
  void EndOfEvent(G4HCofThisEvent* hcote);
  
  G4TouchableHistory* touchable;  // to access informations on hit

  G4int repnum;  // number of the bar which has been hit
  G4String param; // Name of the plane which has been hit
  
  G4int collid;  // pointer on collection ID, 1 per sensitive detector
  ndhitscollection* pointhitcoll;
};

#endif
