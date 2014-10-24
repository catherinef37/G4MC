// ********************************************************************
//
// $Id: HRSStdSD.hh,v 1.0, 2010/12/26 HRS Exp $
// --------------------------------------------------------------
//
#ifndef HRSStdSD_h
#define HRSStdSD_h 1

#include "G4VSensitiveDetector.hh"
#include "HRSStdHit.hh"   //for HRSStdHitsCollection
#include "G4String.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class HRSStdSD : public G4VSensitiveDetector
{
public:
	HRSStdSD(G4String name);
	virtual ~HRSStdSD();

	virtual void Initialize(G4HCofThisEvent *HCE);
	virtual G4bool ProcessHits(G4Step *aStep,G4TouchableHistory *ROhist);
	virtual void EndOfEvent(G4HCofThisEvent*HCE);

private:
	HRSStdHitsCollection *hitsCollection;
	G4int HCID;
};

#endif

