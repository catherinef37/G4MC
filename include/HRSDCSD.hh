// ********************************************************************
//
// $Id: HRSDCSD.hh,v 1.0, 2010/12/26 HRS Exp $
// --------------------------------------------------------------
//
#ifndef HRSDCSD_h
#define HRSDCSD_h 1

#include "G4VSensitiveDetector.hh"
#include "HRSDCHit.hh"
#include "G4String.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class HRSDCSD : public G4VSensitiveDetector
{
public:
	HRSDCSD(G4String name);
	virtual ~HRSDCSD();

	virtual void Initialize(G4HCofThisEvent *HCE);
	virtual G4bool ProcessHits(G4Step *aStep,G4TouchableHistory *ROhist);
	virtual void EndOfEvent(G4HCofThisEvent*HCE);

private:
	HRSDCHitsCollection *hitsCollection;
	G4int HCID;
};

#endif

