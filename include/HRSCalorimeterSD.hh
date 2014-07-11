// ********************************************************************
//
// $Id: HRSCalorimeterSD.hh,v 1.0, 2010/12/26 HRS Exp $
// --------------------------------------------------------------
//
#ifndef HRSCalorimeterSD_h
#define HRSCalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"
#include "HRSCalorimeterHit.hh"   //for HRSCalorimeterHitsCollection
#include "G4String.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class HRSCalorimeterSD : public G4VSensitiveDetector
{
public:
	HRSCalorimeterSD(G4String name);
	virtual ~HRSCalorimeterSD();

	virtual void Initialize(G4HCofThisEvent *HCE);
	virtual G4bool ProcessHits(G4Step *aStep,G4TouchableHistory *ROhist);
	virtual void EndOfEvent(G4HCofThisEvent*HCE);

private:
	HRSCalorimeterHitsCollection *hitsCollection;
	G4int HCID;
};

#endif

