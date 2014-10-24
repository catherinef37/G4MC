// ********************************************************************
//
// $Id: HRSEventAction.hh,v 1.0, 2010/12/26 HRS Exp $
// --------------------------------------------------------------
//
//
#ifndef HRSEventAction_h
#define HRSEventAction_h 1

#include "G4HCofThisEvent.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"

//this number define the maximum number of traj that will be flaged 
//the flag stored  in mStoreTrackIdList[MaxStoredTrjN] determine if the traj 
//will be stored or not. Only traj that fired SD and its parent track will be flaged
//Note that the buffer of the root tree is defined by MaxParticle in HRSRootTree.hh

#define MaxStoredTrjN 20480

class HRSEventActionMessenger;
class G4PrimaryParticle;

class HRSEventAction : public G4UserEventAction
{
public:
	HRSEventAction();
	virtual ~HRSEventAction();

public:
	virtual void BeginOfEventAction(const G4Event*);
	virtual void EndOfEventAction(const G4Event*);
	
	void ProcessDCSDHits(G4HCofThisEvent* HCE);
	void ProcessStdSDHits(G4HCofThisEvent* HCE);
	void ProcessCalorimeterSDHits(G4HCofThisEvent* HCE);

	void ProcessSDHits(G4HCofThisEvent* HCE);


	void PrintPrimary(G4PrimaryParticle* pp,G4int ind);

private:

	HRSEventActionMessenger* messenger;
	G4int verboseLevel;
	//the trajectory flags which will be store into the ntuple
	//only MaxStoredTrjN particles will be considered
	bool mStoreTrackIdList[MaxStoredTrjN];

public:
	inline void SetVerbose(G4int val) { verboseLevel = val; }
	inline G4int GetVerbose() const { return verboseLevel; }
};

#endif
