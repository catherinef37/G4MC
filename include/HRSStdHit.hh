// ********************************************************************
//
// $Id: HRSStdHit.hh,v 1.0, 2010/12/26 HRS Exp $
// --------------------------------------------------------------
//

//
#ifndef HRSStdHit_h
#define HRSStdHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class HRSStdHit : public G4VHit
{
public:

	HRSStdHit(G4int i,G4double t);
	virtual ~HRSStdHit();
	HRSStdHit(const HRSStdHit &right);
	const HRSStdHit& operator=(const HRSStdHit &right);
	int operator==(const HRSStdHit &right) const;


	inline void *operator new(size_t);
	inline void operator delete(void*aHit);

	void Draw();
	void Print();

private:
	G4int         id;		//the detector copy number
	G4int         pdgid;
	G4double      time;
	G4double      edep;
	G4double      edep_NonIon;  //none ionization energy loss
	G4ThreeVector inpos;
	G4ThreeVector outpos;
	G4ThreeVector inmom;
	G4ThreeVector outmom;
	G4int         trackid;
	G4int         parenttrackid;
	const G4VPhysicalVolume* physV;

public:

	inline G4int GetId() const { return id; }

	inline G4double GetTime() const { return time; }
	inline void SetTime(G4double val) { time = val; }

	inline G4int GetPdgid() const { return pdgid; }
	inline void SetPdgid(G4int val) { pdgid = val; }

	inline void SetInPos(G4ThreeVector &xyz) { inpos = xyz; }
	inline G4ThreeVector GetInPos() const { return inpos; }

	inline void SetOutPos(G4ThreeVector &xyz) { outpos = xyz; }
	inline G4ThreeVector GetOutPos() const { return outpos; }

	inline void SetInMom(G4ThreeVector &pxpypz) { inmom = pxpypz; }
	inline G4ThreeVector GetInMom() const { return inmom; }

	inline void SetOutMom(G4ThreeVector &pxpypz) { outmom = pxpypz; }
	inline G4ThreeVector GetOutMom() const { return outmom; } 

	inline void SetEdep(G4double val) { edep = val; }
	inline void AddEdep(G4double val) { edep += val; }
	inline G4double GetEdep() { return edep; }
	
	inline void SetNonIonEdep(G4double val) { edep_NonIon = val; }
	inline void AddNonIonEdep(G4double val) { edep_NonIon += val; }
	inline G4double GetNonIonEdep() { return edep_NonIon; }

	inline void SetTrackId(G4int val) { trackid = val; }
	inline G4int GetTrackId() const { return trackid; }
	
	inline void SetParentTrackId(G4int val) { parenttrackid = val; }
	inline G4int GetParentTrackId() const { return parenttrackid; }
	
	inline void SetPhysV(G4VPhysicalVolume* val) { physV = val; }
	inline const G4VPhysicalVolume* GetPhysV() const { return physV; }

};

typedef G4THitsCollection<HRSStdHit> HRSStdHitsCollection;

extern G4Allocator<HRSStdHit> HRSStdHitAllocator;

inline void* HRSStdHit::operator new(size_t)
{
	void* aHit;
	aHit = (void*)HRSStdHitAllocator.MallocSingle();
	return aHit;
}

inline void HRSStdHit::operator delete(void*aHit)
{
	HRSStdHitAllocator.FreeSingle((HRSStdHit*) aHit);
}

#endif


