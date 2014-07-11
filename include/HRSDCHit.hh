// ********************************************************************
//
// $Id: HRSDCHit.hh,v 1.0, 2010/12/26 HRS Exp $
// --------------------------------------------------------------
//

//
#ifndef HRSDCHit_h
#define HRSDCHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class HRSDCHit : public G4VHit
{
public:

	HRSDCHit(G4int i,G4double t);
	virtual ~HRSDCHit();
	HRSDCHit(const HRSDCHit &right);
	const HRSDCHit& operator=(const HRSDCHit &right);
	int operator==(const HRSDCHit &right) const;


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

	inline void SetInPos(const G4ThreeVector &xyz) { inpos = xyz; }
	inline G4ThreeVector GetInPos() const { return inpos; }

	inline void SetOutPos(const G4ThreeVector &xyz) { outpos = xyz; }
	inline G4ThreeVector GetOutPos() const { return outpos; }

	inline void SetInMom(const G4ThreeVector &pxpypz) { inmom = pxpypz; }
	inline G4ThreeVector GetInMom() const { return inmom; }

	inline void SetOutMom(const G4ThreeVector &pxpypz) { outmom = pxpypz; }
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

typedef G4THitsCollection<HRSDCHit> HRSDCHitsCollection;

extern G4Allocator<HRSDCHit> HRSDCHitAllocator;

inline void* HRSDCHit::operator new(size_t)
{
	void* aHit;
	aHit = (void*)HRSDCHitAllocator.MallocSingle();
	return aHit;
}

inline void HRSDCHit::operator delete(void*aHit)
{
	HRSDCHitAllocator.FreeSingle((HRSDCHit*) aHit);
}

#endif


