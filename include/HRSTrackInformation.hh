#ifndef HRSTrackInformation_h
#define HRSTrackInformation_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Allocator.hh"
#include "G4VUserTrackInformation.hh"

class HRSTrackInformation : public G4VUserTrackInformation
{

public:
  HRSTrackInformation();
  HRSTrackInformation(const G4Track* aTrack);
  HRSTrackInformation(const HRSTrackInformation* aTrackInfo);
  virtual ~HRSTrackInformation();
  void Print() const;
  void Print(const G4Track* aTrack) const;
   
  inline void *operator new(size_t);
  inline void operator delete(void *aTrackInfo);
  inline int operator ==(const HRSTrackInformation& right) const
  {return (this==&right);}


private:
  G4int                 originalTrackID;
  G4ParticleDefinition* originalDefinition;
  G4ThreeVector         originalPosition;
  G4ThreeVector         originalMomentum;
  G4double              originalEnergy;
  G4double              originalTime;
  G4int                 motherTrackID;
						
  G4int iStepOutputStatus; //-1 fStopAndKill;0 fStopAndAlive;1 write the output
    
public:
  inline G4int GetOriginalTrackID() const {return originalTrackID;}
  inline G4ParticleDefinition* GetOriginalDefinition() const {return originalDefinition;}
  inline G4ThreeVector GetOriginalPosition() const {return originalPosition;}
  inline G4ThreeVector GetOriginalMomentum() const {return originalMomentum;}
  inline G4double GetOriginalEnergy() const {return originalEnergy;}
  inline G4double GetOriginalTime() const {return originalTime;}
  inline G4int GetMotherTrackID() const {return motherTrackID;}

  inline void	SetStepOutputStatus(G4int val){iStepOutputStatus=val;}
  inline G4int	GetStepOutputStatus() const {return iStepOutputStatus;}
};

extern G4Allocator<HRSTrackInformation> aTrackInformationAllocator;

inline void* HRSTrackInformation::operator new(size_t)
{ void* aTrackInfo;
 aTrackInfo = (void*)aTrackInformationAllocator.MallocSingle();
 return aTrackInfo;
}

inline void HRSTrackInformation::operator delete(void *aTrackInfo)
{ aTrackInformationAllocator.FreeSingle((HRSTrackInformation*)aTrackInfo);}

#endif
