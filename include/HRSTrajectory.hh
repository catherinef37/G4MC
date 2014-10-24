// ********************************************************************
//
// $Id: HRSTrajectory.hh,v 1.0, 2010/12/26 HRS Exp $
// --------------------------------------------------------------
//

class HRSTrajectory;

#ifndef HRSTrajectory_h
#define HRSTrajectory_h 1

#include "G4VTrajectory.hh"
#include "G4Allocator.hh"
#include <stdlib.h>                 // Include from 'system'
#include "G4ios.hh"                 // Include from 'system'
#include <vector>					//
#include <iostream>					//
#include "globals.hh"               // Include from 'global'
#include "G4ParticleTable.hh"		// Include from 'particle+matter'
#include "G4ParticleTypes.hh"		// Include from 'particle+matter'
#include "G4ParticleDefinition.hh"  // Include from 'particle+matter'
#include "G4TrajectoryPoint.hh"     // Include from 'tracking'
#include "G4Track.hh"
#include "G4Step.hh"

typedef std::vector<G4VTrajectoryPoint*> HRSTrajectoryPointContainer;

class G4Polyline;                   // Forward declaration.

///////////////////
class HRSTrajectory : public G4VTrajectory
///////////////////
{

//--------
   public:
//--------

// Constructor/Destrcutor

   HRSTrajectory();

   HRSTrajectory(const G4Track* aTrack);
   HRSTrajectory(HRSTrajectory &);
   virtual ~HRSTrajectory();

// Operators
   inline void* operator new(size_t);
   inline void  operator delete(void*);
   inline int operator == (const HRSTrajectory& right) const
   {return (this==&right);}

// Get/Set functions
   inline G4int GetTrackID() const { return trackID; }
   inline G4int GetParentID() const  { return parentTrackID; }
   inline G4String GetParticleName() const  { return particleName; }
   inline G4double GetCharge() const  { return PDGCharge; }
   inline G4int GetPDGEncoding() const  { return PDGEncoding; }
   inline G4double GetInitialKineticEnergy() const  { return initialKineticEnergy; }
   inline G4ThreeVector GetInitialMomentum() const  { return initialMomentum; }
   inline G4double GetInitialTime() const  { return initialTime; }


// Other member functions
   virtual void ShowTrajectory() const;
   virtual void ShowTrajectory(std::ostream& o) const;
   virtual void DrawTrajectory(G4int i_mode=0) const;
   virtual void AppendStep(const G4Step* aStep);
   virtual void MergeTrajectory(G4VTrajectory* secondTrajectory);

   virtual int GetPointEntries() const { return positionRecord->size(); }
   virtual G4VTrajectoryPoint* GetPoint(G4int i) const { return (*positionRecord)[i]; }

   G4ParticleDefinition* GetParticleDefinition() { return particleDefinition; }

   virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
   virtual std::vector<G4AttValue>* CreateAttValues() const;

//---------
   private:
//---------

  HRSTrajectoryPointContainer* positionRecord;

  G4int    trackID;
  G4int    parentTrackID;
  G4String particleName;
  G4double PDGCharge;
  G4int    PDGEncoding;
  G4double initialTime;
  G4double initialKineticEnergy;
  G4ThreeVector initialMomentum;

  G4ParticleDefinition* particleDefinition;

};

#if defined G4TRACKING_ALLOC_EXPORT
  extern G4DLLEXPORT G4Allocator<HRSTrajectory> myTrajectoryAllocator;
#else
  extern G4DLLIMPORT G4Allocator<HRSTrajectory> myTrajectoryAllocator;
#endif

inline void* HRSTrajectory::operator new(size_t)
{
  void* aTrajectory;
  aTrajectory = (void*)myTrajectoryAllocator.MallocSingle();
  return aTrajectory;
}

inline void HRSTrajectory::operator delete(void* aTrajectory)
{
  myTrajectoryAllocator.FreeSingle((HRSTrajectory*)aTrajectory);
}

#endif
