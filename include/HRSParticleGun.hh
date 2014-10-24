// ********************************************************************
// G4ParticleGun always print out some warning message while seting the
// energy or momentum of the primary particle, thie really slow down the
// speed and very anoying. I wrote this class to override it
//
// $Id: HRSParticleGun.hh,v 1.0 2011/2/11 HRS Exp $
// GEANT4 tag $Name: geant4-09-04 $
// ********************************************************************

#ifndef HRSParticleGun_h
#define HRSParticleGun_h 1


#include "globals.hh"
#include "G4ParticleGun.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryVertex.hh"
#include "G4ParticleMomentum.hh"

class G4Event;

class HRSParticleGun:public G4ParticleGun
{
  public: // with description
     HRSParticleGun();
     HRSParticleGun(G4int numberofparticles);
     HRSParticleGun(G4ParticleDefinition * particleDef, 
                   G4int numberofparticles = 1);
     // costructors. "numberofparticles" is number of particles to be shoot at one invokation
     // of GeneratePrimaryVertex() method. All paricles are shoot with the same physical
     // quantities.

     virtual ~HRSParticleGun();

  public: // with description
     virtual void GeneratePrimaryVertex(G4Event* evt);
     // Creates a primary vertex at the given point and put primary particles to it.
     // Followings are set methods for the particle properties.
     //   SetParticleDefinition should be called first.  
     //   By using SetParticleMomentum(), both particle_momentum_direction and
     //   particle_energy(Kinetic Energy) are set.
     //   
  
     void SetParticleEnergy(G4double aKineticEnergy);
     void SetParticleMomentum(G4double aMomentum);
     void SetParticleMomentum(G4ParticleMomentum aMomentum);


};

#endif


