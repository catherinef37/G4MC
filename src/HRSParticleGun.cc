// ********************************************************************
// G4ParticleGun always print out some warning message while seting the
// energy or momentum of the primary particle, thie really slow down the
// speed and very anoying. I wrote this class to override it
//
// $Id: HRSParticleGun.cc,v 1.0 2011/2/11 HRS Exp $
// GEANT4 tag $Name: geant4-09-04 $
// ********************************************************************


// HRSParticleGun
#include "HRSParticleGun.hh"
#include "G4PrimaryParticle.hh"
#include "G4Event.hh"
#include "G4ios.hh"

HRSParticleGun::HRSParticleGun()
{
	SetInitialValues();
}

HRSParticleGun::HRSParticleGun(G4int numberofparticles)
{
	SetInitialValues();
	NumberOfParticlesToBeGenerated = numberofparticles;
}

HRSParticleGun::HRSParticleGun(G4ParticleDefinition * particleDef, G4int numberofparticles)
{
	SetInitialValues();
	NumberOfParticlesToBeGenerated = numberofparticles;
	SetParticleDefinition( particleDef );
}

HRSParticleGun::~HRSParticleGun()
{
	;
}


void HRSParticleGun::SetParticleEnergy(G4double aKineticEnergy)
{
	particle_energy = aKineticEnergy;
}

void HRSParticleGun::SetParticleMomentum(G4double aMomentum)
{	
	if(particle_definition==0)
	{
		G4cout <<"Particle Definition not defined yet for HRSParticleGun"<< G4endl;
		G4cout <<"Zero Mass is assumed"<<G4endl;
		particle_momentum = aMomentum;
		particle_energy = aMomentum;
	}
	else
	{
		G4double mass =  particle_definition->GetPDGMass();
		particle_momentum = aMomentum;
		particle_energy = std::sqrt(particle_momentum*particle_momentum+mass*mass)-mass;
	}
}

void HRSParticleGun::SetParticleMomentum(G4ParticleMomentum aMomentum)
{
	if(particle_definition==0)
	{
		G4cout <<"Particle Definition not defined yet for HRSParticleGun"<< G4endl;
		G4cout <<"Zero Mass is assumed"<<G4endl;
		particle_momentum_direction =  aMomentum.unit();
		particle_momentum = aMomentum.mag();
		particle_energy = aMomentum.mag();
	} 
	else 
	{
		G4double mass =  particle_definition->GetPDGMass();
		particle_momentum = aMomentum.mag();
		particle_momentum_direction =  aMomentum.unit();
		particle_energy = std::sqrt(particle_momentum*particle_momentum+mass*mass)-mass;
	}
}

void HRSParticleGun::GeneratePrimaryVertex(G4Event* evt)
{
	if(particle_definition==0) return;

	// create a new vertex
	G4PrimaryVertex* vertex = new G4PrimaryVertex(particle_position,particle_time);

	// create new primaries and set them to the vertex
	G4double mass =  particle_definition->GetPDGMass();
	G4double energy = particle_energy + mass;
	G4double pmom = std::sqrt(energy*energy-mass*mass);
	G4double px = pmom*particle_momentum_direction.x();
	G4double py = pmom*particle_momentum_direction.y();
	G4double pz = pmom*particle_momentum_direction.z();
	for( G4int i=0; i<NumberOfParticlesToBeGenerated; i++ )
	{
		G4PrimaryParticle* particle = new G4PrimaryParticle(particle_definition,px,py,pz);
		particle->SetMass( mass );
		particle->SetCharge( particle_charge );
		particle->SetPolarization(particle_polarization.x(),
			particle_polarization.y(),particle_polarization.z());
		vertex->SetPrimary( particle );
	}

	evt->AddPrimaryVertex( vertex );
}

