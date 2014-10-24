// ********************************************************************
//
// $Id: HRSDCSD.cc,v 1.0, 2010/12/26 HRS Exp $
// --------------------------------------------------------------
//
//By Jixie: The DC SD is just for position detector, good for DC and SC. 
//It can also be used to measure energy with complicate analysis.
//It will do the following: 
//1) It will create one hit for each interaction 
//  (one hit one interaction, therefore multiple hits for each track)
//2) Record the inpos and inmom from the prestep of current interaction, and
//   record the time, outpos and outmom from the poststep of current interaction 
//3) Sign deposited energy to non-ionized energy if it is a transportation or MSC   
//In the event action, these hits will be stored into the root ntuple


#include "HRSDCSD.hh"
#include "HRSDCHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4VProcess.hh"
#include "G4ios.hh"

HRSDCSD::HRSDCSD(G4String name):G4VSensitiveDetector(name)
{
	G4String HCname;
	collectionName.insert(HCname="DCColl");
	HCID = -1;

	//The given 'name' will be stored into G4VSensitiveDetector::SensitiveDetectorName
	//you can give '/RTPCGEM' or 'RTPCGEM', the Geant4 system will both set 'RTPCGEM'
	//as the SensitiveDetectorName and '/RTPCGEM' as fullPathName
	//
	//For current hit colection, the collectionName is "SensitiveDetectorName/HCname"
	//To access this HC in event action, one need to get HCID through
	//G4HCtable::GetCollectionID(G4String collectionName) 

	hitsCollection = 0;
	verboseLevel = 0; //defined in G4VSensitiveDetector.hh
}

HRSDCSD::~HRSDCSD()
{
	;
}

void HRSDCSD::Initialize(G4HCofThisEvent* HCE)
{
	hitsCollection = new HRSDCHitsCollection(SensitiveDetectorName,collectionName[0]);
	if(HCID<0)
	{ 
		HCID = G4SDManager::GetSDMpointer()->GetCollectionID(hitsCollection); 
	}
	HCE->AddHitsCollection(HCID,hitsCollection);
}


G4bool HRSDCSD::ProcessHits(G4Step* aStep,G4TouchableHistory* /*aROHist*/)
{	
	if(!hitsCollection) return false;

	G4double edep = aStep->GetTotalEnergyDeposit();
	if(edep/keV<=0.)  return true;


	G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
	G4TouchableHistory* theTouchable = (G4TouchableHistory*)(preStepPoint->GetTouchable());
	G4VPhysicalVolume* thePhysical = theTouchable->GetVolume();
	
	G4int copyNo = thePhysical->GetCopyNo();
	G4int trackId = aStep->GetTrack()->GetTrackID();
	G4double hitTime = preStepPoint->GetGlobalTime();

	G4StepPoint* postStepPoint = aStep->GetPostStepPoint();
	G4ThreeVector outpos = postStepPoint->GetPosition();       
	G4ThreeVector outmom = postStepPoint->GetMomentum();

	//By Jixie; it turns out that this is always zero due to the fact that G4 physics process never set 
	//a value to it
	G4double edep_NonIon = aStep->GetNonIonizingEnergyDeposit();
	G4String process=postStepPoint->GetProcessDefinedStep()->GetProcessName();
	if(process=="Transportation" || process=="msc") edep_NonIon=edep;	
	//if(edep_NonIon) G4cout<<"edep_NonIon="<<edep_NonIon<<G4endl;

	
	// generally a hit is the total energy deposit in a sensitive detector within a time window
	// it does not matter which particle deposite the energy
	// For drift chamber, I only care the time and position, and also I need to know 
	// which track create this hit, therefore I also store parent trackid

	HRSDCHit* aHit = new HRSDCHit(copyNo,hitTime);

	G4int parentTrackId=aStep->GetTrack()->GetParentID();
	G4ThreeVector inpos = preStepPoint->GetPosition();       
	G4ThreeVector inmom = preStepPoint->GetMomentum();

	G4int pdgid=aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
	if(pdgid>3000)
	{
		//G4cout<<"Wrong PDGCCoding pdgid="<<pdgid
		//	<<" name="<<aStep->GetTrack()->GetParticleDefinition()->GetParticleName()<<G4endl;
		pdgid=(pdgid%1000000)/10;
	}

	aHit->SetPhysV(theTouchable->GetVolume());
	aHit->SetEdep(edep);
	aHit->SetNonIonEdep(edep_NonIon);
	aHit->SetTrackId(trackId);
	aHit->SetOutPos(outpos);
	aHit->SetOutMom(outmom);
	aHit->SetInPos(inpos);
	aHit->SetInMom(inmom);
	aHit->SetParentTrackId(parentTrackId);
	aHit->SetPdgid(pdgid);

	hitsCollection->insert(aHit);

	if(verboseLevel >= 2)
	{
		G4cout<<"<<<Create a new DC Hit of type <" << SensitiveDetectorName << "/"<<collectionName[0]
		<<" in "<<theTouchable->GetVolume()->GetName()<<"["<<copyNo<<"] \n"
			<<" at "<<inpos<<", by track "<< trackId <<", momentum="<<inmom<<G4endl;
	}

	return true;
}

void HRSDCSD::EndOfEvent(G4HCofThisEvent* HCE)
{
	if(!HCE) return;  //no hits collection found 
	//the above line just to avoid warning of not use HCE

	int nhitC = hitsCollection->GetSize();
	if(!nhitC) return;
	
	if(verboseLevel >= 2 && nhitC)
	{
		G4cout<<"<<<End of Hit Collections <" << SensitiveDetectorName << "/"<<collectionName[0]
		<<">: " << nhitC << " hits." << G4endl; 
		for (int i=0; i<nhitC; i++)
		{
			(*hitsCollection)[i]->Print();
		}
	}

	return;
}

