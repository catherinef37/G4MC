// ********************************************************************
//
// $Id: HRSTrackInformation.cc,v 1.0, 2010/12/26 HRS Exp $
// --------------------------------------------------------------
//

#include "HRSTrackInformation.hh"
#include "G4ios.hh"

G4Allocator<HRSTrackInformation> aTrackInformationAllocator;

HRSTrackInformation::HRSTrackInformation()
{
	originalTrackID = 0;
	originalDefinition = 0;
	originalPosition = G4ThreeVector(0.,0.,0.);
	originalMomentum = G4ThreeVector(0.,0.,0.);
	originalEnergy = 0.;
	originalTime = 0.;
	motherTrackID = 0;
	iStepOutputStatus=1;
}

HRSTrackInformation::HRSTrackInformation(const G4Track* aTrack)
{
	originalTrackID = aTrack->GetTrackID();
	originalDefinition = aTrack->GetDefinition();
	originalPosition = aTrack->GetPosition();
	originalMomentum = aTrack->GetMomentum();
	originalEnergy = aTrack->GetTotalEnergy();
	originalTime = aTrack->GetGlobalTime();
	motherTrackID = aTrack->GetParentID();
}

HRSTrackInformation::HRSTrackInformation(const HRSTrackInformation* aTrackInfo)
{
	originalTrackID = aTrackInfo->GetOriginalTrackID();
	originalDefinition = aTrackInfo->GetOriginalDefinition();
	originalPosition = aTrackInfo->GetOriginalPosition();
	originalMomentum = aTrackInfo->GetOriginalMomentum();
	originalEnergy = aTrackInfo->GetOriginalEnergy();
	originalTime = aTrackInfo->GetOriginalTime();
	motherTrackID = aTrackInfo->GetMotherTrackID();
}

HRSTrackInformation::~HRSTrackInformation(){;}

void HRSTrackInformation::Print() const
{
	G4cout <<"At time"<< originalTime
		<< "Original track ID " << originalTrackID
		<< " at " << originalPosition
		<< "with Momentum "<< originalMomentum<<G4endl;
}

void HRSTrackInformation::Print(const G4Track* aTrack) const
{
	//print the following info
	//*********************************************************************************************
	//>>> Simulation truth:   Theta=6.000deg    Phi=-90.000deg    Position=(10.000,0.000,15.000mm)
	//>>> Initial Momentum=(-0.000,-99.302,944.796)=950.000MeV    Particle=e-
	//*********************************************************************************************
	//cout.precision Return Value:  The format flags of the stream before the call.
	G4int prec = G4cout.precision(3);
	G4cout<< "\n*********************************************************************************************\n";
	G4cout<<">>> Simulation truth:"
		<< "   Theta=" << aTrack->GetMomentum().theta()/deg<< "deg "
		<< "   Phi=" << aTrack->GetMomentum().phi()/deg<< "deg "
		<< std::fixed<< "   Position=" << aTrack->GetPosition()<< "mm"
		<< std::fixed<< "\n>>> Initial Momentum=" << aTrack->GetMomentum()
		<< "="<< std::fixed<< aTrack->GetMomentum().mag()<< "MeV "
		<< "   Particle=" << aTrack->GetDefinition()->GetParticleName()
		<< G4endl;
	G4cout<< "***********************************************************************************************\n";
	G4cout.precision(prec);
}
