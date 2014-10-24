// ********************************************************************
//
// $Id: HRSSteppingVerbose.cc,v 1.0, 2010/12/26   HRS Exp $
//
// ********************************************************************
//..............................................................................

#include "HRSSteppingVerbose.hh"

#include "G4SteppingManager.hh"
#include "G4UnitsTable.hh"

//..............................................................................

HRSSteppingVerbose::HRSSteppingVerbose()
{}

//..............................................................................

HRSSteppingVerbose::~HRSSteppingVerbose()
{}

//..............................................................................

void HRSSteppingVerbose::StepInfo()
{
	if(fTrack->GetCurrentStepNumber()>=1024)
	{
		fTrack->SetTrackStatus(fStopAndKill);
		return;
	}
	CopyState();

	// G4int prec = G4cout.precision(3);
	G4int prec = G4cout.precision(5);

	if( verboseLevel >= 1 ){
		if( verboseLevel >= 4 ) VerboseTrack();
		if( verboseLevel >= 3 ){
			G4cout << G4endl;

			G4cout << std::setw( 5) << "Step#"   << " "
				<< std::setw( 8) << "X(mm)"      << " "
				<< std::setw( 8) << "Y(mm)"      << " "
				<< std::setw( 8) << "Z(mm)"      << " "
				<< std::setw( 8) << "KineE(MeV)" << " "
				<< std::setw(11) << "dE(KeV)"    << " "
				<< std::setw(11) << "StepL(mm)"  << " "
				<< std::setw(9)  << "TrackL(mm)" << " "
				<< std::setw(22) << "NextVolume" << " "
				<< std::setw(14) << "ProcessName"<< G4endl;

		}

		G4cout << std::setw(5) << fTrack->GetCurrentStepNumber() << " "
			<< std::setw(8) << fTrack->GetPosition().x()/mm << " "
			<< std::setw(8) << fTrack->GetPosition().y()/mm << " "
			<< std::setw(8) << fTrack->GetPosition().z()/mm << " "
			<< std::setw(8) << fTrack->GetKineticEnergy()/MeV << " "
			<< std::setw(11)<< fStep->GetTotalEnergyDeposit()/keV << " "
			<< std::setw(11)<< fStep->GetStepLength()/mm << " "
			<< std::setw(9) << fTrack->GetTrackLength()/mm << " ";

		if( fTrack->GetNextVolume() != 0 ) {
			G4cout << std::setw(22) << fTrack->GetVolume()->GetName()<< " ";
		} else {
			G4cout << std::setw(22) << "OutOfWorld"<< " ";
		}

		if(fStep->GetPostStepPoint()->GetProcessDefinedStep() != NULL){
			G4cout << std::setw(14) 
				<< fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
		} else {
			G4cout <<  std::setw(14) << "UserLimit";
		}
		G4cout << G4endl;

		if( verboseLevel == 2 ){
			G4int tN2ndariesTot = fN2ndariesAtRestDoIt +
				fN2ndariesAlongStepDoIt +
				fN2ndariesPostStepDoIt;
			if(tN2ndariesTot>0){
				G4cout << "    :----- List of 2ndaries - "
					<< "#SpawnInStep=" << std::setw(3) << tN2ndariesTot
					<< "(Rest="  << std::setw(2) << fN2ndariesAtRestDoIt
					<< ",Along=" << std::setw(2) << fN2ndariesAlongStepDoIt
					<< ",Post="  << std::setw(2) << fN2ndariesPostStepDoIt
					<< "), "
					<< "#SpawnTotal=" << std::setw(3) << (*fSecondary).size()
					<< " ---------------"
					<< G4endl;

				for(size_t lp1=(*fSecondary).size()-tN2ndariesTot;
					lp1<(*fSecondary).size(); lp1++){
						G4cout << "    : "
							<< std::setw(8)
							<< G4BestUnit((*fSecondary)[lp1]->GetPosition().x(),"Length")
							<< std::setw(8)
							<< G4BestUnit((*fSecondary)[lp1]->GetPosition().y(),"Length")
							<< std::setw(8)
							<< G4BestUnit((*fSecondary)[lp1]->GetPosition().z(),"Length")
							<< std::setw(8)
							<< G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(),"Energy")
							<< std::setw(10)
							<< (*fSecondary)[lp1]->GetDefinition()->GetParticleName();
						G4cout << G4endl;
				}

				G4cout << "    :-----------------------------"
					<< "----------------------------------"
					<< "-- EndOf2ndaries Info ---------------"
					<< G4endl;
			}
		}

	}
	G4cout.precision(prec);

	if(fTrack->GetCurrentStepNumber()>=1024)
	{
		fTrack->SetTrackStatus(fStopAndKill);
		return;
	}
}

//..............................................................................

void HRSSteppingVerbose::TrackingStarted()
{
	CopyState();
	G4int prec = G4cout.precision(5);
	if( verboseLevel > 0 ){

		G4cout << std::setw( 5) << "Step#"   << " "
			<< std::setw( 8) << "X(mm)"      << " "
			<< std::setw( 8) << "Y(mm)"      << " "
			<< std::setw( 8) << "Z(mm)"      << " "
			<< std::setw( 8) << "KineE(MeV)" << " "
			<< std::setw(11) << "dE(KeV)"    << " "
			<< std::setw(11) << "StepL(mm)"  << " "
			<< std::setw(9)  << "TrackL(mm)" << " "
			<< std::setw(22) << "NextVolume" << " "
			<< std::setw(14) << "ProcessName"<< G4endl;

		G4cout << std::setw(5) << fTrack->GetCurrentStepNumber() << " "
			<< std::setw(8) << fTrack->GetPosition().x()/mm << " "
			<< std::setw(8) << fTrack->GetPosition().y()/mm << " "
			<< std::setw(8) << fTrack->GetPosition().z()/mm << " "
			<< std::setw(8) << fTrack->GetKineticEnergy()/mm << " "
			<< std::setw(11)<< fStep->GetTotalEnergyDeposit()/keV << " "
			<< std::setw(11)<< fStep->GetStepLength()/mm << " "
			<< std::setw(9) << fTrack->GetTrackLength()/mm << " ";

		if(fTrack->GetNextVolume()){
			G4cout << std::setw(22) << fTrack->GetVolume()->GetName() << " ";
		} else {
			G4cout << std::setw(22) << "OutOfWorld" << " ";
		}
		G4cout  <<  std::setw(14) <<"initStep" << G4endl;
	}
	G4cout.precision(prec);

}

//..............................................................................

