/***********************************************

Used to implement traking in the
scintillators and momentum "measurement".

It the first part we are dealing with
sensitive detectors, and in the second part
with the 'hits collection' system.

***********************************************/

#include "Randomize.hh"
#include "BigBiteHit.hh"
G4Allocator<BigBiteHit> ndhitsallocator;//can comment it--Ramesh

//## 1: Hits collection and utilisation

// Nothing to do here

BigBiteHit::BigBiteHit(){}
BigBiteHit::~BigBiteHit(){}
//void BigBiteHit::draw(G4double xpos,G4double ypos,G4double zpos){}
//void BigBiteHit::Print(){}
//## 2: Sensitive detectors members

BigBiteSD::BigBiteSD(G4String sdname) : G4VSensitiveDetector(sdname) 
{
	param = sdname;
	collectionName.insert(param);  // needed to give a name to collection
	collid = 1;  // so only 1 pointer per detector is created
}

BigBiteSD::~BigBiteSD()
{;}

// ProcessHits takes the values associated to the hits, both for particle
// and the detector, and allow us to save them in an BigBiteHit class.
G4bool BigBiteSD::ProcessHits(G4Step* step,G4TouchableHistory* touchist)
{
	G4double c = 300000.*km/s/1.58;  // in mm/ns

	G4Track* thetrack = step->GetTrack();
	G4double steptime = thetrack->GetGlobalTime();
	// G4double steptime = thetrack->GetGlobalTime()-50.*ns;
	if (steptime < 0.) { return false; }

	// 'creation' of the hits
	BigBiteHit* anotherhit = new BigBiteHit();

	anotherhit->hitime = steptime;
	// G4cout << "\nsteptime is: " <<steptime <<G4endl;
	anotherhit->posofhit = step->GetPreStepPoint()->GetPosition();
	anotherhit->kinenergy = thetrack->GetKineticEnergy();
	// rem : momentum not accessible
	anotherhit->kineloss = step->GetTotalEnergyDeposit();
	// G4cout << "\nKinetic energy is: "<< anotherhit->kineloss<<G4endl;

	// G4double xhpos = anotherhit->posofhit[0];//nowhere used--Ramesh
	G4double yhpos = anotherhit->posofhit[1];
	// G4double zhpos = anotherhit->posofhit[2];//nowhere used--Ramesh


	// Get the replica number
	G4int repnum;
	touchable = (G4TouchableHistory*)(step->GetPreStepPoint()->GetTouchable());
	repnum = touchable->GetReplicaNumber() + 1;
	// if it is not a replica

	//  G4cout << "\n @@ " << repnum << param << "\n" << flush;
	anotherhit->barnumber = repnum;

	// Set the 'measured' time
	// no smearing for the moment
	if (param == "auxbar")
	{anotherhit->ltime = (185.+yhpos)/c + steptime;
	anotherhit->rtime = (185.-yhpos)/c + steptime;}
	if (param == "delta" || param == "mainE" )
	{anotherhit->ltime = (250.+yhpos)/c + steptime;
	anotherhit->rtime = (250.-yhpos)/c + steptime;}
	else 
	{ 
		if (param !="auxbar") 
		{anotherhit->ltime = (500.+yhpos)/c + steptime;
		anotherhit->rtime = (500.-yhpos)/c + steptime;
		// G4cout << "\nltime is: "<<anotherhit->ltime<<G4endl;
		// G4cout << "\nrtime is: "<<anotherhit->rtime<<G4endl;
		}
	}


	anotherhit->lmestime = G4RandGauss::shoot(0.,0.5*ns) + anotherhit->ltime;
	anotherhit->rmestime = G4RandGauss::shoot(0.,0.5*ns) + anotherhit->rtime;



	if ((param == "wirec1") || (param == "wirec2"))
	{ // nothing better for the moment
		anotherhit->ltime = steptime;
		anotherhit->rtime = steptime;
		anotherhit->lmestime = steptime;
		anotherhit->rmestime = steptime;
	}

	//  //Ran.S added accidentals statistics on 03/03
	//   if ( G4UniformRand()>0.5 )
	//   {
	//    anotherhit->lmestime = G4RandGauss::shoot(0.,15.*ns);
	//    anotherhit->rmestime = G4RandGauss::shoot(0.,15.*ns);   
	//   }



	pointhitcoll->insert(anotherhit);

	return true;
}

// Initialize initialize the collections of events for each detector
//void BigBiteSD::Initialize(G4HCofThisEvent* hcote)//Ramesh
void BigBiteSD::Initialize(G4HCofThisEvent* hcte)//Ramesh
{
	pointhitcoll = new ndhitscollection(SensitiveDetectorName, param);
	if (collid == 1)
	{ collid = G4SDManager::GetSDMpointer()->GetCollectionID(param); }
}

// EndOfEvent adds the values collected during event to the collections.

void BigBiteSD::EndOfEvent(G4HCofThisEvent* hcte)
{

	hcte->AddHitsCollection(collid,pointhitcoll);
}

//[Part of the code of J.Annand, not made by me. PM]

// This is a forward declarations of an instantiated G4Allocator<Type> object.
// It has been added in order to make code portable for the GNU g++ 
// (release 2.7.2) compiler. 
// Whenever a new Type is instantiated via G4Allocator, it has to be forward
// declared to make object code (compiled with GNU g++) link successfully. 
// 
#ifdef GNU_GCC//can comment it--Ramesh
template class G4Allocator<BigBiteHit>;//can comment it--Ramesh
#endif//can comment it--Ramesh
