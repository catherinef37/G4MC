// ********************************************************************
//
// $Id: HRSTrajectory.cc,v 1.0, 2010/12/26 HRS Exp $
// --------------------------------------------------------------
//

#include "HRSTrajectory.hh"
#include "G4TrajectoryPoint.hh"
#include "G4ThreeVector.hh"
#include "G4Polyline.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4VVisManager.hh"
#include "G4UnitsTable.hh"
#include "G4UIcommand.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include <iomanip>

G4Allocator<HRSTrajectory> myTrajectoryAllocator;

HRSTrajectory::HRSTrajectory()
{
	particleDefinition = 0;
	particleName = "";
	PDGCharge = 0;
	PDGEncoding = 0;
	trackID = 0;
	parentTrackID = 0;
	positionRecord = 0;
	initialKineticEnergy = 0;
	initialMomentum.set(0.,0.,0.);
	initialTime = 0;
}

HRSTrajectory::HRSTrajectory(const G4Track* aTrack)
{
	particleDefinition = aTrack->GetDefinition();
	particleName = particleDefinition->GetParticleName();
	PDGCharge = particleDefinition->GetPDGCharge();
	PDGEncoding = particleDefinition->GetPDGEncoding();
	trackID = aTrack->GetTrackID();
	parentTrackID = aTrack->GetParentID();
	positionRecord = new HRSTrajectoryPointContainer();
	positionRecord->push_back(new G4TrajectoryPoint(aTrack->GetPosition()));
	initialKineticEnergy = aTrack->GetKineticEnergy(); 
	initialMomentum = aTrack->GetMomentum(); 
	initialTime = aTrack->GetGlobalTime(); 
}

HRSTrajectory::HRSTrajectory(HRSTrajectory & right)
: G4VTrajectory()
{
	particleName = right.particleName;
	particleDefinition = right.particleDefinition;
	PDGCharge = right.PDGCharge;
	PDGEncoding = right.PDGEncoding;
	trackID = right.trackID;
	parentTrackID = right.parentTrackID;
	positionRecord = new HRSTrajectoryPointContainer();
	for(int i=0;i<(int)right.positionRecord->size();i++)
	{
		G4TrajectoryPoint* rightPoint = (G4TrajectoryPoint*)((*(right.positionRecord))[i]);
		positionRecord->push_back(new G4TrajectoryPoint(*rightPoint));
	}
	initialKineticEnergy = right.initialKineticEnergy;
	initialMomentum = right.initialMomentum;
	initialTime = right.initialTime; 
}

HRSTrajectory::~HRSTrajectory()
{
	for(size_t i=0;i<positionRecord->size();i++)
	{
		delete  (*positionRecord)[i];
	}
	positionRecord->clear();

	delete positionRecord;
}

void HRSTrajectory::ShowTrajectory() const
{
	// Invoke the default implementation in G4VTrajectory...
	//G4VTrajectory::ShowTrajectory();
	ShowTrajectory(G4cout);
}

void HRSTrajectory::ShowTrajectory(std::ostream& out) const
{
	out << G4endl << "TrackID = " << trackID << "  ParentID = " << parentTrackID << std::endl;
	out << "ParticleName = " << particleName << "  PDGCode = " << PDGEncoding
		<< "  Charge = " << PDGCharge << std::endl;
	out << "OriginalMomentum = " << initialMomentum.mag()/GeV <<" GeV"<< std::endl;
	out << "Current trajectory has " << positionRecord->size() << " steps." << std::endl;

	size_t precision=out.precision(2);
	out.setf(std::ios_base::fixed); 
	G4ThreeVector aPoint;
	for( size_t i=0 ; i < positionRecord->size() ; i++)
	{
		aPoint = ((G4TrajectoryPoint*)((*positionRecord)[i]))->GetPosition();
		out << "Step " <<std::setw(4)<< i << ": "
			<<std::setw(10)<< aPoint.x()/mm << ", " 
			<<std::setw(10)<< aPoint.y()/mm << ", " 
			<<std::setw(10)<< aPoint.z()/mm << " mm" 
			<< std::endl;
	}
	out.precision(precision);
	out<<std::endl;
}


void HRSTrajectory::DrawTrajectory(G4int /*i_mode*/) const
{
	// Invoke the default implementation in G4VTrajectory...
	// G4VTrajectory::DrawTrajectory();
	// ... or override with your own code here.

	G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
	G4ThreeVector pos;

	G4Polyline pPolyline;
	for (int i = 0; i < (int)positionRecord->size() ; i++) 
	{
		G4TrajectoryPoint* aTrajectoryPoint = (G4TrajectoryPoint*)((*positionRecord)[i]);
		pos = aTrajectoryPoint->GetPosition();
		pPolyline.push_back( pos );
	}

	G4Colour colour(0.75,0.75,0.75);    // LightGray
	if(particleDefinition==G4Gamma::GammaDefinition())
		colour = G4Colour(0.,1.,1.);     // Cyan
	else if(particleDefinition==G4Electron::ElectronDefinition() ||
		particleDefinition==G4Positron::PositronDefinition())
		colour = G4Colour(1.,1.,0.);      // Yellow
	else if(particleDefinition==G4MuonMinus::MuonMinusDefinition() || 
		particleDefinition==G4MuonPlus::MuonPlusDefinition())
		colour = G4Colour(1.,0.,1.);      // Magenta
	else if(particleDefinition->GetParticleType()=="meson")
	{
		if(PDGCharge!=0.)
			colour = G4Colour(1.,0.,0.);   // Red
		else
			colour = G4Colour(0.5,0.,0.);  // HalfRed
	}
	else if(particleDefinition->GetParticleType()=="baryon")
	{
		if(PDGCharge!=0.)
			colour = G4Colour(1.,0.78,0.); // Orange
		else
			colour = G4Colour(0.5,0.39,0.);// HalfOrange
	}

	G4VisAttributes attribs(colour);
	pPolyline.SetVisAttributes(attribs);
	if(pVVisManager) pVVisManager->Draw(pPolyline);
}

const std::map<G4String,G4AttDef>* HRSTrajectory::GetAttDefs() const
{
	G4bool isNew;
	std::map<G4String,G4AttDef>* store
		= G4AttDefStore::GetInstance("HRSTrajectory",isNew);
	if (isNew) 
	{
		G4String ID("ID");
		(*store)[ID] = G4AttDef(ID,"Track ID","Physics","","G4int");

		G4String PID("PID");
		(*store)[PID] = G4AttDef(PID,"Parent Track ID","Physics","","G4int");

		G4String PN("PN");
		(*store)[PN] = G4AttDef(PN,"Particle Name","Physics","","G4String");

		G4String Ch("Ch");
		(*store)[Ch] = G4AttDef(Ch,"Charge","Physics","e+","G4double");

		G4String PDG("PDG");
		(*store)[PDG] = G4AttDef(PDG,"PDG Encoding","Physics","","G4int");

		G4String IKE("IKE");
		(*store)[IKE] = 
			G4AttDef(IKE, "Initial kinetic energy",
			"Physics","G4BestUnit","G4double");

		G4String IMom("IMom");
		(*store)[IMom] = G4AttDef(IMom, "Initial momentum",
			"Physics","G4BestUnit","G4ThreeVector");

		G4String IMag("IMag");
		(*store)[IMag] = 
			G4AttDef(IMag, "Initial momentum magnitude",
			"Physics","G4BestUnit","G4double");

		G4String NTP("NTP");
		(*store)[NTP] = G4AttDef(NTP,"No. of points","Physics","","G4int");
	}
	return store;
}

std::vector<G4AttValue>* HRSTrajectory::CreateAttValues() const
{
	std::vector<G4AttValue>* values = new std::vector<G4AttValue>;

	values->push_back
		(G4AttValue("ID",G4UIcommand::ConvertToString(trackID),""));

	values->push_back
		(G4AttValue("PID",G4UIcommand::ConvertToString(parentTrackID),""));

	values->push_back(G4AttValue("PN",particleName,""));

	values->push_back
		(G4AttValue("Ch",G4UIcommand::ConvertToString(PDGCharge),""));

	values->push_back
		(G4AttValue("PDG",G4UIcommand::ConvertToString(PDGEncoding),""));

	values->push_back
		(G4AttValue("IKE",G4BestUnit(initialKineticEnergy,"Energy"),""));

	values->push_back
		(G4AttValue("IMom",G4BestUnit(initialMomentum,"Energy"),""));

	values->push_back
		(G4AttValue("IMag",G4BestUnit(initialMomentum.mag(),"Energy"),""));

	values->push_back
		(G4AttValue("NTP",G4UIcommand::ConvertToString(GetPointEntries()),""));

#ifdef G4ATTDEBUG
	G4cout << G4AttCheck(values,GetAttDefs());
#endif

	return values;
}

void HRSTrajectory::AppendStep(const G4Step* aStep)
{
	positionRecord->push_back( 
		new G4TrajectoryPoint(aStep->GetPostStepPoint()->GetPosition()) );
}

void HRSTrajectory::MergeTrajectory(G4VTrajectory* secondTrajectory)
{
	if(!secondTrajectory) return;

	HRSTrajectory* seco = (HRSTrajectory*)secondTrajectory;
	G4int ent = seco->GetPointEntries();
	for(int i=1;i<ent;i++) // initial point of the second trajectory should not be merged
	{
		positionRecord->push_back((*(seco->positionRecord))[i]);
	}
	delete (*seco->positionRecord)[0];
	seco->positionRecord->clear();
}


