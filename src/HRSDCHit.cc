// ********************************************************************
//
// $Id: HRSDCHit.cc,v 1.0, 2010/12/26 HRS Exp $
// --------------------------------------------------------------
//DC hit does not care energy deposited(very small), just time and position,
//therefore no need to accumulate the deposited energy in SD

#include "HRSDCHit.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"

#include <iomanip>
#include <iostream>
using namespace std;

G4Allocator<HRSDCHit> HRSDCHitAllocator;

HRSDCHit::HRSDCHit(G4int i,G4double t)
{
	id = i;
	pdgid=0;
	time = t;
	edep=edep_NonIon=0;
	trackid=parenttrackid=0;
	inpos.set(0,0,0); 
	outpos.set(0,0,0); 
	inmom.set(0,0,0); 
	outmom.set(0,0,0); 
	physV = 0;
}

HRSDCHit::~HRSDCHit()
{;}

HRSDCHit::HRSDCHit(const HRSDCHit &right)
: G4VHit() 
{
	id = right.id;
	pdgid=right.pdgid; 
	edep = right.edep;
	edep_NonIon=right.edep_NonIon;
	time = right.time;
	inpos = right.inpos;
	outpos = right.outpos;
	inmom = right.inmom;
	outmom = right.outmom;
	trackid = right.trackid;
	parenttrackid = right.parenttrackid;
	physV = right.physV;
}

const HRSDCHit& HRSDCHit::operator=(const HRSDCHit &right)
{
	id = right.id;
	pdgid=right.pdgid; 
	edep = right.edep;
	edep_NonIon=right.edep_NonIon;
	time = right.time;
	inpos = right.inpos;
	outpos = right.outpos;
	inmom = right.inmom;
	outmom = right.outmom;
	trackid = right.trackid;
	parenttrackid = right.parenttrackid;
	physV = right.physV;
	return *this;
}

int HRSDCHit::operator==(const HRSDCHit &right) const
{
	return ((id == right.id) && (time == right.time) &&  
		(trackid==right.trackid) && (physV == right.physV));
}

void HRSDCHit::Draw()
{
	G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
	if(pVVisManager)
	{
		//just draw a circle with size 3 (mm?)
		G4Circle circle(inpos);
		circle.SetScreenSize(3);
		circle.SetFillStyle(G4Circle::filled);
		G4Colour colour(1.,1.,0.);
		G4VisAttributes attribs(colour);
		circle.SetVisAttributes(attribs);
		pVVisManager->Draw(circle);     
	}
}

void HRSDCHit::Print()
{
	int prec=G4cout.precision(4);
	G4cout.setf(ios::fixed);
	G4cout<<std::setw(20)<<physV->GetName()<<"[" << id << "] was hitted by track "<<trackid 
		<<" (parent_track="<<parenttrackid<<", PdgID="<<pdgid
		<<"), at time="<< time/ns << " (nsec)" <<"\n "
		<<"\t entrance momentum=" << inmom <<inmom.mag()/MeV<<"MeV \n"
		<<"\t outgoing momentum=" << outmom <<outmom.mag()/MeV<<"MeV \n"
		<<"\t Edep=" << edep/MeV  <<"MeV NonIon_Edep=" << edep_NonIon/MeV  
		<< G4endl;
	G4cout.unsetf(ios::fixed);
	G4cout.precision(prec);
}


