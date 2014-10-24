#include "BigBiteEventAction.hh"
#include "BigBiteHit.hh"

#include "math.h"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4Square.hh"
#include "Randomize.hh"

//Includes needed for output file
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;


/************************************************/
// constructor for Tree generation :
BigBiteEventAction:: BigBiteEventAction(TTree* ndtree, char z)
{ a = z;
// here we create the branches to be filled
// The name of each branches explains what is the variable for

tnd =  ndtree;
/*
//these 4 are from primary generator action, I want to change it
tnd->Branch("event_number",&eventn,"en/I");
tnd->Branch("in_plane_angle",&theta,"theta/D");
tnd->Branch("out_of_plane_angle",&phi,"phi/D");
tnd->Branch("momentum",&mestpp,"mestpp/D");
*/

tnd->Branch("aux_tdc_time",auxhtimeta,"ta[56][6]/D");
tnd->Branch("energy_deposita",ekinapta,"kina[56]/D"); 


tnd->Branch("r_tdc_valuede",rmestimede,"rtde[24][6]/D");
tnd->Branch("l_tdc_valuede",lmestimede,"ltde[24][6]/D");
tnd->Branch("energy_depositde",ekinaptde,"ekinde[24]/D");

tnd->Branch("r_tdc_valuee",rmestimee,"rte[24][6]/D");
tnd->Branch("l_tdc_valuee",lmestimee,"lte[24][6]/D");
tnd->Branch("energy_deposite",ekinapte,"ekine[24]/D");


tnd->Branch("r_tdc_valuep0",rmestimep0,"rtp0[30][6]/D");
tnd->Branch("l_tdc_valuep0",lmestimep0,"ltp0[30][6]/D");
tnd->Branch("energy_depositp0",ekinaptp0,"ekinp0[30]/D");


tnd->Branch("r_tdc_valuep1",rmestimep1,"rtp1[30][6]/D");
tnd->Branch("l_tdc_valuep1",lmestimep1,"ltp1[30][6]/D");
tnd->Branch("energy_depositp1",ekinaptp1,"ekinp1[30]/D");

tnd->Branch("r_tdc_valuep2",rmestimep2,"rtp2[24][6]/D");
tnd->Branch("l_tdc_valuep2",lmestimep2,"ltp2[24][6]/D");
tnd->Branch("energy_depositp2",ekinaptp2,"ekinp2[24]/D");

tnd->Branch("r_tdc_valuep3",rmestimep3,"rtp3[6][6]/D");
tnd->Branch("l_tdc_valuep3",lmestimep3,"ltp3[6][6]/D");
tnd->Branch("energy_depositp3",ekinaptp3,"ekinp3[6]/D");

tnd->Branch("r_tdc_valuep4",rmestimep4,"rtp4[4][6]/D");
tnd->Branch("l_tdc_valuep4",lmestimep4,"ltp4[4][6]/D");  
tnd->Branch("energy_depositp4",ekinaptp4,"ekinp4[4]/D");

tnd->Branch("r_tdc_valuep5",rmestimep5,"rtp5[2][6]/D");
tnd->Branch("l_tdc_valuep5",lmestimep5,"ltp5[2][6]/D");
tnd->Branch("energy_depositp5",ekinaptp5,"ekinp5[2]/D");


tnd->Branch("r_tdc_valuep6",rmestimep6,"rtp6[4][6]/D");
tnd->Branch("l_tdc_valuep6",lmestimep6,"ltp6[4][6]/D");
tnd->Branch("energy_depositp6",ekinaptp6,"ekinp6[4]/D");

tnd->Branch("r_tdc_valuep7",rmestimep7,"rtp7[6][6]/D");
tnd->Branch("l_tdc_valuep7",lmestimep7,"ltp7[6][6]/D");
tnd->Branch("energy_depositp7",ekinaptp7,"ekinp7[6]/D");

tnd->Branch("r_tdc_valuep8",rmestimep8,"rtp8[12][6]/D");
tnd->Branch("l_tdc_valuep8",lmestimep8,"ltp8[12][6]/D");
tnd->Branch("energy_depositp8",ekinaptp8,"ekinp8[12]/D");


istree = 1;

}





// Destructor :
BigBiteEventAction::~BigBiteEventAction(){}

// What we want to do at the begining of each event :
// making the collections of variables


void BigBiteEventAction::BeginOfEventAction(const G4Event*)
{
	// pointer needed to get the collection ID

	sdman = G4SDManager::GetSDMpointer();  //used at the end of the event

	collidaux = -1;
	colliddlt = -1;
	collidmain= -1; 
	collidp0  = -1;
	collidp1  = -1;
	collidp2  = -1;
	collidp3  = -1;
	collidp4  = -1;
	collidp5  = -1;
	collidp6  = -1;
	collidp7  = -1;
	collidp8  = -1; 


	collidaux  = sdman->GetCollectionID("auxbar");
	colliddlt  = sdman->GetCollectionID("delta");
	collidmain = sdman->GetCollectionID("mainE");
	collidp0   = sdman->GetCollectionID("barp0");
	collidp1   = sdman->GetCollectionID("barp1");
	collidp2   = sdman->GetCollectionID("barp2");
	collidp3   = sdman->GetCollectionID("barp3");
	collidp4   = sdman->GetCollectionID("barp4");
	collidp5   = sdman->GetCollectionID("barp5");
	collidp6   = sdman->GetCollectionID("barp6");
	collidp7   = sdman->GetCollectionID("barp7");
	collidp8   = sdman->GetCollectionID("barp8");
}

/************************************************/
// Actions we want to be done at the end of the event :
// reconstruction and storage of variables

void BigBiteEventAction::EndOfEventAction(const G4Event* event)
{

	// Before anything we take the needed pointer on collections...

	G4HCofThisEvent* hcte = event->GetHCofThisEvent();

	if (hcte == NULL) { G4cout << "\n pointer NULL" ; return; }

	ndhitscollection* ndhc   = NULL;  //  needed pointer for  hits collection 

	ndhitscollection* ndhca  = NULL;
	ndhitscollection* ndhcde = NULL;  
	ndhitscollection* ndhce  = NULL;  
	ndhitscollection* ndhcp0 = NULL;  
	ndhitscollection* ndhcp1 = NULL;
	ndhitscollection* ndhcp2 = NULL;
	ndhitscollection* ndhcp3 = NULL;
	ndhitscollection* ndhcp4 = NULL;
	ndhitscollection* ndhcp5 = NULL;
	ndhitscollection* ndhcp6 = NULL;
	ndhitscollection* ndhcp7 = NULL;
	ndhitscollection* ndhcp8 = NULL;

	// For the begining, declaration of useful variables and constants

	G4int nbhits,nbhitsp0,nbhitsp1,nbhitsp2,nbhitsp3,nbhitsp4,nbhitsp5,
		nbhitsp6,nbhitsp7,nbhitsp8;//number of hits
	nbhits =0;
	G4int nbhitsa,nbhitsde,nbhitse;
	G4int bn,bn0,bn1,bn2,bn3,bn4,bn5,bn6,bn7,bn8;
	G4int bna,bnd,bne;//bna added for aux-plane bar number--Ramesh

	G4double ndangle = 25.*deg; // inclination of the exit face
	G4double t2 = tan(ndangle); // useful for reconstruction 
	G4double sca = 25.*mm;   // y-size of auxilliary scintillators bars
	G4double scd = 43.*mm;   //     ""    deltaE and mainE    ""
	// x-position of auxilliary plane
	// G4double dist2target = 10000.; 
	G4double dist2target = 4.*m; //This is the dist from front face 
	//of nd to the beam_line_target--Ramesh added

	// rem : this number has to be put in ndPrimaryGeneratorAction
	// dist2target = ndpga->b2b;
	G4double pldist = 752.9-817.4*.5*tan(ndangle)+270.;
	G4double const xce = 752.9-408.7*t2;  // x position for y=z=0 on exit face
	G4double const xcf = (752.9-408.7*t2)*.5; // x-pos of field center

	G4double ndx = 752.9*mm ; 
	G4double ndy = 817.4*mm ;
	G4double ndx2 = ndx-ndy*tan(ndangle);


	// NOTE : do not forget to change one or more of these values if the
	// caracteristics of BigBite apparatus receives some modifications.

	// G4double  mesrpp,  mestpp,mesrpp0,  thetmp, rtofpp;



	//------------------------------------------------
	/** Step 1 : collecting information */

	// First, the auxilliary plane or wire chambers
	// rem : declarations are made here to avoid compilation problems
	// for auxilliary plane
	G4double rxposa, ryposa, rzposa, zposa, xposa, yposa, auxhtime;

	//  // for wire chambers
	//   G4double rxpos1, rypos1, rzpos1, ypos1, xpos1, zpos1;
	//   G4double rxpos2, rypos2, rzpos2, ypos2, xpos2, zpos2;

	// Energy loss in scintillators
	G4double ekinlossa, ekinlossd, ekinlossm;
	G4double zpde,   zpe,   zpp0,   zpp1,   zpp2,   zpp3,   zpp4;
	zpde=0.;zpe=0.;zpp0=0.;zpp1=0.;zpp2=0.;zpp3=0.;zpp4=0.;





	// for Neutron Detector    

	//if (a=='a') { // if we have the auxilliary plane

	ndhca  = (ndhitscollection*)(hcte->GetHC(collidaux));
	ndhcde = (ndhitscollection*)(hcte->GetHC(colliddlt));
	ndhce  = (ndhitscollection*)(hcte->GetHC(collidmain));


	ndhcp0 = (ndhitscollection*)(hcte->GetHC(collidp0));   
	ndhcp1 = (ndhitscollection*)(hcte->GetHC(collidp1));
	ndhcp2 = (ndhitscollection*)(hcte->GetHC(collidp2));
	ndhcp3 = (ndhitscollection*)(hcte->GetHC(collidp3));
	ndhcp4 = (ndhitscollection*)(hcte->GetHC(collidp4));
	ndhcp5 = (ndhitscollection*)(hcte->GetHC(collidp5));
	ndhcp6 = (ndhitscollection*)(hcte->GetHC(collidp6));
	ndhcp7 = (ndhitscollection*)(hcte->GetHC(collidp7));

	ndhcp8 = (ndhitscollection*)(hcte->GetHC(collidp8));


	nbhitsa  = ndhca->entries();
	nbhitsde = ndhcde->entries();
	nbhitse  = ndhce->entries();


	nbhitsp0 = ndhcp0->entries();
	nbhitsp1 = ndhcp1->entries();
	nbhitsp2 = ndhcp2->entries();
	nbhitsp3 = ndhcp3->entries();
	nbhitsp4 = ndhcp4->entries();
	nbhitsp5 = ndhcp5->entries();
	nbhitsp6 = ndhcp6->entries();
	nbhitsp7 = ndhcp7->entries();
	nbhitsp8 = ndhcp8->entries();




	///////////////////////////////////////////////////


	bn0=0;bn1=0;bn2=0;bn3=0;bn4=0;bn5=0;bn6=0;bn7=0;bn8=0;

	if (nbhitsp0 > 0) { bn0 = (*ndhcp0)[0]->barnumber; }
	if (nbhitsp1 > 0) { bn1 = (*ndhcp1)[0]->barnumber; }
	if (nbhitsp2 > 0) { bn2 = (*ndhcp2)[0]->barnumber; }
	if (nbhitsp3 > 0) { bn3 = (*ndhcp3)[0]->barnumber; }
	if (nbhitsp4 > 0) { bn4 = (*ndhcp4)[0]->barnumber; }
	if (nbhitsp5 > 0) { bn5 = (*ndhcp5)[0]->barnumber; }
	if (nbhitsp6 > 0) { bn6 = (*ndhcp6)[0]->barnumber; }
	if (nbhitsp7 > 0) { bn7 = (*ndhcp7)[0]->barnumber; }
	if (nbhitsp8 > 0) { bn8 = (*ndhcp8)[0]->barnumber; }



	G4double NTheta,NZpos,NYpos,NPhi,rtimep0,ltimep0,rtimep8,ltimep8;
	G4double NYpos0,NYpos1,NYpos2,NYpos3,NYpos4;
	G4double NZpos0,NZpos1,NZpos2,NZpos3,NZpos4;




	const G4int kmax = 1000;

	G4double xhpp0[kmax];
	G4double xhpp1[kmax];
	G4double xhpp2[kmax];
	G4double xhpp3[kmax];
	G4double xhpp4[kmax];
	G4double yhpde[kmax];
	G4double yhpe[kmax];
	G4double yhpp0[kmax];
	G4double yhpp1[kmax];
	G4double yhpp2[kmax];
	G4double yhpp3[kmax];
	G4double yhpp4[kmax];
	G4double zhpde[kmax];
	G4double zhpe[kmax];
	G4double zhpp0[kmax];
	G4double zhpp1[kmax];
	G4double zhpp2[kmax];
	G4double zhpp3[kmax];
	G4double zhpp4[kmax];



	G4double xmhpp0,xmhpp1,xmhpp2,xmhpp3,xmhpp4;//x measured hit position at p0-plane,etc.
	G4double ymhpp0,ymhpp1,ymhpp2,ymhpp3,ymhpp4;
	G4double zmhpp0,zmhpp1,zmhpp2,zmhpp3,zmhpp4;

	G4int nhp0,nhp1,nhp2,nhp3,nhp4;//number of hits at p0, etc.
	G4int hitnum,nhitp;

	G4double scdez = 25. ; //z-position in delta E plane
	G4double scez  = 25. ; //z-position in E plane
	G4double sc0z  = 50. ;
	G4double sc1z  = 50. ;
	G4double sc2z  = 62.5;
	G4double sc3z  = 75.0;
	G4double sc4z  = 125.;


	//G4double scsx = 165.; // distance between the planes//nowhere used

	// G4double basedist = 57.5; // usefull later
	G4double de_dist = 1000.+270.+2.5+945.; //target to 
	//bigbite_origin + origin to aux_face + aux_thickness+ aux to delta E.
	G4double dex     = de_dist;
	G4double ex      = dex + 3.;
	G4double basedist= 3000.;// usefull later.This is the dist 
	//of nd w.r.t. the origin of bigbite-axes.
	G4double p0x = basedist;
	G4double p1x = p0x+50+57.5;
	G4double p2x = p1x+165.;
	G4double p3x = p2x+165.;
	G4double p4x = p3x+165.;

	G4double Tleft,Tright;

	// In this part we add all hits per physical plane of the 
	//detector into one array 

	// add up the hits from different bar types,for each plane 

	nhp0=nbhitsp0;// Number of total hits in p0 
	nhp1=nbhitsp1;
	nhp2=nbhitsp2;
	nhp3=nbhitsp3+nbhitsp4+nbhitsp5+nbhitsp6+nbhitsp7;
	nhp4=nbhitsp8;




	// veto plane

	xmhpp0=0;ymhpp0=0;zmhpp0=0;

	for (hitnum=0;hitnum<nbhitsp0;hitnum++) 
	{
		bn = (*ndhcp0)[hitnum]->barnumber;
		zpp0 +=  (*ndhcp0)[0]->kineloss;
		Tleft = (*ndhcp0)[hitnum]->lmestime;
		Tright= (*ndhcp0)[hitnum]->rmestime;
		NYpos = ((Tright-Tleft)/2 )* 300000*km/s/1.58;
		NZpos = -(16.5-bn)*2*sc0z;
		xhpp0[hitnum]= p0x;
		yhpp0[hitnum]= NYpos;
		zhpp0[hitnum]= NZpos;
		xmhpp0+=p0x;ymhpp0+=NYpos;zmhpp0+=NZpos;

	}

	xmhpp0=xmhpp0/nhp0;ymhpp0=ymhpp0/nhp0;zmhpp0=zmhpp0/nhp0;zpp0=zpp0/nhp0;


	// 10x10 plane 
	xmhpp1=0;ymhpp1=0;zmhpp1=0;
	for (hitnum=0;hitnum<nbhitsp1;hitnum++)
	{
		bn = (*ndhcp1)[hitnum]->barnumber;
		zpp1 +=  (*ndhcp1)[0]->kineloss; 
		Tleft = (*ndhcp1)[hitnum]->lmestime;
		Tright= (*ndhcp1)[hitnum]->rmestime; 
		NYpos = ((Tright-Tleft)/2 )* 300000*km/s/1.58;
		NZpos = -(15.5-bn)*2*sc1z;
		xhpp1[hitnum]= p1x;
		yhpp1[hitnum]= NYpos;
		zhpp1[hitnum]= NZpos;
		xmhpp1+=p1x;ymhpp1+=NYpos;zmhpp1+=NZpos;
	}

	xmhpp1=xmhpp1/nhp1;ymhpp1=ymhpp1/nhp1;zmhpp1=zmhpp1/nhp1;zpp1=zpp1/nhp1; zpp1=zpp1-zpp0;



	// 12.5x10 plane
	xmhpp2=0;ymhpp2=0;zmhpp2=0;
	for (hitnum=0;hitnum<nbhitsp2;hitnum++) { 
		bn = (*ndhcp2)[hitnum]->barnumber;
		zpp2 +=  (*ndhcp2)[0]->kineloss; 
		Tleft = (*ndhcp2)[hitnum]->lmestime;
		Tright= (*ndhcp2)[hitnum]->rmestime;
		NYpos = ((Tright-Tleft)/2 )* 300000*km/s/1.58;
		NZpos = -(12.5-bn)*2*sc2z;
		xhpp2[hitnum]= p2x;
		yhpp2[hitnum]= NYpos;
		zhpp2[hitnum]= NZpos;
		xmhpp2+=p2x;ymhpp2+=NYpos;zmhpp2+=NZpos;
	}

	xmhpp2=xmhpp2/nhp2;ymhpp2=ymhpp2/nhp2;zmhpp2=zmhpp2/nhp2;zpp2=zpp2/nhp2;


	// mixed plane
	xmhpp3=0;ymhpp3=0;zmhpp3=0; 
	for (hitnum=0;hitnum<nbhitsp3;hitnum++) {
		bn = (*ndhcp3)[hitnum]->barnumber; 
		zpp3  +=  (*ndhcp3)[0]->kineloss;
		Tleft = (*ndhcp3)[hitnum]->lmestime;
		Tright= (*ndhcp3)[hitnum]->rmestime; 
		NYpos = ((Tright-Tleft)/2 )* 300000*km/s/1.58;
		NZpos = -sc3z*12-sc2z*8-sc1z*2+(bn-.5)*2*sc3z;
		xhpp3[hitnum]= p3x; 
		yhpp3[hitnum]= NYpos;
		zhpp3[hitnum]= NZpos;
		xmhpp3+=p3x;ymhpp3+=NYpos;zmhpp3+=NZpos;
	}


	for (hitnum=0;hitnum<nbhitsp4;hitnum++) {
		bn = (*ndhcp4)[hitnum]->barnumber;
		zpp3 +=  (*ndhcp4)[0]->kineloss; 
		Tleft = (*ndhcp4)[hitnum]->lmestime; 
		Tright= (*ndhcp4)[hitnum]->rmestime; 
		NYpos = ((Tright-Tleft)/2 )* 300000*km/s/1.58; 
		NZpos = -sc2z*8-sc1z*2+(bn-.5)*2*sc2z; 
		xhpp3[hitnum+nbhitsp3]= p3x; 
		yhpp3[hitnum+nbhitsp3]= NYpos;
		zhpp3[hitnum+nbhitsp3]= NZpos;
		xmhpp3+=p3x;ymhpp3+=NYpos;zmhpp3+=NZpos;
	}

	for (hitnum=0;hitnum<nbhitsp5;hitnum++) {
		bn = (*ndhcp5)[hitnum]->barnumber;
		zpp3 +=  (*ndhcp5)[0]->kineloss;
		Tleft = (*ndhcp5)[hitnum]->lmestime;
		Tright= (*ndhcp5)[hitnum]->rmestime;
		NYpos = ((Tright-Tleft)/2 )* 300000*km/s/1.58;
		NZpos = (bn-1.5)*2*sc1z;
		xhpp3[hitnum+nbhitsp3+nbhitsp4]= p3x;
		yhpp3[hitnum+nbhitsp3+nbhitsp4]= NYpos;
		zhpp3[hitnum+nbhitsp3+nbhitsp4]= NZpos;
		xmhpp3+=p3x;ymhpp3+=NYpos;zmhpp3+=NZpos;
	}

	for (hitnum=0;hitnum<nbhitsp6;hitnum++) { 
		bn = (*ndhcp6)[hitnum]->barnumber;
		zpp3 +=  (*ndhcp6)[0]->kineloss; 
		Tleft = (*ndhcp6)[hitnum]->lmestime; 
		Tright= (*ndhcp6)[hitnum]->rmestime; 
		NYpos = ((Tright-Tleft)/2 )* 300000*km/s/1.58; 
		NZpos =  sc1z*2+(bn-.5)*2*sc2z; 
		xhpp3[hitnum+nbhitsp3+nbhitsp4+nbhitsp5]= p3x; 
		yhpp3[hitnum+nbhitsp3+nbhitsp4+nbhitsp5]= NYpos; 
		zhpp3[hitnum+nbhitsp3+nbhitsp4+nbhitsp5]= NZpos;
		xmhpp3+=p3x;ymhpp3+=NYpos;zmhpp3+=NZpos;
	}

	for (hitnum=0;hitnum<nbhitsp7;hitnum++) {
		bn = (*ndhcp7)[hitnum]->barnumber;
		zpp3 +=  (*ndhcp7)[0]->kineloss;
		Tleft = (*ndhcp7)[hitnum]->lmestime;
		Tright= (*ndhcp7)[hitnum]->rmestime;
		NYpos = ((Tright-Tleft)/2 )* 300000*km/s/1.58;
		NZpos =  sc1z*2+sc2z*8+(bn-.5)*2*sc3z;
		xhpp3[hitnum+nbhitsp3+nbhitsp4+nbhitsp5+nbhitsp6]= p3x;
		yhpp3[hitnum+nbhitsp3+nbhitsp4+nbhitsp5+nbhitsp6]= NYpos;
		zhpp3[hitnum+nbhitsp3+nbhitsp4+nbhitsp5+nbhitsp6]= NZpos;
		xmhpp3+=p3x;ymhpp3+=NYpos;zmhpp3+=NZpos;
	}

	xmhpp3=xmhpp3/nhp3;ymhpp3=ymhpp3/nhp3;zmhpp3=zmhpp3/nhp3;zpp3=zpp3/nhp3;

	// 10x25 plane

	xmhpp4=0;ymhpp4=0;zmhpp4=0;
	for (hitnum=0;hitnum<nbhitsp8;hitnum++) {
		bn = (*ndhcp8)[hitnum]->barnumber;
		zpp4  +=  (*ndhcp8)[0]->kineloss;
		Tleft = (*ndhcp8)[hitnum]->lmestime;
		Tright= (*ndhcp8)[hitnum]->rmestime;
		NYpos = ((Tright-Tleft)/2 )* 300000*km/s/1.58;
		NZpos =  -(6.5-bn)*2*sc4z;
		xhpp4[hitnum]= p4x;
		yhpp4[hitnum]= NYpos;
		zhpp4[hitnum]= NZpos;
		xmhpp4+=p4x;ymhpp4+=NYpos;zmhpp4+=NZpos;
	}

	xmhpp4=xmhpp4/nhp4;ymhpp4=ymhpp4/nhp4;zmhpp4=zmhpp4/nhp4;zpp4=zpp4/nhp4;


	theta =0;phi=0;rthe=0;rphi=0;ptof=0;mestpp=0;
	G4double tpa,tpb,tpd,wpa,wpb,wpd;


	// if we have a hit on all planes 

	if (bn0 != 0 & bn1 !=0 & bn2!=0 & bn8!=0 & (bn3!=0 | bn4!=0 | 
		bn5!=0 |bn6!=0 | bn7!=0 )) 
	{
		NPhi = atan((zmhpp0/(xmhpp0+5000.)+zmhpp2/(xmhpp2+5000.))/2.)*180./3.14;

		rtimep0 = (*ndhcp0)[0]->rmestime;
		ltimep0 = (*ndhcp0)[0]->lmestime;
		rtimep8 = (*ndhcp8)[0]->rmestime;
		ltimep8 = (*ndhcp8)[0]->lmestime;

		NYpos = (rtimep8-ltimep8)/2;
		NTheta = -180.*atan((NYpos * 300000*km/s/1.58)/3000.)/(3.14*3.14);

		G4double tof = (rtimep0 + ltimep0)/2. - 500./300000*s/km*1.58  +
			G4RandGauss::shoot(0.,1.*ns);

		mestpp = dist2target/tof;
		mestpp = mestpp*s/m;
		mestpp = mestpp/sqrt(1-mestpp*mestpp/9.*pow(10.,-16.));  //relativistic 
		mestpp = mestpp*1.6726*3./1.6*pow(10.,-6.);  // in MeV/c

		rthe = NTheta;
		rphi=NPhi;
	}

	/*
	//I want to change there 3
	eventn=ndpga->evnb;
	theta = ndpga->theta;
	phi = ndpga->phi;

	*/ 

	// write tdc and adc to simple ascii file - added by Ran 5/2004

	// Ramesh modified it for outputing variables in the tree--8/2004

	/////////////////////////////////////////////////////////////////////
	// tdc value (time) is in ns and adc value (energy) is in MeV, tested with 
	//the following code--Ramesh.

	// G4double c = 300000.*km/s/1.58; // speed of light in plastic 
	//(here c=189.873)
	// G4double c = 300000000.*m/s/1.58; // speed of light in plastic 
	//(here c=189.87342 again)
	// G4cout <<"\nvalue of c is: "<<c<<G4endl;
	// The way geant4 calculates c is in  mm/ns  if c is given in some units 
	//in the code.

	// G4double energy;
	// energy=5.*GeV;//energy=5000
	// energy=5.*MeV;//energy=5
	// energy=5.*eV;//energy=5e-06
	// G4cout <<"\nMy energy is: "<<energy<<G4endl;
	// The way geant4 calculates energy is in  MeV  if energy is given in 
	//some units in the code.

	// G4double time;
	// time = 3.*s;
	// G4cout <<"\nMy Time is: "<<time<<G4endl;//gives time = 3e+09 => 
	//geant4 measures time in ns.

	// G4double length;
	// length = 4.*km;
	// G4cout<<"\nMy Length is: "<<length<<G4endl;//gives length=4000000=>geant4 
	//measures length in mm.

	// In every plane adc and tdc number goes from 1 to the maximum number of 
	//bar in that plane. This is true
	// for p3 to p7 as well (in these planes this number can not be 
	//greater than 6). Number 1 is at the bottom
	// of every stack (I found this when executing the code only with protons 
	//in which case only one ray (which
	// was neutron produced by collision) was able to pass  the nd somewhere 
	//at the top position, when I saw 
	// the output, it was number 22.). Interesting thing was that 
	//any proton was not able to cross the lead wall.
	////////////////////////////////////////////////////////////////////////



	// // Number of entries in the root:  All entries less than or equal to 10 
	// //in the multihit  entries are registered  and also registered the more 8 
	// //entries for extra eight varibles in the tree.Thus total entries will be 
	// //less or equal to 12 X 10 X number of events + 8.


	/////////////////////////////////////////////////////////////////////////////////


	// Tree filling 
	if (istree != 0) 
	{

		G4int i,bn,barnum,hitnum,hitn,ib,ihit,bnta;
		bnta = 0; 

		// 	 //////////////////////////////////// 
		// 	 eventn=ndpga->evnb;
		// 	 theta = ndpga->theta;
		// 	 phi = ndpga->phi;

		// outfile <<"\n\nTheta : "<<theta <<"\tPhi : "<<phi<<G4endl;

		/////////////////////////////////////////

		//Initializations of all variables so that they can be ready to accept 
		//new value starting from the fresh-start.
		// These two loops at the top of each code-fragment is for 
		//Initialization of variables

		for ( ib=0;ib<56;ib++) {
			for ( ihit=0;ihit<6;ihit++) {
				auxhtimeta[ib][ihit]=0.;      
			}
			ekinapta[ib]=0.; 
		}

		ndhca = (ndhitscollection*)(hcte->GetHC(collidaux)); 
		nbhitsa = ndhca->entries();

		if (nbhitsa !=0)
		{

			G4int numberofhits[56];
			for (bn=0; bn<56;bn++) numberofhits[bn]=0;


			for ( hitn=0;hitn<nbhitsa;hitn++) 
			{

				barnum = (*ndhca)[hitn]->barnumber;
				numberofhits[barnum-1] = numberofhits[barnum-1]+1;
				hitnum = numberofhits[barnum-1];
				if (hitnum>6) break;

				if (floor(bnta/2.)-bnta/2. == 0)
				{ auxhtimeta[barnum-1][hitnum-1] = (*ndhca)[hitn]->rmestime;}
				else{auxhtimeta[barnum-1][hitnum-1] = (*ndhca)[hitn]->lmestime;}

				//Kinetic energy is summed over all multihits 
				ekinapta[barnum-1] = ekinapta[barnum-1]+(*ndhca)[hitn]->kineloss*100000;
				//	 G4cout << "adc  ...." <<  ekinapta[barnum-1]<< endl;

			}        

		}
		//    ////////////////////////////
		for ( ib=0;ib<24;ib++) {
			for (ihit=0;ihit<6;ihit++) {
				rmestimede[ib][ihit]=0.;
				lmestimede[ib][ihit]=0.;
			}
			ekinaptde[ib]=0.;
		}

		ndhcde = (ndhitscollection*)(hcte->GetHC(colliddlt));
		nbhitsde = ndhcde->entries();

		if (nbhitsde !=0)
		{

			G4int numberofhits[24];
			for (bn=0; bn<24;bn++) numberofhits[bn]=0;	     

			for (hitn=0;hitn<nbhitsde;hitn++) 
			{

				barnum = (*ndhcde)[hitn]->barnumber;
				numberofhits[barnum-1] = numberofhits[barnum-1]+1;
				G4int hitnum = numberofhits[barnum-1];
				if (hitnum>6) break;
				rmestimede[barnum-1][hitnum-1] = (*ndhcde)[hitn]->rmestime;
				lmestimede[barnum-1][hitnum-1] = (*ndhcde)[hitn]->lmestime;
				ekinaptde[barnum-1] = ekinaptde[barnum-1] + (*ndhcde)[hitn]->kineloss*100000;


			}

		}
		///////////////////////////////
		for ( ib=0;ib<24;ib++) {
			for ( ihit=0;ihit<6;ihit++) {
				rmestimee[ib][ihit]=0.;
				lmestimee[ib][ihit]=0.;      
			}
			ekinapte[ib]=0.;
		}
		ndhce = (ndhitscollection*)(hcte->GetHC(collidmain));
		nbhitse = ndhce->entries();

		if (nbhitse !=0)
		{

			G4int numberofhits[24];
			for (bn=0; bn<24;bn++) numberofhits[bn]=0;


			for (hitn=0;hitn<nbhitse;hitn++) 
			{

				barnum = (*ndhce)[hitn]->barnumber;
				numberofhits[barnum-1] = numberofhits[barnum-1]+1;
				hitnum = numberofhits[barnum-1];
				if (hitnum>6) break;
				rmestimee[barnum-1][hitnum-1] = (*ndhce)[hitn]->rmestime;
				lmestimee[barnum-1][hitnum-1] = (*ndhce)[hitn]->lmestime;
				ekinapte[barnum-1] =ekinapte[barnum-1]+ (*ndhce)[hitn]->kineloss*100000;

			}

		}


		////////////
		for ( ib=0;ib<30;ib++) {
			for ( ihit=0;ihit<6;ihit++) {
				rmestimep0[ib][ihit]=0.;
				lmestimep0[ib][ihit]=0.;
			}
			ekinaptp0[ib]=0.;
		}

		ndhcp0 = (ndhitscollection*)(hcte->GetHC(collidp0));
		nbhitsp0 = ndhcp0->entries();

		if (nbhitsp0 !=0)
		{

			G4int numberofhits[30];
			for (bn=0; bn<30;bn++) numberofhits[bn]=0;


			for (hitn=0;hitn<nbhitsp0;hitn++) 
			{

				barnum = (*ndhcp0)[hitn]->barnumber;
				numberofhits[barnum-1] = numberofhits[barnum-1]+1;
				hitnum = numberofhits[barnum-1];
				if (hitnum>6) break;
				rmestimep0[barnum-1][hitnum-1] = (*ndhcp0)[hitn]->rmestime;
				lmestimep0[barnum-1][hitnum-1] = (*ndhcp0)[hitn]->lmestime;
				ekinaptp0[barnum-1] =ekinaptp0[barnum-1]+ (*ndhcp0)[hitn]->kineloss*100000;

			}

		}

		//////////
		for ( ib=0;ib<30;ib++) {
			for ( ihit=0;ihit<6;ihit++) { 
				rmestimep1[ib][ihit]=0.;
				lmestimep1[ib][ihit]=0.;
			}
			ekinaptp1[ib]=0.;
		}

		ndhcp1 = (ndhitscollection*)(hcte->GetHC(collidp1));
		nbhitsp1 = ndhcp1->entries();

		if (nbhitsp1 !=0)
		{

			G4int numberofhits[30];
			for (bn=0; bn<30;bn++) numberofhits[bn]=0;


			for (hitn=0;hitn<nbhitsp1;hitn++) 
			{

				barnum = (*ndhcp1)[hitn]->barnumber;

				numberofhits[barnum-1] = numberofhits[barnum-1]+1;
				hitnum = numberofhits[barnum-1];

				// G4cout<<"barnum "<<barnum<<" hitnum "<<hitnum<<G4endl<<G4endl;
				if (hitnum>6) break;
				rmestimep1[barnum-1][hitnum-1] = (*ndhcp1)[hitn]->rmestime;
				lmestimep1[barnum-1][hitnum-1] = (*ndhcp1)[hitn]->lmestime;
				ekinaptp1[barnum-1] = ekinaptp1[barnum-1]+ (*ndhcp1)[hitn]->kineloss*100000;


			}
		}

		////////////
		for ( ib=0;ib<24;ib++) {
			for ( ihit=0;ihit<6;ihit++) {
				rmestimep2[ib][ihit]=0.;
				lmestimep2[ib][ihit]=0.;
			}
			ekinaptp2[ib]=0.;
		}

		ndhcp2 = (ndhitscollection*)(hcte->GetHC(collidp2));
		nbhitsp2 = ndhcp2->entries();

		if (nbhitsp2 !=0)
		{
			G4int numberofhits[24];
			for (bn=0; bn<24;bn++) numberofhits[bn]=0;


			for (hitn=0;hitn<nbhitsp2;hitn++) 
			{

				barnum = (*ndhcp2)[hitn]->barnumber;
				numberofhits[barnum-1] = numberofhits[barnum-1]+1;
				hitnum = numberofhits[barnum-1];
				if (hitnum>6) break;
				rmestimep2[barnum-1][hitnum-1] = (*ndhcp2)[hitn]->rmestime;
				lmestimep2[barnum-1][hitnum-1] = (*ndhcp2)[hitn]->lmestime;
				ekinaptp2[barnum-1] =ekinaptp2[barnum-1]+ (*ndhcp2)[hitn]->kineloss*100000;


			}  

		}

		////////////////
		for ( ib=0;ib<6;ib++) {
			for ( ihit=0;ihit<6;ihit++) {
				rmestimep3[ib][ihit]=0.;
				lmestimep3[ib][ihit]=0.;
			}
			ekinaptp3[ib]=0.;
		}
		ndhcp3 = (ndhitscollection*)(hcte->GetHC(collidp3));
		nbhitsp3 = ndhcp3->entries();

		if (nbhitsp3 !=0)
		{


			G4int numberofhits[6];
			for (bn=0; bn<6;bn++) numberofhits[bn]=0;


			for (hitn=0;hitn<nbhitsp3;hitn++) 
			{

				barnum = (*ndhcp3)[hitn]->barnumber;

				numberofhits[barnum-1] = numberofhits[barnum-1]+1;
				hitnum = numberofhits[barnum-1];
				if (hitnum>6) break;
				rmestimep3[barnum-1][hitnum-1] = (*ndhcp3)[hitn]->rmestime;
				lmestimep3[barnum-1][hitnum-1] = (*ndhcp3)[hitn]->lmestime;
				ekinaptp3[barnum-1] = ekinaptp3[barnum-1]+(*ndhcp3)[hitn]->kineloss*100000;


			}       	                                
		}


		///////////////////////
		for ( ib=0;ib<4;ib++) {
			for ( ihit=0;ihit<6;ihit++) {
				rmestimep4[ib][ihit]=0.;
				lmestimep4[ib][ihit]=0.;      
			}
			ekinaptp4[ib]=0.;
		}

		ndhcp4 = (ndhitscollection*)(hcte->GetHC(collidp4));
		nbhitsp4 = ndhcp4->entries();

		if (nbhitsp4 !=0)
		{

			G4int numberofhits[4];
			for (bn=0; bn<4;bn++) numberofhits[bn]=0;


			for (hitn=0;hitn<nbhitsp4;hitn++) 
			{

				barnum = (*ndhcp4)[hitn]->barnumber;
				numberofhits[barnum-1] = numberofhits[barnum-1]+1;
				hitnum = numberofhits[barnum-1];
				if (hitnum>6) break;
				rmestimep4[barnum-1][hitnum-1] = (*ndhcp4)[hitn]->rmestime;
				lmestimep4[barnum-1][hitnum-1] = (*ndhcp4)[hitn]->lmestime;
				ekinaptp4[barnum-1] = ekinaptp4[barnum-1]+ (*ndhcp4)[hitn]->kineloss*100000;


			}      

		}


		////////////////////
		for ( ib=0;ib<2;ib++) {
			for ( ihit=0;ihit<6;ihit++) { 
				rmestimep5[ib][ihit]=0.;
				lmestimep5[ib][ihit]=0.;
			}
			ekinaptp5[ib]=0.;

		}
		ndhcp5 = (ndhitscollection*)(hcte->GetHC(collidp5));
		nbhitsp5 = ndhcp5->entries();

		if (nbhitsp5 !=0)
		{

			G4int numberofhits[2];
			for (bn=0; bn<2;bn++) numberofhits[bn]=0;


			for (hitn=0;hitn<nbhitsp5;hitn++) 
			{

				barnum = (*ndhcp5)[hitn]->barnumber;
				numberofhits[barnum-1] = numberofhits[barnum-1]+1;
				hitnum = numberofhits[barnum-1];
				if (hitnum>6) break;
				rmestimep5[barnum-1][hitnum-1] = (*ndhcp5)[hitn]->rmestime;
				lmestimep5[barnum-1][hitnum-1] = (*ndhcp5)[hitn]->lmestime;
				ekinaptp5[barnum-1] = ekinaptp5[barnum-1]+ (*ndhcp5)[hitn]->kineloss*100000;
			}

		}  

		///////////////////////////
		for ( ib=0;ib<4;ib++) {
			for ( ihit=0;ihit<6;ihit++) { 
				rmestimep6[ib][ihit]=0.;
				lmestimep6[ib][ihit]=0.;

			}
			ekinaptp6[ib]=0.;
		}
		ndhcp6 = (ndhitscollection*)(hcte->GetHC(collidp6));
		nbhitsp6 = ndhcp6->entries();

		if (nbhitsp6 !=0)
		{

			G4int numberofhits[4];
			for (bn=0; bn<4;bn++) numberofhits[bn]=0;


			for (hitn=0;hitn<nbhitsp6;hitn++) 
			{

				barnum = (*ndhcp6)[hitn]->barnumber;
				numberofhits[barnum-1] = numberofhits[barnum-1]+1;
				hitnum = numberofhits[barnum-1];
				if (hitnum>6) break;
				rmestimep6[barnum-1][hitnum-1] = (*ndhcp6)[hitn]->rmestime;
				lmestimep6[barnum-1][hitnum-1] = (*ndhcp6)[hitn]->lmestime;
				ekinaptp6[barnum-1] = ekinaptp6[barnum-1]+(*ndhcp6)[hitn]->kineloss*100000;



			}     

		}
		///////////////////////
		for ( ib=0;ib<6;ib++) {
			for ( ihit=0;ihit<6;ihit++) { 
				rmestimep7[ib][ihit]=0.;
				lmestimep7[ib][ihit]=0.;

			}
			ekinaptp7[ib]=0.;
		}

		ndhcp7 = (ndhitscollection*)(hcte->GetHC(collidp7));
		nbhitsp7 = ndhcp7->entries();

		if (nbhitsp7 !=0)
		{

			G4int numberofhits[6];
			for (bn=0; bn<6;bn++) numberofhits[bn]=0;


			for (hitn=0;hitn<nbhitsp7;hitn++) 
			{

				barnum = (*ndhcp7)[hitn]->barnumber;
				numberofhits[barnum-1] = numberofhits[barnum-1]+1;
				hitnum = numberofhits[barnum-1];
				if (hitnum>6) break;
				rmestimep7[barnum-1][hitnum-1] = (*ndhcp7)[hitn]->rmestime;
				lmestimep7[barnum-1][hitnum-1] = (*ndhcp7)[hitn]->lmestime;
				ekinaptp7[barnum-1] =  ekinaptp7[barnum-1]+(*ndhcp7)[hitn]->kineloss*100000;


			}  

		}

		///////////////////////////
		for ( ib=0;ib<12;ib++) {
			for ( ihit=0;ihit<6;ihit++) { 
				rmestimep8[ib][ihit]=0.;
				lmestimep8[ib][ihit]=0.;

			}
			ekinaptp8[ib]=0.;
		}
		ndhcp8 = (ndhitscollection*)(hcte->GetHC(collidp8));
		nbhitsp8 = ndhcp8->entries();

		if (nbhitsp8 !=0)
		{

			G4int numberofhits[12];
			for (bn=0; bn<12;bn++) numberofhits[bn]=0;


			for (hitn=0;hitn<nbhitsp8;hitn++) 
			{

				barnum = (*ndhcp8)[hitn]->barnumber;
				numberofhits[barnum-1] = numberofhits[barnum-1]+1;
				hitnum = numberofhits[barnum-1];
				if (hitnum>6) break;      
				rmestimep8[barnum-1][hitnum-1] = (*ndhcp8)[hitn]->rmestime;
				lmestimep8[barnum-1][hitnum-1] = (*ndhcp8)[hitn]->lmestime;
				ekinaptp8[barnum-1] = ekinaptp8[barnum-1]+ (*ndhcp8)[hitn]->kineloss*100000;



			}    


		}




		tnd->Fill();

	}//end of tree feeling "if".

	//} // end of if(a = 'a') at the very begining

}// end of the function EndOfEventAction(const G4Event* event).


// Number of entries in the root:  All entries less than or equal to 10 
//in the multihit  entries are registered  and also registered the more 8 
//entries for extra eight varibles in the tree.Thus total entries will be 
//less or equal to 12 X 10 X number of events + 8.





