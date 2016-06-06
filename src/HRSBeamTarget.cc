#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4Material.hh"
#include "G4RunManager.hh"
#include "G4GeometryManager.hh"
#include "UsageManager.hh"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"

#include "HRSBeamTarget.hh"
#include "HRSMultScatt.hh"

#include "HRSDatabase.hh"

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
using namespace std;

#define __MAX_MAT 100
#define Euler 0.5772157

extern UsageManager* gConfig;
HRSBeamTarget *HRSBeamTarget::gSingleton = NULL;

HRSBeamTarget::HRSBeamTarget(){
    gSingleton = this;
    fMother = NULL;
    UpdateInfo();

    fRasterX = fRasterY = 4.0*mm;
    fX0 = fY0 = fTh0 = fPh0 = fdTh = fdPh = 0.0;

    fCorrTh = fCorrPh = 0.0;

    fMS = new HRSMultScatt();

    //fBeamE   = gDefaultBeamE;
    fBeamPol = gDefaultBeamPol;

    fBeamCurr = gDefaultBeamCur;

    fEcut = 1e-6*MeV;

    fDefaultMat = new G4Material("Default_proton"   , 1., 1.0, 1e-19*g/mole);

    fAlreadyWarned = false;

    G4cout << "Access CREX Database" << G4endl;
    fDatabase = new HRSDatabase(1);

}

HRSBeamTarget::~HRSBeamTarget(){
}

HRSBeamTarget *HRSBeamTarget::GetBeamTarget() {
    if( gSingleton == NULL ){
	gSingleton = new HRSBeamTarget();
    }
    return gSingleton;
}

G4double HRSBeamTarget::GetEffLumin(){
    G4double lumin = fEffMatLen*fBeamCurr/(e_SI*coulomb);
    return lumin;
}

void HRSBeamTarget::UpdateInfo(){
    std::vector<G4VPhysicalVolume *>::iterator it;

    fLH2Length   = -1e9;
    fZpos        = -1e9;
    fLH2pos      = -1e9;
    fTotalLength = 0.0;

    // Can't calculate anything without mother
    if( !fMother ) return;
    fZpos = fMother->GetFrameTranslation().z();

    for(it = fTargVols.begin(); it != fTargVols.end(); it++ ){
	// Assume everything is non-nested tubes
	if( !dynamic_cast<G4Tubs *>( (*it)->GetLogicalVolume()->GetSolid() ) ){
	    G4cerr << "ERROR:  " << __PRETTY_FUNCTION__ << " line " << __LINE__ <<
		":  Target volume not made of G4Tubs" << G4endl; 
	    exit(1);
	}

	if( (*it)->GetLogicalVolume()->GetMaterial()->GetName() == "LiquidHydrogen" ){
	    if( fLH2Length >= 0.0 ){
		G4cerr << "ERROR:  " << __PRETTY_FUNCTION__ << " line " << __LINE__ <<
		    ":  Multiply defined LH2 volumes" << G4endl; 
		exit(1);
	    }

	    fLH2Length = ((G4Tubs *) (*it)->GetLogicalVolume()->GetSolid())->GetZHalfLength()*2.0
		*(*it)->GetLogicalVolume()->GetMaterial()->GetDensity();

	    fLH2pos    = (*it)->GetFrameTranslation().z();

	    fTotalLength += ((G4Tubs *) (*it)->GetLogicalVolume()->GetSolid())->GetZHalfLength()*2.0
		*(*it)->GetLogicalVolume()->GetMaterial()->GetDensity();
	}
    }

    return;
}


void HRSBeamTarget::SetScatAngle(G4double pTh){
  fTh = pTh;
}

void HRSBeamTarget::SetTargetLen(G4double z){
    std::vector<G4VPhysicalVolume *>::iterator it;

    for(it = fTargVols.begin(); it != fTargVols.end(); it++ ){
	G4GeometryManager::GetInstance()->OpenGeometry((*it));
	if( (*it)->GetLogicalVolume()->GetMaterial()->GetName() == "LiquidHydrogen" ){
	    // Change the length of the target volume
	    ((G4Tubs *) (*it)->GetLogicalVolume()->GetSolid())->SetZHalfLength(z/2.0);
	} else {

	    G4cerr << "WARNING " << __PRETTY_FUNCTION__ << " line " << __LINE__ <<
		": volume other than cryogen has been specified, but handling not implemented" << G4endl;
	    // Move position of all other volumes based on half length change

	    /*
	    G4ThreeVector pos = (*it)->GetFrameTranslation();

	    if( pos.z() < fLH2pos ){
		pos = pos + G4ThreeVector(0.0, 0.0, (fLH2Length-z)/2.0 );
	    } else {
		pos = pos - G4ThreeVector(0.0, 0.0, (fLH2Length-z)/2.0 );
	    }

	    (*it)->SetTranslation(pos);
	    */
	}
	G4GeometryManager::GetInstance()->CloseGeometry(true, false, (*it));
    }


    G4RunManager* runManager = G4RunManager::GetRunManager();
    runManager->GeometryHasBeenModified();

    UpdateInfo();
}

void HRSBeamTarget::SetTargetPos(G4double z){
    std::vector<G4VPhysicalVolume *>::iterator it;

    //G4double zshift = z-(fZpos+fLH2pos);


    for(it = fTargVols.begin(); it != fTargVols.end(); it++ ){
	G4GeometryManager::GetInstance()->OpenGeometry((*it));
	if( (*it)->GetLogicalVolume()->GetMaterial()->GetName() == "LiquidHydrogen" ){
	    // Change the length of the target volume
	    (*it)->SetTranslation( G4ThreeVector(0.0, 0.0, z-fZpos) );
	} else {
	    G4cerr << "WARNING " << __PRETTY_FUNCTION__ << " line " << __LINE__ <<
		": volume other than cryogen has been specified, but handling not implemented" << G4endl;

	    // Move position of all other volumes based on half length change

	    /*
	    G4ThreeVector prespos = (*it)->GetFrameTranslation();

	    G4ThreeVector pos = prespos + G4ThreeVector(0.0, 0.0, zshift );

	    (*it)->SetTranslation(prespos);
	    */
	}
	G4GeometryManager::GetInstance()->CloseGeometry(true, false, (*it));
    }

    G4RunManager* runManager = G4RunManager::GetRunManager();
    runManager->GeometryHasBeenModified();

    UpdateInfo();
}


////////////////////////////////////////////////////////////////////////////////////////////
//  Sampling functions

HRSVertex HRSBeamTarget::SampleVertex(G4double ztarget){//ztarget is how far is has gone through the lead
    HRSVertex thisvert;
    
    gConfig->GetArgument("SnakeModel",mSnakeModel);
    gConfig->GetArgument("LHRSMomentum",mLHRSMomentum);
    gConfig->GetArgument("RHRSMomentum",mRHRSMomentum);
    gConfig->GetParameter("LSeptumAngle",mLSeptumAngle);
    gConfig->GetParameter("RSeptumAngle",mRSeptumAngle);

    //G4double ztrav, len;

    // Sample where along target weighted by density (which roughly corresponds to A
    // or the number of electrons, which is probably good enough for this

    // Figure out how far along the target we got
    //switch( samp ){
    //case kCryogen: 
    //fSampLen = fLH2Length;
    //break;

    //case kWalls:
    //G4cerr << "ERROR" << __PRETTY_FUNCTION__ << " line " << __LINE__ <<
    //": scattering from cell walls has been specified, but handling not implemented" << G4endl;
    //exit(1);
    //break;

	    /*
	case kWalls:
	    fSampLen = fTotalLength-fLH2Length;
	    break;
	    */
    //case kFullTarget:
    //fSampLen = fTotalLength;
    //break;
    //}

    //ztrav = CLHEP::RandFlat::shoot(0.0, fSampLen);


    //G4bool isLH2;
    //G4bool foundvol = false;
    //G4Material *mat;
    //G4double zinvol;

    //G4double cumz   = 0.0;

    int radiate = 1;

    G4double radsum = 0.0;

    int      nmsmat;
    double   msthick [__MAX_MAT];
    double   msA     [__MAX_MAT];
    double   msZ     [__MAX_MAT];
    double   mradlen [__MAX_MAT];
    double   mdensity[__MAX_MAT];
    double   mtravel [__MAX_MAT];

    //G4cout << mLHRSMomentum << " <--------> " << gDefaultBeamE_PREX << G4endl;
    if( mSnakeModel != 53){
      fBeamE = mLHRSMomentum;//gDefaultBeamE_PREX;
      //fBeamE = gDefaultBeamE_PREX;
      nmsmat = 2;
      //for PREX. make flag for CREX and other studies
      msthick [0] = 0.15;                   //mm
      msA     [0] = 12.;                    //
      msZ     [0] = 6.;                     //
      mradlen [0] = 42.7 / 10. / 10.;       //g / cm ^ 2 converted to g / mm ^ 2
      mdensity[0] =  3.52 / 10. / 10. / 10.;//g / cm ^ 3 converted to g / mm ^ 3
      mtravel [0] = msthick[0];             //travel through the whole thing
      msthick [0]*= mdensity[0];            //
      
      msthick [1] = 0.5;                    //mm
      msA     [1] = 208;                    //
      msZ     [1] = 82;                     //
      mradlen [1] = 6.37 / 10. / 10.;       //g / cm ^ 2 converted to g / mm ^ 2
      mdensity[1] = 11.34 / 10. / 10. / 10.;//g / cm ^ 3 converted to g / mm ^ 3
      mtravel [1] = ztarget + msthick[1] / 2.;                //travel just through a portion of it. GEANT will take care of the second half
      msthick [1]*= mdensity[1];            //
    }else{
      fBeamE = mLHRSMomentum;//gDefaultBeamE_CREX;
      nmsmat = 1;
      //for CREX, but you have to check these numbers still
      msthick [0] = 0.3;                   //mm
      msA     [0] = 48.;                    //
      msZ     [0] = 20.;                     //
      mradlen [0] = 16.14 / 10. / 10.;       //g / cm ^ 2 converted to g / mm ^ 2
      mdensity[0] =  1.55 / 10. / 10. / 10.;//g / cm ^ 3 converted to g / mm ^ 3
      mtravel [0] = msthick[0];             //travel through the whole thing
      msthick [0]*= mdensity[0];            //      
    }
    //test:
    //mdensity [0] *= 10.;
    //mdensity [1] *= 10.;

    for(Int_t i = 0; i < nmsmat; i++){
      //radsum += mtravel[i] * mdensity[i] / mradlen[i];
      radsum += mtravel[i] * mdensity[i] / mradlen[i];
      //G4cout << mtravel[i] << " " <<  mdensity[i] << " " <<  mradlen[i] << G4endl;
      //G4cout << "My radsum: " << radsum << G4endl;
    }


    //radsum * 4.0; //test
    //But there are also internal effects as well!
    // Approximation for Q2, just needs to be order of magnitude
    //G4double effQ2 = 2.0*fBeamE*fBeamE*(1.0-cos(0.5*deg));
    
    G4double effQ2;
    if(mSnakeModel == 53){
      effQ2 = 2.0*fBeamE*fBeamE*(1.0-cos(4.*deg));//5 deg, or 4, for prex/crex
    }else{
      effQ2 = 2.0*fBeamE*fBeamE*(1.0-cos(5.*deg));//5 deg, or 4, for prex/crex
    }
    // About ~1.5%
    G4double m_e    = 0.511;
    G4double int_bt = 0.75*(alpha/pi)*( log( effQ2/(m_e*m_e) ) - 1.0 );
    //Mo and Tsai eq. I.1, it's an unlabelled eq, on the first page.
    //G4cout << "1: " << fBeamE << G4endl;
    //G4cout << fBeamE << " " << m_e << G4endl;
    //G4cout << int_bt << " " << radsum << G4endl;

    G4double bt  = ( 4.0 / 3.0 ) * ( radsum + int_bt); //contributions first from straggling, but also internal
    //Mo and Tsai eq. A.4, with some approximations.
    //G4cout << radsum << " " << int_bt << G4endl;
    //bt*=4.;

    G4double prob, prob_sample, sample, eloss, value, env;
    G4double ref, fSampE;
    fSampE = fBeamE;
    value = 1.0;
    G4double Ekin = fBeamE - m_e;
    prob = 1.- pow(fEcut/Ekin,bt) - bt/(bt+1.)*(1.- pow(fEcut/Ekin,bt+1.))
      + 0.75*bt/(2.+bt)*(1.- pow(fEcut/Ekin,bt+2.));
    prob = prob/(1.- bt*Euler + bt*bt/2.*(Euler*Euler+pi*pi/6.)); /* Gamma function */
    prob_sample = G4UniformRand();        /* Random sampling */
    
    if (prob_sample <= prob) {//Bremsstrahlung has taken place!                                                                                                                                                       
      do {
        sample = G4UniformRand();
        eloss = fEcut*pow(Ekin/fEcut,sample);
        env = 1./eloss;
        value = 1./eloss*(1.-eloss/Ekin+0.75*pow(eloss/Ekin,2))*pow(eloss/Ekin,bt);
	
        sample = G4UniformRand();
        ref = value/env;
      } while (sample > ref);
      
      //fSampE = Ekin-eloss+m_e;
      //G4cout << fBeamE << " " << eloss << " " << fSampE << " " << m_e << " " << fSampE / fBeamE << G4endl;
      if( radiate )
	fSampE -=eloss;

      assert( fSampE > m_e );
    }    

    //G4cout << "2: " << fSampE << G4endl;

    ///////////////////////////////////////////////////////////////////
    // Sample multiple scattering + angles
    
    G4double msth, msph;
    
    if( nmsmat > 0 ){
      //G4cout << "Initializing multiple scattering\n";
	fMS->Init( fSampE, nmsmat, msthick, msA, msZ );
	msth = fMS->GenerateMSPlane();
	msph = fMS->GenerateMSPlane();
    } else {
	msth = 0.0;
	msph = 0.0;
    }
    //msth = 0.; msph = 0.;//activate to turn off multiple scattering
    //G4cout << "theta and phi rad: " << msth << " " << msph << G4endl;
    //G4cout << "theta and phi deg: " << msth * 180. / pi << " " << msph * 180. / pi << G4endl;

    assert( !std::isnan(msth) && !std::isnan(msph) );


    double ZZZ;
    if( mSnakeModel == 53 ){
      ZZZ = 20.;
    }else{
      ZZZ = 82.;
    }
    //double EEE          = pTrack->P0;
    double EEE          = fSampE;
    //double MMM          = 0.938272 * GeV;
    double MMM;
    if( mSnakeModel == 53 ){
      MMM = 48 * 0.931494095 * GeV;
    }else{
      MMM = 207.9766521 * 0.931494095 * GeV;
    }

    //double theta_over_2 = pTrack->Theta0 / 2.;
    double theta_over_2 = ( msth + fTh ) / 2.;
    //G4cout << "msth and th: " << msth << " " << fTh << G4endl;
    //double QQ2          = 4 * EEE * EEE * sin( theta_over_2 );
    double QQ2          = 4 * EEE * EEE * pow( sin( theta_over_2 ), 2) / ( 1 - 2 * EEE / MMM * pow( sin( theta_over_2 ), 2) );
    double ascreen      = 5e5;
    double xscreen      = 38. / ( pow( ascreen , 2 ) );
    double mott         = pow( ZZZ * .197 * cos( theta_over_2 ) / 137. , 2 ) /
      ( xscreen + .0004 * EEE * EEE * pow( sin( theta_over_2 ) , 4 ) );//mott is in barns, and all energy is in MeV
    //  ( xscreen + 400. * EEE * EEE * pow( sin( theta_over_2 ) , 4 ) );
    double FFF;
    if( mSnakeModel == 53){
      //FFF = InterpolateNickieCa48(QQ2 / GeV / GeV);//this is indeed FF squared, already included.
      //G4cout << "About to interpolate" << G4endl;
      FFF = fDatabase->Interpolate(EEE / GeV, fTh * 180. / pi, 0, 0);//0 means non stretched, and XS, whereas 1, 1 means stretched and asym.
      //G4cout << "Interpolated" << G4endl;
    }else{
      FFF = InterpolateNickie(QQ2 / GeV / GeV);//this is indeed FF squared, already included.
    }
    //G4cout << ZZZ << " " << EEE << " " << MMM << " " << theta_over_2 * 180. / 3.141592654<< " " << QQ2 << " " << mott << " " << FFF << G4endl;

    double XXS    = mott * FFF;// * 1000000000.; //convert to nb from barns
    //cout << "Cross section for 208Pb at Q2 = " << QQ2 << " is " << XXS_208Pb << " which is at polar angle " << theta_over_2 * 2. * 180. / pi<< endl;
    //pTrack->rate_208Pb = XXS_208Pb;

    //my @t = (4.00, 4.00, 6.00, 4.50, 5.50, 6.50);
    //my @t1= (3.00, 3.00, 4.00, 3.00, 4.00, 4.00);
    //my @t2= (6.00, 6.00, 8.00, 6.50, 7.50, 8.50);
    
    double minimum_angle;
    //G4cout << mLSeptumAngle << " is the septum angle right now." << G4endl;
    if( mLSeptumAngle == 4. ){
      minimum_angle = 3.;
      //G4cout << minimum_angle << " is the minimum angle right now." << G4endl;
    }else if ( mLSeptumAngle == 6.  ){
      minimum_angle = 4.;
      //G4cout << minimum_angle << " is the minimum angle right now." << G4endl;
    }else if ( mLSeptumAngle == 4.5 ){
      minimum_angle = 3.;
      //G4cout << minimum_angle << " is the minimum angle right now." << G4endl;
    }else if ( mLSeptumAngle == 5.5 ){
      minimum_angle = 4.;
      //G4cout << minimum_angle << " is the minimum angle right now." << G4endl;
    }else if ( mLSeptumAngle == 6.5 ){
      minimum_angle = 4.;
      //G4cout << minimum_angle << " is the minimum angle right now." << G4endl;
    }
    double anglemin_over_2 = minimum_angle / 2 * pi / 180.;//half of min angle
    double QQ2min       = 4 * EEE * EEE * pow( sin( anglemin_over_2 ), 2) / ( 1 - 2 * EEE / MMM * pow( sin( anglemin_over_2 ), 2) );
    double FFFmax       = InterpolateNickie(QQ2min / GeV / GeV);//this is indeed FF squared, already included.
    //    double mottmax = pow( ZZZ * .197 * cos( anglemin_over_2 ) / 137. , 2 ) /
    //( xscreen + .0004 * EEE * EEE * pow( sin( anglemin_over_2 ) , 4 ) );
    double mottmax = pow( ZZZ * .197 * cos( anglemin_over_2 ) / 137. , 2 ) /
      ( xscreen + .0004 * EEE * EEE * pow( sin( anglemin_over_2 ) , 4 ) );
    
    prob_sample = G4UniformRand() * mottmax * FFFmax;
    thisvert.pass = 0;
    //G4cout << thisvert.pass << " " << prob_sample << " " << mott * FFF << ", mott = " << mott << ", FFF = " << FFF << G4endl ;
    
    if (prob_sample <= mott * FFF) {
      thisvert.pass = 1;      
      
      // REradiate////////////////////////////////////////////////////////////////////////////
      // We're going to use the new kinematics for this guy
      
      G4double new_mass_c2;
      if( mSnakeModel == 53 ){
	new_mass_c2 = 48 * 0.931494095 * GeV;
      }else{
	new_mass_c2 = 207.9766521 * 0.931494095 * GeV;
      }
      //G4cout << "Pb mass in MeV? " << new_mass_c2 << G4endl;
      
      double ef  = new_mass_c2*fSampE/(new_mass_c2 + fSampE*(1.0-cos(fTh + msth)));//this should be the actual angle of scattering
      double q2  = 2.0*fSampE*ef*(1.0-cos(fTh + msth));//should be actual angle scattering
      
      //G4cout << ef << " " << q2 << G4endl;
      
      int_bt = (alpha/pi)*( log( q2/(m_e*m_e) ) - 1.0 );
      Ekin = ef - electron_mass_c2;;
      
      //int_bt *= 4.;//test
      //G4cout << int_bt << G4endl;
      
      prob = 1.- pow(fEcut/Ekin, int_bt) - int_bt/(int_bt+1.)*(1.- pow(fEcut/Ekin,int_bt+1.))
	+ 0.75*int_bt/(2.+int_bt)*(1.- pow(fEcut/Ekin,int_bt+2.));
      prob = prob/(1.- int_bt*Euler + int_bt*int_bt/2.*(Euler*Euler+pi*pi/6.)); /* Gamma function */
      prob_sample = G4UniformRand();        /* Random sampling */
      
      if (prob_sample <= prob) {//Bremsstrahlung has taken place!
	do {
	  sample = G4UniformRand();
	  eloss = fEcut*pow(Ekin/fEcut,sample);
	  env = 1./eloss;
	  value = 1./eloss*(1.-eloss/Ekin+0.75*pow(eloss/Ekin,2))*pow(eloss/Ekin,bt);
	  
	  sample = G4UniformRand();
	  ref = value/env;
	} while (sample > ref);
	
	//ef = Ekin-eloss+electron_mass_c2;
	//G4cout << ef << " " << eloss << G4endl;
	if( radiate )
	  ef -= eloss;

	assert( ef > electron_mass_c2 );
      }
      
      //G4cout << "3: " << ef << G4endl;
      
      // Sample beam energy based on radiation
      // We do this so it doesn't affect the weighting
      //
      // This can be ignored and done in a generator by itself
      /*
	fRadLen = radsum;
	
	G4double  Ekin = fBeamE - electron_mass_c2;
	G4double  bt   = fRadLen*4.0/3.0;
	G4double  prob_sample, eloss, sample, env, value, ref;
	
	//G4cout << Ekin << " " << bt << " " << fRadLen << G4endl;
	//G4cout << "Ekin  bt fRadLen" << G4endl;
	
	G4double prob = 1.- pow(fEcut/Ekin,bt) - bt/(bt+1.)*(1.- pow(fEcut/Ekin,bt+1.))
	+ 0.75*bt/(2.+bt)*(1.- pow(fEcut/Ekin,bt+2.));
	prob = prob/(1.- bt*Euler + bt*bt/2.*(Euler*Euler+pi*pi/6.)); */ // Gamma function 
      
      //G4cout << "My prob: " << prob << G4endl;
      /*
	prob_sample = G4UniformRand();
	
	//G4cout << "My prob sample: " << prob_sample << G4endl;
	
	if (prob_sample <= prob) {
	do {
	sample = G4UniformRand();
	eloss = fEcut*pow(Ekin/fEcut,sample);
	env = 1./eloss;
	value = 1./eloss*(1.-eloss/Ekin+0.75*pow(eloss/Ekin,2))*pow(eloss/Ekin,bt);
	
	sample = G4UniformRand();
	ref = value/env;
	} while (sample > ref);
	
	//G4cout << "eloss, fBeamE" << G4endl;
	//G4cout << eloss << " " << fBeamE << G4endl;
	
	fSampE = fBeamE - eloss;
	assert( fSampE > electron_mass_c2 );
	} else {
	fSampE = fBeamE;
	}
      */
      
      
      //thisvert.fBeamE = fSampE;
      //thisvert.fBeamE = fBeamE;
      thisvert.fBeamE = ef;
      thisvert.fmsth = msth;
      thisvert.fmsph = msph;
      thisvert.XS_208Pb = XXS;
      
      assert( thisvert.fBeamE >= electron_mass_c2 );
    }
    //G4cout << "exiting beamtarget" << G4endl;
    return thisvert;
}

double HRSBeamTarget::InterpolateNickie(double Q2_in){
  ifstream INFILE;
  INFILE.open("/home/Nickie/JLab/HallA/G4MC/CREX_Project/mefcal.pb208_1_25deg_fine.out");
  
  string ignore0, ignore1, ignore2, ignore3, ignore4, ignore5, ignore6, ignore7;
  double lastQ2     = 0.;
  double currentQ2  = 0.;
  double lastFF2    = 0.;
  double currentFF2 = 0.;
  string q_data     = "";
  string ff2_data   = "";
  while( (! INFILE.eof()) && Q2_in > currentQ2 ){
    //cout << Q2_in << " versus " << currentQ2 << endl;
    lastQ2  = currentQ2;
    lastFF2 = currentFF2;
    INFILE >> ignore0 >> q_data >> ignore1 >> ignore2 >> ff2_data >> ignore3 >> ignore4 >> ignore5 >> ignore6 >> ignore7;
    currentQ2 = atof( q_data.c_str() ) * .197; //convert table value of inverse fm to GeV.
    currentQ2*= currentQ2;
    currentFF2= atof( ff2_data.c_str() );
    //table is actually q, not q2
    //G4cout << lastQ2 << " " << currentQ2 << " " << lastFF2 << " " << currentFF2 << G4endl;
    //G4cout   << " " << ignore1 << " " << currentQ2 << " " << ignore2 << " " << ignore3 << " " << currentFF2
    //<< " " << ignore4 << " " << ignore5   << " " << ignore6 << " " << ignore7 << G4endl;
    
  }
  INFILE.close();
  //float mmm = ( atof(currentFF2.c_str()) - atof(lastFF2.c_str()) ) / ( atof(currentQ2.c_str()) - atof(lastQ2.c_str()) );
  float mmm = ( currentFF2 - lastFF2 ) / ( currentQ2 - lastQ2 );
  //cout << lastQ2 << " " << currentQ2 << " " << lastFF2 << " " << currentFF2 << endl;
  //float FFF = mmm * ( Q2_in - atof(currentQ2.c_str()) ) + atof(currentFF2.c_str());
  float FFF = mmm * ( Q2_in - currentQ2 ) + currentFF2;
  //G4cout << mmm << " " << FFF << G4endl;
  return FFF;
}

double HRSBeamTarget::InterpolateNickieCa48(double Q2_in){
  ifstream INFILE;
  INFILE.open("/home/Nickie/JLab/HallA/G4MC/CREX_Project/calcium48.txt");
  
  double lastQ2     = 0.;
  double currentQ2  = 0.;
  double lastFF2    = 0.;
  double currentFF2 = 0.;
  string q_data     = "";
  string ff2_data   = "";
  //THE TABLE IS IN Q, FF2, and by the way, Q is in inverse fermi.

  while( (! INFILE.eof()) && Q2_in > currentQ2 ){
    //cout << Q2_in << " versus " << currentQ2 << endl;
    lastQ2  = currentQ2;
    lastFF2 = currentFF2;
    INFILE >> q_data >> ff2_data;
    currentQ2 = atof( q_data.c_str() ) * .197; //convert table value of inverse fm to GeV.
    currentQ2*= currentQ2;//you don't want q, as is given, you want Q^2
    currentFF2= atof( ff2_data.c_str() );
    //table is actually q, not q2
    //G4cout << lastQ2 << " " << currentQ2 << " " << lastFF2 << " " << currentFF2 << G4endl;
    //G4cout   << " " << ignore1 << " " << currentQ2 << " " << ignore2 << " " << ignore3 << " " << currentFF2
    //<< " " << ignore4 << " " << ignore5   << " " << ignore6 << " " << ignore7 << G4endl;
    
  }
  INFILE.close();
  //float mmm = ( atof(currentFF2.c_str()) - atof(lastFF2.c_str()) ) / ( atof(currentQ2.c_str()) - atof(lastQ2.c_str()) );
  float mmm = ( currentFF2 - lastFF2 ) / ( currentQ2 - lastQ2 );
  //cout << lastQ2 << " " << currentQ2 << " " << lastFF2 << " " << currentFF2 << endl;
  //float FFF = mmm * ( Q2_in - atof(currentQ2.c_str()) ) + atof(currentFF2.c_str());
  float FFF = mmm * ( Q2_in - currentQ2 ) + currentFF2;
  //G4cout << mmm << " " << FFF << G4endl;
  return FFF;
}
