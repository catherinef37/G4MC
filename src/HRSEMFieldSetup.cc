// ********************************************************************
//
// $Id: HRSEMField.hh,v 1.0, 2010/12/26 HRS Exp $
// GEANT4 tag $Name: geant4-09-04 $
//
//
//   User Field class Setup implementation.
//
#include "HRSEMFieldSetup.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4EquationOfMotion.hh"
#include "G4EqMagElectricField.hh"
#include "BField_Dipole.hh"
#include "BField_Quad.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ChordFinder.hh"
#include "UsageManager.hh"

#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"
#include "G4ios.hh"
//#include "QuadFringe.hh"

#include <iostream>
using namespace std;

extern UsageManager* gConfig;	

//////////////////////////////////////////////////////////////////////////
HRSEMFieldSetup* HRSEMFieldSetup::fHRSEMFieldSetup=0;
HRSEMFieldSetup* HRSEMFieldSetup::GetHRSEMFieldSetup()
{ 
	if(!fHRSEMFieldSetup)  
	{
		G4cout<<"HRSEMFieldSetup is not initialized yet...exit...\n";
		exit(-99);
	}
	return fHRSEMFieldSetup; 
}

//////////////////////////////////////////////////////////////////////////
//
HRSEMFieldSetup::HRSEMFieldSetup()
: fChordFinder(0), fStepper(0), fIntgrDriver(0)
{
  gConfig->GetArgument("LHRSMomentum",mLHRSMomentum);
  gConfig->GetArgument("RHRSMomentum",mRHRSMomentum);
  mLHRSMomentum /= 1000.;
  mRHRSMomentum /= 1000.;
  gConfig->GetArgument("SnakeModel",mSnakeModel);
  gConfig->GetParameter("LHRSAngle",mLHRSAngle);
  mLHRSAngle*=deg;
  gConfig->GetParameter("RHRSAngle",mRHRSAngle);
  mRHRSAngle*=deg;
  gConfig->GetParameter("LSeptumAngle",mLSeptumAngle);
  gConfig->GetParameter("RSeptumAngle",mRSeptumAngle);

  G4cout << "HRS angles: " << mLHRSAngle << " " << mRHRSAngle << G4endl;

  //G4cout << "Quad fringe?" << G4endl;
  //QuadFringe* fringe = new QuadFringe();
  //G4cout << "Quad fringe!" << G4endl;
  fHRSEMFieldSetup=this;
  
  //global EM field
  fEMfield = new HRSEMField();
  messenger = new HRSEMFieldSetupMessenger(this) ;
  fEquation = new G4EqMagElectricField(fEMfield);
  fMinStep  = 0.00001*mm ; // minimal step of 1 miron, default is 0.01 mm, Nickie finely
  fStepperType = 4 ;     // ClassicalRK4 -- the default stepper
  
  fFieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  
  fChordFinder = 0;   //will be set in UpdateField()
  UpdateField();
  
  //G4double snakemagnumber = 1. / -4.77577734;//unitless, STD - correct, I believe
  //G4double snakemagnumber = 1. / -4.77577734 * ( 4.00 / 1.063 ) ;//unitless, PREX tune B - correct, I believe
  G4double snakemagnumber = -4.77577734 / 1.063;//Nickie's calculation//this is just 4.77577734 / 1.063
  //snakemagnumber *= 0.83756;
  snakemagnumber *= 0.83762;
  //snakemagnumber *= 0.7;//This is tuned with hrstrans #1
  //snakemagnumber *= 1.33;//This is tuned with hrstrans when hrstrans is tuned to JLR
  //snakemagnumber *= 1.24;//This is my new tune, with snake d.dat problems fixed #2
  //snakemagnumber *= 0.965; #3
  
  if( mSnakeModel == 53 || mSnakeModel == 55){
    //snakemagnumber *= 1.063 / 2.2;
    snakemagnumber *= 1.063 / mLHRSMomentum;    
  }
  
  G4int    quads_on = 1;
  G4double KAPPA1 =  0.;
  G4double KAPPA2 =  0.;
  G4double KAPPA3 =  0.;
  G4int sos = 1;
  if( mSnakeModel > 52 ){
    sos = 1;
  }

  double a_sup = 0.1492;                    //Bore radius, in meters
  double l_sup = 0.9413;                    //Length of quad, in meters
  double a_sos = 0.12827;
  double l_sos = 0.70;
  //BELOW ARE SNAKE VALUES, NOT NIM VALUES                                                                                  
  double pTarget =             0.0  * cm;
  double pQ1en   = pTarget + 160.0  * cm;
  double pQ1ex   = pQ1en   +  94.13 * cm;
  double pQ2en   = pQ1ex   + 115.58 * cm;
  double pQ2ex   = pQ2en   + 182.66 * cm;
  //ABOVE ARE SNAKE VALUES, NOT NIM VALUES  
  //Local field  FZB2, Q1
  //double pHallCenter2Q1Face=1.69*m;//NIM
  double pHallCenter2Q1Face=pQ1en;//SNAKE
  //double pQ1Length= 80*cm;//NIM
  double pQ1Length= sos ? 70.    * cm : 94.13*cm;//SNAKE
  double pQ1Radius= sos ? 12.827 * cm : 0.1492 * m;//SNAKE
  double q1shift = sos ? 0.0 * m : 0.0 * m;
  double pQ1Pos_Z=(pHallCenter2Q1Face+94.13*cm/2.0+q1shift);//NIM
  //double pHallCenter2Q2Face=3.74*m;//NIM
  double pHallCenter2Q2Face=pQ2en;//SNAKE
  //double pQ2Length=180*cm;//NIM
  double pQ2Length=182.66*cm;//SNAKE
  double pQ2Radius=  0.3 * m;//SNAKE
  double pQ2Pos_Z=(pHallCenter2Q2Face+pQ2Length/2.0);
  //double pQ3Length=180*cm;//NIM
  double pQ3Length=182.68*cm;//SNAKE
  double pQ3Radius=  pQ2Radius;//SNAKE

  if( quads_on == 1 && ( mSnakeModel == 49 || mSnakeModel > 50 ) ){
    //sos: Pole tip field when particle momentum is 837.56 MeV/c = 0.2804 Tesla
    //G4cout << "The momentum:! " << mRHRSMomentum << G4endl;
    //mRHRSMomentum is in fact in GeV, and not in MeV.
    //need an extra minus for the sos case, for the direction of focusing...
    //KAPPA1 = sos ? ( 0.2804 * tesla) * mRHRSMomentum / 0.83756 : -0.8476  * tesla / snakemagnumber;//diameter is taken into account later, in BField_Quad//This is the current one I use, usually.
    //KAPPA1 = sos ? 1.5554 * -0.8476  * tesla / snakemagnumber : -0.8476  * tesla / snakemagnumber; //test
    //KAPPA1 = sos ? 0.302875 * tesla : -0.8476  * tesla / snakemagnumber; //test2 , this one is great
    if( ( mSnakeModel == 55 && mLSeptumAngle >= 4.9 ) || 
	( mSnakeModel != 53 && mSnakeModel   != 54 && mSnakeModel   != 55 ) ){//prex, 5 degrees
      G4cout << "IN 5 DEGREE MODE!" << G4endl;
      G4cout << mSnakeModel << " " << mLSeptumAngle << " " << mRSeptumAngle << G4endl;
      KAPPA1 = sos ? 0.260387 * tesla : -0.8476  * tesla / snakemagnumber; //test3
      //KAPPA1 = -0.8476  * tesla / snakemagnumber;
      //G4cout << "Ratio " << mRHRSMomentum << " " << 0.83756 << G4endl;
      //KAPPA1 *= 1.1;
      //KAPPA1 *= 0.75;
      KAPPA2 =  0.93528 * tesla / snakemagnumber;//so don't worry!
      KAPPA3 =  1.15762 * tesla / snakemagnumber;
    }else{//CREX, 4 degrees
      G4cout << "IN 4 DEGREE MODE!" << G4endl;
      KAPPA1 = -0.7146 * tesla / snakemagnumber * l_sup / l_sos * a_sos / a_sup;
      KAPPA1 *= 0.85;//this is acceptable as a first edit of Seamus
      //KAPPA1 *= 0.8;//Dustin wants me to tighten the transverse part
      KAPPA2 =  0.8680 * tesla / snakemagnumber;
      KAPPA3 =  1.1748 * tesla / snakemagnumber;
    }
    if( mSnakeModel == 55 && mLSeptumAngle == 6. ){
      G4cout << "Checking the new tune at 6 degrees" << G4endl;
      G4cout << KAPPA1 << " " << KAPPA2 << " " << KAPPA3 << G4endl;
      KAPPA1 = 1.16070 * tesla * pQ1Radius / 1000.;
      KAPPA2 =-1.37117 * tesla * pQ2Radius / 1000.;
      KAPPA3 =-1.72440 * tesla * pQ3Radius / 1000.;
      G4cout << pQ1Radius << " " << pQ2Radius << " " << pQ3Radius << " " << tesla << G4endl;
      G4cout << KAPPA1 << " " << KAPPA2 << " " << KAPPA3 << G4endl;
    }else if(mSnakeModel == 55 && mLSeptumAngle == 4.5){
      G4cout << "Checking the new tune at 4.5 degrees" << G4endl;
      G4cout << KAPPA1 << " " << KAPPA2 << " " << KAPPA3 << G4endl;
      KAPPA1 = 1.02555 * tesla * pQ1Radius / 1000. / 2.2 * mLHRSMomentum;;
      KAPPA2 =-1.33393 * tesla * pQ2Radius / 1000. / 2.2 * mLHRSMomentum;;
      KAPPA3 =-1.70387 * tesla * pQ3Radius / 1000. / 2.2 * mLHRSMomentum;;
      G4cout << pQ1Radius << " " << pQ2Radius << " " << pQ3Radius << " " << tesla << G4endl;
      G4cout << KAPPA1 << " " << KAPPA2 << " " << KAPPA3 << G4endl;
    }else if(mSnakeModel == 55 && mLSeptumAngle == 5.5){
      G4cout << "Checking the new tune at 5.5 degrees" << G4endl;
      G4cout << KAPPA1 << " " << KAPPA2 << " " << KAPPA3 << G4endl;
      KAPPA1 = 1.10973 * tesla * pQ1Radius / 1000. / 2.2 * mLHRSMomentum;;
      KAPPA2 =-1.35688 * tesla * pQ2Radius / 1000. / 2.2 * mLHRSMomentum;;
      KAPPA3 =-1.71687 * tesla * pQ3Radius / 1000. / 2.2 * mLHRSMomentum;;
      G4cout << pQ1Radius << " " << pQ2Radius << " " << pQ3Radius << " " << tesla << G4endl;
      G4cout << KAPPA1 << " " << KAPPA2 << " " << KAPPA3 << G4endl;
    }else if(mSnakeModel == 55 && mLSeptumAngle == 6.5){
      G4cout << "Checking the new tune at 6.5 degrees" << G4endl;
      G4cout << KAPPA1 << " " << KAPPA2 << " " << KAPPA3 << G4endl;
      KAPPA1 = 1.20291 * tesla * pQ1Radius / 1000. / 2.2 * mLHRSMomentum;;
      KAPPA2 =-1.38139 * tesla * pQ2Radius / 1000. / 2.2 * mLHRSMomentum;;
      KAPPA3 =-1.73056 * tesla * pQ3Radius / 1000. / 2.2 * mLHRSMomentum;;
      G4cout << pQ1Radius << " " << pQ2Radius << " " << pQ3Radius << " " << tesla << G4endl;
      G4cout << KAPPA1 << " " << KAPPA2 << " " << KAPPA3 << G4endl;
    }
    //this is just for a test. You better comment this line out.
    //KAPPA3 *= 1.1;

    //G4cout << "gradients: " << KAPPA1 << " " << KAPPA2 << " " << KAPPA3 << G4endl;
    //G4cout << (- 0.2804 * tesla) * mRHRSMomentum / 837.56 << " " <<  -0.8476  * tesla / snakemagnumber << G4endl;
    //G4cout << mRHRSMomentum << " " << snakemagnumber << G4endl;
    //KAPPA1 = -0.8476  * tesla / snakemagnumber / .1492 / m; //PREX tune
    //KAPPA2 =  0.93528 * tesla / snakemagnumber / .300  / m; //PREX tune
    //KAPPA3 =  1.15762 * tesla / snakemagnumber / .300  / m; //PREX tune
    
    //KAPPA1 = -0.2445  * tesla / snakemagnumber / .1492 / m; //STD tune
    //KAPPA2 =  0.1939  * tesla / snakemagnumber / .300  / m; //STD tune
    //KAPPA3 =  0.1794  * tesla / snakemagnumber / .300  / m; //STD tune

  }
  //G4double dipoleField = -0.4192 * tesla; //using dipolert2 from SNAKE, thin vb, this one is reconciled with SNAKE
  G4double dipoleField = -0.4205 * tesla; //matches nicely with DATA (not snake comparison) 
  if( mSnakeModel == 53 || mSnakeModel == 55 ){
    //dipoleField *= 2.2 / 1.063;
    dipoleField *= mLHRSMomentum / 1.063;
  }
  G4RotationMatrix* LROTATED = new G4RotationMatrix;
  G4RotationMatrix* RROTATED = new G4RotationMatrix;
  LROTATED->rotateY( mLHRSAngle );
  RROTATED->rotateY( mRHRSAngle );
  double pDipoleRCenterY=8.4  * m;
  double pDipoleRCenterZ=9.961 * m;//SNAKE

  G4ThreeVector     LORIGIND(pDipoleRCenterZ * sin( mLHRSAngle ), pDipoleRCenterY, pDipoleRCenterZ * cos( mLHRSAngle ));
  G4ThreeVector     RORIGIND(pDipoleRCenterZ * sin( mRHRSAngle ), pDipoleRCenterY, pDipoleRCenterZ * cos( mRHRSAngle ));
    
  
  G4ThreeVector     LORIGINQ1(pQ1Pos_Z * sin(mLHRSAngle), 0., pQ1Pos_Z * cos(mLHRSAngle));
  G4ThreeVector     RORIGINQ1(pQ1Pos_Z * sin(mRHRSAngle), 0., pQ1Pos_Z * cos(mRHRSAngle));
  G4RotationMatrix* LROTATEQ1 = new G4RotationMatrix;
  G4RotationMatrix* RROTATEQ1 = new G4RotationMatrix;
  LROTATEQ1->rotateY( mLHRSAngle);
  RROTATEQ1->rotateY( mRHRSAngle);
  
  G4ThreeVector     LORIGINQ2(pQ2Pos_Z * sin(mLHRSAngle), 0., pQ2Pos_Z * cos(mLHRSAngle));
  G4ThreeVector     RORIGINQ2(pQ2Pos_Z * sin(mRHRSAngle), 0., pQ2Pos_Z * cos(mRHRSAngle));
  G4RotationMatrix* LROTATEQ2 = new G4RotationMatrix;
  G4RotationMatrix* RROTATEQ2 = new G4RotationMatrix;
  
  LROTATEQ2->rotateY( mLHRSAngle);
  RROTATEQ2->rotateY( mRHRSAngle);
  
  //double pFPR       = 8.4  * m;//radius of curvature of dipole
  //double pFPA       = 9.96 * m;//distance from pivot to entrance of dipole//NIM
  //double pFPA       = pDipoleRCenterZ ;//distance from pivot to entrance of dipole//SNAKE
  //double pFPH       = pFPR * tan ( 22.5 * deg );//height
  //double pFPCenterX = ( pFPA + pFPH + ( pFPH + 1.5 * m + 0.9 * m ) / sqrt(2) ) *-sin( mRHRSAngle ); //not including half of Q3//NIM
  //double pFPCenterZ = ( pFPA + pFPH + ( pFPH + 1.5 * m + 0.9 * m ) / sqrt(2) ) * cos( mRHRSAngle ); //not including half of Q3//NIM
  //double pFPCenterY = ( pFPH + 1.5 * m + 0.9 * m ) / sqrt(2); //including half of Q3//NIM
  
  double pLQ3CenterX = (17.0257042 * m  + pQ3Length / sqrt(2.) / 2. ) *  sin( mLHRSAngle );//SNAKE
  double pLQ3CenterZ = (17.0257042 * m  + pQ3Length / sqrt(2.) / 2. ) *  cos( mLHRSAngle );//SNAKE
  double pLQ3CenterY = ( 3.5853101 * m  + pQ3Length / sqrt(2.) / 2. );//SNAKE	
  double pRQ3CenterX = (17.0257042 * m  + pQ3Length / sqrt(2.) / 2. ) *  sin( mRHRSAngle );//SNAKE
  double pRQ3CenterZ = (17.0257042 * m  + pQ3Length / sqrt(2.) / 2. ) *  cos( mRHRSAngle );//SNAKE
  double pRQ3CenterY = ( 3.5853101 * m  + pQ3Length / sqrt(2.) / 2. );//SNAKE	
  
  //cout << pQ3CenterX << " " << pFPCenterX << endl;
  //cout << pQ3CenterY << " " << pFPCenterY << endl;
  //cout << pQ3CenterZ << " " << pFPCenterZ << endl;
  G4ThreeVector     LORIGINQ3(pLQ3CenterX, pLQ3CenterY, pLQ3CenterZ);//try this one
  G4ThreeVector     RORIGINQ3(pRQ3CenterX, pRQ3CenterY, pRQ3CenterZ);//try this one
  G4RotationMatrix* LROTATEQ3 = new G4RotationMatrix;
  G4RotationMatrix* RROTATEQ3 = new G4RotationMatrix;
  
  LROTATEQ3->rotateX(-45.0 * deg);
  LROTATEQ3->rotateY( mLHRSAngle);
  RROTATEQ3->rotateX(-45.0 * deg);
  RROTATEQ3->rotateY( mRHRSAngle);
  
  fMagFieldFZBL1 = new BField_Quad(KAPPA1, LORIGINQ1, LROTATEQ1, pQ1Length, pQ1Radius, 1);
  fEquationFZBL1 = new G4Mag_UsualEqRhs(fMagFieldFZBL1);	
  fStepperFZBL1  = new G4ClassicalRK4(fEquationFZBL1);
  fLocalFieldManagerFZBL1 = new G4FieldManager();
  fChordFinderFZBL1 = 0;
  UpdateFieldFZBL1();
  fMagFieldFZBR1 = new BField_Quad(KAPPA1, RORIGINQ1, RROTATEQ1, pQ1Length, pQ1Radius, 1);
  fEquationFZBR1 = new G4Mag_UsualEqRhs(fMagFieldFZBR1);	
  fStepperFZBR1  = new G4ClassicalRK4(fEquationFZBR1);
  fLocalFieldManagerFZBR1 = new G4FieldManager();
  fChordFinderFZBR1 = 0;
  UpdateFieldFZBR1();
  
  //Local field  FZB2, Q2
  fMagFieldFZBL2 = new BField_Quad(KAPPA2, LORIGINQ2, LROTATEQ2, pQ2Length, pQ2Radius, 2);
  fEquationFZBL2 = new G4Mag_UsualEqRhs(fMagFieldFZBL2);	
  fStepperFZBL2  = new G4ClassicalRK4(fEquationFZBL2);
  fLocalFieldManagerFZBL2 = new G4FieldManager();
  fChordFinderFZBL2 = 0;
  UpdateFieldFZBL2();
  fMagFieldFZBR2 = new BField_Quad(KAPPA2, RORIGINQ2, RROTATEQ2, pQ2Length, pQ2Radius, 2);
  fEquationFZBR2 = new G4Mag_UsualEqRhs(fMagFieldFZBR2);	
  fStepperFZBR2  = new G4ClassicalRK4(fEquationFZBR2);
  fLocalFieldManagerFZBR2 = new G4FieldManager();
  fChordFinderFZBR2 = 0;
  UpdateFieldFZBR2();
  
  fMagFieldFZBL3 = new BField_Dipole( dipoleField, LORIGIND, LROTATED );
  fEquationFZBL3 = new G4Mag_UsualEqRhs(fMagFieldFZBL3);	
  fLocalFieldManagerFZBL3 = new G4FieldManager();
  fChordFinderFZBL3 = 0;
  UpdateFieldFZBL3();
  fMagFieldFZBR3 = new BField_Dipole( dipoleField, RORIGIND, RROTATED );
  fEquationFZBR3 = new G4Mag_UsualEqRhs(fMagFieldFZBR3);	
  fLocalFieldManagerFZBR3 = new G4FieldManager();
  fChordFinderFZBR3 = 0;
  UpdateFieldFZBR3();
  
  //Local field  FZB4, Q3
  fMagFieldFZBL4 = new BField_Quad(KAPPA3, LORIGINQ3, LROTATEQ3, pQ3Length, pQ3Radius, 3);
  //fMagFieldFZBL4 = new BField_Quad(KAPPA3, ORIGINACTUAL, ROTATETEST);
  fEquationFZBL4 = new G4Mag_UsualEqRhs(fMagFieldFZBL4);	
  fStepperFZBL4  = new G4ClassicalRK4(fEquationFZBL4);
  fLocalFieldManagerFZBL4 = new G4FieldManager();
  fChordFinderFZBL4 = 0;
  UpdateFieldFZBL4();
  //Local field  FZB4, Q3
  fMagFieldFZBR4 = new BField_Quad(KAPPA3, RORIGINQ3, RROTATEQ3, pQ3Length, pQ3Radius, 3);
  //fMagFieldFZBR4 = new BField_Quad(KAPPA3, ORIGINACTUAL, ROTATETEST);
  fEquationFZBR4 = new G4Mag_UsualEqRhs(fMagFieldFZBR4);	
  fStepperFZBR4  = new G4ClassicalRK4(fEquationFZBR4);
  fLocalFieldManagerFZBR4 = new G4FieldManager();
  fChordFinderFZBR4 = 0;
  UpdateFieldFZBR4();
    
}

/////////////////////////////////////////////////////////////////////////////////
//

HRSEMFieldSetup::~HRSEMFieldSetup()
{
	if(fChordFinder) delete fChordFinder;
	if(fStepper)     delete fStepper;
	if(fEquation)    delete fEquation;
	if(fEMfield)     delete fEMfield;
	if(messenger)    delete messenger;
}

/////////////////////////////////////////////////////////////////////////////
//
// Register this field to 'global' Field Manager and
// Create Stepper and Chord Finder with predefined type, minstep (resp.)
//

void HRSEMFieldSetup::UpdateField()
{
  //SetStepper();
  fStepper = new G4ClassicalRK4( fEquation, 8 );
  //G4cout<<"HRSEMFieldSetup:: The minimal step is equal to "<<fMinStep/mm<<" mm"<<G4endl ;
  
  fFieldManager->SetDetectorField(fEMfield);
  
  if(fChordFinder) delete fChordFinder;
  fIntgrDriver = new G4MagInt_Driver(fMinStep,fStepper,fStepper->GetNumberOfVariables());
  fChordFinder = new G4ChordFinder(fIntgrDriver);
  fFieldManager->SetChordFinder( fChordFinder );
  
}


/////////////////////////////////////////////////////////////////////////////
void HRSEMFieldSetup::UpdateFieldFZBL1()
{
  //G4cout<<"HRSEMFieldSetup:: The minimal step for FZBL1 is equal to "<<fMinStep/mm<<" mm"<<G4endl ;

	fLocalFieldManagerFZBL1->SetDetectorField(fMagFieldFZBL1);

	if(fChordFinderFZBL1) delete fChordFinderFZBL1;
	fIntgrDriverFZBL1 = new G4MagInt_Driver(fMinStep,fStepperFZBL1,fStepperFZBL1->GetNumberOfVariables());
	fChordFinderFZBL1 = new G4ChordFinder((G4MagneticField*) fMagFieldFZBL1, fMinStep, fStepperFZBL1);
	fLocalFieldManagerFZBL1->SetChordFinder( fChordFinderFZBL1 );
	
}
void HRSEMFieldSetup::UpdateFieldFZBR1()
{
  //G4cout<<"HRSEMFieldSetup:: The minimal step for FZBR1 is equal to "<<fMinStep/mm<<" mm"<<G4endl ;

	fLocalFieldManagerFZBR1->SetDetectorField(fMagFieldFZBR1);

	if(fChordFinderFZBR1) delete fChordFinderFZBR1;
	fIntgrDriverFZBR1 = new G4MagInt_Driver(fMinStep,fStepperFZBR1,fStepperFZBR1->GetNumberOfVariables());
	fChordFinderFZBR1 = new G4ChordFinder((G4MagneticField*) fMagFieldFZBR1, fMinStep, fStepperFZBR1);
	fLocalFieldManagerFZBR1->SetChordFinder( fChordFinderFZBR1 );
	
}

/////////////////////////////////////////////////////////////////////////////
void HRSEMFieldSetup::UpdateFieldFZBL2()
{
  //G4cout<<"HRSEMFieldSetup:: The minimal step for FZBL2 is equal to "<<fMinStep/mm<<" mm"<<G4endl ;

	fLocalFieldManagerFZBL2->SetDetectorField(fMagFieldFZBL2);

	if(fChordFinderFZBL2) delete fChordFinderFZBL2;
	fIntgrDriverFZBL2 = new G4MagInt_Driver(fMinStep,fStepperFZBL2,fStepperFZBL2->GetNumberOfVariables());
	fChordFinderFZBL2 = new G4ChordFinder((G4MagneticField*) fMagFieldFZBL2, fMinStep, fStepperFZBL2);
	fLocalFieldManagerFZBL2->SetChordFinder( fChordFinderFZBL2 );
}
void HRSEMFieldSetup::UpdateFieldFZBR2()
{
  //G4cout<<"HRSEMFieldSetup:: The minimal step for FZBR2 is equal to "<<fMinStep/mm<<" mm"<<G4endl ;

	fLocalFieldManagerFZBR2->SetDetectorField(fMagFieldFZBR2);

	if(fChordFinderFZBR2) delete fChordFinderFZBR2;
	fIntgrDriverFZBR2 = new G4MagInt_Driver(fMinStep,fStepperFZBR2,fStepperFZBR2->GetNumberOfVariables());
	fChordFinderFZBR2 = new G4ChordFinder((G4MagneticField*) fMagFieldFZBR2, fMinStep, fStepperFZBR2);
	fLocalFieldManagerFZBR2->SetChordFinder( fChordFinderFZBR2 );
}



void HRSEMFieldSetup::UpdateFieldFZBL3()
{
  fStepperFZBL3 = new G4ClassicalRK4( fEquationFZBL3, 8 );
  //G4cout<<"HRSEMFieldSetup:: G4ClassicalRK4 (default) is called"<<G4endl;

  //G4cout<<"HRSEMFieldSetup:: The minimal step for FZBL3 is equal to "<<fMinStep/mm<<" mm"<<G4endl ;
  
  fLocalFieldManagerFZBL3->SetDetectorField(fMagFieldFZBL3);
  
  if(fChordFinderFZBL3) delete fChordFinderFZBL3;
  fChordFinderFZBL3 = new G4ChordFinder((G4MagneticField*) fMagFieldFZBL3, fMinStep, fStepperFZBL3);
  fLocalFieldManagerFZBL3->SetChordFinder( fChordFinderFZBL3 );

  fLocalFieldManagerFZBL3->SetMinimumEpsilonStep( 1.0e-5 );
  fLocalFieldManagerFZBL3->SetMaximumEpsilonStep( 1.0e-4 );
  fLocalFieldManagerFZBL3->SetDeltaOneStep( 0.5e-3 * mm );  // 0.5 micrometer
}
void HRSEMFieldSetup::UpdateFieldFZBR3()
{
  fStepperFZBR3 = new G4ClassicalRK4( fEquationFZBR3, 8 );
  //G4cout<<"HRSEMFieldSetup:: G4ClassicalRK4 (default) is called"<<G4endl;

  //G4cout<<"HRSEMFieldSetup:: The minimal step for FZBR3 is equal to "<<fMinStep/mm<<" mm"<<G4endl ;
  
  fLocalFieldManagerFZBR3->SetDetectorField(fMagFieldFZBR3);
  
  if(fChordFinderFZBR3) delete fChordFinderFZBR3;
  fChordFinderFZBR3 = new G4ChordFinder((G4MagneticField*) fMagFieldFZBR3, fMinStep, fStepperFZBR3);
  fLocalFieldManagerFZBR3->SetChordFinder( fChordFinderFZBR3 );

  fLocalFieldManagerFZBR3->SetMinimumEpsilonStep( 1.0e-5 );
  fLocalFieldManagerFZBR3->SetMaximumEpsilonStep( 1.0e-4 );
  fLocalFieldManagerFZBR3->SetDeltaOneStep( 0.5e-3 * mm );  // 0.5 micrometer
}


/////////////////////////////////////////////////////////////////////////////
void HRSEMFieldSetup::UpdateFieldFZBL4()
{
  //G4cout<<"HRSEMFieldSetup:: The minimal step for FZBL4 is equal to "<<fMinStep/mm<<" mm"<<G4endl ;

	fLocalFieldManagerFZBL4->SetDetectorField(fMagFieldFZBL4);

	if(fChordFinderFZBL4) delete fChordFinderFZBL4;
	fIntgrDriverFZBL4 = new G4MagInt_Driver(fMinStep,fStepperFZBL4,fStepperFZBL4->GetNumberOfVariables());
	fChordFinderFZBL4 = new G4ChordFinder((G4MagneticField*) fMagFieldFZBL4, fMinStep, fStepperFZBL4);
	fLocalFieldManagerFZBL4->SetChordFinder( fChordFinderFZBL4 );
}
void HRSEMFieldSetup::UpdateFieldFZBR4()
{
  //G4cout<<"HRSEMFieldSetup:: The minimal step for FZBR4 is equal to "<<fMinStep/mm<<" mm"<<G4endl ;

	fLocalFieldManagerFZBR4->SetDetectorField(fMagFieldFZBR4);

	if(fChordFinderFZBR4) delete fChordFinderFZBR4;
	fIntgrDriverFZBR4 = new G4MagInt_Driver(fMinStep,fStepperFZBR4,fStepperFZBR4->GetNumberOfVariables());
	fChordFinderFZBR4 = new G4ChordFinder((G4MagneticField*) fMagFieldFZBR4, fMinStep, fStepperFZBR4);
	fLocalFieldManagerFZBR4->SetChordFinder( fChordFinderFZBR4 );
}

/////////////////////////////////////////////////////////////////////////////
//
// Set stepper according to the stepper type
//


void HRSEMFieldSetup::SetStepper()
{
	G4int nvar = 8;

	if(fStepper) delete fStepper;

	switch ( fStepperType )
	{
	case 0:
		fStepper = new G4ExplicitEuler( fEquation, nvar );
		//G4cout<<"HRSEMFieldSetup:: G4ExplicitEuler is calledS"<<G4endl;
		break;
	case 1:
		fStepper = new G4ImplicitEuler( fEquation, nvar );
		//G4cout<<"HRSEMFieldSetup:: G4ImplicitEuler is called"<<G4endl;
		break;
	case 2:
		fStepper = new G4SimpleRunge( fEquation, nvar );
		//G4cout<<"HRSEMFieldSetup:: G4SimpleRunge is called"<<G4endl;
		break;
	case 3:
		fStepper = new G4SimpleHeum( fEquation, nvar );
		//G4cout<<"HRSEMFieldSetup:: G4SimpleHeum is called"<<G4endl;
		break;
	case 4:
		fStepper = new G4ClassicalRK4( fEquation, nvar );
		//G4cout<<"HRSEMFieldSetup:: G4ClassicalRK4 (default) is called"<<G4endl;
		break;
	case 5:
		fStepper = new G4CashKarpRKF45( fEquation, nvar );
		//G4cout<<"HRSEMFieldSetup:: G4CashKarpRKF45 is called"<<G4endl;
		break;

	//The following not working for electric field
	case 6:
		fStepper = 0 ; // new G4RKG3_Stepper( fMagEquation );
		//G4cout<<"HRSEMFieldSetup:: G4RKG3_Stepper is not currently working for Electric Field"<<G4endl;
		break;
	case 7:
		fStepper = 0 ; // new G4HelixExplicitEuler( fMagEquation );
		//G4cout<<"HRSEMFieldSetup:: G4HelixExplicitEuler is not valid for Electric Field"<<G4endl;
		break;
	case 8:
		fStepper = 0 ; //  new G4HelixImplicitEuler( fMagEquation );
		//G4cout<<"HRSEMFieldSetup:: G4HelixImplicitEuler is not valid for Electric Field"<<G4endl;
		break;
	case 9:
		fStepper = 0 ; //  new G4HelixSimpleRunge( fMagEquation );
		//G4cout<<"HRSEMFieldSetup:: G4HelixSimpleRunge is not valid for Electric Field"<<G4endl;
		break;
	default: fStepper = 0;
	}
}


///////////////////////////////////////////////////////////////////////////////
