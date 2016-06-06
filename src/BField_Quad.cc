//#include "Python.h"
#include "BField_Quad.hh"
#include "BField_Septum.hh"
#include "G4RotationMatrix.hh"
//#include "QuadFringe.hh"
#include "UsageManager.hh"
//#define DEBUG_BFIELD_QUAD 0

static G4RotationMatrix IdentityMatrix; 

extern UsageManager* gConfig;

/////////////////////////////////////////////////////////////////////////

BField_Quad::BField_Quad(G4double pGradient, G4ThreeVector
			 pOrigin, G4RotationMatrix* pMatrix, G4double pLength, G4double pRadius, G4int pQuadNumber)
//pOrigin, G4RotationMatrix* pMatrix, QuadFringe* pQuadFringe, G4double pLength, G4double pRadius, G4int pQuadNumber)
{
   fGradient    = pGradient ;
   fOrigin      = pOrigin ;
   fpMatrix     = pMatrix ;
   //fQuadFringe  = pQuadFringe;
   fLength      = pLength;
   fRadius      = pRadius;
   fQuadNumber  = pQuadNumber;
   gConfig->GetArgument("SnakeModel",fSnakeModel);

   string pTargetFieldIni=gConfig->GetArgument("TargetFieldIni");
   string pTargetFieldMap=gConfig->GetArgument("TargetFieldMap");
   gConfig->GetArgument("LHRSMomentum",pLHRSMomentum);
   gConfig->GetArgument("RHRSMomentum",pRHRSMomentum);
   string pSeptumFieldIni=gConfig->GetArgument("SeptumFieldIni");
   string pSeptumFieldMap=gConfig->GetArgument("SeptumFieldMap");
   mBField_Septum = new BField_Septum(pLHRSMomentum,pRHRSMomentum,
				      pSeptumFieldIni.c_str(),pSeptumFieldMap.c_str());


}

/////////////////////////////////////////////////////////////////////////

BField_Quad::~BField_Quad()
{
}

////////////////////////////////////////////////////////////////////////

void BField_Quad::GetFieldValue( const G4double y[7],
				 G4double B[3]  ) const  
{
  G4int sos = 1;
  if( fSnakeModel == 52 )
    sos = 1;
  G4ThreeVector r_global= G4ThreeVector
    (y[0] - fOrigin.x(), 
     y[1] - fOrigin.y(), 
     y[2] - fOrigin.z());
  
  G4ThreeVector r_local = G4ThreeVector
    (fpMatrix->colX() * r_global,
     fpMatrix->colY() * r_global,
     fpMatrix->colZ() * r_global);

  
  G4ThreeVector B_local_old = G4ThreeVector
    (fGradient * r_local.y(),
     fGradient * r_local.x(),
     0);
  
  //G4cout << "Fringe Field Subroutine beginning!" << G4endl;
  //G4cout << "The field is " << fGradient << " " << tesla << G4endl;
  //G4ThreeVector B_local = (G4ThreeVector)fQuadFringe->GetFieldValue( r_local.x() / 1000.0 , r_local.y() / 1000.0 , r_local.z() / 1000.0 , fGradient * 1000. , fLength / 1000. , fRadius / 1000. ); //python code is in m, and GEANT4 is in mm//FRINGE FIELDS
  G4ThreeVector B_local;  
  if( sos && ( fQuadNumber == 1 ) ){
    /*The fringe field coeficients for the SOS quad are as follows:
      a1       a2       a3         a4    a5    a6
      0.2486 5.3626 -2.412  0.9874  0.0  0.0
      F(z) = 1 / (1 + exp(a1 +a2*(z/D) + ... + a6*(z/D)^5))
      where F(z) is the enge function, z is the distance to the effective field boundary, and D = 2*radius. (See attached plot).*/
    G4double fringe = 0.;
    G4double D      = fRadius * 2.;
    //G4cout << "The radius is: " << fRadius << G4endl;
    G4double myz    = ( fabs( r_local.z() ) - fLength * 0.5 );
    G4double a[6]   = {0.2486, 5.3626, -2.412,
		       0.9874, 0.0   ,  0.0  };
    for(int i = 0; i < 6; i++){
      fringe += a[i] * pow( myz / D, i );
    }
    fringe = 1. / ( 1. + exp( fringe ) );
    //G4double fringe = 1 / (1 + exp( a1 + a2 * ( z / D ) + ... + a6*(z/D)^5));
    //B_local = G4ThreeVector( r_local.y() * fGradient / fRadius * fringe ,  r_local.x() * fGradient / fRadius * fringe, 0.0 ); //YES FRINGE FIELDS
    if( myz < 0. ){
      B_local = G4ThreeVector( r_local.y() * fGradient / fRadius,  r_local.x() * fGradient / fRadius, 0.0 ); //NO FRINGE FIELDS
    }else{
      B_local = G4ThreeVector( 0., 0., 0. );
    }
    //G4cout << B_local << G4endl;
    //G4cout << fQuadNumber << " " << B_local << G4endl;
  }else{
    if( fabs( r_local.z() ) - fLength * 0.5 < 0 ){
      B_local = G4ThreeVector( r_local.y() * fGradient / fRadius ,  r_local.x() * fGradient / fRadius, 0.0 ); //NO FRINGE FIELDS
    }else{
      B_local = G4ThreeVector( 0.0, 0.0, 0.0 ); //NO FRINGE FIELDS
    }
    //G4cout << fQuadNumber << " " << B_local << G4endl;
  }
  //G4cout << "Fringe Field Subroutine Finished!" << G4endl;
  //G4cout << "r    : " << r_local.x() / 1000.0 << " " << r_local.y() / 1000.0 << " " <<  r_local.z() / 1000.0 << G4endl;
  //G4cout << "B new: " << B_local.x()          << " " << B_local.y()          << " " <<  B_local.z()          << G4endl;
  //G4cout << "B old: " << B_local_old.x()      << " " << B_local_old.y()      << " " <<  B_local_old.z()      << G4endl;
  G4ThreeVector B_global = G4ThreeVector
    (fpMatrix->rowX() * B_local,
     fpMatrix->rowY() * B_local,
     fpMatrix->rowZ() * B_local);
  
  B[0] = B_global.x() ;
  B[1] = B_global.y() ;
  B[2] = B_global.z() ;

  //B[0] = 0.0; //WARNING, THESE LINES TURN OFF THE QUAD FIELD
  //B[1] = 0.0;
  //B[2] = 0.0;

  

#ifdef DEBUG_BFIELD_QUAD

  /*
  G4cout << "SOS: " << sos << G4endl;
  G4cout << "Global coordinates:" << G4endl;
  G4cout << y[0] << " " << y[1] << " " << y[2] << G4endl;
  G4cout << "Origin coordinates:" << G4endl;
  G4cout << fOrigin.x()  << " " << fOrigin.y()  << " " << fOrigin.z()  << G4endl;
  G4cout << "Rotation matrix:" << G4endl;
  G4cout << fpMatrix->xx() << " " << fpMatrix->xy() << " " << fpMatrix->xz() << G4endl;
  G4cout << fpMatrix->yx() << " " << fpMatrix->yy() << " " << fpMatrix->yz() << G4endl;
  G4cout << fpMatrix->zx() << " " << fpMatrix->zy() << " " << fpMatrix->zz() << G4endl;
  //G4cout << "Transpose matrix:" << G4endl;
  //G4cout << fpMatrix->colX().x() << " " << fpMatrix->colX().y() << " "<< fpMatrix->colX().z() << G4endl;
  //G4cout << fpMatrix->colY().x() << " " << fpMatrix->colY().y() << " "<< fpMatrix->colY().z() << G4endl;
  //G4cout << fpMatrix->colZ().x() << " " << fpMatrix->colZ().y() << " "<< fpMatrix->colZ().z() << G4endl;
  //G4cout << "Rotation matrix:" << G4endl;
  //G4cout << fpMatrix->rowX().x() << " " << fpMatrix->rowX().y() << " "<< fpMatrix->rowX().z() << G4endl;
  //G4cout << fpMatrix->rowY().x() << " " << fpMatrix->rowY().y() << " "<< fpMatrix->rowY().z() << G4endl;
  //G4cout << fpMatrix->rowZ().x() << " " << fpMatrix->rowZ().y() << " "<< fpMatrix->rowZ().z() << G4endl;
  */
  G4cout << "Local position:" << G4endl;
  G4cout << r_local .x() << " " << r_local .y() << " " << r_local .z() << G4endl;
  G4cout << "Local field:" << G4endl;
  G4cout << B_local .x() << " " << B_local .y() << " " << B_local .z() << G4endl;
  G4cout << "Global field:" << G4endl;
  G4cout << B_global.x() << " " << B_global.y() << " " << B_global.z() << G4endl;
#endif

}
