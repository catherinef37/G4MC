//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: BField_Dipole.cc 69786 2013-05-15 09:38:51Z gcosmo $
//
// -------------------------------------------------------------------

#include "BField_Dipole.hh"
#include "G4RotationMatrix.hh"
#include <math.h>
//#define DEBUG_BFIELD_DIPOLE 0

extern"C" {
  void dfringe_(float* x,    float* y,  float* z,
                float* fact, float* a,
                float* fx,   float* fy, float* fz);
  void snakedipole_(float* x,  float* y,  float* z,
		    float* bx, float* by, float* bz,
		    float* b0);
}

static G4RotationMatrix IdentityMatrix; 

BField_Dipole::BField_Dipole(G4double pBend)
{
   fBend    = pBend ;
   fOrigin  = G4ThreeVector( 0.0, 0.0, 0.0) ;
   fpMatrix = &IdentityMatrix;
}

/////////////////////////////////////////////////////////////////////////

BField_Dipole::BField_Dipole(G4double pBend, G4ThreeVector
pOrigin, G4RotationMatrix* pMatrix)
{
   fBend    = pBend   ;
   fOrigin  = pOrigin ;
   fpMatrix = pMatrix ;
}

/////////////////////////////////////////////////////////////////////////

BField_Dipole::~BField_Dipole()
{
}

////////////////////////////////////////////////////////////////////////
//  Allow displaced origin and rotation 
//  Extensions by Björn Riese (GSI)

void BField_Dipole::GetFieldValue( const G4double y[7],
				 G4double B[3]  ) const  
{
  //This old way doesn't work, I believe.
  //G4ThreeVector r_global= G4ThreeVector
  //(y[0], 
  //y[1],
  //y[2]);
  
  //G4ThreeVector r_local = G4ThreeVector
  //(fpMatrix->colX() * r_global - fOrigin.x(),
  //fpMatrix->colY() * r_global - fOrigin.y(),
  //fpMatrix->colZ() * r_global - fOrigin.z());

  G4ThreeVector r_global= G4ThreeVector
    (y[0] - fOrigin.x(),
     y[1] - fOrigin.y(),
     y[2] - fOrigin.z());

  G4ThreeVector r_local = G4ThreeVector
    (fpMatrix->colX() * r_global,
     fpMatrix->colY() * r_global,
     fpMatrix->colZ() * r_global);

  //indexed dipole
  //B = ( r / r_0 ) ^ -n  * B_0
  //G4double B_0 = fBend;
  float B_0 = -(float)fBend / tesla;
  //  G4cout << B_0 << " is the B_0 in G4MC" << G4endl;
  //G4double r   = sqrt( pow(r_local.z(), 2) + pow(r_local.y(), 2) );
  //G4double r_0 = 8.4 * m;
  //G4double n   =-1.25;
  //G4double n   =-1.263;
  //G4double indexed_field = -pow( r / r_0, -n ) * B_0;
  // y is -x and x is y
  //G4double indexed_field_x = ( 1. - n * ( r - r_0 ) / r_0  - 0.5 * n * r_local.x() * r_local.x()  / r_0 / r_0 ) * B_0;
  //G4double indexed_field_x = ( 1. - n * ( r - r_0 ) / r_0 ) * B_0;
  //G4double indexed_field_y = ( 0. + n * r_local.x() / r_0 ) * B_0;
  //G4double chi = 0.;
  //if( r_local.z() < 0.){
  //chi = -atan( r_local.z() / abs(r_local.y()) );
  //}else{
  //chi = atan( r_local.z() / abs(r_local.y()) );
  //}

  float snx = r_local.x();
  float sny = r_local.y();
  float snz = r_local.z();
  float snbx,snby,snbz;
  //G4cout << "input " << snx << " " << sny << " " << snz << G4endl;
  snakedipole_(&snx,  &sny,  &snz,
	       &snbx, &snby, &snbz,
	       &B_0);
  //  G4cout << snbx << " " << snby << " " << snbz << G4endl;
  snbx *= -tesla;
  snby *= -tesla;
  snbz *= -tesla;
  //These lines are original without the fringe field, overwrite them later, if you desire fringe field.
  //straight from SNAKE
  //if( mtyp .eq. 3) by=
  //1     bf* ( 1. - ndx*drr1 + bet1*drr2 + gama*drr3 + delt*drr4 )
  //if( mtyp .eq. 4) by= bf/ (1. + ndx*drr1 )
  //G4double indexed_field_x = ( 1. + n * ( r_0 - r ) / r_0 + 0.5 * n * r_local.x() * r_local.x() / r_0 / r_0 ) * B_0;
  //G4double indexed_field_y = ( 0. + n * r_local.x() / r_0 ) * B_0 * cos( chi );
  //G4double indexed_field_z = ( 0. + n * r_local.x() / r_0 ) * B_0 *-sin( chi );
  /*
    float dnr1 = 1. + n * ( r - r_0 ) / r_0;
  float dnr2 = dnr1 * dnr1;
  float dnr3 = dnr2 * dnr1;
  float dnr4 = dnr3 * dnr1;
  float dnr5 = dnr4 * dnr1;
  float yr1 = r_local.x() / r_0;
  float yr2 = yr1 * yr1;
  float yr3 = yr2 * yr1;
  float yr4 = yr3 * yr1;
  float rr1 = r_0 / r;
  float rr2 = rr1 * rr1;
  float rr3 = rr2 * rr1;
  float snake_b = B_0 * ( 1. / dnr1
			  + .5 * yr2 * n * ( -2. * n / dnr3 + rr1 / dnr2 ) +
			  yr4 * n * ( 24. * n * n * n / dnr5 - 12. * n * n * rr1 / dnr4 - 2. * n * rr2 / dnr3 - rr3 / dnr2 ) / 24. );
*/
  //G4double indexed_field_x = ( 1. / ( 1 + n * ( r - r_0 ) / r_0 )  ) * B_0;
  //G4double indexed_field_x = snake_b;
  //G4double indexed_field_y = ( 0. + n * r_local.x() / r_0 ) * ( 1. / ( 1 + n * ( r - r_0 ) / r_0 )  ) * ( 1. / ( 1 + n * ( r - r_0 ) / r_0 )  ) * B_0 * cos( chi );
  //G4double indexed_field_z = ( 0. + n * r_local.x() / r_0 ) * ( 1. / ( 1 + n * ( r - r_0 ) / r_0 )  ) * ( 1. / ( 1 + n * ( r - r_0 ) / r_0 )  ) * B_0 *-sin( chi );
  /*
  G4cout << "r, r_0, x, y, z (local)" << G4endl;
  G4cout << r << " " << r_0 << " " << r_local.x() << " " << r_local.y() << " " << r_local.z() << G4endl;
  G4cout << "chi (rad), chi (deg), Bx, By, Bz (local)" << G4endl;
  G4cout << chi << " " << chi * 180. / 3.141592654 << " " << indexed_field_x << " " << indexed_field_y << " " << indexed_field_z << G4endl;
  */

  /*
  float fact      = (float)B_0 / tesla;
  float aperture = 27;//in cm
  float fx        = 0;
  float fy        = 0;
  float fz        = 0;
  //if(0){// turn it off
  //Note: In the local dipole coordinates, all boundaries are defined completely in terms of y and z. x never matters!
  float z1 = - 1.5 * m; //The leftmost limit, where you say: it's just about 0, so save computing time and call it 0
  float z2 = ( r_local.y() + r_0 ) * tan( pi / 6. );//The entrance of the dipole
  float chi3 = pi / 16.;//When you are inside the magnet far enough that it is "full strength" with no fringe effect. We will just use an angle for ease of calculation
  float chi4 = pi * 3. / 16.;//The same as before, but you are appraoching the exit, still sufficiently inside the magnet, though
  float z5 = ( r_local.y() + r_0 * cos( pi / 4. ) ) / -tan( pi / 12. ) + r_0 * sin( pi / 4. );//The exit of the dipole
  float z6 =   8.0 * m; //The rightmost limit, just call it zero
  if( r_local.z() < z1 || r_local.z() >= z6){
    indexed_field_x = 0.;
    indexed_field_y = 0.;
    indexed_field_z = 0.;
    //do nothing
  }else if( r_local.z() >= z1 && r_local.z() < z2 ){
    float snake_x = - r_local.y() - r_0;
    float snake_y = - r_local.z() + ( r_local.y() + r_0 ) * tan( pi / 6 );
    float snake_z =   r_local.x();
    snake_x /= 10.;
    snake_y /= 10.;
    snake_z /= 10.;
    dfringe_( &snake_x, &snake_y , &snake_z,
              &fact   , &aperture,
              &fx     , &fy      , &fz);
    indexed_field_x =   fz * tesla;//snake
    indexed_field_y = - fx * tesla;//snake
    indexed_field_z = - fy * tesla;//snake
  }else if( r_local.z() >= z2 && chi < chi3  ){//This is where I will cut off the fringe field, for computational reasons.
                      //Also, this is only for the entrance. You have to program the exit separately.
    float theta_i = asin( sin( pi / 6. ) * r_0 / r ) - pi / 6.; 
    float theta_f = chi;
    float d_theta = theta_f - theta_i;
    //float snake_x   = - r_local.y() - r_0;
    float snake_x   = ( r - r_0 );  //-y
    float snake_y   = - r * d_theta;//-z
    float snake_z   =   r_local.x();//+x
    snake_x /= 10.;
    snake_y /= 10.;
    snake_z /= 10.;
    
    dfringe_(&snake_x, &snake_y  , &snake_z,
	     &fact   , &aperture,
	     &fx     , &fy       , &fz    );
    //G4cout << "Fringing fields: " << fx << " " << fy << " " << fz << G4endl;
    //For now, the dipole term of the dipole, NOT the quad (indexed) terms of the dipole, will have fringe field.
    //The quad behavior of the index dipole needs a separate treatment, obviously, but it is supressed by a factor of 1/r or 1/r^2.
    //That makes fringe fields due to those a refinement of a refinement, so... very small.
    //indexed_field_x = (  fz * tesla ) + ( n * ( r_0 - r ) / r_0 + 0.5 * n * r_local.x() * r_local.x() / r_0 / r_0 ) * B_0;
    //indexed_field_y = ( -fx * tesla ) * cos( chi ) + ( n * r_local.x() / r_0 ) * B_0 * cos( chi );
    //indexed_field_z = ( -fy * tesla ) *-sin( chi ) + ( n * r_local.x() / r_0 ) * B_0 *-sin( chi );
    indexed_field_x = (  fz * tesla ) - B_0 + ( 1. / ( 1 + n * ( r - r_0 ) / r_0 )  ) * B_0;
    indexed_field_y = ( -fx * tesla ) * cos( chi ) + ( 0. + n * r_local.x() / r_0 ) * ( 1. / ( 1 + n * ( r - r_0 ) / r_0 )  ) * ( 1. / ( 1 + n * ( r - r_0 ) / r_0 )  ) * B_0 * cos( chi );
    indexed_field_z = ( -fy * tesla ) *-sin( chi ) + ( 0. + n * r_local.x() / r_0 ) * ( 1. / ( 1 + n * ( r - r_0 ) / r_0 )  ) * ( 1. / ( 1 + n * ( r - r_0 ) / r_0 )  ) * B_0 *-sin( chi );
  }else if( chi >= chi3 && chi < chi4 ){
    //do nothing
  }else if( chi >= chi4 && r_local.z() < z5 ){
    float theta_f = pi / 4. - ( asin( sin( pi / 6. ) * r_0 / r ) - pi / 6.); 
    float theta_i = chi;
    float d_theta = theta_f - theta_i;
    //float snake_x   = - r_local.y() - r_0 * m;
    float snake_x   = ( r - r_0 );  //-y
    float snake_y   = - r * d_theta;//+z
    float snake_z   = - r_local.x();//-x
    snake_x /= 10.;
    snake_y /= 10.;
    snake_z /= 10.;
    
    dfringe_(&snake_x, &snake_y  , &snake_z,
	     &fact   , &aperture,
	     &fx     , &fy       , &fz    );
    
    //dipole exit, my x and z signs flip.
    //indexed_field_x = (  fz * tesla ) + ( n * ( r_0 - r ) / r_0 + 0.5 * n * r_local.x() * r_local.x() / r_0 / r_0 ) * B_0;
    //indexed_field_y = ( -fx * tesla ) * cos( chi ) + ( n * r_local.x() / r_0 ) * B_0 * cos( chi );
    //indexed_field_z = ( -fy * tesla ) *-sin( chi ) + ( n * r_local.x() / r_0 ) * B_0 *-sin( chi );
    indexed_field_x = (  fz * tesla ) - B_0 + ( 1. / ( 1 + n * ( r - r_0 ) / r_0 )  ) * B_0;
    indexed_field_y = ( -fx * tesla ) * cos( chi ) + ( 0. + n * r_local.x() / r_0 ) * ( 1. / ( 1 + n * ( r - r_0 ) / r_0 )  ) * ( 1. / ( 1 + n * ( r - r_0 ) / r_0 )  ) * B_0 * cos( chi );
    indexed_field_z = ( -fy * tesla ) *-sin( chi ) + ( 0. + n * r_local.x() / r_0 ) * ( 1. / ( 1 + n * ( r - r_0 ) / r_0 )  ) * ( 1. / ( 1 + n * ( r - r_0 ) / r_0 )  ) * B_0 *-sin( chi );
  }else if( r_local.z() >= z5 && r_local.z() < z6 ){
    float tilt    = pi / 4.;
    float snake_x =   ( r_local.z() * sin( tilt ) - r_local.y() * cos( tilt ) - r_0 ) / cos( pi / 6 );                     //-y
    float snake_y =   r_local.z() * cos( tilt ) + r_local.y() * sin( tilt ) + ( r_0 + r_local.y() * cos( tilt ) - r_local.z() * sin( tilt ) ) * tan( pi / 6 ); //+z
    float snake_z = - r_local.x();                                                                     //-x
    snake_x /= 10.;
    snake_y /= 10.;
    snake_z /= 10.;
    dfringe_( &snake_x, &snake_y , &snake_z,
              &fact   , &aperture,
              &fx     , &fy      , &fz);
    indexed_field_x =   fz * tesla;//snake
    indexed_field_y = - fx * tesla * cos( tilt );//snake
    indexed_field_z = - fy * tesla *-sin( tilt );//snake

  }else{
    indexed_field_x = 0.;
    indexed_field_y = 0.;
    indexed_field_z = 0.;
  }
  indexed_field_y = 0.;
  indexed_field_z = 0.;
  */

  G4ThreeVector B_local = G4ThreeVector( snbx, snby, snbz );
  
  G4ThreeVector B_global = G4ThreeVector
    (fpMatrix->rowX() * B_local,
     fpMatrix->rowY() * B_local,
     fpMatrix->rowZ() * B_local);
  
  B[0] = B_global.x() ;
  B[1] = B_global.y() ;
  B[2] = B_global.z() ;
  /*
    B[0] = 0.;
    B[1] = 0.;
    B[2] = 0.;
  */
#ifdef DEBUG_BFIELD_DIPOLE
  G4cout << "Dipole: global coordinates, local, global field: " << G4endl;
  G4cout << r_global.x() << " " << r_global.y() << " " << r_global.z() << " " << G4endl;
  G4cout << r_local.x() << " " << r_local.y() << " " << r_local.z() << " " << G4endl;
  G4cout << B_local .x() << " " << B_local .y() << " " << B_local .z() << G4endl;
  G4cout << B_global .x() << " " << B_global .y() << " " << B_global .z() << G4endl;
#endif

}
