#ifndef BFIELD_QUAD_HH
#define BFIELD_QUAD_HH

#include "G4MagneticField.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
//#include "QuadFringe.hh"
#include "BField_Septum.hh"
class BField_Quad : public G4MagneticField
{
  public: // with description

  BField_Quad(G4double          pGradient, 
	      G4ThreeVector     pOrigin, 
	      G4RotationMatrix* pMatrix,
	      //QuadFringe*       pFringe,
	      G4double          pLength,
  	      G4double          pRadius,
	      G4int             pQuadNumber);
  
  ~BField_Quad();
  
  void GetFieldValue(const G4double yTrack[],
		     G4double B[]     ) const;
  G4ThreeVector GetFringeField(G4ThreeVector) const;
private:
  
  G4double          fGradient;
  G4ThreeVector     fOrigin;
  G4RotationMatrix* fpMatrix;
  //QuadFringe*       fQuadFringe;
  G4double          fLength;
  G4double          fRadius;
  G4int             fSnakeModel;
  G4int             fQuadNumber;
  BField_Septum*    mBField_Septum; 
  G4double          pLHRSMomentum;
  G4double          pRHRSMomentum;

};
#endif

