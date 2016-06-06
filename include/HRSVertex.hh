#ifndef __HRSVERTEX_HH
#define __HRSVERTEX_HH

#include "G4ThreeVector.hh"

/*!
  Vertex information that only the
  user defined generators will see
*/

class G4Material;

class HRSVertex {
public:
  HRSVertex();
  ~HRSVertex();
  
  G4double    GetBeamE(){ return fBeamE; }
  G4double    GetRadLen(){ return fRadLen; }
  G4Material *GetMaterial(){ return fMaterial; }
  
public:
  G4double fBeamE;
  G4double fRadLen;
  G4double fmsth;
  G4double fmsph;
  G4double XS_208Pb;
  G4Material *fMaterial;
  G4int pass;
};

#endif//__HRSVERTEX_HH
