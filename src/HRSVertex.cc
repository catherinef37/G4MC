#include "HRSVertex.hh"

#include "G4Material.hh"

HRSVertex::HRSVertex(){
    // Some default material
    fMaterial = NULL;
    fBeamE   = 0.0*GeV;
    fRadLen  = 0.0;
    fmsth    = 0.0;
    fmsph    = 0.0;
    XS_208Pb = 0.0;
    pass     = 0;
}

HRSVertex::~HRSVertex(){
}
