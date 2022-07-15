// ----------------------------------------------------------------------------
// nexus | OpticalFibre.h
//
// Cylindrical optical fibre to be placed on the walls to increase light collection,
// with single photosensor.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef OPTICAL_FIBRE_H
#define OPTICAL_FIBRE_H

#include <G4ThreeVector.hh>
#include "GeometryBase.h"

class G4Material;
class G4GenericMessenger;
namespace nexus
{
    class CylinderPointSampler;
    class GenericPhotosensor;
    class OpticalFibre:public GeometryBase
    {
        public:
        //constructor
        OpticalFibre();
        //destructor
        ~OpticalFibre();

        //generate vertex
        G4ThreeVector GenerateVertex(const G4String& region) const;
        void Construct();

        private:
        G4double radius_;    //radius of the cylindrical optical fibre
        G4double length_;    //length of the cylindrical optical fibre
        G4double thickness_; //thickness of the photosensor

        GenericPhotosensor*   sensor_;
        G4GenericMessenger*   msg_;
        CylinderPointSampler* cyl_vertex_gen_;
    };
    
} // namespace nexus
#endif
