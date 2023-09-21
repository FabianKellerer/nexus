#ifndef FIBER_ABSORPTION_H
#define FIBER_ABSORPTION_H

#include <G4ThreeVector.hh>
#include "GeometryBase.h"
#include "PmtR7378A.h"

class G4Material;
class G4GenericMessenger;
namespace nexus
{
    class CylinderPointSampler;
    class GenericPhotosensor;
    class FiberAbsorption:public GeometryBase
    {
        public:
        //constructor
        FiberAbsorption();
        //destructor
        ~FiberAbsorption();

        //generate vertex
        G4ThreeVector GenerateVertex(const G4String& region) const;
        void Construct();

        private:
        G4double lambdawls_;

        PmtR7378A*   sensor_;
        G4GenericMessenger*   msg_;
        CylinderPointSampler* cyl_vertex_gen_;
        G4LogicalVolume* pmt_logic_;
    };
    
} // namespace nexus
#endif