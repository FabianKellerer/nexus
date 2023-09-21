#ifndef PMTCALIB_H
#define PMTCALIB_H

#include <G4ThreeVector.hh>
#include "GeometryBase.h"

class G4Material;
class G4GenericMessenger;
namespace nexus
{
    class CylinderPointSampler;
    class GenericPhotosensor;
    class PMTcalib:public GeometryBase
    {
        public:
        //constructor
        PMTcalib();
        //destructor
        ~PMTcalib();

        //generate vertex
        G4ThreeVector GenerateVertex() const;
        void Construct();

        private:
        G4double length_;     //length of the cylindrical optical fibre

        G4LogicalVolume* pmt_logic_;
    };
    
} // namespace nexus
#endif