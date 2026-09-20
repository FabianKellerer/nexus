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
        bool issquare(G4int n) const;
        G4ThreeVector ComputeFiberPositions();

        private:
        G4double radius_;     //radius of the cylindrical optical fibre
        G4double length_;     //length of the cylindrical optical fibre
        G4double fiber_dist_; //distance between the fibers
        G4bool   al_;         //aluminises fibers
        G4bool   tefl_;       //puts teflon block behind fibers
        G4bool   isround_;    //if false: square fibers
        G4bool   doubleclad_; //if true: double cladding, if false: single cladding
        G4double thickness_;  //thickness of the photosensor
        G4String core_mat_;   //core material of the fibre (EJ280, EJ286 or Y11)
        G4String sensortype_; //photosensor type (PMT at fibre end or SiPM opposite light source)
        G4String sensor_pos_;  //position of the sensor (end or side)
        G4int    num_fibers_; //number of fibres in the bundle
        G4double lamp_size_;  //size of photon generation region in z direction if LAMP region is chosen
        G4double gap_;        //size of gap between fiber and sensor
        G4double rand_wls_;   //number of sigmas of the WLS uncertainty to shift the WLS abs. length by
        G4double rand_att_;   //number of sigmas of the attenuation uncertainty to shift
        G4double rand_sigma_; //number of sigmas of the radial distribution uncertainty to shift the radial distribution by
        G4double qe_;          //quantum efficiency of the WLS

        GenericPhotosensor*   sensor_;
        G4GenericMessenger*   msg_;
        CylinderPointSampler* cyl_vertex_gen_;
        G4LogicalVolume* pmt_logic_;
    };
    
} // namespace nexus
#endif
