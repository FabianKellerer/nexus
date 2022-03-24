// ----------------------------------------------------------------------------
// nexus | XeBox.h
//
// Box filled with xenon, with single photosensor.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef XE_BOX_H
#define XE_BOX_H

#include <G4ThreeVector.hh>
#include "GeometryBase.h"

class G4Material;
class G4GenericMessenger;
namespace nexus { class BoxPointSampler; }


namespace nexus {

  /// Box-shaped chamber filled with xenon (liquid or gas)
  class GenericPhotosensor;

  class XeBox: public GeometryBase
  {
  public:
    /// Constructor
    XeBox();
    /// Destructor
    ~XeBox();

    /// Return vertex within region <region> of the chamber
    G4ThreeVector GenerateVertex(const G4String& region) const;

    void Construct();

  private:
    G4bool   liquid_;   ///< Whether xenon is liquid or not
    G4double pressure_; ///< Pressure (if gaseous state was selected)
    G4double length_;   ///< Length of the box
    G4double width_;    ///< Width  of the box
    G4double height_;   ///< Height of the box
    G4double thickness_;///< thickness of the sensor

    GenericPhotosensor* sensor_;

    /// Vertexes random generator
    BoxPointSampler* box_vertex_gen_;

    /// Messenger for the definition of control commands
    G4GenericMessenger* msg_;
  };

} // end namespace nexus

#endif
