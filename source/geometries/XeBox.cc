// ----------------------------------------------------------------------------
// nexus | XeBox.cc
//
// Box filled with xenon, with single photosensor.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------



#include "XeBox.h"

#include "BoxPointSampler.h"
#include "MaterialsList.h"
#include "IonizationSD.h"
#include "OpticalMaterialProperties.h"
#include "FactoryBase.h"
#include "Visibilities.h"

#include <G4GenericMessenger.hh>
#include <G4Box.hh>
#include <G4NistManager.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4Material.hh>
#include <GenericPhotosensor.h>
#include <G4VisAttributes.hh>
#include <G4SDManager.hh>
#include <G4VUserDetectorConstruction.hh>

#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>


using namespace nexus;
using namespace CLHEP;

REGISTER_CLASS(XeBox, GeometryBase)

namespace nexus {

  XeBox::XeBox():
    GeometryBase(), liquid_(true), pressure_(STP_Pressure),
    length_(1.*m), width_(1.*m), height_(1.*m), box_vertex_gen_(0)
  {
    msg_ = new G4GenericMessenger(this, "/Geometry/XeBox/",
      "Control commands of geometry XeBox.");

    msg_->DeclareProperty("LXe", liquid_,
      "Build the box with liquid xenon.");

    G4GenericMessenger::Command& pressure_cmd =
      msg_->DeclareProperty("pressure", pressure_,
      "Set pressure for gaseous xenon (if selected).");
    pressure_cmd.SetUnitCategory("Pressure");
    pressure_cmd.SetParameterName("pressure", false);
    pressure_cmd.SetRange("pressure>0.");

    G4GenericMessenger::Command& length_cmd =
      msg_->DeclareProperty("length", length_, "Length of the xenon box.");
    length_cmd.SetUnitCategory("Length");
    length_cmd.SetParameterName("length", false);
    length_cmd.SetRange("length>0.");
    
    G4GenericMessenger::Command& width_cmd =
      msg_->DeclareProperty("width", width_, "Width of the xenon box.");
    width_cmd.SetUnitCategory("Length");
    width_cmd.SetParameterName("width", false);
    width_cmd.SetRange("width>0.");
    
    G4GenericMessenger::Command& height_cmd =
      msg_->DeclareProperty("height", height_, "Height of the xenon box.");
    height_cmd.SetUnitCategory("Length");
    height_cmd.SetParameterName("height", false);
    height_cmd.SetRange("height>0.");

    // Create a vertex generator for a box
    box_vertex_gen_ = new BoxPointSampler(length_, width_, height_,  0., G4ThreeVector(0., 0., 0.), 0);

    // hardcoded thickness of the photosensor
    thickness_=2.*mm;
  }



  XeBox::~XeBox()
  {
    delete box_vertex_gen_;
    delete msg_;
  }



  void XeBox::Construct()
  {
    G4String name = "XE_BOX";

    // LAB. This is just a volume of air surrounding the detector
    G4Box* lab_solid = new G4Box("LAB", (length_+1.*cm)/2., (width_+1.*cm)/2., (height_+1.*cm)/2.);

    G4LogicalVolume* lab_logic =
      new G4LogicalVolume(lab_solid,
                          G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"),
                          "LAB");
    lab_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    this->SetLogicalVolume(lab_logic);

    // Define solid volume as a box
    G4Box* Box_solid = new G4Box(name, length_/2., width_/2., height_/2.);

    // Define the material (LXe or GXe) for the box.
    // We use for this the NIST manager or the nexus materials list.
    G4Material* xenon = 0;
    if (liquid_){
      xenon = G4NistManager::Instance()->FindOrBuildMaterial("G4_lXe");
      xenon->SetMaterialPropertiesTable(opticalprops::LXe());}
    else{
      xenon = materials::GXe(pressure_);
      xenon->SetMaterialPropertiesTable(opticalprops::GXe());
    }
    // Define the logical volume of the box using the material
    // and the solid volume defined above
    G4LogicalVolume* box_logic =
    new G4LogicalVolume(Box_solid, xenon, name);
    GeometryBase::SetLogicalVolume(box_logic);
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), box_logic,
		      name, lab_logic, false, 0, false);

    // Set the logical volume of the box as an ionization
    // sensitive detector, i.e. position, time and energy deposition
    // will be stored for each step of any charged particle crossing
    // the volume.
    IonizationSD* ionizsd = new IonizationSD("/XE_BOX");
    G4SDManager::GetSDMpointer()->AddNewDetector(ionizsd);
    box_logic->SetSensitiveDetector(ionizsd);

    //Build the sensor
    sensor_  = new GenericPhotosensor("SENSOR", length_, width_, thickness_);
    sensor_ -> SetVisibility(false);

    //Set the sensor window material
    G4Material* window_mat_ = materials::TPB();
    window_mat_->SetMaterialPropertiesTable(opticalprops::TPB());
    G4MaterialPropertyVector* window_rindex = window_mat_->GetMaterialPropertiesTable()->GetProperty("RINDEX");
    sensor_ -> SetWindowRefractiveIndex(window_rindex);

    //Set the optical properties of the sensor
    G4MaterialPropertiesTable* photosensor_mpt = new G4MaterialPropertiesTable();
    G4double energy[]       = {0.2 * eV, 3.5 * eV, 3.6 * eV, 11.5 * eV};
    G4double reflectivity[] = {0.0     , 0.0     , 0.0     ,  0.0     };
    G4double efficiency[]   = {1.0     , 1.0     , 0.0     ,  0.0     };
    photosensor_mpt->AddProperty("REFLECTIVITY", energy, reflectivity, 4);
    photosensor_mpt->AddProperty("EFFICIENCY",   energy, efficiency,   4);
    sensor_->SetOpticalProperties(photosensor_mpt);

    //Set sensor depth and naming order
    sensor_ ->SetSensorDepth(1);
    sensor_ ->SetMotherDepth(2);
    sensor_ ->SetNamingOrder(1);

    sensor_ -> Construct();

    //Placing the sensor
    G4LogicalVolume* sensor_logic = sensor_ -> GetLogicalVolume();
    G4RotationMatrix sensor_rot;
    //sensor_rot.rotateY(pi);
    G4ThreeVector sensor_pos = G4ThreeVector(0,
                                             0,
                                             height_);

    new G4PVPlacement(G4Transform3D(sensor_rot, sensor_pos), sensor_logic,
                        sensor_logic->GetName(), box_logic, true,
                        0, false);

  }



  G4ThreeVector XeBox::GenerateVertex(const G4String& region) const
  {
    return box_vertex_gen_->GenerateVertex(region);
  }


} // end namespace nexus
