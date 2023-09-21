// ----------------------------------------------------------------------------
// nexus | FiberEfficiency.cc
//
// Box containing optical fibers
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "FiberEfficiency.h"

#include "FactoryBase.h"
#include "OpticalMaterialProperties.h"
#include "Visibilities.h"

#include <G4GenericMessenger.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4NistManager.hh>
#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4VisAttributes.hh>

using namespace nexus;

REGISTER_CLASS(FiberEfficiency, GeometryBase)

namespace nexus {

  FiberEfficiency::FiberEfficiency():
    GeometryBase(),
    fiber_radius_(0.5 * mm),
    length_ (25 * cm),
    coating_ ("TPB"),
    fiber_type_ ("Y11"),
    coated_(true)
  {
    msg_ = new G4GenericMessenger(this, "/Geometry/FiberEfficiency/",
      "Control commands of geometry FiberEfficiency.");

    G4GenericMessenger::Command&  length_cmd =
      msg_->DeclareProperty("length", length_,
                            "Barrel fiber length");
    length_cmd.SetUnitCategory("Length");

    G4GenericMessenger::Command&  fiber_radius_cmd =
      msg_->DeclareProperty("fiber_radius", fiber_radius_,
                            "Fiber radius");
    fiber_radius_cmd.SetUnitCategory("Length");

    msg_->DeclareProperty("coating", coating_, "Fiber coating (TPB or PTH)");
    msg_->DeclareProperty("fiber_type", fiber_type_, "Fiber type (Y11 or B2)");
    msg_->DeclareProperty("coated", coated_, "Coat fibers with WLS coating");

  }



  FiberEfficiency::~FiberEfficiency()
  {
    delete msg_;
  }



  void FiberEfficiency::Construct()
  {

    G4cout << "[FiberEfficiency] *** Fiber efficiency geometry ***" << G4endl;
    G4cout << "[FiberEfficiency] Using " << fiber_type_ << " fibers";
    if (coated_)
      G4cout << " with " << coating_ << " coating";
    G4cout << G4endl;

    world_z_ = 50 * cm;
    world_xy_ = 50 * cm;

    inside_cylinder_ = new CylinderPointSampler2020(0, fiber_radius_, length_/2, 0, 2 * M_PI);

    G4Material *this_fiber = nullptr;
    G4MaterialPropertiesTable *this_fiber_optical = nullptr;
    if (fiber_type_ == "Y11") {
      this_fiber = materials::Y11();
      this_fiber_optical = opticalprops::Y11();
    } else {
      G4Exception("[FiberEfficiency]", "Construct()",
                  FatalException, "Invalid fiber type, must be Y11 or B2");
    }

    G4Material *this_coating = nullptr;
    G4MaterialPropertiesTable *this_coating_optical = nullptr;
    if (coated_) {
      if (coating_ == "TPB") {
        this_coating = materials::TPB();
        this_coating_optical = opticalprops::TPB();
      } else if (coating_ == "TPH") {
        this_coating = materials::TPH();
        this_coating_optical = opticalprops::TPH();
      } else {
        G4Exception("[FiberEfficiency]", "Construct()",
                    FatalException, "Invalid coating, must be TPB or TPH");
      }
    }

    fiber_ = new GenericWLSFiber(fiber_type_, true, fiber_radius_, length_, true, coated_, this_coating, this_fiber, true);

    // WORLD /////////////////////////////////////////////////

    G4String world_name = "WORLD";

    G4Material* world_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");

    world_mat->SetMaterialPropertiesTable(opticalprops::Vacuum());

    G4Box* world_solid_vol =
     new G4Box(world_name, world_xy_/2., world_xy_/2., world_z_/2.);

    G4LogicalVolume* world_logic_vol =
      new G4LogicalVolume(world_solid_vol, world_mat, world_name);
    world_logic_vol->SetVisAttributes(G4VisAttributes::GetInvisible());
    GeometryBase::SetLogicalVolume(world_logic_vol);

    // FIBER ////////////////////////////////////////////////////

    fiber_->SetCoreOpticalProperties(this_fiber_optical);
    fiber_->SetCoatingOpticalProperties(this_coating_optical);

    fiber_->Construct();
    G4LogicalVolume* fiber_logic = fiber_->GetLogicalVolume();
    if (fiber_type_ == "Y11")
      fiber_logic->SetVisAttributes(nexus::LightGreenAlpha());
    else if (fiber_type_ == "B2")
      fiber_logic->SetVisAttributes(nexus::LightBlueAlpha());

    // TEFLON PANEL //////////////////////////////////////////
    G4Box* teflon_panel =
     new G4Box("teflon", 3.6, 0.5/2, 3.6);

    G4Material* teflon =
      G4NistManager::Instance()->FindOrBuildMaterial("G4_TEFLON");

    G4LogicalVolume* teflon_logic_vol =
      new G4LogicalVolume(teflon_panel, teflon, "TEFLON");

    G4OpticalSurface* opsur_teflon =
      new G4OpticalSurface("TEFLON_OPSURF", unified, polished, dielectric_metal);
    opsur_teflon->SetMaterialPropertiesTable(opticalprops::PTFE());

    new G4LogicalSkinSurface("TEFLON_OPSURF", teflon_logic_vol, opsur_teflon);

    // Aluminzed endcap /////////////////////////////////////////////

    G4Material* fiber_end_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");

    G4double fiber_end_z = 0.1 * mm;

    G4Tubs* fiber_end_solid_vol =
      new G4Tubs("fiber_end", 0, fiber_radius_/2, fiber_end_z, 0, 2 * M_PI);

    G4LogicalVolume* fiber_end_logic_vol =
      new G4LogicalVolume(fiber_end_solid_vol, fiber_end_mat, "FIBER_END");
    G4OpticalSurface* opsur_al =
      new G4OpticalSurface("POLISHED_AL_OPSURF", unified, polished, dielectric_metal);
    opsur_al->SetMaterialPropertiesTable(opticalprops::PolishedAl());

    new G4LogicalSkinSurface("POLISHED_AL_OPSURF", fiber_end_logic_vol, opsur_al);

    // PMT /////////////////////////////////////////////

    G4double photosensor_thickness = 5*cm;
    pmt_  = new GenericPhotosensor("FIBER_SENSOR", 1 * mm, 1 * mm, photosensor_thickness);
    pmt_->SetVisibility(true);
    pmt_->SetWindowRefractiveIndex(opticalprops::FusedSilica()->GetProperty("RINDEX"));
    pmt_->SetSensorDepth(1);
    G4MaterialPropertiesTable* photosensor_mpt = new G4MaterialPropertiesTable();

    G4int h_n_entries = 29;
    G4double h_energy[] = {
      h_Planck * c_light / (943.02 * nm), h_Planck * c_light / (925.42 * nm),
      h_Planck * c_light / (911.73 * nm), h_Planck * c_light / (900.00 * nm),
      h_Planck * c_light / (890.22 * nm), h_Planck * c_light / (882.40 * nm),
      h_Planck * c_light / (874.58 * nm), h_Planck * c_light / (862.85 * nm),
      h_Planck * c_light / (851.12 * nm), h_Planck * c_light / (829.61 * nm),
      h_Planck * c_light / (794.41 * nm), h_Planck * c_light / (737.71 * nm),
      h_Planck * c_light / (665.36 * nm), h_Planck * c_light / (612.57 * nm),
      h_Planck * c_light / (571.51 * nm), h_Planck * c_light / (534.36 * nm),
      h_Planck * c_light / (505.03 * nm), h_Planck * c_light / (479.61 * nm),
      h_Planck * c_light / (458.10 * nm), h_Planck * c_light / (432.68 * nm),
      h_Planck * c_light / (411.17 * nm), h_Planck * c_light / (364.25 * nm),
      h_Planck * c_light / (323.18 * nm), h_Planck * c_light / (305.59 * nm),
      h_Planck * c_light / (299.72 * nm), h_Planck * c_light / (293.85 * nm),
      h_Planck * c_light / (289.94 * nm), h_Planck * c_light / (284.08 * nm),
      h_Planck * c_light / (274.30 * nm),
    };

    G4double h_efficiency[] = {
      0.02, 0.05,
      0.12, 0.33,
      0.92, 2.03,
      4.71, 11.14,
      20.79, 38.80,
      47.09, 54.74,
      60.96, 65.03,
      70.87, 75.60,
      75.60, 70.87,
      62.29, 49.16,
      37.97, 23.65,
      12.95, 6.10,
      3.20, 1.61,
      0.72, 0.38,
      0.23
    };

    G4int sipm_entries = 24;
    G4double sipm_energy[] = {
      h_Planck * c_light / (866.20 * nm), h_Planck * c_light / (808.45 * nm),
      h_Planck * c_light / (766.20 * nm), h_Planck * c_light / (721.13 * nm),
      h_Planck * c_light / (685.92 * nm), h_Planck * c_light / (647.89 * nm),
      h_Planck * c_light / (623.94 * nm), h_Planck * c_light / (597.18 * nm),
      h_Planck * c_light / (573.24 * nm), h_Planck * c_light / (545.07 * nm),
      h_Planck * c_light / (518.31 * nm), h_Planck * c_light / (502.82 * nm),
      h_Planck * c_light / (454.93 * nm), h_Planck * c_light / (421.13 * nm),
      h_Planck * c_light / (395.77 * nm), h_Planck * c_light / (378.87 * nm),
      h_Planck * c_light / (367.61 * nm), h_Planck * c_light / (359.15 * nm),
      h_Planck * c_light / (349.30 * nm), h_Planck * c_light / (340.85 * nm),
      h_Planck * c_light / (336.62 * nm), h_Planck * c_light / (332.39 * nm),
      h_Planck * c_light / (326.76 * nm), h_Planck * c_light / (319.72 * nm)
    };

    G4double sipm_efficiency[] = {
      5.75,
      9.38,
      12.48,
      16.37,
      20.18,
      24.34,
      28.67,
      33.01,
      37.26,
      41.50,
      45.93,
      48.32,
      50.00,
      47.61,
      43.54,
      39.03,
      34.25,
      29.91,
      25.40,
      20.62,
      16.02,
      11.50,
      6.81,
      3.36
    };

    for (G4int i=0; i < h_n_entries; i++) {
      h_efficiency[i] = 124 * h_efficiency[i] / (( 1 / h_energy[i] / nm) * (h_Planck * c_light)) / 100;
    }

    for (G4int i=0; i < sipm_entries; i++) {
      sipm_efficiency[i] /= 100;
    }

    G4double energy[]       = {opticalprops::optPhotMinE_, opticalprops::optPhotMaxE_};
    G4double reflectivity[] = {0.0     , 0.0     };
    G4double efficiency[]   = {1.      , 1.      };
    photosensor_mpt->AddProperty("REFLECTIVITY", energy, reflectivity, 2);
    photosensor_mpt->AddProperty("EFFICIENCY",   h_energy, h_efficiency, h_n_entries);
    pmt_->SetOpticalProperties(photosensor_mpt);
    pmt_->Construct();
    G4LogicalVolume* pmt_logic = pmt_->GetLogicalVolume();
    new G4PVPlacement(0, G4ThreeVector(0.,0.,-(length_/2.+photosensor_thickness/2.+0.5*mm)),
		      pmt_logic, "PMT",
		      world_logic_vol, false, 0, true);

    // PLACEMENT /////////////////////////////////////////////
    G4double x = 0;
    G4double y = 0;
    new G4PVPlacement(0, G4ThreeVector(x,y),
                      fiber_logic, "FIBER", world_logic_vol,
                      false, 0, false);
    new G4PVPlacement(0, G4ThreeVector(x,y,length_/2 + fiber_end_z),
                      fiber_end_logic_vol, "FIBER_END", world_logic_vol,
                      false, 0, false);
                      // length_/2*0.7
    new G4PVPlacement(0, G4ThreeVector(x,y-0.5,length_ / 2 * 0.8),
                      teflon_logic_vol, "TEFLON", world_logic_vol,
                      false, 0, false);

  }



  G4ThreeVector FiberEfficiency::GenerateVertex(const G4String& region) const
  {
    G4ThreeVector vertex(0., 150. * mm, length_ / 2 * 0.8);

    // WORLD
    if (region == "INSIDE_FIBER") {
      return inside_cylinder_->GenerateVertex("VOLUME");

    }
    else if (region == "OUTSIDE_FIBER") {
      return vertex;
    }
    else {
      G4Exception("[FiberEfficiency]", "GenerateVertex()", FatalException,
        "Unknown vertex generation region!");
    }
    return vertex;
  }


} // end namespace nexus