// ----------------------------------------------------------------------------
// nexus | BlackBox.cc
//
// Sapphire disk and DICE board with configurable mask thickness parallelly
// placed in a black box.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "BlackBox.h"
#include "BlackBoxSiPMBoard.h"
#include "FactoryBase.h"

#include "MaterialsList.h"
#include "IonizationSD.h"
#include "OpticalMaterialProperties.h"
#include "Visibilities.h"

#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4GenericMessenger.hh>
#include <G4Orb.hh>
#include <G4NistManager.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4Material.hh>
#include <G4ThreeVector.hh>
#include <G4VisAttributes.hh>
#include <G4Colour.hh>
#include <G4SDManager.hh>
#include <G4UserLimits.hh>
#include <Randomize.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4RotationMatrix.hh>
#include <G4SystemOfUnits.hh>

namespace nexus {

  REGISTER_CLASS(BlackBox, GeometryBase)

  BlackBox::BlackBox():
    GeometryBase(),
    world_z_ (3. * m),
    world_xy_ (2. *m),
    box_z_ (2. * m),
    box_xy_ (1. *m),

    visibility_(0)
  {
    msg_ = new G4GenericMessenger(this, "/Geometry/BlackBox/",
				  "Control commands of BlackBox.");
    msg_->DeclareProperty("visibility", visibility_, "Black Box visibility");

    msg_->DeclarePropertyWithUnit("specific_vertex", "mm",  specific_vertex_,
      "Set generation vertex.");

    G4GenericMessenger::Command&  dice_board_x_pos_cmd =
      msg_->DeclareProperty("dice_board_x_pos", dice_board_x_pos_,
                            "X position of the dice");
    dice_board_x_pos_cmd.SetParameterName("dice_board_x_pos", true);
    dice_board_x_pos_cmd.SetUnitCategory("Length");

    G4GenericMessenger::Command&  dice_board_y_pos_cmd =
      msg_->DeclareProperty("dice_board_y_pos", dice_board_y_pos_,
                            "Y position of the dice");
    dice_board_y_pos_cmd.SetParameterName("dice_board_y_pos", true);
    dice_board_y_pos_cmd.SetUnitCategory("Length");

    G4GenericMessenger::Command&  dice_board_z_pos_cmd =
      msg_->DeclareProperty("dice_board_z_pos", dice_board_z_pos_,
                            "Z position of the dice");
    dice_board_z_pos_cmd.SetParameterName("dice_board_z_pos", true);
    dice_board_z_pos_cmd.SetUnitCategory("Length");

    G4GenericMessenger::Command&  rotation_cmd =
      msg_->DeclareProperty("rotation", rotation_,
                            "Angle of rotation of the DICE");
    rotation_cmd.SetParameterName("rotation", true);
    rotation_cmd.SetUnitCategory("Angle");

    G4GenericMessenger::Command&  mask_thickn_cmd =
      msg_->DeclareProperty("mask_thickn", mask_thickn_,
                            "Mask Thickness");
    mask_thickn_cmd.SetParameterName("mask_thickn", true);
    mask_thickn_cmd.SetUnitCategory("Length");

    G4GenericMessenger::Command&  membrane_thickn_cmd =
      msg_->DeclareProperty("membrane_thickn", membrane_thickn_,
                            "Membrane Thickness");
    membrane_thickn_cmd.SetParameterName("membrane_thickn", true);
    membrane_thickn_cmd.SetUnitCategory("Length");

    G4GenericMessenger::Command&  coating_thickn_cmd =
      msg_->DeclareProperty("coating_thickn", coating_thickn_,
                            "Coating Thickness");
    coating_thickn_cmd.SetParameterName("coating_thickn", true);
    coating_thickn_cmd.SetUnitCategory("Length");

    G4GenericMessenger::Command&  hole_diameter_cmd =
      msg_->DeclareProperty("hole_diameter", hole_diameter_,
                            "Hole Diameter");
    hole_diameter_cmd.SetParameterName("hole_diameter", true);
    hole_diameter_cmd.SetUnitCategory("Length");

    G4GenericMessenger::Command&  hole_x_cmd =
      msg_->DeclareProperty("hole_x", hole_x_,
                            "Hole X");
    hole_x_cmd.SetParameterName("hole_x", true);
    hole_x_cmd.SetUnitCategory("Length");

    G4GenericMessenger::Command&  hole_y_cmd =
      msg_->DeclareProperty("hole_y", hole_y_,
                            "Hole Y");
    hole_y_cmd.SetParameterName("hole_y", true);
    hole_y_cmd.SetUnitCategory("Length");

    msg_->DeclareProperty("membrane_hole", membrane_hole_,
                          "True if membrane has holes, false if not.");

    G4GenericMessenger::Command&  membrane_hole_diameter_cmd =
      msg_->DeclareProperty("membrane_hole_diameter", membrane_hole_diam_,
                            "Membrane Hole Diameter");
    membrane_hole_diameter_cmd.SetParameterName("membrane_hole_diameter", true);
    membrane_hole_diameter_cmd.SetUnitCategory("Length");

    G4GenericMessenger::Command&  membrane_hole_x_cmd =
      msg_->DeclareProperty("membrane_hole_x", membrane_hole_x_,
                            "Membrane Hole X");
    membrane_hole_x_cmd.SetParameterName("membrane_hole_x", true);
    membrane_hole_x_cmd.SetUnitCategory("Length");

    G4GenericMessenger::Command&  membrane_hole_y_cmd =
      msg_->DeclareProperty("membrane_hole_y", membrane_hole_y_,
                            "Membrane Hole Y");
    membrane_hole_y_cmd.SetParameterName("membrane_hole_y", true);
    membrane_hole_y_cmd.SetUnitCategory("Length");
  }

  BlackBox::~BlackBox()
  {
    delete msg_;
  }

  void BlackBox::Construct()
  {

  dice_ = new BlackBoxSiPMBoard();

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

  // BLACK BOX /////////////////////////////////////////////////

  G4String box_name = "BLACK BOX";

  G4Box* box_solid_vol =
    new G4Box(box_name, box_xy_/2., box_xy_/2., box_z_/2.);

  G4LogicalVolume* box_logic_vol =
    new G4LogicalVolume(box_solid_vol, world_mat, box_name);
  box_logic_vol->SetVisAttributes(G4VisAttributes::GetInvisible());

  G4VPhysicalVolume* box_phys_vol =
   new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),
                     box_logic_vol, box_name, world_logic_vol,
                     false, 0, false);

  // DICE BOARD //////////////////////////////////////////////

  dice_->SetMotherPhysicalVolume(box_phys_vol);
  dice_->SetMaskThickness    (mask_thickn_);
  dice_->SetMembraneThickness(membrane_thickn_);
  dice_->SetCoatingThickness (coating_thickn_);
  dice_->SetHoleDiameter     (hole_diameter_);
  dice_->SetHoleX            (hole_x_);
  dice_->SetHoleY            (hole_y_);
  dice_->SetHolesMembrane     (membrane_hole_);
  dice_->SetMembraneHoleDiameter     (membrane_hole_diam_);
  dice_->SetMembraneHoleX            (membrane_hole_x_);
  dice_->SetMembraneHoleY            (membrane_hole_y_);
  dice_->Construct();
  G4LogicalVolume* dice_board_logic = dice_->GetLogicalVolume();

  ////Dice Board placement
  G4ThreeVector post(dice_board_x_pos_,dice_board_y_pos_,dice_board_z_pos_);
  G4RotationMatrix* rot = new G4RotationMatrix();
  rot -> rotateY(rotation_);
  new G4PVPlacement(rot, post, dice_board_logic,
            	      dice_board_logic->GetName(), box_logic_vol, false, 1, false);

  // SAPPHIRE //////////////////////////////////////////////////

  G4String sapphire_name = "SAPPHIRE";

  G4Material* sapphire_mat = materials::Sapphire();

  sapphire_mat->SetMaterialPropertiesTable(opticalprops::Sapphire());

  //G4double pedot_thickn  = 1. *um;
  G4double tpb_thickn    = 1. *um;
  G4double sapphire_thickn  = 6. *mm;
  G4double pedot_thickn  = sapphire_thickn/2.;
  G4double window_thickn  = sapphire_thickn + pedot_thickn;
  G4double sapphire_diam = 20 *cm;

  G4Tubs* sapphire_solid_vol =
    new G4Tubs(sapphire_name, 0, sapphire_diam/2.,
                window_thickn/2., 0, twopi);

  G4LogicalVolume* sapphire_logic_vol =
    new G4LogicalVolume(sapphire_solid_vol, sapphire_mat, sapphire_name);

  G4VisAttributes sapphire_col = nexus::LightBlue();
  sapphire_logic_vol->SetVisAttributes(sapphire_col);

  G4ThreeVector sapphire_pos(dice_board_x_pos_,dice_board_y_pos_,-49*cm);

   new G4PVPlacement(rot, sapphire_pos,
                     sapphire_logic_vol, sapphire_name, box_logic_vol,
                     false, 0, false);

  // PEDOT ////////////////////////////////////////////////////////

  G4String pedot_name = "PEDOT";

  G4Material* pedot_mat = materials::TPB();

  sapphire_mat->SetMaterialPropertiesTable(opticalprops::Sapphire());

  G4Tubs* pedot_solid_vol =
    new G4Tubs(pedot_name, 0, sapphire_diam/2., pedot_thickn/2., 0, twopi);

  G4LogicalVolume* pedot_logic_vol =
    new G4LogicalVolume(pedot_solid_vol, pedot_mat, pedot_name);

    /// Add optical surface to PEDOT ///
    //G4OpticalSurface* pedot_coating_opsurf =
    //new G4OpticalSurface("SAPPHIRE_WNDW_PEDOT_OPSURF", glisur, ground,
      //                   dielectric_dielectric, .01);
    //new G4LogicalSkinSurface("SAPPHIRE_WNDW_PEDOT_OPSURF",
      //                       pedot_logic_vol, pedot_coating_opsurf);

  G4VisAttributes pedot_col = nexus::CopperBrown();
  pedot_logic_vol->SetVisAttributes(pedot_col);

  G4ThreeVector pedot_pos(0,0,window_thickn/2.- pedot_thickn/2.);

   new G4PVPlacement(0, pedot_pos,
                     pedot_logic_vol, pedot_name, sapphire_logic_vol,
                     false, 0, false);



  }

    G4ThreeVector BlackBox::GenerateVertex(const G4String& region) const
  {
    G4ThreeVector vertex(0.,0.,0.);

    // WORLD
    if (region == "WORLD") {
      vertex = G4ThreeVector(0.,0.,0.*mm);
    }
    else if (region == "AD_HOC") {
      // AD_HOC does not need to be shifted because it is passed by the user
      return specific_vertex_;
    }
    else {
      G4Exception("[BlackBox]", "GenerateVertex()", FatalException,
		  "Unknown vertex generation region!");
    }

    G4ThreeVector displacement = G4ThreeVector(0., 0., 0.);
    vertex = vertex + displacement;

    return vertex;
  }

} // end namespace nexus
