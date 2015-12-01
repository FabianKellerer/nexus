// ---------------------------------------------------------------------------
//  $Id$
// 
//  Author:  <miquel.nebot@ific.uv.es>
//  Created: 18 Sept 2013
//
//  Copyright (c) 2013 NEXT Collaboration. All rights reserved.
// ---------------------------------------------------------------------------

#include "PetPlainDice.h"

#include "SiPMpetVUV.h"
#include "PmtSD.h"
#include "MaterialsList.h"
#include "OpticalMaterialProperties.h"

#include <G4Box.hh>
#include <G4VisAttributes.hh>
#include <G4NistManager.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4SDManager.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4GenericMessenger.hh>

#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>

#include <sstream>

namespace nexus {

  using namespace CLHEP;

  PetPlainDice::PetPlainDice(): 
    BaseGeometry(),
    rows_(8),
    columns_(8),
    visibility_ (1),
    ysize_(5.*cm),
    xsize_(5.*cm)
  {
    /// Messenger
    msg_ = new G4GenericMessenger(this, "/Geometry/PetalX/", "Control commands of geometry Pet.");
    msg_->DeclareProperty("plain_dice_vis", visibility_, "Kapton Dice Boards Visibility");
    msg_->DeclareProperty("plain_columns", columns_, "Number of rows in SiPMs");
    msg_->DeclareProperty("plain_rows", rows_, "Number of rows in SiPMs");
  }

  PetPlainDice::~PetPlainDice()
  {
  }

  void PetPlainDice::SetSize(G4double xsize, G4double ysize)
  {
  xsize_ = xsize;
  ysize_ = ysize;
  }

  void PetPlainDice::Construct()
  {
   

    //   const G4double sipm_pitch = 6.2 * mm;
     G4double sipm_pitch = ysize_/rows_;
    const G4double coating_thickness = 0.1 * micrometer;
    const G4double board_thickness = 0.3 * mm;
    //const G4double board_side_reduction = .5 * mm;
    const G4double board_side_reduction = 0. * mm;    
    // const G4double db_x = columns_ * sipm_pitch - 2. * board_side_reduction ;  
    // const G4double db_y =    rows_ * sipm_pitch - 2. * board_side_reduction ;
     const G4double db_x = xsize_ ;  
    const G4double db_y =    rows_ * sipm_pitch - 2. * board_side_reduction ;
    const G4double db_z = board_thickness;

    // Outer element volume  /////////////////////////////////////////////////// 
    const G4double border = .5*mm;
    const G4double out_x = db_x ;
    const G4double out_y = db_y;
    const G4double out_z = db_z + 2.*border;

    dimensions_.setX(out_x);
    dimensions_.setY(out_y);
    dimensions_.setZ(out_z);

    G4Material* out_material = G4NistManager::Instance()->FindOrBuildMaterial("G4_lXe");

    G4Box* out_solid = new G4Box("LXE_DICE", out_x/2., out_y/2., out_z/2.);
    G4LogicalVolume* out_logic = 
      new G4LogicalVolume(out_solid, out_material,  "LXE_DICE");
    this->SetLogicalVolume(out_logic);


    // KAPTON BOARD /////////////////////////////////////////////////

    G4Box* board_solid = new G4Box("DICE_BOARD", db_x/2., db_y/2., db_z/2.);
 
    // G4Material* kapton =
    //   G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON");
    G4Material* teflon = 
      G4NistManager::Instance()->FindOrBuildMaterial("G4_TEFLON");
 
    G4LogicalVolume* board_logic = 
      new G4LogicalVolume(board_solid, teflon, "DICE_BOARD");
    new G4PVPlacement(0, G4ThreeVector(0.,0., -border), board_logic,
			"DICE_BOARD", out_logic, false, 0, true);

    // OPTICAL SURFACE FOR REFLECTION
    G4OpticalSurface* db_opsur = new G4OpticalSurface("DICE_BOARD");
    db_opsur->SetType(dielectric_metal);
    db_opsur->SetModel(unified);
    db_opsur->SetFinish(ground);
    db_opsur->SetSigmaAlpha(0.1);
   
    // db_opsur->SetMaterialPropertiesTable(OpticalMaterialProperties::PTFE_with_TPB());
    db_opsur->SetMaterialPropertiesTable(OpticalMaterialProperties::PTFE_LXe());
    
    new G4LogicalSkinSurface("DICE_BOARD", board_logic, db_opsur);

   
    // WLS COATING //////////////////////////////////////////////////
    
    G4Box* coating_solid = 
      new G4Box("DB_WLS_COATING", db_x/2., db_y/2., coating_thickness/2.);

    G4Material* TPB = MaterialsList::TPB();
    TPB->SetMaterialPropertiesTable(OpticalMaterialProperties::TPB_LXe());

    G4LogicalVolume* coating_logic =
      new G4LogicalVolume(coating_solid, TPB, "DB_WLS_COATING");

    G4double pos_z = db_z/2. - coating_thickness / 2.;
    new G4PVPlacement(0, G4ThreeVector(0., 0., pos_z), coating_logic,
    		      "DB_WLS_COATING", board_logic, false, 0, true);

    
   

    // SETTING VISIBILITIES   //////////
    // _visibility  = true;
    if (visibility_) {
      /*
      G4VisAttributes silicon_col(G4Colour(1., 1., 0.));
      silicon_col.SetForceSolid(true);
      sipm_logic->SetVisAttributes(silicon_col);
      */
      G4VisAttributes tpb_col(G4Colour(1., 1., 1.));
      tpb_col.SetForceSolid(true);
      //   coating_logic->SetVisAttributes(tpb_col);
      
    }
    else {
      //   coating_logic->SetVisAttributes(G4VisAttributes::Invisible);
      //     sipm_logic->SetVisAttributes(G4VisAttributes::Invisible);
    }
  }

  G4ThreeVector PetPlainDice::GetDimensions() const
  {
    return dimensions_;
  }
  
  const std::vector<std::pair<int, G4ThreeVector> >& PetPlainDice::GetPositions()
  {
    return positions_;
  }

  // void PetPlainDice::SetMaterial(G4Material& mat)
  // {
  //   _out_mat = mat;
  // }


} // end namespace nexus
