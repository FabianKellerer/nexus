#include "FiberAbsorption.h"

#include "CylinderPointSampler.h"
#include "BoxPointSampler.h"
#include "MaterialsList.h"
#include "IonizationSD.h"
#include "SensorSD.h"
#include "OpticalMaterialProperties.h"
#include "FactoryBase.h"
#include "Visibilities.h"
#include "GenericWLSFiber.h"
#include "PmtR7378A.h"

#include <G4GenericMessenger.hh>
#include <G4Tubs.hh>
#include <G4Box.hh>
#include <G4NistManager.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4Material.hh>
#include <GenericPhotosensor.h>
#include <G4VisAttributes.hh>
#include <G4SDManager.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>

#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>
#include <Randomize.hh>
#include <string>

using namespace nexus;
using namespace CLHEP;

REGISTER_CLASS(FiberAbsorption,GeometryBase)

FiberAbsorption::FiberAbsorption():
    GeometryBase(), lambdawls_(1.3*mm), cyl_vertex_gen_(0)
    {
        msg_=new G4GenericMessenger(this,"/Geometry/FiberAbsorption/","Control commands of geometry FiberAbsorption");

        G4GenericMessenger::Command& lambdawls_cmd = 
            msg_->DeclareProperty("lambdawls",lambdawls_,"Minimal WLS absorption length");
        lambdawls_cmd.SetUnitCategory("Length");
        lambdawls_cmd.SetParameterName("lambdawls",false);
        lambdawls_cmd.SetRange("lambdawls>0.");
    }

FiberAbsorption::~FiberAbsorption()
{
    delete cyl_vertex_gen_;
}

void FiberAbsorption::Construct()
{
    G4Box* lab_solid = new G4Box("LAB", 46.*mm,25.4*mm,11.*cm);

    G4Material* air=G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    air->SetMaterialPropertiesTable(opticalprops::Vacuum());
    G4LogicalVolume* lab_logic =
      new G4LogicalVolume(lab_solid,
                          air,
                          "LAB");
    lab_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    this->SetLogicalVolume(lab_logic);
    G4String name = "FIBER_ABSORPTION";

    G4Material* core_mat = materials::PVT();
    core_mat->SetMaterialPropertiesTable(opticalprops::BCF92(lambdawls_));
    GenericWLSFiber* fiber =
    new GenericWLSFiber("FIBER", true, 2.*mm,
                        11.*cm, true, false, materials::TPB(),
                        core_mat, true);
    fiber->Construct();
    G4LogicalVolume* fiber_logic = fiber->GetLogicalVolume();
    new G4PVPlacement(0,G4ThreeVector(23.*mm,-2.*mm,0),fiber_logic,fiber_logic->GetName(),lab_logic,true,0,true);
    new G4PVPlacement(0,G4ThreeVector(23.*mm, 0.*mm,0),fiber_logic,fiber_logic->GetName(),lab_logic,true,1,true);
    new G4PVPlacement(0,G4ThreeVector(23.*mm, 2.*mm,0),fiber_logic,fiber_logic->GetName(),lab_logic,true,2,true);

    G4RotationMatrix sensor_rot;
    sensor_rot.rotateY(pi/2);
    G4ThreeVector sensor_pos = G4ThreeVector(0,0,0);   
    PmtR7378A pmt;
    pmt.Construct();
    pmt_logic_ = pmt.GetLogicalVolume();
    new G4PVPlacement(G4Transform3D(sensor_rot, sensor_pos),
  		      pmt_logic_, "PMT", lab_logic, true, 3, true);

    G4Box* filter = new G4Box("FILT",0.1*mm,24.5/2*mm,24.5/2*mm);
    G4MaterialPropertiesTable* filter_optprop = new G4MaterialPropertiesTable();
    G4double energy2[]       = {0.2 * eV, 2.694 * eV, 2.695 * eV, 11.5 * eV};
    G4double reflectivity2[] = {1.0     , 1.0     , 0.0     ,  0.0     };
    G4double transmission[]  = {0.0     , 0.0     , 1.0     ,  1.0     };
    G4double refraction[]    = {1.1     , 1.1     , 1.1     ,  1.1     };
    G4double absorption[]    = {0.1*nm  , 0.1*nm  , 10*m    ,  10*m    };
    //filter_optprop->AddProperty("REFLECTIVITY", energy2, reflectivity2, 4);
    //filter_optprop->AddProperty("TRANSMITTANCE", energy2, transmission, 4);
    filter_optprop->AddProperty("RINDEX"       , energy2, refraction,   4);
    filter_optprop->AddProperty("ABSLENGTH"    , energy2, absorption,   4);
    G4Material* filt = G4NistManager::Instance()->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
    filt->SetMaterialPropertiesTable(filter_optprop);
    G4LogicalVolume* filter_log = new G4LogicalVolume(filter,filt,"FILT");
    new G4PVPlacement(0,G4ThreeVector(21.9*mm,0,0),filter_log,filter_log->GetName(),lab_logic,true,4,true);
    // Reflective surface
    //G4MaterialPropertiesTable* refl_surf = new G4MaterialPropertiesTable();
    //refl_surf->AddProperty("REFLECTIVITY", energy2, reflectivity2, 4);
    //refl_surf->AddProperty("TRANSMITTANCE", energy2, transmission, 4);
    //G4OpticalSurface* refl_opsurf =
    //new G4OpticalSurface("Refl_optSurf", unified, polished, dielectric_metal);
    //refl_opsurf->SetMaterialPropertiesTable(refl_surf);
    //new G4LogicalSkinSurface("Refl_optSurf", filter_log, refl_opsurf);
}


G4ThreeVector FiberAbsorption::GenerateVertex(const G4String& region) const
{
    G4ThreeVector position = G4ThreeVector(24.5*mm,0,0);
    BoxPointSampler* cyl_vertex_gen_ = new BoxPointSampler(0.1*mm, 6.*mm, 6.*mm, 0, position, 0);
    return cyl_vertex_gen_->GenerateVertex("WHOLE_VOL");
}
    