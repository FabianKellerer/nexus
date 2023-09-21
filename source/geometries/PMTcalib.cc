#include "PMTcalib.h"

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

REGISTER_CLASS(PMTcalib,GeometryBase)

PMTcalib::PMTcalib():
    GeometryBase()
    {
    }

PMTcalib::~PMTcalib()
{
}

void PMTcalib::Construct()
{   

    // LAB. This is just a volume of air surrounding the detector
    G4double xlab = 3. *cm;
    G4double ylab = 3. *cm;
    G4double zlab = 30.*cm;

    G4Box* lab_solid = new G4Box("LAB", xlab,ylab,zlab);

    G4Material* air=G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    air->SetMaterialPropertiesTable(opticalprops::Vacuum());
    G4LogicalVolume* lab_logic =
      new G4LogicalVolume(lab_solid,
                          air,
                          "LAB");
    lab_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    this->SetLogicalVolume(lab_logic);

    G4ThreeVector sensor_pos = G4ThreeVector(0,0,-25.*cm-21.5*mm);   
    PmtR7378A pmt;
    pmt.Construct();
    pmt_logic_ = pmt.GetLogicalVolume();
    new G4PVPlacement(0,sensor_pos,
  		      pmt_logic_, "PMT", lab_logic, true, 0, true);
}

G4ThreeVector PMTcalib::GenerateVertex() const
{

    G4ThreeVector vertex(0., 0., -30.*cm);

    return vertex;
}