#include "FibreTest.h"
#include "GenericWLSFiber.h"
#include "MaterialsList.h"
#include "OpticalMaterialProperties.h"
#include "Visibilities.h"
#include "CylinderPointSampler.h"
#include "FactoryBase.h"

#include <G4Tubs.hh>
#include <G4Box.hh>
#include <G4Colour.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4Material.hh>
#include <G4GenericMessenger.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4NistManager.hh>

#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>
#include <Randomize.hh>
#include <string>

using namespace nexus;
using namespace CLHEP;

REGISTER_CLASS(FibreTest,GeometryBase)

FibreTest::FibreTest():GeometryBase(), radius_(1.*mm), length_(1.*cm), cyl_vertex_gen_(0) {
    msg_=new G4GenericMessenger(this,"/Geometry/FibreTest/","Control commands of geometry OpticalFibre");

    G4GenericMessenger::Command& radius_cmd = 
            msg_->DeclareProperty("radius",radius_,"Radius of the cylindrical optical fibre");
    radius_cmd.SetUnitCategory("Length");
    radius_cmd.SetParameterName("radius",false);
    radius_cmd.SetRange("radius>0.");

    G4GenericMessenger::Command& length_cmd = 
            msg_->DeclareProperty("length",length_,"Length of the cylindrical optical fibre");
    length_cmd.SetUnitCategory("Length");
    length_cmd.SetParameterName("length",false);
    length_cmd.SetRange("length>0.");


    cyl_vertex_gen_ = new CylinderPointSampler(radius_, 0.5*length_, 0.,  0., G4ThreeVector(0., 0., 0.), 0);
}
FibreTest::~FibreTest() {
    delete cyl_vertex_gen_;
    delete msg_;
}
void FibreTest::Construct(){
    std::vector<G4double> syst_WLSY11 = {
        0.29954168, 0.27814827, 0.24220817, 0.23282155, 0.20920007,
        0.17464791, 0.13866487, 0.09812573, 0.05833202, 0.03885444,
        0.0329607 , 0.03328822, 0.03278581, 0.03170133, 0.0314333 ,
        0.03183566, 0.03327871, 0.03279741, 0.03196641, 0.03212028,
        0.03311408, 0.03479884, 0.03856599, 0.04695476, 0.06639167,
        0.12054313, 0.2819549 , 0.5832129 , 1.10675271, 1.08053982,
        1.56033944, 1.89307361, 1.72522633};
    std::vector<G4double> syst_Y11 = {
       0.03897649, 0.05880983, 0.07513664, 0.09312284, 0.10294135,
       0.10330598, 0.10221196, 0.09634171, 0.08374408, 0.0662849 ,
       0.05914838, 0.06013325, 0.06854581, 0.07668044, 0.08113814,
       0.08112786, 0.07526827, 0.04917664, 0.02943515, 0.03517056,
       0.05099695, 0.0640956 , 0.06888065, 0.06928472, 0.06627835,
       0.06126615, 0.056651  , 0.05352115, 0.05140073, 0.04819723,
       0.04489344, 0.04094288, 0.03628014, 0.03055602, 0.0282274 ,
       0.02671656, 0.0245285 , 0.02159721, 0.01852729, 0.01483524,
       0.01117362, 0.0109576 , 0.00861888, 0.01479128, 0.06207768};
    G4Box* lab_solid = new G4Box("LAB", 2 * mm,2 * mm,1.1*cm);

    G4Material* air=G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    air->SetMaterialPropertiesTable(opticalprops::Vacuum());
    G4LogicalVolume* lab_logic =
      new G4LogicalVolume(lab_solid,
                          air,
                          "LAB");
    lab_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    this->SetLogicalVolume(lab_logic);
    G4Material* ps = materials::PS();
    G4Material* tpb = materials::TPB();
    GenericWLSFiber* fiber_ = new GenericWLSFiber("Y11", true, radius_, length_, true, true, tpb, ps, true);
    fiber_->SetCoreOpticalProperties(opticalprops::Y11(syst_WLSY11, syst_Y11, 1, 1));
    fiber_->SetCoatingOpticalProperties(opticalprops::TPB());
    fiber_->Construct();
    G4LogicalVolume* fiber_logic = fiber_->GetLogicalVolume();
    new G4PVPlacement(0,G4ThreeVector(0,0,0),fiber_logic,
                            fiber_logic->GetName(),lab_logic,true,0,true);
}

G4ThreeVector FibreTest::GenerateVertex(const G4String& region) const {
    return cyl_vertex_gen_->GenerateVertex(region);
}