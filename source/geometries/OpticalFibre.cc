// ----------------------------------------------------------------------------
// nexus | OpticalFibre.cc
//
// Bundle of cylindrical optical fibres with single photosensor.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "OpticalFibre.h"

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

REGISTER_CLASS(OpticalFibre,GeometryBase)

OpticalFibre::OpticalFibre():
    GeometryBase(), radius_(1.*mm), length_(1.*cm), fiber_dist_(0.*mm), isround_(true), core_mat_("EJ280"), num_fibers_(1), lamp_size_(1.*cm),gap_(0.1*mm), cyl_vertex_gen_(0)
    {
        msg_=new G4GenericMessenger(this,"/Geometry/OpticalFibre/","Control commands of geometry OpticalFibre");

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

        G4GenericMessenger::Command& dist_cmd = 
            msg_->DeclareProperty("fiber_dist",fiber_dist_,"Distance between the fibers in the bundle");
        length_cmd.SetUnitCategory("Length");
        length_cmd.SetParameterName("fiber_dist",false);
        length_cmd.SetRange("fiber_dist>0.");

        G4GenericMessenger::Command& shape_cmd =
            msg_->DeclareProperty("shape",isround_,"Shape of the fibers (round or square)");
        shape_cmd.SetParameterName("shape_cmd",false);

        G4GenericMessenger::Command& mat_cmd = 
            msg_->DeclareProperty("core_mat",core_mat_,"Core material (EJ280, EJ286 or Y11)");
        mat_cmd.SetParameterName("core_mat",false);

        G4GenericMessenger::Command& num_cmd =
            msg_->DeclareProperty("num_fibers",num_fibers_,"Number of fibers");
        num_cmd.SetParameterName("num_fibers",false);

        G4GenericMessenger::Command& lamp_cmd =
            msg_->DeclareProperty("lamp_size",lamp_size_,"Size of the lamp");
        lamp_cmd.SetUnitCategory("Length");
        lamp_cmd.SetParameterName("lamp_size",false);
        lamp_cmd.SetRange("lamp_size>0");

        G4GenericMessenger::Command& gap_cmd =
            msg_->DeclareProperty("gap_size",gap_,"Size of the gap between fiber and sensor");
        gap_cmd.SetUnitCategory("Length");
        gap_cmd.SetParameterName("gap_size",false);
        gap_cmd.SetRange("gap_size>0");

        cyl_vertex_gen_ = new CylinderPointSampler(radius_, length_, 0.,  0., G4ThreeVector(0., 0., 0.), 0);

        // hardcoded thickness of the photosensor
        thickness_=2.*mm;
    }

OpticalFibre::~OpticalFibre()
{
    delete cyl_vertex_gen_;
    delete msg_;
}

void OpticalFibre::Construct()
{   

    // LAB. This is just a volume of air surrounding the detector
    G4double xlab;
    G4double ylab;
    G4double puffer = 10*mm;
    if (issquare(num_fibers_)&&num_fibers_>1) {
        xlab = (sqrt(num_fibers_))*(2*radius_+fiber_dist_)-fiber_dist_;
        ylab = (sqrt(num_fibers_))*(2*radius_+fiber_dist_)-fiber_dist_;
    } 
    else {
        xlab = 25.4*mm;//2*radius_;
        ylab = (num_fibers_)*(2*radius_+fiber_dist_)-fiber_dist_+puffer;
    }

G4Box* lab_solid = new G4Box("LAB", xlab,ylab,length_+gap_+1.*cm);

    G4Material* air=G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    air->SetMaterialPropertiesTable(opticalprops::Vacuum());
    G4LogicalVolume* lab_logic =
      new G4LogicalVolume(lab_solid,
                          air,
                          "LAB");
    lab_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
    this->SetLogicalVolume(lab_logic);
    G4String name = "OPTICAL_FIBRE";

    //define materials
    G4Material* FP = materials::FPethylene();
    FP->SetMaterialPropertiesTable(opticalprops::FPethylene());

    G4Material* pmma = materials::PMMA();
    pmma->SetMaterialPropertiesTable(opticalprops::PMMA());

    G4Material* core_mat;
    if (core_mat_=="EJ280") {
        core_mat = materials::PVT();
        core_mat->SetMaterialPropertiesTable(opticalprops::EJ280());
    }
    if (core_mat_=="EJ286") {
        core_mat = materials::PVT();
        core_mat->SetMaterialPropertiesTable(opticalprops::EJ286());
    }
    if (core_mat_=="Y11") {
        core_mat = materials::Y11();
        core_mat->SetMaterialPropertiesTable(opticalprops::Y11());
    }
    if (core_mat_=="BCF92") {
        core_mat = materials::PVT();
        core_mat->SetMaterialPropertiesTable(opticalprops::BCF92(1.3*mm));
    }

    G4Material* tpb = materials::TPB();
    G4bool  coating = false;

    //define logical volume
    GenericWLSFiber* fiber =
    new GenericWLSFiber("FIBER", isround_, 2*radius_,
                        length_, true, coating, tpb,
                        core_mat, true);
    fiber->Construct();
    G4LogicalVolume* fiber_logic = fiber->GetLogicalVolume();


    G4int cntr=0;
    if (issquare(num_fibers_)) {

        //place fibers in a square
        for(G4int i=0; i<sqrt(num_fibers_); i++){
            for(G4int j=0; j<sqrt(num_fibers_); j++){
                G4ThreeVector position = G4ThreeVector((int(sqrt(num_fibers_))-1)*(2*radius_+fiber_dist_)-i*(2*radius_+fiber_dist_),
                                                       (int(sqrt(num_fibers_))-1)*(2*radius_+fiber_dist_)-j*(2*radius_+fiber_dist_),
                                                        0);
                new G4PVPlacement(0,position,fiber_logic,
                                fiber_logic->GetName(),lab_logic,true,cntr,true);
                cntr+=1;
            }
        }
    }
    else {
        //place fibers in a line
        for(G4int i=0; i<num_fibers_; i++){
            G4ThreeVector position = G4ThreeVector((xlab-2*radius_)/2,(int(num_fibers_)-1)*(2*radius_+fiber_dist_)-i*(2*radius_+fiber_dist_)+puffer/2,0);
            new G4PVPlacement(0,position,fiber_logic,
                            fiber_logic->GetName(),lab_logic,true,cntr,true);
            cntr+=1;
        }
    }
            


    //Build the sensor
/*    sensor_  = new GenericPhotosensor("SENSOR", xlab, ylab, thickness_);
    sensor_ -> SetVisibility(true);

    //Set the sensor window material
    G4Material* window_mat_ =
      G4NistManager::Instance()->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
    window_mat_->SetMaterialPropertiesTable(opticalprops::FusedSilica());
    G4MaterialPropertyVector* window_rindex = window_mat_->GetMaterialPropertiesTable()->GetProperty("RINDEX");
    //G4MaterialPropertyVector* window_rindex = air->GetMaterialPropertiesTable()->GetProperty("RINDEX");
    sensor_ -> SetWindowRefractiveIndex(window_rindex);

    //Set the optical properties of the sensor
    G4MaterialPropertiesTable* photosensor_mpt = new G4MaterialPropertiesTable();
    G4double energy[]       = {0.2 * eV, 3.5 * eV, 3.6 * eV, 11.5 * eV};
    G4double reflectivity[] = {0.0     , 0.0     , 0.0     ,  0.0     };
    G4double efficiency[]   = {1.0     , 1.0     , 1.0     ,  1.0     };
    photosensor_mpt->AddProperty("REFLECTIVITY", energy, reflectivity, 4);
    photosensor_mpt->AddProperty("EFFICIENCY",   energy, efficiency,   4);
    sensor_->SetOpticalProperties(photosensor_mpt);
    //sensor_->SetTimeBinning(1*us);

    sensor_->SetWithWLSCoating(false);

    //Set sensor depth and naming order
    sensor_ ->SetSensorDepth(1);
    //sensor_ ->SetMotherDepth(0);
    //sensor_ ->SetNamingOrder(0);

    sensor_ -> Construct();

    //Placing the sensor
    G4LogicalVolume* sensor_logic = sensor_ -> GetLogicalVolume();
    G4RotationMatrix sensor_rot;
    sensor_rot.rotateY(pi);
    G4ThreeVector sensor_pos = G4ThreeVector((xlab-2*radius_)/2,
                                (ylab-2.*radius_)/2,
                                length_/2+thickness_/2+gap_/2);   
    

    new G4PVPlacement(G4Transform3D(sensor_rot, sensor_pos), sensor_logic,
                        sensor_logic->GetName(), lab_logic, true,
                        cntr+1, true);
*/

    G4RotationMatrix sensor_rot;
    sensor_rot.rotateY(pi);
    G4ThreeVector sensor_pos = G4ThreeVector((xlab-2*radius_)/2,
                                (ylab-2.*radius_)/2,
                                length_/2+gap_/2+21.5*mm);   
    PmtR7378A pmt;
    pmt.Construct();
    pmt_logic_ = pmt.GetLogicalVolume();
    new G4PVPlacement(G4Transform3D(sensor_rot, sensor_pos),
  		      pmt_logic_, "PMT", lab_logic, true, cntr+1, true);

    // Debug volume to reflect trapped photons (Teflon)
    G4Box* absorb_box = new G4Box("ABS",xlab/2,ylab/2,0.2*mm);
    G4Material* tefl_mat_ = materials::PolishedAl();
    tefl_mat_->SetMaterialPropertiesTable(opticalprops::PolishedAl());
    G4LogicalVolume* abs_log = new G4LogicalVolume(absorb_box,tefl_mat_,"ABS");
    //new G4PVPlacement(0,G4ThreeVector((xlab-2*radius_)/2,(ylab-2.*radius_)/2,-length_/2-0.2*mm),abs_log,abs_log->GetName(),lab_logic,true,4,true);
    // Reflective surface
    //G4MaterialPropertiesTable* refl_surf = new G4MaterialPropertiesTable();
    //G4double energy2[]       = {0.2 * eV, 3.5 * eV, 3.6 * eV, 11.5 * eV};
    //G4double reflectivity2[] = {0.93     , 0.93     , 0.93     ,  0.93     };
    //refl_surf->AddProperty("REFLECTIVITY", energy2, reflectivity2, 4);
    //G4OpticalSurface* refl_opsurf =
    //new G4OpticalSurface("Refl_optSurf", unified, ground, dielectric_dielectric);
    //refl_opsurf->SetMaterialPropertiesTable(refl_surf);
    //new G4LogicalSkinSurface(name + "_optSurf", abs_log, refl_opsurf);

}

G4ThreeVector OpticalFibre::GenerateVertex(const G4String& region) const
{

    if (region!="LAMP") {
        
        if (isround_) {
            G4int i = G4RandFlat::shootInt((long) 0, sqrt(num_fibers_));
            G4int j = G4RandFlat::shootInt((long) 0, sqrt(num_fibers_));
            G4ThreeVector position = G4ThreeVector((int(sqrt(num_fibers_))-1)*(2*radius_+fiber_dist_)-i*(2*radius_+fiber_dist_),
                                                   (int(sqrt(num_fibers_))-1)*(2*radius_+fiber_dist_)-j*(2*radius_+fiber_dist_),
                                                    0);
            CylinderPointSampler* cyl_vertex_gen_ = new CylinderPointSampler(0., 2*length_, radius_,  0., position, 0);
            return cyl_vertex_gen_->GenerateVertex(region);
        } else {
            G4int i = G4RandFlat::shootInt((long) 0, sqrt(num_fibers_));
            G4int j = G4RandFlat::shootInt((long) 0, sqrt(num_fibers_));
            G4ThreeVector position = G4ThreeVector((int(sqrt(num_fibers_))-1)*(2*radius_+fiber_dist_)-i*(2*radius_+fiber_dist_),
                                                   (int(sqrt(num_fibers_))-1)*(2*radius_+fiber_dist_)-j*(2*radius_+fiber_dist_),
                                                    0);
            BoxPointSampler* cyl_vertex_gen_ = new BoxPointSampler(0., 2*length_, radius_,  0., position, 0);
            return cyl_vertex_gen_->GenerateVertex(region);
        }
    }
    else {
        G4double xlab = 2*radius_+1.*mm;
        G4double ylab;
        bool k = false;
        for (G4int i = 0; i < num_fibers_ / 2 + 2; i++) {
            if (i * i == num_fibers_) {
                k = true;
            }
        }
        if (k) {
            ylab = (sqrt(num_fibers_))*(2*radius_+fiber_dist_);
        } else {
            ylab = (num_fibers_)*(2*radius_+fiber_dist_)-fiber_dist_;
        }
        //BoxPointSampler* cyl_vertex_gen_ = new BoxPointSampler(0.1*mm,ylab,
        //                                lamp_size_,0,G4ThreeVector((xlab-6*radius_)/2-0.1*mm,(ylab-radius_-1.*mm)/2,-lamp_size_/2),0);
        G4RotationMatrix *lamp_rot = new G4RotationMatrix();
        lamp_rot->rotateY(pi/2);
        CylinderPointSampler* cyl_vertex_gen_ = new CylinderPointSampler(0.,0.1*mm,25.4/2.*mm,
                                        0.,G4ThreeVector((xlab-6*radius_)/2-0.1*mm,(int(num_fibers_)-1)*(2*radius_+fiber_dist_)-2*(2*radius_+fiber_dist_)+5.*mm,-5.*mm),lamp_rot);
        return cyl_vertex_gen_->GenerateVertex("WHOLE_VOL");
    }
}

bool OpticalFibre::issquare(G4int n) {
    for (G4int i = 0; i < n / 2 + 2; i++) {
        if (i * i == n) {
      return true;
        }
    }
    return false;
}