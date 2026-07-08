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
    GeometryBase(), radius_(1.*mm), length_(1.*cm), fiber_dist_(0.*mm), al_(false), tefl_(false), isround_(true), core_mat_("EJ280"), sensortype_("PMT"), num_fibers_(1), lamp_size_(1.*cm), gap_(0.1*mm), rand_wls_(1), rand_att_(1), cyl_vertex_gen_(0)
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

        G4GenericMessenger::Command& al_cmd =
            msg_->DeclareProperty("aluminised",al_,"Are fibers aluminised or not");
        al_cmd.SetParameterName("al_cmd",false);

        G4GenericMessenger::Command& tefl_cmd =
            msg_->DeclareProperty("teflon",tefl_,"Teflon block behind fibers yes or no");
        tefl_cmd.SetParameterName("tefl_cmd",false);

        G4GenericMessenger::Command& shape_cmd =
            msg_->DeclareProperty("shape",isround_,"Shape of the fibers (round or square)");
        shape_cmd.SetParameterName("shape_cmd",false);

        G4GenericMessenger::Command& sensor_cmd =
            msg_->DeclareProperty("sensor",sensortype_,"Photosensor type (PMT or SiPM)");
        sensor_cmd.SetParameterName("sensor",false);

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

        G4GenericMessenger::Command& wls_cmd =
            msg_->DeclareProperty("rand_wls",rand_wls_,"Number of sigmas of the WLS uncertainty to shift the WLS abs. length by");
        wls_cmd.SetParameterName("rand_wls",false);

        G4GenericMessenger::Command& att_cmd =
            msg_->DeclareProperty("rand_att",rand_att_,"Number of sigmas of the attenuation uncertainty to shift the attenuation length by");
        att_cmd.SetParameterName("rand_att",false);

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

    // Uncertainty vectors
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
    std::vector<G4double> syst_WLSBCF92 = {
        0.31756971, 0.34330534, 0.49647789, 0.91688226,  1.92173165,
        0.32017544,  0.13471447,  0.08054783,  0.05475556,  0.04266491,
        0.03681777,  0.03377093,  0.03221736,  0.03149993,  0.0310701 ,
        0.03106284,  0.03152224,  0.03248003,  0.03392987,  0.03603926,
        0.0400336 ,  0.04683477,  0.06055102,  0.0867318 ,  0.15747915,
        0.55962362,  0.53056352,  0.958891  ,  1.60406058,   1.81228,
        2.04437714,  1.99898282,  1.7834041 };
    std::vector<G4double> syst_BCF92 = {
        0.07668161, 0.09358926, 0.10533346, 0.11013016, 0.11483066,
        0.12353031, 0.11796668, 0.1168047 , 0.11342559, 0.10970878,
        0.10140249, 0.1055834 , 0.11126771, 0.11538286, 0.11790927,
        0.11775898, 0.11775165, 0.10193637, 0.0816313 , 0.08958056,
        0.10466879, 0.11480426, 0.12031607, 0.12031458, 0.12058542,
        0.11967866, 0.1176351 , 0.11796293, 0.11801438, 0.11783206,
        0.11732908, 0.11689157, 0.11483542, 0.10616123, 0.10540304,
        0.10780882, 0.10568602, 0.10039032, 0.09108796, 0.07671777,
        0.05638106, 0.03654115, 0.0221245 , 0.01866518, 0.01386619};
    // LAB. This is just a volume of air surrounding the detector
    G4double xlab;
    G4double ylab;
    G4double puffer = 10*mm;
    if (issquare(num_fibers_)&&num_fibers_>1) {
        xlab = (sqrt(num_fibers_))*(2*radius_+fiber_dist_)-fiber_dist_;
        ylab = (sqrt(num_fibers_))*(2*radius_+fiber_dist_)-fiber_dist_;
    } 
    else {
        xlab = 25.4*mm;
        ylab = (num_fibers_)*(2*radius_+fiber_dist_)-fiber_dist_+puffer;
        if (ylab<xlab) ylab = xlab;
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
        core_mat->SetMaterialPropertiesTable(opticalprops::Y11(syst_WLSY11, syst_Y11, rand_wls_, rand_att_));
    }
    if (core_mat_=="BCF92") {
        core_mat = materials::PVT();
        core_mat->SetMaterialPropertiesTable(opticalprops::BCF92(0.395*mm, syst_WLSBCF92, syst_BCF92, rand_wls_, rand_att_));
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
            G4ThreeVector position = G4ThreeVector((xlab-2*radius_)/2,
                                    (ylab-2*radius_)/2+(int(num_fibers_)-1)*(fiber_dist_/2+radius_)-i*(2*radius_+fiber_dist_),0);
            new G4PVPlacement(0,position,fiber_logic,
                            fiber_logic->GetName(),lab_logic,true,cntr,true);
            cntr+=1;
        }
    }
            


    //Build the sensor
    if (sensortype_=="SiPM") {
        sensor_  = new GenericPhotosensor("SENSOR", 0.25 * mm, 0.25 * mm, thickness_);
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
        sensor_rot.rotateY(-pi/2);
        G4ThreeVector sensor_pos = G4ThreeVector((xlab+2*radius_)/2+0.1*mm, (ylab-2*radius_)/2, -5.*mm);   
        

        new G4PVPlacement(G4Transform3D(sensor_rot, sensor_pos), sensor_logic,
                            sensor_logic->GetName(), lab_logic, true,
                            cntr+1, true);
    }
    else if (sensortype_=="PMT") {

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
    }

    // Endcap volume to reflect trapped photons (Teflon/Aluminium/perfect absorber)
    if (al_) {
        G4Box* absorb_box = new G4Box("ABS",xlab/2,ylab/2,0.2*mm);
        G4Material* al_mat_ = materials::PolishedAl();
        al_mat_->SetMaterialPropertiesTable(opticalprops::PolishedAl());
        G4LogicalVolume* abs_log = new G4LogicalVolume(absorb_box,al_mat_,"ABS");
        new G4PVPlacement(0,G4ThreeVector((xlab-2*radius_)/2,(ylab-2.*radius_)/2,-length_/2-0.2*mm),abs_log,abs_log->GetName(),lab_logic,true,4,true);
    }

    // Reflective volume behind the fibers to increase efficiency (Teflon block)
    if(tefl_) {
        G4Box* absorb_box = new G4Box("ABS",0.2*mm,ylab/2,lamp_size_);
        G4Material* tefl_mat_ = materials::PVT();
        tefl_mat_->SetMaterialPropertiesTable(opticalprops::PTFE());
        G4LogicalVolume* abs_log = new G4LogicalVolume(absorb_box,tefl_mat_,"ABS");
        new G4PVPlacement(0,G4ThreeVector((xlab-4*radius_)/2,(ylab-2.*radius_)/2,-5*mm),abs_log,abs_log->GetName(),lab_logic,true,4,true);
    }
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
        // Calculate world size
        G4double xlab;
        G4double ylab;
        G4double puffer = 10*mm;
        bool k = issquare(num_fibers_);
        if (k && num_fibers_>1) {
            xlab = (sqrt(num_fibers_))*(2*radius_+fiber_dist_)-fiber_dist_;
            ylab = (sqrt(num_fibers_))*(2*radius_+fiber_dist_)-fiber_dist_;
        } 
        else {
            xlab = 25.4*mm;
            ylab = (num_fibers_)*(2*radius_+fiber_dist_)-fiber_dist_+puffer;
            if (ylab<xlab) ylab = xlab;
        }
        //BoxPointSampler* cyl_vertex_gen_ = new BoxPointSampler(0.1*mm,ylab,
        //                                lamp_size_,0,G4ThreeVector((xlab-6*radius_)/2-0.1*mm,(ylab-radius_-1.*mm)/2,-lamp_size_/2),0);
        //G4RotationMatrix *lamp_rot = new G4RotationMatrix();
        //lamp_rot->rotateY(pi/2);
        //CylinderPointSampler* cyl_vertex_gen_ = new CylinderPointSampler(0.,0.1*mm,lamp_size_,
        //                                0.,G4ThreeVector((xlab-6*radius_)/2-0.1*mm,(ylab-2*radius_)/2,-5.*mm),lamp_rot);
        //return cyl_vertex_gen_->GenerateVertex("WHOLE_VOL");
        // The center point of the lamp cylindrical volume
        G4ThreeVector origin_pos((xlab-6*radius_)/2-0.1*mm, (ylab-2*radius_)/2, -5.*mm);
        
        G4double max_radius = lamp_size_;
        G4double half_length = 0.1*mm; // very thin cylinder in z direction to approximate a disk-shaped lamp
        G4ThreeVector local_pos;

        // Generate points using a spherically symmetric probability distribution
        // Rejection sampling ensures it stays within the cylindrical volume dimensions
        while (true) {
            // Custom spherically symmetric parameter (3D Gaussian here)
            // Adjust 'sigma' for narrower or wider radial distributions
            G4double sigma = (8.64*12.7/7) * mm;  // see Systematics.ipynb fit parameter, scaled to mm in small angle approx.
            
            G4double x = G4RandGauss::shoot(0., sigma);
            G4double y = G4RandGauss::shoot(0., sigma);
            G4double z = 0.05*mm;
            
            // Constrain to the local cylinder bounds (prior to rotation)
            if ((x * x + y * y) <= (max_radius * max_radius) && std::abs(z) <= half_length) {
                local_pos = G4ThreeVector(x, y, z);
                break;
            }
        }

        // Apply rotation (equivalent to rotating Y by pi/2)
        local_pos.rotateY(pi / 2.0);

        return origin_pos + local_pos;
    }
}

bool OpticalFibre::issquare(G4int n) const {
    if (n==1){return false;}
    for (G4int i = 0; i < n / 2 + 2; i++) {
        if (i * i == n) {
      return true;
        }
    }
    return false;
}