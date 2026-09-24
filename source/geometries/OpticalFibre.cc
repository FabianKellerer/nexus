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
    GeometryBase(), radius_(1.*mm), length_(1.*cm), fiber_dist_(0.*mm), doubleclad_(true), al_(true), tefl_(false), spec_(false), isround_(true), core_mat_("EJ280"), sensortype_("PMT"), num_fibers_(1), lamp_size_(1.*cm), gap_(0.1*mm), rand_wls_(0), rand_att_(0), rand_sigma_(0), qe_(1), cyl_vertex_gen_(0)
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

        G4GenericMessenger::Command& spec_cmd =
            msg_->DeclareProperty("spectrometer",spec_,"Spectrometer behind fibers yes or no");
        spec_cmd.SetParameterName("spec_cmd",false);

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

        G4GenericMessenger::Command& sigma_cmd =
            msg_->DeclareProperty("rand_sigma",rand_sigma_,"Number of sigmas of the attenuation uncertainty to shift the attenuation length by");
        att_cmd.SetParameterName("rand_sigma",false);

        G4GenericMessenger::Command& qe_cmd =
            msg_->DeclareProperty("qe",qe_,"Quantum efficiency of the WLS");
        qe_cmd.SetParameterName("qe",false);

        G4GenericMessenger::Command& doubleclad_cmd =
            msg_->DeclareProperty("doubleclad",doubleclad_,"Double cladding or not");
        doubleclad_cmd.SetParameterName("double_clad",false);

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

    // Uncertainty vectors - Total
    std::vector<G4double> syst_WLSY11 = {2.0326059884672274, 
        1.9350266098720144, 1.844026589692775, 1.6977332960106888, 1.583671272564732, 1.5532213155328924, 
        1.2583903815759254, 0.7041010022308332, 0.2358227664808773, 0.08640247001168902, 0.05105627821895477, 
        0.052305012819547844, 0.06640723168334015, 0.07309266491398161, 0.05683510523121253, 0.04320956705914462, 
        0.05085866341968649, 0.06275230376545193, 0.09026151263296603, 0.09507205236156549, 0.11069969210107729, 
        0.13243752202895698, 0.19100960435603911, 0.23891103480995216, 0.29441243325452043, 0.30851916087412656, 
        0.40690977443953824, 0.5036874005140725, 0.46368448491001896, 0.3434358948133334};
    std::vector<G4double> syst_Y11 = {
        3.951850977559591, 
        6.523479217306917, 10.138275603792696, 14.22228924405363, 15.152039831734859, 13.97549707972901, 
        13.668649694642248, 12.145729413748258, 8.566377027282508, 7.623004535459913, 8.15758232548256, 
        10.848687197244342, 13.492160121790265, 16.12292190375665, 15.640088034030192, 13.737027754220081, 
        6.090717834964776, 2.13391696084334, 3.323121621783207, 6.948685986636572, 10.944914048829565, 
        13.261700183578158, 13.282529017184137, 13.099687814276546, 11.524488707485332, 10.14755834699325, 
        9.460607620097704, 9.002647971894174, 8.380006223287438, 7.52551731741322, 6.493902480995727, 
        5.098187638250474, 3.33218377618802, 2.758404010038841, 2.498739459720643, 1.9262697340774957, 
        1.213225591856089, 0.584665774620175, 0.18916233832095639, 0.05234985769054529, 0.022605286447232305,
        0.007464580175701787};
    std::vector<G4double> syst_WLSBCF92 = {1.557111137489016, 
        1.5275946439949935, 1.5100087418647234, 1.4316654373831685, 1.275867602335821, 1.0569717481420842, 
        0.7570795508140168, 0.4463270481650585, 0.25310730366827117, 0.13509815003416398, 0.07916098844846176, 
        0.053521761272318905, 0.041591891206679686, 0.035153926939535715, 0.030640639205775075, 0.027430323595788345, 
        0.024976148228522187, 0.02484833666183279, 0.025901235829404113, 0.029152950202253603, 0.03382354439989134, 
        0.04015294719732863, 0.05553194236514641, 0.08790765810070811, 0.09365081313885587, 0.12928637486786582, 
        0.3591056597577124};
    std::vector<G4double> syst_BCF92 = {
        31.863572757063935, 
        7.729369091977371, 1.4941940384715897, 2.1447210007699153, 5.945655571874951, 11.687694819571842, 
        14.604070313323685, 15.879811117879047, 14.079037179604573, 11.333897360216987, 9.246267682013915, 
        8.486454106431106, 8.201598546988535, 7.706345766410162, 7.054567474015823, 5.900852455445136, 
        4.815131404386551, 2.921835510599543, 2.5792013659104325, 2.779389028640027, 2.517181534443222, 
        1.820991423686093, 1.1086276547131013, 0.5512087383141567, 0.20932099498554863, 0.06274993681972137, 
        0.015364191758086398};
    std::vector<G4double> syst_WLSBCF92_2mm = {22.322064152807084, 
        3.1738433982709613, 1.5151039731349998, 0.9721540251214319, 0.804155121943178, 0.7902766719458124, 
        0.9069793935524967, 0.957337409529103, 0.838215318254885, 0.664525820046821, 0.5341424538310738, 
        0.3657741726547332, 0.2489770741715708, 0.1755154509749875, 0.14573947822702618, 0.13538736024294376, 
        0.15043575804164114, 0.4371360790020529, 7.466986930294742, 76.22489101435404};
    std::vector<G4double> syst_BCF92_2mm = {
        13.475633187827086, 3.0845067530882084, 0.6635737644568879, 0.1372851938232902, 0.02889022616609853
    };
    std::vector<G4double> syst_WLSBCF91 = {1.4804189670184869, 
        1.4245922173493049, 1.412199109797844, 1.2704069670475233, 1.141128992554561, 1.1788611272543892, 
        1.0013955004247612, 0.5813148433023889, 0.24761175421817522, 0.08838634671053519, 0.05568382484614725, 
        0.06336617285504313, 0.0752942857793936, 0.07653810982346956, 0.06412151929066548, 0.05217835896429217, 
        0.05231237823897066, 0.0713134189043534, 0.08267493071032511, 0.09875372854223037, 0.1032615523069587, 
        0.12535552159688595, 0.16481701767864287, 0.21938749805950375, 0.23451803740505084, 0.3279313874781274, 
        0.413849691854618, 0.5697998634412138, 0.5136141240563457, 0.31872845972197444};
    std::vector<G4double> syst_BCF91 = {22.049162978403363, 
        16.445936458845246, 9.523807719011671, 3.3824464097002007, 2.1020408066088527, 1.9516197180490398, 
        1.9524089559025435, 1.4018378087717605, 0.65291910017797, 0.19297860324917632, 0.03952396470633514, 
        0.008702577117923212, 0.0029427092803867563};

    //uncertainty vectors - manufacturing differences
    std::vector<G4double> syst_BCF91_M = {7.442732374080068, 
        6.351464822268071, 4.7686257927590505, 1.9463678847464594, 1.1321188953280765, 1.0910804600070183, 
        1.0563017259739813, 0.744917471386619, 0.3256208103105999, 0.0903225677653937, 0.016316278197959087, 
        0.0024978535955453864, 0.0};
    std::vector<G4double> syst_Y11_M = {2.0202068143622927, 
        3.7334765210753686, 5.60332730097983, 7.7185108054463365, 9.714304038318435, 9.558478519229922, 
        9.241328798908409, 8.93229300246718, 7.421759217393541, 5.219459540145046, 4.814535727019093, 
        5.334057455293114, 6.42019135430776, 7.721229760348157, 8.46819056326845, 8.007287613568712, 
        6.856329618737617, 2.8903868910253045, 1.368680732878047, 2.5339104656206093, 4.826214641945471, 
        6.892656148532167, 7.747421494067968, 7.92248184738345, 7.56446872492184, 6.709973476524898, 
        6.058949773492239, 5.795655949719391, 5.636500122311823, 5.41207193209774, 5.034524515182438, 
        4.355031271061526, 3.4342198984001975, 2.30093779313343, 2.064547596863941, 1.8329679247800283, 
        1.401642959714907, 0.8502551160969012, 0.3793001366260554, 0.11337251513907083, 0.025014542108357615, 
        0.00891885788925397};
    std::vector<G4double> syst_WLSBCF92_M = {0.8796310230472657, 
        0.8549825170662467, 0.861317307396165, 0.7973600808507065, 0.721071073733736, 0.6016111117224587, 
        0.4417465717375406, 0.260656513106415, 0.14854533818087237, 0.08160087463601452, 0.0478291356335527, 
        0.03295158306142624, 0.025496077787160232, 0.021985154357095422, 0.01942328236009304, 0.017498180643784743, 
        0.016107746592005007, 0.016057530047128975, 0.017000427132760815, 0.01921376701638304, 0.023206200661527863, 
        0.026308141347914267, 0.03666499402912999, 0.058140597012252444, 0.05900475227704639, 0.09434755876027955, 
        0.22035787419162275};
    std::vector<G4double> syst_WLSBCF92_2mm_M = {21.274600324142046, 
        3.03611034628292, 1.4474720393363152, 0.9293613297805439, 0.7677300417497889, 0.752688148857639, 
        0.8632207561503051, 0.9093710786847438, 0.7963782127550056, 0.6309367882840236, 0.5072220025252183, 
        0.34719789929087314, 0.23644252782683872, 0.1663158426298309, 0.1377055227564194, 0.1275957421588165, 
        0.13649531575213733, 0.37930955623262, 6.241878076368693, 59.04948782942525};
    /*std::vector<G4double> syst_WLSY11_M = {
        3.4449607689059287, 
        2.0845485163075423, 0.6778577728811916, 0.16641590840746764, 0.043946284234532194, 0.02425885133130992, 
        0.023941045107760347, 0.029185955333199343, 0.0334670725960524, 0.02760599627287378, 0.020463652324858122, 
        0.02127323572925406, 0.02825587907503537, 0.040764958247356015, 0.0498174199688172, 0.054857023805037174, 
        0.07458520996065733, 0.11488546878393072, 0.16227881297081112, 0.2064917974989712, 0.24143658768585094};*/
    std::vector<G4double> syst_WLSY11_M = { 
        0.9568432108254918, 0.905292086467, 0.857262903408, 0.8031960026253675, 0.7273967709515743, 
        0.6831279696779592, 0.4159287465936, 0.1500484193752647, 0.05994408820112083, 0.04574944088201120, 
        0.034307253395308386, 0.036326958739466, 0.04988484741575766, 0.0381765349996225, 0.029310337210823195, 
        0.03646800268972645, 0.042106151249108624, 0.06264114209806373, 0.06484520609139349, 0.07837493901099597, 
        0.08960107935159256, 0.11124206862, 0.1544935038252093, 0.19730573059124693, 0.2121481068234, 
        0.22353214016591683, 0.30364636840448433, 0.31579967611721066, 0.16384515974529076};
    std::vector<G4double> syst_WLSBCF91_M = {1.1441603263529758, 1.159702647210629, 1.0478785092548384, 0.9051250843875618, 0.9589573042642915, 
        0.8288613242407584, 0.49106593782593944, 0.22077412246626912, 0.07986282745231416, 0.0502248005163302, 
        0.057560833750762964, 0.06864589344015042, 0.06978657388325514, 0.058123671054423606, 0.04799687168539973, 
        0.047690625957815924, 0.06503519881425975, 0.0740201925779011, 0.08985705348906284, 0.09295966867574001, 
        0.1120690554472324, 0.14512307733245858, 0.19357109226981778, 0.20181796027435137, 0.2828306020976765, 
        0.35095722625135556, 0.5120258348180301, 0.4481730570715898, 0.26022386307175654};

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

    G4Material* teflon = G4NistManager::Instance()->FindOrBuildMaterial("G4_TEFLON");

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
        core_mat->SetMaterialPropertiesTable(opticalprops::Y11(syst_WLSY11, syst_Y11, rand_wls_, rand_att_, qe_));
    }
    if (core_mat_=="BCF92") {
        core_mat = materials::PVT();
        core_mat->SetMaterialPropertiesTable(opticalprops::BCF92(syst_WLSBCF92, syst_BCF92, rand_wls_, rand_att_, qe_));
    }
    if (core_mat_=="BCF92_2mm") {
        core_mat = materials::PVT();
        core_mat->SetMaterialPropertiesTable(opticalprops::BCF92_2mm(syst_WLSBCF92_2mm, syst_BCF92_2mm, rand_wls_, rand_att_, qe_));
    }
    if (core_mat_=="BCF91") {
        core_mat = materials::PVT();
        core_mat->SetMaterialPropertiesTable(opticalprops::BCF91(syst_WLSBCF91, syst_BCF91, rand_wls_, rand_att_, qe_));
    }

    G4Material* tpb = materials::TPB();
    G4bool  coating = false;

    //define logical volume
    G4double clad_percent_;
    if (core_mat_=="Y11") {
        clad_percent_ = 0.02;
    }
    else {
        clad_percent_ = 0.03;
    }
    GenericWLSFiber* fiber =
    new GenericWLSFiber("FIBER", isround_, 2*radius_,
                        length_, clad_percent_, doubleclad_, coating, tpb,
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
                            cntr, true);
        cntr+=1;
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
                pmt_logic_, "PMT", lab_logic, true, cntr, true);
        cntr+=1;
    }

    // Endcap volume to reflect trapped photons (Teflon/Aluminium/perfect absorber)
    if (al_) {
        G4Box* absorb_box = new G4Box("ABS",xlab/2,ylab/2,0.2*mm);
        G4Material* al_mat_ = materials::PolishedAl();
        al_mat_->SetMaterialPropertiesTable(opticalprops::PerfectAbsorber());
        G4LogicalVolume* abs_log = new G4LogicalVolume(absorb_box,al_mat_,"ABS");
        new G4PVPlacement(0,G4ThreeVector((xlab-2*radius_)/2,(ylab-2.*radius_)/2,-length_/2-0.2*mm),abs_log,abs_log->GetName(),lab_logic,true,cntr,true);
        cntr+=1;
    }

    // Reflective volume behind the fibers to increase efficiency (Teflon block)
    if(tefl_) {
        G4Box* tefl_box = new G4Box("TEFL",0.2*mm,ylab/2,lamp_size_);
        G4LogicalVolume* tefl_log = new G4LogicalVolume(tefl_box,teflon,"TEFL");
        G4OpticalSurface* opsur_teflon = new G4OpticalSurface("TEFLON_OPSURF", unified, ground, dielectric_metal);
        opsur_teflon->SetMaterialPropertiesTable(opticalprops::PTFE());
        new G4LogicalSkinSurface("TEFLON_OPSURF", tefl_log, opsur_teflon);
        new G4PVPlacement(0,G4ThreeVector((xlab)/2+0.21*mm,(ylab-2.*radius_)/2,-5*mm),tefl_log,tefl_log->GetName(),lab_logic,true,cntr,true);
        cntr+=1;
    }
    // Spectrometer volume to measure the WLS absorption of the fibers
    if(spec_) {
        G4Tubs* spectrometer = new G4Tubs("SPEC",0,0.1*mm,0.01*mm,0,2*pi);
        G4LogicalVolume* spec_log = new G4LogicalVolume(spectrometer,teflon,"SPEC");
        G4RotationMatrix spec_rot;
        spec_rot.rotateY(pi/2);
        G4ThreeVector spec_pos = G4ThreeVector((xlab)/2+0.02*mm,(ylab-2.*radius_)/2,-5*mm);
        new G4PVPlacement(G4Transform3D(spec_rot,spec_pos),spec_log,spec_log->GetName(),lab_logic,true,cntr,true);
        cntr+=1;
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
            G4double sigma = (1+rand_sigma_)*(8.47*12.7/7) * mm;  // see Systematics.ipynb fit parameter, scaled to mm in small angle approx.
            
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