// ----------------------------------------------------------------------------
// nexus | PmtR7378A.cc
//
// Geometry of the Hamamatsu R7378A photomultiplier.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "PmtR7378A.h"

#include "SensorSD.h"
#include "OpticalMaterialProperties.h"
#include "MaterialsList.h"
#include "Visibilities.h"

#include <G4NistManager.hh>
#include <G4Tubs.hh>
#include <G4Sphere.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4Material.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4LogicalBorderSurface.hh>
#include <G4VisAttributes.hh>
#include <G4SDManager.hh>

#include <CLHEP/Units/SystemOfUnits.h>

namespace nexus {

  using namespace CLHEP;

  PmtR7378A::PmtR7378A(): GeometryBase()
  {
  }



  PmtR7378A::~PmtR7378A()
  {
  }



  void PmtR7378A::Construct()
  {
    // PMT BODY //////////////////////////////////////////////////////

    pmt_diam_   = 25.4 * mm;
    pmt_length_ = 43.0 * mm;

    G4Tubs* pmt_solid =
      new G4Tubs("PMT_R7378A", 0., pmt_diam_/2., pmt_length_/2., 0., twopi);

    G4Material* aluminum =
      G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");

    G4LogicalVolume* pmt_logic =
      new G4LogicalVolume(pmt_solid, aluminum, "PMT_R7378A");

    this->SetLogicalVolume(pmt_logic);


    // PMT WINDOW ////////////////////////////////////////////////////

    G4double window_diam = pmt_diam_;
    G4double window_length = 6. * mm;

    G4Tubs* window_solid =
      new G4Tubs("PMT_WINDOW", 0., window_diam/2., window_length/2., 0., twopi);

    G4Material* quartz =
      G4NistManager::Instance()->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
    quartz->SetMaterialPropertiesTable(opticalprops::FusedSilica());

    G4LogicalVolume* window_logic =
      new G4LogicalVolume(window_solid, quartz, "PMT_WINDOW");

    new G4PVPlacement(0, G4ThreeVector(0., 0., (pmt_length_-window_length)/2.),
		      window_logic, "PMT_WINDOW", pmt_logic, false, 0, false);

    G4VisAttributes wndw_col = nexus::Blue();
    window_logic->SetVisAttributes(wndw_col);


    // PHOTOCATHODE //////////////////////////////////////////////////

    G4double phcath_diam   = 22.0 * mm;
    G4double phcath_height =  4.0 * mm;
    G4double phcath_thickn =  0.1 * mm;
    G4double phcath_posz   =-16.0 * mm;

    G4double rmax =
      0.5 * (phcath_diam*phcath_diam / (4*phcath_height) + phcath_height);
    G4double rmin = rmax - phcath_thickn;
    G4double theta = asin(phcath_diam/(2.*rmax));

    G4Sphere* phcath_solid =
      new G4Sphere("PHOTOCATHODE", rmin, rmax, 0, twopi, 0, theta);

    G4LogicalVolume* phcath_logic =
      new G4LogicalVolume(phcath_solid, aluminum, "PHOTOCATHODE");

    G4VisAttributes vis_solid;
    vis_solid.SetForceSolid(true);
    phcath_logic->SetVisAttributes(vis_solid);

    //G4PVPlacement* phcath_physi =
      new G4PVPlacement(0, G4ThreeVector(0.,0.,phcath_posz), phcath_logic,
			"PHOTOCATHODE", window_logic, false, 0, false);

    // Sensitive detector
    SensorSD* pmtsd = new SensorSD("/PMT_R7378A/Pmt");
    pmtsd->SetDetectorVolumeDepth(2);
    pmtsd->SetTimeBinning(100.*nanosecond);
    G4SDManager::GetSDMpointer()->AddNewDetector(pmtsd);
    phcath_logic->SetSensitiveDetector(pmtsd);

    // OPTICAL SURFACES //////////////////////////////////////////////

    // The values for the efficiency are chosen in order to match
    // the curve of the quantum efficiency provided by Hamamatsu:
    // http://sales.hamamatsu.com/en/products/electron-tube-division/detectors/photomultiplier-tubes/part-r7378a.php
    // The source of light is point-like, isotropic and it has been placed at a
    // distance of 25 cm from the surface of the PMT window.
    // The quantity to be compared with the Hamamatsu curve is:
    // number of detected photons/ number of photons that reach the PMT window.
    // The total number of generated photons is taken as
    // number of photons that reach the PMT window because
    // light is generated only in that fraction of solid angle that subtends the
    // window of the PMT.

    const G4int entries = 61;

    G4double ENERGIES[entries] =
    { h_Planck * c_light / (710. * nm),  h_Planck * c_light / (700. * nm),
      h_Planck * c_light / (690. * nm),  h_Planck * c_light / (680. * nm),
      h_Planck * c_light / (670. * nm),  h_Planck * c_light / (660. * nm),
      h_Planck * c_light / (650. * nm),  h_Planck * c_light / (640. * nm),
      h_Planck * c_light / (630. * nm),  h_Planck * c_light / (620. * nm),
      h_Planck * c_light / (610. * nm),  h_Planck * c_light / (600. * nm),
      h_Planck * c_light / (590. * nm),  h_Planck * c_light / (580. * nm),
      h_Planck * c_light / (570. * nm),  h_Planck * c_light / (560. * nm),
      h_Planck * c_light / (550. * nm),  h_Planck * c_light / (540. * nm),
      h_Planck * c_light / (530. * nm),  h_Planck * c_light / (520. * nm),
      h_Planck * c_light / (510. * nm),  h_Planck * c_light / (500. * nm),
      h_Planck * c_light / (490. * nm),  h_Planck * c_light / (490. * nm),
      h_Planck * c_light / (470. * nm),  h_Planck * c_light / (460. * nm),
      h_Planck * c_light / (450. * nm),  h_Planck * c_light / (440. * nm),
      h_Planck * c_light / (430. * nm),  h_Planck * c_light / (420. * nm),
      h_Planck * c_light / (410. * nm),  h_Planck * c_light / (400. * nm),
      h_Planck * c_light / (390. * nm),  h_Planck * c_light / (380. * nm),
      h_Planck * c_light / (370. * nm),  h_Planck * c_light / (360. * nm),
      h_Planck * c_light / (350. * nm),  h_Planck * c_light / (340. * nm),
      h_Planck * c_light / (330. * nm),  h_Planck * c_light / (330. * nm),
      h_Planck * c_light / (320. * nm),  h_Planck * c_light / (310. * nm),
      h_Planck * c_light / (300. * nm),  h_Planck * c_light / (290. * nm),
      h_Planck * c_light / (280. * nm),  h_Planck * c_light / (270. * nm),
      h_Planck * c_light / (260. * nm),  h_Planck * c_light / (250. * nm),
      h_Planck * c_light / (240. * nm),  h_Planck * c_light / (230. * nm),
      h_Planck * c_light / (220. * nm),  h_Planck * c_light / (210. * nm),
      h_Planck * c_light / (200. * nm),  h_Planck * c_light / (190. * nm),
      h_Planck * c_light / (188. * nm),  h_Planck * c_light / (182. * nm),
      h_Planck * c_light / (175. * nm),  h_Planck * c_light / (170. * nm),
      h_Planck * c_light / (165. * nm),  h_Planck * c_light / (161. * nm),
      h_Planck * c_light / (150. * nm)};

    /*G4double EFFICIENCY[entries] =
      { 0.0000, 0.0001, 0.0002, 0.0004, 0.0007, 0.0014, 0.0025, 0.0040, 0.0067, 0.0103,
        0.02036148857142857,
        0.023049426020408166,
        0.032076,
        0.041441807250000004,
        0.05257076611111111,
        0.06271214166007906,
        0.07413015176470589,
        0.08247699493917275,
        0.103258155,
        0.13918327500000002,
        0.17217333547677263,
        0.19051794074534162,
        0.205,
        0.21861783530638854,
        0.23402662215477996,
        0.25951520745967743,
        0.27307700723995265,
        0.29805508839779005,
        0.3136967046027678,
        0.3173974841201717,
        0.3332063071143375,
        0.3288931391666667,
        0.3457820935793358,
        0.3520075733226195,
        0.3380934076937407,
        0.32539800473812003,
        0.323261621089109,
        0.3313663558186398,
        0.33501,
        0.33859124922053235,
        0.3257554706939502,
        0.31526222645290586,
        0.30845409437352245,
        0.2839027647936508,
        0.2645969619493671,
        0.24519174026185625,
        0.23669509320066337,
        0.23184893928821176,
        0.22118048119220837,
        0.21390258960509553,
        0.2136632483714445,
        0.21618593770586444,
        0.23137359594506215,
        0.2794295775648949,
        0.2803400049844237,
        0.31889201589595373,
        0.1750, 0.1310, 0.0470,
        0.0000};*/

        G4double EFFICIENCY[entries] = {1,1,1,1,1,1,1,1,1,1,
                                        1,1,1,1,1,1,1,1,1,1,
                                        1,1,1,1,1,1,1,1,1,1,
                                        1,1,1,1,1,1,1,1,1,1,
                                        1,1,1,1,1,1,1,1,1,1,
                                        1,1,1,1,1,1,1,1,1,1,1};

/* RAW EFFICIENCIES (from Hamamatsu, PMT 1666):
    { 0.0000, 0.0001, 0.0002, 0.0004, 0.0007, 0.0014, 0.0025, 0.0040, 0.0067, 0.0103,
    	0.0149, 0.0205, 0.0270, 0.0342, 0.0419, 0.0502, 0.0594, 0.0713, 0.0897, 0.1155,
      0.1393, 0.1534, 0.1635, 0.1744, 0.1880, 0.2041, 0.2191, 0.2310, 0.2397, 0.2466,
      0.2517, 0.2534, 0.2556, 0.2550, 0.2525, 0.2503, 0.2513, 0.2492, 0.2472, 0.2441,
      0.2393, 0.2310, 0.2192, 0.2057, 0.1923, 0.1774, 0.1708, 0.1707, 0.1662, 0.1634,
      0.1669, 0.1777, 0.1879, 0.2055, 0.2020, 0.2050, 0.1890, 0.1750, 0.1310, 0.0470,
      0.0000};*/

    G4double REFLECTIVITY[entries] =
      { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
	      0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
	      0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
	      0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
	      0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
	      0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0.};

    G4MaterialPropertiesTable* phcath_mpt = new G4MaterialPropertiesTable();
    phcath_mpt->AddProperty("EFFICIENCY", ENERGIES, EFFICIENCY, entries);
    phcath_mpt->AddProperty("REFLECTIVITY", ENERGIES, REFLECTIVITY, entries);

    G4OpticalSurface* phcath_opsur =
      new G4OpticalSurface("PHOTOCATHODE", unified, polished, dielectric_metal);
    phcath_opsur->SetMaterialPropertiesTable(phcath_mpt);

    new G4LogicalSkinSurface("PHOTOCATHODE", phcath_logic, phcath_opsur);
  }


} // end namespace nexus
