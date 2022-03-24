//------------------------------------------------------------------
// Nexus | Alpha particle generator
//
// This class produces alpha particles of a specific energy range
//
//------------------------------------------------------------------

#include "AlphaGenerator.h"

#include "DetectorConstruction.h"
#include "GeometryBase.h"
#include "FactoryBase.h"

#include <G4Event.hh>
#include <G4GenericMessenger.hh>
#include <G4RunManager.hh>
#include <G4ParticleTable.hh>
#include <G4RandomDirection.hh>
#include <Randomize.hh>

#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Units/PhysicalConstants.h"

using namespace nexus;
REGISTER_CLASS(AlphaGenerator, G4VPrimaryGenerator)

namespace nexus {
    AlphaGenerator::AlphaGenerator() : geom_(0), energy_low_(1*MeV), energy_high_(2*MeV) {

        msg_ = new G4GenericMessenger(this,"/Generator/AlphaGenerator/","Controls commands of the alpha particle generator.");

        msg_->DeclareProperty("region",region_,"Generation region");

        msg_->DeclareProperty("energy_low",energy_low_,"Lowest possible alpha particle energy");

        msg_->DeclareProperty("energy_high",energy_high_,"Highest possible alpha particle energy");

        particle_defalpha = G4ParticleTable::GetParticleTable()->FindParticle("alpha");

        DetectorConstruction* detconst = (DetectorConstruction*)
            G4RunManager::GetRunManager()->GetUserDetectorConstruction();
        geom_ = detconst->GetGeometry();
    }

    AlphaGenerator::~AlphaGenerator()
    {

    }

    void AlphaGenerator::GeneratePrimaryVertex(G4Event* evt) {

        //ask geometry to generate primary position:
        G4ThreeVector position = geom_->GenerateVertex(region_);
        G4double time=0;

        G4PrimaryVertex *vertex = new G4PrimaryVertex(position,time);

        //generate alpha particle with random energy from the selected energy range and add it to the vertex:
        double ekin = energy_low_+G4UniformRand()*(energy_high_-energy_low_);
        G4double mass = particle_defalpha->GetPDGMass();
        G4double energy = ekin + mass;
        G4double pmod = std::sqrt(energy*energy - mass*mass);
        G4ThreeVector momentum_dir = G4RandomDirection();
        G4double px = pmod * momentum_dir.x();
        G4double py = pmod * momentum_dir.y();
        G4double pz = pmod * momentum_dir.z();
        G4PrimaryParticle* particle =
        new G4PrimaryParticle(particle_defalpha);
        particle->SetMomentum(px, py, pz);
        particle->SetPolarization(0.,0.,0.);
        particle->SetProperTime(time);
        vertex->SetPrimary(particle);

        evt->AddPrimaryVertex(vertex);
    }
}