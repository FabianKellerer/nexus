//------------------------------------------------------------------
// Nexus | Alpha particle generator
//
// This class produces alpha particles of a specific energy range
//
//------------------------------------------------------------------

#ifndef Alpha_GENERATOR_H
#define Alpha_GENERATOR_H

#include <G4VPrimaryGenerator.hh>

class G4Event;
class G4ParticleDefinition;
class G4GenericMessenger;

namespace nexus {
    class GeometryBase;

    class AlphaGenerator : public G4VPrimaryGenerator {
        public:
        //Constructor
        AlphaGenerator();
        //Destructor
        ~AlphaGenerator();

        void GeneratePrimaryVertex(G4Event* evt);

        private:

        G4GenericMessenger* msg_;
        const GeometryBase* geom_;

        G4double energy_low_;
        G4double energy_high_;
        G4String region_;
        G4ParticleDefinition* particle_defalpha;


    };
}// end namespace nexus
#endif