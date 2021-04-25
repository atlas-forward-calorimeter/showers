/*
This is the source code of Primary Generator Action class
Author: Zhaoyuan Cui

*/

#include "FCalPrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"


//// Constructor
FCalPrimaryGeneratorAction::FCalPrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction(), fParticleGun(0)
{
    G4int nofParticles = 1;
    fParticleGun = new G4ParticleGun(nofParticles);

    // Find particle from particle table.
    G4String particleName;

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle \
        = particleTable->FindParticle(particleName = "e-");

    // Initialize particle parameters.
    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticleEnergy(350. * GeV);

    // Straight down z-axis towards -z.
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., -1.));
}


//// Destructor
FCalPrimaryGeneratorAction::~FCalPrimaryGeneratorAction()
{
  delete fParticleGun;
}


void FCalPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    // Control the x and y offset of the electron beam. Set to (0, 0, 0) for a centered beam.
    G4ThreeVector offset = G4ThreeVector(0, 0, 0) * CLHEP::mm;

    G4ThreeVector position = offset + G4ThreeVector(  // 0.3 mm beam spread
        G4RandGauss::shoot(0, 0.3) * CLHEP::mm,
        G4RandGauss::shoot(0, 0.3) * CLHEP::mm,
        200 * CLHEP::mm
    );

    fParticleGun->SetParticlePosition(position);
    fParticleGun->GeneratePrimaryVertex(anEvent);
}
