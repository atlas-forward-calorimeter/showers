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
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., -1.));
}

//// Destructor
FCalPrimaryGeneratorAction::~FCalPrimaryGeneratorAction()
{
  delete fParticleGun;
}

void FCalPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    G4double tubtriside = 7.5 * CLHEP::mm;

    // x, y for the "upper right" (+x and +y) tube.
    G4double offset[2] = { tubtriside / 4 * 3, tubtriside / 4 * sqrt(3) };

    fParticleGun->SetParticlePosition(G4ThreeVector(
        offset[0] + G4RandGauss::shoot(0, 0.03) * CLHEP::cm,
        offset[1] + G4RandGauss::shoot(0, 0.03) * CLHEP::cm,
        20.0 * CLHEP::cm
    ));

    fParticleGun->GeneratePrimaryVertex(anEvent);
}
