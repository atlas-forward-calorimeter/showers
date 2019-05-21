//This is the source code of Primary Generator Action class
//Author: Zhaoyuan Cui

#include"FCalPrimaryGeneratorAction.hh"
#include"G4Event.hh"
#include"G4ParticleGun.hh"
#include"G4ParticleTable.hh"
#include"G4ParticleDefinition.hh"
#include"globals.hh"
#include"G4ThreeVector.hh"
#include "Randomize.hh"


FCalPrimaryGeneratorAction::FCalPrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction(),
    fParticleGun(0)
{
  G4int nofParticles=1;
  fParticleGun=new G4ParticleGun(nofParticles);

  //Find particle from particle table
    G4String particleName;

  G4ParticleTable*particleTable=G4ParticleTable::GetParticleTable();
  G4ParticleDefinition *particle=particleTable->FindParticle(particleName="e-");

  //Initialize particle parameters
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleEnergy(200.*GeV);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
}

FCalPrimaryGeneratorAction::~FCalPrimaryGeneratorAction()
{
  delete fParticleGun;
}

void FCalPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    fParticleGun->SetParticlePosition(G4ThreeVector(G4RandGauss::shoot(0,0.03)*CLHEP::cm,G4RandGauss::shoot(0,0.03)*CLHEP::cm,20.0*CLHEP::cm));

  fParticleGun->GeneratePrimaryVertex(anEvent);
}
