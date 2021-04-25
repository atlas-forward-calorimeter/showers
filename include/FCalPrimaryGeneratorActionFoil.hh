/*
This is the header file for Primary Generator Action class

*/

#ifndef FCalPrimaryGeneratorActionFoil_h
#define FCalPrimaryGeneratorActionFoil_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"


class G4ParticleGun;
class G4Event;


class FCalPrimaryGeneratorActionFoil : public G4VUserPrimaryGeneratorAction
{
    public:
        FCalPrimaryGeneratorActionFoil();
        ~FCalPrimaryGeneratorActionFoil();

        virtual void GeneratePrimaries(G4Event*);

    private:
        G4ParticleGun* fParticleGun;
};

#endif
