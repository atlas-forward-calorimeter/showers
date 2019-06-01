/*
This is the source code of event action of FCal.
Author: Zhaoyuan Cui

*/

#include "FCalEventAction.hh"
#include "FCalEmCalorimeterHit.hh"
#include "FCalAnalysis.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"


FCalEventAction::FCalEventAction() : G4UserEventAction(), fECHCID(-1) {}

FCalEventAction::~FCalEventAction() {}

void FCalEventAction::BeginOfEventAction(const G4Event*) {}

void FCalEventAction::EndOfEventAction(const G4Event*) {}
