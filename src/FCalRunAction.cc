/*
This is the source code of run action of FCal
Author: Zhaoyuan Cui

*/

#include "FCalRunAction.hh"

#include "G4Run.hh"
#include "G4Event.hh"


//// Constructor and Destructor
FCalRunAction::FCalRunAction() : G4UserRunAction() {}
FCalRunAction::~FCalRunAction() {}


G4Run* FCalRunAction::GenerateRun()
{
	return G4UserRunAction::GenerateRun();
}


void FCalRunAction::BeginOfRunAction(const G4Run*) {}


void FCalRunAction::EndOfRunAction(const G4Run* run) {}
