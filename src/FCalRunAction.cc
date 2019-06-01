/*
This is the source code of run action of FCal
Author: Zhaoyuan Cui

*/

#include "FCalRunAction.hh"
#include "FCalAnalysis.hh"
#include "FCalEmCalorimeterSD.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


//// Constructor
FCalRunAction::FCalRunAction() : G4UserRunAction() {
    G4cout << "Hi, I'm a new thread." << G4endl;
}


//// Destructor
FCalRunAction::~FCalRunAction()
{
    delete G4AnalysisManager::Instance();
}


G4Run* FCalRunAction::GenerateRun()
{
	return G4UserRunAction::GenerateRun();
}


void FCalRunAction::BeginOfRunAction(const G4Run*) {}


void FCalRunAction::EndOfRunAction(const G4Run*) {}
