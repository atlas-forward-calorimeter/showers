/*
This is the header file of user initialization class.
Author: Zhaoyuan Cui

*/

#include "FCalActionInitialization.hh"
#include "FCalPrimaryGeneratorAction.hh"
#include "FCalRunAction.hh"
#include "FCalEventAction.hh"


FCalActionInitialization::FCalActionInitialization()
    : G4VUserActionInitialization() {}

FCalActionInitialization::~FCalActionInitialization() {}


void FCalActionInitialization::BuildForMaster() const
{
    SetUserAction(new FCalRunAction);
}


void FCalActionInitialization::Build() const
{
    SetUserAction(new FCalPrimaryGeneratorAction);
    SetUserAction(new FCalRunAction);
    SetUserAction(new FCalEventAction);
}
