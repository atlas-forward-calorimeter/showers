/*
This is the header file of user initialization class.
Author: Zhaoyuan Cui

-

Edited by Anson Kost with the help of Professor John Rutherfoord, March 2020.

*/

#include "FCalActionInitialization.hh"
#include "FCalPrimaryGeneratorAction.hh"
#include "FCalPrimaryGeneratorActionFoil.hh"
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
    SetUserAction(new FCalPrimaryGeneratorAction);  // Set source as single beta or radioactive foil.
    SetUserAction(new FCalRunAction);
    SetUserAction(new FCalEventAction);
}
