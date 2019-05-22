/*
This is the header file of user action class.
Author: Zhaoyuan Cui

*/

#ifndef FCalActionInitialization_hh


#define FCalActionInitialization_hh 1

#include "G4VUserActionInitialization.hh"

class FCalDetectorConstruction;

class FCalActionInitialization : public G4VUserActionInitialization
{
public:
    FCalActionInitialization();
    virtual ~FCalActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;
};


#endif
