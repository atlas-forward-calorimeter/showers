/*
This is the run action header file of FCal
Author: Zhaoyuan Cui

*/

#ifndef FCalRunAction_hh
#define FCalRunAction_hh 1

#include "G4UserRunAction.hh"
#include "globals.hh"


class G4Run;


class FCalRunAction : public G4UserRunAction
{
    public:
        FCalRunAction();
        virtual ~FCalRunAction();

		virtual G4Run* GenerateRun();
        virtual void BeginOfRunAction(const G4Run*);
        virtual void EndOfRunAction(const G4Run*);
};

#endif
