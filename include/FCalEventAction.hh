/*
This is the header file of event action of FCal.
Author: Zhaoyuan Cui

*/

#ifndef FCalEventAction_hh
#define FCalEventAction_hh 1

#include "FCalEmCalorimeterHit.hh"

#include "G4UserEventAction.hh"
#include "globals.hh"

class FCalEventAction : public G4UserEventAction
{
    public:
        FCalEventAction();
        virtual ~FCalEventAction();

        virtual void BeginOfEventAction(const G4Event*);
        virtual void EndOfEventAction(const G4Event*);

        void WriteHits(G4VHitsCollection* hitsCollection, G4int eventID);
        
    private:
        G4int fECHCID;
};

#endif
