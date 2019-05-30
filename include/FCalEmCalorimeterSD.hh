/* This is the header file of FCal Sensitive detector class
Author: Zhaoyuan Cui (Maxwell)

*/

#ifndef FCalEmCalorimeterSD_hh
#define FCalEmCalorimeterSD_hh 1

#include "G4VSensitiveDetector.hh"

#include "FCalEmCalorimeterHit.hh"
#include "G4Accumulable.hh"
#include "G4AccumulableManager.hh"

#include <vector>


class G4Step;
class G4HCofThisEvent;


class FCalEmCalorimeterSD : public G4VSensitiveDetector
{
    public:
        FCalEmCalorimeterSD(
            const G4String& name, const G4String& hitsCollectionName);
        virtual ~FCalEmCalorimeterSD();

        // Methods from base class.
        virtual void Initialize(G4HCofThisEvent* hitCollection);
        virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
        virtual void EndOfEvent(G4HCofThisEvent* hitCollection);
        virtual void GetETot(G4double *ETot_Lar_Tp);

    private:
        FCalEmCalorimeterHitsCollection* fHitsCollection;
        std::vector<std::vector<std::vector<G4double>>>* fDoseRZ;
        std::vector<G4double>* ED_Lar_Tp;
        G4double ETot_Lar_Tp;
        G4double ETot_Lar;
        G4double ETot_Lar_CTub;
        G4Accumulable<G4double> ETot_Lar_CTub_Acc;
};

#endif
