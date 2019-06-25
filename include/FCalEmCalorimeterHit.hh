/*
This is the header file of Hit class of FCal sensitive detector
Author: Zhaoyuan Cui (Maxwell)

Edited by Anson Kost with the help of Professor John Rutherfoord, May 2019.

*/

#ifndef FCalEmCalorimeterHit_hh
#define FCalEmCalorimeterHit_hh 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"


class FCalEmCalorimeterHit : public G4VHit
{
    public:
        // Constructors and Destructor
        FCalEmCalorimeterHit();
        FCalEmCalorimeterHit(const FCalEmCalorimeterHit&);
        virtual ~FCalEmCalorimeterHit();
        
        // Operators
        const FCalEmCalorimeterHit& operator=(const FCalEmCalorimeterHit&);
        G4int operator==(const FCalEmCalorimeterHit&) const;

        inline void *operator new(size_t);
        inline void operator delete(void*);

        // Methods from base class.
        virtual void Draw();
        virtual void Print();

        // Setters
        void SetEdep(G4double de) { fEdep = de; };
        void SetPos(G4ThreeVector xyz) { fPos = xyz; };
        void SetKineticEnergy(G4double energy) { fKineticEnergy = energy; };
        void SetParticleName(G4String name) { fParticleName = name; };
        // void SetMomentum(G4ThreeVector xyz) { fMomentum = xyz; };
        // void SetTotalEnergy(G4double energy) { fTotalEnergy = energy; };
        void SetTrackID(G4int id) { fTrackID = id; };
        // void SetParticleName(G4String& name) { fParticleName = name; };

        // Getters
        G4double GetEdep() const { return fEdep; };
        G4ThreeVector GetPos() const { return fPos; };
        G4double GetKineticEnergy() const { return fKineticEnergy; };
        G4String GetParticleName() const { return fParticleName; };
        // G4ThreeVector GetMomentum() const { return fMomentum; };
        // G4double GetTotalEnergy() const { return fTotalEnergy; };
        G4int GetTrackID() const { return fTrackID; };
        // G4String& GetParticleName() const { return fParticleName; };

    private:
        // Data to Keep
        G4double fEdep;
        G4ThreeVector fPos;
        G4double fKineticEnergy;
        G4String fParticleName;
        // G4ThreeVector fMomentum;
        // G4double fTotalEnergy;
        G4int fTrackID;
        // G4String& fParticleName;
};


typedef G4THitsCollection<FCalEmCalorimeterHit> FCalEmCalorimeterHitsCollection;

extern G4ThreadLocal G4Allocator<FCalEmCalorimeterHit>* \
    FCalEmCalorimeterHitAllocator;


inline void* FCalEmCalorimeterHit::operator new(size_t)
{
    if (!FCalEmCalorimeterHitAllocator)
        FCalEmCalorimeterHitAllocator = new G4Allocator<FCalEmCalorimeterHit>;
    return (void*)FCalEmCalorimeterHitAllocator->MallocSingle();
}


inline void FCalEmCalorimeterHit::operator delete(void *hit)
{
    FCalEmCalorimeterHitAllocator->FreeSingle((FCalEmCalorimeterHit*)hit);
}

#endif
