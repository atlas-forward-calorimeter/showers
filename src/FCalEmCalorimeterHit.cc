/*
This is the source file of Hit class of FCal sensitive detector
Author: Zhaoyuan Cui (Maxwell)

Edited by Anson Kost with the help of Professor John Rutherfoord, May 2019.

*/

#include "FCalEmCalorimeterHit.hh"
#include "G4UnitsTable.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"

#include <iomanip>

G4ThreadLocal G4Allocator<FCalEmCalorimeterHit>* \
    FCalEmCalorimeterHitAllocator = 0;


// Constructor
FCalEmCalorimeterHit::FCalEmCalorimeterHit()
    : G4VHit(), fEdep(0.), fPos(G4ThreeVector()) {}


// Constructor
FCalEmCalorimeterHit::FCalEmCalorimeterHit(const FCalEmCalorimeterHit& right)
    : G4VHit()
{
    fEdep = right.fEdep;
    fPos = right.fPos;
}


// Destructor
FCalEmCalorimeterHit::~FCalEmCalorimeterHit() {}


// =
const FCalEmCalorimeterHit& FCalEmCalorimeterHit::operator \
    =(const FCalEmCalorimeterHit& right)
{
    fEdep = right.fEdep;
    fPos = right.fPos;
    return *this;
}


// ==
G4int FCalEmCalorimeterHit::operator ==(const FCalEmCalorimeterHit& right) const
{
    return (this == &right) ? 1 : 0;
}


void FCalEmCalorimeterHit::Draw() {}


void FCalEmCalorimeterHit::Print()
{
    G4cout << "Edep:\t" << fEdep / MeV << " MeV"
        << " Z: " << fPos.getZ() / mm << " mm"
        << G4endl;
}
