/* This is the source file of SD class of FCal
Author: Zhaoyuan Cui (Maxwell)

Edited by Anson Kost with the help of Professor John Rutherfoord, May 2019.

*/

#include "FCalEmCalorimeterHit.hh"
#include "FCalEmCalorimeterSD.hh"
#include "FCalAnalysis.hh"

#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4Accumulable.hh"
#include "G4AccumulableManager.hh"

#include <fstream>


///|////////////////////////////////////////////////////////////////////////////
//|| (Repeated) Detector Geometry Parameters
///|////////////////////////////////////////////////////////////////////////////
const int numCals = 4;
const int tungPN = 8;
const int tungBPN = 4;


//// Constructor
FCalEmCalorimeterSD::FCalEmCalorimeterSD(
    const G4String& name, const G4String& hitsCollectionName
) : G4VSensitiveDetector(name), fHitsCollection(NULL)
{
    collectionName.insert(hitsCollectionName);
}


//// Destructor
FCalEmCalorimeterSD::~FCalEmCalorimeterSD() {}


void FCalEmCalorimeterSD::Initialize(G4HCofThisEvent* hce)
{
    // Create hits collection.
    fHitsCollection = new FCalEmCalorimeterHitsCollection(
        SensitiveDetectorName, collectionName[0]);

    // Add this collection in hce.
    G4int hcID \
        = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    hce->AddHitsCollection(hcID, fHitsCollection);
}


G4bool FCalEmCalorimeterSD::ProcessHits(G4Step* aStep,
					G4TouchableHistory*)
{
    G4double edep = aStep->GetTotalEnergyDeposit();  // Energy deposit.

    if (edep == 0.) return false;  // Skip if no energy deposit.

    FCalEmCalorimeterHit* newHit = new FCalEmCalorimeterHit();

    newHit->SetEdep(edep);
    newHit->SetPos(aStep->GetPostStepPoint()->GetPosition());
    // newHit->SetMomentum(aStep->GetPostStepPoint()->GetMomentum());
    // newHit->SetTotalEnergy(aStep->GetPostStepPoint()->GetTotalEnergy());
    // newHit->SetTrackID(aStep->GetTrack()->GetTrackID());

    fHitsCollection->insert(newHit);
    return true;
}


void FCalEmCalorimeterSD::EndOfEvent(G4HCofThisEvent*)
{

	///|////////////////////////////////////////////////////////////////////
	//|| Data Output Filename(s)
	///|////////////////////////////////////////////////////////////////////
	std::string hitsFilename = "data/hits.csv";

    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    G4int eID = 0;
    const G4Event* evt = G4RunManager::GetRunManager()->GetCurrentEvent();
    if (evt) eID = evt->GetEventID();
    G4cout << "Finished processing run " << eID << '.' << G4endl;

    G4int nofHits = fHitsCollection->entries();

    std::ofstream file(hitsFilename);
	if (file.is_open())
		G4cout << "Will write data to " << hitsFilename << '.' << G4endl;
	else
		G4cerr << "Couldn't open " << hitsFilename << " for writing data."
			<< G4endl;

    G4double xdep, ydep, zdep, edep;
    for (G4int i = 0; i < nofHits; i++)
    {
        FCalEmCalorimeterHit* hit = (*fHitsCollection)[i];
        G4ThreeVector localPos = hit->GetPos();

        xdep = localPos.x();
        ydep = localPos.y();
        zdep = localPos.z();
        edep = hit->GetEdep();

        // Python!
        file << edep << ", " << xdep << ", " << ydep << ", " << zdep
        /*  << ", " << hit->GetTotalEnergy()
            << ", " << hit->GetMomentum().x()
            << ", " << hit->GetMomentum().y()
            << ", " << hit->GetMomentum().z()
            << ", " << hit->GetTrackID()  */
            << G4endl;
    }

    file.close();
}
