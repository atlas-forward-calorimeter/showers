/*
This is the source code of event action of FCal.
Author: Zhaoyuan Cui

*/

#include "FCalEventAction.hh"
#include "FCalEmCalorimeterHit.hh"

#include "G4Event.hh"
#include "G4SDManager.hh"

#include <fstream>
#include <sstream>
#include <ctime>


//// Constructor and Destructor
FCalEventAction::FCalEventAction() : G4UserEventAction(), fECHCID(-1) {}
FCalEventAction::~FCalEventAction() {}


void FCalEventAction::BeginOfEventAction(const G4Event*) {}


void FCalEventAction::EndOfEventAction(const G4Event* event) {
    // Get hit collection(s) by name.
    auto sdManager = G4SDManager::GetSDMpointer();
    auto collID = sdManager->GetCollectionID("FCalHitsCollection");
    auto hitCollection = event->GetHCofThisEvent()->GetHC(collID);

    G4int eventID = event->GetEventID();

    // Print out useful information.
    G4cout << "Event (" << eventID << ") simulated. (" 
           << hitCollection->GetSize() << ") hits."
           << G4endl;

    // Write data to file.
    WriteHits(hitCollection, eventID);
}


void FCalEventAction::WriteHits(
    G4VHitsCollection* hitsCollection, G4int eventID)
{
    ///|////////////////////////////////////////////////////////////////////
    //|| Data Folder and Filenames
    ///|////////////////////////////////////////////////////////////////////
    const G4String outFolderPath = "data";
    const G4String filePrefix = "hits";
    const G4String fileExtension = ".csv";

    ///|////////////////////////////////////////////////////////////////////
    //|| .csv File Comments
    ///|////////////////////////////////////////////////////////////////////
    time_t now = time(NULL);
    G4String nowString = std::ctime(&now);
    G4String comments = \
        "FCal simulation hits output. " + nowString;
    //

    // Full relative path for this event's data file.
    const G4String filePath = outFolderPath 
                            + '/' 
                            + filePrefix 
                            + '-' 
                            + std::to_string(eventID) 
                            + fileExtension;

    // Hopefully open the file.
    std::ofstream file(filePath);
    if (!file.is_open()) {
        G4cerr << "ERROR! Could not open " << filePath << '.' << G4endl;
        return;
    }

    // Data column headers.
    G4String header = "energy_deposit,x,y,z,particle_name,track_id,kinetic_energy";

    // Write to a stringstream first (for speed).

    std::stringstream fileContents;
    fileContents << comments << header << G4endl;

    FCalEmCalorimeterHit* hit = 0;
    G4ThreeVector position;
    for (size_t i = 0; i < hitsCollection->GetSize(); i++) {
        // Loop over hits.
        hit = static_cast<FCalEmCalorimeterHit*>(hitsCollection->GetHit(i));
        position = hit->GetPos();
        fileContents << hit->GetEdep() / CLHEP::MeV
                     << ',' << position.x() / CLHEP::mm
                     << ',' << position.y() / CLHEP::mm
                     << ',' << position.z() / CLHEP::mm
                     << ',' << hit->GetParticleName()
                     << ',' << hit->GetTrackID()
                     << ',' << hit->GetKineticEnergy() / CLHEP::MeV
                     << G4endl;
    }

    // Then write the stringstream to the file.
    std::ofstream outFile(filePath);
    outFile << fileContents.rdbuf();
    outFile.close();

    G4cout << "Wrote to " << filePath << '.' << G4endl;
}
