/*
`main`: Where it all starts.

Edited by Anson Kost with the help of Professor John Rutherfoord, May 2019.

*/

#include"FCalDetectorConstruction.hh"
#include"FCalActionInitialization.hh"

#ifdef G4MULTITHREADED
#include"G4MTRunManager.hh"
#else
#include"G4RunManager.hh"
#endif

#include"G4UImanager.hh"
#include"FTFP_BERT.hh"

#include"G4VisExecutive.hh"
#include"G4UIExecutive.hh"

#include "Randomize.hh"

#include <chrono>

int main(int argc, char** argv)
{
    // For timing the whole program.
    const time_t startTime = time(NULL);

    // Detect interactive mode and define UI session.
    G4UIExecutive* ui = 0;
    if (argc == 1)
        // Let G4UIExecutive guess the best available UI.
        ui = new G4UIExecutive(argc, argv);
    else if (strcmp(argv[1], "i") == 0)
        // Interactive mode with tcsh.
        ui = new G4UIExecutive(argc, argv, "tcsh");

    // Choose the Random engine.
    G4Random::setTheEngine(new CLHEP::RanecuEngine);

#ifdef G4MULTITHREADED
    G4MTRunManager* runManager = new G4MTRunManager;
    ///|////////////////////////////////////////////////////////////////////
    //|| Number of Cores
    ///|////////////////////////////////////////////////////////////////////
    // Set the default number of threads to be the number of available cores of 
    // the machine.
    runManager->SetNumberOfThreads(32);
#else
    G4RunManager* runManager = new G4RunManager;
#endif

    // Set mandatory user initialization classes:

    // The Geometry
    runManager->SetUserInitialization(new FCalDetectorConstruction());

    // The Physics
    G4VModularPhysicsList* physicsList = new FTFP_BERT;
    runManager->SetUserInitialization(physicsList);

    // User Action Initialization
    runManager->SetUserInitialization(new FCalActionInitialization());

    // Initialize Visualization
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();

    // Get the pointer to the user interface manager.
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    if (!ui)
    {
        // Batch mode.
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command + fileName);
    }
    else
    {
        // Interactive mode.
        UImanager->ApplyCommand("/control/execute init_vis.mac");
        if (ui->IsGUI())
        {
            // UImanager->ApplyCommand("/control/execute gui.mac");
        }
        ui->SessionStart();
        delete ui;
    }

    const time_t endTime = time(NULL);
    G4cout << "This run took " << (endTime - startTime) 
           << " seconds in all. It finished on " << std::ctime(&endTime);

    // Job termination.
    // Free the store: user actions, physics_list and detector_description are 
    // owned and deleted by the run manager, so they should not be deleted
    // in the `main` program!
    delete visManager;
    delete runManager;
}
