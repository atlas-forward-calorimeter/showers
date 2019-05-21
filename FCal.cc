//This is an example main program code 
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

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main(int argc, char** argv)
{
  //Detect interactive mode and define UI session
  G4UIExecutive* ui=0;
  if(argc==1)
    {
      //Let G4UIExecutive guess what is the best available UI
      ui=new G4UIExecutive(argc,argv);
    }

  //Choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

#ifdef G4MULTITHREADED
  G4MTRunManager* runManager=new G4MTRunManager;
  //Set the default number of threads to be the number of available cores of the machine
  runManager->SetNumberOfThreads(4);
#else
  G4RunManager* runManager=new G4RunManager;
#endif

  //Set Mandatory user initialization classes
  //

  //The Geometry
  runManager->SetUserInitialization(new FCalDetectorConstruction());

  //The Physics
  G4VModularPhysicsList* physicsList=new FTFP_BERT;
  runManager->SetUserInitialization(physicsList);

  //User action initialization
  runManager->SetUserInitialization(new FCalActionInitialization());

  //Initialize visualization
  //
  G4VisManager* visManager=new G4VisExecutive;
  visManager->Initialize();

  //Get the pointer to the User Interface manager
  G4UImanager* UImanager=G4UImanager::GetUIpointer();

  if(!ui)
    {
      //batch mode
      G4String command="/control/execute ";
      G4String fileName=argv[1];
      UImanager->ApplyCommand(command+fileName);
    }
    else
      {
	//interactive mode
	UImanager->ApplyCommand("/control/execute init_vis.mac");
	if(ui->IsGUI())
	  { 
	    //UImanager->ApplyCommand("/control/execute gui.mac");
	  }
	ui->SessionStart();
	delete ui;
      }
  
  //Job termination
  //Free the store: user actions, physics_list and detector_description are 
  //owned and deleted by the run mangager, so they should not be deleted
  //in the main() program!

  delete visManager;
  delete runManager;

}
