//This is the source code of run action of FCal
//Author: Zhaoyuan Cui

#include "FCalRunAction.hh"
#include "FCalAnalysis.hh"
#include "FCalEmCalorimeterSD.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

FCalRunAction::FCalRunAction()
  :G4UserRunAction()
{
  G4AnalysisManager* analysisManager=G4AnalysisManager::Instance();
  G4cout<<"Using "<<analysisManager->GetType()<<G4endl;

  //Default settings
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFileName("1TeV_simulation");

  //Book histogram, ntuple
  //

  //Creating 2D histogram
  /*
  analysisManager
    ->CreateH2("FCal_sig_XY","FCal signal X vs Y",
	       450,-450,450,450,-450,450,"mm","mm");
  analysisManager
    ->CreateH2("FCal_sig_YZ","FCal signal Y vs Z",
	       450,0,450,900,-450,450,"mm","mm");
  */
  //bruce setting up histograms
  G4double rMax = 0.2625*CLHEP::cm;
  G4double rMin = 0.2356*CLHEP::cm;
  /*
  analysisManager
    ->CreateH1("Tube1","Tube1 energy deposit",
	       64,rMin,rMax,"cm");
  */
  // G4cout<<"after createH1 in Runaction"<<G4endl;
  /*
  analysisManager
    ->CreateH1("Tube1test","Tube1test energy deposit",
	       10,0,10);
  */
  int runn =10;
  for(int eIDH=0; eIDH<runn; eIDH++)
    {
      char Hisname [100];
      sprintf(Hisname,"event %i: all SD hits distribution in Z", eIDH);
  analysisManager
    ->CreateH1(Hisname,Hisname,
	       1100,-5.5,5.5);
  // sprintf(Hisname,"event %i: all SD energy distribution in Z and rho ", eIDH);

  //  analysisManager
  // ->CreateH2(Hisname,Hisname,
  //	       1100,-5.5,5.5, 300,0,3);
    }
  for(int eIDH=0; eIDH<runn; eIDH++)
    {
            char Hisname [100];
      sprintf(Hisname,"event %i: energy deposit in LAr plate on xy plane", eIDH);
    analysisManager
    ->CreateH2(Hisname, Hisname,
	       300,-1.5,1.5, 300,-1.5,1.5);
    }
  analysisManager
    ->CreateH1("all SD hits distribution in Z","all SD hits distribution in Z",
	       1100,-5.5,5.5); 
   analysisManager
    ->CreateH2("all SD energy distribution in Z and rho ","all SD energy distribution in Z and rho ",
	       1100,-5.5,5.5, 300,0,3);
   

   
}


FCalRunAction::~FCalRunAction()
{
  delete G4AnalysisManager::Instance();
}

void FCalRunAction::BeginOfRunAction(const G4Run*)
{
  //Get analysis manager
  G4AnalysisManager* analysisManager=G4AnalysisManager::Instance();

  //Output file
  analysisManager->OpenFile();
}

void FCalRunAction::EndOfRunAction(const G4Run*)
{
  //Write and close output file
  //save histogram & ntuple
  //
  G4AnalysisManager* analysisManager=G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
  G4double Etot_Lar_Tp;
  //FCalEmCalorimeterSD::GetETot(&Etot_Lar_Tp);
  //G4cout<<"Etot in Lar tungsten plate is "<<Etot_Lar_Tp<<G4endl;
  
}
