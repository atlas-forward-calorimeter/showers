//This is the source file of SD class of FCal
//Author: Zhaoyuan Cui (Maxwell)

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



#include "driftSignal.hh"

FCalEmCalorimeterSD::FCalEmCalorimeterSD( const G4String& name,
					  const G4String& hitsCollectionName)
  : G4VSensitiveDetector(name),
    ETot_Lar_CTub_Acc("ETot_Lar_CTub_Acc", 0.),
    fHitsCollection(NULL)
{
  collectionName.insert(hitsCollectionName);
  const int nbins = 64;
  const int nrings = 3;
  const int ntubs = 4;
  const int nLar_Tp = 10;
  std::vector<G4double> fDose0(nbins,0);
  std::vector< std::vector<G4double> > fDose1(nrings, fDose0);
  fDoseRZ = new std::vector< std::vector< std::vector<G4double> > > (ntubs, fDose1);
  ED_Lar_Tp = new std::vector<G4double> (nLar_Tp ,0);
  ETot_Lar_Tp=0.;
  ETot_Lar=0.;
  ETot_Lar_CTub=0.;

  G4AccumulableManager* AccumulableManager = G4AccumulableManager::Instance();
  AccumulableManager->RegisterAccumulable(ETot_Lar_CTub_Acc);

  

}

FCalEmCalorimeterSD::~FCalEmCalorimeterSD(){
  if ( fDoseRZ ) {
    delete fDoseRZ;
    fDoseRZ = 0;
  }
  if ( ED_Lar_Tp ) {
    delete ED_Lar_Tp;
    ED_Lar_Tp = 0;
  }

}

void FCalEmCalorimeterSD::Initialize(G4HCofThisEvent* hce)
{
  //Create hits collection
  fHitsCollection
    =new FCalEmCalorimeterHitsCollection(SensitiveDetectorName, collectionName[0]);

  //Add this collection in hce
  G4int hcID
    =G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection(hcID, fHitsCollection);
  for(unsigned int k=0;k<fDoseRZ->size();k++)
    {
      for(unsigned int j=0;j<(*fDoseRZ)[k].size();j++)
	{
	  for(unsigned int i=0;i<(*fDoseRZ)[k][j].size();i++) {
	    (*fDoseRZ)[k][j][i] = 0.;
	  }
	}
    }
  for(unsigned int k=0;k<ED_Lar_Tp->size();k++)
    {
      (*ED_Lar_Tp)[k] = 0.;
    }

  
}

G4bool FCalEmCalorimeterSD::ProcessHits(G4Step* aStep,
					G4TouchableHistory*)
{
  //Energy deposit
  G4double edep=aStep->GetTotalEnergyDeposit();
  
  if(edep==0.)
    return false;

  FCalEmCalorimeterHit* newHit=new FCalEmCalorimeterHit();
  
  newHit->SetEdep(edep);
  newHit->SetPos(aStep->GetPostStepPoint()->GetPosition());

  fHitsCollection->insert(newHit);
  //G4cout<<"processing hits"<<G4endl;

  //newHit->Print();

  return true;
}
void FCalEmCalorimeterSD::GetETot(G4double *ETot_Lar_Tp)
{
  *ETot_Lar_Tp=0.;


}

void FCalEmCalorimeterSD::EndOfEvent(G4HCofThisEvent*)
{
  G4AnalysisManager* analysisManager=G4AnalysisManager::Instance();
  G4int eID=0;
  const G4Event* evt= G4RunManager::GetRunManager()->GetCurrentEvent();
  if(evt) eID = evt->GetEventID();
  G4cout<<"finished processing run number "<<eID<<G4endl;
  G4int nofHits=fHitsCollection->entries();
  G4double E_total;
  G4double E_pos;
  G4double E_ave;
  E_total=0;
  E_pos=0;
  E_ave=0;
  G4double ETot_Lar_Tpf=0;
  G4double ETot_Lar_Tpb=0;
    ETot_Lar_CTub=0.;

  

  G4double reading[90000];

  for(G4int i=0;i<90000;i++)
    {
      reading[i]=0;
    }
  //bruce adding 4 tube setup
  G4double zMax = 3.5/2.*CLHEP::cm; 
  G4double rMax = 0.2625*CLHEP::cm;   
  G4double rMin = 0.2356*CLHEP::cm;  

  //bruce getting coords for new shifted tubes
  const int tubcount=4;
  double tubtriside=0.75*CLHEP::cm;
  double tshifty[tubcount];
  double tshiftx[tubcount];
  double tshiftz[tubcount];
  double rde[tubcount];
  double zde[tubcount];
  for(int it=0; it<tubcount;it++)
    {
      tshifty[it]=-tubtriside/4.0*sqrt(3.)*pow(-1,it);
      tshiftx[it]=tubtriside/4.0*(it*2-3);
      tshiftz[it]=0.0*CLHEP::cm;

    }
  //adding Tp setup
  G4double boxhz=4.7/2.;
  G4double tunghzTot=1.;
  G4double larGThz=0.1;
  G4double larGThzvar=0.02*2;
  int tungPN=6;
  const int LarTpn=10;
  double zLarTpc[LarTpn];
  for(int il=0; il<LarTpn; il++)
    {
      if(il<LarTpn/2)
	{
	  zLarTpc[il]=-(boxhz+tunghzTot/(double)tungPN*2*(tungPN-1)+larGThz*2*(tungPN-2)+larGThz)+il*(tunghzTot/tungPN+larGThz)*2;
	  
	  //G4cout<<"in SD define Tungplate series: z ceter for Lar gap is "<<zLarTpc[il]<<G4endl;// test if tungplate coordinates are correct as in geometry

	}
      else
	{
	  zLarTpc[il]=boxhz+tunghzTot/(double)tungPN*2+(il-LarTpn/2)*(tunghzTot/tungPN+larGThz)*2+larGThz;
	  // G4cout<<"in SD define Tungplate series: z ceter for Lar gap is "<<zLarTpc[il]<<G4endl;// test if tungplate coordinates are correct as in geometry

	}

    }


  G4double xdep,  ydep,  zdep, edep;
  for(G4int i=0;i<nofHits;i++)
    {
      //G4cout<<"in nofHits i="<<i<<G4endl;

      FCalEmCalorimeterHit *hit=(*fHitsCollection)[i];
      G4ThreeVector localPos=hit->GetPos();
      //analysisManager->FillH2(0,localPos.x(),localPos.y());
      //analysisManager->FillH2(1,localPos.z(),localPos.y());
      
      // G4cout<<"in nofHits after GetPos"<<G4endl;

      E_total+=hit->GetEdep();
      E_pos+=hit->GetEdep()*localPos.z();
      //bruce tube shifted scoring
      //G4cout<<"in nofHits after GetEdep"<<G4endl;

      xdep=localPos.x();
      ydep=localPos.y();
      zdep=localPos.z();
      edep=hit->GetEdep();
      ETot_Lar+=edep;
     for(int it=0; it<tubcount;it++)
       {
	 rde[it]=sqrt(pow(xdep-tshiftx[it],2)+pow(ydep-tshifty[it],2));
	 zde[it]=zdep-tshiftz[it];	 
       }

     //filling  tube ring bin vector
   for(int itub=0; itub<tubcount;itub++)
     {
       if(rde[itub]<=rMax && rde[itub]>=rMin &&fabs(zde[itub])<zMax)
	 {
	   int ir; //z position in the tube
	   // G4cout<<"zde is "<<zde[itub]<<"0.4*CLHEP::cm is "<<0.4*CLHEP::cm<<G4endl;
	   if(zde[itub]<-0.4*CLHEP::cm)
	     {
	       ir=0;
	       //  G4cout<<"setting ir to 0 "<<G4endl;

	     }
	   else if(zde[itub]<0.4*CLHEP::cm)
	     {
	       ir=1;
	       // G4cout<<"setting ir to 1 and edep is"<<edep<<G4endl;

	     }
	   else
	     {
	       ir=2;
	       // G4cout<<"setting ir to 2 "<<G4endl;

	     }
	    
	   int iBin=(int)((rde[itub]-rMin)/(rMax-rMin)*(*fDoseRZ)[itub][ir].size());
	   if ( iBin > -1 && iBin < (*fDoseRZ)[itub][ir].size() ) {
	     double rBin0 = (iBin+0.0)/(*fDoseRZ)[itub][ir].size()*rMax;
	     double rBin1 = (iBin+1.0)/(*fDoseRZ)[itub][ir].size()*rMax;
	     (*fDoseRZ)[itub][ir][iBin] += edep;
	     ETot_Lar_CTub+=edep;
	     ETot_Lar_CTub_Acc+=edep;


	   }
	 }


     }

   //filling LarGap Tung plate  vector
   for(int il=0; il<LarTpn;il++)
     {
       if(fabs(zdep/CLHEP::cm-zLarTpc[il])<=larGThz)
	 {
	   (*ED_Lar_Tp)[il]+=edep;
	   ETot_Lar_Tp+=edep;
	   if(zdep/CLHEP::cm>0)
	     {
	       //front part (particle incoming)
	       	   ETot_Lar_Tpf+=edep;


	     }
	   else if(zdep/CLHEP::cm<0)
	     {
	       //back part (particle outgoing)
	       ETot_Lar_Tpb+=edep;


	       
	     }
	 }

     }
  
      if(fabs(zdep/CLHEP::cm-zLarTpc[4])<=larGThz)
	 {
	   //H2 for each event in xy plane of one plate order in 10+eID
	   analysisManager->FillH2(eID,xdep/CLHEP::cm,ydep/CLHEP::cm,edep/CLHEP::keV);

	   //storing only plate 4 in LAr in the xy plane histogram
	 }
  
      if(rde[0]<=rMax && rde[0]>=rMin &&fabs(zde[0])<0.4*CLHEP::cm)
      	{
	  // G4cout<<"in condition"<<G4endl;
	  //analysisManager->FillH1(0,rde[0],edep/CLHEP::keV);
	  //1 copper tube histo
	  
      	}
      
      
      //driftSignal *p=new driftSignal(hit->GetEdep(),localPos.z());
      //p->driftFunction(reading);
      
      // delete p;
            //z distribution histogram fill

      if(fabs(zdep/CLHEP::cm)<boxhz-0.1)
	{
	  
      analysisManager->FillH1(eID,zdep/CLHEP::cm,edep*30/CLHEP::keV);

	}
      else
	{

      analysisManager->FillH1(eID,zdep/CLHEP::cm,edep/CLHEP::keV);

	}
      //H2 for each event in z rho plane order in eID
      // analysisManager->FillH2(eID,zdep/CLHEP::cm,sqrt(pow(xdep,2)+pow(ydep,2))/CLHEP::cm,edep/CLHEP::keV);

      analysisManager->FillH1(10,zdep/CLHEP::cm,edep/CLHEP::keV);
      analysisManager->FillH2(10,zdep/CLHEP::cm,sqrt(pow(xdep,2)+pow(ydep,2))/CLHEP::cm,edep/CLHEP::keV);
    } //hits loop end


  for(double ihb=0; ihb<10; ihb=ihb+1.0)
    {
      // analysisManager->FillH1(1,ihb+0.5,(*ED_Lar_Tp)[ihb]/CLHEP::keV);
      //tung plate LArGap series energy

    }
 
  E_ave=E_pos/E_total;

  //analysisManager->FillH1(0,E_ave);

  G4double ave;
  G4double te;
  G4double criticalTime;
  
  ave=0;
  te=0;
  criticalTime=0;
  
  for(G4int i=0;i<90000;i++)
    {
      ave+=((G4double)i+0.5)*10*reading[i];
      te+=reading[i];
    }
  criticalTime=ave/te;
  criticalTime=criticalTime/1000;

  G4cout<<"The time is: "<<criticalTime<<G4endl;
  
  //analysisManager->FillH1(1,criticalTime);
  G4cout<<"This is event number "<<eID<<G4endl;
  G4cout<<G4endl
	<<"----------->Hits Collection: in this event there are "<<nofHits
	<<" hits in the lAr gap: "<<"total energy in Lar tungplate is "<<ETot_Lar_Tp/CLHEP::keV<<" in Lar copper tube is "<<ETot_Lar_CTub/CLHEP::keV<<" tungplate+copper tube is "<<(ETot_Lar_CTub+ETot_Lar_Tp)/CLHEP::keV<<G4endl;
  G4cout<<"ratio of front/back in LAr plate series ="<<ETot_Lar_Tpf/ETot_Lar_Tpb<<" front energy deposit (keV) ="<<ETot_Lar_Tpf/CLHEP::keV<<" back energy deposit (keV) ="<<ETot_Lar_Tpb/CLHEP::keV<<G4endl;
  G4AccumulableManager* parameterManager = G4AccumulableManager::Instance();
  parameterManager->Merge();

    G4double ELCA  = ETot_Lar_CTub_Acc.GetValue();

  G4cout<<"accumulated energy in Lar copper tube over each event is"<<ELCA/CLHEP::keV<<G4endl;
}
