//This is the header file of Hit class of FCal sensitive detector
//Author: Zhaoyuan Cui (Maxwell)

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
  FCalEmCalorimeterHit();
  FCalEmCalorimeterHit(const FCalEmCalorimeterHit&);
  virtual ~FCalEmCalorimeterHit();

  //Operators
  const FCalEmCalorimeterHit& operator=(const FCalEmCalorimeterHit&);
  G4int operator==(const FCalEmCalorimeterHit&) const;

  inline void *operator new(size_t);
  inline void operator delete(void*);

  //Methods from base class
  virtual void Draw();
  virtual void Print();

  //Set methods
  void SetEdep(G4double de){fEdep=de;};
  void SetPos(G4ThreeVector xyz){fPos=xyz;};

  //Get methods
  G4double GetEdep() const {return fEdep;};
  G4ThreeVector GetPos()  const {return fPos;};

private:
  G4double fEdep;
  G4ThreeVector fPos;
};

typedef G4THitsCollection<FCalEmCalorimeterHit> FCalEmCalorimeterHitsCollection;

extern G4ThreadLocal G4Allocator<FCalEmCalorimeterHit>* FCalEmCalorimeterHitAllocator;

inline void* FCalEmCalorimeterHit::operator new(size_t)
{
  if(!FCalEmCalorimeterHitAllocator)
    FCalEmCalorimeterHitAllocator=new G4Allocator<FCalEmCalorimeterHit>;
  return (void*) FCalEmCalorimeterHitAllocator->MallocSingle();
}

inline void FCalEmCalorimeterHit::operator delete(void *hit)
{
  FCalEmCalorimeterHitAllocator->FreeSingle((FCalEmCalorimeterHit*) hit);
}

#endif
