//This is the header file of drift signal
//Author: Zhaoyuan Cui

#ifndef driftSignal_hh
#define driftSignal_hh

#include "globals.hh"

class driftSignal
{
public:

  //Default constructor and destructor
  driftSignal();
  ~driftSignal();

  //Customized constructor
  driftSignal(G4double e, G4double p);

  //The funcion of pulse
  void driftFunction(G4double *res);

  void get_info();

private:

  //Parameter for pulse function
  G4double a;
  G4double b;

  //Total energy of the siganl
  G4double energy;

  //The postition of the signal
  G4double position;

  //The time needed to travle to the position of the signal
  G4double collision_time;

  //Velocity of signal propagation and speed of light
  const double v=240; //Unit: mm/ns
  const double c=300;
 
  //Energy readout method
  void readout(G4double s_E, G4double s_T, G4double *res);
};

#endif
