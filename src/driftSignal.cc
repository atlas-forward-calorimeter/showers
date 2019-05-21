//Source code of driftSignal class
//Author: Zhaoyuan Cui

#include"driftSignal.hh"
#include<iostream>
#include<cmath>

//Default constructor does nothing
driftSignal::driftSignal()
{}

//Default destructor does nothing
driftSignal::~driftSignal(){}

//Customized constructor
driftSignal::driftSignal(G4double e, G4double p)
{
  //Initialize the total energy of the signal and position
  //Calculate the collision time
  //Calculate the parameters for pulse function
  energy=e;
  position=p;
  collision_time=position/c;
  a=-2*e/3721;
  b=-a*(collision_time+61);
}

//Energy readout function
void driftSignal::readout(G4double s_E, G4double s_T, G4double *res)
{
  //Define the split signals' energy
  G4double f_runtime_energy, s_runtime_energy;
  
  //The time for the first signal to travel to readout end of the detector
  G4double individual_time;

  //The time for the signal to finish a round trip
  G4double travel_time;

  //The time difference between sequencial readout
  G4double difference_time;

  //The accumulated time since the signal
  G4double accumulated_time;

  //The total energy readout from the detector
  G4double readout_energy;

  //Counter for time-energy array
  G4double counter,counter_i,counter_d;
  G4int c_n;

  //Initialization of variables
  readout_energy=0.0;
  f_runtime_energy=0.5*s_E;
  s_runtime_energy=0.5*s_E;
  accumulated_time=0.0;
  individual_time=position/v;
  travel_time=900/v;
  difference_time=(900-2*position)/v;

  //Start readout
  while(readout_energy<0.9*s_E)
    {
      if(readout_energy==0.0)
	{
	  //Time and energy for first signal
	  accumulated_time+=individual_time;
	  readout_energy+=0.1*f_runtime_energy;

	  //Store in the time-energy array

	  //|||||||||||||||||||||||||||||
	  counter=(accumulated_time+s_T)*100;
	  counter_d=std::modf(counter,&counter_i);
	  c_n=(int)counter_i;
	  res[c_n]+=0.1*f_runtime_energy;
	  //||||||||||||||||||||||||||||||

	  //Time and energy for second signal
	  accumulated_time+=difference_time;
	  readout_energy+=0.1*s_runtime_energy;

	  //Store in the res array
	  counter=(accumulated_time+s_T)*100;
	  counter_d=std::modf(counter,&counter_i);
	  c_n=(int)counter_i;
	  res[c_n]+=0.1*s_runtime_energy;

	  //Update remaining energy
	  f_runtime_energy=0.9*f_runtime_energy;
	  s_runtime_energy=0.9*s_runtime_energy;
	}
      
      else
	{
	  //First signal
	  accumulated_time+=(travel_time-difference_time);
	  readout_energy+=0.1*f_runtime_energy;

	  counter=(accumulated_time+s_T)*100;
	  counter_d=std::modf(counter,&counter_i);
	  c_n=(int)counter_i;
	  res[c_n]+=0.1*f_runtime_energy;


	  //Check if readout is enough
	  if(readout_energy>=0.9*s_E)
	    {
	      break;
	    }

	  //Update remaining energy
	  f_runtime_energy=0.9*f_runtime_energy;


	  //Second signal
	  accumulated_time+=difference_time;
	  readout_energy+=0.1*s_runtime_energy;
   
	  counter=(accumulated_time+s_T)*100;
	  counter_d=std::modf(counter,&counter_i);
	  c_n=(int)counter_i;
	  res[c_n]+=0.1*s_runtime_energy;
	  
	  //Check if readout is enough
	  if(readout_energy>=0.9*s_E)
	    {
	      break;
	    }

	  //Update remaining energy
	  s_runtime_energy=0.9*s_runtime_energy;
	}
    }
}

//Drift function calculate the the energy from the pulse at a specific time     
void driftSignal::driftFunction(G4double *pt)
{ 
  G4double temp_E=0;
  G4double temp_T=0;
  G4double n=6100;
  G4double time_interval=61/n;

  for(G4int i=0;i<n;i++)
    {
      temp_T=collision_time+i*time_interval;
      temp_E=(a*(2*temp_T+time_interval)+2*b)*time_interval/2;
      readout(temp_E,temp_T+time_interval,pt);
    }
}


//Get information function is used to achieve the information that user want
void driftSignal::get_info()
{
  std::cout<<"--"<<std::endl;
  std::cout<<"a is: "<<a<<" "<<"b is: "<<b<<std::endl;
  std::cout<<"Energy is: "<<energy<<" [MeV]"<<std::endl;
  std::cout<<"Position is: "<<position<<" [mm]"<<std::endl;
  std::cout<<"Collision time is: "<<collision_time<<" [ns]"<<std::endl;
  std::cout<<"--"<<std::endl;
}
