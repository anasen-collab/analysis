////////////////////////////////////////////////////////////////////////////////////////
//   Code: EnergyLoss.cpp

//   Description: Simple class that calculates the energy loss of an ion
//   in a gas target. The input is a three-column text file with the
//   energy of the ion and the electrical and nuclear stopping powers
//   (dE/dx) of the target for that ion energy. The units are assumed
//   to be MeV and MeV/mm for the ion's energy and the stopping powers,
//   respectively. This information (E and dE/dx) can be obtained from
//   SRIM. Also, it is assumed that the first line are three strings
//   describing the columns.

//  Author: Daniel Santiago-Gonzalez //2012-Sep.
//  Edited by: Nabin Rijal // 2013-October..

////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <string.h>
#include <TGraph.h>

#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <cmath>

#include "/home/manasta/Desktop/parker_codes/Include/EnergyLoss.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////
// Get the energy loss of the ion in the gas for a given ion's energy and a differential
// distance through the gas target. Most important function in this class. 
////////////////////////////////////////////////////////////////////////////////////////


Double_t EnergyLoss::GetEnergyLoss(Float_t energy /*MeV*/, Float_t distance /*cm*/)
{
  Int_t i = -1;
  // Look for two points for which the initial energy lays in between.
  // This for-loop should find the points unless there was a big jump from
  // the energy used in the last point and the energy used now.
  for(Int_t p=last_point-1; p<points-1; p++){
    if (last_point>=0)
      if(energy>=IonEnergy[p]  && energy<IonEnergy[p+1]){
	i = p+1;
	last_point = p;
	break;
      }
  }
  // It is probable that if the point wasn't found could have been because of
  // a big jump in the energy (see above), so we need to look in the remaining
  // points.
  if (i==-1) {
    for(Int_t p=0; p<last_point-1; p++){
      if(energy>=IonEnergy[p]  && energy<IonEnergy[p+1]){
	i = p+1;
	last_point = p;
	break;
      }
    }
  }
  // If after the last two for-loops 'i' is still -1 it means the energy was out of range.
  if(i==-1){
    cout << "*** EnergyLoss Error: energy not within range: " << energy << endl;
    Energy_in_range = 0;
    return 0;
  }

  // If the initial energy is within the range of the function, get the stopping power
  // for the initial energy.  
  Double_t E1 = IonEnergy[i-1];        Double_t E2 = IonEnergy[i];
  Double_t dEdx_e1 = dEdx_e[i-1];      Double_t dEdx_e2 = dEdx_e[i];  
  Double_t dEdx_n1 = dEdx_n[i-1];      Double_t dEdx_n2 = dEdx_n[i];  


  // Interpolating the electric stopping power (from point 1 to 'e').
  Double_t dEdx_e1e = dEdx_e1 + (energy - E1)*(dEdx_e2 - dEdx_e1)/(E2 - E1);


  // Interpolating the nuclear stopping power (usually negligable).
  Double_t dEdx_n1e = dEdx_n1 + (energy - E1)*(dEdx_n2 - dEdx_n1)/(E2 - E1);


  // The stopping power units are in MeV/mm so we multiply by 10 to convert to MeV/cm.
  return ( (dEdx_e1e+dEdx_n1e)*10*distance );  
}


////////////////////////////////////////////////////////////////////////////////////////
// This function sets the points of the energy v. distance TGraph object (called EvD)
// for a given initial energy, going from distance 0 to a final distance (depth) in a
// given number of steps.  The EvD object can then be diplayed by 
//   EL->EvD->Draw("l");
// assuming EL is an EnergyLoss pointer.
////////////////////////////////////////////////////////////////////////////////////////

void EnergyLoss::GetEvDCurve(Float_t InitEne/*MeV*/, Float_t FinalDepth/*cm*/, Int_t steps)
{
  Double_t current_ene, current_depth;
  current_ene = InitEne;
  current_depth = 0;
  for(Int_t s=0; s<steps; s++){
    current_ene -= GetEnergyLoss(current_ene,FinalDepth/steps);
    current_depth += FinalDepth/steps;
    EvD->SetPoint(s, current_depth, current_ene);
  }
}


////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////

Double_t EnergyLoss::GetInitialEnergy(Float_t FinalEnergy /*MeV*/, Float_t PathLength /*cm*/, Float_t StepSize/*cm*/)
{
  Double_t Energy = FinalEnergy;
  Int_t Steps = (int)floor(PathLength/StepSize);
  last_point = 0;
  // The function starts by assuming FinalEnergy is within the energy range, but
  // this could be changes in the GetEnergyLoss() function.
  Energy_in_range = 1;

  for (Int_t s=0; s<Steps; s++) {
    Energy = Energy + GetEnergyLoss(Energy,PathLength/Steps);
    if (!Energy_in_range)
      break;
  } 
  Energy = Energy + GetEnergyLoss(Energy,PathLength-Steps*StepSize);
  if (!Energy_in_range)
    Energy = -1000; // Return an unrealistic value.
  
  //  cout << "d: K_lf=" << FinalEnergy << "  K_lr=" << Energy << "  l=" << PathLength <<endl;

  return Energy;
}


////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////

Double_t EnergyLoss::GetFinalEnergy(Float_t InitialEnergy /*MeV*/, Float_t PathLength /*cm*/, Float_t StepSize/*cm*/)
{

  Double_t Energy = InitialEnergy;
  Int_t Steps = (int)floor(PathLength/StepSize);

  // The function starts by assuming InitialEnergy is within the energy range, but
  // this could be changes in the GetEnergyLoss() function.

  Energy_in_range = 1;

  for (Int_t s=0; s<Steps; s++) {
    Energy = Energy - GetEnergyLoss(Energy,PathLength/Steps);
    if (!Energy_in_range)
      break;
  }  

  Energy = Energy - GetEnergyLoss(Energy,PathLength-Steps*StepSize);

  if (!Energy_in_range) 
    Energy = -1000;
  //  cout << "O: K_bw=" << InitialEnergy << "  K_br=" << Energy << "  l=" << PathLength <<endl;
  return Energy;

}


////////////////////////////////////////////////////////////////////////////////////////
// Calulates the ion's path length in cm.
////////////////////////////////////////////////////////////////////////////////////////

Double_t EnergyLoss::GetPathLength(Float_t InitialEnergy /*MeV*/, Float_t FinalEnergy /*MeV*/, Float_t DeltaT /*ns*/)
{
  Double_t L = 0, DeltaX = 0;
  Double_t Kn = InitialEnergy;
  Int_t n=0;
  if (IonMass==0)
    cout << "*** EnergyLoss Error: Path length cannot be calculated for IonMass = 0." << endl;
  else {
    // The path length (L) is proportional to sqrt(Kn). After the sum, L will be multiplied by
    // the proportionality factor.
    while (Kn>FinalEnergy && n<(int)pow(10.0,6)) {
      L += sqrt(Kn);                       // DeltaL going from point n to n+1.
      DeltaX = sqrt(2*Kn/IonMass)*DeltaT*c;// dx = v*dt
      Kn -= GetEnergyLoss(Kn, DeltaX);     // After L is incremented the kinetic energy at n+1 is calculated.
      n++;    
    }
    if (n>=(int)pow(10.0,6)) {
      cout << "*** EnergyLoss Warning: Full path length wasn't reached after 10^6 iterations." << endl;
      L = 0;
    }
    else
      L *= sqrt(2/IonMass)*DeltaT*c;
  }
  return L;
}



////////////////////////////////////////////////////////////////////////////////////////
// Calulates the ion's time of flight in ns.
////////////////////////////////////////////////////////////////////////////////////////

Double_t EnergyLoss::GetTimeOfFlight(Float_t InitialEnergy /*MeV*/, Float_t PathLength /*cm*/, Float_t StepSize /*cm*/)
{
  Double_t TOF = 0;
  Double_t Kn = InitialEnergy;
  Int_t Steps = (Int_t)PathLength/(Int_t)StepSize;
  if (IonMass==0)
    cout << "*** EnergyLoss Error: Time of flight cannot be calculated for IonMass = 0." << endl;
  else {
    // The TOF is proportional to 1/sqrt(Kn). After the sum, TOF will be multiplied by
    // the proportionality factor.
    for (Int_t n=0; n<Steps; n++) {
      TOF += 1/sqrt(Kn);                 // DeltaT going from point n to n+1.
      Kn -= GetEnergyLoss(Kn, StepSize); // After the TOF is added the kinetic energy at point n+1 is calc.
    }
    TOF *= sqrt(IonMass/2)*StepSize/c;
  }
  return TOF;
}



////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////

bool EnergyLoss::LoadSRIMFile(string FileName)
{

  Double_t IonEnergy, dEdx_e, dEdx_n;
  string aux, unit;
  Int_t space_pos = 0, str_counter=0;


  ifstream Read(FileName.c_str());
  last_point = 0;


  if(!Read.is_open()) {
    cout << "*** EnergyLoss Error: File " << FileName << " was not found." << endl;
    GoodELossFile = 0;
  } 

  else {
    GoodELossFile = 1;
    Energy_in_range = 1;        
    // Read all the string until you find "Straggling", then read the next 7 strings.
    // (this method is not elegant at all but there is no time to make it better)
    do 
      Read >> aux;
    while (aux!="Straggling");
    cout << "Found!" << endl;
    for (Int_t i=0; i<7; i++)
      Read >> aux;

    do {
      Read >> aux;
      str_counter++;
    } while (aux!="-----------------------------------------------------------");
    Read.close();
    
    cout << str_counter << endl;
    str_counter--;
    points = str_counter/10;
    cout << points << endl;

    // Create the arrays depending on the number rows in the file.
    this->IonEnergy = new Double_t[points];
    this->dEdx_e = new Double_t[points];
    this->dEdx_n = new Double_t[points]; 

   
    // Go to the begining of the file and read it again to now save the info in the
    // newly created arrays.
    Read.open(FileName.c_str());
    do 
      Read >> aux;
    while (aux!="Straggling");
    cout << "Found!" << endl;
    for (Int_t i=0; i<7; i++)
      Read >> aux;
    

    for (Int_t p=0; p<points; p++) {
      Read >> IonEnergy >> unit >> dEdx_e >> dEdx_n >> aux >> aux >> aux >> aux >> aux >> aux;


      if (unit=="keV")
	IonEnergy *= 0.001;
      //  cout << p << " " << IonEnergy << " " << unit << " " << dEdx_e << " " << dEdx_n << endl;


      this->IonEnergy[p] = IonEnergy;
      this->dEdx_e[p] = dEdx_e*0.008752; // !!!!!
      this->dEdx_n[p] = dEdx_n*0.008752; // !!!!!
    }    
  }

  cout << " You have a modulation on the energies!!!" << endl;
  return GoodELossFile;
}

//////////////////////////////////////////////////////////////////////////////////

void EnergyLoss::SetIonMass(Float_t IonMass)
{
  this->IonMass = IonMass;
  return;  
}

//////////////////////////////////////////////////////////////////////////////////
