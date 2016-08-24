///////////////////////////////////////////////////////////////////////////////////////
//  Header file: EnergyLoss.hpp

//  Description: Simple class that calculates the energy loss of an ion
//  in a gas target. The input is a three-column text file with the
//  energy of the ion and the electrical and nuclear stopping powers
//  (dE/dx) of the target for that ion energy. The units are assumed
//  to be MeV and MeV/mm for the ion's energy and the stopping powers,
//  respectively. This information (E and dE/dx) can be obtained from
//  SRIM. Also, it is assumed that the first line are three strings
//  describing the columns.

//  Author: Daniel Santiago-Gonzalez  //2012-Sep
//  Edited by: Nabin Rijal // 2013-October..
////////////////////////////////////////////////////////////////////////////////////////
  using namespace std;
class EnergyLoss{

public:

  EnergyLoss(){
    c = 29.9792458;           // Speed of light in cm/ns.
    dEdx_e = 0;
    dEdx_n = 0;
    Energy_in_range = 1;
    EvD = new TGraph();
    GoodELossFile = 0;
    IonEnergy = 0;
    IonMass = 0;
    last_point = 0;
    points = 0;
  };
  // Old constructor. Requires special format for the eloss file (two or three columns).
  // In the new version SRIM output files can be read directly.
  ////////////////////////////////////////////////////////////////////////////////////////

  EnergyLoss(string Eloss_file, Float_t IonMass){  /*MeV/c^2*/
    Double_t IonEnergy, dEdx_e, dEdx_n;
    std::string aux;
    Int_t space_pos = 0;


    ifstream Read(Eloss_file.c_str());
    last_point = 0;

    if(!Read.is_open()) {
      cout << "*** EnergyLoss Error: File " << Eloss_file << " was not found." << endl;
      GoodELossFile = 0;
    } 
    else {
      GoodELossFile = 1;        
      // The first line has three strings (columns' description).
      Read >> aux >> aux >> aux;
      // Cout the number of points.
      points = 0;
      
      do{
	Read >> IonEnergy >> dEdx_e >> dEdx_n;
	points++;
      }while(!Read.eof());
      Read.close();
      
      
      // Create the arrays depending on the number rows in the file.
      this->IonEnergy = new Double_t[points];
      this->dEdx_e = new Double_t[points];
      this->dEdx_n = new Double_t[points];    
      
      
      // Go to the begining of the file and read it again to now save the info in the
      // newly created arrays.
      Read.open(Eloss_file.c_str());
      Read >> aux >> aux >> aux;
      
      for(Int_t p=0; p<points; p++){
	Read >> IonEnergy >> dEdx_e >> dEdx_n;
	this->IonEnergy[p] = IonEnergy;
	this->dEdx_e[p] = dEdx_e;
	this->dEdx_n[p] = dEdx_n;
      }    
      
      Energy_in_range = 1;
      this->IonMass = IonMass;  // In MeV/c^2
      c = 29.9792458;           // Speed of light in cm/ns.
      EvD = new TGraph();
      
    }
  };

  Double_t GetEnergyLoss(Float_t initial_energy, Float_t distance);
  void GetEvDCurve(Float_t InitEne, Float_t FinalDepth, Int_t steps);
  Double_t GetInitialEnergy(Float_t FinalEnergy, Float_t PathLength, Float_t StepSize);
  Double_t GetFinalEnergy(Float_t InitialEnergy, Float_t PathLength, Float_t StepSize);  
  Double_t GetPathLength(Float_t InitialEnergy, Float_t FinalEnergy, Float_t DeltaT);
  Double_t GetTimeOfFlight(Float_t InitialEnergy, Float_t PathLength, Float_t StepSize);
  bool LoadSRIMFile(string FileName);
  void SetIonMass(Float_t IonMass);


  bool GoodELossFile;
  TGraph* EvD;

private:

  Double_t c;
  Double_t* IonEnergy;
  Double_t IonMass;
  Double_t* dEdx_e;
  Double_t* dEdx_n;
  Int_t points;
  Int_t last_point;
  bool Energy_in_range;

};

///////////////////////////////////////////////////////////////////////////////////////////
