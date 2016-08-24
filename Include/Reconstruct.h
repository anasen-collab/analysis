//////////////////////////////////////////////////////////////////////////////////////
#include <TROOT.h>
#include <vector>

//#include "/home2/parker/ANASEN/LSU/Include/organizetree.h"
#include "/home/manasta/Desktop/parker_codes/Include/EnergyLoss.h"
using namespace std;

// 12/09/15 Modified for 24 wires in PC
//2016/01/04 Modified for different X, Y, Z for Silicon & PC
//////////////////////////////////////////////////////////////////////////////////////
#define MaxPCHits 24
#define NPCWires  24
#define MaxSiHits 500
#define MaxTracks 100
//////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

#define BeamE 71.816  //Energy of 18Ne beam inside Kapton Window.
#define pcr 3.75    //PC radius.
#define La 53.65    //Length of ANASEN gas volume.
///////////////////Nuclear Masses ///////////////////////////////////////
//nuclear masses //MeV

#define M_Be8 7454.85043438849
#define M_Be7 6534.1836677282 
#define M_D2 1875.61291385342   
#define M_Li5 4667.6163636366931

#define M_Li7 6533.83277448969
#define M_N 939.565413351413   
#define M_He5 4667.67970996292

#define M_Al27 25133.14125
#define M_Na21 20997.655206
#define M_Ne18 16772.20962

#define M_alpha 3727.37929745092
#define M_P     938.27206671856
#define M_Beam  16772.20962    //18Ne
#define M_Heavy 19559.105656 //21Na

/////////////////////////////////////////////////////////////////////////

class Reconstruct {
 public:
  Reconstruct(){};

  void ReconstructHeavy(Track &Tr, Int_t proton);

};


void Reconstruct::ReconstructHeavy(Track &Tr, Int_t proton)

{
  //this reconstructs 21Na from 1 proton
  TLorentzVector P_LV(0.,0.,0.,0.);

  TLorentzVector Beam_LV(0.,0.,0.,0.);
  TLorentzVector Target_LV(0.,0.,0.,0.);
  TLorentzVector Parent_LV(0.,0.,0.,0.); //22Mg

  TLorentzVector Heavy_LV(0.,0.,0.,0.); //21Na

  Float_t E_p_si = 0.0;
  Float_t d_p = 0.0;
  Float_t E_p_loss = 0.0;
  Float_t E_p_rxn = 0.0;
  Float_t theta_p= 0.0;
  Float_t phi_p= 0.0;

  //total energy of proton.
  Float_t E_p_tot =0.0;

  //Total momentum of proton.
  Float_t P_p =0.0;

  //Components of 4-vectors for protons.
  Float_t p_x =0.0;
  Float_t p_y =0.0;
  Float_t p_z =0.0;

  //Kinetic energy and Excitation energy of the 18Ne.

  Float_t Beam_E_tot =0.0;
  Float_t Beam_Pz =0.0;
  Float_t Target_E =0.0;	 

  EnergyLoss *E_Loss_proton = new EnergyLoss("/home/manasta/Desktop/parker_codes/CalParamFiles/H_in_HeCO2(4)_377Torr.txt",M_P);

  E_p_si = Tr.TrEvent.at(proton).SiEnergy;

  d_p = Tr.TrEvent.at(proton).PathLength;

  theta_p = Tr.TrEvent.at(proton).Theta;
  phi_p = Tr.TrEvent.at(proton).SiPhi;
  
  //E_p_loss =E_Loss_proton->GetEnergyLoss(E_p_si,d_p);
  //E_p_rxn = E_p_si + E_p_loss;
  E_p_rxn = E_Loss_proton->GetInitialEnergy(E_p_si,d_p,0.1);
  
  //total energy of protons.
  E_p_tot = E_p_rxn + M_P;
  
  //Total momentum of protons.
  P_p = sqrt(E_p_tot*E_p_tot - M_P*M_P);
  //Components of 4-vectors for proton.
  p_x = P_p*sin(theta_p)*cos(phi_p);
  p_y = P_p*sin(theta_p)*sin(phi_p);
  p_z = P_p*cos(theta_p);

  //Four vectors for proton.
  P_LV.SetPxPyPzE(p_x,p_y,p_z,E_p_tot); 
  
  //Proton reconstruction
  Beam_E_tot = Tr.TrEvent.at(proton).BeamEnergy + M_Beam;
  Beam_Pz = sqrt(Beam_E_tot*Beam_E_tot -M_Beam*M_Beam);
  Target_E = M_alpha;

  //four vectors for beam & Target
  Beam_LV.SetPxPyPzE(0,0,Beam_Pz,Beam_E_tot);
  Target_LV.SetPxPyPzE(0,0,0,Target_E);
  
  //four vectors for parents 
  Parent_LV= Beam_LV + Target_LV;

  //four vectors for Sodium 21Na
  Heavy_LV = Parent_LV-P_LV;

  Tr.AlEvent.IntPoint = Tr.TrEvent.at(proton).IntPoint;
  Tr.AlEvent.BeamEnergy = Tr.TrEvent.at(proton).BeamEnergy;
  Tr.AlEvent.SiEnergy_tot = Tr.TrEvent.at(proton).SiEnergy;
  // Tr.AlEvent.PCEnergy_tot = Tr.TrEvent.at(proton).PCEnergy*Tr.TrEvent.at(proton).PathLength;
  Tr.AlEvent.PCEnergy_tot = Tr.TrEvent.at(proton).PCEnergy*Tr.TrEvent.at(proton).Theta;
  //Tr.AlEvent.Energy_tot = Tr.AlEvent.SiEnergy_tot + Tr.AlEvent.PCEnergy_tot; //kinetic energy of the proton
  Tr.AlEvent.Energy_tot = E_p_rxn;

  Tr.AlEvent.Theta = Heavy_LV.Theta()*180/TMath::Pi();
  Tr.AlEvent.Phi = Heavy_LV.Phi()*180/TMath::Pi();
  //Kinetic energy and Excitation energy of the 8Be.

  Tr.AlEvent.KE = Heavy_LV.E()-Heavy_LV.M(); // .M()=sqrt(E^2-p^2) .E()=4th component of 4-vector, Energy
  Tr.AlEvent.Ex = Heavy_LV.M() - M_Heavy;
  //Tr.AlEvent.Ex = M_Heavy - Heavy_LV.M();
  
 

}
