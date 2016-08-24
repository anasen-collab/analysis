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

#define BeamE 19.6  //Energy of 7Be beam inside Kapton Window.
#define pcr 3.75    //PC radius.
#define La 53.65    //Length of ANASEN gas volume.
///////////////////Nuclear Masses ///////////////////////////////////////
//nuclear masses //MeV
#define M_alpha 3727.37929745092
#define M_Be8 7454.85043438849
#define M_Be7 6534.1836677282 
#define M_D2 1875.61291385342
#define M_P 938.27206671856      
#define M_Li5 4667.6163636366931

#define M_Li7 6533.83277448969
#define M_N 939.565413351413   
#define M_He5 4667.67970996292

#define M_Al27 25133.14125
#define M_Na21 20997.655206
#define M_Ne18 16772.20962
/////////////////////////////////////////////////////////////////////////

class Reconstruct {
 public:
  Reconstruct(){};

  void ReconstructBe(Track &Tr, Int_t alpha1, Int_t alpha2, Int_t proton);
  void ReconstructAl(Track &Tr, Int_t proton);

};

void Reconstruct::ReconstructBe(Track &Tr, Int_t alpha1, Int_t alpha2, Int_t proton){
  //this reconstructs 8Be from 2 alphas
  //if we have two tracks, proton should equal -1
  //if we have three tracks (two alpha), proton will equal the location of the proton track in the array (0, 1, or 2)
  TLorentzVector A1_LV(0.,0.,0.,0.);
  TLorentzVector A2_LV(0.,0.,0.,0.);
  TLorentzVector Be8_LV(0.,0.,0.,0.);

  TLorentzVector Beam_LV(0.,0.,0.,0.);
  TLorentzVector Target_LV(0.,0.,0.,0.);
  TLorentzVector Parent_LV(0.,0.,0.,0.);
  TLorentzVector Proton_LV(0.,0.,0.,0.);
  TLorentzVector Proton_LV_Tr(0.,0.,0.,0.);

  TLorentzVector Exe_LV(0.,0.,0.,0.);

  Float_t E_a1_si = 0.0;
  Float_t E_a2_si = 0.0;
  Float_t d_a1 = 0.0;
  Float_t d_a2 = 0.0;
  Float_t E_a1_loss = 0.0;
  Float_t E_a2_loss = 0.0; 
  Float_t E_a1_rxn = 0.0;
  Float_t E_a2_rxn = 0.0;
  Float_t E_8Be = 0.0;
  Float_t theta_a1= 0.0;
  Float_t theta_a2= 0.0;
  Float_t phi_a1= 0.0;
  Float_t phi_a2= 0.0;

  //total energy of alphas.
  Float_t E_a1_tot =0.0;
  Float_t E_a2_tot =0.0;

  //Total momentum of alphas.
  Float_t P_a1 =0.0;
  Float_t P_a2 =0.0;

  //Components of 4-vectors for alphas.
  Float_t a1x =0.0;
  Float_t a1y =0.0;
  Float_t a1z =0.0;

  Float_t a2x =0.0;
  Float_t a2y =0.0;
  Float_t a2z =0.0;	
	  
  //components of 4-vectors of 8Be.
	  
  Float_t Bex =0.0;
  Float_t Bey =0.0;
  Float_t Bez =0.0;
  Float_t Be_E =0.0;

  Float_t theta_8Be=0.0;
  Float_t phi_8Be=0.0;
      
  //Kinetic energy and Excitation energy of the 8Be. 
  Float_t Be_KE=0.0;
  Float_t Be_exe =0.0;

  //////////////////////Proton Reconstruction from 2 alphas ////////

  Float_t Beam_E_tot =0.0;
  Float_t Beam_Pz =0.0;
  Float_t Target_E =0.0;	 

  Float_t Proton_Px =0.0;
  Float_t Proton_Py =0.0;
  Float_t Proton_Pz =0.0;
  Float_t Proton_E =0.0;	 
  Float_t M_P1	 =0.0;

  Float_t Proton_theta =0.0;
  Float_t Proton_theta_deg =0.0;
  Float_t Proton_phi =0.0;
  Float_t Proton_phi_deg =0.0;
  Float_t Proton_KE =0.0;		 
  Float_t Proton_KE1=0.0;

  EnergyLoss *E_Loss_alpha = new EnergyLoss("/home2/parker/ANASEN/LSU/CalParamFiles/He4_D2_400Torr.eloss",M_alpha);
  EnergyLoss *E_Loss_proton = new EnergyLoss("/home2/parker/ANASEN/LSU/CalParamFiles/p_D2_400Torr.eloss",M_P);

  Tr.BeEvent.DiffIntPoint = Tr.TrEvent.at(alpha1).IntPoint - Tr.TrEvent.at(alpha2).IntPoint;
  //4 vector way
  Tr.BeEvent.SiEnergy_tot = Tr.TrEvent.at(alpha1).SiEnergy + Tr.TrEvent.at(alpha2).SiEnergy;
  Tr.BeEvent.PCEnergy_tot = Tr.TrEvent.at(alpha1).PCEnergy*Tr.TrEvent.at(alpha1).PathLength + Tr.TrEvent.at(alpha2).PCEnergy*Tr.TrEvent.at(alpha2).PathLength;
  Tr.BeEvent.Energy_tot = Tr.BeEvent.SiEnergy_tot + Tr.BeEvent.PCEnergy_tot;

  E_a1_si = Tr.TrEvent.at(alpha1).SiEnergy;
  E_a2_si = Tr.TrEvent.at(alpha2).SiEnergy;

  d_a1 = Tr.TrEvent.at(alpha1).PathLength;
  d_a2 = Tr.TrEvent.at(alpha2).PathLength;

  theta_a1 = Tr.TrEvent.at(alpha1).Theta;
  theta_a2 = Tr.TrEvent.at(alpha2).Theta;
  phi_a1 = Tr.TrEvent.at(alpha1).SiPhi;
  phi_a2 = Tr.TrEvent.at(alpha2).SiPhi;
  
  E_a1_loss =E_Loss_alpha->GetEnergyLoss(E_a1_si,d_a1);
  E_a2_loss =E_Loss_alpha->GetEnergyLoss(E_a2_si,d_a2);
  
  E_a1_rxn = E_a1_si + E_a1_loss;
  E_a2_rxn = E_a2_si + E_a2_loss;
  
  //total energy of alphas.
  E_a1_tot = E_a1_rxn + M_alpha;
  E_a2_tot = E_a2_rxn + M_alpha;
  
  //Total momentum of alphas.
  P_a1 = sqrt(E_a1_tot*E_a1_tot - M_alpha*M_alpha);
  P_a2 = sqrt(E_a2_tot*E_a2_tot - M_alpha*M_alpha);
  
  //Components of 4-vectors for alphas.
  a1x = P_a1*sin(theta_a1)*cos(phi_a1);
  a1y = P_a1*sin(theta_a1)*sin(phi_a1);
  a1z = P_a1*cos(theta_a1);
  
  a2x = P_a2*sin(theta_a2)*cos(phi_a2);
  a2y = P_a2*sin(theta_a2)*sin(phi_a2);
  a2z = P_a2*cos(theta_a2);
  
  //Four vectors for alphas.
  A1_LV.SetPxPyPzE(a1x,a1y,a1z,E_a1_tot);
  A2_LV.SetPxPyPzE(a2x,a2y,a2z,E_a2_tot);
  
  Tr.BeEvent.AlphaAngle = A1_LV.Angle(A2_LV.Vect());
  
  //Four vectors for 8Be
  Be8_LV=A1_LV+A2_LV;   
  
  //accessing components of 4-vectors of 8Be.
  Bex = Be8_LV.Px();
  Bey = Be8_LV.Py();
  Bez = Be8_LV.Pz();
  Be_E = Be8_LV.E();
  
  theta_8Be=Be8_LV.Theta()*180/TMath::Pi();//Theta angle of the 8Be
  Tr.BeEvent.BeTheta = theta_8Be;
  phi_8Be=Be8_LV.Phi()*180/TMath::Pi();//Phi angle of the 8Be
  Tr.BeEvent.BePhi = phi_8Be;
  //Kinetic energy and Excitation energy of the 8Be.

  Be_KE=Be_E-Be8_LV.M();
  Be_exe = Be8_LV.M() - M_Be8;

  Tr.BeEvent.BeKE = Be_E - M_Be8;
  Tr.BeEvent.Ex = Be_exe;
  Tr.BeEvent.AvgBeamEnergy = (Tr.TrEvent.at(alpha1).BeamEnergy + Tr.TrEvent.at(alpha2).BeamEnergy)/2;

  
  //Proton reconstruction
  Beam_E_tot = Tr.BeEvent.AvgBeamEnergy + M_Be7;
  Beam_Pz = sqrt(Beam_E_tot*Beam_E_tot -M_Be7*M_Be7);
  Target_E = M_D2;

  //four vectors for beam & Target
  Beam_LV.SetPxPyPzE(0,0,Beam_Pz,Beam_E_tot);
  Target_LV.SetPxPyPzE(0,0,0,M_D2);
  
  //four vectors for parents 
  Parent_LV= Beam_LV + Target_LV;

  //four vectors for Proton
  Proton_LV_Tr = Parent_LV-Be8_LV;//this is the 4 vector for the expected proton path
  Tr.BeEvent.ProtonEnergy_Rec = Proton_LV_Tr.E() - M_P;
  if ( proton != -1 ){
    //if we have a third track that is not an alpha, assume it is a proton
    Double_t E_p1_si = Tr.TrEvent.at(proton).SiEnergy;
    Double_t d_p1 = Tr.TrEvent.at(proton).PathLength;
    Double_t theta_p1 = Tr.TrEvent.at(proton).Theta;
    Double_t phi_p1 = Tr.TrEvent.at(proton).SiPhi;

    Double_t E_p1_loss = E_Loss_proton->GetEnergyLoss(E_p1_si,d_p1);
    Double_t E_p1_rxn = E_p1_si + E_p1_loss;

    Double_t E_p1_tot = E_p1_rxn + M_P;//+3 makes E-E_Tr centered about zero (ish)
    Double_t P_p1 = sqrt(E_p1_tot*E_p1_tot - M_P*M_P);

    Double_t p1x = P_p1*sin(theta_p1)*cos(phi_p1);
    Double_t p1y = P_p1*sin(theta_p1)*sin(phi_p1);
    Double_t p1z = P_p1*cos(theta_p1);

    Proton_LV.SetPxPyPzE(p1x,p1y,p1z,E_p1_tot);//4-vector for actual proton track

    Tr.BeEvent.ProtonEnergy = Proton_LV.E() - M_P;
    Tr.BeEvent.ProtonDiffEnergy = Proton_LV.E() - Proton_LV_Tr.E();
    Tr.BeEvent.ProtonAngle = Proton_LV.Angle(Proton_LV_Tr.Vect())/TMath::Pi()*180;

  }

}
void Reconstruct::ReconstructAl(Track &Tr, Int_t proton){
  //this reconstructs 21Na from 1 proton
  TLorentzVector P_LV(0.,0.,0.,0.);

  TLorentzVector Beam_LV(0.,0.,0.,0.);
  TLorentzVector Target_LV(0.,0.,0.,0.);
  TLorentzVector Parent_LV(0.,0.,0.,0.);

  TLorentzVector Ne_LV(0.,0.,0.,0.);

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

  EnergyLoss *E_Loss_proton = new EnergyLoss("/home/manasta/Desktop/parker_codes/CalParamFiles/p_alpha_300Torr.eloss",M_P);

  E_p_si = Tr.TrEvent.at(proton).SiEnergy;

  d_p = Tr.TrEvent.at(proton).PathLength;

  theta_p = Tr.TrEvent.at(proton).Theta;
  phi_p = Tr.TrEvent.at(proton).SiPhi;
  
  E_p_loss =E_Loss_proton->GetEnergyLoss(E_p_si,d_p);
  E_p_rxn = E_p_si + E_p_loss;
  
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
  Beam_E_tot = Tr.TrEvent.at(proton).BeamEnergy + M_Ne18;
  Beam_Pz = sqrt(Beam_E_tot*Beam_E_tot -M_Ne18*M_Ne18);
  Target_E = M_alpha;

  //four vectors for beam & Target
  Beam_LV.SetPxPyPzE(0,0,Beam_Pz,Beam_E_tot);
  Target_LV.SetPxPyPzE(0,0,0,Target_E);
  
  //four vectors for parents 
  Parent_LV= Beam_LV + Target_LV;

  //four vectors for Sodium 21Na
  Ne_LV = Parent_LV-P_LV;

  Tr.AlEvent.IntPoint = Tr.TrEvent.at(proton).IntPoint;
  Tr.AlEvent.BeamEnergy = Tr.TrEvent.at(proton).BeamEnergy;
  Tr.AlEvent.SiEnergy_tot = Tr.TrEvent.at(proton).SiEnergy;
  // Tr.AlEvent.PCEnergy_tot = Tr.TrEvent.at(proton).PCEnergy*Tr.TrEvent.at(proton).PathLength;
  Tr.AlEvent.PCEnergy_tot = Tr.TrEvent.at(proton).PCEnergy*Tr.TrEvent.at(proton).Theta;
  //Tr.AlEvent.Energy_tot = Tr.AlEvent.SiEnergy_tot + Tr.AlEvent.PCEnergy_tot; //kinetic energy of the proton
  Tr.AlEvent.Energy_tot = E_p_rxn;

  Tr.AlEvent.Theta = Al_LV.Theta()*180/TMath::Pi();
  Tr.AlEvent.Phi = Al_LV.Phi()*180/TMath::Pi();
  //Kinetic energy and Excitation energy of the 8Be.

  Tr.AlEvent.KE = Al_LV.E()-Al_LV.M(); // .M()=sqrt(E^2-p^2) .E()=4th component of 4-vector, Energy
  Tr.AlEvent.Ex = Al_LV.M() - M_Al27;
 

}
