//////////////////////////////////////////////////////////////////////////////////////
using namespace std;

// 12/09/15 Modified for 24 wires in PC
//2016/01/04 Modified for different X, Y, Z for Silicon & PC
//////////////////////////////////////////////////////////////////////////////////////
#define MaxPCHits 24
#define NPCWires  24
//////////////////////////////////////////////////////////////////////////////////////
class SiHit {
  // This class is for Si hit
 public:
  SiHit(){};
  Int_t NSiHits,DetID[MaxSiHits],UpChNum[MaxSiHits],DownChNum[MaxSiHits],BackChNum[MaxSiHits];
  Int_t HitType[MaxSiHits];
  Double_t X[MaxSiHits],Y[MaxSiHits],Z[MaxSiHits],Energy[MaxSiHits],EnergyBack[MaxSiHits],EnergyUp[MaxSiHits],EnergyDown[MaxSiHits];
  Double_t XW[MaxSiHits],YW[MaxSiHits],ZW[MaxSiHits];
  Double_t RW[MaxSiHits], PhiW[MaxSiHits];//Added by JJPIV 6Jan2015
  Double_t TimeBack[MaxSiHits],TimeUp[MaxSiHits],TimeDown[MaxSiHits];
  Double_t RFSubtract[MaxSiHits];
  Double_t DeltaZ[MaxSiHits];

  Double_t XRand[MaxSiHits],ZRand[MaxSiHits],RRand[MaxSiHits],PhiRand[MaxSiHits];
  
  void zeroSiHit()
  {
    NSiHits = 0;
    for (Int_t i=0; i<MaxSiHits;i++) {
      DetID[i] = -1;
      UpChNum[i] = -1;
      DownChNum[i] = -1;
      BackChNum[i] = -1;
      HitType[i] = 0;
      X[i]  = 0;
      Y[i]  = 0;
      Z[i]  = 0;
      XW[i] = 0;
      YW[i] = 0;
      ZW[i] = 0;
      RW[i] = 0;
      PhiW[i] = 0;
      Energy[i] = 0;
      EnergyBack[i] = 0;
      EnergyUp[i] = 0;
      EnergyDown[i] = 0;
      TimeBack[i] = 0;
      TimeUp[i] = 0;
      TimeDown[i] = 0;
      RFSubtract[i] = 0;
      DeltaZ[i] = 0;
      XRand[i] = 0;
      RRand[i] = 0;
      ZRand[i] = 0;
      PhiRand[i] = 0;
    }
  };
};
//////////////////////////////////////////////////////////////////////////////////////
class PCHit {
  // This class is for Proportional Counter hit
public:
	PCHit(){};
	Int_t NPCHits,WireID[MaxPCHits];
	Int_t Down[MaxPCHits],Up[MaxPCHits];
	Double_t DownVoltage[MaxPCHits],UpVoltage[MaxPCHits],Energy[MaxPCHits];
	Double_t Zp[MaxPCHits];
	Double_t XWp[MaxPCHits],YWp[MaxPCHits],ZWp[MaxPCHits];
	Double_t RWp[MaxPCHits], PhiWp[MaxPCHits];//Added by JJPIV 6Jan2015
	
	Double_t Xrec[MaxPCHits],Yrec[MaxPCHits], Zrec[MaxPCHits];

	void zeroPCHit()
	{
		NPCHits = 0;
		for (Int_t i=0; i<MaxPCHits; i++) {
			WireID[i] = -1;
			Down[i] = 0;
			DownVoltage[i] = 0;
			Up[i]   = 0;
			UpVoltage[i] = 0;
			Energy[i] =0;
			Zp[i] = 0;
			XWp[i] = 0;
			YWp[i] = 0;
			ZWp[i] = 0;
			RWp[i] = 0;
			PhiWp[i] = 0;

			Xrec[i] = 0;
			Yrec[i] = 0; 
			Zrec[i] = 0;

		}
	};
};
/////////////////////////////////////////////////////////////////////////////////////
class Track {
 public:
  Track(){};
  Int_t NTracks;
  Int_t NTracks1;
  Int_t NTracks2;
  Int_t NTracks3;
  Int_t TrackType[MaxTracks];
  Double_t SiZ[MaxTracks];
  Double_t SiR[MaxTracks];
  Double_t SiPhi[MaxTracks];
  Double_t SiEnergy[MaxTracks];
  Double_t PCZ[MaxTracks];
  Double_t PCZ_Ref[MaxTracks];
  Double_t PCR[MaxTracks];
  Double_t PCPhi[MaxTracks];
  Double_t PCEnergy[MaxTracks];
  Double_t PCUpVoltage[MaxTracks];
  Double_t PCDownVoltage[MaxTracks];
  Double_t IntPoint[MaxTracks];
  Double_t Theta[MaxTracks];
  Double_t PathLength[MaxTracks];
  Double_t BeamEnergy[MaxTracks];
  Double_t EnergyLoss[MaxTracks];
  Double_t DiffIntPoint;
  Double_t DiffIntPoint_Proton;
  Double_t SiEnergy_tot;
  Double_t PCEnergy_tot;
  Double_t Energy_tot;
  Double_t Ex;
  Double_t Angle;
  Int_t WireID[MaxTracks];
  Int_t DetID[MaxTracks];
  Int_t ChNum[MaxTracks];
  Double_t DiffPhi[MaxTracks];
  Double_t DiffZ[MaxTracks];
  Double_t DiffR[MaxTracks];
  Double_t BeTheta;
  Double_t BeKE;
  Double_t AvgBeamEnergy;
  Double_t ProtonEnergy;
  Double_t ProtonEnergy_Rec;
  Double_t ProtonDiffEnergy;
  Double_t ProtonAngle;

  void zeroTrack(){
    NTracks = 0;
    NTracks1 = 0;
    NTracks2 = 0;
    NTracks3 = 0;
    DiffIntPoint = -10;
    DiffIntPoint_Proton = -10;
    SiEnergy_tot = 0;
    PCEnergy_tot = 0;
    Energy_tot = 0;
    Ex = -10;
    Angle = -10;
    BeTheta = -10;
    BeKE = -10;
    AvgBeamEnergy = -10;
    ProtonEnergy = -10;
    ProtonEnergy_Rec = -10;
    ProtonDiffEnergy = -1000;
    ProtonAngle = -1000;

    for (Int_t i=0; i<MaxTracks; i++){
      TrackType[i] = -1;
      SiZ[i] = 0;
      SiR[i] = 0;
      SiPhi[i] = 0;
      SiEnergy[i] = 0;
      PCZ[i] = 0;
      PCZ_Ref[i] = 0;
      PCR[i] = 0;
      PCPhi[i] = 0;
      PCEnergy[i] = 0;
      PCUpVoltage[i] = 0;
      PCDownVoltage[i] = 0;
      IntPoint[i] = 0;
      Theta[i] = 0;
      PathLength[i] = 0;
      BeamEnergy[i] = 0;
      EnergyLoss[i] = 0;
      WireID[i] = -1;
      DetID[i] = -1;
      ChNum[i] = -1;
      DiffPhi[i] = -10;
      DiffZ[i] = -10;
      DiffR[i] = -10;
    }

  };

};
/////////////////////////////////////////////////////////////////////////////////////
