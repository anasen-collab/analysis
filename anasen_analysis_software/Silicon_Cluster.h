////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Description: Sorts silicon signals and assigns them as particle hit and makes Cluster of hits.
//
// Different Methods to sort QQQ3 && SX3 detectors
//
// Speciality: Multiple particles per detectors are sorted & Front-Back matching.
// 
//
// Usage: Include in Main & call for SortQQQ() && SortSX3().
//
// Author: Nabin Rijal, 2016 August
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <TROOT.h>
#include <TRandom3.h>
#include <TMath.h>
#include <algorithm>
#include <vector>

#include "ChannelMap.h"
#include "tree_structure.h"

using namespace std;

#define Q3_Qdiff 0.05    // 0.15
#define SX3_Qdiff 0.2  // 0.5
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Silicon_Cluster{

 public:
  Silicon_Cluster(){  

  };

  int Initialize(){
    Random = new TRandom3();
    Random->SetSeed();  
  };

  struct data {
    Double_t Energy;
    Double_t Time;
    Int_t Channel;
    Int_t Channel_Up;
    Int_t Channel_Down;
    Int_t Channel_Front;
    Double_t Energy_Up;
    Double_t Energy_Down;
    Double_t Energy_Front;
  };

  TRandom3 *Random;

  data data_obj;
  vector<data> front;
  vector<data> back;

  Double_t QQQR;
  Double_t QQQPhi;

  //Declaring Variables to access SX3 Z info to the Main for Filling Histos.  
  static bool Esort_method(struct data a,struct data b);
  static bool Csort_method(struct data a,struct data b);
 
  void SortQ3(SiHit *Si, ChannelMap *CMAP);
  void SortSX3(SiHit *Si, ChannelMap *CMAP);
  
  ~Silicon_Cluster(){   
    delete Random;
    cout<<"========Random deleted==========="<<endl;
  };  
  
};
/////////////////////////////// Sorting by Energy or Channels ///////////////////////////////////////////////

bool Silicon_Cluster::Esort_method(struct data a,struct data b){
  if(a.Energy>b.Energy)
    return 1;
  return 0;
};

bool Silicon_Cluster::Csort_method(struct data a,struct data b){
  if(a.Channel<b.Channel)
    return 1;
  return 0;
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Silicon_Cluster::SortQ3(SiHit *Si, ChannelMap *CMAP){

  const Double_t InnerRadius = 5.01; //cm
  const Double_t OuterRadius = 10.1; //cm 
  
  const Double_t RingPitch = (OuterRadius-InnerRadius)/16;
  //const Double_t ArcLength = 0.96; //cm (distance between strips at outer radius) 

  //////// Calculating strip width at outer radius /////
  Double_t Gap =2.863636; //in degrees  //gap between Q3 from center
  Double_t Nab =((360-4*Gap)/360)*2*TMath::Pi()/(4*16); //No. of Q3 ==4; No. of Channels per Q3 ==16;
  Double_t ArcLength =Nab*OuterRadius;
  Double_t StripAngle = 2*TMath::ASin(ArcLength/(2*OuterRadius)); // angle spanned by a single strip  
  

  front.clear();
  back.clear();
 
  //copy contents of det_obj into new arrays, 
  for ( Int_t i=0; i<Si->det_obj.FrontMult; i++ ){
    data_obj.Channel = Si->det_obj.FrontChNum[i];
    data_obj.Energy = Si->det_obj.EFront_Cal[i];
    data_obj.Time = Si->det_obj.TFront[i];
    if (data_obj.Energy>0){//if the energy is less than or equal to zero,discard it
      front.push_back(data_obj);
    }
  }
  for ( Int_t i=0; i<Si->det_obj.BackMult; i++ ){
    data_obj.Channel = Si->det_obj.BackChNum[i];
    data_obj.Energy = Si->det_obj.EBack_Cal[i];
    data_obj.Time = Si->det_obj.TBack[i];
    if (data_obj.Energy>0){//if the energy is less than or equal to zero, discard it
      back.push_back(data_obj);
    }
  }
  if (front.size()==0 || back.size()==0){
    return;
  }  
  //---------------------------------------------------------------------------------
  //sort energies, largest to smallest
  sort(front.begin(), front.end(), Esort_method );
  sort(back.begin(), back.end(), Esort_method );

  ////////////////////////////////////////////////////////////////////////////////////////////////////////// 
  Float_t fEnergy[16] = {0};
  Float_t bEnergy[16] = {0};
  Int_t fChannel[16] = {0};
  Int_t bChannel[16] = {0};
  Float_t bTime[16]= {0};
  int i=0,j=0,k=0;

  while ((i< front.size()) && (j< back.size()) && (k<16)){   

    if (fabs(back[j].Energy - front[i].Energy) < Q3_Qdiff*back[j].Energy + 0.2 ){
    
      fEnergy[k] = front[i].Energy;
      bEnergy[k] = back[j].Energy;
      fChannel[k] = front[i].Channel;
      bChannel[k] = back[j].Channel;
      bTime[k] = back[j].Time;		
      i++;
      j++;
      k++;	

    }
    else{      
      // test for two rings / one back-wedge
      if (front[i].Energy < back[j].Energy){
	// is there another front-hit that matches as a sum 
	if (((i+1)< front.size()) && (fabs(front[i].Energy + front[i+1].Energy - back[j].Energy)< Q3_Qdiff*back[j].Energy+0.2 )){

	  // first hit 	
	  fEnergy[k] = front[i].Energy;
	  bEnergy[k] = front[i].Energy;
	  fChannel[k] = front[i].Channel;
	  bChannel[k] = back[j].Channel;
	  bTime[k] = back[j].Time;	
	  k++;
	  i++;

	  // second hit 
	  fEnergy[k] = front[i].Energy;
	  bEnergy[k] = front[i].Energy;
	  fChannel[k] = front[i].Channel;
	  bChannel[k] = back[j].Channel;
	  bTime[k] = back[j].Time;	
	  i++;
	  j++;
	  k++;	  

	}else{
	  break;
	}
      }     
      else if (front[i].Energy > back[j].Energy){
	// is there another back-hit that matches as a sum 
	if (((j+1)< back.size()) && (fabs(back[j].Energy + back[j+1].Energy - front[i].Energy)< Q3_Qdiff*front[i].Energy+0.2 )){

	  // first hit //in the back
	  fEnergy[k] = back[j].Energy;
	  bEnergy[k] = back[j].Energy;
	  fChannel[k] = front[i].Channel;
	  bChannel[k] = back[j].Channel;
	  bTime[k] = back[j].Time;	
	  k++;
	  j++;

	  // second hit //in the back
	  fEnergy[k] = back[j].Energy;
	  bEnergy[k] = back[j].Energy;
	  fChannel[k] = front[i].Channel;
	  bChannel[k] = back[j].Channel;
	  bTime[k] = back[j].Time;	
	  i++;
	  j++;
	  k++;

	}else{
	  break;
	}
      }else{
	// this might have been noise, keep what we already had 
	break;
      }
    }
     
  }//end of while 
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //-----------------------------------------------------------------------------------------
  for(int r=0;r<k;r++){
  
    //determine radius and angle
    QQQR = OuterRadius - (fChannel[r]+Random->Rndm() )*RingPitch;
    QQQPhi = (bChannel[r] + Random->Rndm() )*StripAngle;
   
    //fill tree
    Si->hit_obj.NHitsInDet = k;  
    Si->hit_obj.DetID = Si->det_obj.DetID;
    Si->hit_obj.HitType = back.size()*10 + front.size();
    //cout<<"Q3  DetNum = "<<Si->hit_obj.DetID<<"  HitinDet = "<<m<<"  HitType  = "<<Si->hit_obj.HitType<<endl;

    Si->hit_obj.FrontChannel = fChannel[r];
    Si->hit_obj.BackChannel = bChannel[r];
    Si->hit_obj.EnergyFront = fEnergy[r];
    Si->hit_obj.EnergyBack = bEnergy[r];
    Si->hit_obj.Energy = bEnergy[r];
    Si->hit_obj.Time = bTime[r]; 
 
    Si->hit_obj.X = (QQQR*TMath::Cos(QQQPhi));
    Si->hit_obj.Y = (QQQR*TMath::Sin(QQQPhi));

    Double_t xw=0, yw=0, rw=0, phiw=0;
    
    CMAP->GetQ3WorldCoordinates(Si->hit_obj.DetID,Si->hit_obj.X,Si->hit_obj.Y,xw,yw,rw,phiw);
    Si->hit_obj.XW = (xw);
    Si->hit_obj.YW = (yw);
    Si->hit_obj.RW = (rw);
    Si->hit_obj.PhiW = (phiw);
  
    //Since top 2 Q3 are at zero of our Z-reference and bottom 2 Q3 are  Upstream by 0.762 cm
    //Z-positions for Q3 are set accordingly.

    if ( Si->hit_obj.DetID==1 || Si->hit_obj.DetID==2 ){
      Si->hit_obj.Z = (0.762);
      Si->hit_obj.ZW = (0.762);
    }else if (Si->hit_obj.DetID==0 || Si->hit_obj.DetID==3){
      Si->hit_obj.Z = (0);
      Si->hit_obj.ZW = (0);
    }else{
      //cout<<" Who is this ?"<<endl;   
    }
    
    Si->Hit.push_back(Si->hit_obj);   
    Si->NSiHits++;
  } 
};
//////////////////////////////////////////////////////////////////////////////////////////////////
void Silicon_Cluster::SortSX3(SiHit *Si, ChannelMap *CMAP){

  front.clear();
  back.clear();  
  //------------------------------------------------------------------
  Double_t up[4][3];
  Double_t down[4][3];
  ////Initialization can be done all -1 //since channels are 0,1,2,3
  for ( Int_t i=0; i<4; i++ ){
    for( Int_t j=0; j<3; j++ ){
      up[i][j] = -1;
      down[i][j] = -1;
    }
  }   
  
  for ( Int_t i=0; i<Si->det_obj.UpMult; i++ ){
    up[Si->det_obj.UpChNum[i]][0] = Si->det_obj.UpChNum[i];
    up[Si->det_obj.UpChNum[i]][1] = Si->det_obj.EUp_Cal[i];
    up[Si->det_obj.UpChNum[i]][2] = Si->det_obj.TUp[i];
  }
  
  for ( Int_t i=0; i<Si->det_obj.DownMult; i++ ){
    down[Si->det_obj.DownChNum[i]][0] = Si->det_obj.DownChNum[i];
    down[Si->det_obj.DownChNum[i]][1] = Si->det_obj.EDown_Cal[i];
    down[Si->det_obj.DownChNum[i]][2] = Si->det_obj.TDown[i];
  }
  //------------------------------------------------------------------
  //copy contents of det_obj into new arrays, get rid of bad energies
  for ( Int_t i=0; i<4; i++ ){

    if ( up[i][1]>0 || down[i][1]>0 ){ 
     
      data_obj.Channel_Up =(Int_t)up[i][0];
      data_obj.Channel_Down =(Int_t)down[i][0];
      
      //if(data_obj.Channel_Up == data_obj.Channel_Down){
      //data_obj.Energy = up[i][1] + down[i][1];

      // assign channel from the active side of the detector if only one is present 
	if ( up[i][1]>0){ 
	  data_obj.Channel = (Int_t)up[i][0];//
	}else if(down[i][1]>0 ){ 
	  data_obj.Channel = (Int_t)down[i][0];//added 20160927
	}
	//}
      
      data_obj.Energy_Up = up[i][1];
      data_obj.Energy_Down = down[i][1];

      data_obj.Energy=0.;//added 20160927
      if (up[i][1]>0){
	data_obj.Energy += up[i][1];
      }
      if (down[i][1]>0){
	data_obj.Energy += down[i][1];
      }
      
      // assign time from the active side of the detector if only one is present 
      if ( up[i][1]>0){ 
	data_obj.Time = up[i][2];
      }else if( down[i][1]>0){ 
	data_obj.Time = down[i][2];
      }
      front.push_back(data_obj);
    }
  }
  
  for ( Int_t i=0; i<Si->det_obj.BackMult; i++ ){
    data_obj.Channel = Si->det_obj.BackChNum[i];
    data_obj.Energy = Si->det_obj.EBack_Cal[i];
    data_obj.Time = Si->det_obj.TBack[i];
    
    //if the energy is less than or equal to zero, discard it
    if (data_obj.Energy>0){
      back.push_back(data_obj);
    }
  } 
  
  if (front.size()<=0 || back.size()<=0){
    return;
  }  
  //-----------------------------------------------------------------
  //sort energies, largest to smallest
  sort(front.begin(), front.end(), Esort_method );
  sort(back.begin(), back.end(), Esort_method );
  //-----------------------------------------------------------------
  ////////////////////////////////////////////////////////////////////////
  Float_t fEn[4] = {0};
  Float_t fEn_Up[4] = {0};
  Float_t fEn_Down[4] = {0};
  Float_t bEn[4] = {0};
  Int_t fCh[4] = {0};
  Int_t bCh[4] = {0};
  Float_t bTi[4]= {0};
  int p=0,q=0,m=0;

  while ((p< front.size()) && (q< back.size()) && (m<4)){
    if (fabs(back[q].Energy - front[p].Energy) < SX3_Qdiff*back[q].Energy + 0.2 ){
      
      fEn[m] = front[p].Energy;
      fEn_Up[m] = front[p].Energy_Up;
      fEn_Down[m] = front[p].Energy_Down;
      bEn[m] = back[q].Energy;
      fCh[m] = front[p].Channel;
      bCh[m] = back[q].Channel;
      bTi[m] = back[q].Time;		
      p++;
      q++;
      m++;      
    }else{      
      // for two front / one back firing
      if (front[p].Energy < back[q].Energy){
	// is there another front-hit that matches as a sum 
	if (((p+1)< front.size()) && (fabs(front[p].Energy + front[p+1].Energy - back[q].Energy)< SX3_Qdiff*back[q].Energy+0.2 )){
	  
	  // first hit 	
	  fEn[m] = front[p].Energy;
	  fEn_Up[m] = front[p].Energy_Up;
	  fEn_Down[m] = front[p].Energy_Down;
	  bEn[m] = front[p].Energy;
	  fCh[m] = front[p].Channel;
	  bCh[m] = back[q].Channel;
	  bTi[m] = back[q].Time;	
	  m++;
	  p++;
	 
	  if(m==4)break;

	  // second hit 
	  fEn[m] = front[p].Energy;
	  fEn_Up[m] = front[p].Energy_Up;
	  fEn_Down[m] = front[p].Energy_Down;
	  bEn[m] = front[p].Energy;
	  fCh[m] = front[p].Channel;
	  bCh[m] = back[q].Channel;
	  bTi[m] = back[q].Time;	
	  p++;
	  q++;
	  m++;	 
 
	}else{
	  break;
	}
      } //   for one front / two back  firing
      else if (front[p].Energy > back[q].Energy){
	// is there another back-hit that matches as a sum 
	if (((q+1)< back.size()) && (fabs(back[q].Energy + back[q+1].Energy - front[p].Energy)< SX3_Qdiff*front[p].Energy+0.2 )){
	  
	  // first hit //in the back
	  fEn[m] = back[q].Energy;
	  fEn_Up[m] = back[q].Energy*front[p].Energy_Up/front[p].Energy;
	  fEn_Down[m] = back[q].Energy*front[p].Energy_Down/front[p].Energy;
	  bEn[m] = back[q].Energy;
	  fCh[m] = front[p].Channel;
	  bCh[m] = back[q].Channel;
	  bTi[m] = back[q].Time;	
	  m++;
	  q++;

	  if(m==4)break;

	  // second hit //in the back	
	  fEn[m] = back[q].Energy;
	  fEn_Up[m] = back[q].Energy*front[p].Energy_Up/front[p].Energy;
	  fEn_Down[m] = back[q].Energy*front[p].Energy_Down/front[p].Energy;
	  bEn[m] = back[q].Energy;
	  fCh[m] = front[p].Channel;
	  bCh[m] = back[q].Channel;
	  bTi[m] = back[q].Time;	
	  p++;
	  q++;
	  m++;

	}else{
	  break;
	}
      }     
      else{
	// this might have been noise, keep what we already had 
	break;
      }
    }
  }//end of while 
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  for(int s=0;s<m;s++){   
    //fill tree

    Si->hit_obj.NHitsInDet = m; 
    Si->hit_obj.DetID = Si->det_obj.DetID;
    Si->hit_obj.HitType = Si->det_obj.BackMult*100 + Si->det_obj.UpMult*10 + Si->det_obj.DownMult;
    //cout<<" DetNum = "<<Si->hit_obj.DetID<<"  HitinDet = "<<m<<"  HitType  = "<<Si->hit_obj.HitType<<endl;

    Si->hit_obj.FrontChannel = fCh[s];
    Si->hit_obj.BackChannel = bCh[s];
    Si->hit_obj.EnergyFront = fEn[s];
    Si->hit_obj.EnergyBack = bEn[s];
    Si->hit_obj.Energy = bEn[s]; 
    Si->hit_obj.Time = bTi[s];
   
    Double_t ZUp=0,ZUpCal=-10, ZDown=0, ZDownCal=-10;
    Double_t xw=0, yw=0, zw=0, rw=0, phiw=0;      

    if(fEn_Down[s]>0 && (fEn_Down[s]>=fEn_Up[s])) {     
      ZDown = ((2*fEn_Down[s]/bEn[s])-1);
      Si->hit_obj.ZDown_Dummy = ZDown;
      //CMAP->PosCal(Si->hit_obj.DetID,Si->hit_obj.FrontChannel,Si->hit_obj.BackChannel,ZDown,ZDownCal);   
      CMAP->PosCal(Si->hit_obj.DetID,fCh[s],bCh[s],ZDown,ZDownCal);  

    }else if(fEn_Up[s]>0 && (fEn_Down[s]<fEn_Up[s])) {
      ZUp = (1-(2*fEn_Up[s]/bEn[s]));    
      Si->hit_obj.ZUp_Dummy = ZUp;  
      //CMAP->PosCal(Si->hit_obj.DetID,Si->hit_obj.FrontChannel,Si->hit_obj.BackChannel,ZUp,ZUpCal);
      CMAP->PosCal(Si->hit_obj.DetID,fCh[s],bCh[s],ZUp,ZUpCal);
    }
      //--------------------------------------------------------------
    //determine z positions using back channels but it is too wide i.e. bad resolution
    //Si->hit_obj.Z = (3-bCh[s] + ZRandom->Rndm())*1.875; 

    //Si->hit_obj.X = XRandom->Rndm()+(3-Si->hit_obj.FrontChannel);
    Si->hit_obj.X = Random->Rndm()+(3-fCh[s]); 
    //-------------------------------------------------------------

    if ((ZDownCal > -1) && (fEn_Down[s] > fEn_Up[s])) {  
      Si->hit_obj.Z = ZDownCal;
      CMAP->GetSX3WorldCoordinates(Si->hit_obj.DetID,Si->hit_obj.X,Si->hit_obj.Z,xw,yw,zw,rw,phiw);
      //if(Si->hit_obj.HitType != 111){ cout<<"Up:  hit_obj.Z ="<<Si->hit_obj.Z<<"hit_obj.X ="<<Si->hit_obj.X<<"  HitType  = "<<Si->hit_obj.HitType<<endl;}
    }else if((ZUpCal > -1) && (fEn_Down[s] <= fEn_Up[s])) {       
      Si->hit_obj.Z = ZUpCal;
      CMAP->GetSX3WorldCoordinates(Si->hit_obj.DetID,Si->hit_obj.X,Si->hit_obj.Z,xw,yw,zw,rw,phiw);
      //if(Si->hit_obj.HitType != 111){ cout<<"Up:  hit_obj.Z ="<<Si->hit_obj.Z<<"hit_obj.X ="<<Si->hit_obj.X<<"  HitType  = "<<Si->hit_obj.HitType<<endl;}
    }  else if(ZUpCal == -10 && ZDownCal == -10) {  
      continue;  
    } else{    
      //continue;
      //cout<<"=============== ========= ==========="<<endl;
    }      

    Si->hit_obj.XW = (xw);
    Si->hit_obj.YW = (yw);
    Si->hit_obj.ZW = (zw);
    Si->hit_obj.RW = (rw);
    Si->hit_obj.PhiW = (phiw);

    Si->Hit.push_back(Si->hit_obj);
    Si->NSiHits++;
    //}
  }
};
////////////////////////////////////////////////////////////////////////////////////////////
