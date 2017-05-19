 /////////////////////////////////////////////////////////////////////////////////////////////
////////////// checking for the heavy hit &/or cross talk in the wire ///////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////    
#ifdef PCPlots
Int_t pct;
for(Int_t pc=0; pc<Tr.NTracks; pc++){
  //cout<<"   Check 10 " <<Tr.NTracks<<endl;

  //All PCWire vs Energy
  MyFill("WireID_vs_PCEnegy",25,0,24,Tr.TrEvent[pc].WireID,500,0,2,Tr.TrEvent[pc].PCEnergy);

  for(Int_t pca=0; pca<Tr.NTracks1; pca++){

    //PCWire with track in channel 12 & other tracks & non-tracks in other channels
    MyFill("WireID_mod1_vs_PCEnegy",25,0,24,(((Int_t)Tr.TrEvent[pc].WireID-(Int_t)Tr.TrEvent[pca].WireID +12)%24),500,0,2,Tr.TrEvent[pc].PCEnergy);
  }

 }
#endif
/////////////////////////////////////////////////////////////////////////////////////////////
