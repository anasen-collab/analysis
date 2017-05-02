///////////////////////////////////////////////////////////////////////////////////////////
////Do alpha calibration
#include <TGraph.h>
#include <TF1.h>

void SimpleFit(void) {
  const Int_t size = 6;
  Double_t energy[size] = { 5.42315, 5.68537, 6.05, 6.28808, 6.77803, 8.74886 };
  Double_t channel[size] = { 2474, 2605, 2781, 2890, 3122, 4074 };//for detector 16
  
  TGraph *graph = new TGraph(size, channel, energy);
  graph->Draw("A*");

  TF1 *myfit = new TF1("myfit","[0]+[1]*x",0,6000);
  graph->Fit("myfit");
}
