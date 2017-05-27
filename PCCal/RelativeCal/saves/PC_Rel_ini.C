// create calibration file with trivial parameters
#include <fstream>
void PC_Rel_ini (void) {
  Double_t gain = 1;
  ofstream outfile; 
  outfile.open("PCWire_RelGain_init.dat");
  outfile << "Wire\tGain\n";
  for (Int_t WireNum=0; WireNum<24; WireNum++){
    outfile << WireNum << "\t"<< gain <<endl;
  }
}
