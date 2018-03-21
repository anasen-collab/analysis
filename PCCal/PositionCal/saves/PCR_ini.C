// create calibration file with trivial parameters
#include <fstream>
void PCR_ini (void) {
  Double_t gain = 3.85;
  ofstream outfile; 
  outfile.open("PCR_init.dat");
  outfile << "WireID\tRadius\n";
  for (Int_t WireNum=0; WireNum<24; WireNum++){
    outfile << WireNum << "\t"<< gain <<endl;
  }
}
