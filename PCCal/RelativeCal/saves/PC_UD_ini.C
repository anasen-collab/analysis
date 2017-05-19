// create calibration file with trivial parameters
#include <fstream>
void PC_UD_ini (void) {
  Double_t slope = 0;
  Double_t shift = 1;
  ofstream outfile; 
  outfile.open("PC_UD_RelCal_init.dat");
  outfile << "Wire\tSlope\tShift\n";
  for (Int_t WireNum=0; WireNum<24; WireNum++){
    outfile << WireNum << "\t"<< shift << "\t" << slope <<endl;
  }
}
