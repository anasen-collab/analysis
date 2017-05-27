// create calibration file with trivial parameters
#include <fstream>
void PC_UD_ini (void) {
  Double_t slope = 1;
  Double_t offset = 0;
  ofstream outfile; 
  outfile.open("PC_UD_RelCal_init.dat");
  outfile << "Wire\tSlope\tOffset\n";
  for (Int_t WireNum=0; WireNum<24; WireNum++){
    outfile << WireNum << "\t"<< slope << "\t" << offset <<endl;
  }
}
