// create calibration file with trivial parameters
#include <fstream>
void PCPulser_ini (void) {
  Double_t zeroshift = 0;
  Double_t vperch = 1;
  ofstream outfile; 
  outfile.open("PCpulser_ini.dat");
  outfile << "ID\tChan\tOffset\tVolts/Chan" << endl;
  for (Int_t id=2; id<4; id++) {  
    for (Int_t chan=0; chan<32; chan++) {
      if(id==3&&chan>15)
	continue;
      outfile << id << "\t" << chan << "\t"<< zeroshift << "\t" << vperch <<endl;
    }
  }
}
