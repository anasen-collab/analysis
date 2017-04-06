 ANASEN proportional counter (PC) pulser calibration.
 Developed by : Jon Lighthall Nov 2016
 
## Instructions
 Usage: root -l PCPulser_All.C 

### Output files
The output files are saved in the `saves` directory.
The output file (e.g."PCpulserCal.dat") has the following columns:
 ID, Chan, zero shift (fit offset), voltage per can (fit slope)

The first line of dat files is a dummy line.
The PC channels correspond to ADCs 2-3.
ADC 2  has 32 channels (0-31) and ADC 3 has 16 channels (0-15).

The `.dat` files are included in the repository as an example. The run-to-run `.dat` file changes are excuded by `.gitignore`. The force a save of updated files use the command `git add -f file.dat`.
