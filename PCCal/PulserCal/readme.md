# ANASEN proportional counter (PC) pulser calibration.
This program is used to calibrate the raw PC data from ADC channel number into voltages. During calibration, a variety of pulser signals of known voltages is supplied to the PC preamplifers. By comparing the peak positions in the PC spectrum with the known voltages from the pulser, a calibration is developed.
 
## Instructions
The program is intened to be run on "raw" `.root` files that have been processed by the `evt2root` progam and assumes that the voltage values from the pulser are known. Voltages and file locations are entered into the program and run using the following command.
````
root -l PCPulser_All.C
````

## Output files
The output files are saved in the `saves` directory.
The output file (e.g.`PCpulserCal.dat`) has the following columns:
 ID, Chan, zero shift (fit offset), voltage per can (fit slope)

The first line of dat files is a dummy line.
The PC channels correspond to ADCs 2-3.
ADC 2  has 32 channels (0-31) and ADC 3 has 16 channels (0-15).

The `.dat` files are included in the repository as an example. The run-to-run `.dat` file changes are excuded by `.gitignore`. The force a save of updated files use the command `git add -f file.dat`.

## Next steps
Once the calibration file has been generated, modify the program `Main.cpp` to read in the calibration file. Once processed through `Main.cpp` the resultant `.root` files will contain PC spectra that have been calibrated to voltages.
