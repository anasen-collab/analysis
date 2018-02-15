# PC pulser calibration.
This program is used to calibrate the raw PC data from ADC channel number into voltages.
During calibration, a variety of pulser signals of known voltages is supplied to the PC preamplifers.
By comparing the peak positions in the PC spectrum with the known voltages from the pulser, a calibration is developed.

The PC is calibrated by applying negative polarity tail pulses to the Mesytec preamplifers located inside the ANASEN chamber.
Typical votage amplitudes used for calibration are 0.1V down to zero in 5 steps. 
 
## Instructions
 Voltages and file locations are entered into the code.
### Standard
The program is run using the following command.
````
root -l PCPulser_All.C
````
### Zero
To test the inclusion of the pedistal peak at "zero" use the command
````
root -l PCPulser_zero.C
````
It is not recommended to use the files generated which include the zero peak. Doing so will optimize for the zero peak and not for gain matching for position calibration.

## Input files
The program is intened to be run on "raw" `.root` files that have been processed by the `evt2root` progam and assumes that the voltage values from the pulser are known.

## Output files
The output files are saved in the `saves` directory.
The output file (e.g.`PCpulserCal.dat`) has the following columns:
 ID, Chan, zero shift (fit offset), voltage per can (fit slope)

The first line of dat files is a dummy line.
The PC channels correspond to ADCs 2-3.
ADC 2  has 32 channels (0-31) and ADC 3 has 16 channels (0-15).
Each calibration file should have 49 lines, including the header line.

The `.dat` files are included in the repository as an example. The run-to-run `.dat` file changes are excuded by `.gitignore`. The force a save of updated files use the command `git add -f file.dat`.

For reference, trivial calibration files can be generated with the file `PCPulser_ini.C` in the [saves](saves) directory.
This program is run using the command `root -l PCPulser_ini.C`

### PCpulserCal.dat
Fit of peaks from `xpeaks`, center of bin of heighest peak from TSpectrum peakfit
Linear fit with ROB, fit function `fit`
outliers are excluded
### PCpulserCal_full.dat
Fit of peaks from `xpeaks`
Linear fit, without ROB, `fit3`
in other words, the full distribution is fit with no exclusion
### PCpulserCal_centroid.dat
Fit of peaks from `gpeaks`, centroids of gaussian fits
Linear fit with ROB `fit4`
## Next steps
Once the calibration file has been generated, modify the program `Main.cpp` to read in the calibration file.
Once processed through `Main.cpp` the resultant `.root` files will contain PC spectra that have been calibrated to voltages.

To shift the "zero" peak to zero use the Main definition rezero and set rezero value to the mean of the zero peak in the voltage sum distribution.

After voltage calibration, the PC signals should be calibrated for gain matching using the macros located in [RelativeCal](../RelativeCal).
