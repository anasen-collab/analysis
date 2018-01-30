# Silicon Pulser Calibration

## Instructions
Use `SiPulser_All.C` to perform the silicon pulser calibration.

Usage: (from the same directory).
````
root -l SiPulser_All.C 
````

## Output files
Output file (e.g.`Sipulser_2015Dec13.dat`) has the following columns:
 MBID, CBID, ASICs_Channel, ZeroShift(offset), Voltage_per_Ch(slope)
 
There are 416 total silicon channels to calibrate: 4 QQQ detectors with 32 channels each (128 channels) and 24 SX3 detectors with 12 channels each (288 channels).

## Next steps
Once the calibration file has been generated, modify the program `Main.cpp` to read in the calibration file. Once processed through `Main.cpp` the resultant `.root` files will contain silicon spectra that have been calibrated to voltages.

 `SiPulser_test.C` may be used to test the results of the silicon pulser calibration.

After voltage calibration, the silicon signals should be calibrated for gain matching using the macros located in [RelativeCal](../RelativeCal).
