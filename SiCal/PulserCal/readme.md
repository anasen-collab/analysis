# Silicon Pulser Calibration

Use `SiPulser_All.C` to perform the silicon pulser calibration.
Use `SiPulser_test.C` to test the results of the silicon pulser calibration.

Usage: (from the same directory).
````
root -l SiPulser_All.C 
````
 Output file (e.g.`Sipulser_2015Dec13.dat`) has the following columns:
 MBID, CBID, ASICs_Channel, ZeroShift(offset), Voltage_per_Ch(slope)
