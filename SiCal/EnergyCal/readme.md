# Silicon detector energy calibration
The energy calibration macro calibrates the back detector signals only.
Therefore, the energy calibration is not dependent on the position calibration.
## Instructions
### General Usage
Usage: `root -l Alpha_All.C+` (from the same directory).
### Output files (.dat files)
Each of the 27 detectors has one calibration.
The output file should have either 28 lines, including the header; or one more line than the detectors with good energy signals.
### Input files
A root file with multiple known eneryg peaks should be used.
## Details
### Histograms
The energy is read out from the data tree.
No pre-existing histogram is used.
### Fiting Method
The energy calibration is determined using a peakfit routine.
### Gains
Since the pedistal, offset, etc. have been carefully determined in the previous calibration steps, only one parameter is used for the energy calibration: the gain.
The ChannelMap methods Get MeVPerChannel are used.
These methods are called in Main.
## Next steps



