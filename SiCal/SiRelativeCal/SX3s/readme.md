# Super-X3 Relative Calibration
Relative calibration of Si gains

## General Usage
### Instructions
In all three steps

If you loop over a subset of detectors, only the paramaters for those detectors will be written into the new `.dat` file
1. First, run `Main.cpp` with a given `X3cal.dat` file
  * Input `organize.root` and the `.dat` file into this code
  * This code will output another .dat file that you should define
  * `root -l SiRelativeGains_Step1.C++`
2. Second, run Main.cpp again with the new .dat file
  * Input the new organize.root and new .dat file into Step2.C
  * `root -l SiRelativeGains_Step2.C++`
3. Finally, run Main.cpp again with the new .dat file
  * Input the new organize.root and new .dat file into Step3.C
  * `root -l SiRelativeGains_Step3.C++`

### .dat files
Output file (e.g.`X3RelativeGains_Slope1.dat`) has the following columns:
Detector number, Front channel, Slope
The first line of dat files is a dummy line.
The SX3 detectors correspond to detectors number 4--27. Each detector has 12 front channels (0--11).
The .dat files are included as an example. The run-to-run .dat file changes are excuded by .gitignore. The force a save of updated files use the command `git add -f file.dat`.

## Step 1
This program fixes the relative gains for the up and down on the SX3
How to use: Create Histograms in the `Main.cpp`:
* Plot down vs up energy for each front channel of each detector
* Code reads in histogram of name `down_vs_up%i_front%i`--e.g. `down_vs_up4_front0` (det4, channel0)
* It can loop over any range of detectors you want, but if the histogram does not exist, the code will crash and your work will not be saved

The program reads in an `X3RelativeGains.dat` file and outputs a new file with updated coefficients
Before running, make sure that the `.root` file you are reading in has the right histograms and has been created using the `X3RelativeGains.dat` file that you are inputting into this code

This file changes the relative gains on the down relative to the up
Once you have completed this program, rerun `Main.cpp` with this new `X3RelativeGains_Step1.dat`
You will input this new root file with this new relative gains file into step 2.

Edited by : John Parker , 2016Jan22
Edited by : Maria Anastasiou, 2016Sept20

### Options
 if it is necessary to limit the area of your data that you want to fit see in our Canvas and 
  input below on the CUT the coordinates of the points that surround this area.

 if you are just starting the calibration you can use an .dat input file where ALL SLOPES ARE ONE(1) apart from the MASKED CHANNELS which are ZERO(0)

  input the root file that was created in the Main(Organize) using the .dat file where ALL SLOPES are ONE.

You can do the calibration detector by detector looping one detector at a time or loop from DetNum=4 to 27.
  if you do so change the loop below.

  histo "down_vs_up_divideBack%i_front_%i" is a new extra histo that was created in the OrganizeIC.cpp for data that are normalized with the BackEnergy
  if your data is not normalized use histo "down_vs_up%i_front_%i" for this code. 

## Step 2


