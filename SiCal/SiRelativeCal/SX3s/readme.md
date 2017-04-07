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
The SX3 detectors correspond to detectors number 4-27. Each detector has 12 front channels (0-11).
The `.dat` files are included in the repository as an example. The run-to-run changes in the `.dat` files are excuded by `.gitignore`. To force the updated files to be saved to the repository, use the command `git add -f file.dat`.

## Step 1
This program fixes the relative gains for the up and down on the SX3
How to use: Create Histograms in the `Main.cpp`:
* Plot down vs up energy for each front channel of each detector
* Code reads in histogram of name `down_vs_up%i_front%i`--e.g. `down_vs_up4_front0` (det4, channel0) histo "down_vs_up_divideBack%i_front_%i" is a new extra histo that was created in the OrganizeIC.cpp for data that are normalized with the BackEnergy
  if your data is not normalized use histo "down_vs_up%i_front_%i" for this code. 
* It can loop over any range of detectors you want You can do the calibration detector by detector looping one detector at a time or loop from DetNum=4 to 27.
  if you do so change the loop below.

The input `.root` file should be generated from runs with fixed particle energy; either alpha calibration or proton scattering. The calibration assumes that up+down=constant. This is only true for a fixed energy.

The program reads in an `X3RelativeGains.dat` file and outputs a new file with updated coefficients
Before running, make sure that the `.root` file you are reading in has the right histograms and has been created using the `X3RelativeGains.dat` file that you are inputting into this code

if you are just starting the calibration you can use an `.dat` input file where ALL SLOPES ARE ONE(1) apart from the MASKED CHANNELS which are ZERO(0)

  input the root file that was created in the Main(Organize) using the `.dat` file where ALL SLOPES are ONE.

This file changes the relative gains on the down relative to the up
Once you have completed this program, rerun `Main.cpp` with this new `X3RelativeGains_Step1.dat`
You will input this new root file with this new relative gains file into step 2.

### Options
 

### Methods 
1. Method 1 - calculates slope of points wihtin pre-defined cut using TGraph. developed from Step 1 in 'Old' directory. if it is necessary to limit the area of your data that you want to fit see in our Canvas and 
  input below on the CUT the coordinates of the points that surround this area.
2. Method 2 - calculates slope of points wihtin user-defined cut using TGraph manual cut, from FrontFirst directory
   This method works very similar to method 1, except that instead of a predefined cut, the user must
   draw their own graphical cut. To do this, in the canvas: View->Toolbar and click on the scissors on the top
   right hand corner. Then, select the region around your data by clicking. Double click to close the cut.
   A best fit line should appear through your data
3. Automatic calculation of cut.

## Step 2


