# QQQ Relative Calibration
Relative calibration of Si gains for QQQ. Essentially the same progam as that for the SX3 (or vice-versa).

## Instructions
### General Usage
1. First, run `Main.cpp` with a given `.dat` file all slopes equal to 1
   * Edit `Main.cpp` to initialize the channel map with a `.dat` file with all 1's
   * Compile and run `Main.cpp` to generate a new `.root` file
   * Input the new `.root` and the `.dat` file into this code. Make sure you are using the same `.dat` file that was reference by `Main.cpp`.
   * Run `root -l SiRelativeGains_Step1.C+`
   * This code will output a new `.dat` file with updated gains
2. Second, run Main.cpp again with the new .dat file
   * Edit `Main.cpp` to initialize the channel map with the new `.dat` file generated in Step 1
   * Compile and run `Main.cpp` to generate a new `.root` file
   * Input the new `.root` and the `.dat` file from Step 1 into this code. Make sure you are using the same `.dat` file that was reference by `Main.cpp`.
   * Run  `root -l SiRelativeGains_Step2.C+`
   
### Input files (`.root` files)
The data for the QQQ relaive gains calibration should have high statistics over a range of energies; such as no-target in-gas runs.
This is equivalent to gas calibration runs with the target fully removed, e.g., Spacer 7.

### Output files (`.dat` files)
Output file (e.g.`QQQRelativeGains_Step1.dat`) has the following columns:
Detector number, Front channel, Slope
The first line of `.dat` files is a dummy line.
The QQQ detectors correspond to detectors number 0-3. Each detector has 32 channels (0-31).
Each calibration file should have 129 lines, including the header line.

The `.dat` files are included in the repository as an example. The run-to-run changes in the `.dat` files are excuded by `.gitignore`. To force the updated files to be saved to the repository, use the command `git add -f file.dat`.

For reference, a trivial calibration file can be generated with the file `QQQRelativeGain_ini.C` in the [saves](saves) directory. This program is run using the command `root -l QQQRelativeGain_ini.C`

### Procedure
The QQQ detectors are double-sidded silicon strip detectors divided into rings and wedges. The detector is gain-matched section-by-section. First, the front sections are gain-matched to a particular back section. Then, the remaining back sections are gain-matched to the back section used in the first step. This fitting procedure is refered to "front first" throughout the code.
1. [Step 1](#step-1) `SiRelativeGains_Step1.C`
2. [Step 2](#step-2) `SiRelativeGains_Step2.C`
3. [Final Fix](#final-fix)

In each step, a variety of [Fitting Methods](#fitting-methods) may be specified. The fitting methods are defined in `.h` and are shared by each `SiRelativeGains_Step1.C` and `SiRelativeGains_Step2.C`. The fitting method is selected by modifying the appropriate line in each `.C` file with the following line given as an example.
````
Double_t gain = gainmatch.Fit7(DetNum,FrontChNum,BackChNum);
````
In this example Fit7 is specified.

## Step 1
Loop over front channels.
For each detector (`DetNum=0-3`), each front channel (`FrontChNum=0-15`) is gain-matched with a particular back channel (typically `BackChNum=0`).
### Histograms
`Q3_back_vs_front%i_%i_%i`
### Files
Channels 16-31 in the `.dat` file will be filled in; this is written in the code as '[FrontChNum+16]'.
### Gains
The measured slope is multiplied by the previous gain (typically 1) to produce the new gain. Using this method, that is, using fits on histotrams to generate the slope parameters, the results are good to within 0.1 or better. That is the final slopes will be within 1.000 +/- 0.003 or better.
### Next steps
Once you have completed this program, edit, recompile, and re-run `Main.cpp` using the new `QQQRelativeGains_Step1.dat` file.
Once complet, verify that the slope of each histogram measured in Step 1 is now equal to 1.00.
You will input this new root file with this new relative gains file into [Step 2](step-2).

## Step 2
Loop over back channels. For each detector (`DetNum=0-3`), each new back channel (`BackChNum=1-15`) is gain-matched with a particular front channel (typically `FrontChNum=0`). Since `FrontChNum=0` was already gain matched to `BackChNum=0` in Step 1, the loop runs over the range 1-15.
### Histograms
`Q3_back_vs_front%i_b%i`
### Files
Channels 1-15 in the `.dat` file will be filled in. 
### Gains
Using this method, that is, using fits on histotrams to generate the slope parameters, the results are good to within 0.3%. That is the final slopes will be within 1.000 +/- 0.003 or better.
### Next steps
Check historams `Q3_back_vs_front%i`

## Final Fix
The "final fix" step in the calibration corrects the offset between the front and back signals. This performed automatically in Step 1 and Step 2.
### Histograms
`q3offset_back_vs_front%i`
### Files

### Offsets
The offset for each detector segment should be on the order of +/-30 channels or smaller.
### Next steps
Check that the back vs. front offset for each detector is consistient with zero.

## Fitting Methods 
The file `SiRelativeGains.h` provides the following fit methods.

1. Method 1 - calculates slope of points wihtin pre-defined cut using TGraph
2. Method 2 - calculates slope of points wihtin user-defined cut using TGraph
   This method works very similar to method 1, except that instead of a predefined cut, the user must
   draw their own graphical cut. To do this, in the canvas: View->Toolbar and click on the scissors on the top
   right hand corner. Then, select the region around your data by clicking. Double click to close the cut.
   A best fit line should appear through your data
3. Expanded Method 1 - three slopes are calculated:		
   1) slope from un-gated fit of 2D histogram		
   2) slope from gated fit of TGraph (same as Method 1)		
   3) slope from un-gated fit of profile of 2D histogram (uses the average y-position for each x-position)
4. Automated cut generation based on TProfile slope. A parrallelpiped cut is generated.
5. Expansion of Method 4. Used for Det 2. Uses upper and lower limits for the lower bound of the gate; i.e., the gate "backs in" towards the origin over the steps.
6. Automated cut generation based on TProfile slope; cone-shaped cut. Calculates slopes within 0.2% of those of Mehtod 4 with the advantage of beign able to better detect the slope for detectors with multiple loci.

## Images
In the `images` folder there are two macros. The image in the folder was generated using the macro `Draw2.C`, however, that macro is not compatible with `Main.cpp`.

```
root -l Draw.C

```
