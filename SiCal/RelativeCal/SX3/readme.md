# Super-X3 Relative Calibration
Relative calibration of silicon detector gains for SX3 detectors.
The gains are calibrated in four steps (three steps plus a "final fix").
First, the two position signals are gain-matched for each strip on the front of the detectors.
Then, the front signals are gain-matched to a particular back signal, thus gain-matching the front the the detector.
Next, the back signals are gain-matched to each other.
Last, the a linear offset is appplied to align the sum of the position signals with the energy signal.

## Instructions
### General Usage
1. First, run `Main.cpp` with a given `X3cal.dat` file
   * Input the `.root` file and the `.dat` file into `Step1.C`
   * Run `root -l SiRelativeGains_Step1.C+`
   * This code will output another `.dat` file.
2. Second, run `Main.cpp` again with the new `.dat` file from Step 1.
   * Input the new `.root` and new `.dat` file into `Step2.C`
   * `root -l SiRelativeGains_Step2.C+`
3. Finally, run `Main.cpp` again with the new `.dat` file from Step 2.
   * Input the new organize.root and new .dat file into `Step3.C`
   * `root -l SiRelativeGains_Step3.C+`

### Input files (`.root` files)
The input `.root` file should be generated from runs with fixed particle energy; either alpha calibration or proton scattering.
The calibration assumes that up+down=constant.
This is only true for a fixed energy.

### Output files (`.dat` files)
The output file (e.g.`X3RelativeGains_Slope1.dat`) has the following columns:
Detector number, Front channel number, Slope.
The first line of dat files is a header line.
The SX3 detectors correspond to detectors number 4-27. Each detector has 12 channels (0-11). The channels are enumerated as follows
* 0- 3 back  channel                 (total energy)
* 4- 7 front channel, up? position   (partial energy)
* 8-11 front channel, down? position (partial energy)
Each calibration file should have 289 lines, including the header line.

Note, that if you loop over a subset of detectors, only the paramaters for those detectors will be written into the new `.dat` file.

The `.dat` files  included in the repository are an example. The run-to-run changes in the `.dat` files are excuded by `.gitignore`. To force the updated files to be saved to the repository, use the command `git add -f file.dat`.

For reference, trivial calibration files can be generated with the files `X3RelativeGain_ini.C` and `X3FinalFix_ini.C` in the [saves](saves) directory. These programs are  run, for example, using the command `root -l X3RelativeGain_ini.C`

### A note on history
In the previous development of the code, contained in the `Old` folder, the backs were gain-matched first. That is, in Step 2, each back channel was gain-matched to a particular front channel. This is true for both the standard and "clickable" sets of code.

while looping over the front, both new and old sets of code use the histogram 
`hist = (TH2F*)f1->Get(Form("back_vs_front%i_%i_%i",DetNum,FrontChNum,BackChNum));` 
This same histogram is used for looping over the back in Step 2 old. In the new step which loops over the back, Step 3, the following histogram is used.
`hist = (TH2F*)f1->Get(Form("back_vs_front%i_back%i",DetNum,BackChNum))` (New)

The fit function used in `Old` was a one-parameter scaling coefficient. In all other versions, the fit function is a linear function.

## Step 1
This program fixes the relative gains for the up and down on the SX3.
### Histograms
Plot down vs up energy for each front channel of each detector.
These histograms are created in in `Main.cpp` with the following code.
````
//Step 0 RelCal/U-D
	  MyFill(Form("down_vs_up%i_f%i",Si.det_obj.DetID,Si.det_obj.UpChNum[0]),512,0,16384,Si.det_obj.EUp_Pulser[0],512,0,16384,Si.det_obj.EDown_Pulser[0]);
	  MyFill(Form("down_vs_up_divideBack%i_front%i",Si.det_obj.DetID,Si.det_obj.UpChNum[0]),100,0,1,(Si.det_obj.EUp_Cal[0]/Si.det_obj.EBack_Cal[0]),100,0,1,(Si.det_obj.EDown_Cal[0]/Si.det_obj.EBack_Cal[0]));
````

Histograms named `down_vs_up_divideBack%i_front_%i` are an extra histogram, also created in `Main.cpp` for data that are normalized with the BackEnergy. If your data is not normalized use the histograms previously mentioned.
`hist = (TH2F*)f1->Get(Form("down_vs_up_divideBack%i_front%i",DetNum,FrontChNum));`

* Code reads in histogram of name `down_vs_up%i_front%i`--e.g. `down_vs_up4_front0` (`DetNum=4`, `FrontChNum=0`) 

* It can loop over any range of detectors you want You can do the calibration detector by detector looping one detector at a time or loop from DetNum=4 to 27.
  if you do so change the loop
### Files
The program reads in an `X3RelativeGains.dat` file and outputs a new file with updated coefficients.
Before running, make sure that the `.root` file you are reading in has the right histograms and has been 
created using the same `X3RelativeGains.dat` file that you are inputting into this code.
If you are just starting the calibration, you can use a `.dat` input file where all sopes are one (1) apart from the masked channels which are zero (0).
A trivial file can be generated using the program `X3RelativeGain_ini.C`

### Gains
This file changes the relative gains on the down relative to the up.
The gains are applied as -[old]/[new] where old is the previous gain read in by the file and new is the measured slope of the histogram. Channel numbers 8-11 will be updated.

### Next steps
Once you have completed this program, rerun `Main.cpp` with this new `X3RelativeGains_Step1.dat`
You will input this new root file with this new relative gains file into Step 2.

## Step 2
Loop over front (same as "clickable" Step 3)
This step fixes the relative gains of the 4 front strips with respect to a single back strip
### Histograms
First, histograms should be created for events in which one front strip (up and down) and one back strip fired. Using more complicated multiplicities here will confuse things (if only the up fired, then front != back, which defeats the initial assumption that front == back).

The root file should have histograms that are plotting `back_vs_front`

 * Code reads in histogram of name `back_vs_front%i_0_%i`--e.g. `down_vs_up4_0_1` (`det4`, `front channel0`, `back channel 1`)
qy`hist = (TH2F*)f1->Get(Form("back_vs_front%i_%i_%i",DetNum,FrontChNum,BackChNum));`
 * If the front channel 0 does not exist for a given detector, choose a different front channel for everything to be relative to.
 * It can loop over any range of detectors you want, 
```` 
//Step 2 RelCal//F-B //Condition: RelGain Cal from Up-Down is applied
	`MyFill(Form("back_vs_front%i_%i_%i",Si.det_obj.DetID,Si.det_obj.UpChNum[0],Si.det_obj.BackChNum[0]),512,0,16384,Si.det_obj.EUp_Rel[0]+Si.det_obj.EDown_Rel[0],512,0,16384,Si.det_obj.EBack_Rel[0]);`
````
### Files
The program reads in an `X3RelativeGains.dat` file and outputs a new `.dat` file with updated coefficients.
Before running, make sure that the `.root` file you are reading in has the right histograms and has been created using the `X3RelativeGains.dat` file that you are inputting into this code. That is, whichever `.dat` file was used in the program that creates the histograms (`Main.cpp`) should be the input file to this program.

### Gains
The gains are applied as [old]*[new] where old is the previous gain read in by the file and new is the measured slope of the histogram.
The two corresponding front gains (up and down) are then multiplied by the slope to get the new gain.
this is why the output indices are `[DetNum-4][FrontChNum+4]` and `[DetNum-4][FrontChNum+8]`.
Channel numbers 4-11 will be updated.
 
#### Next steps
Once you have completed this program, rerun `Main.cpp` with this new `X3RelativeGains_Step2.dat`
You will input this new root file with this new relative gains file into Step 3. 

## Step 3
Loop over back; (same as "old" and "clickable" Step 2).
This program fixes the back gains in the SX3 detectors relative to a front channel.
This step fixes the relative gains of the 3 remaining back channels to that of a designated front channel.
It works in the same was as Step 2.
### Histograms
The root file should have histograms that are plotting `back_vs_front`
Code reads in histogram of name `back_vs_front%i_0_%i`--e.g. `down_vs_up4_0_1` (`det4`, `front channel0`, `back channel 1`)
`hist = (TH2F*)f1->Get(Form("back_vs_front%i_%i_%i",DetNum,FrontChNum,BackChannelNumber));` (Old)
`hist = (TH2F*)f1->Get(Form("back_vs_front%i_back%i",DetNum,BackChNum))` (New)
Step 3 RelCal//F-B //RelGain Cal from Step 1 is applied
`MyFill(Form("back_vs_front%i_b%i",Si.det_obj.DetID,Si.det_obj.BackChNum[0]),512,0,16384,Si.det_obj.EUp_Rel[0]+Si.det_obj.EDown_Rel[0],512,0,16384,Si.det_obj.EBack_Rel[0]);`

First, histograms should be created for events in which one front strip (up and down) and one back strip fired. Using more complicated multiplicities here will confuse things (if only the up fired, then front != back, which defeats the initial assumption that front == back).

 * If step 2 fixed the gains relative to back (front) channel 0, in this code, you should loop over bakc (front) channels 1,2,3
 * If the back channel 0 does not exist for a given detector, choose a different front channel for everything to be relative too (in this code, det 11 uses back channel 1).
 * If the front channel 0 does not exist for a given detector, choose a different front channel for everything to be relative to.

### Files
The program reads in an `X3RelativeGains_Step2.dat` file and outputs a new `.dat` file with updated coefficients.
Before running, make sure that the `.root` file you are reading in has the right histograms and has been created using the same `X3RelativeGains.dat` file that you are inputting into this code (presumably, the `.dat` file generated by Step 2).
Whichever gains file was used in the program that creates the histograms should be the input file to this program.

### Gains
The back gain is equal to the inverse of the slope of the best fit line.
The gains are applied as [old]/[new] where old is the previous gain read in by the file and new is the measured slope of the histogram.

### Next steps
Once you have completed this program, rerun `Main.C` with this new `X3RelativeGains_Step3.dat`
Everything should now be calibrated properly, but check the histograms to make sure.
If the histograms are not as good as desired, you can repeat this process for any individual detector (all three programs). Just be sure to input the correct `RelativeGains.dat` file and root file.

## Final fix (Step 4)
### Hisotgrams
````
//check offset
  	  MyFill(Form("sx3offset_back_vs_front%i",Si.det_obj.DetID),512,0,16384,Si.det_obj.EUp_Rel[0]+Si.det_obj.EDown_Rel[0],200,-400,400,(Si.det_obj.EBack_Rel[0]-(Si.det_obj.EUp_Rel[0]+Si.det_obj.EDown_Rel[0])));
	  
	  MyFill(Form("sx3offset_back_vs_front_normback%i",Si.det_obj.DetID),512,0,16384,Si.det_obj.EBack_Rel[0],300,-0.2,0.2,((Si.det_obj.EBack_Rel[0]-(Si.det_obj.EUp_Rel[0]+Si.det_obj.EDown_Rel[0]))/Si.det_obj.EBack_Rel[0]));

	  MyFill(Form("offset_check%i",Si.det_obj.DetID),512,0,16384,Si.det_obj.EBack_Rel[0],100,-1,1,((Si.det_obj.EDown_Rel[0]-Si.det_obj.EUp_Rel[0])/Si.det_obj.EBack_Rel[0]));
	  MyFill(Form("offset2_check%i",Si.det_obj.DetID),512,0,16384,Si.det_obj.EBack_Rel[0],100,-1,1,((Si.det_obj.EDown_Rel[0]-Si.det_obj.EUp_Rel[0])/(Si.det_obj.EUp_Rel[0]+Si.det_obj.EDown_Rel[0])));
````
### Files
### Gains
### Next steps
Once you have completed this program, rerun `Main.cpp` with this new `X3FinalFix.dat`


## Fiting Methods 
The fitting methods have been collected in the file `SiRelativeGains.h` and are shared by each step.

1. `Fit1`
   Convert the histogram bin-by-bin to a TGraph and fit.
   It works by dumping the x-y coordinates of a 2D histogram into a `TGraph`.
   The coordinates are weighted by the bin content of the 2D hist
   It then fits the `TGraph` with a function of the form m*x+b.
   
   The do-cut flag is added to method 1 to turn on or off the use of a gate.
    * `docut=kFALSE` Method 0 - Ignore the graphical cut and convert the entire histogram bin-by-bin to a `TGraph` and fit.
   Used by 'SX3s/Step2,3' and 'Old/Step2,3'
    This method shows that a `TGraph` should be generated directly from the data tree.
    * `docut=kTRUE` Method 1 - Calculate the slope of points wihtin pre-defined cut using `TGraph`. 
	If it is necessary to limit the area of your data that you want to fit see in our Canvas and input below on the CUT the coordinates of the points that surround this area.
  Only the coordinates that fall within a predefined cut are accepted. This gets rid of noise that can affect the fit.
  However, on first iteration, not all of the data will fit neatly into the cut, so it may need to be adjusted for a few detectors. 
  Developed from Step 1 in 'Old' directory and  Step 3 in QQQ 'Method 2'
  Note: Fit4 can be used to calulate the gate automaticly.
2. `Fit2`
   Method 2 - calculates slope of points wihtin user-defined cut using TGraph manual cut.
   This method works very similar to Method 1, except that instead of a predefined cut, the user must
   draw their own graphical cut. To do this, in the canvas: View->Toolbar and click on the scissors on the top
   right hand corner. Then, select the region around your data by clicking. Double click to close the cut.
   A best fit line should appear through your data.
   Developed from 'FrontFirst' directory.
3. `Fit3`
   Use a `TCutG` to draw a line over the data. in canvas: View->Toolbar->GraphicalCut (pair of scissors on right).
   The verticies of the cut will be used to generate a linear fit. The last point, the point which is double-cliked to close the cut, is not included.
   All of the segments for a detector should have the same slope; the plot should look like a straight line. 
   Click along the straight line. When you are done, double click in canvas and a best fit line will appear. The best fit line should follow the data very well, if not you are doing something wrong.
   If the back gains are not set properly (or when doing front-first calibration), the line may appear segmented. If so, choose your favorite segment and get a best fit line for that. Do not click in each segment as that will throw your best fit line off. 
After the best fit line appears, double click to move onto the next channel.
   Used by 'SX3s/Clickable' and 'Old/Clickable_Step2,3' 
3. `Fit4`
   Automatic calculation of cut. The user provides an estimate of the slope. If the given slope is less than 0, it is assumed the data is of the form `up_vs_down` and the offset is estimated by the highest bin with a non-zero count. The shape of the grahical cut is a parrallelpiped centered on a line of a given slope and offset. Once the fit is determined, the width of the graphical cut is reduced and position of the cut is refined using the slope and offset from the previous fit.
   
   Since the graphical cut is calculated automatically, no user input is required. Use of the commands `can->Update();` and `can->WaitPrimitive();' can be used to monitor the performance of the program, but dramatically increase the amount of time the program takes to run.
