ANASEN Analysis Software
=

## Getting the code
Make sure you have Git installed on your local computer. Use the following command to download a copy of the git repository.
''''
git clone https://github.com/jonlighthall/analysis.git
''''
Use Git to track changes in your local copy of the code and to receive updates.
## Procedure 
1. [Conversion](#Conversion)
2. [Calibration](#Calibration)
3. [Analysis](#Analysis)
## Conversion
Before the data files can be calibrated and analysed, they must fist be converted. The programs used for converting `.evt` data files into `.root` files is located in the `evt2root` directory. The output of these programs are un-calibrated and un-analyzed. 
## Calibration
In the folder `analysis_software`, the program `Main.cpp` is used for converting the "raw" `.root` files produced by `evt2root` into `.root` files used for calibration.
### Silicon detectors
The macros for calibrating the silicon detectors are located in `SiCal`.
### Proportional counter
The macros for calibrating the proportional counter detector are located in `PCCal`.
## Analysis
