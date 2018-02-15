This program fixes the relative gains for the up and down on the PC

# Usage: 
````
root -l PosCal.C
````

# Requirements

The histograms referenced in this program are defined in Analyzer.cpp
The calibration constants are calculated by comparing the projected PC Z position, called PCZ_Ref to 

histograms can be gated on the silicon enegy. For example, for a 10MeV beam, the silicon energy can be gated on the range 9.4-9.7 Mev
