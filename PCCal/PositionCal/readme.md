# PC Position Calibration

This program calculates the gain and offset parameters to scale the calculated PC position to match the projected PC position.

# Usage: 
````
root -l PosCal.C
````

# Requirements

The histograms referenced in this program are defined in Analyzer.cpp
The calibration constants are calculated by comparing the projected PC Z-position, called PCZ_Ref, to the calculated PC Z-position.

histograms can be gated on the silicon enegy. For example, for a 10MeV beam, the silicon energy can be gated on the range 9.4-9.7 Mev
