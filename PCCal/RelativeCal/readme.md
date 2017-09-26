script to do the PC UP vs Down calibration. This needs to be done right after the PC pulser calibration. I have only attached the changes in the files.

Here is how you can use it: you first have to create 2D histos in the main for up vs down (you can look in the Main file attached), then use the script to calculate the slope and offset. Then load the file with slope and offset in the main (add these changes files  to your Main.cpp & ChannelMap.h).
