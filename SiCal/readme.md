# ANASEN silicon detector array calibration
The silicon detector array is comprised of two structures in three groups: (1) A forward endcap formed by four Micron QQQ detectors and (2) two rings of Micron SX3 detectors, with 12 detectors in each ring. Both detector types follow the same basic calibration procedure, but the SX3 detectors require additional position calibration.
 
## Order of calibration
 1. [Pulser calibration](#pulser-calibration)
 2. [Relative calibration](#relative-calibration)
 3. [Position calibration](#position-calibration)
 3. [Energy calibration](#energy-calibration)

## Pulser calibration
The macros for calibrating the voltage of the silicon detectors are located in [PulserCal](PulserCal).

## Relative calibration
The macros for gain matching the silicon detector signals relative to one another are located in [RelativeCal](RelativeCal). The calibration methods for the QQQ and SX3 detectors are similar, but they are handeled differently.

### QQQ
The macros for calibrating the QQQ detectors are located in [RelativeCal/QQQ](RelativeCal/QQQ). There are four QQQ detectors in total; detector numbers 0-3.
### SX3
The macros for calibrating the SX3 detectors are located in [RelativeCal/SX3](RelativeCal/SX3). There are 24 SX3 detectors in two rings of 12; detector numbers 4-15 in Ring 1 and detectors 16-27.

## Position calibration
The macros for calibrating the relative position of the SX3 detectors are located in [PositionCal](PositionCal). This calibration is based on the physical geometry of the silicon wafers and is refered to as the geometry calibration throughout the code.

## Energy calibration
The macros for calibrating the energy of the silicon detectors are located in [EnergyCal](EnergyCal).
The energy calibration macro calibrates the back detector signals only.
Therefore, the energy calibration is not dependent on the position calibration.

---

## Directory names
The names of a number of directories have chaned for clarity during the development of the code. For reference, the following commands were used. 
The command `git filter-branch` was used instead of `git mv` to preserve the fold-specific history.

````bash
REWRITE_FROM='SiCal/SiPulserCal'
REWRITE_TO='SiCal/PulserCal'
git filter-branch -f --index-filter "
git read-tree --prefix='$REWRITE_TO'/ \$GIT_COMMIT:'$REWRITE_FROM'
git rm -r --cached --ignore-unmatch '$REWRITE_FROM'
"
````

from [Stack Overflow]: http://stackoverflow.com/questions/42355621/git-filter-branch-moving-a-folder-with-index-filter-does-not-work

The following changes have been made
````bash
REWRITE_FROM='SiCal/SiPulserCal'
REWRITE_TO='SiCal/PulserCal'

REWRITE_FROM='SiCal/AlphaCalibration'
REWRITE_TO='SiCal/EnergyCal'

REWRITE_FROM='SiCal/X3GeometryCal'
REWRITE_TO='SiCal/PositionCal'

REWRITE_FROM='PCCal/PosCal'
REWRITE_TO='PCCal/PositionCal'

REWRITE_FROM='SiCal/SiRelativeCal'
REWRITE_TO='SiCal/RelativeCal'

REWRITE_FROM='SiCal/RelativeCal/QQQ3s'
REWRITE_TO='SiCal/RelativeCal/QQQ'

REWRITE_FROM='SiCal/RelativeCal/SX3s'
REWRITE_TO='SiCal/RelativeCal/SX3'

REWRITE_FROM='anasen_analysis_software'
REWRITE_TO='analysis_software'
````
After each instance of git filter-branch called, the following commands must be used to push the changes to the repository and to pull them to remote clones.
`git push -f` (fork) and `git push -f upstream` (original). Please note that his technique only works if there is no one else developing the repsoitory in parrallel.

To undo the actions of `git filter-branch` use the command `git reset --hard refs/original/refs/heads/master`.

On the remote clones of the repository, the following command must be used to bring the clone up to date without duplicating every commit. 
check the difference with 
`git diff origin` or `git diff origin --diff-filter=r`
ensure that the only changes are renames then 

`git reset --hard origin/master`
