# Directory structure
The following chanages have been made using the following commands. `git filter-branch` was used instead of `git mv` to preserve the fold-specific history.

````bash
REWRITE_FROM='SiCal/SiPulserCal'
REWRITE_TO='SiCal/PulserCal'

REWRITE_FROM='SiCal/AlphaCalibration'
REWRITE_TO='SiCal/EnergyCal'
git filter-branch -f --index-filter "
git read-tree --prefix='$REWRITE_TO'/ \$GIT_COMMIT:'$REWRITE_FROM'
git rm -r --cached '$REWRITE_FROM'
"

````
