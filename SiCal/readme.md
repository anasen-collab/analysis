# Directory names
The names of a number of directories have chaned using the following commands. 
The command `git filter-branch` was used instead of `git mv` to preserve the fold-specific history.

````bash
REWRITE_FROM='SiCal/SiPulserCal'
REWRITE_TO='SiCal/PulserCal'
git filter-branch -f --index-filter "
git read-tree --prefix='$REWRITE_TO'/ \$GIT_COMMIT:'$REWRITE_FROM'
git rm -r --cached '$REWRITE_FROM'
"
````

The following changes have been made
````bash
REWRITE_FROM='SiCal/SiPulserCal'
REWRITE_TO='SiCal/PulserCal'

REWRITE_FROM='SiCal/AlphaCalibration'
REWRITE_TO='SiCal/EnergyCal'

REWRITE_FROM='SiCal/SiRelativeCal'
REWRITE_TO='SiCal/RelativeCal'
````
After each instance of git filter-branch called, the following commands must be used to push the changes to the repository and to pull them to remote clones.
`git push -f` (fork) and `git push -f upstream` (original). Please note that his technique only works if there is no one else developing the repsoitory in parrallel. On the remote clones of the repository, the following command must be used to bring the clone up to date without duplicating every commit. `git reset --hard origin/master`
