# Directory names
The names of a number of directories have chaned using the following commands. 
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
````
After each instance of git filter-branch called, the following commands must be used to push the changes to the repository and to pull them to remote clones.
`git push -f` (fork) and `git push -f upstream` (original). Please note that his technique only works if there is no one else developing the repsoitory in parrallel.

To undo the actions of `git filter-branch` use the command `git reset --hard refs/original/refs/heads/master`.

On the remote clones of the repository, the following command must be used to bring the clone up to date without duplicating every commit. `git reset --hard origin/master`
