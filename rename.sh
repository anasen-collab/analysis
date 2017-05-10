REWRITE_FROM='SiCal/SiRelativeCal'
REWRITE_TO='SiCal/RelativeCal'
echo $REWRITE_FROM
echo $REWRITE_TO

git filter-branch -f --index-filter "
        git read-tree --prefix='$REWRITE_TO'/ \$GIT_COMMIT:'$REWRITE_FROM'
        git rm -r --cached '$REWRITE_FROM'
"
