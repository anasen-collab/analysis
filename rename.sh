REWRITE_FROM='PCCal/PosCal'
REWRITE_TO='PCCal/PositionCal'
echo $REWRITE_FROM
echo $REWRITE_TO

git filter-branch -f --index-filter "
        git read-tree --prefix='$REWRITE_TO'/ \$GIT_COMMIT:'$REWRITE_FROM'
        git rm -r --cached '$REWRITE_FROM'
"
