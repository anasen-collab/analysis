REWRITE_FROM='SiCal/RelativeCal/QQQ3s'
REWRITE_TO='SiCal/RelativeCal/QQQ'
echo $REWRITE_FROM
echo $REWRITE_TO

git filter-branch -f --index-filter "
        git read-tree --prefix='$REWRITE_TO'/ \$GIT_COMMIT:'$REWRITE_FROM'
        git rm -r --cached --ignore-unmatch '$REWRITE_FROM'
"
