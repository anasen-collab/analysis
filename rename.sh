REWRITE_FROM='anasen_analysis_software'
REWRITE_TO='analysis_software'
echo $REWRITE_FROM
echo $REWRITE_TO

git filter-branch -f --index-filter "
        git read-tree --prefix='$REWRITE_TO'/ \$GIT_COMMIT:'$REWRITE_FROM'
        git rm -r --cached --ignore-unmatch '$REWRITE_FROM'
"
