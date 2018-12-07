# usage csh t.com
cd $WORKDIR
find . -exec touch {} \;
#

