#/bin/sh

cd lj
perl testsw.pl ../$1
cd ../intracheck
perl testsw.pl ../$1
