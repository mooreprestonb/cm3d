#!/usr/bin/perl -W

#Check for significant differences between lines!

# read in diff file
$name = $ARGV[0];

open(FILE,$name) || die "Couldn't open file: $name\n";
@inp = <FILE>;
close(FILE);

$nsig = 0; # number of significant difference 

#read in lines

@line1 = split(" ",$inp[1]);
@line2 = split(" ",$inp[3]);

#print @line1,@line2;
# check for numbers

$icol = 0;
foreach $col1 (@line1){
  if ($col1 =~ /[0-9]/ ){ # number
    $diff = abs($col-$line2[$icol]);
    print "diff = $diff\n";
  # is absolute difference small?

  #yes then skip.

  #no then count sigdiffse
  }
  ++$icol;
}

exit($nsig);



