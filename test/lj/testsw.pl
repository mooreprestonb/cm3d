#!/usr/bin/perl

# perl script to test code

$DEBUG=1;

use Cwd;

sub Help{
  print"perl script to test cm3d lj\n";
  exit(1);
}

sub sysran {
  my $rc = $_[0];
  if($rc==0) {
    # print "successful system call \"",$_[1],"\"\n";
  } elsif ($rc == 256) {
    die "command failed: exit($rc)\n";
  } else {
    print "\"$_[1]\" ran, but exited with signal ",$rc,"\n";
    if($Error){die;}
  }
}

sub stripls {  # strip lines with $_[0] out of $_[1]
  my $pattern = $_[0];
  my $file = $_[1];
  my $ofile = $file."tmp";
  open(OFILE,">$ofile");
  open(FILE,"<$file");
  while (<FILE>) {
    if ( /$pattern/ ) {} #print "Matched $_";
    else {print OFILE $_ ;}
  }
  close(OFILE);
  close(FILE);
  rename $ofile,$file;
}

sub striphl {  # strip the first line out of $_[0]
  my $num = $_[0];
  my $file = $_[1];
  my $ofile = $file."tmp";

  open(OFILE,">$ofile");
  open(FILE,"<$file");
  while (<FILE>) {
    if ($. > $num) {print OFILE $_;}
    # else {print $_;}
  }
  close(OFILE);
  close(FILE);
  rename $ofile,$file;
}

sub stripll {  # strip the last $_[0] lines out of $_[1]
  my $num = $_[0];
  my $file = $_[1];
  my $ofile = $file."tmp";
  my $tnum = 1;
  open(FILE,"<$file");
  while (<FILE>) {$tnum++;} # count lines
  close(FILE);
  open(FILE,"<$file");
  open(OFILE,">$ofile");
  $num = $tnum-$num;
  while (<FILE>) {
    if ($. < $num) {print OFILE $_;}
  } 
  close(OFILE);
  close(FILE);
  rename $ofile,$file;
}

sub sysdiff {
  my $command = "diff $_[0] $_[1] > $_[2]";

  if($DEBUG) {print "executing \"$command\"\n";}
  my $rc = system($command);
  if(!($rc==0)) {
    $nerrors++;
    print "\nThere are differences between \"$_[0]\" and \"$_[1]\"\n";
    if($Error){die;}
  } else {
    unlink("$_[2]");
  # if (eof tmp){print "EOF!\n"; }else{print "junk\n";} #test for empty file
  }
}

# this is the start of the main routines

if($#ARGV==-1){&Help;} # help the user, they forgot the syntax, otherwise...

if(($#ARGV==0) && ($ARGV[0] !~ /^-/)) {
  $code  = shift(@ARGV);
} else {
  foreach $arg (@ARGV)  {
    if ($arg =~ /^-/) {
      if ($arg =~ /h/)  {
	&Help;
      }
    } 
  }
  shift;
  $code  = shift(@ARGV);
}

local $dir = cwd();

$nerrs = 0;
print "Testing script, directory = $dir\n";

$| = 1;  # fflush buffer

#  system("cp -f $dir/../$code $code");
print "\nTesting code=$code\n";

$nerrs += test_ljnve2();
$nerrs += test_ljnve();
$nerrs += test_ljnvt();
$nerrs += test_ljnpt_i();

if ($nerrs > 1){
  die "There were $nerrs errors\n";
} else {
  if ($nerrs == 1){
    die "There was a single error detected!\n";
  }
}

print "There were NO errors detected!\n ALL test PASSED\n";

exit(0);

# testing
sub test_ljnve2() {
  $nerrors = 0;
  print "Testing lj (nve) 2\n";

  $ifile = "lj.input_ljnve.2";
  $ofile = "lj.out"; unlink "$ofile";
  $tfile = "lj.out_nve2.t";
  $dfile = "diff.lj.out";

  $nerrors += &test($ifile, $ofile);
  sysdiff("$ofile","$tfile","$dfile");

  if($nerrors==0) {
    print "Test lj NVE 2 particles.... Passed\n";
    unlink "$dfile"; unlink "$ofile";
  } else {
    die "Test failed with $nerrors \n";
  }
  return $nerrors;
}

# testing
sub test_ljnve() {
  $nerrors = 0;
  print "Testing lj (nve) \n";

  $ifile = "lj.input_ljnve";
  $ofile = "lj.out"; unlink "$ofile";
  $tfile = "lj.out_nve.t";
  $dfile = "diff.lj.out";

  $nerrors += &test($ifile, $ofile);
  sysdiff("$ofile","$tfile","$dfile");

  if($nerrors==0) {
    print "Test lj NVE.... Passed\n";
    unlink "$dfile"; unlink "$ofile";
  } else {
    die "Test failed with $nerrors \n";
  }
  return $nerrors;
}

sub test_ljnvt() {
  $nerrors = 0;
  print "Testing lj (nvt) \n";

  $ifile = "lj.input_ljnvt";
  $ofile = "lj.out"; unlink "$ofile";
  $tfile = "lj.out_nvt.t";
  $dfile = "diff.lj.out";

  $nerrors += &test($ifile, $ofile);
  sysdiff("$ofile","$tfile","$dfile");

  if($nerrors==0) {
    print "Test lj NVT.... Passed\n";
    unlink "$dfile";  unlink "$ofile";
  } else {
    die "Test failed with $nerrors errors\n";
  }
  return $nerrors;
}

sub test_ljnpt_i() {
  $nerrors = 0;
  print "Testing lj (npt_i) \n";

  $ifile = "lj.input_ljnpt_i";
  $ofile = "lj.out"; unlink "$ofile";
  $tfile = "lj.out_npt_i.t";
  $dfile = "diff.lj.out";

  $nerrors += &test($ifile, $ofile);
  sysdiff("$ofile","$tfile","$dfile");

  if($nerrors==0) {
    print "Test lj NPT_I.... Passed\n";
    unlink "$dfile"; unlink "$ofile";
  } else {
    die "Test failed with $nerrors errors\n";
  }
  return $nerrors;
}

sub test($$) {
  my $command;
  my $rc;
  my ($ifile, $ofile) = @_;
  my $efile = "lj.err";

  $command = "$code $ifile $ofile 2> $efile";
  if($DEBUG) {print "executing \"$command\"\n";}
  $rc = system($command);
	
  sysran($rc,$command);

  stripls("^time steps=","$ofile");
  stripls("^Bytes","$ofile");
  stripls("completed","$ofile");
  striphl(1,"$ofile"); # strip head
  #  stripll(1,"$tmpf"); # strip tail

  return $nerrors;
}
