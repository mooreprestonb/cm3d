#!/usr/bin/perl -W

# perl script to test code of intra molecular interactions 

$DEBUG=1;

use Cwd;

sub Help{
  print"perl script to test cm3d intramolecular\n";
  print"argv[0] <cm3d>\n";
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

#testing routine runs command and strips different lines out of the file

sub test($$) {
  my $command;
  my $rc;
  my ($ifile, $ofile) = @_;
  my $efile = "out.err";

  $command = "$code $ifile > $ofile 2> $efile";
  if($DEBUG) {print "executing \"$command\"\n";}
  $rc = system($command);
	
  sysran($rc,$command);

  stripls("^time steps=","$ofile");
  stripls("^Bytes","$ofile");
  stripls("completed","$ofile");
  stripls("Searching","$ofile");
  striphl(1,"$ofile"); # strip head
  #  stripll(1,"$tmpf"); # strip tail

  return $nerrors;
}

#test bonding interactions

# testing diatomic harmonic
sub test_diharm() {
  $nerrors = 0;
  print "Testing di 2\n";

  $ifile = "di.input_harm.2";
  $ofile = "di.out"; unlink "$ofile";
  $tfile = "di.out_diharm.t";
  $dfile = "diff.lj.out";

  $nerrors += &test($ifile, $ofile);
  sysdiff("$ofile","$tfile","$dfile");

  if($nerrors==0) {
    print "Test di NVE 2 particles.... Passed\n";
    unlink "$dfile"; unlink "$ofile";
  } else {
    die "Test failed with $nerrors \n";
  }
  return $nerrors;
}

# testing diatomic morse
sub test_dimorse() {
  $nerrors = 0;
  print "Testing di 2\n";

  $ifile = "di.input_morse.2";
  $ofile = "di.out"; unlink "$ofile";
  $tfile = "di.out_dimorse.t";
  $dfile = "diff.lj.out";

  $nerrors += &test($ifile, $ofile);
  sysdiff("$ofile","$tfile","$dfile");

  if($nerrors==0) {
    print "Test di NVE 2 connect with morse .... Passed\n";
    unlink "$dfile"; unlink "$ofile";
  } else {
    die "Test failed with $nerrors \n";
  }
  return $nerrors;
}
# testing diatomic quartic and harm
sub test_diquartich() {
  $nerrors = 0;
  print "Testing di 2\n";

  $ifile = "di.input_quartic.2.h";
  $ofile = "di.out"; unlink "$ofile";
  $tfile = "di.out_diharm.t";
  $dfile = "diff.lj.out";

  $nerrors += &test($ifile, $ofile);

  sysdiff("$ofile","$tfile","$dfile");

  if($nerrors==0) {
    print "Test di NVE 2 particles.... Passed\n";
    unlink "$dfile"; unlink "$ofile";
  } else {
    die "Test failed with $nerrors \n";
  }
  return $nerrors;
}

# testing diatomic quartic
sub test_diquartic() {
  $nerrors = 0;
  print "Testing di 2\n";

  $ifile = "di.input_quartic.2";
  $ofile = "di.out"; unlink "$ofile";
  $tfile = "di.out_diquartic.t";
  $dfile = "diff.lj.out";

  $nerrors += &test($ifile, $ofile);
  sysdiff("$ofile","$tfile","$dfile");

  if($nerrors==0) {
    print "Test di NVE 2 particles.... Passed\n";
    unlink "$dfile"; unlink "$ofile";
  } else {
    die "Test failed with $nerrors \n";
  }
  return $nerrors;
}

# Testing Bending routines....

# testing triatomic harmonic
sub test_triharm() {
  $nerrors = 0;
  print "Testing tri 3\n";

  $ifile = "tri.input_harm.3";
  $ofile = "tri.out"; unlink "$ofile";
  $tfile = "tri.out_triharm.t";
  $dfile = "diff.lj.out";

  $nerrors += &test($ifile, $ofile);
  sysdiff("$ofile","$tfile","$dfile");

  if($nerrors==0) {
    print "Test tri NVE harmoic bonds+bends 3 particles.... Passed\n";
    unlink "$dfile"; unlink "$ofile";
  } else {
    die "Test failed with $nerrors \n";
  }
  return $nerrors;
}

# testing triatomic quartic (same as harm) bend
sub test_triquartich() {
  $nerrors = 0;
  print "Testing tri 3\n";

  $ifile = "tri.input_quartich.3";
  $ofile = "tri.out"; unlink "$ofile";
  $tfile = "tri.out_triharm.t";
  $dfile = "diff.lj.out";

  $nerrors += &test($ifile, $ofile);
  sysdiff("$ofile","$tfile","$dfile");

  if($nerrors==0) {
    print "Test tri NVE harm bonds+ quartic 3 particles.... Passed\n";
    unlink "$dfile"; unlink "$ofile";
  } else {
    die "Test failed with $nerrors \n";
  }
  return $nerrors;
}

# testing triatomic cosine bend
sub test_tricosine() {
  $nerrors = 0;
  print "Testing tri 3\n";

  $ifile = "tri.input_cosine.3";
  $ofile = "tri.out"; unlink "$ofile";
  $tfile = "tri.out_tricosine.t";
  $dfile = "diff.lj.out";

  $nerrors += &test($ifile, $ofile);
  sysdiff("$ofile","$tfile","$dfile");

  if($nerrors==0) {
    print "Test tri NVE harm bonds+ cosinebends 3 particles.... Passed\n";
    unlink "$dfile"; unlink "$ofile";
  } else {
    die "Test failed with $nerrors \n";
  }
  return $nerrors;
}

# testing triatomic quartic (same as harm) bend
sub test_triquartic() {
  $nerrors = 0;
  print "Testing tri 3\n";

  $ifile = "tri.input_quartic.3";
  $ofile = "tri.out"; unlink "$ofile";
  $tfile = "tri.out_triquartic.t";
  $dfile = "diff.lj.out";

  $nerrors += &test($ifile, $ofile);
  sysdiff("$ofile","$tfile","$dfile");

  if($nerrors==0) {
    print "Test tri NVE harm bonds+ quartic 3 particles.... Passed\n";
    unlink "$dfile"; unlink "$ofile";
  } else {
    die "Test failed with $nerrors \n";
  }
  return $nerrors;
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

$nerrs += test_diharm();
$nerrs += test_dimorse();
$nerrs += test_diquartich();
$nerrs += test_diquartic();
$nerrs += test_triharm();
$nerrs += test_tricosine();
$nerrs += test_triquartich();
$nerrs += test_triquartic();

if ($nerrs > 1){
  die "There were $nerrs errors\n";
} else {
  if ($nerrs == 1){
    die "There was a single error detected!\n";
  }
}

print "There were NO errors detected!\n ALL test PASSED\n";

exit(0);
