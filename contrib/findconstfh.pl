#
# $Id$
#
# perl script searches for un-terminated character constants in include files
#
# The script is envoked with the command:
#    perl findconstfh.pl file1.f [file2.f ...]
#
# 
# Written:  3/19/97
# By:       Ricky A. Kendall
#           High Performance Computational Chemistry Group
#           Theory Modeling and Simulation Program
#           Environmental Molecular Sciences Laboratory
#           Pacific Northwest National Laboratory
#           P.O. Box 999
#           Richland, WA 99352-0999
#           email: ra_kendall@pnl.gov
#
$debug = 0;
foreach $file (@ARGV){
    if ($debug){print "file        : $file\n";}
    open(FIXEDFILE,$file) || die "Could not open file: $file\n";
    $found = 0;
    $lines = 0;
  FOUNDIT: {  
      while (<FIXEDFILE>) {
	  $lines ++;
	  if (/^[cC\*]/){
	      if ($debug) {print "$file:input:$lines: $_";} 
	      $line = $_;
	      if ($debug) {print "$file:vline:$lines: $line";} 
	      $quotes = 0;
	      while (length($line)) {
		  $c = chop($line);
		  if ($c eq "'") {$quotes++;}
		  if ($debug) {print "c: < $c > [$quotes] \n";}
	      }
# 	      $times = 10;
# 	      while ($times) {
# 		  $lenc = length($line);
# 		  if ($c = chop($line))  {
# 		      print "true  character < $c >[$lenc]  $times \n"
# 		      }
# 		  else {
# 		      print "false character < $c >[$lenc]  $times \n"
# 		      }
# 		  $times--;
# 	      }
	      if ($debug && $quotes) {print "number of quotes $quotes\n";}
	      if ($quotes % 2) {
		  print "file: $file at line: $lines is a problem || ";
		  print "text: $_";
	      }
	  }
      }
  }
}
