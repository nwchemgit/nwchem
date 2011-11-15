#
# $Id$
#
# perl script searches for both "double" and "single" values of the 
# blas and lapack routines reporting only those files that have a 
# recognized routine.  This script also prints the matching lines
# so you may see the text that must be transformed for "single" 
# precision computers.
#
# The script is envoked with the command:
#    perl showblas.pl file1.f [file2.f ...]
#
# This script uses $NWCHEM_TOP/src/config/data.dbl2sngl which is a datafile 
# that containes the specific transliterations for "double" to "single" in 
# simple ascii format.  If you need to add a conversion see the comments in 
# the data file.
#
# 
# Written:  3/14/97
# By:       Ricky A. Kendall
#           High Performance Computational Chemistry Group
#           Theory Modeling and Simulation Program
#           Environmental Molecular Sciences Laboratory
#           Pacific Northwest National Laboratory
#           P.O. Box 999
#           Richland, WA 99352-0999
#
$debug = 0;
@tokens = ();
$data_path = $ENV{'NWCHEM_TOP'} ;
if ($data_path eq "") {
    print "Error: environment variable NWCHEM_TOP is not set\n";
    print "dbl2sngl: Fatal error\n" ;
    exit 1;
}
if($debug) {print "{$data_path} \n";}
$data_path = $data_path . "/src/config/data.dbl2sngl";
if($debug) {print "{$data_path} \n";}
open (DATA,$data_path) || die " unable to open: $data_path \n";
while (<DATA>)
{
    if (/^[^\#]/) {
	if($debug) {print $_;}
	@newtokens = split(' ');      
	$num_tokens = @newtokens ;
	if($debug){print "tokens: @newtokens $#newtokens $num_tokens \n";}
	push(tokens,@newtokens);
    }
}
close (DATA);
$num_tokens = @tokens;
if ($debug){
    print "tokens array @tokens \n";
    print "number of tokens: $num_tokens\n";
}
if ($debug) { print "arguments: @ARGV\n";}
@found_files = ();
foreach $file (@ARGV){
    if ($debug){print "file        : $file\n";}
    open(FIXEDFILE,$file) || die "Could not open file: $file\n";
    $found = 0;
    $lines = 0;
  FOUNDIT: {  
      while (<FIXEDFILE>) {
	  $lines ++;
	  if (/^[ \d]/){
	      $itok = 0;
	      while ($itok < $num_tokens )
	      {
		  if (/(\W{1})$tokens[$itok](\W{1})/i) {
		      $found++;
		      if ($debug) {print "token: $tokens[$itok]\n";}
		      print "$file: $_";
		  }
		  $itok++;
	      }
	  }
      }
  }
    if ($found) {push(found_files,$file);}
    close(FIXEDFILE);
}
$num_found_files = @found_files;
if ($num_found_files) {
    print "showblas: found $num_found_files files: @found_files\n";
}
