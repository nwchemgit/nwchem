#
# $Id
#
#
#
# perl script to check data for  transliteration from "double" values to "single" values
#
# The script is envoked with the command:
#    perl data.check.pl
#
# This script uses/checks $NWCHEM_TOP/src/config/data.dbl2sngl which is a datafile 
# that containes the specific transliterations for "double" to "single" in 
# simple ascii format.  If you need to add a conversion see the comments in 
# the data file.
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
@from = ();
@to   = ();
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
	@tokens = split(' ');      
	$num_tokens = @tokens ;
	if($debug){print "tokens: @tokens $#tokens $num_tokens \n";}
	if ($tokens[0] ne $tokens[1])  {
	    push(from,$tokens[0]);
	    push(to,  $tokens[1]);
	}
	else {
	    print "to and from tokens are the same:\n";
	    print "from: $tokens[0]\n";
	    print "to  : $tokens[1]\n";
	}
    }
}
close (DATA);
if($debug){
    print "from array @from \n";
    print "to   array @to \n";
}
$num_from = @from;
$num_to   = @to;
if ($num_from != $num_to) {
    print "To and From token count not identical\n";
    print "Number of From tokens: $num_from\n";
    print "number of To   tokens: $num_to\n";
    die "Fatal dbl2sngl error\n";
}
else
{
    if ($debug){
	print "Number of From tokens: $num_from\n";
	print "number of To   tokens: $num_to\n";
    }
    $num_compare = $num_from;
}
