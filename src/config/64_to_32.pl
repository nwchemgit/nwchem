#
# $Id$
#
#
#
# perl script to do transliteration from "double" values to "8wrap" values
#
# The script is envoked with the command:
#    perl 64_to_32.pl file1.f [file2.f ...]
#
# This script uses $NWCHEM_TOP/src/config/data.64_to_32 which is a datafile 
# that containes the specific transliterations for "double" to "8wrap" in 
# simple ascii format.  If you need to add a conversion see the comments in 
# the data file.
# Tthis conversion is needed for using 32-bit integer libraries (blas,lapack
# and scalapack) when the rest of the code uses 64-bit integers
# 8wrap files convert from 64-bit to 32-bit
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

sub copy_case() {   # Take case from "from" and apply it to "to" and return that new string
    my ($from, $to) = @_;
    my ($lf, $lt) = map length, @_;
    if ($lt < $lf) { $from = substr $from, 0, $lt }
    else { $from .= substr $to, $lf }
    return uc $to | ($from ^ uc $from);
}                   
 
$debug = 0;
@from = ();
@to   = ();
use File::stat;
use File::Basename;
$data_path = dirname(__FILE__);
if($debug) {print "{$data_path} \n";}

#$data_path = $data_path . "/src/config/data.64_to_32";
$data_path = $data_path . "/data.64_to_32";
if($debug) {print "{$data_path} \n";}
open (DATA,$data_path) || die " unable to open: $data_path \n";
while (<DATA>)
{
    if (/^[^\#]/) {
	if($debug) {print $_;}
	@tokens = split(' ');      
	$num_tokens = @tokens ;
	if($debug){print "tokens: @tokens $#tokens $num_tokens \n";}
	push(@from,$tokens[0]);
	push(@to,  $tokens[1]);
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
if ($debug) { print "arguments: @ARGV\n";}
foreach my $file (@ARGV){
    if ($debug){print "file        : $file\n";}
    my $orgfile = $file;
    $filebak = $file . ".$$" ;
    if ($debug){print "backup file : $filebak\n";}
    rename($file,$filebak);
    $file = '>' . $file;
    my $my_pid = $pid;
    if (open(FILETOFIX,$filebak)) {
#	print " open OK\n";
    }else{
	print "Exiting: Could not open file: $filebak\n";
	kill 15, -$my_pid;
    }
    open(FIXEDFILE,$file) || die "Could not open file: $file\n";
    while (<FILETOFIX>) {
	if ( /^c/ || /^C/ || /^\*/ || /^$/){
	    print FIXEDFILE $_;
	}
	else	{
	    for ($compare = 0; $compare < $num_compare ; $compare++)
	    {
		if (/$from[$compare]/i){

		    $froom = index(uc $_,uc $from[$compare]); # Find where the string starts ie the "from"
		    $toot = &copy_case(substr($_,$froom),$to[$compare]) ; # Generate a "to" that matches the case

		    if (/^[ ]{5}[^\s]/) {
			s/([ ]{5}.)$from[$compare](\W{1})/$1${toot}$2/gi ;
		    }
#this takes care of the tab chars in C
		    if (/^[ \t]/){
			s/(\W{1})$from[$compare](\W{1})/$1${toot}$2/gi ;
		    }
#this takes care of declarations in C
		    if (/^[ \S]/){
			s/(\W{1})$from[$compare](\W{1})/$1${toot}$2/gi ;
		    }
		    if (/^[ \d]/){
			s/(\W{1})$from[$compare](\W{1})/$1${toot}$2/gi ;
		    }
		}
	    }
	    print FIXEDFILE $_;
	}
    }
    close(FIXEDFILE);
    close(FILETOFIX);
    unlink($filebak);
    my $filesize = stat("$orgfile")->size;
    if($filesize eq "0") {
	print("Exiting: converted file $orgfile has size $filesize\n");
	kill 15, -$my_pid;
    }
}


	
	
