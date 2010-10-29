#!/bin/env perl
# -*-Perl-*-
#
# $Id$
#
# look at NWCHEM.BSLIB and identify lines where 
# END
# END
# exits.
# also checks to make sure that all BASIS blocks end with END
# 

$debug = 0;

$num_argv = @ARGV;
if ($num_argv == 0) {
    &Usage;
    die "fatal error: no file to check \n";
}

foreach $libraryfile (@ARGV) {
# open input file 
    open (FILEIN, $libraryfile) || die "fatal error: Could not open file: $libraryfile\n";
    $count_lines = 0;
    $inblock = 0;
    while (<FILEIN>) {
	$count_lines++;
	if (/^BASIS/) {
	    if ($inblock) {
		print "At line $count_lines BASIS was found inside a block\n";
	    }
	    else {
		$inblock++;
	    }
	}
	if (/^END/) {
	    if (!($inblock)) {
		print "At line $count_lines END was found while not in a block\n";
	    }
	    else {
		$inblock--;
	    }
	}
    }
    close(FILEIN);
    print "checked file $libraryfile with $count_lines lines\n";
}

sub Usage
{
    print "\n\n\n";
    print "USAGE: end2end.pl filename [filename2] [filename3] [...]\n";
    print "  where filename is a NWCHEM.BSLIB file from DFFeller that\n";
    print "  needs to be converted to an NWChem library file\n";
    print "\n\n\n";
}
