#!/bin/env perl
# -*-Perl-*-
#
# $Id$
#
# get the basis family names and atom list from an nwchem formatted library
#
# 

$debug = 0;

$num_argv = @ARGV;
if ($num_argv == 0) {
    &Usage;
    die "fatal error: no file to extract from \n";
}

$tmp_stub = ".tmp.gen_liblist";
foreach $libraryfile (@ARGV) {
# open input file 
    open (FILEIN, $libraryfile) || die "fatal error: Could not open file: $libraryfile\n";
    $temp_file = $tmp_stub . $libraryfile;
    $otemp_file = ">" . $temp_file;
    open (FILETMP, $otemp_file) || die "fatal error: Could not open file: $temp_file for output\n";
    $count_family = 0;
    while (<FILEIN>) {
	if (/^basis/) {
	    if ($debug) {print $_;}
	    $count_family++;
            $string = $_;
	    $string =~ s/\n//;
	    $string =~ s/basis//;
	    $string =~ s/CARTESIAN//;
	    $string =~ s/SPHERICAL//;
	    $string =~ s/\"//g;  #" for hilight in emacs
	    $string =~ s/ //g;
	    if ($debug) {print "string is <$string>\n";}
	    @mystring = split(/_/,$string);
	    $num = @mystring;
	    if ($debug) {print "num is $num mystring is < @mystring >\n";}
	    $atom = @mystring[0];
	    $family = @mystring[1];
	    for ($i=2;$i<$num;$i++){
		$family = $family . " " . @mystring[$i];
	    }
	    if ($debug) {print "Atom::: $atom\n";}
	    if ($debug) {print "Family: $family\n";}
	    print FILETMP "Basis Set:$family:$atom\n";
	}
	if (/^ecp/) {
	    if ($debug) {print $_;}
	    $count_family++;
            $string = $_;
	    $string =~ s/\n//;
	    $string =~ s/ecp//;
	    $string =~ s/CARTESIAN//;
	    $string =~ s/SPHERICAL//;
	    $string =~ s/\"//g;  #" for hilight in emacs
	    $string =~ s/ //g;
	    if ($debug) {print "string is <$string>\n";}
	    @mystring = split(/_/,$string);
	    $num = @mystring;
	    if ($debug) {print "num is $num mystring is < @mystring >\n";}
	    $atom = @mystring[0];
	    $family = @mystring[1];
	    for ($i=2;$i<$num;$i++){
		$family = $family . " " . @mystring[$i];
	    }
	    if ($debug) {print "Atom::: $atom\n";}
	    if ($debug) {print "Family: $family\n";}
	    print FILETMP "ECP:$family:$atom\n";
	}
    }
    close(FILEIN);
    close(FILETMP);
    $debug=0;
    open (FILETMP, $temp_file) || die "fatal error: Could not open file: $temp_file for input\n";
    $count_lines = 0;
    $fams = 0;
    $ats  = 0;
    @famatcnt = (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,00,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,00,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
    @families = ();
    while(<FILETMP>) {
	$count_lines++;
	$line = $_;
	$line =~ s/\n//g;
	@tokens = split(/:/,$line);
	$num = @tokens;
	if ($debug) {print "tokens:$num: @tokens\n";}
	$family = @tokens[0] . " " . "\"". @tokens[1] . "\"";
	$atom = @tokens [2];
#	if ($count_lines>1){
	    $fams = @families;
	    $foundit = 0;
	    for ($i=0;$i<$fams;$i++){
		if ($family eq $families[$i]) {$foundit++;}
	    }
	    if (!($foundit)){
		push(@families,$family);
	    }
#	}
	$fams = @families;
	$fams--;
	$famatcnt[$fams] = $famatcnt[$fams] + 1;
	$ats = $famatcnt[$fams];
	$atoms{$family}[$ats] = $atom;
#	push @{atoms{$family}}, $atom;
#	$ats++;
#	push(@families,$family);
#	%atoms = ($family,$atom);
#	push(@atoms,($family,$atom));
#
#	if ($count_lines > 50) {die "tmp death\n";}
	
	if ($debug) {
	    $fams = @families;
	    if ($debug) {print "fams $fams\n";}
	    for ($i=0;$i<$fams;$i++) {
		$ats = $famatcnt[$i];
		print "$families[$i] \n";
		$count = 0;
		for ($j=0;$j<$ats;$j++) {
		    print " $atoms{$families[$i]}[$j]";
		    $count++;
		    if ($count == 25){
			$count=0;
			print "\n";
		    }
		}
	    }
	}
#	print "families @families\n";
#	print "atoms @atoms\n";
    }
    close(FILETMP);
    unlink($temp_file) or die " could not unlink: $temp_file\n";
    $listfile = ">" . $libraryfile . ".list";
    $famfile  = ">" . $libraryfile . ".fam";
    $texfile  = ">" . $libraryfile . ".tex";
    open (LISTFILE, $listfile) || die "fatal error: Could not open file: $listfile\n";
    open (FAMFILE, $famfile) || die "fatal error: Could not open file: $famfile\n";
    open (TEXFILE, $texfile) || die "fatal error: Could not open file: $texfile\n";
    $fams = @families;
    for ($i=0;$i<$fams;$i++) {
	$atsp = $famatcnt[$i];
	$ats  = $atsp + 1;
	print LISTFILE "\n\n";
	print LISTFILE "$families[$i] (number of atoms $atsp)\n";
	$texline = "\\item " . "$families[$i] (number of atoms $atsp)" . "  \\newline"  ;
	$texline =~ s/Basis Set /Basis Set \\verb\#/;
	$texline =~ s/ECP /ECP \\verb\#/;
	$texline =~ s/\" \(/\" \# \(/ ;   #"
        $texline =~ s/ \#/\#/;
        $texline =~ s/ \#/\#/;
        $texline =~ s/ \#/\#/;
	print TEXFILE  "$texline \n";
	print FAMFILE "$families[$i] (number of atoms $atsp)\n";
	$count = 0;
	for ($j=0;$j<$ats;$j++) {
	    print LISTFILE " $atoms{$families[$i]}[$j]";
	    print TEXFILE " $atoms{$families[$i]}[$j]";
	    $count++;
	    if ($count > 25){
		$count=1;
		print LISTFILE "\n";
		print TEXFILE "\n";
	    }
	}
	print LISTFILE "\n";
	print TEXFILE "\n\n\n";
    }
    close(LISTFILE);
    close(FAMFILE);
    close(TEXFILE);
    $listfile = $libraryfile . ".list";
    $famfile  = $libraryfile . ".fam";
    $texfile  = $libraryfile . ".tex";
    print "files generated: $listfile $famfile $texfile \n"
}

sub Usage
{
    print "gen_liblist.pl filename [filename2] [filename3] [...]\n";
    print "   where filename is an NWChem library file\n";
}
