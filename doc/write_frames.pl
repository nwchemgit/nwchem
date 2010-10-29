#!/usr/bin/perl
#########!/msrc/apps/bin/perl
#
# New frames code to eat latex2html output and write a frames document
#
# Ricky A. Kendall
#
# 3/19/98
#
# $Id$
#
# remove nwchem banner stuff 3/23/98
#
###@INC = ("/msrc/apps/perl-5.005/lib/5.00502","/usr/lib/perl5");
require File::Copy;

if (($ARGV[0] eq "") || ($ARGV[1] eq "")) {
    print "Usage: write_frames.pl document_name title_string\n";
    exit(1);
}
else {
    $document = $ARGV[0];
    $title    = $ARGV[1];
}
# figure_out where directory is:
# 
$okay = 0;
if (-d "$document") {$okay++;}
if (-e "$document/index.html"){$okay++;}
if (-e "$document/$document.html"){$okay++;}
#print "okay is $okay \n";
if ($okay != 3) {
    print "write_frames.pl: you must be running in the latex source directory\n";
    print "                 where the latex2html directory $document exists\n\n";
    exit(1);
}
if ($document eq "user") {
    $bodystring = "<BODY BGCOLOR=\"#FFFFFF\">\n";
}
elsif ($document eq "prog") {
    $bodystring = "<BODY BGCOLOR=\"#FFFFFF\">\n";
}
else{
    $bodystring = "<BODY>\n";
}
# write index.html (frames document)
if (!(open(FINDEX,">$document/index.html"))){
    die "write_frames.pl: could not open index.html for writing\n";
}
print FINDEX "<HTML>\n<HEAD>\n<TITLE> $title </TITLE>\n</HEAD>\n";
print FINDEX "<FRAMESET COLS=\"2*,5*\">\n";
print FINDEX "    <FRAMESET ROWS=\"3*,5*\">\n";
print FINDEX "        <FRAME SRC=\"$document.search.html\">\n";
print FINDEX "        <FRAME SRC=\"contents.html\">\n";
print FINDEX "    </FRAMESET>\n";
#-noban#print FINDEX "    <FRAMESET ROWS=\"2*,10*\">\n";
#-noban#print FINDEX "        <FRAME SRC=\"banner.html\">\n";
print FINDEX "    <FRAMESET>\n";
print FINDEX "        <FRAME SRC=\"$document.html\" NAME=\"main\">\n";
print FINDEX "    </FRAMESET>\n";
print FINDEX "</FRAMESET>\n";
print FINDEX "</HTML>\n";
close(FINDEX);
#-noban## write banner.html
#-noban#if (!(open(FBAN,">$document/banner.html"))){
#-noban#    die "write_frames.pl: could not open banner.html for writing\n";
#-noban#}
#-noban#print FBAN "<HTML>\n<HEAD>\n<TITLE> NWChem logo </TITLE>\n</HEAD>\n";
#-noban#print FBAN "<BODY BGCOLOR=\"aqua\">\n<CENTER>\n<P>\n";
#-noban#print FBAN "<IMG SRC=\"/docs/nwchem/nwchem_logo.gif\" ";
#-noban#print FBAN " ALT=\"NWChem - computational chemistry on parallel computers\"> \n";
#-noban#print FBAN "</P>\n<HR>\n</CENTER>\n</BODY>\n</HTML>\n";
#-noban#close(FBAN);
#
# write $document_search.html
if (!(open(FHTMLSEARCH,">$document/$document.search.html"))){
    die "write_frames.pl: could not open $document/$document.search.html for writing\n";
}
if (!(open(FSEARCH,"$document.search"))){
    die "write_frames.pl: could not open $document.search for reading(1)\n";
}
print FHTMLSEARCH "<HTML>\n<HEAD>\n<TITLE> $title Search Window</TITLE>\n";
print FHTMLSEARCH "<BASE TARGET=\"main\">\n";
print FHTMLSEARCH "</HEAD>\n";
print FHTMLSEARCH "$bodystring";
while (<FSEARCH>){
    print FHTMLSEARCH $_;
}
print FHTMLSEARCH "<A TARGET=\"_TOP\" HREF=\"$document.html\"> No Frames Version </A><HR>";
print FHTMLSEARCH "\n</BODY>\n</HTML>\n";
close(FHTMLSEARCH);
close(FSEARCH);
# write contents.html  (Contents are in a block that 
#  starts with: <!--Table of Contents--> and 
#  ends with  : <!--End of Table of Contents-->
# this is usually in node2.html.  This is assumed
#
if (!(open(FHTMLCONT,">$document/contents.html"))){
    die "write_frames.pl: could not open $document/contents.html for writing\n";
}
if ($document eq "user") {
  if (!(open(FCONT,"$document/node2.html"))){
    die "write_frames.pl: could not open $document/node2.html for reading(2)\n";
  }
}
if ($document eq "prog") {
  if (!(open(FCONT,"$document/node4.html"))){
    die "write_frames.pl: could not open $document/node2.html for reading(2)\n";
  }
}
print FHTMLCONT "<HTML>\n<HEAD>\n<TITLE><B>Contents of $title </B></TITLE>\n";
print FHTMLCONT "<BASE TARGET=\"main\">\n</HEAD>\n";
print FHTMLCONT "<BR>\n<B>Contents of<br>&nbsp;&nbsp;&nbsp; $title </B><BR>\n";
print FHTMLCONT "$bodystring";
$printit = 0;
$printdone = 0;
while(<FCONT>){
    if (!($printdone)){
	if (/<!--Table of Contents-->/) {$printit = 1;}
	if (/<!--End of Table of Contents-->/) {$printit = 0;$printdone=1;}
	if ($printit) {print FHTMLCONT $_;}
    }
}
print FHTMLCONT "</BODY>\n</HTML>\n";
close(FHTMLCONT);
close(FCONT);
# fix No Title
chdir $document;
$tmp = "tmp_$$";
#print "tmp is $tmp \n";
foreach $file (`ls -1 *.html`) {
    chop($file);
    if (-e "$tmp") {unlink $tmp || die "could not delete $tmp\n";}
#    print "file: <: $file :>\n";
#    print "tmp is <: $tmp :> \n";
    $result_cp = ` /bin/cp $file $tmp `;
# perl copy seems to be broke ??
#    copy($file,$tmp);
#    $result_ls = `ls -l $tmp `;
#    print $result_ls;
#    exit;
    $Ofile = '>' . $file;
    open(FO,$Ofile) || die " could not open $Ofile for writting \n";
    open(FI,$tmp)   || die " could not open $tmp for reading(3) \n";
    while (<FI>) {
	if (/No Title/) {
	    $line = chop();
	    $line =~ s/No Title/$title/; 
	    print FO "$line \n";
	}
	else{
	    print FO $_;
	}
    }
    close(FO);close(FI);
}
if (-e "$tmp") {unlink $tmp || die "could not delete $tmp\n";}
exit(0);
