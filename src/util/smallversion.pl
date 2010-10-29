#!/bin/env perl
# -*-Perl-*-
#
# $Id$
# 
# smallversion.pl is a perl script that breaks up a large util_version.F file
# into several subroutines for compilers that run out of internal memory
# trying to compile a large subroutine.
#
# The goal is to preserve the order of the write statements. 
#
# Ricky Kendall, May 1998
# 
@INC = ("/dfs/apps/perl/lib","/usr/lib/perl5", "/usr/local/lib/perl", "/usr/local/perl/lib");
use File::Copy;

$debug = 0;

if (!(-e 'util_version.F')) {die "util_version.F does not exist";}

copy("util_version.F","big_util_version.F");
unlink "util_version.F" || die "could not delete util_version.F";
open(FILE_BIG,"big_util_version.F");
open(FILE,">util_version.F");
$sub_cnt = 0;
$sub_lines = 0;
$sub_max_lines = 200;
$sub_open = 0;
while (<FILE_BIG>){
    if (/write/){
	if ($sub_lines > $sub_max_lines) {
	    $sub_lines = 0;
	    print FILE "      endif\n";
	    print FILE "      end ! subroutine $sub_name\n";
            print FILE "*----------------------------------------------------------------------\n";
	    print FILE "\n\n";
	    $sub_open = 0;
	}
	if (!($sub_open)) {
	    $sub_open = 1;
	    $sub_cnt++;
	    print " $sub_cnt .. ";
	    $sub_name = "util_ver_" . $sub_cnt;
            print FILE "*----------------------------------------------------------------------\n";
            print FILE "*                $sub_name \n";
            print FILE "*----------------------------------------------------------------------\n";
	    print FILE "      subroutine $sub_name()\n";
	    print FILE "      implicit none\n";
	    print FILE "#include \"global.fh\"\n";
	    print FILE "      if (ga_nodeid().eq.0) then\n";
	}
	$sub_lines++;
	print FILE $_;
    }
}
if ($sub_open) {
    print FILE "      endif\n";
    print FILE "      end ! subroutine $sub_name\n";
    print FILE "*----------------------------------------------------------------------\n";
    print FILE "\n\n";
}
print FILE "      subroutine util_version()\n";
print FILE "      implicit none\n";
print FILE "#include \"global.fh\"\n";
print FILE "      if (ga_nodeid().eq.0) then\n";
$count = 0;
while ($count < $sub_cnt){
    $count++;
    $sub_name = "util_ver_" . $count;
    print FILE "      call $sub_name()\n";
}
print FILE "      call util_flush(6)\n";
print FILE "      endif\n";
print FILE "      call ga_sync()\n";
print FILE "      end\n";
print "\n";
unlink "big_util_version.F" || die "could not delete big_util_version.F";
exit 0;
