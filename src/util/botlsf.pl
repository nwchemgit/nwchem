# $Id$

$jobid = $ENV{"LSB_JOBID"};
if ($jobid eq "") {
    exit(1); 
#die("LSB_JOBID not defined\n");
}
#exit(1) unless open(LL,"cat /msrc/home/edo/cr1|");
exit(1) unless open(LL,"$bjobs -l $jobid|");
    $nproc=0;
    $used=0;
while (<LL>) {
    chop;
    ($field, $value) = split(/ :/);
    if ($field =~ /(\d+\.\d+) min of */) {
	$walllimit = $1*60;
    }    elsif ($field =~ /Started on (\d+) Hosts*/) {
	$nproc = $1;
    }    elsif ($field =~ / *The CPU time used is (\d+) seconds./) {
	$used = int($1/$nproc*1.05);}
}

close(LL);

exit 1 unless (defined($used) &&defined($walllimit));


$left = $walllimit - $used;
#print "wsec = $walllimit used = $used\n";

#print "The job has been running for $used seconds and has $left seconds remaining.\n";

print "$left\n"

