# $Id$

$jobid = $ENV{"PBS_JOBID"};
if ($jobid eq "") {
    exit(1); 
#die("PBS_JOBID not defined\n");
}
exit(1) unless open(LL,"$qstat -f $jobid|");
    $nproc=0;
    $used=0;
while (<LL>) {
    chop;
    ($field, $value) = split(/ :/);
    if ($field =~ /Resource_List.walltime = (\d\d):(\d\d):(\d\d)/){
	$walllimit = $1*3600 + $2*60 + $3;
    } elsif ($field =~ /resources_used.walltime = (\d\d):(\d\d):(\d\d)/){
	$used = $1*3600 + $2*60 + $3;}
}

close(LL);

exit 1 unless (defined($used) &&defined($walllimit));


$left = $walllimit - $used;
#print "wsec = $walllimit used = $used\n";

#print "The job has been running for $used seconds and has $left seconds remaining.\n";

print "$left\n"

