#!/usr/local/bin/perl5
#
# $Id$
#

# ON THE IBM SP DETERMINE THE TIME LEFT TO A LL BATCH JOB

# PRINT THE TIME LEFT IN SEC ON STDOUT AND EXIT(0)
# ON FAILURE EXIT(1) AND PRINT NON-INTEGER MESSAGE.


sub walltosec {

    use integer;
 $date = shift(@_); 
    @days = split(/[ +:]/,$date);
    if($days[0] eq '00') {
    @fields = split(/[ +:]/,$date);
    } else{
    @fields = split(/[ +:]/,$days[1]);
    }
#    print "days |$days[0]| fld1 $fields[0]  fld2 $fields[1] fld3 $fields[2]\n";
    $walltosec = $days[0]*24*3600+$fields[0]*3600+$fields[1]*60+$fields[2]; 
}

sub datetosec {

    use integer;

    my $datetosec, $day, $month, $daynum, $hour, $min, $sec, $zone, $year, $dayinyear, $secinyear, $date, $dayssince1997;

    $date = shift(@_);

    $precdays{Jan} = 0;
    $precdays{Feb} = 31;
    $precdays{Mar} = 59;
    $precdays{Apr} = 90;
    $precdays{May} = 120;
    $precdays{Jun} = 151;
    $precdays{Jul} = 181;
    $precdays{Aug} = 212;
    $precdays{Sep} = 243;
    $precdays{Oct} = 273;
    $precdays{Nov} = 304;
    $precdays{Dec} = 334;
    
    @fields = split(/[ +:]/,$date);
#    ($day, $month, $daynum, $hour, $min, $sec, $zone, $year) = split(/[ :]/,$date);

    for ($i=0; $i<@fields; $i++) {
#       print "field $i = $fields[$i]\n";
       if ($fields[$i] =~ /[ :]/ || $fields[$i] eq "") {
#         print "shifting\n";
         for ($j=($i+1); $j<@fields; $j++) {
           $fields[$j-1] = $fields[$j];
         }
         $fields[@fields-1] = ' ';
       }
    }
        
   ($day, $month, $daynum, $hour, $min, $sec, $zone, $year) = @fields;
# Handle case when no zone info is included (from LL output)

    if ($zone =~ /\d+/) {
	$year = $zone;
    }
    
# print "day=$day month=$month day#=$daynum hour=$hour min=$min sec=$sec zone=$zone year=$year\n";

    $dayssince1997 = ($year - 1997)*365 + ($year - 1997)/4 ;

    $dayinyear = $precdays{$month} + $daynum;
    
    if (($dayinyear >= 59) && (($year % 4) > 0)) {$dayinyear += 1};
    
    $secinyear = $sec + 60*($min + 60*($hour + 24*($dayssince1997 + $dayinyear)));
    
#    print "seconds = $secinyear\n";

    $datetosec = $secinyear;
}


$jobid = $ENV{"LOADL_JOB_NAME"};

if ($jobid eq "") {
    exit(1); 
#die("LOADL_JOB_NAME not defined\n");
}

#die("llq failed\n") unless open(LL,"llq -l $jobid|");
exit(1) unless open(LL,"llq -l $jobid|");
#exit(1) unless open(LL,"cat /u/edo/crap|");

while (<LL>) {
    chop;
    ($field, $value) = split(/: /);
    if ($field =~ /[ ]*Dispatch Time/) {
	$dispatch = $value;
    }
    elsif ($field =~ /Wall Clk Hard Limit/) {
        if($value =~ /(\d+) seconds/){
        $walllimit=$1;
        }
    }
}

close(LL);

#die("Unable to determine dispatch/limit\n") unless (defined($dispatch) &&defined($walllimit));
exit 1 unless (defined($dispatch) &&defined($walllimit));

#print "Dispatch date = '$dispatch'; walllimit = $walllimit\n";

$now = `date`;
#$now = "Thu Oct 24 00:11:01 2002           ";
chop($now);

#print("now = $now\n");
$used = &datetosec($now) - &datetosec($dispatch);
$left = $walllimit - $used;
#print "  walllimit $walllimit \n";
#print "used = $used\n";

#print "The job has been running for $used seconds and has $left seconds remaining.\n";

print "$left\n"

