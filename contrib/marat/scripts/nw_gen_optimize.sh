#!/bin/csh -f
set sname = `echo ${0}| awk '{print $(NF)}' FS="/"`
# Check that there are some arguments present. 
set args = ($argv[1-])
#if ( $#args == 0 ) then
#goto help
#exit 
#endif
# Process each of the arguments to the script. 
echo $args
while ( $#args > 0 )
  switch ($args[1]) 
  case -charge:	
     shift args 
     set charge = $args[1] 
     shift args 
     breaksw
   case -sim
     shift args 
     set sim = $args[1] 
     shift args 
     breaksw
   case -structure:	# output file
     shift args 
     set structure = $args[1] 
     shift args 
     breaksw
   case -mult:	# output file
     shift args 
     set mult = $args[1] 
     shift args 
     breaksw
   case -o  : # output file
     shift args 
     set out = $args[1] 
     shift args 
     breaksw
   case *:
     goto help
   endsw
end


if(($?out) == 0) then
  echo "please provide output file"
  goto help
endif

if(($?structure) == 0) then
  echo "please provide structure file"
  goto help
endif

if(($?sim) == 0) then
  set sim = `echo $structure | cut -d "." -f1`
endif


if(($?charge) == 0) set charge = 0
if(($?mult) == 0) set mult = 1
echo $charge
echo $mult

cat >> $out  <<INPFILE2
memory total 2000 Mb

start ${sim}


permanent_dir ./perm

scratch_dir   ./data

geometry nocenter noautoz noautosym
load ${structure}
end


basis "ao basis"
 * library "6-31+G*"
 end

dft
 mult ${mult}
 XC b3lyp
 iterations 5000
end

driver
 xyz ${sim}.out
 maxiter 100
end

task dft optimize

INPFILE2
exit
help:
echo "NAME"
echo "    ${sname} -- generates nwchem input file for optimization"
echo "SYNOPSIS"
echo "     -structure xyz or pdb structure (required)"
echo "     -o output file (required)"
echo "     -charge charge (optional defaults to 0) "
echo "     -mult multiplicity (optional defaults to 1)"
echo "     -sim simulation name (optional defaults to structure name)"
exit
