#!/msrc/apps/bin/perl
#
#  view_steps.pl
#  Hugh L Taylor, Battelle / Pacific Northwest Laboratory 
#  <pg511@pnl.gov> 
#  Wed May  3 1995
#  $Log: not supported by cvs2svn $
#Revision 1.3  1995/07/07  19:41:44  pg511
#typo
#
#Revision 1.2  1995/07/07  19:32:44  pg511
#Works for SCF and DFT.
#
#Revision 1.1  1995/05/12  18:23:38  pg511
#Initial revision
#
# Revision 1.1  1995/05/11  04:45:10  pg511
# Initial revision
#
#
#  try constructs like:
#while ( <> ) {
#  $converged = 1 if ( /Requested convergence criteria/);
#  $energy = $1 if ( $converged && //\s*Total energy:\s*(\S+)/);
#}

#
# open files
#
{
  open (INPUT,$ARGV[0]) || die "can't open $ARGV[0]";

  $outname = $ARGV[0].".xyz";
  open (OUTPUT,">$outname") || die "can't open $outname";

  $stepname = $ARGV[0].".steps";
  open (STEPFILE,">$stepname") || die "can't open $outname";

  $Title = 'Your advertising here';
#
# loop over lines of input file
#
  while (<INPUT>) {
    if ( /^\s*No\. of atoms\s*:\s*(\d+)/ ) {
#     found DFT # of atoms
      $natoms = $1;
      print "DFT Calculation Output  natoms: $natoms\n";
      last;
    }
    if ( /^\s*atoms\s*=\s*(\d+)/ ) {
#     found SCF # of atoms
      $natoms = $1;
      print "SCF Calculation Output  natoms: $natoms\n";
      last;
    }
  }

  $ioff_steps = 0;
  $istep = 0;

SEARCHLOOP:  
  while (<INPUT>) {
    if ( /.*Requested convergence criteria on energy and density met.*/ ) {
#
#     converged; looking for total energy
#
      print "found DFT energy convergence flag\n";
    }

    if (/\s*Total energy:\s*(\S+)/) {
      $energy = $1;
      $Title = "Total DFT Energy: $energy Hartrees";
      next SEARCHLOOP;
    }

    if (/\s*Total SCF energy =\s*(\S+)/) {
      $energy = $1;
      $Title = "Total SCF Energy: $energy Hartrees";
      next SEARCHLOOP;
    }

    if ( /^\s+ENERGY GRADIENTS/ ) {
#     skip 3 lines
      $temp = <INPUT>;
      $temp = <INPUT>;
      $temp = <INPUT>;
      print OUTPUT "$natoms\n";
      print OUTPUT "$Title\n";
#     loop over natoms lines of coordinates
      for ($ii=1; $ii<=$natoms; $ii++) {
	$temp = <INPUT>;
#
#                               symbol   x        y       z
	if ($temp =~ /^\s*\d+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+/) {
	  $atomlabel[$ii] = $1;
	  $x[$ii] = $2;
	  $y[$ii] = $3;
	  $z[$ii] = $4;
	  $atomlabel[$ii] =~ s/Ru/Fe/; # xmol does not know what Ru is
#
#         convert to angstroms for xmol
#
	  $x[$ii] = $x[$ii] * .529177249;
	  $y[$ii] = $y[$ii] * .529177249;
	  $z[$ii] = $z[$ii] * .529177249;

#	  print "$atomlabel[$ii] $x[$ii]  $y[$ii]  $z[$ii]\n";
	}
      }
#
#     write coordinates to .xyz file
#
      for ($ii=1; $ii<=$natoms; $ii++) {
	print OUTPUT "$atomlabel[$ii] $x[$ii]  $y[$ii]  $z[$ii]\n";
      }
    } else {
#
#     found step vector
#
      if (/^\s*----* Current step/ ) {
#
#       get step vector; change to include last N steps or all steps,
#       so that direction cosines can be generated
#
	$istep++;
	$temp = <INPUT>;		# read line
	for ($ii=1; $ii<=$natoms; $ii++) {
	  $temp = <INPUT>;		# read line
	  if ($temp =~  /^\s*\d+\s+(\S+)\s+(\S+)\s+(\S+)/) {
	    $stepx[$ii + $ioff_steps] = $1;
	    $stepy[$ii + $ioff_steps] = $2;
	    $stepz[$ii + $ioff_steps] = $3;
	  }
	}
#
#       output step coords to step file
#
	for ($ii=1; $ii<=$natoms; $ii++) {
	  print STEPFILE 
	  "$atomlabel[$ii+$ioff_steps] $stepx[$ii+$ioff_steps]  $stepy[$ii+$ioff_steps]  $stepz[$ii]\n";
	}
	$ioff_steps += $natoms;
	print STEPFILE "\n";
#
#       output norms and direction cosines of the last few steps
#
	$stepnorm = &norm($istep,$natoms,*stepx,*stepy,*stepz);

	$nless = $istep-1;
	for ($iless=1; $iless<=$nless && $iless<$istep; $iless++) {
	  $jstep = $istep - $iless;
	  $dir_cosine[$iless] =
	      &direction_cosine($natoms,$istep,$jstep,*stepx,*stepy,*stepz);
	} 
	print "step: $istep norm: $stepnorm\n";
	for ($ii=1; $ii<=$nless; $ii++) {
	  print "     dot(n-$ii): $dir_cosine[$ii]\n";
	}
	print "\n";
      }	
    }
  }


  close(OUTPUT);
  close(STEPFILE);
#  system("xmol -readFormat xyz $outname &");
#  print "xmol call disabled\n";
  exit;
}
#
#
#


sub norm {			#($istep,$natoms,$stepx,$stepy,$stepz);
  local($istep,$natoms,*stepx,*stepy,*stepz) = @_;
  local($ioff,$norm,$ii);
#
# norm of step $istep
#
  $ioff = ($istep-1)*$natoms;
  $norm = 0.0;
  for ($ii=1; $ii<=$natoms; $ii++) {
    $norm += $stepx[$ioff+$ii]*$stepx[$ioff+$ii];
    $norm += $stepy[$ioff+$ii]*$stepy[$ioff+$ii];
    $norm += $stepz[$ioff+$ii]*$stepz[$ioff+$ii];
  }
  $norm = sqrt($norm);
}

sub dot_product {		# ($N,*vec1,*vec2)
  local($N,*vec1,*vec2);
  local($dotval,$ii);

  $dotval = 0.0;
  for ($ii=1; $ii<=$N; $ii++) {
    $dotval += $vec1[$ii] * $vec2[$ii]
  }
  $dotval;
}

sub direction_cosine {		#  ($natoms,$istep,$jstep,stepx,$stepy,$stepz);
  local($natoms,$istep,$jstep,*stepx,*stepy,*stepz) = @_;
  local($ioff,$joff,$temp,$normij);

  $ioff = ($istep-1)*$natoms;
  $joff = ($jstep-1)*$natoms;
  $temp = 0.0;

  for ($ii=1; $ii<=$natoms; $ii++) {
    $temp += $stepx[$ioff+$ii]*$stepx[$joff+$ii];
    $temp += $stepy[$ioff+$ii]*$stepy[$joff+$ii];
    $temp += $stepz[$ioff+$ii]*$stepz[$joff+$ii];
  }

  $normij = &norm($istep,$natoms,*stepx,*stepy,*stepz) *
      &norm($jstep,$natoms,*stepx,*stepy,*stepz);

  if ($normij != 0.0) {
    $temp = $temp/($normij);
  } else {
    $temp;
  }
}










sub init_options {
  @gen_opts = ('h','help','v','verbose');
  @all_opts = ();
  push(@all_opts,@gen_opts);
}
#
#
#
sub handle_options {
  &help && exit if ( $opt_h || $opt_help );
  $verbose = $opt_v || $opt_verbose ? 1 : 0;
}
#
#
#
sub help {
  print "Program:    stepper_to_xyz.pl   A Perl script\n";
  print "Purpose:    \n";
  print "Usage:      stepper_to_xyz.pl infile [options] > outfile\n";
  print "Options:    -h        this help\n";
  print "            -v        verbose\n";
  print "            -verbose  verbose\n";
}


#
# search down to desire
#
#
sub downto {			# might make this a logical function
  local(*FILE,$string) = @_;
  if (/$string/) {
    $result = 0;
    return;
  }
  while(<FILE>) {  
    if (/$string/) {
      $result = 0;
      return;
    }
  }
  $result = 1;			# nonzero for failed search
}
