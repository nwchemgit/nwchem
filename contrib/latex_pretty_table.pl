#!/msrc/apps/bin/perl
#
#  latex_pretty_table.pl
#  Hugh L Taylor, Battelle / Pacific Northwest Laboratory 
#  hl_taylor@pnl.gov 
#  Sun Jun 25 1995
#  $Log: not supported by cvs2svn $
# Revision 1.7  1995/06/27  05:04:39  pg511
# Cleanup.  Comments, variable names, help option.
#;;
#;; the following is an elisp keyboard macro to call this function from emacs
#;; place it in your .emacs file
#;;
#;; temporary keyboard macro form of latex-pretty-table
#;;
#(fset 'latex-pretty-table
#   [19 101 110 100 123 116 97 98 117 108 97 114 1 4194336 18 103 105 110 123 116 97 98 117 108 14 1 24 24 21 escape 124 108 97 116 101 120 95 112 114 101 116 116 121 95 116 97 98 108 101 46 112 108 return])

#
#
{ require 'newgetopt.pl';
#
#  Handle the options (not used for the moment)
#
  &init_options;
  &NGetOpt(@all_opts);
  &handle_options;
#
  $iline = -1;

#
# read lines into array, saving the maximum length field for each column
# of the table
#
  while (<>) {
    chop;
    $lines[++$iline] = " ".$_;	  # add lead space
    $lines[$iline] =~ s/^\s+/ /g; # exactly one lead space
    $lines[$iline] =~ s/&/ & /g;  # spaces around ampersands

    @fields = split(/&/,$_);
#
#   last field might include \\ or that and following \hline(s)
#   last field is then split in two
#
    @end_fields = split(/\\\\/,$fields[$#fields]);
    pop(@fields);
    push(@fields,($end_fields[0]));
    
    $multicolumn_skip = 0;
    for ( $jfield=0; $jfield <=$#fields; $jfield++) { # 
      $fields[$jfield] =~ s/^ +/ /g;	# 
      $fields[$jfield] =~ s/ +$/ /g;	# <= 1 leading or trailing space
      $temp_len = length($fields[$jfield]) 
	  + &number_braces_in_string($fields[$jfield]);

      if ($fields[$jfield] =~ /\\multicolumn\{(\d+)/) {	# multicolumn entry
	$mc_columns = $1;
	push(@multicolumn_list,($jfield+$multicolumn_skip,
				$jfield+$multicolumn_skip+$mc_columns-1,
				$temp_len));
	$multicolumn_skip += $mc_columns-1;
      } else {
	$max_field_lengths[$jfield] = $temp_len
	    if $max_field_lengths[$jfield] < $temp_len;
      }
    }
  }
#
#  if a \multicolumn is longer than the combined lengths of the columns
#  it covers, add to the length of the rightmost column
#
  $num_mc = ($#multicolumn_list+1)/3 - 1 ;

  for ( $jfield=0; $jfield <=$#max_field_lengths; $jfield++) {
    $max_field_lengths[$jfield] += 4; # a little extra room
  }

  for ($ii = 0; $ii <= $num_mc; $ii++) {
    $left_column =  shift(@multicolumn_list);
    $right_column = shift(@multicolumn_list);
    $mc_length =  shift(@multicolumn_list);
    $tot_length = 0;
    for ($jj = $left_column; $jj <= $right_column; $jj++) {
      $tot_length += $max_field_lengths[$jj];
    }
    if ($mc_length > $tot_length ) {
      $max_field_lengths[$right_column] += $mc_length - $tot_length;
    }
  }

#
# loop over lines, rewriting them
#
  for ( $iline=0; $iline <= $#lines; $iline++ ) { # last line starting from 0
    $#fields = 0;
    @fields = split(/&/,$lines[$iline]);
#
#   last field might include \\ or that and following \hline
#   field is then split in two
#
    @end_fields = split(/\\\\/,$fields[$#fields]);
    pop(@fields);
    push(@fields,($end_fields[0]));

    $delim = " & ";
    $multicolumn_skip = 0;	# running total of reduction in no. of columns
    $mc_width = 0;

    for ( $jfield=0; $jfield <=$#max_field_lengths; $jfield++) {
      $fields[$jfield] =~ s/^ +/ /;	# 
      $fields[$jfield] =~ s/ +$/ /;	# 1 leading or trailing space

      $width = $max_field_lengths[$jfield+$multicolumn_skip];
      if ($fields[$jfield] =~ /\\multicolumn\{(\d+)/) { # multicolumn field
	$mc_width = $1;
	for ($ii = 1; $ii < $mc_width; $ii++) {
	  $width = $width + $max_field_lengths[$ii+$multicolumn_skip+$jfield];
	}
	$multicolumn_skip += $mc_width-1;
      }

      $format_string = "%-".$width."s";
      $delim = "\\\\" if ($jfield == $#max_field_lengths-$multicolumn_skip);
      $delim = "  "   if ($jfield  > $#max_field_lengths-$multicolumn_skip);

      $outfield = sprintf("$format_string %s", $fields[$jfield], $delim);
#
#     sprintf sometimes yields too many chars (e.g. braces aren't counted)
#
      $extrachar = length($outfield) - $width;

      while ($extrachar-- > 0) {$outfield =~ s/  / /};
      printf("$format_string",$outfield);
    }
    printf(" %s",$end_fields[1]) if $#end_fields == 1;
    print " \n";
  }

  exit;
}
#
#
#
sub init_options {
  @gen_opts = ('h','help');
  @all_opts = ();
  push(@all_opts,@gen_opts);
}
#
#
#
sub handle_options {
  &help && exit if ( $opt_h || $opt_help );
#  $verbose = $opt_v || $opt_verbose ? 1 : 0;
}
#
#
#
sub help {
  print "Program:    latex_pretty_table.pl   A Perl script\n";
  print "Purpose:    Reformats latex table contents into nice columns\n";
  print "            for easy editing and reading.\n";
  print "Usage:      latex_pretty_table.pl < infile > outfile\n";
  print "            Meant to be called from emacs C-u M-|; \n";
  print "            shell command on region, with output replacing region.\n";
}

sub number_braces_in_string {    # $string
  local ($string,@temp) = @_;
  @temp = split(/[\{\}]/,$string);
  $number_braces_in_string = $#temp;
}		
  

