#!/msrc/apps/bin/perl5
###########################
# latex2HTML wrapper script
# Pawel Wolinski
# pwolinsk@comp.uark.edu
# 5/21/97
###########################

#
# This script takes the output of htmlize (munged latex2html)
# and adds the frame based indexing
#


$input_dir =".";
$output_dir=".";
$title="NWChem User Documentation";

if ($ARGV[0] eq "") {&print_usage;exit(0);} else {$start_file=$ARGV[0];}
if ($ARGV[1] ne "") {$title=$ARGV[1];}

&make_contents_cover($input_dir,$output_dir,$start_file);
&make_manuals($output_dir,$title);





#-----------------------subroutines-------------------
sub make_contents_cover {
local ($in_dir,$out_dir,$start_f)=@_;
$r1=<<EOT;
.*<IMG ALIGN=BOTTOM ALT="next".*
<B>.*Next:</B>.*
<B>.*Up:</B>.*
<B>.*Previous:</B>.*</A>*.
EOT

$r2=<<EOT2;
<script>
var ontop=true;
function focus_me() {
   if (ontop) {self.focus();}
}
function toggle_ontop() {
   if (ontop==true) ontop=false;
   else ontop=true;
   focus_me();
}
function index() {
   if (self.name=="contentsWin") {
      self.opener.location="index.html";
      self.close();
   }
}
function windows_man() {
   if (self.name=="contents_frame") {
      top.location="windows_man.html";
   }
}
</script>
EOT2

$r3=<<EOT5;
<table>
<tr>
<td><font size=-1><i><input type=\"button\" value=\"frames\" onClick=\"index()\"></i></font></td>
<td><font size=-1><i><input type=\"button\" value=\"windows\" onClick=\"windows_man()\"></i></font></td>
<td width=100%></td>
<td><font size=-1><i><input type=\"button\" value=\"on top\" onClick=\"toggle_ontop()\"></i></font></td>
</tr>
</table>
EOT5
   if (!(open(F,"$in_dir/node2.html"))) {
      print stderr "make_manual: Cannot find the contents file \"node2.html\" - terminating!\n\n";
      exit(0);
   }
   close(F);

   if (!(open(F,"$in_dir/$start_f"))) {
      print stderr "make_manual: Cannot find the start file \"$start_f\" - terminating!\n\n";
      exit(0);
   }
   close(F);

   $contents=`cat $in_dir/node2.html`;
   $cover   =`cat $in_dir/$start_f`;

   $contents=~s/$r1//g;
   $contents=~s/<HR>//g;
   $contents=~s/<BR>  <P>//ig;

   $cover   =~s/<H2>*.<\/H2>//g;
   $cover   =~s/<UL>.*<\/UL>//gs;
   $cover   =~s/<HR>//g;

   if ($cover=~/<A HREF.*>\s*Compressed postscript<\/A>/i) 
      {$comp_postscript=$&;}
   if ($cover=~/<A HREF.*>\s*Postscript<\/A>/i) 
      {$postscript=$&;}
   if ($cover=~/<ISINDEX.*>/i)
      {$search=$&;}

   $cover   =~s/<P>.*Compressed Postscript<\/A>//ig;
   $cover   =~s/<P>.*Postscript<\/A>//ig;
   $cover   =~s/<ISINDEX.*>//gi;

   $contents=~s/<\/HEAD>/<base target=\"main\">\n$r2<\/HEAD>/g;
   $contents=~s/<BODY/<BODY onBlur=\"focus_me()\"/;
   $contents=~s/<H2>/<font size=-1><i><form>$r3$comp_postscript $postscript<\/form><\/i><\/font>\n$search\n<H2>/g;
  
   open(F,">$out_dir/contents.html");
   print F $contents;
   close(F);

   open(F,">$out_dir/cover.html");
   print F $cover;
   close(F);
}




sub print_usage {
print <<EOT2;
usage: make_manual.pl [start file] ["title"]

    start file:     root file of the latex2html tree
    title     :     (optional) title displayed in browser's 
                    header bar; default="NWChem User Documentation"
EOT2
}



sub make_manuals {
local($out_dir,$tit)=@_;

$frames_index=<<EOT3;
<html>
<head>
<title>$tit</title>
<script>
<!--
window.name="something";
//-->
</script>
</head>
<frameset cols="300,*" border=1>
   <frame name="contents_frame" scrolling="yes" src="contents.html">
   <frame name="main" src="cover.html">
</frameset>
</html>
EOT3


$cover=`cat cover.html`;
if ($cover=~/<body.*<\/body>/is) {
   $cover=$&;
   $cover=~s/<body/<body onLoad=\"show_contents()\"/i;
}

$windows_index=<<EOT4;
<html>
<head>
<title>$tit</title>
<script>
<!--
var toolkit=java.awt.Toolkit.getDefaultToolkit();
var screen_size=toolkit.getScreenSize();
screen_h=screen_size.height;
function show_contents() {
   win_features="toolbar=no,scrollbars=yes,status=no,dircetories=no,resizable=yes,width=300,height="+screen_h;
   contents=window.open
      ("contents.html","contentsWin",win_features);
}
window.name="main";
//-->
</script>
</head>
$cover
</html>
EOT4

   open(F,">$out_dir/index.html");
   print F $frames_index;
   close(F);

   open(F,">$out_dir/windows_man.html");
   print F $windows_index;
   close(F);
}
