#!/bin/sh
#	Gadget to take two LaTeX files and produce a third which
#	has changebars highlighting the difference between them.
#
# Version 1.2
# Author:
#	Don Ward, Careful Computing (don@careful.co.uk)
# v1.0	April 1989
# v1.1  Feb 93	Amended to use changebar.sty (v3.0) and dvips
# v1.2  Aug 95  Added support for LaTeX209/LaTeX2e
#		Added RCS support to retrive old files

CMD=`basename $0`

SED=sed
RM="rm -f"
DIFF=diff
ED=ed
AWK=awk
GREP=grep
MV=mv
CAT=cat
MKDIR=mkdir
CO="co"

TMPFILE=/tmp/$CMD.$$
SED_CMD_FILE=$TMPFILE.sed

usage()
{
$CAT << _END_
Usage:
  $CMD [-hgG] [-d dir] old new [output]
	default output is stdout

  $CMD [-hgG] [-d dir] old
	new file on stdin, output on stdout

  $CMD [-hgG] -d dir -r rev files
	old file retrieved using RCS

  Gadget to take two LaTeX files and produce a third which
  has changebars highlighting the difference between them.
  Changebars are inserted for differences after '\begin{document}'.

  Feature: \`new' can not be named \`-'.

  Options are:
  -d dir : Write the output to file  \`dir/new', if \`new' is given or
	   to file \`dir/old'.
	   If \`dir' does not exist, it is created.
	   If \`output' is given, it is discarded.

  -r rev : If the LaTeX \`files' are kept under control of the
	   Revision Control System RCS, the old files of
	   the revision \`rev' can be retrived.
	   \`rev' is specified using the RCS conventions.
	   This option must be used together with the \`-d dir' option.
	   \`files' must be a nonempty list of files.

  -h	 : Print this info text.
  -g	 : Print some debugging info.
  -G	 : Even more debug info.

  Version 1.2: August 3. 1995
_END_
exit 1
}

# parse options and arguments
DEBUG="NO"
DIR=
REV=
# process options
while getopts d:r:gGh i $*
do
  case $i in
	d ) DIR=$OPTARG;;
	r ) REV=$OPTARG;;
	g ) DEBUG="YES" ;;
	G ) set -x; DEBUG="YES";;
	h | \
        * ) usage ;;
  esac
done

shift `expr $OPTIND - 1`

case $# in
  1 ) OLD=$1; NEW="-"; OUT=""   ;;
  2 ) OLD=$1; NEW=$2;  OUT=""   ;;
  3 ) OLD=$1; NEW=$2;  OUT="$3" ;;
  * ) usage ;;
esac

# check correct options
if [ ! -z "$DIR" ]
then
  [ -d $DIR ] || $MKDIR $DIR
fi

if [ ! -z "$REV" ]
then
  [ -z "$DIR" ] && usage
  FILES=$*
else
  FILES=$NEW
fi

# do the work
for NEW in $FILES
do
  if [ ! -z "$DIR" ]
  then
    if [ $NEW = "-" ]
    then
      OUT=$DIR/$OLD
    else
      OUT=$DIR/$NEW
    fi
  fi
  if [ ! -z "$REV" ]
  then
    OLD=${TMPFILE}.old
    $CO -p"$REV" -q $NEW > $OLD
  fi

  [ $DEBUG = "YES" ] && echo "OLD=\`$OLD' NEW=\`$NEW' OUT=\`$OUT'"

  # gather some info about the file
  # Since we have for sure only the name of the OLD file, ...
  $GREP "^\\\\begin{document}" $OLD > /dev/null
  if [ $? -eq 0 ]
  then
    [ $DEBUG = "YES" ] && echo "contains a \\begin{document}"
    HAS_BEGIN_DOC="YES"
  else
    [ $DEBUG = "YES" ] && echo "contains no \\begin{document}"
    HAS_BEGIN_DOC="NO"
  fi

  # Method to do the work:
  # 1	Use diff to get an ed script to go from file1 to file2.
  # 2	Breath on it a bit (with sed) to insert changebar commands.
  # 3	Apply modified ed script to produce (nearly) the output.
  # 4	Use awk to insert the changebars option into the \documentstyle
  #	and to handle changebar commands inside verbatim environments.
  # 5     Remove changebars before \begin{document} with sed

  #	SED commands to edit ED commands to edit old file
  $CAT > $SED_CMD_FILE <<\_END_
/^\.$/i\
\\cbend{}%
/^[0-9][0-9]*[ac]$/a\
\\cbstart{}%
/^[0-9][0-9]*,[0-9][0-9]*[ac]$/a\
\\cbstart{}%
/^[0-9][0-9]*d$/a\
i\
\\cbdelete{}%\
.
/^[0-9][0-9]*,[0-9][0-9]*d$/a\
i\
\\cbdelete{}%\
.
_END_

  # note DIFF accepts `-' as stdin
  $DIFF -b -e $OLD $NEW | \
	( $SED -f $SED_CMD_FILE ; echo w ${TMPFILE}.1 ; echo q ) | \
	$ED - $OLD

  #	AWK commands to insert Changebars style and to protect
  #	changebar commands in verbatim environments
  #       and to tell what driver is in use; we assume the `dvips' driver

  $AWK '
  BEGIN {kind=""; # we saw now \documentXXX[]{}
  }
  /^\\documentstyle/{
    kind = "209";
    if (index($0, "changebar") == 0 ) {
      opts = index($0, "[")
      if (opts > 0)
	printf "%schangebar,%s\n",substr($0,1,opts),substr($0,opts+1)
      else
	printf "\\documentstyle[changebar]%s\n", substr($0,15)
      next
    }
  }
  /^\\documentclass/{
    kind = "2e";
    printf "%s\n", $0
    printf "\\usepackage[dvips]{changebar}\n"
    next
  }
  /\\begin{document}/ {if (kind == "209" ) {print "\\driver{dvips}"}}
  /\\begin{verbatim}/{++nesting}
  /\\end{verbatim}/{--nesting}
  /\\cbstart{}%|\\cbend{}%|\cbdelete{}%/ {
    if ( nesting > 0) {
  #	changebar command in a verbatim environment: Temporarily exit,
  #	do the changebar command and reenter.
  #
  #	The obvious ( printf "\\end{verbatim}%s\\begin{verbatim} , $0 )
  #	leaves too much vertical space around the changed line(s).
  #	The following magic seeems to work
  #
	print  "\\end{verbatim}\\nointerlineskip"
	print  "\\vskip -\\ht\\strutbox\\vskip -\\ht\\strutbox"
	printf "\\vbox to 0pt{\\vskip \\ht\\strutbox%s\\vss}\n", $0
	print  "\\begin{verbatim}"
	next
	}
  }
  { print $0 }
  '  ${TMPFILE}.1 > ${TMPFILE}.2

  # if a \begin{document} is contained in the file,
  # remove the changebar commands before them

  if [ $HAS_BEGIN_DOC = "YES" ]
  then
    SED_CMD="1,/\\\\begin{document}/s/\(\\\\cb[sed][tne][adl][^{}]*{}%\)$/%%\1/"
    $SED "$SED_CMD" ${TMPFILE}.2 > ${TMPFILE}.3
  else
    $CAT ${TMPFILE}.2 > ${TMPFILE}.3
  fi
  if [ -z "$OUT" ]
  then
    $CAT ${TMPFILE}.3
 else
    $MV ${TMPFILE}.3 $OUT
  fi

done

[ $DEBUG =  "NO" ] && $RM ${TMPFILE}.*

###############################################################
