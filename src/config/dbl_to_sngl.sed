#
# $Id$
#
# This script performs blas keyword substitutions using the sed editor.
# It is divided into four separate sed commands because of buffer overflow 
# on some machines, including Cray unicos.  
# (1) The first command substitutes regular embedded keywords in noncomment
#     lines.
# (2) The second version substitutes keywords that occur at the end of 
#     noncomment lines.
# (3) The third version substitutes embedded keywords that are jammed next
#     to continuation characters.
# (4) The fourth version substitutes jammed keywords at the end of
#     continuation lines.
# This is not bulletproof, but it catches almost all keyword occurences.  It
# is recommended that your coding practices be consistent with this script.
#
# 08-feb-90 four-part division. -rls
# 14-dec-88 written by ron shepard. based on a previous script written
#           by eric stahlberg.

/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Aa][Xx][Pp][Yy]\([^A-Za-z0-9_]\)/\1saxpy\2/g
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Dd][Oo][Tt]\([^A-Za-z0-9_]\)/\1sdot\2/g
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Ss][Cc][Aa][Ll]\([^A-Za-z0-9_]\)/\1sscal\2/g
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Cc][Oo][Pp][Yy]\([^A-Za-z0-9_]\)/\1scopy\2/g
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Nn][Rr][Mm]2\([^A-Za-z0-9_]\)/\1snrm2\2/g
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Ii][Dd][Aa][Mm][Aa][Xx]\([^A-Za-z0-9_]\)/\1isamax\2/g
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Aa][Ss][Uu][Mm]\([^A-Za-z0-9_]\)/\1sasum\2/g
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Rr][Oo][Tt]\([^A-Za-z0-9_]\)/\1srot\2/g
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Rr][Oo][Tt][Gg]\([^A-Za-z0-9_]\)/\1srotg\2/g
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Ss][Ww][Aa][Pp]\([^A-Za-z0-9_]\)/\1sswap\2/g
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Mm][Aa][Cc][Hh]\([^A-Za-z0-9_]\)/\1smach\2/g
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Gg][Ee][Mm][Mm]\([^A-Za-z0-9_]\)/\1sgemm\2/g
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Gg][Ee][Mm][Vv]\([^A-Za-z0-9_]\)/\1sgemv\2/g
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Gg][Ee][Ss][Vv]\([^A-Za-z0-9_]\)/\1sgesv\2/g
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Ss][Pp][Mm][Vv]\([^A-Za-z0-9_]\)/\1sspmv\2/g
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Gg][Ee][Rr]\([^A-Za-z0-9_]\)/\1sger\2/g
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Gg][Ee][Ss][Vv][Dd]\([^A-Za-z0-9_]\)/\1sgesvd\2/g
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Ss][Pp][Ss][Vv][Xx]\([^A-Za-z0-9_]\)/\1sspsvx\2/g
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Gg][Ee][Tt][Rr][Ss]\([^A-Za-z0-9_]\)/\1sgetrs\2/g
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Gg][Ee][Tt][Rr][Ff]\([^A-Za-z0-9_]\)/\1sgetrf\2/g
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Pp][Oo][Tt][Rr][Ff]\([^A-Za-z0-9_]\)/\1spotrf\2/g
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Pp][Oo][Tt][Rr][Ii]\([^A-Za-z0-9_]\)/\1spotri\2/g
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Pp][Oo][Ss][Vv]\([^A-Za-z0-9_]\)/\1sposv\2/g
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Ll][Aa][Ss][Ee][Tt]\([^A-Za-z0-9_]\)/\1slaset\2/g
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Ll][Aa][Mm][Cc][Hh]\([^A-Za-z0-9_]\)/\1slamch\2/g
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Ss][Yy][Ee][Vv]\([^A-Za-z0-9_]\)/\1ssyev\2/g
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Ss][Yy][Gg][Vv]\([^A-Za-z0-9_]\)/\1ssygv\2/g
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Ii][Zz][Aa][Mm][Aa][Xx]\([^A-Za-z0-9_]\)/\1icamax\2/g

/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Aa][Xx][Pp][Yy]$/\1saxpy/
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Dd][Oo][Tt]$/\1sdot/
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Ss][Cc][Aa][Ll]$/\1sscal/
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Cc][Oo][Pp][Yy]$/\1scopy/
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Nn][Rr][Mm]2$/\1snrm2/
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Ii][Dd][Aa][Mm][Aa][Xx]$/\1isamax/
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Aa][Ss][Uu][Mm]$/\1sasum/
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Rr][Oo][Tt]$/\1srot/
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Rr][Oo][Tt][Gg]$/\1srotg/
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Ss][Ww][Aa][Pp]$/\1sswap/
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Mm][Aa][Cc][Hh]$/\1smach/
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Gg][Ee][Mm][Mm]$/\1sgemm/
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Gg][Ee][Mm][Vv]$/\1sgemv/
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Gg][Ee][Ss][Vv]$/\1sgesv/
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Ss][Pp][Mm][Vv]$/\1sspmv/
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Gg][Ee][Rr]$/\1sger/
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Gg][Ee][Ss][Vv][Dd]$/\1sgesvd/
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Ss][Pp][Ss][Vv][Xx]$/\1sspsvx/
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Gg][Ee][Tt][Rr][Ss]$/\1sgetrs/
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Gg][Ee][Tt][Rr][Ff]$/\1sgetrf/
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Pp][Oo][Tt][Rr][Ff]$/\1spotrf/
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Pp][Oo][Tt][Rr][Ii]$/\1spotri/
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Ll][Aa][Ss][Ee][Tt]$/\1slaset/
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Ll][Aa][Mm][Cc][Ht]$/\1slamch/
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Ss][Yy][Ee][Vv]$/\1ssyev/
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Ss][Yy][Gg][Vv]$/\1ssygv/
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Ii][Zz][Aa][Mm][Aa][Xx]$/\1icamax/
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Gg][Ee][Bb][Aa][Kk]$/\1sgebak/
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Gg][Ee][Bb][Aa][Ll]$/\1sgebal/
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Gg][Ee][Hh][Rr][Dd]$/\1sgehrd/
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Hh][SS][Ee][Qq][Rr]$/\1shseqr/
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Ll][Aa][Bb][Aa][Dd]$/\1slabad/
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Oo][Rr][Gg][Hh][Rr]$/\1sorghr/
/^[ 0-9]/ s/\([^A-Za-z0-9_]\)[Dd][Tt][Rr][Ee][Vv][Cc]$/\1strevc/
