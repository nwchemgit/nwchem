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

s/^\(     .\)[Dd][Aa][Xx][Pp][Yy]\([^A-Za-z0-9_]\)/\1saxpy\2/
s/^\(     .\)[Dd][Dd][Oo][Tt]\([^A-Za-z0-9_]\)/\1sdot\2/
s/^\(     .\)[Dd][Ss][Cc][Aa][Ll]\([^A-Za-z0-9_]\)/\1sscal\2/
s/^\(     .\)[Dd][Cc][Oo][Pp][Yy]\([^A-Za-z0-9_]\)/\1scopy\2/
s/^\(     .\)[Dd][Nn][Rr][Mm]2\([^A-Za-z0-9_]\)/\1snrm2\2/
s/^\(     .\)[Ii][Dd][Aa][Mm][Aa][Xx]\([^A-Za-z0-9_]\)/\1isamax\2/
s/^\(     .\)[Dd][Aa][Ss][Uu][Mm]\([^A-Za-z0-9_]\)/\1sasum\2/
s/^\(     .\)[Dd][Rr][Oo][Tt]\([^A-Za-z0-9_]\)/\1srot\2/
s/^\(     .\)[Dd][Rr][Oo][Tt][Gg]\([^A-Za-z0-9_]\)/\1srotg\2/
s/^\(     .\)[Dd][Ss][Ww][Aa][Pp]\([^A-Za-z0-9_]\)/\1sswap\2/
s/^\(     .\)[Dd][Mm][Aa][Cc][Hh]\([^A-Za-z0-9_]\)/\1smach\2/
s/^\(     .\)[Dd][Gg][Ee][Mm][Mm]\([^A-Za-z0-9_]\)/\1sgemm\2/
s/^\(     .\)[Dd][Gg][Ee][Mm][Vv]\([^A-Za-z0-9_]\)/\1sgemv\2/
s/^\(     .\)[Dd][Gg][Ee][Ss][Vv]\([^A-Za-z0-9_]\)/\1sgesv\2/
s/^\(     .\)[Dd][Ss][Pp][Mm][Vv]\([^A-Za-z0-9_]\)/\1sspmv\2/
s/^\(     .\)[Dd][Gg][Ee][Rr]\([^A-Za-z0-9_]\)/\1sger\2/
s/^\(     .\)[Dd][Gg][Ee][Ss][Vv][Dd]\([^A-Za-z0-9_]\)/\1sgesvd\2/
s/^\(     .\)[Dd][Ss][Pp][Ss][Vv][Xx]\([^A-Za-z0-9_]\)/\1sspsvx\2/
s/^\(     .\)[Dd][Gg][Ee][Tt][Rr][Ss]\([^A-Za-z0-9_]\)/\1sgetrs\2/
s/^\(     .\)[Dd][Gg][Ee][Tt][Rr][Ff]\([^A-Za-z0-9_]\)/\1sgetrf\2/
s/^\(     .\)[Dd][Pp][Oo][Tt][Rr][Ff]\([^A-Za-z0-9_]\)/\1spotrf\2/
s/^\(     .\)[Dd][Pp][Oo][Tt][Rr][Ii]\([^A-Za-z0-9_]\)/\1spotri\2/
s/^\(     .\)[Dd][Pp][Oo][Ss][Vv]\([^A-Za-z0-9_]\)/\1sposv\2/
s/^\(     .\)[Dd][Ll][Aa][Ss][Ee][Tt]\([^A-Za-z0-9_]\)/\1slaset\2/
s/^\(     .\)[Dd][Ll][Aa][Mm][Cc][Hh]\([^A-Za-z0-9_]\)/\1slamch\2/
s/^\(     .\)[Dd][Ss][Yy][Ee][Vv]\([^A-Za-z0-9_]\)/\1ssyev\2/
s/^\(     .\)[Dd][Ss][Yy][Gg][Vv]\([^A-Za-z0-9_]\)/\1ssygv\2/
s/^\(     .\)[Ii][Zz][Aa][Mm][Aa][Xx]\([^A-Za-z0-9_]\)/\1icamax\2/

s/^\(     .\)[Dd][Aa][Xx][Pp][Yy]$/\1saxpy/
s/^\(     .\)[Dd][Dd][Oo][Tt]$/\1sdot/
s/^\(     .\)[Dd][Ss][Cc][Aa][Ll]$/\1sscal/
s/^\(     .\)[Dd][Cc][Oo][Pp][Yy]$/\1scopy/
s/^\(     .\)[Dd][Nn][Rr][Mm]2$/\1snrm2/
s/^\(     .\)[Ii][Dd][Aa][Mm][Aa][Xx]$/\1isamax/
s/^\(     .\)[Dd][Aa][Ss][Uu][Mm]$/\1sasum/
s/^\(     .\)[Dd][Rr][Oo][Tt]$/\1srot/
s/^\(     .\)[Dd][Rr][Oo][Tt][Gg]$/\1srotg/
s/^\(     .\)[Dd][Ss][Ww][Aa][Pp]$/\1sswap/
s/^\(     .\)[Dd][Mm][Aa][Cc][Hh]$/\1smach/
s/^\(     .\)[Dd][Gg][Ee][Mm][Mm]$/\1sgemm/
s/^\(     .\)[Dd][Gg][Ee][Mm][Vv]$/\1sgemv/
s/^\(     .\)[Dd][Gg][Ee][Ss][Vv]$/\1sgesv/
s/^\(     .\)[Dd][Ss][Pp][Mm][Vv]$/\1sspmv/
s/^\(     .\)[Dd][Gg][Ee][Rr]$/\1sger/
s/^\(     .\)[Dd][Gg][Ee][Ss][Vv][Dd]$/\1sgesvd/
s/^\(     .\)[Dd][Ss][Pp][Ss][Vv][Xx]$/\1sspsvx/
s/^\(     .\)[Dd][Gg][Ee][Tt][Rr][Ss]$/\1sgetrs/
s/^\(     .\)[Dd][Gg][Ee][Tt][Rr][Ff]$/\1sgetrf/
s/^\(     .\)[Dd][Pp][Oo][Tt][Rr][Ff]$/\1spotrf/
s/^\(     .\)[Dd][Pp][Oo][Tt][Rr][Ii]$/\1spotri/
s/^\(     .\)[Dd][Ll][Aa][Ss][Ee][Tt]$/\1slaset/
s/^\(     .\)[Dd][Ll][Aa][Mm][Cc][Hh]$/\1slamch/
s/^\(     .\)[Dd][Ss][Yy][Ee][Vv]$/\1ssyev/
s/^\(     .\)[Dd][Ss][Yy][Gg][Vv]$/\1ssygv/
s/^\(     .\)[Ii][Zz][Aa][Mm][Aa][Xx]$/\1icamax/
s/^\(     .\)[Dd][Gg][Ee][Bb][Aa][Kk]$/\1sgebak/
s/^\(     .\)[Dd][Gg][Ee][Bb][Aa][Ll]$/\1sgebal/
s/^\(     .\)[Dd][Gg][Ee][Hh][Rr][Dd]$/\1sgehrd/
s/^\(     .\)[Dd][Hh][SS][Ee][Qq][Rr]$/\1shseqr/
s/^\(     .\)[Dd][Ll][Aa][Bb][Aa][Dd]$/\1slabad/
s/^\(     .\)[Dd][Oo][Rr][Gg][Hh][Rr]$/\1sorghr/
s/^\(     .\)[Dd][Tt][Rr][Ee][Vv][Cc]$/\1strevc/

