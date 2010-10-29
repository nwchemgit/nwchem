#
# $Id$
#
s/^\(     .\)[Ss][Aa][Xx][Pp][Yy]\([^A-Za-z0-9_]\)/\1daxpy\2/
s/^\(     .\)[Ss][Dd][Oo][Tt]\([^A-Za-z0-9_]\)/\1ddot\2/
s/^\(     .\)[Ss][Ss][Cc][Aa][Ll]\([^A-Za-z0-9_]\)/\1dscal\2/
s/^\(     .\)[Ss][Cc][Oo][Pp][Yy]\([^A-Za-z0-9_]\)/\1dcopy\2/
s/^\(     .\)[Ss][Nn][Rr][Mm]2\([^A-Za-z0-9_]\)/\1dnrm2\2/
s/^\(     .\)[Ii][Ss][Aa][Mm][Aa][Xx]\([^A-Za-z0-9_]\)/\1idamax\2/
s/^\(     .\)[Ss][Aa][Ss][Uu][Mm]\([^A-Za-z0-9_]\)/\1dasum\2/
s/^\(     .\)[Ss][Rr][Oo][Tt]\([^A-Za-z0-9_]\)/\1srot\2/
s/^\(     .\)[Ss][Rr][Oo][Tt][Gg]\([^A-Za-z0-9_]\)/\1drotg\2/
s/^\(     .\)[Ss][Ss][Ww][Aa][Pp]\([^A-Za-z0-9_]\)/\1dswap\2/
s/^\(     .\)[Ss][Mm][Aa][Cc][Hh]\([^A-Za-z0-9_]\)/\1dmach\2/
s/^\(     .\)[Ss][Gg][Ee][Mm][Mm]\([^A-Za-z0-9_]\)/\1dgemm\2/
s/^\(     .\)[Ss][Gg][Ee][Mm][Vv]\([^A-Za-z0-9_]\)/\1dgemv\2/
s/^\(     .\)[Ss][Gg][Ee][Ss][Vv]\([^A-Za-z0-9_]\)/\1dgesv\2/
s/^\(     .\)[Ss][Ss][Pp][Mm][Vv]\([^A-Za-z0-9_]\)/\1dspmv\2/
s/^\(     .\)[Ss][Gg][Ee][Rr]\([^A-Za-z0-9_]\)/\1dger\2/
s/^\(     .\)[Ss][Gg][Ee][Ss][Vv][Dd]\([^A-Za-z0-9_]\)/\1dgesvd\2/
s/^\(     .\)[Ss][Ss][Pp][Ss][Vv][Xx]\([^A-Za-z0-9_]\)/\1dspsvx\2/
s/^\(     .\)[Ss][Gg][Ee][Tt][Rr][Ss]\([^A-Za-z0-9_]\)/\1dgetrs\2/
s/^\(     .\)[Ss][Pp][Oo][Tt][Rr][fF]\([^A-Za-z0-9_]\)/\1dpotrf\2/
s/^\(     .\)[Ss][Pp][Oo][Tt][Rr][Ii]\([^A-Za-z0-9_]\)/\1dpotri\2/
s/^\(     .\)[Ss][Pp][Oo][Ss][Vv]\([^A-Za-z0-9_]\)/\1dposv\2/
s/^\(     .\)[Ss][Ll][Aa][Ss][Ee][Tt]\([^A-Za-z0-9_]\)/\1dlaset\2/
s/^\(     .\)[Ss][Ll][Aa][Mm][Cc][Hh]\([^A-Za-z0-9_]\)/\1dlamch\2/
s/^\(     .\)[Ss][Ss][Yy][Ee][Vv]\([^A-Za-z0-9_]\)/\1dsyev\2/
s/^\(     .\)[Ss][Ss][Yy][Gg][Vv]\([^A-Za-z0-9_]\)/\1dsygv\2/
s/^\(     .\)[Ii][Cc][Aa][Mm][Aa][Xx]\([^A-Za-z0-9_]\)//\1izamax\2/

s/^\(     .\)[Ss][Aa][Xx][Pp][Yy]$/\1daxpy/
s/^\(     .\)[Ss][Dd][Oo][Tt]$/\1ddot/
s/^\(     .\)[Ss][Ss][Cc][Aa][Ll]$/\1dscal/
s/^\(     .\)[Ss][Cc][Oo][Pp][Yy]$/\1dcopy/
s/^\(     .\)[Ss][Nn][Rr][Mm]2$/\1dnrm2/
s/^\(     .\)[Ii][Ss][Aa][Mm][Aa][Xx]$/\1idamax/
s/^\(     .\)[Ss][Aa][Ss][Uu][Mm]$/\1dasum/
s/^\(     .\)[Ss][Rr][Oo][Tt]$/\1drot/
s/^\(     .\)[Ss][Rr][Oo][Tt][Gg]$/\1drotg/
s/^\(     .\)[Ss][Ss][Ww][Aa][Pp]$/\1dswap/
s/^\(     .\)[Ss][Mm][Aa][Cc][Hh]$/\1dmach/
s/^\(     .\)[Ss][Gg][Ee][Mm][Mm]$/\1dgemm/
s/^\(     .\)[Ss][Gg][Ee][Mm][Vv]$/\1dgemv/
s/^\(     .\)[Ss][Gg][Ee][Ss][Vv]$/\1dgesv/
s/^\(     .\)[Ss][Ss][Pp][Mm][Vv]$/\1dspmv/
s/^\(     .\)[Ss][Gg][Ee][Rr]$/\1dger/
s/^\(     .\)[Ss][Gg][Ee][Ss][Vv][Dd]$/\1dgesvd/
s/^\(     .\)[Ss][Ss][Pp][Ss][Vv][Xx]$/\1dspsvx/
s/^\(     .\)[Ss][Gg][Ee][Tt][Rr][Ss]$/\1dgetrs/
s/^\(     .\)[Ss][Pp][Oo][Tt][Rr][Ff]$/\1dpotrf/
s/^\(     .\)[Ss][Pp][Oo][Tt][Rr][Ii]$/\1dpotri/
s/^\(     .\)[Ss][Ll][Aa][Ss][Ee][Tt]$/\1dlaset/
s/^\(     .\)[Ss][Ll][Aa][Mm][Cc][Hh]$/\1dlamch/
s/^\(     .\)[Ss][Ss][Yy][Ee][Vv]$/\1dsyev/
s/^\(     .\)[Ss][Ss][Yy][Gg][Vv]$/\1dsygv/
s/^\(     .\)[Ii][Cc][Aa][Mm][Aa][Xx]$/\1izamax/
s/^\(     .\)[Ss][Gg][Ee][Bb][Aa][Kk]$/\1dgebak/
s/^\(     .\)[Ss][Gg][Ee][Bb][Aa][Ll]$/\1dgebal/
s/^\(     .\)[Ss][Gg][Ee][Hh][Rr][Dd]$/\1dgehrd/
s/^\(     .\)[Ss][Hh][SS][Ee][Qq][Rr]$/\1dhseqr/
s/^\(     .\)[Ss][Ll][Aa][Bb][Aa][Dd]$/\1dlabad/
s/^\(     .\)[Ss][Oo][Rr][Gg][Hh][Rr]$/\1dorghr/
s/^\(     .\)[Ss][Tt][Rr][Ee][Vv][Cc]$/\1dtrevc/





