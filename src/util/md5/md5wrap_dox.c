/* This file contains only documentation for Doxygen.
   It does not contain any code and therefore should not be compiled.
*/

/**
\file md5wrap.c
\defgroup checksum Check sums
\ingroup checksum
@{
\brief Interface to check sum capabilities

Checksums are useful for rapid comparison and validation of data,
such as digital signatures for verification of important messages, or,
more relevant to us, to determine if input and disk resident restart
data are still consistent.  The checksum routines provided here are
wrappers around the RSA implementation of the RSA Data Security, Inc.\
MD5 Message-Digest Algorithm.  It is the reference implementation for
internet RFC 1321, The MD5 Message-Digest Algorithm, and as such has
been extensively tested and there are no restrictions placed upon its
distribution or export.  License is granted by RSA to make and use
derivative works provided that such works are identified as "derived
from the RSA Data Security, Inc. MD5 Message-Digest Algorithm" in all
material mentioning or referencing the derived work.  Consider this
done.  The unmodified network posting is included in `md5.txt` for
reference.

> MD5 is probably the strongest checksum algorithm most people will need
> for common use.  It is conjectured that the difficulty of coming up
> with two messages having the same message digest is on the order of
> \f$2^{64}\f$ operations, and that the difficulty of coming up with any
> message having a given message digest is on the order of \f$2^{128}\f$
> operations.

The checksums are returned (through the NWChem interface) as character
strings containing a 32 character hexadecimal representation of the
128 bit binary checksum.  This form loses no information, may be
readily compared with single statements of standard C/F77, is easily
printed, and does not suffer from byte ordering problems.  The
checksum depends on both the value and order of data, and thus
differing numerical representations, floating-point rounding
behaviour, and byte ordering, make the checksum of all but simple text
data usually machine dependent unless great care is taken when moving
data between machines.  The Fortran test program merely tests the
Fortran interface.  For a more definitive test of MD5 make
`mddriver` and execute it with the `-x` option, comparing
output with that in `md5.txt`.

C routines should include `checksum.h` for prototypes.
There is no Fortran header file since there are no functions.

The checksum of a contiguous block of data may be generated with 
~~~~
call checksum_simple(len, data, sum)
~~~~
to get more sophisticated see below and have a look at `ftest.F`.

@}
*/
/* $Id$ */
