/** 
$Id$

\brief Return the pathname of NWCHEM_SRCDIR

NWChem needs to know the pathname of NWCHEM_SRCDIR in order to find the
basis libraries (both the Gaussian basis sets as well as the pseudopotentials
in the plane wave code). Originally this was handled in the Fortran77 source
code but the line length limitations generated ridiculous restrictions. 
Attempts to move this to Fortran90 failed as some compilers adopted the
silly Fortran77 line length restrictions even for free format Fortran90 code.
Hence we have to do this the nasty way in C.

The value of SRCDIR has been carried over to the string NWCHEM_SRCDIR
and the C-preprocessor inserts this value into the code. C copies the constant
value into a Fortran character string, the remainder of the Fortran string
is filled with spaces, and the resulting string is returned to NWChem. 

To use this routine the Fortran routine must include the interface block
defined in util_nwchem_srcdir.fh. This of course requires a proper Fortran
compiler that understands the ISO_C_BINDING standard. Currently this 
standard is not universally supported therefore the use of the environment
variable NWCHEM_LONG_PATHS is recommended. This variable should be set to
"Y" to indicate that these code constructs can be used with the compiler
you have.
*/
#ifdef XLFLINUX
void utilc_nwchem_srcdir_(
#else
void utilc_nwchem_srcdir(
#endif
    char * pathname, //! Pointer to the Fortran character string
    int  * length    //! The length of the Fortran character string, i.e. len(string)
)
{

    char * nwchem_srcdir = NWCHEM_SRCDIR;
    int i;
    for (i = 0; (i < *length) && (nwchem_srcdir[i] != '\0'); i++) {
        pathname[i] = nwchem_srcdir[i];
    }
    for ( ; i < *length; i++) {
        pathname[i] = ' ';
    }
}
