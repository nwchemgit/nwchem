*
* $Id$
*
*======================================================================
*
* DISCLAIMER
*
* This material was prepared as an account of work sponsored by an
* agency of the United States Government.  Neither the United States
* Government nor the United States Department of Energy, nor Battelle,
* nor any of their employees, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
* ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY,
* COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT,
* SOFTWARE, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT
* INFRINGE PRIVATELY OWNED RIGHTS.
*
* ACKNOWLEDGMENT
*
* This software and its documentation were produced with Government
* support under Contract Number DE-AC06-76RLO-1830 awarded by the United
* States Department of Energy.  The Government retains a paid-up
* non-exclusive, irrevocable worldwide license to reproduce, prepare
* derivative works, perform publicly and display publicly by or for the
* Government, including the right to distribute to other Government
* contractors.
*
*======================================================================
*
*  -- PFFT routine (version 1.0) --
*     Pacific Northwest Laboratory
*     April 5, 1995
*
*======================================================================
      SUBROUTINE PFFT3D( ISIGN, N1, N2, N3, X, MAP1, MAP2,
     $                   LWORK, WORK, LIWORK, IWORK, INFO )

* ARGUMENTS
* =========
      INTEGER           ISIGN, N1, N2, N3, LWORK, LIWORK, INFO
      INTEGER           MAP1(N3), MAP2(N2), IWORK(LIWORK)
      DOUBLE PRECISION  X(*), WORK(LWORK)

*-----------------------------------------------------------------------
*
*  DESCRIPTION
*  ===========
*
*  This routine applies a parallel 2- or 3-dimensional
*  complex-to-complex fft to X.
*
*  Assume that Z and Y are stored in FORTRAN arrays dimensioned
*
*  DOUBLE COMPLEX Y( 1:N1, 1:N2, 1:N3 )
*  DOUBLE COMPLEX Z( 1:N1, 1:N2, 1:N3 )
*
*  then this routine computes:
*
*  Y_{k1+1,k2+1,k3+1} = SUM_{j3=0}^{N3-1}
*                        [ SUM_{j2=0}^{N2-1}
*                          [ SUM_{j1=0}^{N1-1} [ Z_{j1+1,j2+1,j3+1} *
*                                                \omega_1^{j1 k1}  ] *
*                                                \omega_2^{j2 k2}  ] *
*                                                \omega_3^{j3 k3}  ] 
*  for 
*
*    k1 = 0,...,N1-1,
*    k2 = 0,...,N2-1,
*    k3 = 0,...,N3-1,
*
*  where
*
*    ISIGN = +- 1
*
*    omega_1 = exp( ISIGN * 2 * PI * i / N1 ), where i = SQRT( -1 )
*    omega_2 = exp( ISIGN * 2 * PI * i / N2 ),
*    omega_3 = exp( ISIGN * 2 * PI * i / N3 ),
*
*  Calling PFFT3D with ISIGN = +1 followed by ISIGN = -1 will scale
*  the original Z array by a factor of N1*N2*N3.  Thus, except for a
*  scale factor, ISIGN = +1 and ISIGN = -1 are inverses of each other.
*
*
*  Z represents ALL of the data on which the fft is to be applied.  The
*  array X passed to this routine is the portion of Z owned by this
*  processor.  X is assume to be distributed to processors by complete
*  planes of data in the x1-x2 plane.  In particular, processor MAP1(I)
*  is assumed to own Z( 1:N1, 1:N2, I ), I = 1 to N3.   Also, in this
*  routine X is declared DOUBLE PRECISION, rather than DOUBLE COMPLEX
*  like Z above.  See the descriptions of X and MAP1 below for more
*  details on X.
*
*
*  ARGUMENTS
*  =========
*
*   ISIGN ...... (Input) INTEGER
*                Integer determining sign of exponent in transform.
*                Allowed values are +1 and -1.
*                We refer to ISIGN = +1 as the forward transform, and
*                            ISIGN = -1 as the inverse transform.
*                 
*   N1 ......... (Input) INTEGER
*                Number of elements in first dimension of transform.
*                N1 > 0.  If N1 = 1, then this routine really does
*                a 2-dimennsional fft on N2 * N3 points.
*
*   N2 ......... (Input) INTEGER
*                Number of elements in second dimension of transform.
*                N2 > 1.
*
*   N3 ......... (Input) INTEGER
*                Number of elements in third dimension of transform.
*                N3 > 1.
*
*   X .......... (Input/Output) DOUBLE PRECISION, dim( 2*N1*N2*NVECS1 ),
*                            or DOUBLE COMPLEX,   dim(   N1*N2*NVECS1 )
*
*                where NVECS1 is the number of entries in MAP1(1:N3)
*                which are equal to this processor's id, ME.
*
*                On Entry
*                ========
*                Let Z be the DOUBLE COMPLEX Z( 1:N1, 1:N2, 1:N3 )
*                Fortran array described above and MAP1 be as described
*                below.  Then X is a DOUBLE PRECISION array containing
*                the part of Z owned by this processor (as determine by
*                MAP1) stored in packed form.  In particular, let I_m
*                be the M-th value of I such that MAP1(I) = ME, then 
*
*                X( 1+(M-1)*2*N1*N2: M*2*N1*N2 ) = Z( 1:N1, 1:N2, I_m )
*
*                stored according to:
*
*                K = 1
*                DO 30 I3 = 1, N3
*                  IF( MAP1( I3 ) .EQ. ME ) THEN
*
*                     DO 20 I2 = 1, N2
*                      DO 10 I1 = 1, N1
*                         X( K   ) = REAL( Z( I1, I2, I_m ) )
*                         X( K+1 ) = IMAG( Z( I1, I2, I_m ) )
*                         K = K + 2
*            10        CONTINUE
*            20      CONTINUE
*
*                  ENDIF
*            30  CONTINUE
*
*                If X is DOUBLE PRECISION in the calling routine, then
*                the code above is appropriate.
*
*                If X is declared DOUBLE COMPLEX in the calling routine,
*                then X will be corrent in this routine if
*   
*                X( (M-1)*N1*N2 + (I2-1)*N1 + I1 ) = Z( I1, I2, I_m ).
*
*                On Exit
*                =======
*                Let Y be the array define above, i.e., a
*                DOUBLE COMPLEX array containing the transformed Z.
*                Then, on exit X contains the part of Y owned by this
*                processor, which is the same as the part of Z owned
*                by this processor.  The storage format on exit is
*                the same as it was on entry, except now X contains
*                the transformed data Y, rather than the input data Z.
*                 
*   MAP1 ....... (Input) INTEGER array, dim( N3 )
*                MAP1(I) is the id of the processor which owns the
*                I-th x1-x2 plane of data to be transformed, i.e., owns
*                Z( 1:N1, 1:N2, I) 
*
*   MAP2 ....... (Input) INTEGER array, dim( N2 )
*                As part of the fft an inter-processor transpose is done
*                in the x2-x3 plane.  I.e., after this transpose each
*                processor owns contiguous sets of data in the x1-x3
*                plane, rather data in the x1-x2 plane.  MAP2 determines
*                how this x2-x3 data should be distributed to processors
*                following the transpose.  In particular, 
*                MAP2(I) is the id of the processor which should own the
*                I-th X1-X3 plane of X data, i.e., owns
*                Z( 1:N1, I, 1:N3).
*
*                The only thins that is really important about MAP2 is
*                that this transposed data is stored in the work array
*                WORK, and MAP2 is used to specify how much of X
*                transpose that each processor has room to store.
*
*   LWORK ...... (Input) INTEGER
*                The size of the DOUBLE PRECISION work array WORK.
*                LWORK >= ?
*
*   WORK ....... (Workspace) DOUBLE PRECISION array, dim( >= LWORK )
*                DOUBLE PRECISION workspace
*
*   LIWORK ..... (Input) INTEGER
*                The size of the INTEGER work array WORK.
*                LIWORK >= ?
*
*   IWORK ...... (Workspace) INTEGER array, dim( >= LIWORK )
*                INTEGER workspace
*
*   INFO ....... (Output) INTEGER
*                Error flag.
*
*                A negative INFO on a processors causes that processor
*                to STOP program execution via a call to XSTOP.
*
*                PFFT3D attempts to return consistent values of INFO
*                on all processors in MAP1 and/or MAP2, but this is not
*                always possible.  In particular, if a processor finds
*                that N2 or N3 is < 1, or MAP1 and/or MAP2 contains
*                improper data, or WORK or IWORK is not big enough,
*                or the processor's id is not in MAP1 or MAP2,
*                then that processor just exits with a negative INFO,
*                but without communicating with any other processors.
*
*                = 0, then computation completed successfully,
*
*                = -1 to -12, then the INFO-th input argument had
*                             an illegal value on this processor.
*
*                = -50, then this processor is not in MAP1 or MAP2,
*                       hence should not have called this routine.
*
*                = -51, then an input argument to one of the processors
*                       had an illegal value (that processor returns
*                       INFO in the range -1 to -12), but all input
*                       arguments were ok on this processor
*                       
*                = -52, then ISIGN, N1, N2, N3, MAP1 and/or MAP2
*                       are not the same on all processors in MAP1
*                       and/or MAP2.
*
* ----------------------------------------------------------------------

      EXTERNAL          PFFT1I, PFFT1, MXINIT
	character*6 junk

************************************************************************

*     Initialize MX communication library
      CALL MXINIT()

*     Check input data.
c
	junk = 'pfft3d'
      CALL CHK3D( junk, 1, ISIGN, N1, N2, N3, X, MAP1, MAP2,
     $            LWORK, WORK, LIWORK, IWORK, INFO )
c
      IF( INFO .NE. 0 )
     $  RETURN

************************************************************************

      CALL PFFT( PFFT1, PFFT1I, PFFT1, PFFT1I, 3, 1,
     $           ISIGN, N1, N2, N3, X, MAP1, MAP2,
     $           LWORK, WORK, LIWORK, IWORK, INFO )

      RETURN
      END
