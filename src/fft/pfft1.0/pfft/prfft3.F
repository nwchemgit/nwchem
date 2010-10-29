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
      SUBROUTINE PRFFT3( IDIR, ISIGN, N1, N2, N3, X, MAP1, MAP2,
     $                   LWORK, WORK, LIWORK, IWORK, INFO )

* ARGUMENTS
* =========
      INTEGER           IDIR, ISIGN, N1, N2, N3, LWORK, LIWORK, INFO
      INTEGER           MAP1(N3), MAP2(N2), IWORK(LIWORK)
      DOUBLE PRECISION  X(*), WORK(LWORK)

*-----------------------------------------------------------------------
*
* DESCRIPTION
* ===========
*
*  This routine applies a parallel 2- or 3-dimensional
*  real-to-complex (IDIR = +1) or complex-to-real (IDIR = -1 ) 
*  discrete fourier transform (DFT) to X.
*
*  Assume that Z and Y are three dimensional arrays of data:
*
*  Y_{ 1:N1, 1:N2, 1:N3 }, Z_{ 1:N1, 1:N2, 1:N3 }
*
*  If IDIR = +1 (real-to-complex DFT of Z), then this routine computes:
*
*  Y_{ k1+1,k2+1,k3+1} = SUM_{j3=0}^{N3-1}
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
*  If IDIR = -1 (complex-to-real DFT of Z), then the same computation
*   is done but with the summations done in reversed order.
*
*  Calling PRFFT3 with ISIGN = +1 or -1, IDIR = +1 followed by
*  ISIGN = -ISIGN, IDIR = -1 will scale the original Z array by a factor
*  of N1*N2*N3.  Thus, except for a scale factor the complex-to-real
*  and real-to-complex DFTs are inverses of each other.
*
*
*  For a real-to-complex DFT, IDIR = +1, Z is real and Y is complex,
*  hence in principle Y requires twice as much storage as Z.  In
*  practice, however, only 
*
*  Y_{ 1:K1, 1:K2, 1:N3},    where  K1 = N1/2+1, K2 = N2, when N1 > 1,
*                            or     K1 = 1, K2 = N2/2+1  when N1 = 1,
*
*  needs to be computed and stored.  The rest of Y,
*  Y_{ M+1:N1, 1:N2, 1:N3), can be derive, if needed, from symmetries
*  in Y.  Thus, Y only requires only slight more room for storage
*  than Z.
*
*  For a complex-to-real DFT, IDIR = -1, Y is real and Z is complex.
*  This DFT can only be applied to Z if it has the special symmetry
*  associated with the result of a real-to-complex DFT.  Furthermore,
*  this symmetry in Z means that only
*
*  Z_{ 1:K1, 1:K2, 1:N3},    where  K1 = N1/2+1, K2 = N2, when N1 > 1,
*                            or     K1 = 1, K2 = N2/2+1  when N1 = 1,
*
*  is need or used when computing the complex-to-real DFT of Z.
*
*  Z represents ALL of the data on which the DFT is to be applied.  The
*  array X passed to this routine is the portion of Z owned by this
*  processor.  X is assume to be distributed to processors by complete
*  planes of data in the x1-x2 plane.  In particular, processor MAP1(I)
*  is assumed to own Z_{ 1:N1, 1:N2, I ), I = 1 to N3.   Also, in this
*  routine X is always declared DOUBLE PRECISION.
*  See the descriptions of X and MAP1 below for more
*  details on X.
*
*
*  ARGUMENTS
*  =========
*   IDIR ....... (Input) INTEGER
*                = +1 or -1.
*                IDIR = +1, then do real-to-complex DFT.
*                IDIR = -1, then do complex-to-real DFT.
*                 
*   ISIGN ...... (Input) INTEGER
*                Integer determining sign of exponent in transform.
*                Allowed values are +1 and -1.
*                 
*   N1 ......... (Input) INTEGER
*                Number of elements in first dimension of transform.
*                N1 > 0.  If N1 = 1, then this routine really does
*                a 2-dimennsional DFT on N2 * N3 points.
*
*   N2 ......... (Input) INTEGER
*                Number of elements in second dimension of transform.
*                N2 > 1.
*
*   N3 ......... (Input) INTEGER
*                Number of elements in third dimension of transform.
*                N3 > 1.
*
*   X .......... (Input/Output) DOUBLE PRECISION, dim( 2*K1*K2*NVECS1 ),
*                            or DOUBLE COMPLEX,   dim(   K1*K2*NVECS1 )
*
*                where NVECS1 is the number of entries in MAP1(1:N3)
*                             which are equal to this processor's id,
*                             ME,
*
*                      K1 ... = 1        if N1 = 1
*                             = N1/2 + 1 if N1 > 1
*
*                      K2 ... = N2       if N1 > 1
*                             = N2/2 + 1 if N1 = 1
*
*                IF IDIR = +1
*                ============
*                  On Entry
*                  ========
*
*                  Let Z_{ 1:N1, 1:N2, 1:N3 } be the array of real
*                  data to which the real-to-complex DFT is to be
*                  applied.  Let processor MAP1(I) own the
*                  Ith N1-N2 plane of Z data, Z_{ 1:N1, 1:N2, I}.
*                  Let ME be this processor's id.
*
*                  Then X is a DOUBLE PRECISION array containing
*                  the part of Z owned by this processor stored in
*                  a modified packed form.  In particular, let I_m
*                  be the M-th value of I such that MAP1(I) = ME, then 
*
*                  X( 1+(M-1)*2*K1*K2: M*2*K1*K2 ) = Z_{ 1:N1, 1:N2, I_m}
*
*                  stored according to:
*
*                  DOUBLE PRECISION X(*)
*
*                  K = 1+(M-1)*2*K1*K2
*                  DO 20 I2 = 1, N2
*                    DO 10 I1 = 1, N1
*                       X( K   ) = Z_{ I1, I2, I_m }
*                       K = K + 1
*            10      CONTINUE
*                    IF( N1 .GT. 1 ) K = K + (2*K1 - N1 )
*            20    CONTINUE
*                  IF( N1 .EQ. 1 ) K = K + (2*K2 - N2 )
*
*                  This is equivalent to:
*
*                  DOUBLE PRECISION X(2*K1, K2, N3)
*
*                  DO 20 I2 = 1, N2
*                    DO 10 I1 = 1, N1
*                       X( I1, I2, M ) = Z_{ I1, I2, I_m }
*            10      CONTINUE
*            20    CONTINUE
*
*
*                  Note that X has extra storage space allocated to
*                  allow the complex vector Y (the DFT of Z) to
*                  overwrite X.
*
*
*                  On Exit
*                  =======
*                  Let Y be the array define above, i.e., an array of
*                  complex data containing the real-to-complex
*                  transform of  Z.
*                  Then on exit  X is a DOUBLE PRECISION array
*                  containing
*                  the part of Z owned by this processor stored in
*                  a packed form.  In particular, let I_m
*                  be the M-th value of I such that MAP1(I) = ME, then 
*
*                  X( 1+(M-1)*2*K1*K2: M*2*K1*K2 ) = Y_{ 1:K1, 1:K2, I_m }
*
*                  stored according to:
*
*                  DOUBLE PRECISION X(*)
*
*                  K = 1+(M-1)*2*K1*K2
*                  DO 20 I2 = 1, N2
*                    DO 10 I1 = 1, N1
*                       X( K   ) = REAL( Y_{ I1, I2, I_m } )
*                       X( K+1 ) = IMAG( Y_{ I1, I2, I_m } )
*                       K = K + 2
*            10      CONTINUE
*            20    CONTINUE
*
*                  This is equivalent to:
*
*                  DOUBLE PRECISION X(2*K1, K2, N3)
*
*                  DO 20 I2 = 1, N2
*                    DO 10 I1 = 1, N1, 2
*                       X( I1,   I2, M ) = REAL( Y_{ I1, I2, I_m } )
*                       X( I1+1, I2, M ) = IMAG( Y_{ I1, I2, I_m } )
*            10      CONTINUE
*            20    CONTINUE
*
*                  and is also equivalent to:
*
*                  DOUBLE COMPLEX X(K1, K2, N3)
*
*                  DO 20 I2 = 1, N2
*                    DO 10 I1 = 1, N1
*                       X( I1, I2, M ) = Y_{ I1, I2, I_m }
*            10      CONTINUE
*            20    CONTINUE
*
*
*                IF IDIR = -1
*                ============
*
*                  On Entry
*                  ========
*                  X is as describe above for IDIR = +1 "On Exit"
*                 
*                  On Exit
*                  =======
*                  X is as describe above for IDIR = +1 "On Entry"
*                 
*   MAP1 ....... (Input) INTEGER array, dim( N3 )
*                MAP1(I) is the id of the processor which owns the
*                I-th x1-x2 plane of data to be transformed, i.e., owns
*                Z_{ 1:N1, 1:N2, I) 
*
*   MAP2 ....... (Input) INTEGER array, dim( N2 )
*                As part of the DFT an inter-processor transpose is done
*                in the x2-x3 plane.  I.e., after this transpose each
*                processor owns contiguous sets of data in the x1-x3
*                plane, rather data in the x1-x2 plane.  MAP2 determines
*                how this x2-x3 data should be distributed to processors
*                following the transpose.  In particular, 
*                MAP2(I) is the id of the processor which should own the
*                I-th X1-X3 plane of X data, i.e., owns
*                Z_{ 1:N1, I, 1:N3).
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
*                = -1 to -13, then the INFO-th input argument had
*                             an illegal value on this processor.
*
*                = -50, then this processor's id is not in MAP1 or MAP2,
*                       hence should not have called this routine.
*
*                = -51, then an input argument to one of the processors
*                       had an illegal value (that processor returns
*                       INFO in the range -1 to -13), but all input
*                       arguments were ok on this processor
*                       
*                = -52, then IDIR, ISIGN, N1, N2, N3, MAP1 and/or MAP2
*                       are not the same on all processors in MAP1
*                       and/or MAP2.
*
* ----------------------------------------------------------------------

      EXTERNAL          PRFFT, PRFFTI, PFFT1, PFFT1I, MXINIT

************************************************************************

*     Initialize MX communication library
      CALL MXINIT()

*     Check input data.

      CALL CHK3D( 'PRFFT3', IDIR, ISIGN, N1, N2, N3, X, MAP1, MAP2,
     $            LWORK, WORK, LIWORK, IWORK, INFO )

      IF( INFO .NE. 0 )
     $  RETURN

************************************************************************

      CALL PFFT( PRFFT, PRFFTI, PFFT1, PFFT1I, 2, IDIR,
     $           ISIGN, N1, N2, N3, X, MAP1, MAP2,
     $           LWORK, WORK, LIWORK, IWORK, INFO )

      RETURN
      END
