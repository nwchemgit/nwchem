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
      SUBROUTINE PFFT( FFT1, FFT1I, FFT2, FFT2I,
     $                 ITYPE, IDIR, ISIGN, N1, N2, N3, X, MAP1, MAP2,
     $                 LWORK, WORK, LIWORK, IWORK, INFO )

      INTEGER           ITYPE, IDIR, ISIGN, N1, N2, N3, LWORK, LIWORK,
     $                  INFO
      INTEGER           MAP1(N3), MAP2(N2), IWORK(LIWORK)
      DOUBLE PRECISION  X(*), WORK(LWORK)

      EXTERNAL          FFT1I, FFT1, FFT2I, FFT2

*-----------------------------------------------------------------------
*
* DESCRIPTION
* ===========
*  FFT1 ...... (Input) EXTERNAL SUBROUTINE
*              If ITYPE = 1 or 3, then:
*                 The name of external routine used to to do all of the
*                 1d ffts
*
*              If ITYPE = 2, then,
*                 The name of the external routine used to to do the fft
*                 in the first first, N1 (N2 if N1 = 1), dimension.
*                 In particular, FFT1 is the name of the routine to
*                 do the 1d real-to-complex fft, IDIR = +1, or the
*                 complex-to-real fft, IDIR = -1.
*
*  FFT1I ..... (Input) EXTERNAL SUBROUTINE
*              Name of external routine used to initialize FFT1.
*
*  FFT2 ...... (Input) EXTERNAL SUBROUTINE
*              If ITYPE = 1 or 3, then:
*                 The name of external routine used to to do all of the
*                 1d ffts, hence FFT2 is the same as FFT1.
*
*              If ITYPE = 2, then:
*                 The name of external routine used to do 1d ffts in
*                 the second and, if present, the third dimension.  In
*                 particular, this is the name of the routine for
*                 doing 1d complex-to-complex ffts.
*
*  FFT2I ..... (Input) EXTERNAL SUBROUTINE
*              Name of external routine used to initialize FFT2
*
*  ITYPE ..... = 1, then real-to-real (i.e., sin or cos) transform,
*                   !! THIS OPTION HAS NOT BEEN COMPLETED AND TESTED !!!
*
*              = 2, then real-to-complex fft, or its inverse 
*                        (complex-to real fft).
*              = 3, then complex-to-complex fft
*
*  IDIR  ..... (Input) INTEGER
*              Only used when ITYPE = 2.
*              = +1 then, real-to-complex fft
*              = -1 then, complex-to-real fft
*
*  ISIGN ..... (Input) INTEGER
*              Sign of exponent in fft = +-1.
*              Only used when ITYPE = 2 or 3.
*
*-----------------------------------------------------------------------

      INTEGER           ME, NELEBK, IDWORK, NVECS1, NVECS2,
     $                  LLWORK, M, NROW, NELEVC, MROW, MSGTYP, MELEVC,
     $                  ICMPLX

      INTEGER           MXNPRC, MXMYND, LSTCNT
      EXTERNAL          MXNPRC, MXMYND, LSTCNT, MXINIT

#ifdef TIMING
      INTEGER           MXSYNC, MXNPRC, IVAL
      DOUBLE PRECISION  MXCLOCK, T0, T1, T2, TIMES(10)
      EXTERNAL          MXCLOCK, MXSYNC, MXNPRC
#endif
************************************************************************

*     Initialize MX communication library
      CALL MXINIT()

#ifdef TIMING
                 IVAL = MXSYNC()
                 T0 = MXCLOCK()
                 T0 = MXCLOCK()
                 T1 = T0
#endif

      INFO = 0

      ME = MXMYND()

*     Determine number of stored data points in the
*     1st dimension, NROW, and the number of DOUBLE PRECISION elements
*     each vector occupies in the first dimesion.

      IF( N1 .GT. 1 ) THEN
        M = N1
      ELSE
        M = N2
      ENDIF
      IF( ITYPE .EQ. 1 ) THEN

*       Sin or Cos transform, input and output both real.

        NROW   = M
        NELEVC = NROW
        ICMPLX = 1

      ELSE IF( ITYPE .EQ. 2 ) THEN

*       real-to-complex or complex-to-real transform, 
*       one of input/output real and the other complex.
*       Assume partially packed storage of complex data.
*       Input and output data of 1d fft actually differ in 
*       size by one or two items, but ignore this.

        NROW   = M / 2 + 1
        NELEVC = 2 * NROW
        ICMPLX = 2

      ELSE IF( ITYPE .EQ. 3 ) THEN

*       Complex-to-complex transform, input and output both complex.

        NROW   = M
        NELEVC = 2 * NROW
        ICMPLX = 2

      ENDIF

*     Determine number of n1-n2 and n1-n3 planes of data
*     that this processor owns.
* 
      NVECS1 = LSTCNT( ME, MAP1, N3 )

      IF( N1 .EQ. 1 ) THEN
        NVECS2 = LSTCNT( ME, MAP2, NROW )
      ELSE
        NVECS2 = LSTCNT( ME, MAP2, N2 )
      ENDIF

#ifdef TESTTRN
         IF( NVECS1 .GT. 0 ) THEN
           CALL SETAC( ITYPE, 1, ME, N1, N2, N3, WORK, MAP1 )
           CALL CHKA12( ITYPE, ME, N1, N2, N3, MAP1, X, WORK, 0 )
         ENDIF
#endif

************************************************************************

*     Skip next part if doing complex-to-real transform.

************************************************************************

      IF( ITYPE .NE. 2  .OR.  IDIR .EQ. 1 ) THEN


      IF( N1 .GT. 1  .AND. NVECS1 .GT. 0 ) THEN

*        Do NVECS1 2-d ffts on X.
*        A 2-d fft is done via column and then row 1-d ffts)


*       N2 * NVECS1 column 1d ffts

        CALL MFFTC1( FFT1, FFT1I, IDIR, ISIGN, NELEVC, N1, N2 * NVECS1,
     $               X, WORK )

#ifdef TESTFFT
                 CALL CHKFFT( ITYPE, ISIGN, ME, N1, N2, N3, X, MAP1,
     $                        1, WORK )
#endif
#ifdef TESTTRN
                 CALL SETAC( ITYPE, 1, ME, N1, N2, N3, WORK, MAP1 )
                 CALL CHKA12( ITYPE, ME, N1, N2, N3, MAP1, X, WORK, 1 )
#endif

*       NROW * NVECS1 row 1d ffts (on NVECS1 NROW-by-N2 2-d arrays
*                                  of data)

        CALL MFFTR1( FFT2, FFT2I, ICMPLX, IDIR, ISIGN, NROW, N2, NVECS1,
     $               X, WORK )

#ifdef TESTFFT
                 CALL CHKFFT( ITYPE, ISIGN, ME, N1, N2, N3, X, MAP1,
     $                        2, WORK )
#endif
      ELSE IF( N1 .EQ. 1  .AND. NVECS1 .GT. 0 ) THEN

*       do NVECS1 column 1-d ffts

        CALL MFFTC1( FFT1, FFT1I, IDIR, ISIGN, NELEVC, N2, NVECS1, X,
     $               WORK )

#ifdef TESTFFT
                 CALL CHKFFT( ITYPE, ISIGN, ME, N1, N2, N3, X, MAP1,
     $                        20, WORK )
#endif

      ENDIF

      ENDIF

************************************************************************

*     End of section to skip when doing complex-to-real transform.

************************************************************************

*     Do block multi-processor transpose.

*     Determine number of rows, MROW, in the block matrix which is to
*     be transposed and the number of elements  in each block, NELEBK,
*     of the block matrix which is to be transposed.

      IF( N1 .EQ. 1 ) THEN
        MROW   = NROW
        NELEBK = 1
      ELSE
        MROW   = N2
        NELEBK = NROW
      ENDIF

      IF( ITYPE .NE. 1 )
     $  NELEBK = 2 * NELEBK

      MSGTYP = 100

*     Store X Transpose in WORK(1:IDWORK-1) and use WORK(IDWORK:) as
*     double precision workspace.

      IDWORK = NELEBK * N3 * NVECS2 + 1

#ifdef TESTTRN
                 CALL SETAC( ITYPE, 1, ME, N1, N2, N3, WORK, MAP1 )
                 CALL CHKA12( ITYPE, ME, N1, N2, N3, MAP1, X, WORK, 2 )
#endif
#ifdef TIMING
                 IVAL = MXSYNC()
                 T2 = MXCLOCK()
                 TIMES(1) = T2 - T1
                 T1 = T2
#endif

      CALL BLKTRN( MROW, N3, NELEBK, X, MAP1, MAP2,
     $             MSGTYP, WORK, WORK( IDWORK ), IWORK, INFO )

#ifdef TIMING
                 IVAL = MXSYNC()
                 T2 = MXCLOCK()
                 TIMES(2) = T2 - T1
                 T1 = T2
#endif

      IF( INFO .NE. 0 ) THEN
        INFO = 2
        RETURN
      ENDIF

#ifdef TESTTRN
      CALL CHKTRN( ITYPE, ME, N1, N2, N3, WORK, MAP2, 1 )
#endif

************************************************************************

*     do MIN( N1, NROW ) * NVECS2 1d ffts on X "transpose"
*     (stored in WORK(1:IDWORK-1)

      LLWORK = LWORK - IDWORK + 1

      IF( NVECS2 .GT. 0  .AND. N1 .GT. 1 ) THEN
 
        CALL MFFTR1( FFT2, FFT2I, ICMPLX, IDIR, ISIGN, NROW, N3, NVECS2,
     $               WORK, WORK( IDWORK ) )

      ELSE IF( NVECS2 .GT. 0  .AND. N1 .EQ. 1 ) THEN

        IF( ITYPE .EQ. 1 ) THEN
          MELEVC = N3
        ELSE
          MELEVC = 2*N3
        ENDIF

        CALL MFFTC1( FFT2, FFT2I, IDIR, ISIGN, MELEVC, N3, NVECS2, WORK,
     $               WORK( IDWORK ) )

      ENDIF

************************************************************************

*     Do inverse block multi-processor transpose.

      MSGTYP = 200

#ifdef TESTTRN
      CALL CHKTRN( ITYPE, ME, N1, N2, N3, WORK, MAP2, 2 )
#endif
#ifdef TIMING
                 IVAL = MXSYNC()
                 T2 = MXCLOCK()
                 TIMES(3) = T2 - T1
                 T1 = T2
#endif

      CALL BLKTRN( N3, MROW, NELEBK, WORK, MAP2, MAP1,
     $             MSGTYP, X, WORK( IDWORK ), IWORK, INFO )

#ifdef TIMING
                 IVAL = MXSYNC()
                 T2 = MXCLOCK()
                 TIMES(4) = T2 - T1
                 T1 = T2
#endif

      IF( INFO .NE. 0 ) THEN
        INFO = 4
        RETURN
      ENDIF

#ifdef TESTFFT
                 IF( NVECS1 .GT. 0 )
     $           CALL CHKFFT( ITYPE, ISIGN, ME, N1, N2, N3, X, MAP1,
     $                        33, WORK )
#endif
#ifdef TESTTRN
         IF( NVECS1 .GT. 0 ) THEN
           CALL SETAC( ITYPE, 1, ME, N1, N2, N3, WORK, MAP1 )
           CALL CHKA12( ITYPE, ME, N1, N2, N3, MAP1, X, WORK, 5 )
         ENDIF
#endif

************************************************************************

*     Do next part if doing complex-to-real transform.
*     Need to do n1 (n2 if n1 = 1) direction complex-to-real fft last

************************************************************************

      IF( ITYPE .EQ. 2  .AND.  IDIR .EQ. -1 ) THEN

      IF( N1 .GT. 1  .AND. NVECS1 .GT. 0 ) THEN

*        Do NVECS1 2-d ffts on X.
*        A 2-d fft is done via row and then column 1-d ffts)

*       NROW * NVECS1 row 1d ffts (on NVECS1 NROW-by-N2 2-d arrays
*                                  of data)

        CALL MFFTR1( FFT2, FFT2I, ICMPLX, IDIR, ISIGN, NROW, N2, NVECS1,
     $               X, WORK )

#ifdef TESTFFT
                 CALL CHKFFT( ITYPE, ISIGN, ME, N1, N2, N3, X, MAP1,
     $                        2, WORK )
#endif
#ifdef TESTTRN
                 CALL SETAC( ITYPE, 1, ME, N1, N2, N3, WORK, MAP1 )
                 CALL CHKA12( ITYPE, ME, N1, N2, N3, MAP1, X, WORK, 6 )
#endif

*       N2 * NVECS1 column 1d ffts

        CALL MFFTC1( FFT1, FFT1I, IDIR, ISIGN, NELEVC, N1, N2 * NVECS1,
     $               X, WORK )

#ifdef TESTFFT
                 CALL CHKFFT( ITYPE, ISIGN, ME, N1, N2, N3, X, MAP1,
     $                        1, WORK )
#endif

      ELSE IF( N1 .EQ. 1  .AND. NVECS1 .GT. 0 ) THEN

*       do NVECS1 column 1-d ffts

        CALL MFFTC1( FFT1, FFT1I, IDIR, ISIGN, NELEVC, N2, NVECS1, X,
     $               WORK )

#ifdef TESTFFT
                 CALL CHKFFT( ITYPE, ISIGN, ME, N1, N2, N3, X, MAP1,
     $                        20, WORK )
#endif

      ENDIF

#ifdef TESTTRN
                 CALL SETAC( ITYPE, 1, ME, N1, N2, N3, WORK, MAP1 )
                 CALL CHKA12( ITYPE, ME, N1, N2, N3, MAP1, X, WORK, 7 )
#endif

      ENDIF

************************************************************************

*     End of section to do when doing complex-to-real transform.

************************************************************************

#ifdef TIMING

                 IVAL = MXSYNC()
                 T2 = MXCLOCK()
                 TIMES(5) = T2 - T1
                 TIMES(6) = T2 - T0
                 IF( ME .EQ. MAP1(1) ) THEN
                   WRITE(*,'(1X,A,3I5, A, I4, A, 3I3)')
     $                ' pfft times, N1,N2,N3=', N1,N2,N3,
     $                ' NPROC=',MXNPRC(), 
     $                ' ISIGN,IDIR,ITYPE= ', ISIGN, IDIR, ITYPE
                   WRITE(*,'(A,F10.4)') ' Pre 1st trans. = ', TIMES(1)
                   WRITE(*,'(A,F10.4)') ' 1st transpose  = ', TIMES(2)
                   WRITE(*,'(A,F10.4)') ' n3  1d fft     = ', TIMES(3)
                   WRITE(*,'(A,F10.4)') ' 2nd transpose  = ', TIMES(4)
                   WRITE(*,'(A,F10.4)') ' Post 2nd trans = ', TIMES(5)
                   WRITE(*,'(A,F10.4)') ' pfft total     = ', TIMES(6)
                   WRITE(*,'(A,F10.4)') ' %time trans    = ',
     $                100.D0*(TIMES(2)+TIMES(4))/TIMES(6)
                 ENDIF
#endif

      RETURN
      END
