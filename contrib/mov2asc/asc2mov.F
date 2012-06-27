      Program asc2mov
      Implicit NONE
C
C     Reads formatted (ASCII) movecs file, writes binary version
C
C     Usage: asc2mov guess_nbf infile outfile
C
C     guess_nbf must be at least as large as the NBF in the file to be
C     translated.
C
C     The format of the ASCII file is exactly as the companion program
C     mov2asc writes it. 
C
C     $Id$
C
#include "mafdecls.fh"
c
c     Temporary routine
c
      character*(255) data, bfilename, afilename
      character*(255) title       ! Returns title of job that created vectors
      character*(255) basis_name  ! Returns name of basis set
      integer nbf               ! Returns no. of functions in basis
      integer nsets             ! Returns no. of functions in each set
      integer ldnmo             ! Inputs size of nmo
      parameter(ldnmo = 10)
      integer nmo(ldnmo)        ! Returns no. of vectors in each set
c
      integer binlu,asclu      ! Unit no. for writeing
      parameter (binlu = 67,asclu = 68)   ! These need to be managed !!!
      integer lentit
      integer lenbas
      integer ok, i, guess_nbf, NArgs, jset, j, k_vecs, l_vecs
      character*32 geomsum, basissum, bqsum
      character*26 date
      character*20 scftype20
C
#if defined(CRAY)
      integer*4 length,ierror,iarg
#else
      Integer IArgc
#if (__GNUC__ < 4)
      External IArgc
#endif
#endif
C
      double precision energy, enrep
C     Interpret command line
C
#ifdef CRAY
      integer*4 ipxfargc
      external ipxfargc
      nargs =ipxfargc() 
#else
      NArgs = IArgc()
#endif
      If ( NArgs .ne. 3) then
         Write (0, *) 'Usage: asc2mov guess_nbf infile outfile'
         Call Exit(1)
      EndIf
#ifdef CRAY
      iarg=1
      Call pxfGetArg(iarg, data, length, ierror)
#else
      Call GetARg(1, data)
#endif
      Read(Data, '(I10)') Guess_NBF
#ifdef CRAY
      iarg=2
      Call pxfGetArg(iarg, afilename, 255, ok)
      iarg=3
      Call pxfGetArg(iarg, bfilename, 255, ok)
#else
      Call GetARg(2, afilename)
      Call GetARg(3, bfilename)
#endif
C
c$$$      Write (6, 9000) Guess_NBF, Bfilename, Afilename
c$$$ 9000 Format('Arguments ', I10, 5X, '+', A, '+', '+', A, '+')
c$$$C
      If ( .NOT. MA_Init(Mt_Dbl, Guess_NBF, 0)) Call ErrQuit(
     $   'Unable to initialize memory', Guess_NBF)
c
      open(binlu, status='unknown', form='unformatted', file=bfilename)
      open(asclu, status='old', form='formatted', file=afilename)
C
      Read(AscLu, *)
c
      read(asclu,'(a32)') basissum
      read(asclu,'(a32)') geomsum
      read(asclu,'(a32)') bqsum
      read(asclu,'(a20)') scftype20
      read(asclu,'(a26)') date
      write(binlu) basissum, geomsum, bqsum, scftype20, date
c
      read(asclu,'(a20)') scftype20
      write(binlu) scftype20
C
      read(asclu, '(I10)') LenTit
      write(binlu) lentit
C
      title = ' '
      Read(AscLu, '(A)') Title(1:LenTit)
      write(binlu) title(1:lentit)
C
      Read (AscLu, '(I10)') LenBas
      if (len(basis_name) .lt. lenbas) call errquit
     $   ('movecs_write_header: basis_name too short', lenbas)
      write(binlu) lenbas
C
      basis_name = ' '
      Read( AscLu, '(A)') Basis_Name(1:LenBas)
      write(binlu) basis_name(1:lenbas)
C
      Read( AscLu, '(I10)') NSets
      write(binlu) nsets
C
      Read( AscLu, '(I10)') NBF
      write(binlu) nbf
C
      if (nsets .gt. ldnmo) then
         Write(6,*) ' movecs_write_header: ldnmo too small ',
     $      nsets, ldnmo
         close(binlu)
         close(asclu)
         Stop
      endif
      Read( AscLu, '(7(I10,X))') (nmo(i),i=1,nsets)
      write(binlu) (nmo(i),i=1,nsets)
c
      If ( Guess_NBF .lt. NBF ) Call ErrQuit(
     $   'Guessed too small for NBF.  Actual is', NBF)
      If ( .NOT. MA_Push_Get(MT_Dbl, NBF, 'temporary vector', l_vecs,
     $   k_vecs)) Call ErrQuit('Failed to allocate temp vector', NBF)
C
      do jset = 1, Nsets
         Read(Asclu, '(3E25.15)') (dbl_mb(k_vecs+j),j=0,nbf-1)
         write(binlu) (dbl_mb(k_vecs+j),j=0,nbf-1) ! Occupation numbers
C
         Read(Asclu, '(3E25.15)') (dbl_mb(k_vecs+j),j=0,nbf-1)
         write(binlu) (dbl_mb(k_vecs+j),j=0,nbf-1) ! Eigenvalues
C
         do i = 1, nmo(jset)
            Read(Asclu, '(3E25.15)') (dbl_mb(k_vecs+j),j=0,nbf-1)
            write(binlu) (dbl_mb(k_vecs+j),j=0,nbf-1) ! An eigenvector
         enddo
      enddo
C
      read(Asclu,'(2E25.15)',err=901,end=901) energy, enrep
      write(binlu) energy, enrep
      goto 902
 901  write(binlu) 0.0,0.0
 902  continue
c
      close(binlu)
      Close(AscLu)
C
      if (.not. ma_pop_stack(l_vecs)) call errquit
     $   ('movecs_write: pop failed', l_vecs)
C
      Stop
      End
      Subroutine ErrQuit(Msg, Num)
      Implicit NONE
      Character*(*) Msg
      Integer Num
      Write (0, '(A,1X, I10)') Msg, Num
      Call Exit(1)
      Return
      End
