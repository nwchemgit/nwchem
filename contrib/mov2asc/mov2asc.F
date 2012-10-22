      Program mov2asc
      Implicit NONE
C
C     Reads binary movecs file, writes formatted (ASCII) version
C
C     Usage: mov2asc guess_nbf infile outfile
C
C     guess_nbf must be at least as large as the NBF in the file to be
C     translated.
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
      integer binlu,asclu      ! Unit no. for reading
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
c
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
         Write (0, *) 'Usage: mov2asc guess_nbf infile outfile'
         Call Exit(1)
      EndIf
#ifdef CRAY
      iarg=1
      Call pxfGetArg(iarg, data, 255, ierror)
#else
      Call GetARg(1, data)
#endif
      Read(Data, '(I10)') Guess_NBF
#ifdef CRAY
      iarg=2
      Call pxfGetArg(iarg, bfilename, 255, ok)
      iarg=3
      Call pxfGetArg(iarg, afilename, 255, ok)
#else
      Call GetARg(2, bfilename)
      Call GetARg(3, afilename)
#endif
C
c$$$      Write (6, 9000) Guess_NBF, Bfilename, Afilename
c$$$ 9000 Format('Arguments ', I10, 5X, '+', A, '+', '+', A, '+')
c$$$C
      If ( .NOT. MA_Init(Mt_Dbl, Guess_NBF, 0)) Call ErrQuit(
     $   'Unable to initialize memory', Guess_NBF)
c
      open(binlu, status='old', form='unformatted', file=bfilename)
      open(asclu, status='unknown', form='formatted', file=afilename)
C
      write (AscLu, 9010)
 9010 Format('# This is an NWChem movecs file translated by mov2asc')
C
C Record on convergence and scftype
C
      read(binlu,err=1001) basissum, geomsum, bqsum, scftype20, date
      goto 1003
 1001 continue
      rewind(binlu,err=1002)
 1002 continue
      bqsum = ' '
      read(binlu,err=1003) basissum, geomsum, scftype20, date
 1003 continue
      write(asclu,'(a32)') basissum
      write(asclu,'(a32)') geomsum
      write(asclu,'(a32)') bqsum
      write(asclu,'(a20)') scftype20
      write(asclu,'(a26)') date
c
c
      read(binlu) scftype20
      write(asclu,'(a20)') scftype20
C
C
      read(binlu) lentit
      if (len(title) .lt. lentit) call errquit
     $   ('movecs_read_header: title too short', lentit)
      write(asclu, '(I10)') LenTit
C
      title = ' '
      read(binlu) title(1:lentit)
      Write(AscLu, '(A)') Title(1:LenTit)
C
      read(binlu) lenbas
      if (len(basis_name) .lt. lenbas) call errquit
     $   ('movecs_read_header: basis_name too short', lenbas)
      Write (AscLu, '(I10)') LenBas
C
      basis_name = ' '
      read(binlu) basis_name(1:lenbas)
      Write( AscLu, '(A)') Basis_Name(1:LenBas)
C
      read(binlu) nsets
      Write( AscLu, '(I10)') NSets
C
      read(binlu) nbf
      Write( AscLu, '(I10)') NBF
C
      if (nsets .gt. ldnmo) then
         write(6,*) ' movecs_read_header: ldnmo too small ',
     $      nsets, ldnmo
         close(binlu)
         close(asclu)
         Stop
      endif
      read(binlu) (nmo(i),i=1,nsets)
      Write( AscLu, '(7(I10,X))') (nmo(i),i=1,nsets)
c
      If ( Guess_NBF .lt. NBF ) Call ErrQuit(
     $   'Guessed too small for NBF.  Actual is', NBF)
      If ( .NOT. MA_Push_Get(MT_Dbl, NBF, 'temporary vector', l_vecs,
     $   k_vecs)) Call ErrQuit('Failed to allocate temp vector', NBF)
C
      do jset = 1, Nsets
         read(binlu) (dbl_mb(k_vecs+j),j=0,nbf-1) ! Occupation numbers
         Write(Asclu, '(3E25.15)') (dbl_mb(k_vecs+j),j=0,nbf-1)
C
         read(binlu) (dbl_mb(k_vecs+j),j=0,nbf-1) ! Eigenvalues
         Write(Asclu, '(3E25.15)') (dbl_mb(k_vecs+j),j=0,nbf-1)
C
         do i = 1, nmo(jset)
            read(binlu) (dbl_mb(k_vecs+j),j=0,nbf-1) ! An eigenvector
            Write(Asclu, '(3E25.15)') (dbl_mb(k_vecs+j),j=0,nbf-1)
         enddo
      enddo
C
      read(binlu,err=901,end=901) energy, enrep
      write(Asclu,'(2E25.15)') energy, enrep
      goto 902
 901  write(Asclu,'(2E25.15)') 0.0,0.0
 902  continue
c
      close(binlu)
      Close(AscLu)
C
      if (.not. ma_pop_stack(l_vecs)) call errquit
     $   ('movecs_read: pop failed', l_vecs)
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
