      Program Stupid_Energy
C$Id$
      Implicit NONE
C      
      Integer NBF, NOcc, NVir, NFit
      Character*(*) Vectors
C     -------------------- biphenyl ------------------------------
      Parameter (NBF = 204, NOcc = 73, NVir = NBF-NOcc, NFit = 204)
      Parameter (Vectors = 'bip.movecs')
c$$$C     -------------------------- He16 ----------------------------
c$$$      Parameter (NBF = 32, NOcc = 16, NVir = NBF-NOcc, NFit = 272)
c$$$      Parameter (Vectors = 'movecs')
C
      Double precision IntI(NVir, NFit), IntJ(NVir, NFit),
     $   A(NVir*NVIr), B(NVir*NVir), Evals(NBF)
C
      Integer I, J, Junk, AB, AA, BB
      double precision ChkSum, DenIJ, Den, EPair
C
      External SASUM, ISAMAX, SDot
      Double precision SASUM, SDot
      Integer ISAMAX
C
#include "numerical_constants.fh"
C
      Call OpenRec(51, NVir*NFit)
      Call MOVecs_Read(Vectors, Evals)
      Write (6, 8000) EVals
 8000 Format((5F13.6))
C
      Do I = 1, NOcc
         Do J = 1, I
C
            Call GetRec(51, I, NVir*NFit, IntI )
            Call GetRec(51, J, NVir*NFit, IntJ )
C
            Call SGEMM('N', 'T', NVir, NVir, NFit,
     $         FP_One, IntI, NVir, IntJ, NVir,
     $         FP_Zero, A, NVir)
C
            ChkSum = SASUM(NVir*NVir, A, 1)
            Junk = ISAMAX(NVir*NVir, A, 1)
            Write (6, 9000) 'A', I, J, ChkSum, A(Junk)
 9000       Format(1X, A, ' block ', I5, I5, ' checksum = ', F30.15,
     $         ' max = ', F30.15)
C
            Call DGECOP('T', NVir, NVir, A, NVir, B, NVir)
            Call SScal(NVir*NVir, -FP_One, B, 1)
            Call SAXPY(NVir*NVir, FP_Two, A, 1, B, 1)
C
            ChkSum = SASUM(NVir*NVir, B, 1)
            Junk = ISAMAX(NVir*NVir, B, 1)
            Write (6, 9000) 'B', I, J, ChkSum, B(Junk)
C
            DenIJ = Evals(I) + Evals(J)
            AB = 1
            Do BB = 1, NVir
               Do AA = 1, NVir
                  Den = DenIj - Evals(NOcc+AA) - EVals(NOcc+BB)
                  B(AB) = B(AB) / Den
                  AB = AB + 1
               EndDo
            EndDo
C
            ChkSum = SASUM(NVir*NVir, B, 1)
            Junk = ISAMAX(NVir*NVir, B, 1)
            Write (6, 9000) 'Denom B', I, J, ChkSum, B(Junk)
C
            EPair = SDot(NVir*NVir, A, 1, B, 1)
            Write (6, 9010) I, J, EPair
 9010       Format(1X, 'Pair energy ', I5, I5, F30.15)
C
         EndDo
      EndDo
C
C
C
      Stop
      End
C
C
C
C
      Subroutine OpenRec(LU, NElem)
      Implicit NONE
      Integer LU, NElem
C
      Open (UNIT=LU, ACCESS='DIRECT', RECL=NElem*8)
      Return
      End
C
C
C
      Subroutine GetRec(LU, RecNr, NElem, A)
      Implicit NONE
      Integer LU, RecNr, NElem
      Double Precision A(NElem)
C
      Integer Junk
      Double precision ChkSum
C
      External SASUM, ISAMAX
      Double precision SASUM
      Integer IASMAX
C
      Read (LU, REC=RecNr) A
C
      ChkSum = SASUM(NElem, A, 1)
      Junk = ISAMAX(NElem, A, 1)
      Write (6, 9000) RecNr, ChkSum, A(Junk)
 9000 Format(1X, 'Record ', I5, ' checksum = ', F30.15,
     $   ' max = ', F30.15)
      Return
      End
C
C
C
      Subroutine movecs_read(filename, evals)
      implicit none
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "mafdecls.fh"
c
      character*(*) filename
      double precision evals(*)
c
      integer nsets             ! No. of sets of vectors
      integer nbf               ! No. of functions in basis
      integer nmo(2)            ! No. of vectors in each set
      integer ok, jset, i, j
      integer l_vecs, k_vecs
      integer unitno
      parameter (unitno = 67)
c     
      ok = 0
      open(unitno, status='old', form='unformatted', file=filename,
     $   err=1000)
c
c     Skip over uninteresting bits of the header
c
      read(unitno, err=1001, end=1001) ! lentit
      read(unitno, err=1001, end=1001) ! title
      read(unitno, err=1001, end=1001) ! lenbas
      read(unitno, err=1001, end=1001) ! basis_name
      read(unitno, err=1001, end=1001) nsets
      read(unitno, err=1001, end=1001) nbf
      read(unitno, err=1001, end=1001) (nmo(i),i=1,nsets)
C
      read(unitno, err=1001, end=1001) ! (occ(j),j=1,nbf)
      read(unitno, err=1001, end=1001) (evals(j),j=1,nbf)
 10   close(unitno,err=1002)
      ok = 1
c
      return
c
 1000 write(6,*) ' movecs_read: failed to open ', filename
      ok = 0
      goto 10
c
 1001 write(6,*) ' movecs_read: failing reading from ', filename
      ok = 0
      close(unitno,err=1002)
      goto 10
c
 1002 write(6,*) ' movecs_read: failed to close', filename
      ok = 0
      goto 10
c
      end
