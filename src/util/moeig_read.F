c
c $Id$
c
c This routine was essentially stolen from movecs_read, but only read
c the eigenvalues.
c
      logical function moeig_read(filename, iset, evals )
      implicit none
#include "global.fh"
#include "msgtypesf.h"
#include "mafdecls.fh"
#include "msgids.fh"
#include "inp.fh"
#include "util.fh"
c
      character*(*) filename
      integer iset              ! No. (1,2) of set of vectors to read
      double precision evals(*) ! Must be at least nbf long (not nmo)
c
      integer nsets             ! No. of sets of vectors
      integer nbf               ! No. of functions in basis
      integer nmo(2)            ! No. of vectors in each set
      integer ok, jset, i, j
      integer mitob1
      integer unitno
      parameter (unitno = 67)
c
      mitob1=MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_sync()
      ok = 0
      if (ga_nodeid() .eq. 0) then
         open(unitno, status='old', form='unformatted', file=filename,
     $        err=1000)
c
c     Skip over uninteresting bits of the header
c
         read(unitno, err=1001, end=1001) ! convergence info
         read(unitno, err=1001, end=1001) ! scftype
         read(unitno, err=1001, end=1001) ! lentit
         read(unitno, err=1001, end=1001) ! title
         read(unitno, err=1001, end=1001) ! lenbas
         read(unitno, err=1001, end=1001) ! basis_name
         read(unitno, err=1001, end=1001) nsets
         read(unitno, err=1001, end=1001) nbf
         read(unitno, err=1001, end=1001) (nmo(i),i=1,nsets)
c
c     Skip over unwanted sets
c
         do jset = 1, iset-1
            read(unitno, err=1001, end=1001)
            read(unitno, err=1001, end=1001)
            do i = 1, nmo(jset)
               read(unitno, err=1001, end=1001)
            enddo
         enddo
         read(unitno, err=1001, end=1001) ! occ
         read(unitno, err=1001, end=1001) (evals(j),j=1,nbf)
         do i = 1, nmo(iset)
             read(unitno, err=1001, end=1001)  ! movecs
         enddo
 9       close(unitno,err=1002)
         ok = 1
      endif
c
 10   continue
      call ga_sync()
      call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, mitob1, 0) ! Propagate status
      if (ok .eq. 1) then
         call ga_brdcst(Msg_Vec_EVal+MSGDBL, evals,
     $        MA_sizeof(MT_DBL,nbf,MT_BYTE), 0)
      endif
c
      moeig_read = ok .eq. 1
      if (ga_nodeid() .eq. 0 .and. moeig_read .and.
     $     util_print('vectors i/o', print_default)) then
         write(6,22) filename(1:inp_strlen(filename))
 22      format(/' Read molecular orbitals from ',a/)
         call util_flush(6)
      endif
      return
c
 1000 write(6,*) ' moeig_read: failed to open ',
     $     filename(1:inp_strlen(filename))
      call util_flush(6)
      ok = 0
      goto 10
c
 1001 write(6,*) ' moeig_read: failing reading from ',
     $     filename(1:inp_strlen(filename))
      call util_flush(6)
      ok = 0
      close(unitno,err=1002)
      goto 10
c
 1002 write(6,*) ' moeig_read: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(6)
      ok = 0
      goto 10
c
      end
