      logical function movecs_read_header(filename, title, basis_name, 
     $     nbf, nsets, nmo, ldnmo)
      implicit none
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
c
c     Temporary routine
c
      character*(*) filename    ! File to read header from
      character*(*) title       ! Returns title of job that created vectors
      character*(*) basis_name  ! Returns name of basis set
      integer nbf               ! Returns no. of functions in basis
      integer nsets             ! Returns no. of functions in each set
      integer ldnmo             ! Inputs size of nmo
      integer nmo(ldnmo)        ! Returns no. of vectors in each set
c
      integer unitno            ! Unit no. for reading
      parameter (unitno = 67)   ! These need to be managed !!!
      integer lentit
      integer lenbas
      integer ok, i
c
      ok = 0
      if (ga_nodeid() .eq. 0) then
         open(unitno, status='old', form='unformatted', file=filename,
     $        err=1000)
         read(unitno, err=1001, end=1001) lentit
         if (len(title) .lt. lentit) call errquit
     $        ('movecs_read_header: title too short', lentit)
         title = ' '
         read(unitno, err=1001, end=1001) title(1:lentit)
         read(unitno, err=1001, end=1001) lenbas
         if (len(basis_name) .lt. lenbas) call errquit
     $        ('movecs_read_header: basis_name too short', lenbas)
         basis_name = ' '
         read(unitno, err=1001, end=1001) basis_name(1:lenbas)
         
         read(unitno, err=1001, end=1001) nsets
         read(unitno, err=1001, end=1001) nbf
         if (nsets .gt. ldnmo) then
            write(6,*) ' movecs_read_header: ldnmo too small ',
     $           nsets, ldnmo
            close(unitno, err=1002)
            goto 10
         endif
         read(unitno, err=1001, end=1001) (nmo(i),i=1,nsets)
         close(unitno, err=1002)
         ok = 1
      endif
c
 10   call ga_brdcst(1+MSGINT, ok, mitob(1), 0) ! Propagate status
      movecs_read_header = ok .eq. 1
      if (.not. movecs_read_header) return
c
*     call cbrdcst(1, title, 0)
*     call cbrdcst(2, basis_name, 0)
      call ga_brdcst(4+MSGINT, nbf, mitob(1), 0)
      call ga_brdcst(5+MSGINT, nsets, mitob(1), 0)
      call ga_brdcst(6+MSGINT, nmo, mitob(nsets), 0)
c
      return
c
 1000 write(6,*) ' movecs_read_header: failed to open ', filename
      call util_flush(6)
      ok = 0
      goto 10
c
 1001 write(6,*) ' movecs_read_header: failing reading from ', filename
      call util_flush(6)
      ok = 0
      close(unitno, err=1002)
      goto 10
c
 1002 write(6,*) ' movecs_read_header: failed to close', filename
      call util_flush(6)
      ok = 0
      goto 10
c
      end
      logical function movecs_read(filename, iset, occ, evals, g_vecs)
      implicit none
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "mafdecls.fh"
c
      character*(*) filename
      integer iset              ! No. (1,2) of set of vectors to read
      double precision occ(*)
      double precision evals(*)
      integer g_vecs
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
      if (ga_nodeid() .eq. 0) then
         open(unitno, status='old', form='unformatted', file=filename,
     $        err=1000)
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
c
         if (.not. ma_push_get(mt_dbl,nbf,'movecs_read',l_vecs,k_vecs))
     $        call errquit('movecs_read: ma failed', nbf)
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
         read(unitno, err=1001, end=1001) (occ(j),j=1,nbf)
         read(unitno, err=1001, end=1001) (evals(j),j=1,nbf)
         do i = 1, nmo(iset)
            read(unitno, err=1001, end=1001)
     $           (dbl_mb(k_vecs+j), j=0,nbf-1)
            call ga_put(g_vecs, 1, nbf, i, i, dbl_mb(k_vecs), 1)
         enddo
	 if (.not. ma_pop_stack(l_vecs)) call errquit
     $      ('movecs_read: pop failed', l_vecs)
 9       close(unitno,err=1002)
         ok = 1
      endif
c
 10   call ga_brdcst(1+MSGINT, ok, mitob(1), 0) ! Propagate status
      if (ok .eq. 1) then
	 call ga_brdcst(2+MSGINT, nbf, mitob(1), 0)
	 call ga_brdcst(3+MSGDBL, evals, mdtob(nbf), 0)
	 call ga_brdcst(4+MSGDBL, occ,   mdtob(nbf), 0)
      endif
c
      movecs_read = ok .eq. 1
      return
c
 1000 write(6,*) ' movecs_read: failed to open ', filename
      call util_flush(6)
      ok = 0
      goto 10
c
 1001 write(6,*) ' movecs_read: failing reading from ', filename
      call util_flush(6)
      ok = 0
      close(unitno,err=1002)
      goto 10
c
 1002 write(6,*) ' movecs_read: failed to close', filename
      call util_flush(6)
      ok = 0
      goto 10
c
      end
      logical function movecs_write(filename, title, basis_name, 
     $     nbf, nsets, nmo, occ, ldocc, evals, ldevals, g_vecs)
      implicit none
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "mafdecls.fh"
#include "inp.fh"
c
c     Temporary routine
c
      character*(*) filename    ! File to write to
      character*(*) title       ! Title of job that created vectors
      character*(*) basis_name  ! Name of basis set
      integer nbf               ! No. of functions in basis
      integer nsets             ! No. of sets of vectors
      integer nmo(nsets)        ! No. of vectors in each set
      integer ldocc
      integer ldevals
      double precision occ(ldocc, nsets) ! Occupation numbers of each set
      double precision evals(ldevals, nsets) ! Eigenvalues for each set
      integer g_vecs(nsets)     ! Global array with eigen-vectors
c
      integer unitno            ! Unit no. for writing
      parameter (unitno = 67)   ! These need to be managed !!!
      integer lentit
      integer lenbas
      integer l_vecs, k_vecs
      integer ok, iset, i, j
c
      ok = 0
      if (ga_nodeid() .eq. 0) then
         if (nsets .gt. 2) then
            write(6,*) ' movecs_write: nsets > 2 ', nsets
            goto 10
         endif
         open(unitno, status='unknown', form='unformatted',
     $        file=filename, err=1000)
c
c     Check that read routines are both consistent with this
c
         lentit = inp_strlen(title)
         write(unitno, err=1001) lentit
         write(unitno, err=1001) title(1:lentit)
         lenbas = inp_strlen(basis_name)
         write(unitno, err=1001) lenbas
         write(unitno, err=1001) basis_name(1:lenbas)
         write(unitno, err=1001) nsets
         write(unitno, err=1001) nbf
         write(unitno, err=1001) nmo
         if (.not. ma_push_get(mt_dbl,nbf,'movecs_write',l_vecs,k_vecs))
     $        call errquit('movecs_write: ma failed', nbf)
         do iset = 1, nsets
            write(unitno, err=1001) (occ(j,iset),j=1,nbf)
            write(unitno, err=1001) (evals(j,iset),j=1,nbf)
            do i = 1, nmo(iset)
               call ga_get(g_vecs(iset), 1, nbf, i, i, dbl_mb(k_vecs),1)
               write(unitno, err=1001) (dbl_mb(k_vecs+j), j=0,nbf-1)
            enddo
         enddo
         if (.not. ma_pop_stack(l_vecs))
     $        call errquit('movecs_write: ma pop failed', l_vecs)
         close(unitno,err=1002)
         ok = 1
      endif
c
 10   call ga_brdcst(1+MSGINT, ok, mitob(1), 0) ! Propagate status
c
      movecs_write = ok .eq. 1
      return
c
 1000 write(6,*) ' movecs_write: failed to open ', filename
      call util_flush(6)
      ok = 0
      goto 10
c
 1001 write(6,*) ' movecs_write: failing writing to ', filename
      call util_flush(6)
      ok = 0
      close(unitno,err=1002)
      goto 10
c
 1002 write(6,*) ' movecs_write: failed to close', filename
      call util_flush(6)
      ok = 0
      goto 10
c
      end

      
      
