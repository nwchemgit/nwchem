      logical function movecs_read_header(filename, title, basis_name, 
     $     nbf, nsets, nmo, ldnmo)
C$Id: vectors.F,v 1.15 1995-04-19 22:28:41 gg502 Exp $
      implicit none
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "msgids.fh"
#include "cscfps.fh"
#include "inp.fh"
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
      if (oscfps) call pstat_on(ps_vecio)
      ok = 0
      if (ga_nodeid() .eq. 0) then
         open(unitno, status='old', form='unformatted', file=filename,
     $        err=1000)
         read(unitno, err=1001, end=2001) lentit
         if (len(title) .lt. lentit) call errquit
     $        ('movecs_read_header: title too short', lentit)
         title = ' '
         read(unitno, err=1001, end=2001) title(1:lentit)
         read(unitno, err=1001, end=2001) lenbas
         if (len(basis_name) .lt. lenbas) call errquit
     $        ('movecs_read_header: basis_name too short', lenbas)
         basis_name = ' '
         read(unitno, err=1001, end=2001) basis_name(1:lenbas)
         read(unitno, err=1001, end=2001) nsets
         read(unitno, err=1001, end=2001) nbf
         if (nsets .gt. ldnmo) then
            write(6,*) ' movecs_read_header: ldnmo too small ',
     $           nsets, ldnmo
            close(unitno, err=1002)
            goto 10
         endif
         read(unitno, err=1001, end=2001) (nmo(i),i=1,nsets)
         close(unitno, err=1002)
         ok = 1
      endif
c
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, mitob(1), 0) ! Propagate status
      movecs_read_header = ok .eq. 1
      if (.not. movecs_read_header) return
c
      call util_char_ga_brdcst(Msg_Vec_Ttl, title)
      call util_char_ga_brdcst(Msg_Vec_Nam, basis_name)
      call ga_brdcst(Msg_Vec_NBF+MSGINT, nbf, mitob(1), 0)
      call ga_brdcst(Msg_Vec_Set+MSGINT, nsets, mitob(1), 0)
      call ga_brdcst(Msg_Vec_NMO+MSGINT, nmo, mitob(nsets), 0)
c
      if (oscfps) call pstat_off(ps_vecio)
c
      return
c
 1000 write(6,*) ' movecs_read_header: failed to open ',
     $     filename(1:inp_strlen(filename))
      call util_flush(6)
      ok = 0
      goto 10
c
 1001 write(6,*) ' movecs_read_header: failing reading from ',
     $     filename(1:inp_strlen(filename))
      call util_flush(6)
      ok = 0
      close(unitno, err=1002)
      goto 10
c
 2001 write(6,*) ' movecs_read_header:eof: failing reading from '
     &    , filename(1:inp_strlen(filename))
      call util_flush(6)
      ok = 0
      close(unitno, err=1002)
      goto 10
c
 1002 write(6,*) ' movecs_read_header: failed to close',
     $     filename(1:inp_strlen(filename))
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
#include "msgids.fh"
#include "cscfps.fh"
#include "inp.fh"
#include "util.fh"
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
      l_vecs = -1               ! An invalid MA handle
c
      if (oscfps) call pstat_on(ps_vecio)
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
 9       close(unitno,err=1002)
         ok = 1
      endif
c
 10   continue
      if (l_vecs .ne. -1) then
	 if (.not. ma_pop_stack(l_vecs)) call errquit
     $      ('movecs_read: pop failed', l_vecs)
      endif
      call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, mitob(1), 0) ! Propagate status
      if (ok .eq. 1) then
	 call ga_brdcst(Msg_Vec_NBF+MSGINT, nbf, mitob(1), 0)
	 call ga_brdcst(Msg_Vec_EVal+MSGDBL, evals, mdtob(nbf), 0)
	 call ga_brdcst(Msg_Vec_Occ+MSGDBL, occ,   mdtob(nbf), 0)
      endif
c
      movecs_read = ok .eq. 1
      if (ga_nodeid() .eq. 0 .and. movecs_read .and.
     $     util_print('vectors i/o', print_default)) then
         write(6,22) filename(1:inp_strlen(filename))
 22      format(/' Read molecular orbitals from ',a/)
         call util_flush(6)
      endif
      if (oscfps) call pstat_off(ps_vecio)
      return
c
 1000 write(6,*) ' movecs_read: failed to open ',
     $     filename(1:inp_strlen(filename))
      call util_flush(6)
      ok = 0
      goto 10
c
 1001 write(6,*) ' movecs_read: failing reading from ',
     $     filename(1:inp_strlen(filename))
      call util_flush(6)
      ok = 0
      close(unitno,err=1002)
      goto 10
c
 1002 write(6,*) ' movecs_read: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(6)
      ok = 0
      goto 10
c
      end
      logical function movecs_gather(filename, iset, NWant, Want,
     $   occ, evals, g_vecs)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
c
      character*(*) filename        ! [in] Name of movecs file
      integer iset                  ! [in] No. (1,2) of set of vectors to read
      integer NWant                 ! [in] length of gather vector
      integer Want(NWant)           ! [in] list of MOs to be gathered
      double precision occ(NWant)   ! [out] Occupation numbers
      double precision evals(NWant) ! [out] Eigenvalues
      integer g_vecs                ! [in] Handle for GA to hold results
c
      integer nsets                 ! No. of sets of vectors
      integer nbf                   ! No. of functions in basis
      integer ldnmo                 ! Size of nmo
      parameter (ldnmo = 2)
      integer nmo(ldnmo)            ! No. of vectors in each set
      integer i, l_FullOcc, k_FullOcc, l_FullEval, k_FullEval, G_Tmp
      Character*32 Title, Basis
      integer MinChunk              ! For data distirbution in GA
      parameter (MinChunk = -1)     ! Let GA determine even distrib.
C
      Logical MOVecs_Read_Header, MOVecs_Read
      External MOVecs_Read_Header, MOVecs_Read
c     
      If ( .NOT. MOVecs_Read_Header(Filename, title, basis, nbf,
     $   nsets, nmo, ldnmo) ) Call ErrQuit(
     $   'MOVecs_Gather: failed to read movecs file', 1)
c
C     Allocate memory to be used while reading info
C
      if (.not. ma_push_get(mt_dbl,nmo(iset),'movecs_gather',
     $   l_FullOcc, k_FullOcc)) call errquit(
     $   'movecs_gather: ma failed', nmo(iset))
C
      if (.not. ma_push_get(mt_dbl,nmo(iset),'movecs_gather',
     $   l_FullEval, k_FullEval)) call errquit(
     $   'movecs_gather: ma failed', nmo(iset))
C
      if ( .NOT. GA_Create( MT_Dbl, NBF, NMO(ISet), 'Full vectors',
     $   MinChunk, MinChunk, G_Tmp) ) Call ErrQuit(
     $   'MOVecs_Gather: ga_create failed', NBF*NMO(ISet) )
C
C     Read the full data
C
      If ( .NOT. MOVecs_Read(Filename, ISet, Dbl_MB(k_FullOcc),
     $   Dbl_MB(k_FullEval), G_Tmp) ) Call ErrQUit(
     $   'MOVecs_Gather: failed reading movecs file', 2)
C
C     Gather the data into the requested subset
C
      Do i = 1, NWant
         Occ(i)   = Dbl_MB( k_FullOcc  + Want(i) - 1)
         Evals(i) = Dbl_MB( k_FullEval + Want(i) - 1)
         Call GA_Copy_Patch('N', G_Tmp, 1, NBF, Want(I), Want(I),
     $      G_Vecs, 1, NBF, I, I)
      EndDo
C
C     Free temporaries
C
      If ( .NOT. GA_Destroy(G_Tmp)) Call ErrQuit(
     $   'MOVecs_Gather: ga_destroy failed', G_Tmp)
      If (.NOT. MA_Pop_Stack(l_FullEVal)) Call ErrQuit(
     $   'movecs_gather: pop failed', l_FullEval)
      If (.NOT. MA_Pop_Stack(l_FullOcc)) Call ErrQuit(
     $   'movecs_gather: pop failed', l_FullOcc)
C
      MOVecs_Gather = .TRUE.
C
      Return
      End
      logical function movecs_write(filename, title, basis_name, 
     $     nbf, nsets, nmo, occ, ldocc, evals, ldevals, g_vecs)
      implicit none
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "mafdecls.fh"
#include "inp.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "util.fh"
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
      if (oscfps) call pstat_on(ps_vecio)
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
         write(unitno, err=1001) (nmo(i),i=1,nsets)
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
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, mitob(1), 0) ! Propagate status
c
      movecs_write = ok .eq. 1
      if (ga_nodeid() .eq. 0 .and.
     $     util_print('vectors i/o', print_default)) then
         write(6,22) filename(1:inp_strlen(filename))
 22      format(/' Wrote molecular orbitals to ',a/)
         call util_flush(6)
      endif
      if (oscfps) call pstat_off(ps_vecio)
      return
c
 1000 write(6,*) ' movecs_write: failed to open ',
     $     filename(1:inp_strlen(filename))
      call util_flush(6)
      ok = 0
      goto 10
c
 1001 write(6,*) ' movecs_write: failing writing to ',
     $     filename(1:inp_strlen(filename))
      call util_flush(6)
      ok = 0
      close(unitno,err=1002)
      goto 10
c
 1002 write(6,*) ' movecs_write: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(6)
      ok = 0
      goto 10
c
      end
      subroutine movecs_swap(rtdb,g_vecs,occ,eval)
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
c
      integer rtdb, g_vecs
      double precision occ(*), eval(*)
c
      double precision tmp
      integer type, nelem, handle, index
      integer l_veci, k_veci, l_vecj, k_vecj, dim1, dim2
      integer pair, i, j
      logical oprint
c
      oprint = util_print('mo guess', print_default)
      if (.not. rtdb_ma_get(rtdb, 'scf:swap', type, nelem, 
     $     handle)) return
      if (.not. ma_get_index(handle,index)) call errquit
     $     ('movecs_swap: get index failed', handle)
      if (mod(nelem,2).ne.0 .or. type.ne.mt_int) call errquit
     $     ('movecs_swap: invalid list of pairs',nelem)
c
      call ga_sync()
      if (ga_nodeid() .eq. 0) then
         if (oprint) write(6,*)
         call ga_inquire(g_vecs, type, dim1, dim2)
         if (.not. ma_push_get(mt_dbl, dim1, 'swap', l_veci, k_veci))
     $        call errquit('movecs_swap: no scratch space', dim1)
         if (.not. ma_push_get(mt_dbl, dim1, 'swap', l_vecj, k_vecj))
     $        call errquit('movecs_swap: no scratch space', dim1)
         do pair = 1, nelem, 2
            i = int_mb(index+pair-1)
            j = int_mb(index+pair  )
            if (oprint) write(6,1) i, j
 1          format(' Swapping orbitals ', 2i5)
            if (i.lt.0 .or. i.gt.dim2) call errquit
     $           ('movecs_swap: invalid vector ', i)
            if (j.lt.0 .or. j.gt.dim2) call errquit
     $           ('movecs_swap: invalid vector ', j)
c     
            call ga_get(g_vecs, 1, dim1, i, i, dbl_mb(k_veci), 1)
            call ga_get(g_vecs, 1, dim1, j, j, dbl_mb(k_vecj), 1)
            call ga_put(g_vecs, 1, dim1, i, i, dbl_mb(k_vecj), 1)
            call ga_put(g_vecs, 1, dim1, j, j, dbl_mb(k_veci), 1)
c
            tmp = occ(i)
            occ(i) = occ(j)
            occ(j) = tmp
            tmp = eval(i)
            eval(i) = eval(j)
            eval(j) = tmp

c     
         enddo
         if (oprint) call util_flush(6)
         if (.not. ma_pop_stack(l_vecj)) call errquit('ms:ma?', 0)
         if (.not. ma_pop_stack(l_veci)) call errquit('ms:ma?', 0)
      endif
c
      call ga_sync
c
      if (.not. rtdb_delete(rtdb, 'scf:swap')) call errquit
     $     ('movecs_swap: delete of swap entry failed', 0)
c
      if (.not. ma_free_heap(handle)) call errquit('ms:ma?', 0)
c
      end
