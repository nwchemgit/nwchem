      logical function movecs_read_header(filename, title, basis_name, 
     $     scftype, nbf, nsets, nmo, ldnmo)
C$Id: vectors.F,v 1.36 2001-06-29 01:52:41 d3g681 Exp $
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
      character*(*) scftype     ! Returns the SCF type of the vectors
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
      character*20 scftype20
c
      if (oscfps) call pstat_on(ps_vecio)
      ok = 0
      if (ga_nodeid() .eq. 0) then
         open(unitno, status='old', form='unformatted', file=filename,
     $        err=1000)
         read(unitno, err=1001, end=2001) ! SKIP convergence info
         read(unitno, err=1001, end=2001) scftype20
         scftype = scftype20
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
      call util_char_ga_brdcst(Msg_Vec_Ttl, title, 0)
      call util_char_ga_brdcst(Msg_Vec_Nam, basis_name, 0)
      call util_char_ga_brdcst(Msg_Vec_Type, scftype, 0)
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
      double precision occ(*)   ! Must be at least nbf long (not nmo)
      double precision evals(*) ! Must be at least nbf long (not nmo)
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
*            read(unitno, err=1001, end=1001)
*     $           (dbl_mb(k_vecs+j), j=0,nbf-1)
            call sread(unitno, dbl_mb(k_vecs), nbf)
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
      call ga_sync()
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
      logical function movecs_gather(filename, iset, nwant, want,
     $   occ, evals, g_vecs)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
c
      character*(*) filename        ! [in] name of movecs file
      integer iset                  ! [in] no. (1,2) of set of vectors to read
      integer nwant                 ! [in] length of gather vector
      integer want(nwant)           ! [in] list of mos to be gathered
      double precision occ(nwant)   ! [out] occupation numbers
      double precision evals(nwant) ! [out] eigenvalues
      integer g_vecs                ! [in] handle for ga to hold results
c
      integer nsets                 ! no. of sets of vectors
      integer nbf                   ! no. of functions in basis
      integer ldnmo                 ! size of nmo
      parameter (ldnmo = 2)
      integer nmo(ldnmo)            ! no. of vectors in each set
      integer i, l_fullocc, k_fullocc, l_fulleval, k_fulleval, g_tmp
      character*255 title, basis
      integer minchunk              ! for data distirbution in ga
      parameter (minchunk = -1)     ! let ga determine even distrib.
      character*20 scftype
c
      logical movecs_read_header, movecs_read
      external movecs_read_header, movecs_read
c     
      if ( .not. movecs_read_header(filename, title, basis,
     $     scftype, nbf, nsets, nmo, ldnmo) ) call errquit(
     $   'movecs_gather: failed to read movecs file', 1)
c
c     allocate memory to be used while reading info
C     movecs_read expects fullocc and fulleval arrays to be nbf,
C     not nmo(iset)!
c
      if (.not. ma_push_get(mt_dbl,nbf,'movecs_gather',
     $   l_fullocc, k_fullocc)) call errquit(
     $   'movecs_gather: ma failed', nmo(iset))
c
      if (.not. ma_push_get(mt_dbl,nbf,'movecs_gather',
     $   l_fulleval, k_fulleval)) call errquit(
     $   'movecs_gather: ma failed', nmo(iset))
c
      if ( .not. ga_create( mt_dbl, nbf, nmo(iset), 'full vectors',
     $   minchunk, minchunk, g_tmp) ) call errquit(
     $   'movecs_gather: ga_create failed', nbf*nmo(iset) )
c
c     read the full data
c
      if ( .not. movecs_read(filename, iset, dbl_mb(k_fullocc),
     $   dbl_mb(k_fulleval), g_tmp) ) call errquit(
     $   'movecs_gather: failed reading movecs file', 2)
c
c     gather the data into the requested subset
c
      do i = 1, nwant
         occ(i)   = dbl_mb( k_fullocc  + want(i) - 1)
         evals(i) = dbl_mb( k_fulleval + want(i) - 1)
         call ga_copy_patch('n', g_tmp, 1, nbf, want(i), want(i),
     $      g_vecs, 1, nbf, i, i)
      enddo
c
c     free temporaries
c
      if ( .not. ga_destroy(g_tmp)) call errquit(
     $   'movecs_gather: ga_destroy failed', g_tmp)
      if (.not. ma_pop_stack(l_fulleval)) call errquit(
     $   'movecs_gather: pop failed', l_fulleval)
      if (.not. ma_pop_stack(l_fullocc)) call errquit(
     $   'movecs_gather: pop failed', l_fullocc)
c
      movecs_gather = .true.
c
      end
      logical function movecs_write(rtdb, basis, filename, 
     $     scftype, title,
     $     nbf, nsets, nmo, occ, ldocc, evals, ldevals, g_vecs)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "inp.fh"
#include "msgids.fh"
#include "cscfps.fh"
#include "util.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
c
c     Temporary routine
c
      integer rtdb              ! [input] RTDB handle (-1 if not accessible)
      integer basis             ! [input] Basis handle(-1 if not accessible)
      character*(*) filename    ! [input] File to write to
      character*(*) scftype     ! [input] Type of SCF (dft, scf, mcscf)
      character*(*) title       ! [input] Title of job that created vectors
      integer nbf               ! [input] No. of functions in basis
      integer nsets             ! [input] No. of sets of vectors
      integer nmo(nsets)        ! [input] No. of vectors in each set
      integer ldocc             ! [input]
      integer ldevals           ! [input]
      double precision occ(ldocc, nsets) ! [input] Occupation numbers 
      double precision evals(ldevals, nsets) ! [input] Eigenvalues
      integer g_vecs(nsets)     ! Global array with eigen-vectors
c
      integer unitno            ! Unit no. for writing
      parameter (unitno = 67)   ! These need to be managed !!!
      integer lentit
      integer lenbas
      integer l_vecs, k_vecs
      integer ok, iset, i, j
      integer geom, ma_type, nelem
      character*26 date
      character*32 geomsum, basissum, key
      character*20 scftype20    ! DECL MUST match movecs_converged routine
      character*128 basis_name, trans_name
c
      if (oscfps) call pstat_on(ps_vecio)
      call ga_sync()
c
c     Generate info to go into convergence and restart line
c
      if (basis.ne.-1 .and. rtdb.ne.-1) then
         if (.not. bas_checksum(basis, basissum))
     $        call errquit('movecs_write: bad basis',0)
         if (.not. bas_geom(basis, geom))
     $        call errquit('movecs_write: bad basis',0)
         if (.not. geom_checksum(geom, geomsum))
     $        call errquit('movecs_write: bad geometry',0)
         if (.not. bas_name(basis, basis_name, trans_name))
     $        call errquit('movecs_write: bad basis',0)
         if (trans_name .ne. ' ') basis_name = trans_name
c         
         key = ' '
         write(key,'(a,'':converged'')') scftype(1:inp_strlen(scftype))
         if (.not. rtdb_get_info(rtdb, key, ma_type, nelem, date))
     $        date = ' '
      else
         basissum = ' '
         geomsum  = ' '
         date     = ' '
      endif
      scftype20 = scftype
c
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
c     Information about convergence
c
         write(unitno, err=1001) basissum, geomsum, scftype20, date
c
c     Check that read routines are both consistent with this
c
         write(unitno, err=1001) scftype20
         lentit = max(1,inp_strlen(title)) ! 0 length record confuses f2c
         write(unitno, err=1001) lentit
         write(unitno, err=1001) title(1:lentit)
         lenbas = max(1,inp_strlen(basis_name))
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
               call swrite(unitno, dbl_mb(k_vecs), nbf)
*               write(unitno, err=1001) (dbl_mb(k_vecs+j), j=0,nbf-1)
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
     $     util_print('vectors i/o', print_high)) then
         write(6,22) filename(1:inp_strlen(filename))
 22      format(/' Wrote molecular orbitals to ',a/)
         call util_flush(6)
      endif
      call ga_sync()
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
      subroutine movecs_swap(rtdb, module, scftype,g_vecs,occ,eval)
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "inp.fh"
c
      character*(*) scftype
      integer rtdb, g_vecs(*)
      double precision occ(*), eval(*)
      character*(*) module
c
c     In the database may be an entry '<module>:swap' (with
c     alpha/beta appended for UHF) which is a list of pairs
c     of vectors to swap.  Apply these to the input vectors, 
c     occupation numbers and eigenvalues.  In the case of
c     UHF g_vecs is assumed to be an array with 2 handles.
c     When finished delete the databse entries so that
c     a restart does not re-swap the vectors.
c
      double precision tmp
      integer type, nelem, handle, index
      integer l_veci, k_veci, l_vecj, k_vecj, dim1, dim2
      integer pair, i, j, ioff, joff, set, nset, mlen
      logical oprint
      character*40 name
      character*8 text
c
      oprint = util_print('mo guess', print_default)
      mlen = inp_strlen(module)
c
      if (scftype .eq. 'UHF') then
         nset = 2
      else
         nset = 1
      endif
c
      do set = 1, nset
         name = ' '
         if (scftype .eq. 'UHF') then
            if (set.eq.1) then
               write(name,'(a,a)') module(1:mlen),':swap alpha'
               text = 'alpha'
            else
               write(name,'(a,a)') module(1:mlen),':swap beta'
               text = 'beta'
            endif
         else
            write(name,'(a,a)') module(1:mlen),':swap'
            text = scftype
         endif
         if (.not. rtdb_ma_get(rtdb, name, type, nelem, 
     $        handle)) goto 1000
         if (.not. ma_get_index(handle,index)) call errquit
     $        ('movecs_swap: get index failed', handle)
         if (mod(nelem,2).ne.0 .or. type.ne.mt_int) call errquit
     $        ('movecs_swap: invalid list of pairs',nelem)
c
         call ga_sync()
         call ga_inquire(g_vecs, type, dim1, dim2)
         if (ga_nodeid() .eq. 0) then
            if (oprint) write(6,*)
            if (.not. ma_push_get(mt_dbl, dim1, 'swap', l_veci, k_veci))
     $           call errquit('movecs_swap: no scratch space', dim1)
            if (.not. ma_push_get(mt_dbl, dim1, 'swap', l_vecj, k_vecj))
     $           call errquit('movecs_swap: no scratch space', dim1)
            do pair = 1, nelem, 2
               i = int_mb(index+pair-1)
               j = int_mb(index+pair  )
               if (oprint) write(6,1) text(1:inp_strlen(text)), i, j
 1             format(' Swapping ',a,' orbitals ', 2i5)
               if (i.lt.0 .or. i.gt.dim2) call errquit
     $              ('movecs_swap: invalid vector ', i)
               if (j.lt.0 .or. j.gt.dim2) call errquit
     $              ('movecs_swap: invalid vector ', j)
c     
               call ga_get(g_vecs(set), 1, dim1, i, i, dbl_mb(k_veci),1)
               call ga_get(g_vecs(set), 1, dim1, j, j, dbl_mb(k_vecj),1)
               call ga_put(g_vecs(set), 1, dim1, i, i, dbl_mb(k_vecj),1)
               call ga_put(g_vecs(set), 1, dim1, j, j, dbl_mb(k_veci),1)
c
            enddo
            if (oprint) call util_flush(6)
            if (.not. ma_pop_stack(l_vecj)) call errquit('ms:ma?', 0)
            if (.not. ma_pop_stack(l_veci)) call errquit('ms:ma?', 0)
         endif
c
         do pair = 1, nelem, 2
           i = int_mb(index+pair-1)
           j = int_mb(index+pair  )
           ioff = i + (set-1)*dim1
           joff = j + (set-1)*dim1
           tmp = occ(ioff)
           occ(ioff) = occ(joff)
           occ(joff) = tmp
           tmp = eval(ioff)
           eval(ioff) = eval(joff)
           eval(joff) = tmp
         enddo
c
         call ga_sync
c
         if (.not. rtdb_delete(rtdb, name)) call errquit
     $        ('movecs_swap: delete of swap entry failed', 0)
c
         if (.not. ma_free_heap(handle)) call errquit('ms:ma?', 0)
 1000    continue
      enddo
c
      end
      logical function file_write_ga(fname, g_a)
      implicit none
c
c     Generic routine for sequential write of global array to file
c
      character*(*) fname            ! [input] Name of file to write to
      integer g_a                    ! [input] Array to write
c
      integer nrow, ncol, type
      logical file_write_ga_patch
      external file_write_ga_patch
c
      call ga_inquire(g_a, type, nrow, ncol)
c
      file_write_ga = file_write_ga_patch(fname, g_a, 1, nrow, 1, ncol)
c
      end
      logical function file_read_ga(fname, g_a)
      implicit none
c
c     Generic routine for sequential read of global array from file
c
      character*(*) fname            ! [input] Name of file to read from
      integer g_a                    ! [input] Array to read
c
      integer nrow, ncol, type
      logical file_read_ga_patch
      external file_read_ga_patch
c
      call ga_inquire(g_a, type, nrow, ncol)
c
      file_read_ga = file_read_ga_patch(fname, g_a, 1, nrow, 1, ncol)
c
      end
      logical function movecs_converged(rtdb, basis, scftype, filename)
      implicit none
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "bas.fh"
#include "geom.fh"
#include "tcgmsg.fh"
      integer rtdb              ! [input]
      integer basis             ! [input]
      character*(*) scftype     ! [input] Current SCF type
      character*(*) filename    ! [input] MOvecs filename
c
c     Return true/false if the MO vectors in the file correspond
c     to converged vectors of the current SCF type with the
c     current set of options.
c
c     1) The movecs file contains basis checksum, geometry checksum,
c     .  scftype, date on convergence flag in the database at the
c     .  time the vectors were written.
c
c     2) The database contains a convergence flag that is deleted
c     .  whenever any input is changed or options automatically 
c     .  changed.
c
c     The calculation is converged if all of the following are true
c     
c     a) convergence flag is present in the database and the date 
c     .  matches the date in the movecs file
c     b) scftypes match
c     c) basis and geometry checksums match
c
      character*20 filescftype  ! MUST MATCH WRITE ROUTINE DECL
      character*32 filebasissum, filegeomsum, basissum, geomsum
      logical oconverged
      character*40 key
      character*26 date, filedate
      integer geom, ma_type, nelem
      integer unitno            ! Unit no. for reading
      parameter (unitno = 67)   ! These need to be managed !!!
      logical status, oldmode
c
      status = .false.
c
      oldmode =  rtdb_parallel(.false.)
      if (ga_nodeid() .eq. 0) then
         open(unitno, status='old', form='unformatted',
     $        file=filename, err=1000)
         read(unitno, end=1001, err=1001)
     $        filebasissum, filegeomsum, filescftype, filedate
         close(unitno, err=1002)
c     
c     a)
c     
         key = ' '
         write(key,'(a,'':converged'')') scftype(1:inp_strlen(scftype))
         if (.not. rtdb_get(rtdb, key, mt_log, 1, oconverged)) goto 10
         if (.not. oconverged) goto 10
         if (.not. rtdb_get_info(rtdb, key, ma_type, nelem, date))
     $        call errquit('movecs_converged: rtdb corrupt?',0)
         if (date .ne. filedate) goto 10
c     
c     b) 
c     
         if (.not. inp_compare(.false., scftype, filescftype)) goto 10
c     
c     c) 
c     
         if (.not. bas_geom(basis, geom)) 
     $        call errquit('movecs_converged: bad basis',0)
         if (.not. bas_checksum(basis, basissum))
     $        call errquit('movecs_converged: bad basis',0)
         if (.not. geom_checksum(geom, geomsum))
     $        call errquit('movecs_converged: bad geometry',0)
         if (basissum.ne.filebasissum .or. geomsum.ne.filegeomsum) 
     $        goto 10
c     
c     Cool ... it's converged!
c
         status = .true.
      endif
c
 10   oldmode =  rtdb_parallel(oldmode)
      call ga_brdcst(1311, status, mitob(1), 0)
      movecs_converged = status
      return
c
 1000 call errquit('movecs_converged: error opening file',0)
 1001 call errquit('movecs_converged: error reading file',0)
 1002 call errquit('movecs_converged: error closing file',0)
c
      end
      logical function file_read_ga_info(fname,
     $     title, nrow, ncol)
      implicit none
#include "global.fh"
#include "tcgmsg.fh"
#include "inp.fh"
c
      character*(*) fname
      character*(*) title
      integer nrow, ncol
c
c     Read the header info from the dumb file_ga file
c
      integer unitno                 ! Unit no. for writing
      parameter (unitno = 86)
      integer len1, ok
c
      ok = 1
      if (ga_nodeid() .eq. 0) then
         open(unitno, status='old', form='unformatted',
     $        file=fname, err=1000)
         read(unitno, err=1001) len1
         read(unitno, err=1001) title(1:len1)
         read(unitno, err=1001) nrow, ncol
         close(unitno,err=1002)
         goto 10
c
 1000    write(6,*) ' file_read_ga_info: failed to open ',
     $        fname(1:inp_strlen(fname))
         call util_flush(6)
         ok = 0
         goto 10
c     
 1001    write(6,*) ' file_read_ga_info: failing reading from ',
     $        fname(1:inp_strlen(fname))
         call util_flush(6)
         ok = 0
         close(unitno,err=1002)
         goto 10
c     
 1002    write(6,*) ' file_read_ga_info: failed to close',
     $        fname(1:inp_strlen(fname))
         call util_flush(6)
         ok = 0
         goto 10
      end if
c
 10   call ga_sync()
      call ga_brdcst(1, ok, mitob(1), 0) ! Propagate status
      call ga_brdcst(2, nrow, mitob(1), 0)
      call ga_brdcst(3, ncol, mitob(1), 0)
      call util_char_ga_brdcst(4, title, 0)
c
      file_read_ga_info = ok .eq. 1
c
      end
      logical function file_write_ga_patch(fname, g_a, 
     $     ilo, ihi, jlo, jhi)
      implicit none
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "mafdecls.fh"
#include "inp.fh"
#include "msgids.fh"
#include "util.fh"
c
c     Generic routine for sequential write 
c     of global array to file
c
      character*(*) fname            ! [input] Name of file to write to
      integer g_a                    ! [input] Array to write
      integer ilo, ihi, jlo, jhi     ! [input] Patch to write
c
      integer unitno                 ! Unit no. for writing
      parameter (unitno = 86)
      character*255 gtitle
      integer l_tmp, k_tmp
      integer len1, nrow, ncol
      integer ok, j
c
      call ga_sync
c
      ok = 0
      nrow = ihi-ilo+1
      ncol = jhi-jlo+1
      call ga_inquire_name(g_a, gtitle)
      len1 = inp_strlen(gtitle)
      if (ga_nodeid() .eq. 0) then
        open(unitno, status='unknown', form='unformatted',
     $       file=fname, err=1000)
        write(unitno, err=1001) len1
        write(unitno, err=1001) gtitle(1:len1)
        write(unitno, err=1001) nrow, ncol
        if (.not. ma_push_get(mt_dbl,nrow,'movecs_write',l_tmp,k_tmp))
     $    call errquit('matrix_fwrite: ma failed', nrow)
        do j=jlo,jhi
          call ga_get(g_a, ilo, ihi, j, j, dbl_mb(k_tmp), 1)
          call swrite(unitno, dbl_mb(k_tmp), nrow)
****          write(unitno, err=1001) (dbl_mb(k_tmp+i), i=0,nrow-1)
        enddo
        if (.not. ma_pop_stack(l_tmp))
     $    call errquit('file_write_ga: ma pop failed', l_tmp)
        close(unitno,err=1002)
        ok = 1
      endif
c
 10   call ga_sync
      call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, mitob(1), 0)        ! Propagate status
c
      file_write_ga_patch = ok .eq. 1
      if (ga_nodeid() .eq. 0 .and.
     $     util_print('vectors i/o', print_high)) then
         write(6,22) gtitle(1:len1), fname(1:inp_strlen(fname))
 22      format(/' Wrote ',a,' to ',a/)
         call util_flush(6)
      endif
      return
c
 1000 write(6,*) ' file_write_ga: failed to open ',
     $     fname(1:inp_strlen(fname))
      call util_flush(6)
      ok = 0
      goto 10
c
 1001 write(6,*) ' file_write_ga: failing writing to ',
     $     fname(1:inp_strlen(fname))
      call util_flush(6)
      ok = 0
      close(unitno,err=1002)
      goto 10
c
 1002 write(6,*) ' file_write_ga: failed to close',
     $     fname(1:inp_strlen(fname))
      call util_flush(6)
      ok = 0
      goto 10
c
      end
      logical function file_read_ga_patch
     $     (fname, g_a, ilo, ihi, jlo, jhi)
      implicit none
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "mafdecls.fh"
#include "inp.fh"
#include "msgids.fh"
#include "util.fh"
c
c     Generic routine for sequential read 
c     of file to global array
c
      character*(*) fname            ! [input] Name of file to read from
      integer g_a                    ! [output] Array
      integer ilo, ihi, jlo, jhi     ! [input]
c
      integer unitno                 ! Unit no. for writing
      parameter (unitno = 86)
      character*255 gtitle
      integer l_tmp, k_tmp
      integer len1, nrow, ncol
      integer grow, gcol
      integer ok, j
c
      call ga_sync()
c
      ok = 0
      grow = ihi-ilo+1
      gcol = jhi-jlo+1
      call ga_inquire_name(g_a, gtitle)
      len1 = inp_strlen(gtitle)
      if (ga_nodeid() .eq. 0) then
        open(unitno, status='old', form='unformatted',
     $       file=fname, err=1000)
        read(unitno, err=1001) len1
        read(unitno, err=1001) gtitle(1:len1)
        read(unitno, err=1001) nrow, ncol
        if ((nrow.ne.grow).and.(ncol.ne.gcol)) goto 2000
        if (.not. ma_push_get(mt_dbl,nrow,'movecs_write',l_tmp,k_tmp))
     $    call errquit('matrix_fwrite: ma failed', nrow)
        do j=jlo,jhi
          call sread(unitno, dbl_mb(k_tmp), nrow)
***          read(unitno, err=1001) (dbl_mb(k_tmp+i), i=0,nrow-1)
          call ga_put(g_a, ilo, ihi, j, j, dbl_mb(k_tmp), 1)
        enddo
        if (.not. ma_pop_stack(l_tmp))
     $    call errquit('file_read_ga: ma pop failed', l_tmp)
        close(unitno,err=1002)
        ok = 1
      endif
c
 10   call ga_sync()
      call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, mitob(1), 0)        ! Propagate status
c
      file_read_ga_patch = ok .eq. 1
      if (ga_nodeid() .eq. 0 .and.
     $     util_print('vectors i/o', print_high)) then
         write(6,22) gtitle(1:len1), fname(1:inp_strlen(fname))
 22      format(/' Read ',a,' from ',a/)
         call util_flush(6)
      endif
      return
c
 1000 write(6,*) ' file_read_ga: failed to open ',
     $     fname(1:inp_strlen(fname))
      call util_flush(6)
      ok = 0
      goto 10
c
 1001 write(6,*) ' file_read_ga: failing writing to ',
     $     fname(1:inp_strlen(fname))
      call util_flush(6)
      ok = 0
      close(unitno,err=1002)
      goto 10
c
 1002 write(6,*) ' file_read_ga: failed to close',
     $     fname(1:inp_strlen(fname))
      call util_flush(6)
      ok = 0
      goto 10
c
 2000 write(6,*) ' file_read_ga: GA and file ',
     $           'contents mismatch dimensions',
     $           nrow, ncol, grow, gcol
      call util_flush(6)
      ok = 0
      goto 10
c
      end
