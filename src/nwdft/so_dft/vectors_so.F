      logical function movecs_read_header_so
     $    (filename, title, basis_name, 
     $     scftype, nbf)
C$Id$
      implicit none
#include "errquit.fh"
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
c
      integer unitno            ! Unit no. for reading
      parameter (unitno = 67)   ! These need to be managed !!!
      integer lentit
      integer lenbas
      integer ok
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
     $        ('movecs_read_header_so: title too short', lentit,
     &       DISK_ERR)
         title = ' '
         read(unitno, err=1001, end=2001) title(1:lentit)
         read(unitno, err=1001, end=2001) lenbas
         if (len(basis_name) .lt. lenbas) call errquit
     $        ('movecs_read_header_so: basis_name too short', lenbas,
     &       DISK_ERR)
         basis_name = ' '
         read(unitno, err=1001, end=2001) basis_name(1:lenbas)
         read(unitno, err=1001, end=2001) nbf
         close(unitno, err=1002)
         ok = 1
      endif
c
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, mitob(1), 0) ! Propagate status
      movecs_read_header_so = ok .eq. 1
      if (.not. movecs_read_header_so) return
c
      call util_char_ga_brdcst(Msg_Vec_Ttl, title, 0)
      call util_char_ga_brdcst(Msg_Vec_Nam, basis_name, 0)
      call ga_brdcst(Msg_Vec_NBF+MSGINT, nbf, mitob(1), 0)
c
      if (oscfps) call pstat_off(ps_vecio)
c
      return
c
 1000 write(6,*) ' movecs_read_header_so: failed to open ',
     $     filename(1:inp_strlen(filename))
      call util_flush(6)
      ok = 0
      goto 10
c
 1001 write(6,*) ' movecs_read_header_so: failing reading from ',
     $     filename(1:inp_strlen(filename))
      call util_flush(6)
      ok = 0
      close(unitno, err=1002)
      goto 10
c
 2001 write(6,*) ' movecs_read_header_so:eof: failing reading from '
     &    , filename(1:inp_strlen(filename))
      call util_flush(6)
      ok = 0
      close(unitno, err=1002)
      goto 10
c
 1002 write(6,*) ' movecs_read_header_so: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(6)
      ok = 0
      goto 10
c
      end
      logical function movecs_read_so
     $     (filename, occ, evals, g_vecs)
      implicit none
#include "errquit.fh"
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
      double precision occ(*)   ! Must be at least nbf long (not nmo)
      double precision evals(*) ! Must be at least nbf long (not nmo)
      integer g_vecs(2)
c
      integer nbf               ! No. of functions in basis
      integer ok, i, j
      integer l_vecs, k_vecs
      integer iset
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
         read(unitno, err=1001, end=1001) nbf
c
         if (.not. ma_push_get(mt_dbl,nbf,'movecs_read_so',
     $        l_vecs,k_vecs))
     $        call errquit('movecs_read_so: ma failed', nbf, MA_ERR)
c
c     Skip over unwanted sets
c
         read(unitno, err=1001, end=1001) (occ(j),j=1,nbf)
         read(unitno, err=1001, end=1001) (evals(j),j=1,nbf)
         do iset = 1, 2
          do i = 1, nbf
            call sread(unitno, dbl_mb(k_vecs), nbf)
            call ga_put(g_vecs(iset), 1, nbf, i, i, dbl_mb(k_vecs), 1)
          enddo
         enddo
 9       close(unitno,err=1002)
         ok = 1
      endif
c
 10   continue
      if (l_vecs .ne. -1) then
	 if (.not. ma_pop_stack(l_vecs)) call errquit
     $      ('movecs_read_so: pop failed', l_vecs, MA_ERR)
      endif
      call ga_sync()
      call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, mitob(1), 0) ! Propagate status
      if (ok .eq. 1) then
	 call ga_brdcst(Msg_Vec_NBF+MSGINT, nbf, mitob(1), 0)
	 call ga_brdcst(Msg_Vec_EVal+MSGDBL, evals, mdtob(nbf), 0)
	 call ga_brdcst(Msg_Vec_Occ+MSGDBL, occ,   mdtob(nbf), 0)
      endif
c
      movecs_read_so = ok .eq. 1
      if (ga_nodeid() .eq. 0 .and. movecs_read_so .and.
     $     util_print('vectors i/o', print_default)) then
         write(6,22) filename(1:inp_strlen(filename))
 22      format(/' Read molecular orbitals from ',a/)
         call util_flush(6)
      endif
      if (oscfps) call pstat_off(ps_vecio)
      return
c
 1000 write(6,*) ' movecs_read_so: failed to open ',
     $     filename(1:inp_strlen(filename))
      call util_flush(6)
      ok = 0
      goto 10
c
 1001 write(6,*) ' movecs_read_so: failing reading from ',
     $     filename(1:inp_strlen(filename))
      call util_flush(6)
      ok = 0
      close(unitno,err=1002)
      goto 10
c
 1002 write(6,*) ' movecs_read_so: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(6)
      ok = 0
      goto 10
c
      end
      logical function movecs_gather_so(filename, iset, nwant, want,
     $   occ, evals, g_vecs)
      implicit none
#include "errquit.fh"
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
      integer nbf                   ! no. of functions in basis
      integer ldnmo                 ! size of nmo
      parameter (ldnmo = 2)
      integer nmo(ldnmo)            ! no. of vectors in each set
      integer i, l_fullocc, k_fullocc, l_fulleval, k_fulleval, g_tmp(2)
      character*255 title, basis
      integer minchunk              ! for data distirbution in ga
      parameter (minchunk = -1)     ! let ga determine even distrib.
      character*20 scftype
c
      logical movecs_read_header_so, movecs_read_so
      external movecs_read_header_so, movecs_read_so
c     
      if ( .not. movecs_read_header_so(filename, title, basis,
     $     scftype, nbf) ) call errquit(
     $   'movecs_gather_so: failed to read movecs file', 1, DISK_ERR)
c
c     allocate memory to be used while reading info
C     movecs_read_so expects fullocc and fulleval arrays to be nbf,
C     not nmo(iset)!
c
      if (.not. ma_push_get(mt_dbl,nbf,'movecs_gather_so',
     $   l_fullocc, k_fullocc)) call errquit(
     $   'movecs_gather_so: ma failed', nmo(iset), MA_ERR)
c
      if (.not. ma_push_get(mt_dbl,nbf,'movecs_gather_so',
     $   l_fulleval, k_fulleval)) call errquit(
     $   'movecs_gather_so: ma failed', nmo(iset), MA_ERR)
c
      if ( .not. ga_create( mt_dbl, nbf, nmo(iset), 'full vectors',
     $   minchunk, minchunk, g_tmp(iset)) ) call errquit(
     $   'movecs_gather_so: ga_create failed', nbf*nmo(iset) , GA_ERR)
c
c     read the full data
c
       if ( .not. movecs_read_so(filename, dbl_mb(k_fullocc),
     $   dbl_mb(k_fulleval), g_tmp) ) call errquit(
     $   'movecs_gather_so: failed reading movecs file', 2, DISK_ERR)
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
     $   'movecs_gather_so: ga_destroy failed', g_tmp, GA_ERR)
      if (.not. ma_pop_stack(l_fulleval)) call errquit(
     $   'movecs_gather_so: pop failed', l_fulleval, MA_ERR)
      if (.not. ma_pop_stack(l_fullocc)) call errquit(
     $   'movecs_gather_so: pop failed', l_fullocc, MA_ERR)
c
      movecs_gather_so = .true.
c
      end
      logical function movecs_write_so(rtdb, basis, filename, 
     $     scftype, title, nbf, occ, evals, g_vecs)
      implicit none
#include "errquit.fh"
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
#include "bq.fh"
#include "stdio.fh"
c
c     Temporary routine
c
      integer rtdb              ! [input] RTDB handle (-1 if not accessible)
      integer basis             ! [input] Basis handle(-1 if not accessible)
      character*(*) filename    ! [input] File to write to
      character*(*) scftype     ! [input] Type of SCF (dft, scf, mcscf)
      character*(*) title       ! [input] Title of job that created vectors
      integer nbf               ! [input] No. of functions in basis
      double precision occ(nbf) ! [input] Occupation numbers 
      double precision evals(nbf) ! [input] Eigenvalues
      integer g_vecs(2)     ! Global array with eigen-vectors
c
      integer unitno            ! Unit no. for writing
      parameter (unitno = 67)   ! These need to be managed !!!
      integer bq_handle         ! Handle of the BQ fragment
      integer lentit
      integer lenbas
      integer l_vecs, k_vecs
      integer ok, iset, i, j
      integer geom, ma_type, nelem
      character*26 date
      character*32 geomsum, basissum, bqsum, key
      character*20 scftype20    ! DECL MUST match movecs_converged routine
      character*128 basis_name, trans_name
c
      logical converged
c
      logical util_io_unit
      external util_io_unit
c
      if (oscfps) call pstat_on(ps_vecio)
      call ga_sync()
c
c     Generate info to go into convergence and restart line
c
      if (basis.ne.-1 .and. rtdb.ne.-1) then
         if (.not. bas_checksum(basis, basissum))
     $        call errquit('movecs_write_so: bad basis',0, BASIS_ERR)
         if (.not. bas_geom(basis, geom))
     $        call errquit('movecs_write_so: bad basis',0, BASIS_ERR)
         if (.not. geom_checksum(geom, geomsum))
     $        call errquit('movecs_write_so: bad geometry',0, GEOM_ERR)
         if (.not. bas_name(basis, basis_name, trans_name))
     $        call errquit('movecs_write_so: bad basis',0, BASIS_ERR)
         if (trans_name .ne. ' ') basis_name = trans_name
         if (.not. bq_get_active(bq_handle)) then
              bqsum = ' '
         else
           if (.not. bq_namespace(bq_handle,bqsum)) then
              bqsum = ' '
           endif
         endif
c         
         key = ' '
         write(key,'(a,'':converged'')') scftype(1:inp_strlen(scftype))
         if (.not. rtdb_get_info(rtdb, key, ma_type, nelem, date))
     $        date = ' '
      else
         basissum = ' '
         geomsum  = ' '
         bqsum    = ' '
         date     = ' '
      endif
      scftype20 = scftype
c
c     Check convergence flag
c
      converged = .false.
      if(.not. rtdb_get(rtdb, 'sodft:converged',MT_LOG,1,converged))
     &          converged = .false.
c
      ok = 0
      if (ga_nodeid() .eq. 0) then
         open(unitno, status='unknown', form='unformatted',
     $        file=filename, err=1000)
c
c        Information about convergence
c
         write(unitno, err=1001) basissum, geomsum, bqsum, scftype20,
     $                           date
c
c        Check that read routines are both consistent with this
c
         write(unitno, err=1001) scftype20
         lentit = max(1,inp_strlen(title)) ! 0 length record confuses f2c
         write(unitno, err=1001) lentit
         write(unitno, err=1001) title(1:lentit)
         lenbas = max(1,inp_strlen(basis_name))
         write(unitno, err=1001) lenbas
         write(unitno, err=1001) basis_name(1:lenbas)
         write(unitno, err=1001) nbf
         if (.not. ma_push_get(mt_dbl,nbf,'movecs_write_so',
     $        l_vecs,k_vecs))
     $        call errquit('movecs_write_so: ma failed', nbf, MA_ERR)
         write(unitno, err=1001) (occ(j),j=1,nbf)
         write(unitno, err=1001) (evals(j),j=1,nbf)
c
         do iset = 1, 2
            do i = 1, nbf
               call ga_get(g_vecs(iset), 1, 
     $              nbf, i, i, dbl_mb(k_vecs),1)
               call swrite(unitno, dbl_mb(k_vecs), nbf)
            enddo
         enddo
         if (.not. ma_pop_stack(l_vecs))
     $        call errquit('movecs_write_so: ma pop failed', l_vecs,
     &       MA_ERR)
         close(unitno,err=1002)
         ok = 1
      endif ! ga_nodeid() = 0
c
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, mitob(1), 0) ! Propagate status
c
      movecs_write_so = ok .eq. 1
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
 1000 write(6,*) ' movecs_write_so: failed to open ',
     $     filename(1:inp_strlen(filename))
      call util_flush(6)
      ok = 0
      goto 10
c
 1001 write(6,*) ' movecs_write_so: failing writing to ',
     $     filename(1:inp_strlen(filename))
      call util_flush(6)
      ok = 0
      close(unitno,err=1002)
      goto 10
c
 1002 write(6,*) ' movecs_write_so: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(6)
      ok = 0
      goto 10
c
      end
      subroutine movecs_swap_so(rtdb, module, scftype,g_vecs,occ,eval)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "inp.fh"
c
      character*(*) scftype
      integer rtdb, g_vecs(2)
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
      integer pair, i, j, ioff, joff, mlen
      logical oprint
      character*40 name
      character*8 text
c
      oprint = util_print('mo guess', print_default)
      mlen = inp_strlen(module)
c
      name = ' '
      write(name,'(a,a)') 'sodft',':swap'
      text = 'sodft'
      if (.not. rtdb_ma_get(rtdb, name, type, nelem, 
     $     handle)) goto 1000
      if (.not. ma_get_index(handle,index)) call errquit
     $     ('movecs_swap_so: get index failed', handle, MA_ERR)
      if (mod(nelem,2).ne.0 .or. type.ne.mt_int) call errquit
     $     ('movecs_swap_so: invalid list of pairs',nelem, INPUT_ERR)
c
      call ga_sync()
      call ga_inquire(g_vecs(1), type, dim1, dim2)
      if (ga_nodeid() .eq. 0) then
         if (oprint) write(6,*)
         if (.not. ma_push_get(mt_dbl, dim1, 'swap', l_veci, k_veci))
     $        call errquit('movecs_swap_so: no scratch space', dim1,
     &       MA_ERR)
         if (.not. ma_push_get(mt_dbl, dim1, 'swap', l_vecj, k_vecj))
     $        call errquit('movecs_swap_so: no scratch space', dim1,
     &       MA_ERR)
         do pair = 1, nelem, 2
            i = int_mb(index+pair-1)
            j = int_mb(index+pair  )
            if (oprint) write(6,1) text(1:inp_strlen(text)), i, j
 1          format(' Swapping ',a,' orbitals ', 2i5)
            if (i.lt.0 .or. i.gt.dim2) call errquit
     $           ('movecs_swap_so: invalid vector ', i, INPUT_ERR)
            if (j.lt.0 .or. j.gt.dim2) call errquit
     $           ('movecs_swap_so: invalid vector ', j, INPUT_ERR)
c     
            call ga_get(g_vecs(1), 1, dim1, i, i, dbl_mb(k_veci),1)
            call ga_get(g_vecs(1), 1, dim1, j, j, dbl_mb(k_vecj),1)
            call ga_put(g_vecs(1), 1, dim1, i, i, dbl_mb(k_vecj),1)
            call ga_put(g_vecs(1), 1, dim1, j, j, dbl_mb(k_veci),1)
c     
            call ga_get(g_vecs(2), 1, dim1, i, i, dbl_mb(k_veci),1)
            call ga_get(g_vecs(2), 1, dim1, j, j, dbl_mb(k_vecj),1)
            call ga_put(g_vecs(2), 1, dim1, i, i, dbl_mb(k_vecj),1)
            call ga_put(g_vecs(2), 1, dim1, j, j, dbl_mb(k_veci),1)
c     
         enddo  ! pair
c
         if (oprint) call util_flush(6)
         if (.not. ma_pop_stack(l_vecj)) call errquit('ms:ma?', 0,
     &       MA_ERR)
         if (.not. ma_pop_stack(l_veci)) call errquit('ms:ma?', 0,
     &       MA_ERR)
      endif  ! ga_nodeid()
c     
      do pair = 1, nelem, 2
         i = int_mb(index+pair-1)
         j = int_mb(index+pair  )
         ioff = i 
         joff = j 
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
     $     ('movecs_swap_so: delete of swap entry failed', 0, RTDB_ERR)
c     
      if (.not. ma_free_heap(handle)) call errquit('ms:ma?', 0, MA_ERR)
 1000 continue
c     
      end
      logical function file_write_ga_so(fname, g_a)
      implicit none
c
c     Generic routine for sequential write of global array to file
c
      character*(*) fname            ! [input] Name of file to write to
      integer g_a                    ! [input] Array to write
c
      integer nrow, ncol, type
      logical file_write_ga_patch_so
      external file_write_ga_patch_so
c
      call ga_inquire(g_a, type, nrow, ncol)
c
      file_write_ga_so = 
     $     file_write_ga_patch_so(fname, g_a, 1, nrow, 1, ncol)
c
      end
      logical function file_read_ga_so(fname, g_a)
      implicit none
c
c     Generic routine for sequential read of global array from file
c
      character*(*) fname            ! [input] Name of file to read from
      integer g_a                    ! [input] Array to read
c
      integer nrow, ncol, type
      logical file_read_ga_patch_so
      external file_read_ga_patch_so
c
      call ga_inquire(g_a, type, nrow, ncol)
c
      file_read_ga_so = 
     $     file_read_ga_patch_so(fname, g_a, 1, nrow, 1, ncol)
c
      end
      logical function movecs_converged_so
     $     (rtdb, basis, scftype, filename)
      implicit none
#include "errquit.fh"
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
     $        call errquit('movecs_converged_so: rtdb corrupt?',0,
     &       RTDB_ERR)
         if (date .ne. filedate) goto 10
c     
c     b) 
c     
         if (.not. inp_compare(.false., scftype, filescftype)) goto 10
c     
c     c) 
c     
         if (.not. bas_geom(basis, geom)) 
     $        call errquit('movecs_converged_so: bad basis',0,
     &       BASIS_ERR)
         if (.not. bas_checksum(basis, basissum))
     $        call errquit('movecs_converged_so: bad basis',0,
     &       BASIS_ERR)
         if (.not. geom_checksum(geom, geomsum))
     $        call errquit('movecs_converged_so: bad geometry',0,
     &       GEOM_ERR)
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
      movecs_converged_so = status
      return
c
 1000 call errquit('movecs_converged_so: error opening file',0,
     &       DISK_ERR)
 1001 call errquit('movecs_converged_so: error reading file',0,
     &       DISK_ERR)
 1002 call errquit('movecs_converged_so: error closing file',0,
     &       DISK_ERR)
c
      end
      
      logical function file_write_ga_patch_so(fname, g_a, 
     $     ilo, ihi, jlo, jhi)
      implicit none
#include "errquit.fh"
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
        if (.not. ma_push_get
     $       (mt_dbl,nrow,'movecs_write_so',l_tmp,k_tmp))
     $    call errquit('matrix_fwrite: ma failed', nrow,
     &       MA_ERR)
        do j=jlo,jhi
          call ga_get(g_a, ilo, ihi, j, j, dbl_mb(k_tmp), 1)
          call swrite(unitno, dbl_mb(k_tmp), nrow)
****          write(unitno, err=1001) (dbl_mb(k_tmp+i), i=0,nrow-1)
        enddo
        if (.not. ma_pop_stack(l_tmp))
     $    call errquit('file_write_ga_so: ma pop failed', l_tmp,
     &       MA_ERR)
        close(unitno,err=1002)
        ok = 1
      endif
c
 10   call ga_sync
      call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, mitob(1), 0)        ! Propagate status
c
      file_write_ga_patch_so = ok .eq. 1
      if (ga_nodeid() .eq. 0 .and.
     $     util_print('vectors i/o', print_high)) then
         write(6,22) gtitle(1:len1), fname(1:inp_strlen(fname))
 22      format(/' Wrote ',a,' to ',a/)
         call util_flush(6)
      endif
      return
c
 1000 write(6,*) ' file_write_ga_so: failed to open ',
     $     fname(1:inp_strlen(fname))
      call util_flush(6)
      ok = 0
      goto 10
c
 1001 write(6,*) ' file_write_ga_so: failing writing to ',
     $     fname(1:inp_strlen(fname))
      call util_flush(6)
      ok = 0
      close(unitno,err=1002)
      goto 10
c
 1002 write(6,*) ' file_write_ga_so: failed to close',
     $     fname(1:inp_strlen(fname))
      call util_flush(6)
      ok = 0
      goto 10
c
      end
      logical function file_read_ga_patch_so 
     $     (fname, g_a, ilo, ihi, jlo, jhi)
      implicit none
#include "errquit.fh"
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
        if (.not. ma_push_get
     $       (mt_dbl,nrow,'movecs_write_so',l_tmp,k_tmp))
     $    call errquit('matrix_fwrite: ma failed', nrow, MA_ERR)
        do j=jlo,jhi
          call sread(unitno, dbl_mb(k_tmp), nrow)
***          read(unitno, err=1001) (dbl_mb(k_tmp+i), i=0,nrow-1)
          call ga_put(g_a, ilo, ihi, j, j, dbl_mb(k_tmp), 1)
        enddo
        if (.not. ma_pop_stack(l_tmp))
     $    call errquit('file_read_ga_so: ma pop failed', l_tmp, MA_ERR)
        close(unitno,err=1002)
        ok = 1
      endif
c
 10   call ga_sync()
      call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, mitob(1), 0)        ! Propagate status
c
      file_read_ga_patch_so = ok .eq. 1
      if (ga_nodeid() .eq. 0 .and.
     $     util_print('vectors i/o', print_high)) then
         write(6,22) gtitle(1:len1), fname(1:inp_strlen(fname))
 22      format(/' Read ',a,' from ',a/)
         call util_flush(6)
      endif
      return
c
 1000 write(6,*) ' file_read_ga_so: failed to open ',
     $     fname(1:inp_strlen(fname))
      call util_flush(6)
      ok = 0
      goto 10
c
 1001 write(6,*) ' file_read_ga_so: failing writing to ',
     $     fname(1:inp_strlen(fname))
      call util_flush(6)
      ok = 0
      close(unitno,err=1002)
      goto 10
c
 1002 write(6,*) ' file_read_ga_so: failed to close',
     $     fname(1:inp_strlen(fname))
      call util_flush(6)
      ok = 0
      goto 10
c
 2000 write(6,*) ' file_read_ga_so: GA and file ',
     $           'contents mismatch dimensions',
     $           nrow, ncol, grow, gcol
      call util_flush(6)
      ok = 0
      goto 10
c
      end
      logical function spinor_guess(filename)
C$Id$
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
c
      integer unitno            ! Unit no. for reading
      parameter (unitno = 67)   ! These need to be managed !!!
      integer ok
      character*20 scftype20
c
      ok = 0
      if (ga_nodeid() .eq. 0) then
         if(filename.ne.'atomic'.and.
     $      filename.ne.'fragment'.and.filename.ne.'hcore') then
            open(unitno, status='old', form='unformatted', 
     $           file=filename,err=1000)
            read(unitno, err=1001, end=2001) ! SKIP convergence info
            read(unitno, err=1001, end=2001) scftype20
            close(unitno, err=1002)
            if(scftype20.eq.'sodft')ok = 1
          endif
      endif
c
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, mitob(1), 0) ! Propagate status
      spinor_guess = ok .eq. 1
      return
 1000 write(6,*) ' spinor_guess: failed to open ',
     $     filename(1:inp_strlen(filename))
      call util_flush(6)
      ok = 0
      goto 10
c
 1001 write(6,*) ' spinor_guess: failing reading from ',
     $     filename(1:inp_strlen(filename))
      call util_flush(6)
      ok = 0
      close(unitno, err=1002)
      goto 10
c
 2001 write(6,*) ' spinor_guess:eof: failing reading from '
     &    , filename(1:inp_strlen(filename))
      call util_flush(6)
      ok = 0
      close(unitno, err=1002)
      goto 10
 1002 write(6,*) ' spinor_guess: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(6)
      ok = 0
      end
c
c     Write out the alpha/beta split eigenvalues
c
      subroutine print_evals_ab(nbf_ao, occ, evals)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "inp.fh"
#include "msgids.fh"
#include "util.fh"
#include "stdio.fh"
#include "cscfps.fh"
c
      integer nbf_ao                   ! [input] No. of basis functions
      double precision occ(2*nbf_ao)   ! [input] Occupation numbers 
      double precision evals(2*nbf_ao) ! [input] Eigenvalues
c
      integer i,ok
      character*255 evalfile
      integer uniteval
c
      logical util_io_unit
      external util_io_unit
c
      call util_file_name('evals',.false.,.false.,evalfile)
      if(.not.util_io_unit(80,90,uniteval))
     &  call errquit('print_evals_ab:cannot get free unit', 0,0)
c
      call ga_sync()
      ok = 0
      if (ga_nodeid() .eq. 0) then
       open(uniteval, status='unknown', form='formatted',file=evalfile)
       write(uniteval,*) nbf_ao
       write(uniteval,*) 
       do i = 1,nbf_ao
        write(uniteval,*) i, evals(2*i-1), evals(2*i)
       end do
       close(uniteval)
       ok = 1
      endif ! ga_nodeid() = 0
c
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, mitob(1), 0) ! Propagate status
c
      call ga_sync()
      if (oscfps) call pstat_off(ps_vecio)
c
      return
      end
