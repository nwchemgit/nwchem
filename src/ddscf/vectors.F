      logical function movecs_read_header(filename, title, basis_name, 
     $     scftype, nbf, nsets, nmo, ldnmo)
C$Id$
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "msgids.fh"
#include "cscfps.fh"
#include "inp.fh"
#include "stdio.fh"
#include "mafdecls.fh"
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
      integer ioserr,inntsize
      character*20 scftype20
c
      if (oscfps) call pstat_on(ps_vecio)
      ok = 0
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      if (ga_nodeid() .eq. 0) then
         open(unitno, status='old', form='unformatted', file=filename,
     $        iostat=ioserr,err=1000)
         read(unitno, err=1001, end=2001) ! SKIP convergence info
         read(unitno, err=1001, end=2001) scftype20
         scftype = scftype20
         read(unitno, err=1001, end=2001) lentit
         if (len(title) .lt. lentit) call errquit
     $        ('movecs_read_header: title too short', lentit, INPUT_ERR)
         title = ' '
         read(unitno, err=1001, end=2001) title(1:lentit)
         read(unitno, err=1001, end=2001) lenbas
         if (len(basis_name) .lt. lenbas) call errquit
     $        ('movecs_read_header: basis_name too short', lenbas,
     &       INPUT_ERR)
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
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      movecs_read_header = ok .eq. 1
      if (.not. movecs_read_header) return
c
      call util_char_ga_brdcst(Msg_Vec_Ttl, title, 0)
      call util_char_ga_brdcst(Msg_Vec_Nam, basis_name, 0)
      call util_char_ga_brdcst(Msg_Vec_Type, scftype, 0)
      call ga_brdcst(Msg_Vec_NBF+MSGINT, nbf, inntsize, 0)
      call ga_brdcst(Msg_Vec_Set+MSGINT, nsets, inntsize, 0)
      call ga_brdcst(Msg_Vec_NMO+MSGINT, nmo, inntsize*nsets, 0)
c
      if (oscfps) call pstat_off(ps_vecio)
c
      return
c
 1000 write(6,*) ' movecs_read_header: failed to open ',
     $     filename(1:inp_strlen(filename)),
     A ' IERR = ', ioserr
      call util_flush(luout)
      ok = 0
      goto 10
c
 1001 write(6,*) ' movecs_read_header: failing reading from ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno, err=1002)
      goto 10
c
 2001 write(6,*) ' movecs_read_header:eof: failing reading from '
     &    , filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno, err=1002)
      goto 10
c
 1002 write(6,*) ' movecs_read_header: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
c
      end
      logical function movecs_read(filename, iset, occ, evals, g_vecs)
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
#include "stdio.fh"
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
      integer inntsize,ddblsize
c     
      l_vecs = -1               ! An invalid MA handle
c
      if (oscfps) call pstat_on(ps_vecio)
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      ddblsize=MA_sizeof(MT_DBL,1,MT_BYTE)
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
     $        call errquit('movecs_read: ma failed', nbf, MA_ERR)
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
            call ga_put(g_vecs, 1, nbf, i, i, dbl_mb(k_vecs), nbf)
         enddo
 9       close(unitno,err=1002)
         ok = 1
      endif
c
 10   continue
      if (l_vecs .ne. -1) then
	 if (.not. ma_pop_stack(l_vecs)) call errquit
     $      ('movecs_read: pop failed', l_vecs, MA_ERR)
      endif
      call ga_sync()
      call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      if (ok .eq. 1) then
	 call ga_brdcst(Msg_Vec_NBF+MSGINT, nbf, inntsize, 0)
	 call ga_brdcst(Msg_Vec_EVal+MSGDBL, evals, ddblsize*nbf, 0)
	 call ga_brdcst(Msg_Vec_Occ+MSGDBL, occ, ddblsize*nbf, 0)
      endif
c
      movecs_read = ok .eq. 1
      if (ga_nodeid() .eq. 0 .and. movecs_read .and.
     $     util_print('vectors i/o', print_high)) then
         write(6,22) filename(1:inp_strlen(filename))
 22      format(/' Read molecular orbitals from ',a/)
         call util_flush(luout)
      endif
      if (oscfps) call pstat_off(ps_vecio)
      return
c
 1000 write(6,*) ' movecs_read: failed to open ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
c
 1001 write(6,*) ' movecs_read: failing reading from ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
c
 1002 write(6,*) ' movecs_read: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
c
      end
      logical function movecs_gather(filename, iset, nwant, want,
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
     $   'movecs_gather: failed to read movecs file', 1, INPUT_ERR)
c
c     allocate memory to be used while reading info
C     movecs_read expects fullocc and fulleval arrays to be nbf,
C     not nmo(iset)!
c
      if (.not. ma_push_get(mt_dbl,nbf,'movecs_gather',
     $   l_fullocc, k_fullocc)) call errquit(
     $   'movecs_gather: ma failed', nmo(iset), MA_ERR)
c
      if (.not. ma_push_get(mt_dbl,nbf,'movecs_gather',
     $   l_fulleval, k_fulleval)) call errquit(
     $   'movecs_gather: ma failed', nmo(iset), MA_ERR)
c
      if ( .not. ga_create( mt_dbl, nbf, nmo(iset), 'full vectors',
     $   minchunk, minchunk, g_tmp) ) call errquit(
     $   'movecs_gather: ga_create failed', nbf*nmo(iset) , GA_ERR)
c
c     read the full data
c
      if ( .not. movecs_read(filename, iset, dbl_mb(k_fullocc),
     $   dbl_mb(k_fulleval), g_tmp) ) call errquit(
     $   'movecs_gather: failed reading movecs file', 2, DISK_ERR)
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
     $   'movecs_gather: ga_destroy failed', g_tmp, GA_ERR)
      if (.not. ma_pop_stack(l_fulleval)) call errquit(
     $   'movecs_gather: pop failed', l_fulleval, MA_ERR)
      if (.not. ma_pop_stack(l_fullocc)) call errquit(
     $   'movecs_gather: pop failed', l_fullocc,
     &       MA_ERR)
c
      movecs_gather = .true.
c
      end
      logical function movecs_write(rtdb, basis, filename, 
     $     scftype, title,
     $     nbf, nsets, nmo, occ, ldocc, evals, ldevals, g_vecs)
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
#include "stdio.fh"
#include "bq.fh"
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
      integer bq_handle         ! Handle of the BQ fragment
      integer lentit
      integer lenbas
      integer l_vecs, k_vecs
      integer ok, iset, i, j
      integer geom, ma_type, nelem
      integer nmem     ! the number of real*8 words available
      integer nvec     ! the number of vectors that fit in local memory
      integer nmo_top  ! the top MO currently considered
      character*26 date
      character*32 geomsum, basissum, bqsum, key
      character*20 scftype20    ! DECL MUST match movecs_converged routine
      character*128 basis_name, trans_name
      double precision energy, enrep
      integer inntsize,kcount
c
      if (oscfps) call pstat_on(ps_vecio)
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_sync()
c
c     Generate info to go into convergence and restart line
c
      if (basis.ne.-1 .and. rtdb.ne.-1) then
         if (.not. bas_checksum(basis, basissum))
     $        call errquit('movecs_write: bad basis',0, BASIS_ERR)
         if (.not. bas_geom(basis, geom))
     $        call errquit('movecs_write: bad basis',0, BASIS_ERR)
         if (.not. geom_checksum(geom, geomsum))
     $        call errquit('movecs_write: bad geometry',0, GEOM_ERR)
         if (.not. bas_name(basis, basis_name, trans_name))
     $        call errquit('movecs_write: bad basis',0, BASIS_ERR)
         if (trans_name .ne. ' ') basis_name = trans_name
         if (.not. bq_get_active(bq_handle)) then
              bqsum = ' '
         else
           if (.not. bq_checksum(bq_handle,bqsum)) then
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
c
      ok = 0
c
c Write scf energy at the end of the movecs file. If energy
c is not in rtdb, it is given a value of zero.
c
      write(key,'(a,'':energy'')') scftype(1:inp_strlen(scftype))
      if (.not. rtdb_get(rtdb, key, mt_dbl, 1, energy)) then
        energy=0.0d0
      endif
c
c get nuclear repulsion energy
      if (.not. geom_nuc_rep_energy(geom, enrep))
     $  call errquit('unable to get nuclear rep energy',geom, GEOM_ERR)
c
      if (ga_nodeid() .eq. 0) then
         if (nsets .gt. 2) then
            write(6,*) ' movecs_write: nsets > 2 ', nsets
            goto 10
         endif
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
         write(unitno, err=1001) nsets
         write(unitno, err=1001) nbf
         write(unitno, err=1001) (nmo(i),i=1,nsets)
         nmem = ma_inquire_avail(MT_DBL)
         nvec = max(1,min(nbf,nmem/nbf))
         if (.not. ma_push_get(mt_dbl,nbf*nvec,'movecs_write',l_vecs,
     K        k_vecs))
     $        call errquit('movecs_write: ma failed', nbf, MA_ERR)
         do iset = 1, nsets
            write(unitno, err=1001) (occ(j,iset),j=1,nbf)
            write(unitno, err=1001) (evals(j,iset),j=1,nbf)
            do i = 1, nmo(iset), nvec
               nmo_top = min(nmo(iset),i+nvec-1)
               call ga_get(g_vecs(iset), 1, nbf, i, nmo_top,
     D              dbl_mb(k_vecs),nbf)
               kcount = k_vecs
               do j = i, nmo_top
                 call swrite(unitno, dbl_mb(kcount), nbf)
                 kcount = kcount + nbf
               enddo
            enddo

         enddo
         if (.not. ma_pop_stack(l_vecs))
     $        call errquit('movecs_write: ma pop failed', l_vecs,
     &       MA_ERR)
c
c Write scf energy at the end of the movecs file. If energy
c is not in rtdb, it is given a value of zero.
c
         write(unitno, err=1001) energy, enrep
c         
         close(unitno,err=1002)
         ok = 1
      endif
c
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
c
      movecs_write = ok .eq. 1
      if (ga_nodeid() .eq. 0 .and.
     $     util_print('vectors i/o', print_high)) then
         write(6,22) filename(1:inp_strlen(filename))
 22      format(/' Wrote molecular orbitals to ',a/)
         call util_flush(luout)
      endif
      call ga_sync()
      if (oscfps) call pstat_off(ps_vecio)
      return
c
 1000 write(6,*) ' movecs_write: failed to open ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
c
 1001 write(6,*) ' movecs_write: failing writing to ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
c
 1002 write(6,*) ' movecs_write: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
c
      end
      subroutine movecs_swap(rtdb, module, scftype,g_vecs,occ,eval)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "inp.fh"
#include "stdio.fh"
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
      integer pair, i, j, ioff, joff, off, set, nset, mlen
      logical oprint,dorot
      double precision mixang,cosa,sina,pi
      integer k_tj,l_tj,k_ti,l_ti
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
ccheck if we want the mix rotation
      if (rtdb_get(rtdb, 'mixang', mt_dbl, 1, mixang)) then
         dorot=.true.
         pi=acos(-1d0)
         cosa=cos(mixang*pi/180d0)
         sina=sin(mixang*pi/180d0)
c         write(6,*) ' nseeeet ',nset
      else
         dorot=.false.
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
     $        ('movecs_swap: get index failed', handle, MA_ERR)
         if (type.ne.mt_int) call errquit
     $        ('movecs_swap: invalid list of pairs',nelem, MA_ERR)
c
         call ga_sync()
         call ga_inquire(g_vecs, type, dim1, dim2)
         if (ga_nodeid() .eq. 0) then
            if (oprint) write(6,*)
            if (.not. ma_push_get(mt_dbl, dim1, 'swap', l_veci, k_veci))
     $           call errquit('movecs_swap: no scratch space', dim1,
     &       MA_ERR)
            if (.not. ma_push_get(mt_dbl, dim1, 'swap', l_vecj, k_vecj))
     $           call errquit('movecs_swap: no scratch space', dim1,
     &       MA_ERR)
            do pair = 1, nelem, 2
               i = int_mb(index+pair-1)
               j = int_mb(index+pair  )
               if (oprint) write(6,1) text(1:inp_strlen(text)), i, j
 1             format(' Swapping ',a,' orbitals ', 2i5)
               if (i.lt.0 .or. i.gt.dim2) call errquit
     $              ('movecs_swap: invalid vector ', i, INPUT_ERR)
               if (j.lt.0 .or. j.gt.dim2) call errquit
     $              ('movecs_swap: invalid vector ', j, INPUT_ERR)
c     
               call ga_get(g_vecs(set), 1, dim1, i, i, 
     D              dbl_mb(k_veci),dim1)
               call ga_get(g_vecs(set), 1, dim1, j, j, 
     D              dbl_mb(k_vecj),dim1)
c
               if(dorot) then 
                  if (.not. ma_push_get(mt_dbl, dim1, 'swap', 
     J                 l_ti, k_ti))
     $                 call errquit('moswap: pushget fail', dim1, 0)
                  if (.not. ma_push_get(mt_dbl, dim1, 'swap', 
     J                 l_tj, k_tj))
     $                 call errquit('moswap: pushget fail', dim1, 0)
                  if(set.eq.1) then
c alpha
c psi_1 =  cos*psi_1 + sin*psi_2
                 call dcopy(dim1, 0.d0, 0, dbl_mb(k_ti), 1)
                 call daxpy(dim1,cosa,dbl_mb(k_veci),1,dbl_mb(k_ti),1)
                 call daxpy(dim1,sina,dbl_mb(k_vecj),1,dbl_mb(k_ti),1)
c psi_2 = -sin*psi_1 + cos*psi_2
                 call dcopy(dim1, 0.d0, 0, dbl_mb(k_tj), 1)
                 call daxpy(dim1,-sina,dbl_mb(k_veci),1,dbl_mb(k_tj),1)
                 call daxpy(dim1,cosa,dbl_mb(k_vecj),1,dbl_mb(k_tj),1)
c
               call ga_put(g_vecs(set), 1, dim1, i, i, dbl_mb(k_ti),1)
               call ga_put(g_vecs(set), 1, dim1, j, j, dbl_mb(k_tj),1)
               else
c
c     for beta we swap 1 and 2
c
               call ga_get(g_vecs(2), 1, dim1, i, i, dbl_mb(k_veci),1)
               call ga_get(g_vecs(2), 1, dim1, j, j, dbl_mb(k_vecj),1)
c psi_1 =  -sin*psi_1 + cos*psi_2
                 call dcopy(dim1, 0.d0, 0, dbl_mb(k_ti), 1)
                 call daxpy(dim1,-sina,dbl_mb(k_veci),1,dbl_mb(k_ti),1)
                 call daxpy(dim1,cosa,dbl_mb(k_vecj),1,dbl_mb(k_ti),1)
c psi_2 = cos*psi_1 + sin*psi_2
                 call dcopy(dim1, 0.d0, 0, dbl_mb(k_tj), 1)
                 call daxpy(dim1,cosa,dbl_mb(k_veci),1,dbl_mb(k_tj),1)
                 call daxpy(dim1,sina,dbl_mb(k_vecj),1,dbl_mb(k_tj),1)
c
               call ga_put(g_vecs(2), 1, dim1, i, i, dbl_mb(k_ti),1)
               call ga_put(g_vecs(2), 1, dim1, j, j, dbl_mb(k_tj),1)
                  endif
                  if (.not. ma_chop_stack(l_ti)) call errquit(
     M                 'moswap: chop stack?', 0,MA_ERR)
            else
               call ga_put(g_vecs(set), 1, dim1, i, i, 
     D              dbl_mb(k_vecj),dim1)
               call ga_put(g_vecs(set), 1, dim1, j, j, 
     D              dbl_mb(k_veci),dim1)
            endif
c
            enddo
            if (oprint) call util_flush(luout)
            if (.not. ma_pop_stack(l_vecj)) call errquit('ms:ma?', 0,
     &       MA_ERR)
            if (.not. ma_pop_stack(l_veci)) call errquit('ms:ma?', 0,
     &       MA_ERR)
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
     $        ('movecs_swap: delete of swap entry failed', 0, RTDB_ERR)
c
         if (.not. ma_free_heap(handle)) call errquit('ms:ma?', 0,
     &       MA_ERR)
 1000    continue
      enddo
c
      end


      subroutine movecs_swap_nuclei(rtdb,basis,geom,scftype,g_vecs)
      implicit none
#include "errquit.fh"
c
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "inp.fh"
#include "bas.fh"
#include "geom.fh"
#include "stdio.fh"
c
c     In the database may be an entry 'reorder' which is a list
c     containing the new ordering of the atoms. The MO blocks over each atoms
c     range of basis functions will be reordered according to this list.
c
c     Apply these to the input vectors only (eigenvalues and
c     occupation numbers are unchanged since the MO ordering is
c     not changed). In the case of UHF, g_vecs is assumed to be
c     an array with 2 handles.
c
c     When finished delete the databse entries so that
c     a restart does not re-reorder the vectors.
c
      character*(*) scftype
      integer rtdb, basis,g_vecs(*), nbf, nat, geom
      integer type, nelem, handle, index, g_temp
      integer i, j, k, set, nset
      integer l_atm, k_atm, l_batm, k_batm
      integer jlo, jhi
      integer dim1, dim2, start_from, end_from, start_to, end_to
      logical oprint
c
      oprint = util_print('mo guess',print_default)  
c
      if (scftype .eq. 'UHF') then
        nset = 2
      else
        nset = 1
      endif
c
      if (.not.bas_numbf(basis,nbf))
     $     call errquit('movecs_swap_nuclei: bas_numbf failed',0,0)
c
      if (.not. geom_ncent(geom, nat))
     &     call errquit('movecs_swap_nuclei: geom_ncent failed',0,0)
c
c  set up the list of atoms to reorder
c  -----------------------------------
c
      if (.not.ma_push_get(mt_int,nat,'atom list',l_atm, k_atm))
     &    call errquit('movecs_swap_nuclei: ma_push failed',0,0)
c
      if (.not.ma_push_get(mt_int,nat,'nbf_at list',l_batm, k_batm))
     &    call errquit('movecs_swap_nuclei: ma_push failed',0,0)
c
      if (.not. rtdb_ma_get(rtdb, 'reorder', type, nelem,handle)) 
     $   call errquit('movecs_swap_nuclei: invalid list of pairs',
     $   nelem, 0)
      if (.not. ma_get_index(handle,index)) call errquit
     $    ('movecs_swap_nuclei: get index failed', handle , 0)
      if (type.ne.mt_int) call errquit
     $    ('movecs_swap_nuclei: invalid list of pairs',nelem , 0)
      if (nelem.ne.nat) call errquit
     $    ('movecs_swap_nuclei: reorder list.ne.#atoms',0 , 0)
c
      do i = 1, nelem
         j = int_mb(index+i-1)
           if (i.lt.0.or.i.gt.nat) 
     &       call errquit('movecs_swap_nucei: invalid vector',0,0)
           if (j.lt.0.or.j.gt.nat) 
     &       call errquit('movecs_swap_nucei: invalid vector',0,0)
         int_mb(k_atm+i-1) = j
      enddo
c
c  find nbf for each atom 
c  ----------------------
c
      if (ga_nodeid().eq.0) then
      call util_print_centered(LuOut,'Reordering MO coeffs on atoms',25,
     $        .true.)
      write(LuOut,*) 
      write(LuOut,101) 'old','new','old','new'
      write(LuOut,102) 'atom#','atom#','basis range','basis range'
      write(LuOut,*) '-------------------------------------------------'
      endif
c
      int_mb(k_batm) = 0
      do i = 1,nat
        j = int_mb(k_atm+i-1)
        if (.not.bas_ce2bfr(basis,i,jlo,jhi)) 
     &    call errquit('movecs_swap_nuclei: bas_ce2bfr',0,0)
        int_mb(k_batm+j-1) = jhi-jlo+1
      enddo
c
c  reorder atoms in each set
c  -------------------------
c
      do set = 1,nset
      call ga_inquire(g_vecs(set),type,dim1,dim2)
c 
      if ( .not. ga_create( mt_dbl, dim1, dim2, 'temp vectors',
     $    0,0,g_temp) ) 
     $    call errquit('movecs_swap_nuclei: ga_create failed',0, 0)
      call ga_zero(g_temp)
c
      start_to = 1
c
      do i= 1, nat
      j = int_mb(k_atm+i-1)
c
        start_from = 0
        do k = 1,j-1
          start_from = start_from + int_mb(k_batm+k-1)
        enddo
        start_from = start_from + 1
c
        end_from = start_from + int_mb(k_batm+j-1) - 1
c
        end_to = start_to + int_mb(k_batm+j-1) - 1

        if (ga_nodeid().eq.0.and.set.eq.1) then
        write(LuOut,103) j,i,start_from,end_from,start_to,end_to
        endif
c
        call ga_copy_patch('N',g_vecs(set),start_from,end_from,1,dim2,
     &       g_temp,start_to,end_to,1,dim2)
c
        start_to = start_to + int_mb(k_batm+j-1)
       enddo ! nat
c
      call ga_copy(g_temp, g_vecs(set))
c
      if ( .not. ga_destroy(g_temp))
     $    call errquit('movecs_swap_nuclei: ga_destroy failed',0, 0)
c
      enddo  ! set
c
      call ga_sync()
c
  101 format(2x,a3,4x,a3,11x,a3,16x,a3)
  102 format(1x,a5,2x,a5,7x,a11,8x,a11)
  103 format(i4,2x,i4,5x,i6,2x,i6,5x,i6,2x,i6)
        if (.not. ma_pop_stack(l_batm))
     $    call errquit('movecs_swap_nuclei: ma pop failed',0, 0)
c
        if (.not. ma_pop_stack(l_atm))
     $    call errquit('movecs_swap_nuclei: ma pop failed',0, 0)
c
       if (.not.rtdb_delete(rtdb,'reorder') )
     &   call errquit('movecs_swap_nuclei: rtdb delete failed',0,0)
c
       if (.not. ma_free_heap(handle)) 
     &   call errquit('movecs_swap_nuclei: cant free heap',0,0)
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
#include "errquit.fh"
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "bas.fh"
#include "geom.fh"
#include "tcgmsg.fh"
#include "bq.fh"
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
      character*32 filebqsum, bqsum
      logical oconverged
      character*40 key
      character*26 date, filedate
      integer geom, ma_type, nelem
      integer bq_handle         ! Handle of the BQ fragment
      integer unitno            ! Unit no. for reading
      integer inntsize
      parameter (unitno = 67)   ! These need to be managed !!!
      logical status, oldmode
c
      status = .false.
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
c
      oldmode =  rtdb_parallel(.false.)
      if (ga_nodeid() .eq. 0) then
         open(unitno, status='old', form='unformatted',
     $        file=filename, err=1000)
         read(unitno, end=1001, err=2001)
     $        filebasissum, filegeomsum, filebqsum, filescftype,
     $        filedate
         goto 3001
c
c        If the previous read statement failed we may have a movecs
c        file from an old version of the code. So try reading the 
c        old format (but this means that the vectors cannot be
c        considered as converged). Non-convergence of the vectors
c        is detected from the filebqsum value which we set to something
c        unusual so it will not match anything obvious.
c
 2001    continue
         filebqsum = " "
         rewind(unitno, err=1001)
         read(unitno, end=1001, err=1001)
     $     filebasissum, filegeomsum, filescftype, filedate
 3001    continue
         close(unitno, err=1002)
c     
c     a)
c     
         key = ' '
         write(key,'(a,'':converged'')') scftype(1:inp_strlen(scftype))
         if (.not. rtdb_get(rtdb, key, mt_log, 1, oconverged)) goto 10
         if (.not. oconverged) goto 10
         if (.not. rtdb_get_info(rtdb, key, ma_type, nelem, date))
     $        call errquit('movecs_converged: rtdb corrupt?',0,
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
     $        call errquit('movecs_converged: bad basis',0, BASIS_ERR)
         if (.not. bas_checksum(basis, basissum))
     $        call errquit('movecs_converged: bad basis',0, BASIS_ERR)
         if (.not. geom_checksum(geom, geomsum))
     $        call errquit('movecs_converged: bad geometry',0,
     &       BASIS_ERR)
         if (.not. bq_get_active(bq_handle)) then
              bqsum = ' '
         else
           if (.not. bq_checksum(bq_handle,bqsum)) then
              bqsum = ' '
           endif
         endif
         if (basissum.ne.filebasissum .or. geomsum.ne.filegeomsum .or.
     $       bqsum.ne.filebqsum) 
     $        goto 10
c     
c     Cool ... it is converged!
c
         status = .true.
      endif
c
 10   oldmode =  rtdb_parallel(oldmode)
      call ga_brdcst(1311, status, inntsize, 0)
      movecs_converged = status
      return
c
 1000 call errquit('movecs_converged: error opening file',0, DISK_ERR)
 1001 call errquit('movecs_converged: error reading file',0, DISK_ERR)
 1002 call errquit('movecs_converged: error closing file',0, DISK_ERR)
c
      end
      logical function file_read_ga_info(fname,
     $     title, nrow, ncol)
      implicit none
#include "global.fh"
#include "tcgmsg.fh"
#include "inp.fh"
#include "eaf.fh"
#include "stdio.fh"
#include "mafdecls.fh"
c
      character*(*) fname
      character*(*) title
      integer nrow, ncol
c
c     Read the header info from the dumb file_ga file
c
      integer unitno                 ! Unit no. for writing
      integer len1, ok
      integer ierr
      double precision offset
      integer lgth
      character*255 errmsg
      integer inntsize
c
      ok = 1
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      if (ga_nodeid() .eq. 0) then
         offset=0d0
         ierr=eaf_open(fname(1:inp_strlen(fname)),
     ,        eaf_rw, unitno)
        if (ierr .ne. 0) then
           call eaf_errmsg(ierr,errmsg)
           write(LuOut,*) ga_nodeid(),errmsg
           goto 1000
        endif
        lgth=ma_sizeof(MT_INT,1,MT_BYTE)
        ierr = eaf_read(unitno, offset, len1,lgth)
        offset=offset+lgth
        lgth=len1
#ifdef WIN32        
        title=' '
#else
        ierr = eaf_read(unitno, offset, title(1:len1),lgth)
#endif
        offset=offset+lgth
        lgth=ma_sizeof(MT_INT,1,MT_BYTE)
        ierr = eaf_read(unitno, offset, nrow,lgth)
        offset=offset+lgth
        lgth=ma_sizeof(MT_INT,1,MT_BYTE)
        ierr = eaf_read(unitno, offset, ncol,lgth)
        if (ierr .ne. 0) then
           call eaf_errmsg(ierr, errmsg)
           write(LuOut,*) ' IO offset ', offset
           write(LuOut,*) ' IO error message ',
     $          errmsg(1:inp_strlen(errmsg))
           goto 1001
        endif
        ierr=eaf_close(unitno)
        if (ierr .ne. 0) then
           write(luout,*) ga_nodeid(),' closing FD =',unitno,
     ,          fname
            call eaf_errmsg(ierr, errmsg)
            write(LuOut,*) ' IO error message ',
     $           errmsg(1:inp_strlen(errmsg))
            call util_flush(luout)
            goto 1002
         endif
         goto 10
c
 1000    write(6,*) ' file_read_ga_info: failed to open ',
     $        fname(1:inp_strlen(fname))
         call util_flush(luout)
         ok = 0
         goto 10
c     
 1001    write(6,*) ' file_read_ga_info: failing reading from ',
     $        fname(1:inp_strlen(fname))
         call util_flush(luout)
         ok = 0
        ierr=eaf_close(unitno)
        if (ierr .ne. 0) then
           write(luout,*) ga_nodeid(),' closing FD =',unitno,
     ,          fname
            call eaf_errmsg(ierr, errmsg)
            write(LuOut,*) ' IO error message ',
     $           errmsg(1:inp_strlen(errmsg))
            call util_flush(luout)
            goto 1002
         endif
         goto 10
c     
 1002    write(6,*) ' file_read_ga_info: failed to close',
     $        fname(1:inp_strlen(fname))
         call util_flush(luout)
         ok = 0
         goto 10
      end if
c
 10   call ga_sync()
      call ga_brdcst(1, ok, inntsize, 0) ! Propagate status
      call ga_brdcst(2, nrow, inntsize, 0)
      call ga_brdcst(3, ncol, inntsize, 0)
      call ga_brdcst(31, len1, inntsize, 0)
      call util_char_ga_brdcst(4, title(1:len1), 0)
cedo      call util_char_ga_brdcst(4, title, 0)
c
      file_read_ga_info = ok .eq. 1
c
      end
      logical function file_write_ga_patch(fname, g_a, 
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
#include "eaf.fh"
#include "stdio.fh"
c
c     Generic routine for sequential write 
c     of global array to file
c
      character*(*) fname            ! [input] Name of file to write to
      integer g_a                    ! [input] Array to write
      integer ilo, ihi, jlo, jhi     ! [input] Patch to write
c
      integer unitno                 ! Unit no. for writing
      character*255 gtitle
      integer l_tmp, k_tmp
      integer len1, nrow, ncol
      integer ok, j
      integer ierr
      double precision offset
      integer lgth
      integer inntsize
      character*255 errmsg
c
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_sync
c
      ok = 0
      nrow = ihi-ilo+1
      ncol = jhi-jlo+1
      call ga_inquire_name(g_a, gtitle)
      len1 = inp_strlen(gtitle)
      if (ga_nodeid() .eq. 0) then
         offset=0d0
         ierr=eaf_open(fname(1:inp_strlen(fname)),eaf_rw, unitno)
        if (ierr .ne. 0) then
           call eaf_errmsg(ierr,errmsg)
           write(LuOut,*) ga_nodeid(),errmsg
           call errquit('vectors: eaf_open failed',0, ierr)
        endif
        lgth=ma_sizeof(MT_INT,1,MT_BYTE)
        ierr = eaf_write(unitno, offset, len1,lgth)
        offset=offset+lgth
        lgth=len1
#ifdef WIN32        
        gtitle=' '
#else     
        ierr = eaf_write(unitno, offset, gtitle(1:len1),lgth)
#endif
        offset=offset+lgth
        lgth=ma_sizeof(MT_INT,1,MT_BYTE)
        ierr = eaf_write(unitno, offset, nrow,lgth)
        offset=offset+lgth
        lgth=ma_sizeof(MT_INT,1,MT_BYTE)
        ierr = eaf_write(unitno, offset, ncol,lgth)
        if (ierr .ne. 0) then
           call eaf_errmsg(ierr, errmsg)
           write(LuOut,*) ' IO offset ', offset
           write(LuOut,*) ' IO error message ',
     $          errmsg(1:inp_strlen(errmsg))
           goto 1001
        endif
        offset=offset+lgth
         
        if (.not. ma_push_get(mt_dbl,nrow,'movecs_write',l_tmp,k_tmp))
     $    call errquit('matrix_fwrite: ma failed', nrow, MA_ERR)
        lgth=8d0*nrow
        do j=jlo,jhi
          call ga_get(g_a, ilo, ihi, j, j, dbl_mb(k_tmp), ihi-ilo+1)
          ierr = eaf_write(unitno, offset, dbl_mb(k_tmp),lgth)
          if (ierr .ne. 0) then
             call eaf_errmsg(ierr, errmsg)
             write(LuOut,*) ' IO offset ', offset
             write(LuOut,*) ' IO error message ',
     $            errmsg(1:inp_strlen(errmsg))
             goto 1001
        endif
        offset=offset+lgth
        enddo
        if (.not. ma_pop_stack(l_tmp))
     $    call errquit('file_write_ga: ma pop failed', l_tmp, MA_ERR)

      ierr=eaf_close(unitno)
      if (ierr .ne. 0) then
         write(luout,*) ga_nodeid(),' closing FD =',unitno,
     ,        fname
         call eaf_errmsg(ierr, errmsg)
         write(LuOut,*) ' IO error message ',
     $        errmsg(1:inp_strlen(errmsg))
         call util_flush(luout)
         goto 1002
      endif
        ok = 1
      endif
c
 10   call ga_sync
      call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0)        ! Propagate status
c
      file_write_ga_patch = ok .eq. 1
      if (ga_nodeid() .eq. 0 .and.
     $     util_print('vectors i/o', print_high)) then
         write(6,22) gtitle(1:len1), fname(1:inp_strlen(fname))
 22      format(/' Wrote ',a,' to ',a/)
         call util_flush(luout)
      endif
      return
c
 1000 write(6,*) ' file_write_ga: failed to open ',
     $     fname(1:inp_strlen(fname))
      call util_flush(luout)
      ok = 0
      goto 10
c
 1001 write(6,*) ' file_write_ga: failing writing to ',
     $     fname(1:inp_strlen(fname))
      call util_flush(luout)
      ok = 0
      ierr=eaf_close(unitno)
      if (ierr .ne. 0) then
         write(luout,*) ga_nodeid(),' closing FD =',unitno,
     ,        fname
         call eaf_errmsg(ierr, errmsg)
         write(LuOut,*) ' IO error message ',
     $        errmsg(1:inp_strlen(errmsg))
         call util_flush(luout)
         goto 1002
      endif
      goto 10
c
 1002 write(6,*) ' file_write_ga: failed to close',
     $     fname(1:inp_strlen(fname))
      call util_flush(luout)
      ok = 0
      goto 10
c
      end
      logical function file_read_ga_patch
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
#include "eaf.fh"
#include "stdio.fh"
c
c     Generic routine for sequential read 
c     of file to global array
c
      character*(*) fname            ! [input] Name of file to read from
      integer g_a                    ! [output] Array
      integer ilo, ihi, jlo, jhi     ! [input]
c
      integer unitno                 ! Unit no. for writing
      character*255 gtitle
      integer l_tmp, k_tmp
      integer len1, nrow, ncol
      integer grow, gcol
      integer ok, j
      integer ierr
      double precision offset
      integer lgth
      integer inntsize
      character*255 errmsg
c
      inntsize=ma_sizeof(MT_INT,1,MT_BYTE)
      call ga_sync()
c
      ok = 0
      grow = ihi-ilo+1
      gcol = jhi-jlo+1
      call ga_inquire_name(g_a, gtitle)
      len1 = inp_strlen(gtitle)
      if (ga_nodeid() .eq. 0) then
         offset=0d0
         ierr=eaf_open(fname(1:inp_strlen(fname)),
     ,        eaf_rw, unitno)
         if (ierr .ne. 0) then
           call eaf_errmsg(ierr,errmsg)
           write(LuOut,*) ga_nodeid(),errmsg
           goto 1000
        endif
        lgth=inntsize
        ierr = eaf_read(unitno, offset, len1,lgth)
        offset=offset+lgth
        lgth=len1
#ifdef WIN32        
        gtitle=' '
#else
        ierr = eaf_read(unitno, offset, gtitle(1:len1),lgth)
#endif
        offset=offset+lgth
        lgth=inntsize
        ierr = eaf_read(unitno, offset, nrow,lgth)
        offset=offset+lgth
        lgth=inntsize
        ierr = eaf_read(unitno, offset, ncol,lgth)
        if (ierr .ne. 0) then
           call eaf_errmsg(ierr, errmsg)
           write(LuOut,*) ' IO offset ', offset
           write(LuOut,*) ' IO error message ',
     $          errmsg(1:inp_strlen(errmsg))
           goto 1001
        endif
        offset=offset+lgth
        if ((nrow.ne.grow).and.(ncol.ne.gcol)) goto 2000
        if (.not. ma_push_get(mt_dbl,nrow,'movecs_write',l_tmp,k_tmp))
     $    call errquit('matrix_fwrite: ma failed', nrow, MA_ERR)
        lgth=8d0*nrow
        do j=jlo,jhi
          ierr = eaf_read(unitno, offset, dbl_mb(k_tmp),lgth)
          if (ierr .ne. 0) then
             call eaf_errmsg(ierr, errmsg)
             write(LuOut,*) ' IO offset ', offset
             write(LuOut,*) ' IO error message ',
     $            errmsg(1:inp_strlen(errmsg))
             if (.not. ma_pop_stack(l_tmp))
     $        call errquit('file_read_ga: ma pop failed', l_tmp, MA_ERR)
          goto 1001
        endif
        offset=offset+lgth
        call ga_put(g_a, ilo, ihi, j, j, dbl_mb(k_tmp), ihi-ilo+1)
        enddo
        if (.not. ma_pop_stack(l_tmp))
     $    call errquit('file_read_ga: ma pop failed', l_tmp, MA_ERR)
        ierr=eaf_close(unitno)
        if (ierr .ne. 0) then
           write(luout,*) ga_nodeid(),' closing FD =',unitno,
     ,          fname
            call eaf_errmsg(ierr, errmsg)
            write(LuOut,*) ' IO error message ',
     $           errmsg(1:inp_strlen(errmsg))
            call util_flush(luout)
            goto 1002
         endif
        ok = 1
      endif
c
 10   call ga_sync()
      call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0)        ! Propagate status
c
      file_read_ga_patch = ok .eq. 1
      if (ga_nodeid() .eq. 0 .and.
     $     util_print('vectors i/o', print_high)) then
         write(6,22) gtitle(1:len1), fname(1:inp_strlen(fname))
 22      format(/' Read ',a,' from ',a/)
         call util_flush(luout)
      endif
      return
c
 1000 write(6,*) ' file_read_ga: failed to open ',
     $     fname(1:inp_strlen(fname))
      call util_flush(luout)
      ok = 0
      goto 10
c
 1001 write(6,*) ' file_read_ga: failing reading from ',
     $     fname(1:inp_strlen(fname))
      call util_flush(luout)
      ok = 0
      ierr=eaf_close(unitno)
      if (ierr .ne. 0) then
         write(luout,*) ga_nodeid(),' closing FD =',unitno,
     ,        fname
         call eaf_errmsg(ierr, errmsg)
         write(LuOut,*) ' IO error message ',
     $        errmsg(1:inp_strlen(errmsg))
         call util_flush(luout)
         goto 1002
      endif
      goto 10
c
 1002 write(6,*) ' file_read_ga: failed to close',
     $     fname(1:inp_strlen(fname))
      call util_flush(luout)
      ok = 0
      goto 10
c
 2000 write(6,*) ' file_read_ga: GA and file ',
     $           'contents mismatch dimensions',
     $           nrow, ncol, grow, gcol
      call util_flush(luout)
      ok = 0
      goto 10
c
      end
