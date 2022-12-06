      subroutine mp2_nonsep_uhf( rtdb, geom,
     $     basis, oskel,
     $     nbf, 
     $     nir, 
     $     noa, nva, nva_lo, nva_hi, num_va, num_oa, 
     $     nob, nvb, nvb_lo, nvb_hi, num_vb, num_ob,
     $     sym_lo_oa, sym_hi_oa, sym_lo_va, sym_hi_va,
     $     sym_lo_ob, sym_hi_ob, sym_lo_vb, sym_hi_vb,
     $     oseg_lo, oseg_hi, 
     $     irs_a, irs_b,
     $     c_a, c_b,            ! Better left in a GA to conserve memory
     $     nva_lo_local, nva_hi_local, 
     $     tunita, tunitb, grad,
     P     p_file_size)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "bas.fh"
#include "eaf.fh"
#include "cmp2ps.fh"
#include "geom.fh"
#include "rtdb.fh"
      integer basis, rtdb, geom
      logical oskel
      integer nbf               ! No. of basis functions
      integer nir               ! No. of irreducible representations
      integer noa, nob          ! No. of occupied orbitals
      integer nva, nva_lo, nva_hi ! Number and ranges of virtual orbitals
      integer nvb, nvb_lo, nvb_hi
      integer num_va(0:nir-1), num_vb(0:nir-1) ! No. of vir of each symmetry
      integer num_oa(0:nir-1), num_ob(0:nir-1) ! No. of occ of each symmetry
      integer sym_lo_oa(0:nir-1), sym_hi_oa(0:nir-1),
     $     sym_lo_va(0:nir-1), sym_hi_va(0:nir-1) ! Ranges of each symmetry
      integer sym_lo_ob(0:nir-1), sym_hi_ob(0:nir-1),
     $     sym_lo_vb(0:nir-1), sym_hi_vb(0:nir-1)
      integer oseg_lo, oseg_hi  ! Range of occupied for this pass
      integer irs_a(nbf), irs_b(nbf) ! Orbital symmetries
      double precision c_a(nbf,*), c_b(nbf,*) ! MO coefficients
      integer nva_lo_local, nva_hi_local ! Range of virtuals on this node
      integer tunita, tunitb    ! Unit no.s for pure and mixed spin T
      double precision grad(3,*)
c
c     Allocate memory for back transformation routine
c
      integer g_buf,g_buf_trans
      integer nshpair, nshpairlocal, nbfpair, nbfpairlocal
      integer l_shpairs, k_shpairs, l_shpairslocal, k_shpairslocal
      integer l_shdim, k_shdim, l_shlo, k_shlo
      integer l_t, k_t, l_tmp, k_tmp, l_iauv, k_iauv, l_map, k_map
      integer l_c_t, k_c_t
      integer l_act, k_act, l_actsh, k_actsh
      integer twopdmunit
      integer junk(2), ninseg, ierr, i, j, ish, ishlo, ishhi, shmax,
     t     tdim
      integer nsh, natoms, nactive, nblock
c
      integer k_scr, l_scr, k_lab, l_lab, k_eri, l_eri, leneri, lenscr
      integer k_bftosh, l_bftosh, k_bftoce, l_bftoce, k_pdm, l_pdm
      integer eaftype,eaf_size_in_mb,inntsize
      integer mp2_eaftype
      external mp2_eaftype
      double precision p_file_size
c
      character*(nw_max_path_len) fname
c
      double precision tol2e
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
c
      if (.not. rtdb_get(rtdb, 'mp2:backtol', mt_dbl, 1, tol2e))
     $     tol2e = 1d-9
c
      if (.not. bas_numcont(basis, nsh)) call errquit
     $     ('mp2: backt bad basis handle ', basis, BASIS_ERR)
      ninseg = oseg_hi - oseg_lo + 1
c
c     Make mapping of shell pairs to processors
c
      call int_init(rtdb,1,basis)     ! intd_init overwrite mem estimates
      call schwarz_init(geom,basis)
      call int_terminate()
c
      call mp2_backt_info(basis, tol2e, oskel,
     $     .false., nshpair, nshpairlocal,
     $     nbfpair, nbfpairlocal, shmax, junk, junk, junk, junk, junk)
      if (.not. ma_push_get(mt_int, nshpair*2, 'mp2: shpairs',
     $     l_shpairs, k_shpairs)) call errquit
     $     ('mp2: insufficient memory : shpairs ', 2*nshpair, MA_ERR)
      if (.not. ma_push_get(mt_int, nshpairlocal*2, 'mp2:shpairslocal',
     $     l_shpairslocal, k_shpairslocal)) call errquit
     $     ('mp2: insufficient memory : shpairs ', 2*nshpairlocal,
     &       MA_ERR)
      if (.not. ma_push_get(mt_int, ga_nnodes(), 'mp2:map',
     $     l_map, k_map)) call errquit
     $     ('mp2: insufficient memory : map ', ga_nnodes(), MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'mp2:shlo',
     $     l_shlo, k_shlo)) call errquit
     $     ('mp2: insufficient memory : shlo ', nsh, MA_ERR)
      if (.not. ma_push_get(mt_int, nsh, 'mp2:shdim',
     $     l_shdim, k_shdim)) call errquit
     $     ('mp2: insufficient memory : shdim ', nsh, MA_ERR)
      call mp2_backt_info(basis, tol2e, oskel,
     $     .true., nshpair, nshpairlocal,
     $     nbfpair, nbfpairlocal, shmax,
     $     int_mb(k_shdim), int_mb(k_shlo),
     $     int_mb(k_shpairs), int_mb(k_shpairslocal), int_mb(k_map))
c
c     Allocate remaining memory and open file to hold the 
c     3-parts transformed density matrix.
c
      tdim = max(nbf*nbf, nbf*shmax*shmax)
c         
      if (.not. ma_push_get(mt_dbl, tdim, 'mp2: backt t',
     $     l_t, k_t)) call errquit
     $     ('mp2: insufficient memory : t ', tdim, MA_ERR)
      if (.not. ma_push_get(mt_dbl, tdim, 'mp2: backt tmp',
     $     l_tmp, k_tmp)) call errquit
     $     ('mp2: insufficient memory : tmp ', tdim, MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbf*nbf, 'mp2: backt iauv',
     $     l_iauv, k_iauv)) call errquit
     $     ('mp2: insufficient memory : iauv ', nbf*nbf, MA_ERR)
c
c     Right now we have an N**3 array in GA.  This is easily
c     replaced with an N*B array with B=P*S*S (P procs, S=max shell size)
c     by blocking the transformation.
c
c     For small problems some nodes may not have data in which case
c     map(proc) will be set to a value larger then nbfpair (which is the 
c     actual maximum dimension) ... discard these before creating
c     the array.
c
      nblock = 0
      do i = 1, ga_nnodes()
         if (int_mb(k_map+i-1) .le. nbfpair) nblock = nblock + 1
      end do
c
      if (.not. ga_create_irreg(mt_dbl, nva, nbfpair,
     $     'mp2: backt', 1, 1, int_mb(k_map), nblock, g_buf))
     $     call errquit('mp2: backt: ga_create failed', nva*nbfpair,
     &       GA_ERR)
      if (.not.(ga_create( mt_dbl, nva, nbfpair,'mp2: backt2',
     $                     0, nbfpair, g_buf_trans )))
     $   call errquit('mp2:backt: ga_create2 failed', nva*nbfpair,
     &       GA_ERR)
c     
      call util_file_name('2pdm', .true., .true., fname)
#ifdef NOIO
       eaftype=mp2_eaftype(rtdb,p_file_size)
      if(ga_nodeid().eq.0) write(6,*) ' mp2_eaf for 2pdm  ',eaftype
#else
       eaftype=eaf_rw
#endif
      if (eaf_open(fname,eaftype,twopdmunit) .ne. 0)
     $     call errquit('mp2: backt: failed to open file',0, DISK_ERR)
c     
c     Now do the back transformation
c
      call ga_sync()
      call pstat_on(ps_backt)
      call mp2_back_transform_uhf(
     $     rtdb, 
     $     nbf, 
     $     nir, 
     $     nva, nva_lo, nva_hi, num_va, num_oa, 
     $     nvb, nvb_lo, nvb_hi, num_vb, num_ob,
     $     sym_lo_oa, sym_hi_oa, sym_lo_va, sym_hi_va,
     $     sym_lo_ob, sym_hi_ob, sym_lo_vb, sym_hi_vb,
     $     oseg_lo, oseg_hi, 
     $     irs_a, irs_b,
     $     c_a, c_b, 
     $     nva_lo_local, nva_hi_local, 
     $     tunita, tunitb, twopdmunit,
     $     nshpair, nbfpair, int_mb(k_shpairs),
     $     nshpairlocal, nbfpairlocal, int_mb(k_shpairslocal),
     $     int_mb(k_shdim), int_mb(k_shlo),
     $     dbl_mb(k_t), dbl_mb(k_iauv), dbl_mb(k_tmp),
     $     int_mb(k_map), g_buf, g_buf_trans)
      call pstat_off(ps_backt)
      call ga_sync()
c
c     Free up some memory and then contract with gradient integrals.
c
      if (.not. ga_destroy(g_buf)) call errquit
     $     ('mp2: backt: ga destroy failed', 0, GA_ERR)
      if (.not. ga_destroy(g_buf_trans)) call errquit
     $     ('mp2: backt: ga destroy failed', 1, GA_ERR)
      if (.not. ma_chop_stack(l_t)) call errquit
     $     ('mp2:backt: first chop stack failed', 0, MA_ERR)
c
      tdim = shmax*shmax*nbf*ninseg
      if (.not. ma_push_get(mt_dbl, tdim, 'mp2:back tt',
     $     l_t, k_t)) call errquit
     $     ('mp2: backt: failed ma for nonsep test', tdim, MA_ERR)
      tdim = shmax*shmax*nbf
      if (.not. ma_push_get(mt_dbl, tdim*ninseg, 'mp2:back tbuf',
     $     l_tmp, k_tmp)) call errquit
     $     ('mp2: backt: failed ma for nonsep test', tdim, MA_ERR)
      if (.not. ma_push_get(mt_dbl, ninseg*nbf, 'mp2:back c_t',
     $     l_c_t, k_c_t)) call errquit
     $     ('mp2: backt: failed ma for nonsep test', ninseg*nbf, MA_ERR)
      do i = oseg_lo, oseg_hi
         call dcopy(nbf, c_a(1,i), 1, dbl_mb(k_c_t+i-oseg_lo), ninseg)
      end do
c
c     Determine list of active centers
c
      if (.not. geom_ncent(geom, natoms)) call errquit
     $     ('mp2_backt: geom ?',0, GEOM_ERR)
      if (.not. ma_push_get(mt_log, natoms, 'mp2:back act',
     $     l_act, k_act)) call errquit('mp2:back ma ', natoms, MA_ERR)
      call grad_active_atoms(rtdb, natoms, log_mb(k_act), nactive)
c
c     Turn this into a list of active shells
c
      if (.not. ma_push_get(mt_log, nsh, 'mp2:back actsh',
     $     l_actsh, k_actsh)) call errquit('mp2:back ma ', nsh, MA_ERR)
      do i = 1, natoms
         if (.not. bas_ce2cnr(basis, i, ishlo, ishhi))
     $        call errquit('mp2:backt basis?',0, BASIS_ERR)
         do ish = ishlo,ishhi
            log_mb(k_actsh+ish-1) = log_mb(k_act+i-1)
         end do
      end do
c
      call intd_init(rtdb, 1, basis)
      call intb_mem_2e4c(leneri, lenscr) ! blocking algorithm
      leneri = leneri/12
      leneri = max(leneri,15**4) ! 1 G quartet = 39 D quartets ... SEE MP2_MEMORY
      if (.not. ma_push_get(mt_dbl,12*leneri,'deriv buffer',
     $     l_eri,k_eri)) call errquit
     $     ('mp2:backt could not allocate buffer',12*leneri, MA_ERR)
      if (.not. ma_push_get(mt_dbl,shmax**4,'pdm buffer',
     $     l_pdm,k_pdm)) call errquit
     $     ('mp2:backt could not allocate pdm',shmax*4, MA_ERR)
      if (.not. ma_push_get(mt_int,4*leneri,'deriv labels',l_lab,k_lab))
     $     call errquit('mp2:backt could not allocate labels',leneri,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl,lenscr,'deriv scratch',
     $     l_scr,k_scr))
     $     call errquit('mp2:backt could not allocate scr',lenscr,
     &       MA_ERR)
      if (.not. ma_push_get(mt_int,nbf,'shmap',l_bftosh,k_bftosh))
     $     call errquit('mp2:backt could not allocate bftosh',nbf,
     &       MA_ERR)
      if (.not. ma_push_get(mt_int,nbf,'cemap',l_bftoce,k_bftoce))
     $     call errquit('mp2:backt could not allocate bftoce',nbf,
     &       MA_ERR)
c      
      call ga_sync()
      call pstat_on(ps_nonsep)
      call mp2_nonsep(rtdb, basis,
     $     twopdmunit,
     $     dbl_mb(k_c_t),
     $     nsh, nbf, ninseg,
     $     nshpairlocal, int_mb(k_shpairslocal),
     $     int_mb(k_shdim), int_mb(k_shlo),
     $     dbl_mb(k_t), dbl_mb(k_tmp), grad,
     $     tol2e, oskel,
     $     log_mb(k_actsh), leneri, dbl_mb(k_eri), int_mb(k_lab),
     $     lenscr, dbl_mb(k_scr), int_mb(k_bftosh), int_mb(k_bftoce),
     $     dbl_mb(k_pdm))
      call pstat_off(ps_nonsep)
      call ga_sync()
      call intd_terminate()
c
c     Zero out gradients on inactive atoms
c
      do i = 1, natoms
         if (.not. log_mb(k_act+i-1)) then
            do j = 1, 3
               grad(j,i) = 0.0d0
            end do
         end if
      end do
c
      call schwarz_tidy()
c
c     Done.
c     
      if (.not. ma_chop_stack(l_shpairs)) call errquit
     $     ('mp2: backt: failed chopping stack',0, MA_ERR)
c
      if (util_print('iostats', print_high) .and. 
     $     ga_nodeid().eq.0) call eaf_print_stats(twopdmunit)
      if (eaf_close(twopdmunit) .ne. 0)
     $     call errquit('mp2: backt: closing 2pdm',ierr, DISK_ERR)
      call util_file_unlink(fname)
c
      end
      subroutine mp2_back_transform_uhf(
     $     rtdb, 
     $     nbf, 
     $     nir, 
     $     nva, nva_lo, nva_hi, num_va, num_oa, 
     $     nvb, nvb_lo, nvb_hi, num_vb, num_ob,
     $     sym_lo_oa, sym_hi_oa, sym_lo_va, sym_hi_va,
     $     sym_lo_ob, sym_hi_ob, sym_lo_vb, sym_hi_vb,
     $     oseg_lo, oseg_hi, 
     $     irs_a, irs_b,
     $     c_a, c_b, 
     $     nva_lo_local, nva_hi_local, 
     $     tunita, tunitb, twopdmunit,
     $     nshpair, nbfpair, shpairs,
     $     nshpairlocal, nbfpairlocal, shpairslocal,
     $     shdim, shlo,
     $     t, ia_uv, tmp,
     $     map, g_buf, g_buf_trans)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "eaf.fh"
#include "util.fh"
#include "rtdb.fh"
      integer rtdb
      integer nbf               ! No. of basis functions
      integer nir               ! No. of irreducible representations
      integer nva, nva_lo, nva_hi ! Number and ranges of virtual orbitals
      integer nvb, nvb_lo, nvb_hi
      integer num_va(0:nir-1), num_vb(0:nir-1) !  No. of vir of each symmetry
      integer num_oa(0:nir-1), num_ob(0:nir-1) !  No. of occ of each symmetry
      integer sym_lo_oa(0:nir-1), sym_hi_oa(0:nir-1),
     $     sym_lo_va(0:nir-1), sym_hi_va(0:nir-1) ! Ranges of each symmetry
      integer sym_lo_ob(0:nir-1), sym_hi_ob(0:nir-1),
     $     sym_lo_vb(0:nir-1), sym_hi_vb(0:nir-1)
      integer oseg_lo, oseg_hi  ! Range of occupied for this pass
      integer irs_a(nbf), irs_b(nbf) ! Orbital symmetries
      double precision c_a(nbf,*), c_b(nbf,*) ! MO coefficients
      integer nva_lo_local, nva_hi_local ! Range of virtuals on this node
      integer tunita, tunitb    ! Unit no.s for pure and mixed spin T
      integer nshpair           ! total no. of non-zero shell pairs
      integer nshpairlocal      ! no. of pairs assigned to this processor
      integer nbfpair           ! dimension of bf pairs
      integer nbfpairlocal      ! dimension of local bf pairs
      integer shpairs(2,nshpair) ! (u>=v) shell pair indices
      integer shpairslocal(2,*) ! local shell pair indices
      integer shdim(*)          ! no. of bf in shell
      integer shlo(*)           ! first bf in shell
      double precision t(*)     ! Scratch max(nbf*nbf,nbf*S*S)
      double precision tmp(*)   ! Scratch max(nbf*nbf,nbf*S*S)
      double precision ia_uv(nbf,nbf) ! Scratch
      integer twopdmunit        ! Unit no. for part transformed two-PDM
      integer map(0:*)          ! First index of g_buf on each node
      integer g_buf             ! Global array handle for transpose buffer
      integer g_buf_trans       ! Global array handle for tr-transpose buffer
c     
*     double precision tsum

      double precision zero, one, scale
      double precision tunitptra, tunitptrb ! Pointers into files
      double precision twopdmunitptr
      integer count, tcount, ind,  ptr
      integer symi, symj, syma, symb, symia
      integer i, a, u, v, ush, vsh, udim, vdim, ishpair
c
      logical mp2cpybck
      logical scs  ! flag for SCS-MP2
      double precision fss,fos ! scaling factors for SCS-MP2
c
#include "bitops.fh"
c     
      zero=0.0d0
      one=1.0d0
c     
      tunitptra=1
      tunitptrb=1
*     tsum = 0
c     flag and scaling factors for SCS-MP2

      if (.not. rtdb_get(rtdb, 'mp2:scs', mt_log, 1, scs))
     $  scs = .false.       ! flag for SCS gradient
      if(scs)then
        if (.not. rtdb_get(rtdb, 'mp2:fss', mt_dbl, 1, fss))
     $    fss = 1.0d0/3.0d0   ! same spin scaling factor for SCS-MP2
        if (.not. rtdb_get(rtdb, 'mp2:fos', mt_dbl, 1, fos))
     $    fos = 1.2d0         ! opposite spin scaling factor for SCS-MP2
      else
        if (.not. rtdb_get(rtdb, 'mp2:fss', mt_dbl, 1, fss))
     $    fss = one
        if (.not. rtdb_get(rtdb, 'mp2:fos', mt_dbl, 1, fos))
     $    fos = one
        if ((fss.ne.one).or.(fos.ne.one)) then
          if (ga_nodeid().eq.0) then
            write(6,*) "You set fss or fos but not scs; setting scs=T"
          endif
          scs = .true.
        endif
      endif
c     
      do i=oseg_lo,oseg_hi 
cedo	 call ga_fill_patch(g_buf,1,nva,1,nbfpair,99.0d0) ! For debug
         symi=irs_a(i)
         do a=nva_lo_local,nva_hi_local
            syma=irs_a(a)
            symia=ieor(syma,symi)
            call dfill((nbf*nbf),zero,ia_uv,1) ! Will add pure and mixed in
c     
            call mp2_read_tijab(nva_lo, nva_hi, irs_a, symia,
     $           num_oa, sym_hi_oa, sym_lo_oa, tunita, tunitptra, t)
c     
            tcount=1
            do symb=0,nir-1
               symj=ieor(symia,symb)
               if(num_va(symb).gt.0.and.num_oa(symj).gt.0) then
                  call dgemm('n','t', !             t(j,b)Cbv -> t(j,v)
     $                 num_oa(symj), nbf, num_va(symb),
     $                 fss, t(tcount), num_oa(symj),
     $                 c_a(1,sym_lo_va(symb)), nbf,
     $                 zero, tmp, num_oa(symj))
                  call dgemm('n','n', !             Cuj t(j,v) -> t(u,v)
     $                 nbf, nbf, num_oa(symj),
     $                 one, c_a(1,sym_lo_oa(symj)), nbf,
     $                 tmp, num_oa(symj),
     $                 one, ia_uv, nbf)
                  tcount=tcount+num_oa(symj)*num_va(symb)
               end if
            end do
c     
            call mp2_read_tijab(nvb_lo, nvb_hi, irs_b, symia,
     $           num_ob, sym_hi_ob, sym_lo_ob, tunitb, tunitptrb, t)
c     
            tcount=1
            do symb=0,nir-1
               symj=ieor(symia,symb)
               if(num_vb(symb).gt.0.and.num_ob(symj).gt.0) then
                  call dgemm('n','t', !             t(j,b)Cbv -> t(j,v)
     $                 num_ob(symj), nbf, num_vb(symb),
     $                 fos, t(tcount), num_ob(symj),
     $                 c_b(1,sym_lo_vb(symb)), nbf,
     $                 zero, tmp, num_ob(symj))
                  call dgemm('n','n', !             Cuj t(j,v) -> t(u,v)
     $                 nbf, nbf, num_ob(symj),
     $                 one, c_b(1,sym_lo_ob(symj)), nbf,
     $                 tmp, num_ob(symj),
     $                 one, ia_uv, nbf)
                  tcount=tcount+num_ob(symj)*num_vb(symb)
               end if
            end do
c     
c     Pack down into sparse, symmetry unique list and symmetrize over uv
c     Note we put in the whole square of diagonal shell blocks and divide
c     it by two.
c     
            ind = 0
            do ishpair = 1, nshpair
               ush = shpairs(1,ishpair)
               vsh = shpairs(2,ishpair)
               scale = 1.0d0
               if (ush.eq.vsh) scale = 0.5d0
               do u = shlo(ush), shlo(ush)+shdim(ush)-1
                  do v = shlo(vsh), shlo(vsh)+shdim(vsh)-1
                     ind = ind + 1
                     tmp(ind) = (ia_uv(u,v)+ia_uv(v,u))*scale
                  end do
               end do
            end do
            if (ind .ne. nbfpair) call errquit('mp2bt: ind?', ind,
     &       UNKNOWN_ERR)
c     
c     shpairs(1:2,1:nshpair)    (u>=v) shell pair indices
c     shdim(1:nsh)              dimension of this sh (s=1,p=2,...)
c     shlo(1:nsh)               first basis function in shell
c     nbfpair                   total no. of non-zero bf pairs
c     nshpair                   total no. of non-zero shell pairs
c     nshpairlocal              no. of shellpairs assigned to this processor
c     nbfpairlocal              no. of bfpairs assigned to this processor
c     shpairslocal(1:2,1:nshpairlocal) (u>=v) local shell pair indices
c     
            call ga_put(g_buf_trans, a-nva_lo+1, a-nva_lo+1,
     $           1, nbfpair, tmp, 1)
c     
         end do                 ! End of a
         call ga_sync()
        if (.not. rtdb_get(rtdb, 'mp2:copyback', mt_log, 1, mp2cpybck))
     $    mp2cpybck=.false.

         if(mp2cpybck) then
            call mp2_copyback(g_buf,g_buf_trans)
         else
            call ga_copy(g_buf_trans,g_buf)
         endif
         call ga_sync()
c     
c     Now have locally (a,1:nbfpairlocal).  Loop thru local shell pairs
c     transform the a index and write to disk.
c     
         ind = map(ga_nodeid())
         ptr = 0
         do ishpair = 1, nshpairlocal
            ush = shpairslocal(1,ishpair)
            vsh = shpairslocal(2,ishpair)
            udim = shdim(ush)
            vdim = shdim(vsh)
c     
            call ga_get(g_buf,1,nva,ind,ind+udim*vdim-1,tmp,nva)
	    ind = ind + udim*vdim
            call dgemm('n','n', ! C(s,a)*F(a,uv) -> F(s,uv)
     $           nbf, udim*vdim, nva,
     $           one, c_a(1,nva_lo), nbf,
     $           tmp, nva,
     $           zero, t, nbf)
c     
c     On disk is D(1:nbf,1:vdim(v),1:udim(u),oseg_lo:oseg_hi,1:nshpairlocal)
c     with indices D(s,vf,uf,i,ishpair)
c     
c     Note - this presumes that the routine reading the data back in will
c     be using an algorithm with o*n*s*s memory requirement (n=#bf, s=#bf 
c     in shell, o=occupied).
c     
            count = nbf*udim*vdim ! Amount of data to write
            twopdmunitptr = 8.0d0*(ptr + (i-oseg_lo)*count)
c     
            if (eaf_write(twopdmunit, twopdmunitptr, t, count*8).ne.0)
     $           call errquit('mp2_bt: write of two pdm?', 0,
     &       DISK_ERR)
c     
            ptr = ptr + count*(oseg_hi-oseg_lo+1)
         end do
c     
         call ga_sync
c     
      end do                    ! End of i
c     
      end
      subroutine mp2_backt_info(
     $     basis,
     $     tol2e, oskel,
     $     omakearrays,
     $     nshpair, nshpairlocal, 
     $     nbfpair, nbfpairlocal, shmax,
     $     shdim, shlo, shpairs, shpairslocal, map)
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "schwarz.fh"
#include "global.fh"
#include "util.fh"
#include "sym.fh"
      integer basis             ! [input] basis set handle
      double precision tol2e    ! [input] screening threshold
      logical oskel             ! [input] if true use skeleton symm
      logical omakearrays       ! [input] if true then make the arrays
      integer nshpair           ! [output] no. of non-zero shell pairs
      integer nshpairlocal      ! [output] no. of local shell pairs
      integer shmax             ! [output] Max AO shell dimension
      integer nbfpair           ! [output] sum of all non-zero pair dims
      integer nbfpairlocal      ! [output] sum of local non-zero pair dims
      integer shdim(*)          ! [output] dimension of each shell
      integer shlo(*)           ! [output] first bf in each shell
      integer shpairs(2,*)      ! [output] 1->u, 2->v u>=v shells in pair
      integer shpairslocal(2,*) ! [output] ditto but only for local pairs
      integer map(0:*)          ! [output] Map for g_create_irreg for g_buf
c     
c     Form a list of interacting shell pairs (u>=v) and assign them
c     to processors.  Ideally do this so as to optimize the efficiency
c     of derivative integral evaluation and to provide good load balance 
c     for the back transformation.  Right now we just do round-robin.
c     
c     If (omakearrays) then 
c     actually make the arrays
c     else
c     just return the scalar results
c     end if
c     
c     The Schwarz package must be initialized before entry
c     
      integer u, v, nbf, nsh, lo, hi, udim, vdim, me, nproc, owner
      integer count, tmp, ishpair, first, last, iproc, uvdim
      logical odoit
      double precision q2, sss
c     
      integer npairblock        ! Crude attempt to load balance the pairs
      parameter (npairblock = 7)
      integer pairblock(npairblock), ipairb, npairmin, npairmax
      integer logmin, logmax
      intrinsic log10
      data pairblock /1,3,6,9,18,36,1000000/ ! Selected from s*s,s*p,p*p ...
c
      if (.not. bas_numbf(basis, nbf))
     $     call errquit('mp2_backt_info: bad basis handle', 0,
     &       BASIS_ERR)
      if (.not. bas_numcont(basis, nsh))
     $     call errquit('mp2_backt_info: bad basis handle', 0,
     &       BASIS_ERR)
c     
      me = ga_nodeid()
      nproc = ga_nnodes()
      shmax = 0
c     
      if (omakearrays) then
c     
c     First compute offset to first shell pair on each processor
c     
         nshpair = 0
         call ifill(nproc, 0, map, 1)
         logmax = 1000000
         do logmin = -1,-13,-4  ! NOTE THAT SCREENING IS IMPLIED IN THIS LOOP
            npairmin = 0
            do ipairb = 1, npairblock ! Loop thru types of pairs
               npairmax = pairblock(ipairb)
               do u = nsh, 1, -1
                  if (.not. bas_cn2bfr(basis, u, lo, hi))
     $                 call errquit('mp2_backt_info: bas range', u,
     &       BASIS_ERR)
                  udim = hi - lo + 1
                  do v = 1, u
                     if (.not. bas_cn2bfr(basis, v, lo, hi))
     $                    call errquit('mp2_backt_info:bas range',v,
     &       BASIS_ERR)
                     vdim = hi - lo + 1
                     uvdim = udim*vdim
                     sss = schwarz_shell(u,v)*schwarz_max() + 1d-100
                     if (uvdim.gt.npairmin .and. uvdim.le.npairmax .and.
     $                    log10(sss).lt.logmax .and. 
     $                    log10(sss).ge.logmin) then
c     
                        odoit = sss .gt. tol2e*0.1d0
                        if (odoit .and. oskel) 
     $                       odoit = sym_shell_pair(basis,u,v,q2)
                        if (odoit) then
                           owner = mod(nshpair,nproc)
                           map(owner) = map(owner) + 1
                           nshpair = nshpair + 1
                        end if
c     
                     end if
                  end do
               end do
               npairmin = npairmax
            end do
            logmax = logmin
         end do
         count = 1
         do u = 0, nproc-1
            tmp = map(u)
            map(u) = count
            count = count + tmp
         end do
      end if
c     
c     Make list of all shell pairs ordered by their assigned processor,
c     the list of local shell pairs, and also count the other information.
c     
      nshpairlocal = 0
      nbfpair = 0
      nbfpairlocal = 0
      nshpair = 0
      logmax = 1000000
      do logmin = -1,-13,-4  ! NOTE THAT SCREENING IS IMPLIED IN THIS LOOP
         npairmin = 0
         do ipairb = 1, npairblock ! Loop thru types of pairs
            npairmax = pairblock(ipairb)
            do u = nsh, 1, -1
               if (.not. bas_cn2bfr(basis, u, lo, hi))
     $              call errquit('mp2_backt_info: bas range', u,
     &       BASIS_ERR)
               shmax = max(shmax,hi-lo+1)
               if (omakearrays) then
                  shdim(u) = hi - lo + 1
                  shlo(u) = lo
               end if
               udim = hi - lo + 1
               do v = 1, u
                  if (.not. bas_cn2bfr(basis, v, lo, hi))
     $                 call errquit('mp2_backt_info: bas range', v,
     &       BASIS_ERR)
                  vdim = hi - lo + 1
                  uvdim = udim*vdim
                  sss = schwarz_shell(u,v)*schwarz_max() + 1d-100
                  if (uvdim.gt.npairmin .and. uvdim.le.npairmax .and.
     $                 log10(sss).lt.logmax .and. 
     $                 log10(sss).ge.logmin) then
c     
                     odoit = sss.gt.tol2e*0.1d0
                     if (odoit .and. oskel) 
     $                    odoit = sym_shell_pair(basis,u,v,q2)
                     if (odoit) then
                        owner = mod(nshpair,nproc)
                        if (owner.eq.me) then
                           nshpairlocal = nshpairlocal + 1
                           nbfpairlocal = nbfpairlocal + udim*vdim
                           if (omakearrays) then
                              shpairslocal(1,nshpairlocal) = u
                              shpairslocal(2,nshpairlocal) = v
                           end if
                        end if
                        nshpair = nshpair + 1
                        nbfpair = nbfpair + udim*vdim
                        if (omakearrays) then
                           shpairs(1,map(owner)) = u
                           shpairs(2,map(owner)) = v
                           map(owner) = map(owner) + 1
                        end if
                     end if
                  end if
               end do
            end do
            npairmin = npairmax
         end do
         logmax = logmin
      end do
c     
c     Now replace map with the map creating the global array for transposing
c     
      if (omakearrays) then
         count = 1
         first = 1
         do iproc = 0, nproc-1
            last = map(iproc) - 1
            map(iproc) = count
            do ishpair = first, last
               u = shpairs(1,ishpair)
               v = shpairs(2,ishpair)
               count = count + shdim(u)*shdim(v)
            end do
            first = last + 1
         end do
      end if
c     
      if (util_print('mp2_backt', print_debug)) then
         write(6,*) me, ' nshpair ', nshpair
         write(6,*) me, ' nshpairlocal ', nshpairlocal
         write(6,*) me, ' nbfpair ', nbfpair
         write(6,*) me, ' nbfpairlocal ', nbfpairlocal
         if (omakearrays) then
            do u = 1, nsh
               write(6,*) me, ' u shlo shdim ', u, shlo(u), shdim(u)
            end do
            do u = 1, nshpairlocal
               write(6,*) me, ' pair shpairslocal ', u,
     $              shpairslocal(1,u), shpairslocal(2,u)
            end do
            do u = 0, nproc-1
               write(6,*) me, ' proc map ', u, map(u)
            end do
         end if
      end if
c     
      end
      subroutine mp2_nonsep(
     $     rtdb,
     $     basis,
     $     twopdmunit,
     $     c_t,
     $     nsh, nbf, ninseg,
     $     nshpairlocal, shpairslocal, shdim, shlo,
     $     t, tbuf, grad,
     $     tol2e, oskel, oactive, leneri, eri, labels, lenscr, 
     $     scratch, shmap, cemap, pdm)
      implicit none
#include "errquit.fh"
#include "schwarz.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "bas.fh"
#include "geom.fh"
#include "global.fh"
#include "eaf.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "stdio.fh"
c     
      integer rtdb              ! [input]
      integer basis             ! [input]
      integer twopdmunit        ! [input]
      integer nsh, nbf, ninseg  ! [input] 
      double precision c_t(ninseg,nbf) ! [input] Transposed MOs in segment
      integer nshpairlocal, shpairslocal(2,*), shdim(nsh), shlo(nsh) ! [input]
      double precision t(*)     ! [scratch] ninseg*nbf*S*S
      double precision tbuf(*)  ! [scratch] nbf*S*S
      double precision grad(3,*) ! [input/output]
      double precision tol2e    ! [input]
      logical oskel             ! [input]
      logical oactive(*)        ! [input] oactive(ish)=true if shell is active
      integer leneri, lenscr    ! [input]
      integer labels(leneri,4), nq, nint ! [scratch]
      double precision eri(3,4,leneri), scratch(lenscr) ! [scratch]
      integer shmap(*), cemap(*) ! [scratch]
      double precision pdm(*)   ! [scratch]
c     
      double precision fileptr, energy, scale, psum, q4, block_eff
      integer i, count, ishpair
      integer ush, vsh, xsh, ysh, u, v, x, y
      integer udim, vdim, xdim, ydim, xlo, ylo
      integer ush_cur, vsh_cur, xsh_cur, ysh_cur
      integer ush_prev, vsh_prev, xsh_prev, ysh_prev
      integer xcent, ycent, ucent, vcent, ijkl, xsh_start, ysh_start
      logical omore
c
      integer maxq
      parameter (maxq=10000)
      integer sh_list(maxq,4)
      double precision q4_list(maxq)
c
      logical intbd_init4c, intbd_2e4c
      external intbd_init4c, intbd_2e4c
c     
      integer nat, geom,
     $     ulo_cur, uhi_cur, vlo_cur, vhi_cur, xlo_cur, xhi_cur, 
     $     ylo_cur, yhi_cur, udim_cur, vdim_cur, xdim_cur, ydim_cur
      logical status, odebug, odoit, sym_shell_quartet, oenergy
      external sym_shell_quartet
      double precision ddot
      external ddot
c     
      odebug = util_print('mp2_backt', print_debug)
      oenergy = util_print('backtenergy', print_debug)
c     
      energy = 0.0d0
c     
      status = bas_geom ( basis, geom )
      status = geom_ncent ( geom, nat )
c     
      do u = 1, nbf
         if (.not. bas_bf2ce(basis,u,ucent))
     $        call errquit('mp2g: bad something?',0, BASIS_ERR)
         if (.not. bas_bf2cn(basis,u,ush))
     $        call errquit('mp2g: bad something?',0, BASIS_ERR)
         cemap(u) = ucent
         shmap(u) = ush
      enddo
c     
      if (odebug) then
         write(6,*) ' Transposed occupied MOS in segment '
         call output(c_t, 1, ninseg, 1, nbf, ninseg, nbf, 1)
      end if
c     
      fileptr = 0.0d0
      do ishpair = 1, nshpairlocal
         ush = shpairslocal(1,ishpair)
         vsh = shpairslocal(2,ishpair)
         udim = shdim(ush)
         vdim = shdim(vsh)
         count = nbf*udim*vdim
c     
c     Part transformed density is stored as 
c     t(1:nbf,1:vdim,1:udim,1:ninseg,1:nshpairlocal) -> 
c     t(x,v,u,i,ush,vsh)
c     
c     Read this in and transpose to in core structure t(i,x,v,u)
c     
        if (eaf_read(twopdmunit, fileptr, tbuf, ninseg*count*8).ne.0)
     $       call errquit('mp2: ao test: failed reading density', 0,
     &     DISK_ERR)
         do i = 1, ninseg
            call dcopy(count, tbuf(count*(i-1)+1), 1, t(i), ninseg)
         end do
            fileptr = fileptr + ninseg*count*8.0d0
c     
         xsh_start = 1          ! For braindead multipassing
         ysh_start = 1
c
 333     nq = 0
         do xsh = xsh_start, nsh
            do ysh = ysh_start, xsh
               odoit = schwarz_shell(xsh,ysh)*schwarz_shell(ush,vsh)
     $              .gt. tol2e*0.1d0
               odoit = odoit .and. 
     $              (oactive(ush).or.oactive(vsh).or.
     $              oactive(xsh).or.oactive(ysh))
               q4 = 1.0d0
               if (odoit .and. oskel) odoit = 
     $              sym_shell_quartet(basis, ush, vsh, xsh, ysh, q4)
c     
               odoit = .true.
               if (odoit) then
c     
                  scale = q4
                  if (xsh .eq. ysh) scale = scale*0.5d0
                  xdim = shdim(xsh)
                  ydim = shdim(ysh)
                  xlo  = shlo(xsh)
                  ylo  = shlo(ysh)

                  call mp2_make_two_particle_density(psum,
     $                 udim, vdim, xdim, ydim, xlo, ylo, 
     $                 ninseg, nbf, 
     $                 scale, c_t, t, pdm)
c     
c     use density to screen evaluation of gradient integrals
                     if (oenergy) then
                        call int_2e4c(basis, ush, vsh, basis, xsh, ysh, 
     $                       lenscr, scratch, leneri, eri)
                        energy = energy + 
     $                       ddot(udim*vdim*xdim*ydim,eri,1,pdm,1)
                     endif
c
                  if (schwarz_shell(xsh,ysh)*schwarz_shell(ush,vsh)*psum
     $                 .gt. tol2e) then
c
                     nq = nq + 1
                     sh_list(nq,1) = ysh
                     sh_list(nq,2) = xsh
                     sh_list(nq,3) = vsh
                     sh_list(nq,4) = ush
                     q4_list(nq)   = scale
                     if (nq .eq. maxq) goto 666 ! UGLY 
                  endif         ! schwarz and density
               endif            ! schwarz
            enddo               ! ysh
            ysh_start = 1
         enddo                  ! xsh
c
 666     xsh_start = xsh
         ysh_start = ysh+1
c     
         omore = .false.
c     
         if (.not. intbd_init4c(
     $        basis, sh_list(1,1), sh_list(1,2),
     $        basis, sh_list(1,3), sh_list(1,4),
     $        nq, q4_list, .true., lenscr, scratch, leneri, 
     $        block_eff)) call errquit('mp2:backt: txs init?',nq,
     &       INT_ERR)
c     
 1000    omore = intbd_2e4c(
     $        basis, sh_list(1,1), sh_list(1,2),
     $        basis, sh_list(1,3), sh_list(1,4),
     $        nq, q4_list, .true., tol2e, .false.,
     $        labels(1,1),labels(1,2), 
     $        labels(1,3), labels(1,4), 
     $        eri, leneri, nint, lenscr, scratch)

         if (nint .gt. leneri) call errquit('mp2_nonsep: nint', nint,
     &       INT_ERR)
c     
         ush_prev = -1
         vsh_prev = -1
         ysh_prev = -1
         xsh_prev = -1
         do ijkl = 1, nint
            u = labels(ijkl,4)
            v = labels(ijkl,3)
            x = labels(ijkl,2)
            y = labels(ijkl,1)
            ush_cur = shmap(u)
            vsh_cur = shmap(v)
            xsh_cur = shmap(x)
            ysh_cur = shmap(y)
c
            if ( ush_cur.ne.ush_prev .or.
     $           vsh_cur.ne.vsh_prev .or.
     $           xsh_cur.ne.xsh_prev .or.
     $           ysh_cur.ne.ysh_prev ) then
c     
               ush_prev = ush_cur
               vsh_prev = vsh_cur
               xsh_prev = xsh_cur
               ysh_prev = ysh_cur
c     
               ulo_cur  = shlo(ush_cur)
               vlo_cur  = shlo(vsh_cur)
               xlo_cur  = shlo(xsh_cur)
               ylo_cur  = shlo(ysh_cur)
c     
               udim_cur = shdim(ush_cur)
               vdim_cur = shdim(vsh_cur)
               xdim_cur = shdim(xsh_cur)
               ydim_cur = shdim(ysh_cur)
c     
               uhi_cur  = ulo_cur + udim_cur - 1
               vhi_cur  = vlo_cur + vdim_cur - 1
               xhi_cur  = xlo_cur + xdim_cur - 1
               yhi_cur  = ylo_cur + ydim_cur - 1
c     
               ucent = cemap(u)
               vcent = cemap(v)
               xcent = cemap(x)
               ycent = cemap(y)
c     
               call mp2_make_two_particle_density( psum,
     $              udim_cur, vdim_cur, xdim_cur, ydim_cur, 
     $              xlo_cur, ylo_cur, 
     $              ninseg, nbf, 
     $              1.0d0, c_t, t, pdm)
            endif
c     
            call make_mp2grad(ucent,vcent,xcent,ycent, u,v,x,y,
     $           ulo_cur, uhi_cur, vlo_cur, vhi_cur,
     $           xlo_cur, xhi_cur, ylo_cur, yhi_cur,
     $           pdm, eri(1,1,ijkl), grad )
c     
         enddo
c     
         if (omore) goto 1000   ! Texas split the request
c
         if (xsh_start .le. nsh) goto 333
c     
      end do
c     
      if (oenergy) then
         call ga_dgop(1,energy,1,'+')
         if (ga_nodeid().eq.0) 
     +     write(LuOut,*) ' The energy from the two PDM is ', energy
      end if
c     
      end
      subroutine mp2_make_two_particle_density(
     $     psum, udim, vdim, xdim, ydim, xlo, ylo, ninseg, nbf, 
     $     scale, c_t, t, density)
      implicit none
c     
      double precision psum     ! [output] Norm of density
      integer udim, vdim, xdim, ydim ! [input] Shell dimensions
      integer xlo, ylo          ! [input] Start of shells
      integer ninseg            ! [input] No. of i in the batch
      integer nbf               ! [input]
      double precision scale    ! [input] Q4 and other factors
      double precision c_t(ninseg,nbf) ! [input] Transposed MO vectors in batch
      double precision t(ninseg,nbf,vdim,udim) ! [input] Partially transformed amplitudes
      double precision density(ydim,xdim,vdim,udim) ! [output] 
c
      integer u, v, x, y, i
      double precision p
c
      psum = 0.0d0
      do u = 1, udim
         do v = 1, vdim
            do x = 1, xdim
               do y = 1, ydim
                  p = 0.0d0
                  do i = 1, ninseg
                     p = p + t(i,x+xlo-1,v,u)*c_t(i,y+ylo-1)
     $                     + t(i,y+ylo-1,v,u)*c_t(i,x+xlo-1)
                  end do
                  p = p*scale
                  psum = psum + p*p
                  density(y,x,v,u) = p
               end do
            end do
         end do
      end do
c
      psum = sqrt(psum)
c
      end
      subroutine make_mp2grad(ucent,vcent,xcent,ycent, u,v,x,y,
     $                        ulo, uhi, vlo, vhi, xlo, xhi, ylo, yhi,
     $                        density, eri, grad )
      implicit none
      integer ucent,vcent,xcent,ycent, u,v,x,y 
      integer ulo, uhi, vlo, vhi, xlo, xhi, ylo, yhi 
      integer icart
      double precision density(ylo:yhi,xlo:xhi,vlo:vhi,ulo:uhi)
      double precision eri(3,4) , grad(3,*)
c
c ycent - eri(*,1)
c xcent - eri(*,2)
c vcent - eri(*,3)
c ucent - eri(*,4)
c
      do icart=1,3
         grad(icart,ucent)=grad(icart,ucent)+
     $                     density(y,x,v,u)*eri(icart,4)
         grad(icart,vcent)=grad(icart,vcent)+
     $                     density(y,x,v,u)*eri(icart,3)
         grad(icart,xcent)=grad(icart,xcent)+
     $                     density(y,x,v,u)*eri(icart,2)
         grad(icart,ycent)=grad(icart,ycent)+
     $                     density(y,x,v,u)*eri(icart,1)
      enddo
c
      end
c--------------------------------------------------
      subroutine mp2_copyback(g_a,g_a_trans)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
c
c     re-arrange (a la transponse) from one GA distrib. to the opposite
c
      integer g_a,g_a_trans
c
      integer dim1,dim2,type
      integer ilo,ihi,jlo,jhi
      integer i_length,j_length
      integer k_buff,l_buff
      integer i_chunk
      logical dowork
c
      call ga_distribution(g_a, ga_nodeid(), ilo, ihi, jlo, jhi)
      dowork=.true.
      if ((ilo.eq.0).and.(ihi.eq.-1))  dowork=.false.
      if(dowork) then
         i_length=ihi-ilo+1
         j_length=jhi-jlo+1
         if (.not. ma_push_get(mt_dbl,i_length*j_length,
     M        'mp2 copyb',l_buff,k_buff))
     $        call errquit('mp2:copyb could not allocate buff',
     L        i_length*j_length,MA_ERR)
      endif
      call mp2_copyback0(dowork,g_a,g_a_trans,
     C     dbl_mb(k_buff),ilo,ihi,jlo,jhi)
      if(dowork) then
         if (.not. ma_chop_stack(l_buff)) call errquit
     $        ('mp2: backt: failed chopping stack',0, MA_ERR)
      endif

      return
      end
      subroutine mp2_copyback0(dowork,g_a,g_a_trans,
     B     buff,ilo,ihi,jlo,jhi)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
      logical dowork
      integer g_a,g_a_trans
      integer ilo,ihi,jlo,jhi
      double precision buff(ilo:ihi,jlo:jhi)
c
      integer dim1,dim2,type
      integer i_length,j_length
      integer i_chunk,iseg,i_ld,i_0,i_1
      integer t_ilo, t_ihi, t_jlo, t_jhi
      integer nsegs,iseg_in
      double precision tr,tr_tr
      double precision ddot
      external ddot
c
      call ga_distribution(g_a_trans, ga_nodeid(), 
     T     t_ilo, t_ihi, t_jlo, t_jhi)
      i_chunk=t_ihi-t_ilo+1
 123  format(i4,a,2i9)
      call ga_igop(123,i_chunk,1,'max')
cmin      call ga_igop(123,i_chunk,1,'min')
      if(i_chunk.lt.1) then
         write(6,123) ga_nodeid(),'1trans ilo ihi ',t_ilo, t_ihi
         call errquit(' mp2_copyback: bad i_chunk',i_chunk,0)
      endif

      if(dowork) then
         call ga_inquire(g_a, type, dim1, dim2)
c     scatter first value
c     check for remainder
         nsegs=dim1/i_chunk
         if(mod(dim1,i_chunk).ne.0) nsegs=nsegs+1
cstaggered
         do iseg_in=ga_nodeid()+1,ga_nodeid()+nsegs
            iseg=mod(iseg_in,nsegs)+1
c            i_ld=i_chunk
            i_ld=ihi-ilo+1
c     check for remainder
            i_0=(iseg-1)*i_chunk+1
            i_1=iseg*i_chunk
            if(iseg.eq.nsegs) i_1=dim1
c     check for boundaries of buff()
            if((i_1-i_0).gt.(ihi-ilo)) then
               write(6,*) ga_nodeid(),' i_1.gt.ihi ',i_1,ihi
               call errquit(' mp2copybck: i boundary passed ',0,0)
            endif
            call ga_get(g_a_trans,i_0,i_1,
     A           jlo,jhi,buff,i_ld)
            call ga_put(g_a, i_0, i_1, jlo, jhi, buff, i_ld)
         enddo
      endif
c     everybody needs to have the buffer filled
      call ga_sync()
      return
      end
