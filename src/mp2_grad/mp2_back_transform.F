      subroutine mp2_nonsep_uhf( rtdb, geom,
     $     basis,
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
     $     tunita, tunitb, frc)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "bas.fh"
      integer basis, rtdb, geom
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
      double precision frc(3,*)
c
c     Allocate memory for back transformation routine
c
      integer g_buf
      integer nshpair, nshpairlocal, nbfpair, nbfpairlocal
      integer l_shpairs, k_shpairs, l_shpairslocal, k_shpairslocal
      integer l_shdim, k_shdim, l_shlo, k_shlo
      integer l_t, k_t, l_tmp, k_tmp, l_iauv, k_iauv, l_map, k_map
      integer l_c_t, k_c_t
      integer twopdmunit
      integer junk, ninseg, ierr, i, shmax, tdim
      integer nsh
      character*(nw_max_path_len) fname
c
      double precision tol2e
c
      tol2e = 1d-12
      if (.not. bas_numcont(basis, nsh)) call errquit
     $     ('mp2: backt bad basis handle ', basis)
      ninseg = oseg_hi - oseg_lo + 1
c
c     Make mapping of shell pairs to processors
c
      call intd_init(rtdb,1,basis)
      call schwarz_init(geom,basis)

      call mp2_backt_info(basis, tol2e, .false., nshpair, nshpairlocal,
     $     nbfpair, nbfpairlocal, junk, junk, junk, junk, junk)
      if (.not. ma_push_get(mt_int, nshpair*2, 'mp2: shpairs',
     $     l_shpairs, k_shpairs)) call errquit
     $     ('mp2: insufficient memory : shpairs ', 2*nshpair)
      if (.not. ma_push_get(mt_int, nshpairlocal*2, 'mp2:shpairslocal',
     $     l_shpairslocal, k_shpairslocal)) call errquit
     $     ('mp2: insufficient memory : shpairs ', 2*nshpairlocal)
      if (.not. ma_push_get(mt_int, ga_nnodes(), 'mp2:map',
     $     l_map, k_map)) call errquit
     $     ('mp2: insufficient memory : map ', ga_nnodes())
      if (.not. ma_push_get(mt_int, nsh, 'mp2:shlo',
     $     l_shlo, k_shlo)) call errquit
     $     ('mp2: insufficient memory : shlo ', nsh)
      if (.not. ma_push_get(mt_int, nsh, 'mp2:shdim',
     $     l_shdim, k_shdim)) call errquit
     $     ('mp2: insufficient memory : shdim ', nsh)
      call mp2_backt_info(basis, tol2e, .true., nshpair, nshpairlocal,
     $     nbfpair, nbfpairlocal, int_mb(k_shdim), int_mb(k_shlo),
     $     int_mb(k_shpairs), int_mb(k_shpairslocal), int_mb(k_map))
c
c     Allocate remaining memory and open file to hold the 
c     3-parts transformed density matrix.
c
      shmax = 0
      do i = 1, nsh
         shmax = max(shmax,int_mb(k_shdim+i-1))
      enddo
      tdim = max(nbf*nbf, nbf*shmax*shmax)
c         
      if (.not. ma_push_get(mt_dbl, tdim, 'mp2: backt t',
     $     l_t, k_t)) call errquit
     $     ('mp2: insufficient memory : t ', tdim)
      if (.not. ma_push_get(mt_dbl, tdim, 'mp2: backt tmp',
     $     l_tmp, k_tmp)) call errquit
     $     ('mp2: insufficient memory : tmp ', tdim)
      if (.not. ma_push_get(mt_dbl, nbf*nbf, 'mp2: backt iauv',
     $     l_iauv, k_iauv)) call errquit
     $     ('mp2: insufficient memory : iauv ', nbf*nbf)
c
c     Right now we have an N**3 array in GA.  This is easily
c     replaced with an N*B array with B=P*S*S (P procs, S=max shell size)
c     by blocking the transformation.
c
      if (.not. ga_create_irreg(mt_dbl, nva, nbfpair,
     $     'mp2: backt', 1, 1, int_mb(k_map), ga_nnodes(), g_buf))
     $     call errquit('mp2: backt: ga_create failed', nva*nbfpair)
c     
      twopdmunit = 33
      call util_file_name('2pdm', .true., .true., fname)
      call waopen(twopdmunit,fname,1,0,ierr)
      if (ierr .ne. 0) call errquit('mp2: backt: failed to open file',0)

c     
c     Now do the back transformation
c
      call mp2_back_transform_uhf(
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
     $     g_buf)
c
c     Free up some memory and then contract with gradient integrals.
c
      if (.not. ga_destroy(g_buf)) call errquit
     $     ('mp2: backt: ga destroy failed', 0)
      if (.not. ma_chop_stack(l_t)) call errquit
     $     ('mp2:backt: first chop stack failed', 0)
c
      tdim = shmax*shmax*nbf*(oseg_hi-oseg_lo+1)
      if (.not. ma_push_get(mt_dbl, tdim, 'mp2:back tt',
     $     l_t, k_t)) call errquit
     $     ('mp2: backt: failed ma for nonsep test', tdim)
      tdim = shmax*shmax*nbf
      if (.not. ma_push_get(mt_dbl, tdim, 'mp2:back tbuf',
     $     l_tmp, k_tmp)) call errquit
     $     ('mp2: backt: failed ma for nonsep test', tdim)
      if (.not. ma_push_get(mt_dbl, ninseg*nbf, 'mp2:back c_t',
     $     l_c_t, k_c_t)) call errquit
     $     ('mp2: backt: failed ma for nonsep test', ninseg*nbf)
      do i = oseg_lo, oseg_hi
         call dcopy(nbf, c_a(1,i), 1, dbl_mb(k_c_t+i-oseg_lo), ninseg)
      enddo
c      
      call mp2_nonsep_test(basis,
     $     twopdmunit,
     $     dbl_mb(k_c_t),
     $     nsh, nbf, ninseg,
     $     nshpairlocal, int_mb(k_shpairslocal),
     $     int_mb(k_shdim), int_mb(k_shlo),
     $     dbl_mb(k_t), dbl_mb(k_tmp), frc,
     $     tol2e)

      call intd_terminate()
      call schwarz_tidy()
c
c     Done.
c     
      if (.not. ma_chop_stack(l_shpairs)) call errquit
     $     ('mp2: backt: failed chopping stack',0)
c
      call waclose(twopdmunit, ierr)
      if (ierr .ne. 0) call errquit('mp2: backt: closing 2pdm',ierr)
      call util_file_unlink(fname)
c
      end
      subroutine mp2_back_transform_uhf(
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
     $     g_buf)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
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
      integer g_buf             ! Global array handle for transpose buffer
c
      double precision zero, one, scale
      integer tunitptra, tunitptrb, twopdmunitptr  ! Pointers into files
      integer count, tcount, ind, ierr, ptr
      integer symi, symj, syma, symb, symia
      integer i, a, b, u, v, ush, vsh, udim, vdim, ishpair
c     
      zero=0.0d0
      one=1.0d0
c
      tunitptra=1
      tunitptrb=1
c
      do i=oseg_lo,oseg_hi 
	 call ga_fill_patch(g_buf,1,nva,1,nbfpair,99.0d0) ! For debug
         symi=irs_a(i)
         do a=nva_lo_local,nva_hi_local
            syma=irs_a(a)
            symia=ieor(syma,symi)
            call dfill((nbf*nbf),zero,ia_uv,1) ! Will add pure and mixed in
c     
            tcount=1            ! Read pure spin t(ia,**) from disk 
            do b=nva_lo,nva_hi
               symb=irs_a(b)
               symj=ieor(symb,symia)
               if(num_oa(symj).ne.0)then
                  count=sym_hi_oa(symj)-sym_lo_oa(symj)+1
                  call getwa(tunita,t(tcount),
     $                 tunitptra,count,ierr)
                  if(ierr.ne.0)call errquit('mp2_back: reading ta',ierr)
                  tunitptra=tunitptra+count
                  tcount=tcount+count
               endif
            enddo
c     
            tcount=1
            do symb=0,nir-1
               symj=ieor(symia,symb)
               if(num_va(symb).gt.0.and.num_oa(symj).gt.0) then
                  call dgemm('n','t', !             t(j,b)Cbv -> t(j,v)
     $                 num_oa(symj), nbf, num_va(symb),
     $                 one, t(tcount), num_oa(symj),
     $                 c_a(1,sym_lo_va(symb)), nbf,
     $                 zero, tmp, num_oa(symj))
                  call dgemm('n','n', !             Cuj t(j,v) -> t(u,v)
     $                 nbf, nbf, num_oa(symj),
     $                 one, c_a(1,sym_lo_oa(symj)), nbf,
     $                 tmp, num_oa(symj),
     $                 one, ia_uv, nbf)
                  tcount=tcount+num_oa(symj)*num_va(symb)
               endif
            enddo
c
            tcount=1           ! Read mixed spin t(ia,**) from disk 
            do b=nvb_lo,nvb_hi
               symb=irs_b(b)
               symj=ieor(symb,symia)
               if(num_ob(symj).ne.0)then
                  count=sym_hi_ob(symj)-sym_lo_ob(symj)+1
                  call getwa(tunitb,t(tcount),
     $                 tunitptrb,count,ierr)
                  if(ierr.ne.0)call errquit('mp2_back: reading tb',ierr)
                  tunitptrb=tunitptrb+count
                  tcount=tcount+count
               endif
            enddo
c
            tcount=1
            do symb=0,nir-1
               symj=ieor(symia,symb)
               if(num_vb(symb).gt.0.and.num_ob(symj).gt.0) then
                  call dgemm('n','t', !             t(j,b)Cbv -> t(j,v)
     $                 num_ob(symj), nbf, num_vb(symb),
     $                 one, t(tcount), num_ob(symj),
     $                 c_b(1,sym_lo_vb(symb)), nbf,
     $                 zero, tmp, num_ob(symj))
                  call dgemm('n','n', !             Cuj t(j,v) -> t(u,v)
     $                 nbf, nbf, num_ob(symj),
     $                 one, c_b(1,sym_lo_ob(symj)), nbf,
     $                 tmp, num_ob(symj),
     $                 one, ia_uv, nbf)
                  tcount=tcount+num_ob(symj)*num_vb(symb)
               endif
            enddo
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
                  enddo
               enddo
            enddo
            if (ind .ne. nbfpair) call errquit('mp2bt: ind?', ind)
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
            call ga_put(g_buf, a-nva_lo+1, a-nva_lo+1,
     $           1, nbfpair, tmp, 1)
c
         enddo                  ! End of a
         call ga_sync
c
c     Now have locally (a,1:nbfpairlocal).  Loop thru local shell pairs
c     transform the a index and write to disk.
c
         ind = 1
         ptr = 1
         do ishpair = 1, nshpairlocal
            ush = shpairs(1,ishpair)
            vsh = shpairs(2,ishpair)
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
            twopdmunitptr = ptr + (i-oseg_lo)*count
*            write(6,*) ' i = ', i, ' ishpair = ', ishpair
*            call output(t, 1, nbf, 1, udim*vdim, nbf, udim*vdim, 1)
c
            call putwa(twopdmunit, t, twopdmunitptr, count, ierr)
            if (ierr.ne.0) call errquit
     $           ('mp2_bt: write of two pdm?', ierr)
c
            ptr = ptr + count*(oseg_hi-oseg_lo+1)
         enddo
c
         call ga_sync
c            
      enddo                     ! End of i
c
      end
      subroutine mp2_backt_info(
     $     basis,
     $     tol2e,
     $     omakearrays,
     $     nshpair, nshpairlocal, 
     $     nbfpair, nbfpairlocal,
     $     shdim, shlo, shpairs, shpairslocal, map)
      implicit none
#include "bas.fh"
#include "schwarz.fh"
#include "global.fh"
#include "util.fh"
      integer basis             ! [input] basis set handle
      double precision tol2e    ! [input] screening threshold
      logical omakearrays       ! [input] if true then make the arrays
      integer nshpair           ! [output] no. of non-zero shell pairs
      integer nshpairlocal      ! [output] no. of local shell pairs
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
c        actually make the arrays
c     else
c        just return the scalar results
c     endif
c
c     The Schwarz package must be initialized before entry
c
c     Symmetry needs to be added.
c     
      integer u, v, nbf, nsh, lo, hi, udim, vdim, me, nproc
c
      if (.not. bas_numbf(basis, nbf))
     $     call errquit('mp2_backt_info: bad basis handle', 0)
      if (.not. bas_numcont(basis, nsh))
     $     call errquit('mp2_backt_info: bad basis handle', 0)
c
      nshpair = 0
      nshpairlocal = 0
      nbfpair = 0
      nbfpairlocal = 0
      me = ga_nodeid()
      nproc = ga_nnodes()
      if (omakearrays) call ifill(nproc, 0, map, 1)
c
      do u = 1, nsh
         if (.not. bas_cn2bfr(basis, u, lo, hi))
     $        call errquit('mp2_backt_info: bas range failed', u)
         if (omakearrays) then
            shdim(u) = hi - lo + 1
            shlo(u) = lo
         endif
         udim = hi - lo + 1
         do v = 1, u
            if (.not. bas_cn2bfr(basis, v, lo, hi))
     $           call errquit('mp2_backt_info: bas range failed', v)
            vdim = hi - lo + 1
            if (schwarz_shell(u,v)*schwarz_max().gt.tol2e) then
               if (mod(nshpair,nproc).eq.me) then
                  nshpairlocal = nshpairlocal + 1
                  nbfpairlocal = nbfpairlocal + udim*vdim
                  if (omakearrays) then
                     map(me) = map(me) + udim*vdim
                     shpairslocal(1,nshpairlocal) = u
                     shpairslocal(2,nshpairlocal) = v
                  endif
               endif
               nshpair = nshpair + 1
               nbfpair = nbfpair + udim*vdim
               if (omakearrays) then
                  shpairs(1,nshpair) = u
                  shpairs(2,nshpair) = v
               endif
            endif
         enddo
      enddo
c
      if (omakearrays) then
         do u = 1, nproc-1
            map(u) = map(u) + map(u-1)
         enddo
         do u = 1, nproc-1
            map(u) = map(u-1)
         enddo
         map(0) = 1
      endif
c
      if (util_print('mp2_backt', print_debug)) then
         write(6,*) ' nshpair ', nshpair
         write(6,*) ' nshpairlocal ', nshpairlocal
         write(6,*) ' nbfpair ', nbfpair
         write(6,*) ' nbfpairlocal ', nbfpairlocal
         if (omakearrays) then
            do u = 1, nsh
               write(6,*) ' u shlo shdim ', u, shlo(u), shdim(u)
            enddo
            do u = 1, nshpairlocal
               write(6,*) ' pair shpairslocal ', u,
     $              shpairslocal(1,u), shpairslocal(2,u)
            enddo
            do u = 0, nproc-1
               write(6,*) ' proc map ', u, map(u)
            enddo
         endif
      endif
c
      end
      subroutine mp2_nonsep_test(
     $     basis,
     $     twopdmunit,
     $     c_t,
     $     nsh, nbf, ninseg,
     $     nshpairlocal, shpairslocal, shdim, shlo,
     $     t, tbuf, frc,
     $     tol2e)
      implicit none
#include "schwarz.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "bas.fh"
#include "geom.fh"
c     
      integer basis
      integer twopdmunit
      integer nsh, nbf, ninseg
      double precision c_t(ninseg,nbf) ! Transposed MOs in segment
      integer nshpairlocal, shpairslocal(2,*), shdim(nsh), shlo(nsh)
      double precision t(*)     ! ninseg*nbf*S*S
      double precision tbuf(*)  ! nbf*S*S
      double precision tol2e, p

      double precision frc
      dimension frc(3,*)
C     
c     
      integer fileptr, ierr, i, leng, lenscr, count
      integer ishpair
      integer ush, vsh, xsh, ysh, u, v, x, y
      integer udim, vdim, xdim, ydim, xlo, ylo
      double precision energy, scale
      integer l_g, k_g, l_scr, k_scr
      integer uvxy, vuoff, xvuoff, yvuoff
c     
      double precision g
c     
      integer nat, idatom, lengd, lenscrd, l_gd, k_gd, iat, icart, 
     $     k_scrd, l_scrd, idum, geom, idum2, myat
      dimension idatom(4)
      double precision dE
      logical status
c     
      energy = 0.0d0
c     
      call int_mem_2e4c(leng, lenscr)
      call int_mem ( idum, lengd, idum2, lenscrd )
      lenscr = max(lenscr, leng)
      lenscrd = max(lenscrd, lengd)
      status = bas_geom ( basis, geom )
      status = geom_ncent ( geom, nat )
      if (.not. ma_push_get(mt_dbl, lenscr, 'mp2:test', l_scr, k_scr))
     $     call errquit('mp2:test ma', lenscr)
      if (.not. ma_push_get(mt_dbl, leng, 'mp2:test', l_g, k_g))
     $     call errquit('mp2:test ma', leng)
      if (.not. ma_push_get(mt_dbl, lenscrd, 'mp2:test', l_scrd,k_scrd))
     $     call errquit('mp2:test ma', lenscrd)
      if (.not. ma_push_get(mt_dbl, lengd, 'mp2:test', l_gd, k_gd))
     $     call errquit('mp2:test ma', lengd)
c     
      if (util_print('mp2_backt', print_debug)) then
         write(6,*) ' Transposed occupied MOS in segment '
         call output(c_t, 1, ninseg, 1, nbf, ninseg, nbf, 1)
      endif
c     
      fileptr = 1
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
         do i = 1, ninseg
            call getwa(twopdmunit, tbuf, fileptr, count, ierr)
            if (ierr .ne. 0) call errquit
     $           ('mp2: ao test: failed reading density', 0)
            fileptr = fileptr + count
            call dcopy(count, tbuf, 1, t(i), ninseg)
         enddo
c     
         do xsh = 1, nsh
            do ysh = 1, xsh
               if (schwarz_shell(xsh,ysh)*schwarz_shell(ush,vsh)
     $              .gt. tol2e) then
                  scale = 1.0d0
                  if (xsh .eq. ysh) scale = 0.5d0
                  xdim = shdim(xsh)
                  ydim = shdim(ysh)
                  xlo  = shlo(xsh)
                  ylo  = shlo(ysh)

                  call int_2e4c(basis, ush, vsh, basis, xsh, ysh,
     $                 lenscr, dbl_mb(k_scr), leng, dbl_mb(k_g))

                  call intd_2e4c(basis, ush, vsh, basis, xsh, ysh,
     $                 lenscrd, dbl_mb(k_scrd), lengd, dbl_mb(k_gd), 
     $                 idatom)

c     first the energy

                  uvxy = k_g
                  do u = 1, udim
                     do v = 1, vdim
                        vuoff = nbf*(v-1 + vdim*(u-1))
                        do x = xlo, xlo+xdim-1
                           xvuoff = ninseg*(x-1 + vuoff)
                           do y = ylo, ylo+ydim-1
                              yvuoff = ninseg*(y-1 + vuoff)
                              g = dbl_mb(uvxy)
                              if (abs(g).gt.tol2e) then
                                 p = 0.0d0
                                 do i = 1, ninseg
                                    p = p +
     $                                   t(i+xvuoff)*c_t(i,y)+
     $                                   t(i+yvuoff)*c_t(i,x)
                                 enddo
                                 energy = energy + p * g * scale
                              endif
                              uvxy = uvxy + 1
                           enddo
                        enddo
                     enddo
                  enddo
C     and then the gradient
                  uvxy = k_gd
                  do iat = 1, 4
                     if ( idatom(iat) .gt. 0 ) then
                        do icart = 1, 3
                           dE = 0.D0
                           do u = 1, udim
                              do v = 1, vdim
                                 vuoff = nbf*(v-1 + vdim*(u-1))
                                 do x = xlo, xlo+xdim-1
                                    xvuoff = ninseg*(x-1 + vuoff)
                                    do y = ylo, ylo+ydim-1
                                       yvuoff = ninseg*(y-1 + vuoff)
                                       g = dbl_mb(uvxy)
                                       if (abs(g).gt.tol2e) then
                                          p = 0.0d0
                                          do i = 1, ninseg
                                             p = p +
     $                                            t(i+xvuoff)*c_t(i,y)+
     $                                            t(i+yvuoff)*c_t(i,x)
                                          enddo
                                          dE = dE + p * g * scale
                                       endif
                                       uvxy = uvxy + 1
                                    enddo
                                 enddo
                              enddo
                           enddo
                           frc(icart,idatom(iat)) =
     $                          frc(icart,idatom(iat)) + dE
                        enddo   ! icart
                     else
                        uvxy = uvxy + 3*udim*vdim*xdim*ydim
                     endif      ! idatom
                  enddo         ! iat
               endif
            enddo
         enddo
      enddo
c     
      call ga_dgop(1,frc,3*nat,'+')
      call ga_dgop(1,energy,1,'+')
      write(6,*) ' The energy from the two PDM is ', energy
      write(6,*) ' The gradient (xyz down column, atoms across rows)'
      call output(frc, 1, 3, 1, nat, 3, nat, 1)
c     
      if (.not. ma_pop_stack(l_gd))
     $     call errquit('mp2:test: ma pop failed', 0)
      if (.not. ma_pop_stack(l_scrd))
     $     call errquit('mp2:test: ma pop failed', 0)
      if (.not. ma_pop_stack(l_g))
     $     call errquit('mp2:test: ma pop failed', 0)
      if (.not. ma_pop_stack(l_scr))
     $     call errquit('mp2:test: ma pop failed', 0)
c     
      end
