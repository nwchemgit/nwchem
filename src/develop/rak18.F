*
* $Id$
*
C> \ingroup task
C> @{
      logical function task_computeSld(rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "bas.fh"
#include "util.fh"
#include "stdio.fh"
c::functions
      logical int_normalize
      external int_normalize
c::passed
      integer rtdb ! [input] RTDB handle
c::local
      integer geom     ! geometry handle
      integer basis    ! basis set handle
      integer nbf      ! number of basis functions
      integer nshells  ! number of contractions / shells
      integer max1e    ! maximum size of a 1-e integral buffer
      integer mscr     ! estimated size of scratch buffer for computing integrals
      integer g_S      ! GA handle to overlap array
      integer g_V      ! GA handle to eigenvectors
      integer nproc    ! number of processors
      integer mynode   ! lexical index of "my node"
      integer compute_node ! lexical counter for parallel tasks
      integer ish      ! loop index for i shells
      integer jsh      ! loop index for j shells
      integer i        ! loop indes for eigenvalues
      integer ndep     ! number of linear dependent eigenvalues
      integer h_buff, k_buff ! ma handle/index for integral buffer
      integer h_scr, k_scr ! ma handle/index for integral scratch
      integer h_eig, k_eig ! ma handle/index for overlap eigenvalues
      integer h_eigvec, k_eigvec ! ma handle/index for overlap eigenvectors local
      integer ilo, ihi, nbfi, jlo, jhi, nbfj ! bf ranges/sizes
      logical status  ! dummy status variable
      double precision sld_thresh ! linear dependence threshold
c

* get threshold from RTDB
      if (.not.rtdb_get(rtdb,'sld_thresh',mt_dbl,1,sld_thresh))
     &    sld_thresh = 1.0d-05
* create geometry handle
      if (.not.geom_create(geom,'geometry'))
     &      call errquit('task_computeSld: geom_create failed?',911,
     &       GEOM_ERR)
* load geometry from rtdb
      if (.not.geom_rtdb_load(rtdb,geom,'geometry'))
     &      call errquit
     &      ('task_computeSld: geom_rtdb_load failed?',911, RTDB_ERR)
* create basis set handle
      if (.not.bas_create(basis,'ao basis'))
     &      call errquit('task_computeSld: bas_create failed?',911,
     &       BASIS_ERR)
* load basis set
      if (.not.bas_rtdb_load(rtdb,geom,basis,'ao basis'))
     &      call errquit
     &      ('task_computeSld: bas_rtdb_load failed?',911, BASIS_ERR)
c
* query basis set object for the number of functions and shells      
c
      if (.not.bas_numbf(basis,nbf)) call errquit
     &      ('task_computeSld: bas_numbf failed?',911, BASIS_ERR)
      if (.not.bas_numcont(basis,nshells)) call errquit
     &      ('task_computeSld: bas_numcont failed?',911, BASIS_ERR)
*
*normalize basis set
      if (.not.int_normalize(rtdb,basis)) call errquit
     &      ('task_computeSld: int_normalize failed?',911, INT_ERR)
*
* initialize the integral API
      call int_init(rtdb,1,basis)
*
* query integral API for memory requirements.
      call int_mem_1e(max1e, mscr)
*
* create an evenly distributed global array for S
      if (.not.ga_create(mt_dbl,nbf,nbf,'overlap matrix',0,0,g_S))
     &      call errquit('task_computeSld: ga_create failed',911,
     &       GA_ERR)
      call ga_zero(g_S)         ! zero GA memory segment
*
* allocate local memory buffers for integrals and scratch space for 
* the integral API
      if (.not. ma_push_get(mt_dbl,max1e,'integral buffer',
     &      h_buff,k_buff)) call errquit
     &      ('task_computeSld: could not allocat integral buffer',911,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl,mscr,'integral scratch',
     &      h_scr,k_scr)) call errquit
     &      ('task_computeSld: could not allocat integral buffer',911,
     &       MA_ERR)
c start parallel computational loop
      nproc  = ga_nnodes()
      mynode = ga_nodeid()
      compute_node = (mynode-1)
      do ish = 1,nshells
        if (.not.bas_cn2bfr(basis,ish,ilo,ihi)) call errquit
     &      ('task_computeSld: bas_cn2bfr failed',911, BASIS_ERR)
        nbfi = ihi - ilo + 1
        do jsh = 1,nshells
          if (.not.bas_cn2bfr(basis,jsh,jlo,jhi)) call errquit
     &        ('task_computeSld: bas_cn2bfr failed',911, BASIS_ERR)
          nbfj = jhi - jlo + 1
          compute_node = compute_node + 1
          if (mod(compute_node,nproc).eq.0) then
*...compute it: overlap(ilo:ihi,jlo:jhi)
            call int_1eov(basis,jsh,basis,ish,
     &            mscr,dbl_mb(k_scr),
     &            max1e,dbl_mb(k_buff))
            call ga_put(g_S,ilo,ihi,jlo,jhi,dbl_mb(k_buff),nbfi)
          endif
        enddo
      enddo
c
* free local memory (pop stack) 
      status = ma_pop_stack(h_scr)
      status = status .and. ma_pop_stack(h_buff)
      if (.not.status) call errquit
     &          ('task_computeSld: poping stack failed',911, MA_ERR)
c
      call ga_sync()
      if (util_print('sld overlap',print_high)) then
        if (mynode.eq.0) then
          write(luout,*)' overlap matrix'
        endif
        call ga_print(g_S)
      endif
*
* create an evenly distributed global array for eigenvectors
      call ga_sync()
      if (.not.ga_create(mt_dbl,nbf,nbf,'S eigenvectors',0,0,g_V))
     &      call errquit('task_computeSld: ga_create failed',911,
     &       GA_ERR)
      call ga_zero(g_V)         ! zero GA memory segment
      call ga_sync()
*
* allocate memory for eigenvalues
      if (.not. ma_push_get(mt_dbl,nbf,'S eigenvalues',
     &      h_eig,k_eig)) call errquit
     &      ('task_computeSld: could not allocat integral buffer',911,
     &       MA_ERR)
*
* allocate memory for eigenvalues
      if (.not. ma_push_get(mt_dbl,nbf*nbf,'S eigenvectors',
     &      h_eigvec,k_eigvec)) call errquit
     &      ('task_computeSld: could not allocat integral buffer',911,
     &       MA_ERR)
      call dfill(nbf,0.0d00,dbl_mb(k_eig),1)
      call dfill(nbf,0.0d00,dbl_mb(k_eigvec),1)
c
      call ga_diag_std(g_S,g_V,dbl_mb(k_eig))
      call ga_get(g_V,1,nbf,1,nbf,dbl_mb(k_eigvec))
      ndep = 0
      do i = 0,(nbf-1)
        if (dbl_mb(k_eig+i).lt.sld_thresh) ndep = ndep + 1
      enddo
      if (mynode.eq.0) then
        if (ndep.gt.0) then
          write(luout,*)' ',ndep,
     &        ' eigenvalues in S are linear dependent'
          write(luout,*)' threshold:',sld_thresh
          write(luout,*)' '
        endif
        write(luout,*)' eigenvalues in S'
*        call output(dbl_mb(k_eig),1,nbf,1,1,nbf,1,1)
        call print_Sld(dbl_mb(k_eig),dbl_mb(k_eigvec),nbf,luout)
      endif
      call int_terminate()
      if (.not.bas_destroy(basis)) call errquit
     &          ('task_computeSld: bas_destroy failed',911, BASIS_ERR)
      if (.not.geom_destroy(geom)) call errquit
     &          ('task_computeSld: geom_destroy failed',911, GEOM_ERR)
* free memory
      status = ma_pop_stack(h_eigvec)
      status = ma_pop_stack(h_eig)
      if (.not.status) call errquit
     &          ('task_computeSld: poping stack failed',911, MA_ERR)
* free global arrays
      if (.not.ga_destroy(g_S)) call errquit
     &          ('task_computeSld: ga_destroy of g_S failed',911,
     &       GA_ERR)
      if (.not.ga_destroy(g_V)) call errquit
     &          ('task_computeSld: ga_destroy of g_V failed',911,
     &       GA_ERR)
      end
C> @}
      subroutine print_Sld(eigs,eigvec,n,lunit)
      implicit none
      integer n
      double precision eigs(n)
      double precision eigvec(n,n)
      integer lunit
c
      integer ivallo, ivalhi, ival, ivec, blocksize
      write(6,*)' n '
      blocksize = 6
      do ivallo = 1,n,blocksize
        ivalhi = min(n,(ivallo+blocksize-1))
        write(6,*)'ivallo ', ivallo
        write(6,*)'ivalhi ', ivalhi
        write(lunit,10000)(eigs(ival),ival=ivallo,ivalhi)
        do ivec = 1,n
          write(lunit,10001)ivec,(eigvec(ivec,ival),ival=ivallo,ivalhi)
        enddo
      enddo
10000 format(5x,6(1pd15.8,2x),/,1x,(6(15('-'),2x)))
10001 format(1x,i4,6(1pd15.8,2x))
      end
