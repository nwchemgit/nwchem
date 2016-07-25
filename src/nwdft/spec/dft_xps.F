c
c     Calculates xps intensity between the initial and final states
c
      subroutine dft_xps(rtdb,ao_bas_han,geom)
c
      implicit none
c
#include "bas.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "msgtypesf.h"
#include "tcgmsg.fh"
#include "util.fh"
#include "inp.fh"
c
      integer rtdb               ! [input]
      integer ao_bas_han         ! [input]
      integer geom
c
      integer ipol 
      integer nmo(2)
      integer nbf
c
      integer nbf_verify  
      integer ipol_verify  
c
      integer g_mos_ini(2)
      integer g_mos_fin(2)
c
      integer l_evl,k_evl      ! MA handles for eigenvalues 
      integer l_occ,k_occ      ! MA handles for occupancies 
c
      integer g_dpx            ! dipole x-component ao basis
      integer g_dpy            ! dipole y-component ao basis
      integer g_dpz            ! dipole z-component ao basis
c
      integer g_tdx(2)         ! transition dipole x
      integer g_tdy(2)         ! transition dipole y
      integer g_tdz(2)         ! transition dipole z
c
      integer g_tmp
c
      double precision tdx,tdy,tdz      ! transition moments
c
      logical nodezero
      integer i
c
      character*255 title      ! Title
      character*255 basisname  ! Basis set name
      character*255 mos_ini    ! initial state movecs
      character*255 mos_fin    ! final state movecs
      character*255 scftype    ! type of the scf: RHF OR UHF
      character*32 pname
      double precision xps_scale
      double precision xps_wif
      logical oskel
c
      double precision two, three
      parameter (two=2.d0,three=3.d0)
c
      logical movecs_read,movecs_read_header
      external movecs_read,movecs_read_header
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
      double precision ga_trace_diag
      external ga_trace_diag
c
c     Preliminaries
      nodezero=(ga_nodeid().eq.0)
      pname="dft_xps: "
c
c     Get oscillator strength scale factor
      xps_scale = 1.d0
      if (.not. rtdb_get(rtdb, 'xps:scale', mt_dbl, 1, xps_scale))
     &     xps_scale = 1.d0
c
c     Get scftype
      if (.not. rtdb_cget(rtdb, 'dft:scftype', 1, scftype))
     $  call errquit(pname//'DFT scftype not defined',0,RTDB_ERR)
      if (scftype.eq.'RHF') ipol = 1
      if (scftype.eq.'UHF') ipol = 2
c
c     Get the number of basis functions
      if (.not. bas_numbf(ao_bas_han,nbf)) 
     &    call errquit(pname//'failed to get nbf',1,BASIS_ERR)
c
c     Get the file names
      if (.not. rtdb_cget(rtdb,'dft:mos_ini',1,mos_ini)) 
     &  call errquit(pname//'failed to get initial vectors',0,RTDB_ERR)
      if (.not. rtdb_cget(rtdb,'dft:mos_fin',1,mos_fin)) 
     &  call errquit(pname//'failed to get final vectors',0,RTDB_ERR)
c
c      if (nodezero) then
c            write(luout,*) "dft:scftype:",scftype
c            write(luout,*) "dft:mos_ini:",mos_ini
c            write(luout,*) "dft:mos_fin:",mos_fin 
c      end if 
c
c     Allocate global arrays for MO vectors
      if (.not.ma_push_get(mt_dbl,nbf*ipol,'evl',l_evl,k_evl)) 
     &  call errquit(pname//'failed to allocate evl',0, MA_ERR)
      if (.not.ma_push_get(mt_dbl,nbf*ipol,'occ',l_occ,k_occ)) 
     &  call errquit(pname//'failed to allocate occ',0, MA_ERR)
c
c    Read movecs of the initial state
      if (.not.movecs_read_header(mos_ini,title,basisname,
     1 scftype,nbf_verify,ipol_verify,nmo,2)) 
     &  call errquit(pname//'failed to read MO file header',0,DISK_ERR)
      if (nbf.ne.nbf_verify) 
     &  call errquit(pname//'corrupted MO vectors',0,DISK_ERR)
      if (ipol.ne.ipol_verify) 
     &  call errquit(pname//'corrupted MO vectors',0,DISK_ERR)
      do i=1,ipol
        g_mos_ini(i)=ga_create_atom_blocked(geom,ao_bas_han,'ini vecs')
        if (.not.movecs_read(mos_ini,i,dbl_mb(k_occ+(i-1)*nbf),
     2    dbl_mb(k_evl+(i-1)*nbf),g_mos_ini(i))) 
     &    call errquit(pname//'failed to read ini vecs',0, DISK_ERR)
      enddo
c
c    Read movecs of the final state
      if (.not.movecs_read_header(mos_fin,title,basisname,
     1 scftype,nbf_verify,ipol_verify,nmo,2)) 
     &  call errquit(pname//'failed to read MO file header',0,DISK_ERR)
      if (nbf.ne.nbf_verify) 
     &  call errquit(pname//'corrupted MO vectors',0,DISK_ERR)
      if (ipol.ne.ipol_verify) 
     &  call errquit(pname//'corrupted MO vectors',0,DISK_ERR)
      do i=1,ipol
        g_mos_fin(i)=ga_create_atom_blocked(geom,ao_bas_han,'fin vecs')
        if (.not.movecs_read(mos_fin,i,dbl_mb(k_occ+(i-1)*nbf),
     2    dbl_mb(k_evl+(i-1)*nbf),g_mos_fin(i))) 
     &    call errquit(pname//'failed to read ini vecs',0, DISK_ERR)
      enddo
c
c     Allocate ga memory for the AO dipole & td components
      if (.not.ga_create(mt_dbl,nbf,nbf,'Dip-X',nbf,0,g_dpx)) 
     &   call errquit(pname//'error creating g_dpx',0,GA_ERR)
      if (.not.ga_create(mt_dbl,nbf,nbf,'Dip-Y',nbf,0,g_dpy)) 
     &   call errquit(pname//'error creating g_dpy',0,GA_ERR)
      if (.not.ga_create(mt_dbl,nbf,nbf,'Dip-Z',nbf,0,g_dpz)) 
     &   call errquit(pname//'error creating g_dpz',0,GA_ERR)
      if (.not.ga_create(mt_dbl,nbf,nbf,'tmp',nbf,0,g_tmp)) 
     &   call errquit(pname//'error creating g_tmp',0,GA_ERR)
c
      do i=1,ipol
       if (.not.ga_create(mt_dbl,nbf,nbf,'TD-X',nbf,0,g_tdx(i))) 
     &     call errquit(pname//'error creating g_tdx',0,GA_ERR)
           call ga_zero(g_tdx(i))
       if (.not.ga_create(mt_dbl,nbf,nbf,'TD-Y',nbf,0,g_tdy(i))) 
     &     call errquit(pname//'error creating g_tdy',0,GA_ERR)
           call ga_zero(g_tdy(i))
       if (.not.ga_create(mt_dbl,nbf,nbf,'TD-Z',nbf,0,g_tdz(i))) 
     &     call errquit(pname//'error creating g_tdz',0,GA_ERR)
           call ga_zero(g_tdz(i))
      end do
c
c     Compute dipole integrals in the AO basis
      call ga_zero(g_dpx)
      call ga_zero(g_dpy)
      call ga_zero(g_dpz)
      call int_dip_ga(ao_bas_han,ao_bas_han,g_dpx,g_dpy,g_dpz)
c
ctest
c      oskel = .false.
c      call int_1e_ga(ao_bas_han, ao_bas_han,g_dpx,'overlap',oskel)
ctest
c
c     Compute (Cf)T * D * (Ci)
      do i=1,ipol
c      x-comp
       call ga_zero(g_tmp)
       call ga_dgemm('T','N',nbf,nbf,nbf,1.d0,g_mos_fin(i),g_dpx,
     &        0.d0,g_tmp)
       call ga_dgemm('N','N',nbf,nbf,nbf,1.d0,g_tmp,g_mos_ini(i),0.d0,
     &                 g_tdx(i))
c      y-comp
       call ga_zero(g_tmp)
       call ga_dgemm('T','N',nbf,nbf,nbf,1.d0,g_mos_fin(i),g_dpy,
     &        0.d0,g_tmp)
       call ga_dgemm('N','N',nbf,nbf,nbf,1.d0,g_tmp,g_mos_ini(i),0.d0,
     &                 g_tdy(i))
c      z-comp
       call ga_zero(g_tmp)
       call ga_dgemm('T','N',nbf,nbf,nbf,1.d0,g_mos_fin(i),g_dpz,
     &        0.d0,g_tmp)
       call ga_dgemm('N','N',nbf,nbf,nbf,1.d0,g_tmp,g_mos_ini(i),0.d0,
     &                 g_tdz(i))
      end do ! ipol
c
c     Sum up the transition dipole up and down
      if (ipol.eq.2) then
       call ga_dadd(1.d0,g_tdx(1),1.d0,g_tdx(2),g_tdx(1))
       call ga_dadd(1.d0,g_tdy(1),1.d0,g_tdy(2),g_tdy(1))
       call ga_dadd(1.d0,g_tdz(1),1.d0,g_tdz(2),g_tdz(1))
      end if ! ipol
c
c     Calculate transition dipole
      tdx = ga_trace_diag(g_tdx(1))
      tdy = ga_trace_diag(g_tdy(1))
      tdz = ga_trace_diag(g_tdz(1))
c
c     Calculate the normalized intensity
      if (ipol.eq.1) then
        xps_wif = (tdx*tdx + tdy*tdy + tdz*tdz)/(nbf*nbf)
      else
        xps_wif = (tdx*tdx + tdy*tdy + tdz*tdz)/(4.d0*nbf*nbf)
      end if ! ipol
      xps_wif = 2.d0/3.d0 * xps_wif * xps_scale
c
      if (ga_nodeid().eq.0) 
     &      write(luout,*) "XPS Intensity: ",xps_wif
c
c     Clean up
      if (.not.ma_pop_stack(l_occ)) call errquit
     &  (pname//'failed to deallocate occ',0,MA_ERR)
      if (.not.ma_pop_stack(l_evl)) call errquit
     &  (pname//'failed to deallocate evl',0,MA_ERR)
c
      if (.not.ga_destroy(g_tmp)) 
     &  call errquit(pname//'failed to destroy g_tmp',0,GA_ERR)
      if (.not.ga_destroy(g_dpx)) 
     &  call errquit(pname//'failed to destroy g_dpx',0,GA_ERR)
      if (.not.ga_destroy(g_dpy)) 
     &  call errquit(pname//'failed to destroy g_dpy',0,GA_ERR)
      if (.not.ga_destroy(g_dpz)) 
     &  call errquit(pname//'failed to destroy g_dpz',0,GA_ERR)
c
      do i=1,ipol
        if (.not.ga_destroy(g_mos_ini(i))) 
     &    call errquit(pname//'failed to destroy g_mos_ini',0,GA_ERR)
        if (.not.ga_destroy(g_mos_fin(i))) 
     &    call errquit(pname//'failed to destroy g_mos_fin',0,GA_ERR)
        if (.not.ga_destroy(g_tdx(i))) 
     &    call errquit(pname//'failed to destroy g_tdx',0,GA_ERR)
        if (.not.ga_destroy(g_tdy(i))) 
     &    call errquit(pname//'failed to destroy g_tdy',0,GA_ERR)
        if (.not.ga_destroy(g_tdz(i))) 
     &    call errquit(pname//'failed to destroy g_tdz',0,GA_ERR)
      enddo
c
      end 
c $Id$
