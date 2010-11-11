c Subject : Implementation of Nuclear Magnetic Resonance
c           Chemical Shieldings (NMR-CS) using 
c           Scalar Relativistic ZORA (SRZ4) 
c Author  : Fredy Aquino
c Date    : 11-09-10
      subroutine get_NMRCS_SRZORA(rtdb,
     &                            g_dens_at,
     &                            nexc,
     &                            geom, 
     &                            ao_bas_han,
     &                            nbf,
     &                            noc,
     &                            ipol)
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh" 
#include "geom.fh" 
#include "zora.fh" 
c ------ defined in zora.fh ---- START
c  ga_dia     ! OUTPUT
c  ga_para1   ! OUTPUT
c  ga_h01_num ! OUTPUT
c  ga_Fji     ! OUTPUT
c ------ defined in zora.fh ---- END
      integer ga_dia1    
      integer ga_dia2    
      integer rtdb    
      integer g_dens_at(2)   
      integer geom
      integer ao_bas_han    
      integer nbf,noc(2),ipol,nexc
      integer nat_slc
      integer alo(3),ahi(3),ld(2)
      logical status  
      double precision ppm
      integer i,j
      external get_tensor_SR,get_Fji
      data ppm     /26.62566914d+00/   

c === allocate arrays to store diamagnetic tensor === START
c ++++++++ get-nat_slc++++++++++ START
      if (.not.rtdb_get(rtdb, 'giao:natoms', MT_INT, 1,nat_slc)) then  
         if (.not.geom_ncent(geom, nat_slc)) call                   
     &       errquit('get_NMRCS_SRZORA: geom_ncent failed ',0, GEOM_ERR)
      endif
      alo(1) = 3
      alo(2) = -1
      alo(3) = -1
      ahi(1) = 3
      ahi(2) = 3
      ahi(3) = nat_slc
      if (.not.nga_create(MT_DBL,3,ahi,'g_dia matrix',alo,ga_dia)) call 
     &    errquit('get_NMRCS_SRZORA: nga_create failed g_dia',0,GA_ERR)
      if (do_prntNMRCS) then  ! printNMRCS-if-requested-START
       if (ga_nodeid().eq.0) then
        write(*,*) "nat_slc=",nat_slc
       endif
      endif                               ! printNMRCS-if-requested-END
c ++++++++ get-nat_slc ++++++++++ END   
c ------- get_NMRCSdia1 --------- START
      if (do_prntNMRCS) then  ! printNMRCS-if-requested-START
       if (ga_nodeid().eq.0) 
     &     write(*,*) "==== get_NMRCSdia1 ======= START" 
      endif                               ! printNMRCS-if-requested-END
      call get_tensor_SR(ga_dia1,   ! OUTPUT
     &                   ga_dia2,   ! OUTPUT
     &                   ga_para1,  ! OUTPUT
     &                   ga_h01_num,! OUTPUT
     &                   rtdb,g_dens_at,nexc,
     &                   geom,ao_bas_han,nbf,
     &                   noc,ipol)
      call ga_add(1.0d0,ga_dia1,1.0d0,ga_dia2,ga_dia)

      if (do_prntNMRCS) then  ! printNMRCS-if-requested-START
      if (ga_nodeid().eq.0) then
       write(*,*) '----OUT:ga_dia-ppm----- START'
      endif
      call ga_print(ga_dia)
      if (ga_nodeid().eq.0) then
       write(*,*) '----OUT:ga_dia-ppm----- END'
      endif
      if (ga_nodeid().eq.0) then
       write(*,*) '----OUT:ga_dia2-ppm----- START'
      endif
      call ga_print(ga_dia2)
      if (ga_nodeid().eq.0) then
       write(*,*) '----OUT:ga_dia2-ppm----- END'
      endif
c +++++++++ using ga_dia2 as scratch-1 ++++++++ START
      if (ga_nodeid().eq.0) then
       write(*,*) '----OUT:ga_dia-1-NWChem-ppm----- START'
      endif
      call ga_add(1.0d0,ga_dia2,1.0d0,ga_para1,ga_dia2)
      call ga_print(ga_dia2)
      if (ga_nodeid().eq.0) then
       write(*,*) '----OUT:ga_dia-1-NWChem-ppm----- END'
      endif
c +++++++++ using ga_dia2 as scratch-1 ++++++++ END
c +++++++++ using ga_dia2 as scratch ++++++++ START
      if (ga_nodeid().eq.0) then
       write(*,*) '----OUT:ga_dia-NWChem-ppm----- START'
      endif
      call ga_add(1.0d0,ga_dia,1.0d0,ga_para1,ga_dia2)
      call ga_print(ga_dia2)
      if (ga_nodeid().eq.0) then
       write(*,*) '----OUT:ga_dia-NWChem-ppm----- END'
      endif
c +++++++++ using ga_dia2 as scratch ++++++++ END
      if (ga_nodeid().eq.0) then
       write(*,*) '----OUT:ga_para1-ppm----- START'
      endif
      call ga_print(ga_para1)
      if (ga_nodeid().eq.0) then
       write(*,*) '----OUT:ga_para1-ppm----- END'
      endif
      if (ga_nodeid().eq.0) then
       write(*,*) '----OUT:ga_h01----- START'
      endif
      call ga_print(ga_h01_num)
      if (ga_nodeid().eq.0) then
       write(*,*) '----OUT:ga_h01----- END'
      endif
      if (ga_nodeid().eq.0) 
     &     write(*,*) "==== get_NMRCSdia2par1 ======= END"
      endif                               ! printNMRCS-if-requested-END
c ------- get_NMRCSdia1 --------- END
      call get_Fji(ga_Fji, ! OUT: munu-mat-Fji
     &             rtdb,g_dens_at,nexc,
     &             geom,ao_bas_han,nbf)
      if (do_prntNMRCS) then  ! printNMRCS-if-requested-START
         if (ga_nodeid().eq.0) then
           write(*,*) '----OUT:ga_Fji------ START'
         endif
         call ga_print(ga_Fji)
         if (ga_nodeid().eq.0) then
           write(*,*) '----OUT:ga_Fji------ END'
         endif
      endif                               ! printNMRCS-if-requested-END  
      return
      end

      subroutine get_tensor_SR(gdia1, ! OUTPUT
     &                         gdia2, ! OUTPUT
     &                         gpar1, ! OUTPUT
     &                         g_h01, ! OUTPUT
     &                         rtdb,g_dens_at,nexc,
     &                         geom, 
     &                         ao_bas_han,
     &                         nbf,
     &                         noc,
     &                         ipol)
c    Author : Fredy Aquino
c    Date   : 11-06-10
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh" 
#include "geom.fh" 
#include "zora.fh" 
      integer gdia1,gdia2,gpar1, ! out
     &        g_h01              ! out
      integer rtdb    
      integer g_dens_at(2),g_densZ4(3)    
      integer g_v(3)
      integer noc(2),noc1,ipol
      integer geom,ao_bas_han                      
      integer ispin,nexc,iat,iat1,nat
      integer l_xyzpt,k_xyzpt,     
     &        l_zanpt,k_zanpt    
      integer l_AtNr,k_AtNr 
      logical status               
      character*16 at_tag  
      integer stat_read,read_SLCTD_NMR_Atoms
      integer alo(3),ahi(3),ld(2)
      integer i,j,k,t,a,nbf
      integer dims(3),chunk(3)
      double precision scf_dbl,scale2munu
      double precision xyz_NMRQcoords(3)
      integer g_dia1,g_dia2(3),g_dia3(3,3),
     &        g_nmr1(6),g_nmr2(18)
      external zora_getv_NMRCS_SR12,
     &         get_dia1,get_dia2par1h01,
     &         get_densZ4,read_SLCTD_NMR_Atoms,
     &         get_chi_centers_ga
      double precision ppm
      data ppm     /26.62566914d+00/ 
      scale2munu=0.5d0 ! /0.25d0/(clight_au*clight_au) ! 1/(4c^2)  
      if (ipol.eq.1) then
        scf_dbl=2.0d00
        noc1=noc(1)
      else if (ipol.eq.2) then 
        scf_dbl=1.0d00   
        noc1=noc(1)+noc(2)    
      endif
      status=geom_ncent(geom,nat) ! Get nat, # of atoms
c----- Allocate memory - FA
      if (.not. ma_alloc_get(mt_dbl,3*nat,'xyz pnt',l_xyzpt,k_xyzpt))
     &    call errquit('gCSSR: ma failed',911,MA_ERR)
      if (.not. ma_alloc_get(mt_dbl,nat,'zan pnt',l_zanpt,k_zanpt))
     &    call errquit('gCSSR: ma failed',911,MA_ERR)
c +++++++++++++++creating ga_arrays ++++++++START
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &      'gCSSR: g_dia1',0,0,g_dia1))
     $    call errquit('gCSSR: g_dia1',0,GA_ERR)
         call ga_zero(g_dia1)
      do i=1,3
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &      'gCSSR: g_dia2',0,0,g_dia2(i)))
     $    call errquit('gCSSR: g_dia2',0,GA_ERR)
         call ga_zero(g_dia2(i))
         do j=1,3
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &      'gCSSR: g_dia3',0,0,g_dia3(i,j)))
     $    call errquit('gCSSR: g_dia3',0,GA_ERR)
         call ga_zero(g_dia3(i,j))    
         enddo ! end-loop-j
      enddo ! end-loop-i
      do i=1,6
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &     'gCSSR: g_nmr1',0,0,g_nmr1(i)))
     $    call errquit('gCSSR: g_nmr1',0,GA_ERR)
         call ga_zero(g_nmr1(i))
      enddo ! end-loop-i
      do i=1,18
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &     'gCSSR: g_nmr2',0,0,g_nmr2(i)))
     $    call errquit('gCSSR: g_nmr2',0,GA_ERR)
         call ga_zero(g_nmr2(i))
      enddo ! end-loop-i
c +++++++++++++++creating ga_arrays ++++++++END
       dims(1) =nbf
       chunk(1)=nbf 
       do i=1,3
        if (.not. nga_create(mt_dbl,1,dims,"Array v",chunk,g_v(i)))
     $    call errquit('gCSSR: g_v', 0,GA_ERR)
       enddo
       call get_chi_centers_ga(g_v,ao_bas_han,nbf,geom,nat)
c +++++ Read Atom Nr for NMR calc ++
        if (.not. ga_create(mt_dbl,1,nat,
     &   'gCSSR: g_AtNr',0,0,g_AtNr)) 
     $   call errquit('gCSSR: g_AtNr',0,GA_ERR)
        call ga_zero(g_AtNr)
       stat_read=read_SLCTD_NMR_Atoms(rtdb,nat,nlist,g_AtNr)
c === allocate arrays to store diamagnetic tensor === START
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 3*nlist
      if (.not.nga_create(MT_DBL,3,ahi,'H01 matrix',alo,g_h01)) call 
     &    errquit('get_d2p1: nga_create failed g_h01_num',0,GA_ERR)
      call ga_zero(g_h01)
      alo(1) =  3
      alo(2) = -1
      alo(3) = -1
      ahi(1) =  3
      ahi(2) =  3
      ahi(3) =  nlist ! Total nr. atoms requested
      if (.not.nga_create(MT_DBL,3,ahi,'gdia1 matrix',
     &                    alo,gdia1)) 
     &  call errquit('gCSSR: nga_create failed gdia1 all',
     &               0,GA_ERR)
      if (.not.nga_create(MT_DBL,3,ahi,'gdia2 matrix',
     &                    alo,gdia2)) 
     &    call errquit('NMRCSdia2par1: nga_create failed gdia2',
     &            0,GA_ERR)
      if (.not.nga_create(MT_DBL,3,ahi,'gpar1 matrix',
     &                    alo,gpar1)) 
     &    call errquit('NMRCSdia2par1: nga_create failed gpar1',
     &            0,GA_ERR)
c === allocate arrays to store diamagnetic tensor === END
c  Allocate memory for l_AtNr,k_AtNr 
      if (.not.ma_alloc_get(mt_dbl,nat,
     &  'AtNr',l_AtNr,k_AtNr))
     &  call errquit('gCSSR: ma failed',0,MA_ERR)
      call ga_get(g_AtNr,1,1,1,nat,dbl_mb(k_AtNr),1)
      call get_densZ4(rtdb,ao_bas_han,geom,g_densZ4)
      do iat1=1,nlist  ! nlist <= nat
       iat=dbl_mb(k_AtNr+iat1-1)
       status=geom_cent_get(geom,iat,at_tag,
     &                      dbl_mb(k_xyzpt+3*(iat-1)),
     &                      dbl_mb(k_zanpt+iat-1))
       xyz_NMRQcoords(1)= dbl_mb(k_xyzpt  +3*(iat-1))
       xyz_NMRQcoords(2)= dbl_mb(k_xyzpt+1+3*(iat-1))
       xyz_NMRQcoords(3)= dbl_mb(k_xyzpt+2+3*(iat-1))
       call zora_getv_NMRCS_SR12(rtdb,g_dens_at, 
     &                           xyz_NMRQcoords,
     &                           g_dia1,g_dia2,g_dia3,! out: munu matrix
     &                           g_nmr1,g_nmr2,       ! out: munu matrix
     &                           nexc) 
       call get_dia1(gdia1, ! out
     &               g_dia1,g_dia2,g_dia3,g_densZ4, 
     &               iat1,nbf,ipol,xyz_NMRQcoords,g_v)
       call get_dia2par1h01(gdia2, ! out
     &                      gpar1, ! out
     &                      g_h01, ! out
     &                      iat1,ipol,nbf, 
     &                      g_nmr1,g_nmr2,
     &                      g_v,g_densZ4)
      end do ! iat loop
c ---Destroying ga arrays ------- START
      do i=1,6
        if (.not. ga_destroy(g_nmr1(i))) call errquit(
     &    'gCSSR: ga_destroy failed ',0, GA_ERR) 
      enddo ! end-loop-i
      do i=1,18
        if (.not. ga_destroy(g_nmr2(i))) call errquit(
     &    'gCSSR: ga_destroy failed ',0, GA_ERR) 
      enddo ! end-loop-i
        if (.not. ga_destroy(g_dia1)) call errquit(
     &    'gCSSR: ga_destroy failed ',0, GA_ERR) 
       do i=1,3
        if (.not. ga_destroy(g_dia2(i))) call errquit(
     &    'gCSSR: ga_destroy failed ',0, GA_ERR) 
        do j=1,3
        if (.not. ga_destroy(g_dia3(i,j))) call errquit(
     &    'gCSSR: ga_destroy failed ',0, GA_ERR) 
        enddo
       enddo
       do i=1,3
        if (.not. ga_destroy(g_v(i))) call errquit(
     &    'NMRCSdia2par1: ga_destroy failed ',0, GA_ERR)
        if (.not. ga_destroy(g_densZ4(i))) call errquit(
     &    'NMRCSdia2par1: ga_destroy failed ',0, GA_ERR)
       enddo
c ---Destroying ga arrays ------- END
c----deallocate memory 
      if (.not.ma_free_heap(l_zanpt)) call errquit
     &   ('gCSSR:, ma_free_heap of l_zanpt failed',911,MA_ERR)
      if (.not.ma_free_heap(l_xyzpt)) call errquit
     &   ('gCSSR:, ma_free_heap of l_xyzpt failed',911,MA_ERR)
      if (.not.ma_free_heap(l_AtNr)) call
     &    errquit('gCSSR:: ma_free_heap l_AtNr',0, MA_ERR)
      return
      end

      subroutine get_dia1(gdia1,         ! out
     &                    g_dia1,        ! in:
     &                    g_dia2,        ! in:
     &                    g_dia3,        ! in: 
     &                    g_densZ4,      ! in:
     &                    iat,           ! in: atom nr.
     &                    nbf,           ! in: nr. basis functions
     &                    ipol,          ! in: nr. polarizations         
     &                    xyz_NMRQcoords,! in: nuclear coord.
     &                    g_R)           ! in: R_{mu}
c    Purpose: Evaluating dia1
c             1st component of diamagnetic tensor.
c             =P_{munu} g_dia3
c    Info about g_dia3 (after modification) 
c    k ne t: -I_{munu,kt}
c    k eq t: \sum_{p ne k} I_{munu,pp}
c            = I_{munu,22}+I_{munu,33}  k=t=1
c            = I_{munu,33}+I_{munu,11}  k=t=2
c            = I_{munu,11}+I_{munu,22}  k=t=3
c
c    I_{munu,kt}=2 <chi_mu|K r_k r_t/r_Q^3|chi_nu> +
c                R_{Q,k}(R_{mu,t}+R_{nu,t}) <chi_mu|K/r_Q^3|chi_nu> -
c                (R_{mu,t}+R_{nu,t})    
c  Output: g_dia3 (=I_{munu,kt})
c  Input: 
c  g_dia1, <chi_mu|K        /r_Q^3|chi_nu>
c  g_dia2, <chi_mu|K     r_k/r_Q^3|chi_nu>    t=1,2,3
c  g_dia3, <chi_mu|K r_k r_t/r_Q^3|chi_nu>  k,t=1,2,3 
c  FA-11-05-10
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "zora.fh"
      integer gdia1 ! 1st term of dia tensor
      integer g_dia1,g_dia2(3),
     &        g_dia3(3,3) ! output/input
      integer g_densZ4(3)
      integer g_R(3)
      integer k,t,cbuf,ispin,ipol,iat,nbf
      integer g_t1,g_t2   
      integer alo(3),ahi(3),ld(2)
      integer l_buf,k_buf
      double precision ac,n1,scale2munu
      double precision xyz_NMRQcoords(3),scl
      double precision ppm
      data ppm     /26.62566914d+00/ 
c ---- Creating scratch ga-arrays --------- START
       if (.not. ga_create(mt_dbl,nbf,nbf,
     &     'gd1: g_t1',0,0,g_t1))
     $ call errquit('gd1: g_t1',0,GA_ERR)
       call ga_zero(g_t1)
       if (.not. ga_create(mt_dbl,nbf,nbf,
     &     'gd1: g_t2',0,0,g_t2))
     $ call errquit('gd1: g_t2',0,GA_ERR)
       call ga_zero(g_t2)
c ---- Creating scratch ga-arrays --------- END
      if(.not.ma_alloc_get(mt_dbl,3*3,'gd1:buf',
     &                    l_buf,k_buf))
     $  call errquit('gd1: ma failed',3*3, MA_ERR)
      scale2munu=0.5d0
      do k=1,3
       do t=1,3
       scl=xyz_NMRQcoords(k)
       call ga_scale(g_dia3(k,t),2.0d0)        !  +2.0 <chi_mu|r_k r_t/R_Q^3|chi_nu>
       call ga_add(scl,g_dia1,-1.0d0,g_dia2(k),g_t1) 
       call ga_copy(g_t1,g_t2)
       call ga_scale_rows(g_t1,g_R(t))         ! R_{nu,t} g_z
       call ga_scale_cols(g_t2,g_R(t))         ! R_{nu,t} g_z
       call ga_add(1.0d0,g_t1,1.0d0,g_t2,g_t1)       
       call ga_add(1.0d0,g_dia3(k,t),1.0d0,g_t1,g_dia3(k,t))
       call ga_add(1.0d0,g_dia3(k,t),-2.0d0*scl,g_dia2(t),
     &             g_dia3(k,t))
       if (k.ne.t)
     &   call ga_scale(g_dia3(k,t),-1.0d0)
       enddo ! end-loop-t     
      enddo ! end-loop-k
      call ga_add(1.0d0,g_dia3(2,2),1.0d0,g_dia3(3,3),g_t1)
      call ga_add(1.0d0,g_dia3(3,3),1.0d0,g_dia3(1,1),g_t2)
      call ga_add(1.0d0,g_dia3(1,1),1.0d0,g_dia3(2,2),g_dia3(3,3))
      call ga_copy(g_t1,g_dia3(1,1))
      call ga_copy(g_t2,g_dia3(2,2))
       do k=1,3  ! k=x,y,z
        do t=1,3 ! t=x,y,z    
         ac=0.0d0
         do ispin=1,ipol
          n1=ga_ddot(g_densZ4(ispin),g_dia3(k,t)) 
          ac=ac+n1
          if (do_prntNMRCS) then  ! printNMRCS-if-requested-START
           if (ga_nodeid().eq.0) then
           write(*,15) k,t,iat,ispin,
     &                 n1*scale2munu*ppm, 
     &                 ac*scale2munu*ppm       
 15        format('dia1(k,t,iat,ispin)=(',
     &            i3,',',i3,',',i3,',',i3,')=(',
     &            f15.8,',',f15.8,') ppm')  
           endif
          endif
         enddo ! ispin-loop
         cbuf=k_buf+(k-1)*3+t-1
         dbl_mb(cbuf)=ac*scale2munu*ppm
         enddo ! t loop
       enddo ! k loop
       alo(1)=1
       ahi(1)=3
       alo(2)=1
       ahi(2)=3
       alo(3)=iat
       ahi(3)=iat
       ld(1)=3
       ld(2)=3
       call nga_put(gdia1,alo,ahi,dbl_mb(k_buf),ld)
c ---Destroying ga arrays 
       if (.not. ga_destroy(g_t1)) call errquit(
     &    'get_d1: ga_destroy failed ',0, GA_ERR)
       if (.not. ga_destroy(g_t2)) call errquit(
     &    'get_d1: ga_destroy failed ',0, GA_ERR)
      if (.not.ma_free_heap(l_buf)) call errquit
     &   ('get_d1: ma_free_heap of l_buf failed',
     &   911,MA_ERR)
      return
      end

      subroutine get_dia2par1h01(gdia2,   ! out
     &                           gpar1,   ! out
     &                           g_h01,   ! out
     &                           iat,     ! in: atom index
     &                           ipol,    ! in: nr. polarizations
     &                           nbf,     ! in: nr. basis functions
     &                           g_nmr1,  ! in:
     &                           g_nmr2,  ! in:
     &                           g_R,     ! in: R_{mu}
     &                           g_densZ4)! in: e-density
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "zora.fh"
      integer gdia2,gpar1
      integer g_h01
      integer g_nmr1(6),g_nmr2(18),g_R(3)
      integer g_densZ4(3)
      integer g_h1,g_t1,g_t2,g_t3,g_M,g_Nmn(3,3)
      integer iat,nbf
      integer j,k,t,a,b,cbuf,ipol
      double precision ac,ad,n1,n2
      double precision scl2munu
      integer alo(3),ahi(3),ld(2)
      integer l_buf,k_buf,
     &        l_buf1,k_buf1,
     &        l_buf2,k_buf2
      double precision ppm
      data ppm     /26.62566914d+00/   
      integer ind_tmn(2,3)
      data ind_tmn / 2, 3,  ! tmn=123
     &               3, 1,  ! tmn=231
     &               1, 2 / ! tmn=312
      integer indxs(3,3),ind
      data indxs / 1, 2, 3,  ! a=1
     &             7, 8, 9,  ! a=2
     &            13,14,15 / ! a=3
      external get_scld_A
      scl2munu =0.5d0  ! 1/(2c)
c ---- Creating scratch ga-arrays --------- START
       if (.not. ga_create(mt_dbl,nbf,nbf,
     &     'gd2p1: g_t1',0,0,g_t1))
     $ call errquit('get_d2p1: g_t1',0,GA_ERR)
       call ga_zero(g_t1)
       if (.not. ga_create(mt_dbl,nbf,nbf,
     &     'gd2p1: g_t1',0,0,g_t2))
     $ call errquit('get_d2p1: g_t2',0,GA_ERR)
       call ga_zero(g_t2)
       if (.not. ga_create(mt_dbl,nbf,nbf,
     &     'gd2p1: g_t1',0,0,g_t3))
     $ call errquit('get_d2p1: g_t3',0,GA_ERR)
       call ga_zero(g_t3)
       if (.not. ga_create(mt_dbl,nbf,nbf,
     &     'gd2p1: g_t1',0,0,g_h1))
     $ call errquit('get_d2p1: g_h1',0,GA_ERR)
       call ga_zero(g_h1)
       if (.not. ga_create(mt_dbl,nbf,nbf,
     &     'gd2p1: g_t1',0,0,g_M))
     $ call errquit('get_d2p1: g_M',0,GA_ERR)
       call ga_zero(g_M)
       do k=1,3
        do t=1,3
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &     'gd2p1: g_Nmn',0,0,g_Nmn(k,t)))
     $   call errquit('get_d2p1: g_Nmn',0,GA_ERR)
         call ga_zero(g_Nmn(k,t))
        enddo
       enddo
c ---- Creating scratch ga-arrays --------- END
      if(.not.ma_alloc_get(MT_DBL,nbf*nbf,'get_d2p1',
     &                    l_buf,k_buf))
     $ call errquit('get_d2p1: ma failed',911, MA_ERR)
      if(.not.ma_alloc_get(mt_dbl,3*3,'gd2p1:buf1',
     &                    l_buf1,k_buf1))
     $ call errquit('get_d2p1: ma failed',3*3, MA_ERR)
      if(.not.ma_alloc_get(mt_dbl,3*3,'gd2p1:buf2',
     &                    l_buf2,k_buf2))
     $ call errquit('get_d2p1: ma failed',3*3, MA_ERR)
       do t=1,3
        do a=1,3
         ind=indxs(t,a) 
         call ga_add(1.0d0,g_nmr2(ind)  ,
     &              -1.0d0,g_nmr2(ind+3),g_t1)
         call ga_transpose(g_t1,g_t2)
         call ga_add(1.0d0,g_t1,-1.0d0,g_t2,g_Nmn(t,a))                    
        enddo ! end-loop-a
       enddo ! end-loop-t
       do t=1,3
c ----- g_h01 = < chi_{mu} | h_t^{01}| chi_{nu} >
c       h_t^{01}=K/(2c) (\vec{r} x \vec{p})_t/r_Q^3 +
c                       (\vec{r} x \vec{p})_t/r_Q^3 K/(2c)
        call ga_add(1.0d0,g_nmr1(t)  ,
     &             -1.0d0,g_nmr1(t+3),g_t1)
        call ga_transpose(g_t1,g_t2)
        call ga_add(1.0d0,g_t2,-1.0d0,g_t1,g_h1) ! g_h1=g_h01
        do k=1,3
         a=ind_tmn(1,k)
         b=ind_tmn(2,k)
         call ga_zero(g_M)
         call ga_copy(g_Nmn(t,a),g_t1)
         call get_scld_A(g_t1,g_R(b),g_t2) ! output : g_t1
         call ga_add(1.0d0,g_M, 1.0d0,g_t1,g_M)
         call ga_copy(g_Nmn(t,b),g_t1)
         call get_scld_A(g_t1,g_R(a),g_t2) ! output : g_t1
         call ga_add(1.0d0,g_M,-1.0d0,g_t1,g_M)
c g_M= R_{numu,b'}[N_{munu,mn,a'}-N_{numu,mn,a'}^*-
c                 (N_{munu,nm,a'}-N_{numu,nm,a'}^*)]-
c      R_{numu,a'}[N_{munu,mn,b'}-N_{numu,mn,b'}^*-
c                 (N_{munu,nm,b'}-N_{numu,nm,b'}^*)]
         call ga_copy(g_h1,g_t1)
         call ga_scale_cols(g_t1,g_R(b))          ! R_{nu,b} g_N
         call ga_scale_rows(g_t1,g_R(a))          ! R_{mu,a} [R_{nu,b} g_N] -> g_t1
         call ga_copy(g_h1,g_t2)
         call ga_scale_cols(g_t2,g_R(a))          ! R_{nu,b} g_N
         call ga_scale_rows(g_t2,g_R(b))          ! R_{mu,a} [R_{nu,b} g_N] -> g_t2
         call ga_add(1.0d0,g_t1,-1.0d0,g_t2,g_t3) ! g_t3=-4c^2 <chi_mu|i/(2c)(R_mu x R_nu)_k h_t^{01}|chi_nu>
         call ga_add(1.0d0,g_M,1.0d0,g_t3,g_M)
         ac=0.0d0
         ad=0.0d0
         do j=1,ipol    
          n1=ga_ddot(g_densZ4(j),g_t3) ! for par1
          ac=ac+n1
          n2=ga_ddot(g_densZ4(j),g_M)  ! for dia2
          ad=ad+n2
          if (do_prntNMRCS) then  ! printNMRCS-if-requested-START
           if (ga_nodeid().eq.0) then
           write(*,25) k,t,iat,j,
     &                 n2*scl2munu*ppm, 
     &                 ad*scl2munu*ppm,
     &                 n1*scl2munu*ppm, 
     &                 ac*scl2munu*ppm    
 25        format('(dia2,par1)(k,t,iat,spin)=(',
     &            i3,',',i3,',',i3,',',i3,')=(',
     &            f15.8,',',f15.8,',',
     &            f15.8,',',f15.8,') ppm')  
           endif
          endif
         enddo ! end-loop-j  
         cbuf=k_buf1+(k-1)*3+t-1
         dbl_mb(cbuf)=-ac*scl2munu*ppm ! par1
         cbuf=k_buf2+(k-1)*3+t-1
         dbl_mb(cbuf)= ad*scl2munu*ppm ! dia2
        enddo ! end-loop-k
        alo(1)=1
        ahi(1)=nbf
        alo(2)=1
        ahi(2)=nbf
        alo(3)=3*(iat-1)+t
        ahi(3)=3*(iat-1)+t
        ld(1)=nbf
        ld(2)=nbf
        call ga_scale(g_h1,-0.5d0)
        call ga_get(g_h1,1,nbf,1,nbf,dbl_mb(k_buf),nbf)
        call nga_put(g_h01,alo,ahi,dbl_mb(k_buf),ld) ! store g_h01
       enddo ! end-loop-t
       alo(1)=1
       ahi(1)=3
       alo(2)=1
       ahi(2)=3
       alo(3)=iat
       ahi(3)=iat
       ld(1)=3
       ld(2)=3
       call nga_put(gpar1,alo,ahi,dbl_mb(k_buf1),ld)
       call nga_put(gdia2,alo,ahi,dbl_mb(k_buf2),ld)      
c----deallocate memory
       if (.not.ma_free_heap(l_buf)) call errquit
     &    ('get_d2p1, ma_free_heap of l_buf failed',
     &      911,MA_ERR)
       if (.not.ma_free_heap(l_buf1)) call errquit
     &    ('get_d2p1, ma_free_heap of l_buf failed',
     &      911,MA_ERR)
       if (.not.ma_free_heap(l_buf2)) call errquit
     &    ('get_d2p1, ma_free_heap of l_buf failed',
     &      911,MA_ERR)
c ---Destroying ga arrays 
       if (.not. ga_destroy(g_t1)) call errquit(
     &    'get_d2p1: ga_destroy failed ',0, GA_ERR)
       if (.not. ga_destroy(g_t2)) call errquit(
     &    'get_d2p1: ga_destroy failed ',0, GA_ERR)
       if (.not. ga_destroy(g_t3)) call errquit(
     &    'get_d2p1: ga_destroy failed ',0, GA_ERR)
       if (.not. ga_destroy(g_h1)) call errquit(
     &    'get_d2p1: ga_destroy failed ',0, GA_ERR)
       if (.not. ga_destroy(g_M)) call errquit(
     &    'get_d2p1: ga_destroy failed ',0, GA_ERR)
       do k=1,3
        do t=1,3
         if (.not. ga_destroy(g_Nmn(k,t))) 
     &    call errquit(
     &    'get_d2p1: ga_destroy failed ',0, GA_ERR)
        enddo
       enddo
      return
      end

      subroutine get_Fji(g_Fji, ! OUT: munu-mat-Fji
     &                   rtdb,g_dens_at,nexc,
     &                   geom,ao_bas_han,nbf)
c    Purpose: Calculates 
c             F_{ji}^{1k}(1st and 3rd terms in write-up)=
c             < chi_{mu}| K/(4c) [\vect{r}_{mu} x \vect{p}]_k +
c                                [\vect{r}_{nu} x \vect{p}]_k K/(4c) | chi_{nu} > +
c             < chi_{mu}| \vect{p} (iK/(4c)) 
c                        [\vect{r} x \vect{R}_{numu}]_k \vect{p} | chi_{nu} > 
c    Author: Fredy Aquino
c    Date  : 11-08-10
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh" 
#include "geom.fh" 
      integer g_Fji ! OUTPUT
      logical status        
      integer rtdb,g_dens_at(2),geom,ao_bas_han       
      integer nexc            
      integer i,j,k,m,n    
      integer nat,nbf
      double precision scl2munu
      integer g_v(3)
      integer g_nmr(3),g_nmr3(3),g_nmr4(6)
      integer g_t1,g_t2,g_t3,g_Lk
      integer l_buf,k_buf
      integer alo(3),ahi(3),ld(2)
      integer chunk(3),dims(3)
      integer ind_kmn(2,3)
      data ind_kmn / 2, 3,  ! kmn=123
     &               3, 1,  ! kmn=231
     &               1, 2 / ! kmn=312
      external get_chi_centers_ga,get_scld_A,
     &         zora_getv_NMRCS_SR34
c +++++++++++++++creating ga_arrays ++++++++START
      do i=1,3
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &      'gFji: g_nmr',0,0,g_nmr(i)))
     $    call errquit('gFji: g_nmr',0,GA_ERR)
         call ga_zero(g_nmr(i))
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &      'gFji: g_nmr3',0,0,g_nmr3(i)))
     $    call errquit('gFji: g_nmr3',0,GA_ERR)
         call ga_zero(g_nmr3(i))
      enddo ! end-loop-i
      do i=1,6
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &     'gFji: g_nmr4',0,0,g_nmr4(i)))
     $    call errquit('gFji: g_nmr4',0,GA_ERR)
         call ga_zero(g_nmr4(i))
      enddo ! end-loop-i
c +++++++++++++++creating ga_arrays ++++++++END
       call zora_getv_NMRCS_SR34(rtdb,g_dens_at, 
     &                           g_nmr,         ! out: munu matrix
     &                           g_nmr3,g_nmr4, ! out: munu matrix
     &                           nexc) 
      if(.not.ma_alloc_get(mt_dbl,nbf*nbf,'gFji:buf',
     &                    l_buf,k_buf))
     $     call errquit('gFji: ma failed',911, MA_ERR)
      dims(1) =nbf
      chunk(1)=nbf 
      do i=1,3
        if (.not. nga_create(mt_dbl,1,dims,"Array v",chunk,g_v(i)))
     $       call errquit('gFji: g_v', 0,GA_ERR)
      enddo
      status=geom_ncent(geom,nat)
      call get_chi_centers_ga(g_v,ao_bas_han,nbf,geom,nat)
      scl2munu=0.5d0 ! 0.25d0 !/clight_au
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 3
      if (.not.nga_create(MT_DBL,3,ahi,'Fji matrix',alo,g_Fji)) 
     &    call 
     &    errquit('gFji: nga_create failed g_Fji',
     &            0,GA_ERR)
      call ga_zero(g_Fji)
c ---- Creating scratch ga-arrays --------- START
       if (.not. ga_create(mt_dbl,nbf,nbf,
     &     'gFji: g_t1',0,0,g_t1))
     $ call errquit('gFij: g_t1',0,GA_ERR)
       call ga_zero(g_t1)
       if (.not. ga_create(mt_dbl,nbf,nbf,
     &     'gFji: g_t2',0,0,g_t2))
     $ call errquit('gFij: g_t2',0,GA_ERR)
       call ga_zero(g_t2)
       if (.not. ga_create(mt_dbl,nbf,nbf,
     &     'gFji: g_t3',0,0,g_t3))
     $ call errquit('gFij: g_t3',0,GA_ERR)
       call ga_zero(g_t3)
       if (.not. ga_create(mt_dbl,nbf,nbf,
     &     'gFji: g_Lk',0,0,g_Lk))
     $ call errquit('gFij: g_Lk',0,GA_ERR)
c ---- Creating scratch ga-arrays --------- END
      do k=1,3   
       call ga_zero(g_Lk)    
       m=ind_kmn(1,k)
       n=ind_kmn(2,k)
       call ga_copy(g_nmr3(n),g_t2)
       call ga_scale_rows(g_t2,g_v(m))  ! R_{mu,m} g_z
       call ga_add(1.0d0,g_nmr4(k  ),-1.0d0,g_t2,g_nmr4(k))   ! g_nmr4(k)  =L_{numu,mn}^*
       call ga_copy(g_nmr3(m),g_t1)
       call ga_scale_rows(g_t1,g_v(n))  ! R_{mu,n} g_z
       call ga_add(1.0d0,g_nmr4(k+3),-1.0d0,g_t1,g_nmr4(k+3)) ! g_nmr4(k+3)=L_{numu,nm}^*
       call ga_add(1.0d0,g_nmr4(k),-1.0d0,g_nmr4(k+3),g_Lk)
       call ga_transpose(g_Lk,g_t1)
       call ga_add(1.0d0,g_Lk,-1.0d0,g_t1,g_Lk)
       call ga_copy(g_nmr(m),g_t1)
       call get_scld_A(g_t1,g_v(n),g_t2) ! output : g_t1
       call ga_zero(g_t3)
       call ga_add(1.0d0,g_t3, 1.0d0,g_t1,g_t3)
       call ga_copy(g_nmr(n),g_t2)
       call get_scld_A(g_t2,g_v(m),g_t1) ! output : g_t2
       call ga_add(1.0d0,g_t3,-1.0d0,g_t2,g_t3)
       call ga_add(-scl2munu,g_Lk,
     &              scl2munu,g_t3,g_Lk)  ! update g_Lmn
       alo(1)=1
       ahi(1)=nbf
       alo(2)=1
       ahi(2)=nbf
       alo(3)=k
       ahi(3)=k
       ld(1) =nbf
       ld(2) =nbf
       call ga_get(g_Lk,1,nbf,1,nbf,dbl_mb(k_buf),nbf)
       call nga_put(g_Fji,alo,ahi,dbl_mb(k_buf),ld)
      enddo ! end-loop-k
c ---Destroying ga arrays ----- START
       if (.not. ga_destroy(g_t1)) call errquit(
     &    'gFij: ga_destroy failed ',0, GA_ERR)
       if (.not. ga_destroy(g_t2)) call errquit(
     &    'gFij: ga_destroy failed ',0, GA_ERR)
       if (.not. ga_destroy(g_t3)) call errquit(
     &    'gFij: ga_destroy failed ',0, GA_ERR)
       if (.not. ga_destroy(g_Lk)) call errquit(
     &    'gFij: ga_destroy failed ',0, GA_ERR)
       do i=1,3
        if (.not. ga_destroy(g_nmr(i))) call errquit(
     &    'gFij: ga_destroy failed ',0, GA_ERR) 
        if (.not. ga_destroy(g_nmr3(i))) call errquit(
     &    'gFij: ga_destroy failed ',0, GA_ERR) 
       enddo
      do i=1,6
        if (.not. ga_destroy(g_nmr4(i))) call errquit(
     &    'gFij: ga_destroy failed ',0, GA_ERR) 
      enddo ! end-loop-i
c ---Destroying ga arrays ----- END
c----deallocate memory
       if (.not.ma_free_heap(l_buf)) call errquit
     &    ('gFij, ma_free_heap of l_buf failed',911,MA_ERR)
      return
      end

      integer function read_SLCTD_NMR_Atoms
     &                 (rtdb,nat,nlist,g_AtNr)
c---- GA output: g_AtNr   
c

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
#include "rtdb.fh" 
#include "context.fh"

      integer rtdb,ii,nlist,nat,g_AtNr
      integer atomnr(nat)
      integer giaoatoms
      double precision AtNr_dbl

      if (.not. rtdb_get(rtdb, 'giao:natoms',mt_int,
     $                     1,giaoatoms))
     &  giaoatoms=0 ! reset       
         if (giaoatoms.eq.0) then
          giaoatoms=nat
           nlist=giaoatoms
           do ii=1,giaoatoms
            AtNr_dbl=ii
            call ga_put(g_AtNr,1,1,ii,ii,AtNr_dbl,1)
           enddo     
         else
          if (.not. rtdb_get(rtdb, 'giao:atom list',mt_int,
     $                     giaoatoms,atomnr))
     $      call errquit('prop_input-giao: rtdb_get failed',
     $                   555, RTDB_ERR)
           nlist=giaoatoms
           do ii=1,giaoatoms
            AtNr_dbl=atomnr(ii)
            call ga_put(g_AtNr,1,1,ii,ii,AtNr_dbl,1)
           enddo
         endif 
       read_SLCTD_NMR_Atoms = 1
       return
       end   

      subroutine print_NMRCS_SRZORA_version()
       implicit none
       write(*,*) "+++++++++++++++++++++++++++++++++++"
       write(*,*) "+ NMR Shieldings (11-09-10 12:00) +"
       write(*,*) "+ Scalar Relativistic ZORA        +"
       write(*,*) "+ Author: Fredy W. Aquino         +"
       write(*,*) "+++++++++++++++++++++++++++++++++++"
       return
       end
