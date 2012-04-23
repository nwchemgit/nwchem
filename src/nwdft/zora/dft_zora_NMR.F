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
     &                            ipol,
     &                            ga_Fji,
     &                            g_densZ4)
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh" 
#include "geom.fh" 
#include "zora.fh"
c  #include "zora.fh" 
      integer g_densZ4(3)
      integer ga_dia1    
      integer ga_dia2    
      integer rtdb    
      integer g_dens_at(2)   
      integer geom
      integer ao_bas_han    
      integer nbf,noc(2),ipol,nexc
      integer nat_slc,type_nmrdata
      integer alo(3),ahi(3),ld(2)
      logical status
      logical do_prntNMRCS
      character*255 zorafilename
      double precision ppm
      integer i,j
      logical dft_zoraNMR_write
      integer ga_dia,    ! OUTPUT
     &        ga_para1,  ! OUTPUT
     &        ga_h01_num,! OUTPUT
     &        ga_Fji     ! OUTPUT
 
      external get_tensor_SR,get_Fji,
     &         dft_zoraNMR_write
      data ppm     /26.62566914d+00/   

c === allocate arrays to store diamagnetic tensor === START
c ++++++++ get-nat_slc++++++++++ START
      if (.not.rtdb_get(rtdb, 'giao:natoms', MT_INT, 1,nat_slc)) then  
         if (.not.geom_ncent(geom, nat_slc)) call                   
     &       errquit('get_NMRCS_SRZORA: geom_ncent failed ',0, GEOM_ERR)
      endif
      alo(1) =  3
      alo(2) = -1
      alo(3) = -1
      ahi(1) =  3
      ahi(2) =  3
      ahi(3) = nat_slc
      if (.not.nga_create(MT_DBL,3,ahi,'g_dia matrix',alo,ga_dia)) call 
     &    errquit('get_NMRCS_SRZORA: nga_create failed g_dia',0,GA_ERR)
      if(.not.rtdb_get(rtdb,'zora:do_prntNMRCS',        ! FA
     &                 mt_log,1,do_prntNMRCS))          ! FA
     &  do_prntNMRCS= .false.      
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
     &                   noc,ipol,g_densZ4)
      call ga_add(1.0d0,ga_dia1,1.0d0,ga_dia2,ga_dia)
c     if (ga_nodeid().eq.0) write(*,*) '----OUT:ga_dia-ppm--- START'
c     call ga_print(ga_dia)
c     if (ga_nodeid().eq.0) write(*,*) '----OUT:ga_dia-ppm--- END'
c     if (ga_nodeid().eq.0) write(*,*) '----OUT:ga_para1-ppm- START'
c     call ga_print(ga_para1)
c     if (ga_nodeid().eq.0) write(*,*) '----OUT:ga_para1-ppm- END'
      if (do_prntNMRCS) then  ! printNMRCS-if-requested-START
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
     &               rtdb,g_dens_at,nexc,
     &               geom,ao_bas_han,nbf)
      if (do_prntNMRCS) then  ! printNMRCS-if-requested-START
         if (ga_nodeid().eq.0) then
           write(*,*) '----OUT:ga_Fji------ START'
         endif
         call ga_print(ga_Fji)
         if (ga_nodeid().eq.0) then
           write(*,*) '----OUT:ga_Fji------ END'
         endif
      endif                               ! printNMRCS-if-requested-END  
c ----- Store nmrcs data in a file ------- START
c       Note.- lbl_nmrcs defined in zora.fh
        call util_file_name(lbl_nmrcs,.false.,.false.,zorafilename)
        type_nmrdata=1 ! =1,2,3=shieldings,hyperfine,gshift
        if (.not.dft_zoraNMR_write(
     &         zorafilename,
     &         type_nmrdata, ! =1,2,3=shieldings,hyperfine,gshift
     &         nbf,
     &         nat_slc,
     &         g_AtNr,
     &         ga_dia,
     &         ga_para1,
     &         ga_h01_num,
     &         ga_Fji))
     &     call errquit('get_NMRCS_SRZORA: dft_zoraNMR_write failed',
     &                  0,DISK_ERR)
c ----- Store nmrcs data in a file ------- END
c ---- Destroy stored ga arrays ------ START
           if (.not. ga_destroy(ga_dia)) call errquit(
     &    'get_NMRCS_SRZORA: ga_destroy failed ',0, GA_ERR)  
           if (.not. ga_destroy(ga_dia1)) call errquit(
     &    'get_NMRCS_SRZORA: ga_destroy failed ',0, GA_ERR)  
           if (.not. ga_destroy(ga_dia2)) call errquit(
     &    'get_NMRCS_SRZORA: ga_destroy failed ',0, GA_ERR)  
        if (.not. ga_destroy(ga_para1)) call errquit(
     &    'get_NMRCS_SRZORA: ga_destroy failed ',0, GA_ERR)  
        if (.not. ga_destroy(ga_h01_num)) call errquit(
     &    'get_NMRCS_SRZORA: ga_destroy failed ',0, GA_ERR)        
c ---- Destroy stored ga arrays ------ END
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
     &                         ipol,g_densZ4)
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
      integer g_v(3),vectors(2)
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
      double precision scf_dbl
      double precision xyz_NMRQcoords(3)
      integer g_dia1,g_dia2(3),g_dia3(3,3),
     &        g_nmr1(6),g_nmr2(18)
c --------- definitions for shieldings NLMO analysis --- START
      character*255 zorafilename
      integer g_dia2kt(3,3),g_par1kt(3,3)
      integer g_munudia,g_munupar1,g_munu_h01,g_c1,ndata
      integer count,count1,count2,n_munu,ndir,ndir1,shldfile
      logical dft_zoraShield_NLMOAnalysis_write
      external fill_munuFCSD,dft_zoraShield_NLMOAnalysis_write,
     &         util_file_name
c --------- definitions for shieldings NLMO analysis --- END
      external zora_getv_NMRCS_SR12,
     &         get_dia1,get_dia2par1h01,
     &         get_densZ4,read_SLCTD_NMR_Atoms,
     &         get_chi_centers_ga
      double precision ppm,scale2munu
      data ppm,scale2munu /26.62566914d+00,0.5d0/ ! /0.25d0/(clight_au*clight_au) ! 1/(4c^2)  

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
     &    call errquit('get_tensor_SR: ma failed',911,MA_ERR)
      if (.not. ma_alloc_get(mt_dbl,nat,'zan pnt',l_zanpt,k_zanpt))
     &    call errquit('get_tensor_SR: ma failed',911,MA_ERR)
c +++++++++++++++creating ga_arrays ++++++++START
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &      'get_tensor_SR: g_dia1',0,0,g_dia1))
     $    call errquit('get_tensor_SR: g_dia1',0,GA_ERR)
         call ga_zero(g_dia1)
      do i=1,3
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &      'get_tensor_SR: g_dia2',0,0,g_dia2(i)))
     $    call errquit('get_tensor_SR: g_dia2',0,GA_ERR)
         call ga_zero(g_dia2(i))
         do j=1,3
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &      'get_tensor_SR: g_dia3',0,0,g_dia3(i,j)))
     $    call errquit('get_tensor_SR: g_dia3',0,GA_ERR)
         call ga_zero(g_dia3(i,j))    
c ------------ NLMO analysis (dia2,par1) ---------- START
c        g_M --> g_dia2kt(), g_t3 --> g_par1kt()
          shldfile=0 ! not doing NLMO analysis by default
          status=rtdb_get(rtdb,'prop:shldfile',mt_int,1,shldfile) ! for NLMO analysis
          if (shldfile.eq.1) then
           if (.not. ga_create(mt_dbl,nbf,nbf,
     &        'get_tensor_SR: g_dia3',0,0,g_dia2kt(i,j)))
     $      call errquit('get_tensor_SR: g_dia2kt',0,GA_ERR)
           call ga_zero(g_dia2kt(i,j))    
           if (.not. ga_create(mt_dbl,nbf,nbf,
     &        'get_tensor_SR: g_dia3',0,0,g_par1kt(i,j)))
     $      call errquit('get_tensor_SR: g_par1kt',0,GA_ERR)
           call ga_zero(g_par1kt(i,j))    
          endif
c ------------ NLMO analysis (dia2,par1) ---------- END
         enddo ! end-loop-j
      enddo ! end-loop-i
      do i=1,6
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &     'get_tensor_SR: g_nmr1',0,0,g_nmr1(i)))
     $    call errquit('get_tensor_SR: g_nmr1',0,GA_ERR)
         call ga_zero(g_nmr1(i))
      enddo ! end-loop-i
      do i=1,18
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &     'get_tensor_SR: g_nmr2',0,0,g_nmr2(i)))
     $    call errquit('get_tensor_SR: g_nmr2',0,GA_ERR)
         call ga_zero(g_nmr2(i))
      enddo ! end-loop-i
c +++++++++++++++creating ga_arrays ++++++++END
       dims(1) =nbf
       chunk(1)=nbf 
       do i=1,3
        if (.not. nga_create(mt_dbl,1,dims,"Array v",chunk,g_v(i)))
     $    call errquit('get_tensor_SR: g_v', 0,GA_ERR)
       enddo
       call get_chi_centers_ga(g_v,ao_bas_han,nbf,geom,nat)
c +++++ Read Atom Nr for NMR calc ++
        if (.not. ga_create(mt_dbl,1,nat,
     &   'get_tensor_SR: g_AtNr',0,0,g_AtNr)) 
     $   call errquit('get_tensor_SR: g_AtNr',0,GA_ERR)
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
     &  call errquit('get_tensor_SR: nga_create failed gdia1 all',
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
     &  call errquit('get_tensor_SR: ma failed',0,MA_ERR)
      call ga_get(g_AtNr,1,1,1,nat,dbl_mb(k_AtNr),1)
c ---------- for shieldings NLMO analysis --------- START
        if (shldfile.eq.1) then
         ndir=6 ! xx,yy,zz,xy,xz,yz
         n_munu=nbf*(nbf+1)/2*ndir*nlist
         if (.not. ga_create(mt_dbl,1,n_munu,
     &        'get_tensor_SR: g_munu',0,0,g_munudia))
     $    call errquit('get_tensor_SR:',0,GA_ERR)
         call ga_zero(g_munudia)
         if (.not. ga_create(mt_dbl,1,n_munu,
     &        'get_tensor_SR: g_munu',0,0,g_munupar1))
     $    call errquit('get_tensor_SR:',0,GA_ERR)
         call ga_zero(g_munupar1)
         ndir=3 ! x,y,z
         n_munu=nbf*(nbf+1)/2*ndir*nlist
         if (.not. ga_create(mt_dbl,1,n_munu,
     &        'get_tensor_SR: g_munu',0,0,g_munu_h01))
     $    call errquit('get_tensor_SR:',0,GA_ERR)
         call ga_zero(g_munu_h01)
        endif
c ---------- for shieldings NLMO analysis --------- END
      count =1 ! for shieldings NLMO analysis
      count1=1 ! for shieldings NLMO analysis     
      count2=1 ! for shieldings NLMO analysis     
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
     &               g_dia1,g_dia2,
     &               g_dia3,   ! g_dia1kt useful for shieldings-NLMO analysis
     &               g_densZ4, 
     &               iat1,iat,
     &               nbf,ipol,xyz_NMRQcoords,g_v)
       call get_dia2par1h01(gdia2, ! out
     &                      gpar1, ! out
     &                      g_h01, ! out
     &                      iat1,iat,
     &                      ipol,nbf, 
     &                      g_nmr1,g_nmr2,
     &                      g_v,g_densZ4,
     &                      count2,    ! out: for shieldings-NLMO analysis
     &                      g_munu_h01,! out: for shieldings-NLMO analysis
     &                      g_dia2kt,  ! out: for shieldings-NLMO analysis
     &                      g_par1kt,  ! out: for shieldings-NLMO analysis
     &                      rtdb)      ! in : rtdb handle
c -- collect dia: (g_dia3,g_dia2kt) par1: g_par1kt for shieldings NLMO ---- START
       shldfile=0 ! not doing NLMO analysis by default
       status=rtdb_get(rtdb,'prop:shldfile',mt_int,1,shldfile) ! for NLMO analysis
       if (shldfile.eq.1) then
        do k=1,3
         do t=1,3
          call ga_add(scale2munu*ppm,g_dia3(k,t),
     &                scale2munu*ppm,g_dia2kt(k,t),g_dia3(k,t))
          call ga_scale(g_par1kt(k,t),-scale2munu*ppm)
         enddo
        enddo
         call fill_munuFCSD(g_munudia,  !out   : array with matrices (g_dia3--> g_munudia)
     &                      count,      !in/out: counting data stored in g_munuFCSD
     &                      g_dia3,     ! in   : nbf
     &                      nbf) 
         call fill_munuFCSD(g_munupar1, !out   : array with matrices (g_par1kt-->g_munupar1)
     &                      count1,     !in/out: counting data stored in g_munuFCSD
     &                      g_par1kt,   ! in   : nbf
     &                      nbf) 
       endif
c -- collect dia: (g_dia3,g_dia2kt) par1: g_par1kt for shieldings NLMO ---- END
      end do ! iat loop
c ---Destroying ga arrays ------- START
c ------------ NLMO analysis (dia2,par1) ---------- START
       shldfile=0 ! not doing NLMO analysis by default
       status=rtdb_get(rtdb,'prop:shldfile',mt_int,1,shldfile) ! for NLMO analysis
       if (shldfile.eq.1) then
        call ga_scale(g_munu_h01,ppm)
        ndir =6 ! xx,yy,zz,xy,xz,yz
        ndir1=3 ! x,y,z
        ndata=1 ! =1 (dia,par1,h01,dens) =2 (g_c1, perturbed MO)
        call util_file_name(lbl_nlmoshield,.false.,.false.,zorafilename)
         if (.not.dft_zoraShield_NLMOAnalysis_write(
     &       zorafilename,  ! in: filename
     &                nbf,  ! in: nr basis functions
     &               ndir,  ! in: nr of directions: 6 = xx yy zz xy xz yz for g_munudia,g_munupar1
     &              ndir1,  ! in: nr of directions: 3 = x y z             for g_munu_h01
     &              nlist,  ! in: list of selected atoms 
     &                noc,  ! in: dummy not used yet here
     &              ndata,  ! in: =1 write dia,par1,h01,dens =2 write g_c1
     &          g_munudia,  ! in: munu dia 
     &         g_munupar1,  ! in: 1st term of para
     &         g_munu_h01,  ! in: perturbed AO: H01
     &               ipol,  ! in: nr. of polarizations
     &            vectors,  ! in: dummy not used yet here
     &               g_c1,  ! in: dummy not used yet here
     &         g_densZ4(3)))! in: total density
     &   call errquit('get_tensor_SR: dft_zoraShieldNLMO_write failed',
     &                0,DISK_ERR)
          if (.not. ga_destroy(g_munu_h01)) call errquit(
     &      'gCSSR: ga_destroy failed ',0, GA_ERR) 
          if (.not. ga_destroy(g_munudia)) call errquit(
     &      'gCSSR: ga_destroy failed ',0, GA_ERR) 
          if (.not. ga_destroy(g_munupar1)) call errquit(
     &      'gCSSR: ga_destroy failed ',0, GA_ERR) 
        do i=1,3
         do j=1,3
          if (.not. ga_destroy(g_dia2kt(i,j))) call errquit(
     &      'gCSSR: ga_destroy failed ',0, GA_ERR) 
          if (.not. ga_destroy(g_par1kt(i,j))) call errquit(
     &      'gCSSR: ga_destroy failed ',0, GA_ERR) 
         enddo
        enddo
       endif
c ------------ NLMO analysis (dia2,par1) ---------- END
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
     &                    g_dia3,        ! in: ! used for NLMO analysis also
     &                    g_densZ4,      ! in:
     &                    iat,           ! in: atom nr.
     &                    iat1,          ! in: atom index in geometry
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
      integer rtdb
      integer gdia1 ! 1st term of dia tensor
      integer g_dia1,g_dia2(3),
     &        g_dia3(3,3) ! output/input
      integer g_densZ4(3),g_dia1kt(3,3)
      integer g_R(3)
      integer k,t,cbuf,ispin,ipol,iat,iat1,nbf
      integer cbuf1,cbuf2
      integer g_t1,g_t2   
      integer alo(3),ahi(3),ld(2)
      integer l_buf,k_buf
      double precision ac,n1,scale2munu,val
      double precision xyz_NMRQcoords(3),scl
      double precision ppm
      data ppm     /26.62566914d+00/ 
c ---- Creating scratch ga-arrays --------- START
       if (.not. ga_create(mt_dbl,nbf,nbf,
     &     'get_dia1: g_t1',0,0,g_t1))
     $ call errquit('get_dia1: g_t1',0,GA_ERR)
       call ga_zero(g_t1)
       if (.not. ga_create(mt_dbl,nbf,nbf,
     &     'get_dia1: g_t2',0,0,g_t2))
     $ call errquit('get_dia1: g_t2',0,GA_ERR)
       call ga_zero(g_t2)
c ---- Creating scratch ga-arrays --------- END
      if(.not.ma_alloc_get(mt_dbl,3*3,'get_dia1:buf',
     &                    l_buf,k_buf))
     $  call errquit('get_dia1: ma failed',3*3, MA_ERR)
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
c          if (do_prntNMRCS) then  ! printNMRCS-if-requested-START
           if (ga_nodeid().eq.0) then
           write(*,15) k,t,iat1,ispin,
     &                 n1*scale2munu*ppm, 
     &                 ac*scale2munu*ppm       
 15        format('dia1(k,t,iat1,ispin)=(',
     &            i3,',',i3,',',i3,',',i3,')=(',
     &            f15.8,',',f15.8,') ppm')  
           endif
c          endif
         enddo ! ispin-loop
         cbuf=k_buf+(k-1)*3+t-1
         dbl_mb(cbuf)=ac*scale2munu*ppm
         enddo ! t loop
       enddo ! k loop
c ----------- symmetrize dia1 --------- START
       do k=1,3    ! k=x,y,z
        do t=k+1,3 ! t=x,y,z    
         call ga_add(0.5d0,g_dia3(k,t),
     &               0.5d0,g_dia3(t,k),g_t1)
         call ga_copy(g_t1,g_dia3(k,t))
         call ga_copy(g_t1,g_dia3(t,k))
        enddo ! t loop
       enddo ! k loop
c ----------- symmetrize dia1 --------- END
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
     &    'get_dia1: ga_destroy failed ',0, GA_ERR)
       if (.not. ga_destroy(g_t2)) call errquit(
     &    'get_dia1: ga_destroy failed ',0, GA_ERR)
      if (.not.ma_free_heap(l_buf)) call errquit
     &   ('get_dia1: ma_free_heap of l_buf failed',
     &   911,MA_ERR)
      return
      end
      subroutine get_dia2par1h01(gdia2,     ! out
     &                           gpar1,     ! out
     &                           g_h01,     ! out
     &                           iat,       ! in: atom nr
     &                           iat1,      ! in: atom index in geometry
     &                           ipol,      ! in: nr. polarizations
     &                           nbf,       ! in: nr. basis functions
     &                           g_nmr1,    ! in:
     &                           g_nmr2,    ! in:
     &                           g_R,       ! in: R_{mu}
     &                           g_densZ4,  ! in: e-density
     &                           count,     ! io  : counter for shieldings-NLMO
     &                           g_munu_h01,! out : for shieldings-NLMO analysis
     &                           g_dia2kt,  ! out : for shieldings-NLMO analysis
     &                           g_par1kt,  ! out : for shieldings-NLMO analysis
     &                           rtdb)      ! in  : rtdb handle                         
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh" 
#include "msgids.fh"
#include "zora.fh"
      integer rtdb
      integer gdia2,gpar1
      integer g_h01
      integer g_nmr1(6),g_nmr2(18),g_R(3)
      integer g_densZ4(3)
      integer g_h1,g_t1,g_t2,g_t3,g_M,g_Nmn(3,3)
      integer iat,iat1,nbf
      integer j,k,t,a,b,cbuf,ipol
      integer cbuf1,cbuf2
      double precision ac,ad,n1,n2
      double precision scl2munu,val
      integer alo(3),ahi(3),ld(2)
      integer l_buf,k_buf,
     &        l_buf1,k_buf1,
     &        l_buf2,k_buf2
c ----- definitions for shieldings NLMO analysis ---- START
      logical status  
      integer shldfile,count,g_munu_h01
      integer g_dia2kt(3,3),g_par1kt(3,3)
      external fill_munuPSOSO
c ----- definitions for shieldings NLMO analysis ---- END
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
     &             -1.0d0,g_nmr1(t+3),g_t1) ! i.e. t=1  mn-nm=23-32, etc.
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
c ------------ NLMO analysis (dia2,par1) ---------- START
c        g_M --> g_dia2kt(), g_t3 --> g_par1kt()
        shldfile=0 ! not doing NLMO analysis by default
        status=rtdb_get(rtdb,'prop:shldfile',mt_int,1,shldfile) ! for NLMO analysis
        if (shldfile.eq.1) then
         call ga_copy(g_M ,g_dia2kt(k,t))
         call ga_copy(g_t3,g_par1kt(k,t))
        endif
c ------------ NLMO analysis (dia2,par1) ---------- END
         ac=0.0d0
         ad=0.0d0
c ---- WARNING: changing ipol for testing --- 02-17-11
          do j=1,ipol    ! original
           n1=ga_ddot(g_densZ4(j),g_t3) ! for par1
           ac=ac+n1
           n2=ga_ddot(g_densZ4(j),g_M)  ! for dia2
           ad=ad+n2
           if (ga_nodeid().eq.0) then
           write(*,25) k,t,iat1,j,
     &                 n2*scl2munu*ppm, 
     &                 ad*scl2munu*ppm,
     &                 n1*scl2munu*ppm, 
     &                 ac*scl2munu*ppm    
 25        format('(dia2,par1)(k,t,iat1,spin)=(',
     &            i3,',',i3,',',i3,',',i3,')=(',
     &            f15.8,',',f15.8,',',
     &            f15.8,',',f15.8,') ppm')  
           endif
         enddo ! end-loop-j  
         cbuf=k_buf1+(k-1)*3+t-1
         dbl_mb(cbuf)= ac*(-scl2munu*ppm) ! par1
         cbuf=k_buf2+(k-1)*3+t-1
         dbl_mb(cbuf)= ad*scl2munu*ppm    ! dia2
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

c          if (ga_nodeid().eq.0)
c     &     write(*,15) iat1,t
c 15        format('-------h01(',i3,',',i3,')--------- START')
c           call ga_print(g_h1)
c          if (ga_nodeid().eq.0)
c     &     write(*,16) iat1,t
c 16        format('-------h01(',i3,',',i3,')--------- END')

c ------ store g_h01 for NLMO analysis -------- START
c       g_h1 --> g_munu_h01
         shldfile=0 ! not doing NLMO analysis by default
         status=rtdb_get(rtdb,'prop:shldfile',mt_int,1,shldfile) ! for NLMO analysis
         if (shldfile.eq.1) then
          call fill_munuPSOSO(g_munu_h01,  !out   : array with matrices (g_t3--> g_munuPSOSO)
     &                        count,       !in/out: counting data stored in g_munuPSOSO
     &                        g_h1,        ! in:  nbf
     &                        nbf)
         endif
c ------ store g_h01 for NLMO analysis -------- END
       enddo ! end-loop-t
c ---- symmetrize par1,dia2 for NLMO analysis ----- START
        shldfile=0 ! not doing NLMO analysis by default
        status=rtdb_get(rtdb,'prop:shldfile',mt_int,1,shldfile) ! for NLMO analysis
        if (shldfile.eq.1) then
         do k=1,3    ! k=x,y,z
          do t=k+1,3 ! t=x,y,z    
           call ga_add(0.5d0,g_dia2kt(k,t),
     &                 0.5d0,g_dia2kt(t,k),g_t1)
           call ga_copy(g_t1,g_dia2kt(k,t))
           call ga_copy(g_t1,g_dia2kt(t,k))
           call ga_add(0.5d0,g_par1kt(k,t),
     &                 0.5d0,g_par1kt(t,k),g_t1)  
           call ga_copy(g_t1,g_par1kt(k,t))
           call ga_copy(g_t1,g_par1kt(t,k))
          enddo ! t loop
         enddo ! k loop
        endif
c ---- symmetrize par1,dia2 for NLMO analysis ----- END
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
        if (.not. ga_destroy(g_v(i))) call errquit(
     &    'gFij: ga_destroy failed g_v',0, GA_ERR) 
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
c -------- routines for shieldings NLMO analysis ---- START
      logical function dft_zoraShield_NLMOAnalysis_write(
     &           filename, ! in: filename
     &                nbf, ! in: nr basis functions
     &               ndir, ! in: nr of directions: 6 = xx yy zz xy xz yz for g_munudia
     &              ndir1, ! in: nr of directions: 3 = x y z for g_munupar1
     &              nlist, ! in: list of selected atoms 
     &               nocc, ! in: nocc(i) i=1,2 nr. occupations
     &              ndata, ! in: writing order =1,2
     &          g_munudia, ! in: munu dia
     &         g_munupar1, ! in: munu 1st term para
     &          g_munuh01, ! in: munu para for perturbed AO operator H10
     &               npol, ! in: nr. of polarizations
     &            vectors, ! in: MOs
     &               g_c1, ! in: perturbed MO coeffs
     &             g_dens) ! in: total density
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
#include "stdio.fh"
      character*(*) filename    ! [input] File to write to
      integer nmat,   ! = nlst*ndir  (ndir=6)
     &        nmat2,  ! = nlst*ndir1 (ndir1=3)
     &        nmat1,  ! = nbf*nbf
     &        npol
      integer g_munudia,g_munupar1,
     &        g_munuh01,g_dens,
     &        vectors(npol),g_c1      
      integer ndir,ndir1,nbf,nlst,nlist,ndata,ntot,nocc(2)
      integer unitno
      parameter (unitno = 77)
      integer l_mat ,k_mat,
     &        l_mat1,k_mat1,
     &        l_mat2,k_mat2,
     &        l_c1,k_c1,
     &        l_mo,k_mo
      integer ok,iset,i,j,alo(3),ahi(3),ld(2)
      integer inntsize

      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_sync()
      ok = 0
c     Read routines should be consistent with this
c     Write out the atomic zora corrections
      if (ga_nodeid() .eq. 0) then
c     Allocate the temporary buffer
      if (ndata.eq.1) then ! First time writing
c ++++++++ using ma_alloc_get +++++++++++++++++++ START
       nlst=nbf*(nbf+1)/2
       nmat=nlst*ndir*nlist
       if (.not. ma_alloc_get(
     &        mt_dbl,nmat,'dft_zoraNLMO_writeshield',
     &        l_mat,k_mat))
     $  call errquit('dft_zoraNLMO_writeshield: ma failed', 
     &               nmat, MA_ERR)
       nmat2=nlst*ndir1*nlist
       if (.not. ma_alloc_get(
     &        mt_dbl,nmat2,'dft_zoraNLMO_writeshield',
     &        l_mat2,k_mat2))
     $  call errquit('dft_zoraNLMO_writeshield: ma failed', 
     &               nmat2, MA_ERR)
       nmat1=nbf*nbf
       if (.not. ma_alloc_get(
     &        mt_dbl,nmat,'dft_zoraNLMO_writeshield',
     &        l_mat1,k_mat1))
     $  call errquit('dft_zoraNLMO_writeshield: ma failed', 
     &               nmat1, MA_ERR)
c ++++++++ using ma_alloc_get +++++++++++++++++++ END
c     Open the file - 1st time
        open(unitno, status='unknown', form='unformatted',
     $      file=filename, err=1000)
c     Write out the number of sets and basis functions
        write(unitno, err=1001) nbf
        write(unitno, err=1001) nlst
        write(unitno, err=1001) ndir
        write(unitno, err=1001) ndir1
        write(unitno, err=1001) nlist
        call ga_get(g_munudia,1,1,1,nmat,
     &              dbl_mb(k_mat),1)
        call swrite(unitno,dbl_mb(k_mat),nmat)
        call ga_get(g_munupar1,1,1,1,nmat,
     &              dbl_mb(k_mat),1)
        call swrite(unitno,dbl_mb(k_mat),nmat)
        call ga_get(g_munuh01,1,1,1,nmat2,
     &              dbl_mb(k_mat2),1)
        call swrite(unitno,dbl_mb(k_mat2),nmat2)
        call ga_get(g_dens,1,nbf,1,nbf,
     &              dbl_mb(k_mat1),nbf)
        call swrite(unitno,dbl_mb(k_mat1),nmat1)
c     Deallocate the temporary buffer
c ----- Using ma_free_heap ------------ START
       if (.not. ma_free_heap(l_mat))
     $  call errquit('dft_zoraNLMO_writeshield: ma free_heap failed', 
     &               911, MA_ERR)
       if (.not. ma_free_heap(l_mat1))
     $  call errquit('dft_zoraNLMO_writeshield: ma free_heap failed', 
     &               911, MA_ERR)
       if (.not. ma_free_heap(l_mat2))
     $  call errquit('dft_zoraNLMO_writeshield: ma free_heap failed', 
     &               911, MA_ERR)
c ----- Using ma_free_heap ------------ END
      else if (ndata.eq.2) then
       if (.not. ma_alloc_get(
     &        mt_dbl,nbf,'dft_zoraNLMO_writehyp',
     &        l_mo,k_mo))
     $  call errquit('dft_zoraNLMO_writehyp: l_mo failed', 
     &               nbf,MA_ERR)
        ndir = 3 ! x,y,z
        ntot = nocc(1)+nocc(2)    
        nmat = nbf*ndir*ntot
        if (.not. ma_alloc_get(
     &        mt_dbl,nmat,'dft_zoraNLMO_writeshield',l_c1,k_c1))
     $   call errquit('dft_zoraNLMO_writeshield: ma failed', 
     &                nmat, MA_ERR)
c     Open the file - 2nd time
        open(unitno, status='unknown', form='unformatted',
     $       file=filename, err=1000,position='append')
        write(unitno, err=1001) nocc(1)
        write(unitno, err=1001) nocc(2)  
        write(unitno, err=1001) ntot
        write(unitno, err=1001) nmat
        write(unitno, err=1001) npol
        write(unitno, err=1001) nbf
c ----- Add MOs in file ----- START
        do i=1,npol
         do j=1,nbf
         call dcopy(nbf,0.0d0,0,dbl_mb(k_mo),1) ! reset
         call ga_get(vectors(i),1,nbf,j,j,dbl_mb(k_mo),1)
         call swrite(unitno,dbl_mb(k_mo),nbf)         
         enddo ! end-loop-j
        enddo ! end-loop-i
c ----- Add MOs in file ----- END
        alo(1)=1
        ahi(1)=nbf
        alo(2)=1
        ahi(2)=ntot
        alo(3)=1
        ahi(3)=3
        ld(1)=nbf
        ld(2)=ntot
        call nga_get(g_c1,alo,ahi,dbl_mb(k_c1),ld)
        call swrite(unitno,dbl_mb(k_c1),nmat)
        if (.not. ma_free_heap(l_mo))
     $   call errquit('dft_zoraNLMO_writehyp: ma free_heap failed', 
     &               911, MA_ERR)
        if (.not. ma_free_heap(l_c1))
     $   call errquit('dft_zoraNLMO_writeshield: ma free_heap failed', 
     &               911, MA_ERR)
      endif
c     Close the file
       close(unitno,err=1002)
       ok = 1
      end if
c     Broadcast status to other nodes
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      call ga_sync()
      dft_zoraShield_NLMOAnalysis_write = (ok .eq. 1)
      if (ga_nodeid() .eq. 0) then
         write(6,22) filename(1:inp_strlen(filename))
 22      format(/' Wrote ZORA NLMO Gshift data to ',a/)
         call util_flush(luout)
      endif
      return
 1000 write(6,*) 'dft_zoraNLMO_writegshift: failed to open ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
 1001 write(6,*) 'dft_zoraNLMO_writegshift: failed to write ',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
 1002 write(6,*) 'dft_zoraNLMO_writegshift: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
      end
      logical function dft_zoraShield_NLMOAnalysis_read(
     &           filename, ! in: filename
     &                nbf, ! in: nr basis functions
     &               ndir, ! in: nr of directions: 6 = xx yy zz xy xz yz for g_munuFCSD
     &              ndir1, ! in: nr od directions: 3 = x y z for g_munuPSOSO
     &              nlist, ! in: list of selected atoms 
     &               nocc, ! in: nocc(i) i=1,2
     &               npol, ! in: nr polarizations
     &          g_munudia, ! in: munu dia
     &         g_munupar1, ! in: munu 1st term para
     &          g_munuh01, ! in: munu para for perturbed AO operator H10
     &            vectors, ! out: MOs
     &               g_c1, ! in: perturbed MO coeffs
     &             g_dens) ! in: total density
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

      character*(*) filename    ! [input] File to write to
      integer nmat,nmat1,nmat2,npol
      integer g_munudia,g_munupar1,
     &        g_munuh01,g_dens,
     &        vectors(npol),g_c1       
      integer nbf,nbf_read,nlst,nlst_read,
     &        ndir,ndir_read,
     &        ndir1,ndir1_read,   
     &        nlist,nlist_read,
     &        ntot,ntot_read,
     &        nocc(2),nocc_read(2),
     &        n_c1,n_c1_read,
     &        npol_read,
     &        alo(3),ahi(3),ld(2)
      integer unitno
      parameter (unitno = 77)
      integer l_mat ,k_mat,
     &        l_mat1,k_mat1,
     &        l_mat2,k_mat2,
     &        l_c1,k_c1,
     &        l_mo,k_mo
      integer ok,iset,i,j
      integer inntsize
c     Initialise to invalid MA handle
      nlst=nbf*(nbf+1)/2
      inntsize=MA_sizeof(MT_INT,1,MT_BYTE)
      call ga_sync()
      ok = 0
      nmat=nlst*ndir*nlist
       if (.not. ga_create(mt_dbl,1,nmat,
     &   'dft_zoraNLMO_read: g_munudia',0,0,g_munudia)) 
     $   call errquit('dft_zoraNLMO_read: g_munuFCSD',0,GA_ERR)
        call ga_zero(g_munudia)  
       if (.not. ga_create(mt_dbl,1,nmat,
     &   'dft_zoraNLMO_read: g_munudia',0,0,g_munupar1)) 
     $   call errquit('dft_zoraNLMO_read: g_munupar1',0,GA_ERR)
        call ga_zero(g_munupar1)  
      nmat2=nlst*ndir1*nlist
       if (.not. ga_create(mt_dbl,1,nmat2,
     &   'dft_zoraNLMO_read: g_munuPSOSO',0,0,g_munuh01)) 
     $   call errquit('dft_zoraNLMO_read: g_munuh01',0,GA_ERR)
        call ga_zero(g_munuh01)   
       if (.not. ga_create(mt_dbl,nbf,nbf,
     &   'dft_zoraNLMO_read: g_dens',0,0,g_dens)) 
     $   call errquit('dft_zoraNLMO_read: g_dens',0,GA_ERR)
        call ga_zero(g_dens)  
       do i=1,npol
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &   'dft_zoraNLMO_read: g_sdens',0,0,vectors(i))) 
     $   call errquit('dft_zoraNLMO_read: vectors',0,GA_ERR)
        call ga_zero(vectors(i))  
       enddo
       ntot=nocc(1)+nocc(2)
       n_c1=nbf*3*ntot
       alo(1) = nbf
       alo(2) = -1
       alo(3) = -1
       ahi(1) = nbf
       ahi(2) = ntot
       ahi(3) = 3
       if (.not.nga_create(MT_DBL,3,ahi,'c1 matrix',alo,g_c1)) call 
     &     errquit('dft_zoraNMR_read: nga_create failed g_c1',
     &             0,GA_ERR)
       call ga_zero(g_c1)
 
      if (ga_nodeid() .eq. 0) then
c      Print a message indicating the file being read
       write(6,22) filename(1:inp_strlen(filename))
 22    format(/' Read ZORA Shield NLMO data from ',a/)
       call util_flush(luout)
c      Open the file
       open(unitno, status='old', form='unformatted', file=filename,
     $        err=1000)
c      Read in some basics to check if they are consistent with the calculation
       read(unitno, err=1001, end=1001) nbf_read
       read(unitno, err=1001, end=1001) nlst_read
       read(unitno, err=1001, end=1001) ndir_read
       read(unitno, err=1001, end=1001) ndir1_read
       read(unitno, err=1001, end=1001) nlist_read
c      Error checks
       if ((nbf_read      .ne. nbf)   .or.
     &     (nlst_read     .ne. nlst)  .or.
     &     (ndir_read     .ne. ndir)  .or.
     &     (ndir1_read    .ne. ndir1) .or.
     &     (nlist_read    .ne. nlist)) goto 1003
c ---> ma_alloc_get: to allocate memory
c ---> ma_free_heap: to release allocated memory
       nmat=nlst*ndir*nlist
       if (.not. ma_alloc_get(mt_dbl,nmat, ! allocate memory
     &    'dft_zoraNLMO_read',l_mat,k_mat))
     $  call errquit('dft_zoraNLMO_read: ma failed', 
     &               nmat, MA_ERR)
       nmat2=nlst*ndir1*nlist
       if (.not. ma_alloc_get(mt_dbl,nmat2, ! allocate memory
     &    'dft_zoraNLMO_read',l_mat2,k_mat2))
     $  call errquit('dft_zoraNLMO_read: ma failed', 
     &               nmat2, MA_ERR)
       nmat1=nbf*nbf
       if (.not. ma_alloc_get(mt_dbl,nmat1, ! allocate memory
     &    'dft_zoraNLMO_read',l_mat1,k_mat1))
     $  call errquit('dft_zoraNLMO_read: ma failed', 
     &               nmat1, MA_ERR)
       call sread(unitno,dbl_mb(k_mat),nmat)
       call ga_put(g_munudia,1,1,1,nmat,dbl_mb(k_mat),1)
       call sread(unitno,dbl_mb(k_mat),nmat)
       call ga_put(g_munupar1,1,1,1,nmat,dbl_mb(k_mat),1)
       call sread(unitno,dbl_mb(k_mat2),nmat2)
       call ga_put(g_munuh01,1,1,1,nmat2,dbl_mb(k_mat2),1)
       call sread(unitno,dbl_mb(k_mat1),nmat1)
       call ga_put(g_dens,1,nbf,1,nbf,dbl_mb(k_mat1),nbf)
      if (.not. ma_free_heap(l_mat))       ! deallocate memory
     $  call errquit('dft_zoraNLMO_read: ma free_heap failed', 
     &               911, MA_ERR)
      if (.not. ma_free_heap(l_mat1))      ! deallocate memory
     $  call errquit('dft_zoraNLMO_read: ma free_heap failed', 
     &               911, MA_ERR)
      if (.not. ma_free_heap(l_mat2))      ! deallocate memory
     $  call errquit('dft_zoraNLMO_read: ma free_heap failed', 
     &               911, MA_ERR)
c -------- Read Perturbed MO coeffs ------------- START
       read(unitno, err=1001, end=1001) nocc_read(1)
       read(unitno, err=1001, end=1001) nocc_read(2)
       read(unitno, err=1001, end=1001) ntot_read
       read(unitno, err=1001, end=1001) n_c1_read
       read(unitno, err=1001, end=1001) npol_read
       read(unitno, err=1001, end=1001) nbf_read
c      Error checks
       if ((nocc_read(1) .ne. nocc(1)) .or.
     &     (nocc_read(2) .ne. nocc(2)) .or.
     &     (ntot_read    .ne. ntot)    .or.
     &     (nbf_read     .ne. nbf)     .or.
     &     (npol_read    .ne. npol)    .or.
     &     (n_c1_read    .ne. n_c1)) goto 1003
c ----- Read MOs ----- START
       if (.not. ma_alloc_get(
     &        mt_dbl,nbf,'dft_zoraNLMO_readhyp',
     &        l_mo,k_mo))
     $  call errquit('dft_zoraNLMO_readhyp: ma failed', 
     &               nbf,MA_ERR)
        do i=1,npol
         do j=1,nbf
          call dcopy(nbf,0.0d0,0,dbl_mb(k_mo),1) ! reset
          call sread(unitno,dbl_mb(k_mo),nbf)   
          call ga_put(vectors(i),1,nbf,j,j,dbl_mb(k_mo),1)
         enddo ! end-loop-j
        enddo ! end-loop-i
c ----- Read MOs ----- END
       if (.not. ma_alloc_get(mt_dbl,n_c1, ! allocate memory
     &    'dft_zoraNLMO_read',l_c1,k_c1))
     $  call errquit('dft_zoraNLMO_read: ma failed', 
     &               n_c1, MA_ERR)
! Note: n_c1 = nbf*3*ntot    [ ntot=nocc(1)+nocc(2) ]
        alo(1)=1
        ahi(1)=nbf
        alo(2)=1
        ahi(2)=ntot
        alo(3)=1
        ahi(3)=3
        ld(1) =nbf
        ld(2) =ntot
        call sread(unitno,dbl_mb(k_c1),n_c1)
        call nga_put(g_c1,alo,ahi,dbl_mb(k_c1),ld)
      if (.not. ma_free_heap(l_mo))       ! deallocate memory
     $  call errquit('dft_zoraNLMO_read: ma free_heap failed', 
     &               911, MA_ERR)
      if (.not. ma_free_heap(l_c1))       ! deallocate memory
     $  call errquit('dft_zoraNLMO_read: ma free_heap failed', 
     &               911, MA_ERR)
c -------- Read Perturbed MO coeffs ------------- END
c      Close the file
       close(unitno,err=1002)
       ok = 1
      end if
c
c     Broadcast status to other nodes
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      call ga_sync()
      dft_zoraShield_NLMOAnalysis_read = ok .eq. 1
      return
 1000 write(6,*) 'dft_zoraNLMO_read: failed to open',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
 1001 write(6,*) 'dft_zoraNLMO_read: failed to read',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
 1003 write(6,*) 'dft_zoraNLMO_read: file inconsistent',
     &           ' with calculation',
     $           filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      close(unitno,err=1002)
      goto 10
 1002 write(6,*) 'dft_zoraNLMO_read: failed to close',
     $     filename(1:inp_strlen(filename))
      call util_flush(luout)
      ok = 0
      goto 10
      end
c
c -------- routines for shieldings NLMO analysis ---- END
      subroutine print_NMRCS_SRZORA_version()
c
      implicit none
c
#include "stdio.fh"
c
      write(LuOut,*)
      call util_print_centered(LuOut,
     $   'Scalar Relativistic ZORA NMR Shielding', 23, .true.)
      write(LuOut,*)
c
      return
      end
c $Id$
