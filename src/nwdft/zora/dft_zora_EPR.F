      subroutine get_EPRg(rtdb,
     &                    g_dens_at,
     &                    nexc,
     &                    geom, 
     &                    ao_bas_han,
     &                    nbf,
     &                    focc,
     &                    noc,
     &                    ipol,
     &                    done_Fji,
     &                    ga_Fji,
     &                    g_densZ4)
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh" 
#include "geom.fh" 
#include "zora.fh" 
      integer ga_dia_epr,ga_para1_epr,ga_h01_epr,ga_Fji
      character*255 zorafilename
      integer type_nmrdata,nat_slc
      logical dft_zoraNMR_write,dft_zoraNMR_write_AB
      integer g_densZ4(3),g_sdens,g_dens_at(2) 
      integer rtdb    
      integer geom
      integer ao_bas_han    
      integer nbf,noc(2),ipol,nexc
      double precision focc(nbf*ipol) ! contains occupation values if modified
      integer i,j
      integer alo(3),ahi(3),ld(2)
      logical status, Knucl
      integer g_epr1(3),g_epr3(3),
     &        g_epr2(3,3),g_epr4(3,3)
      integer g_dia_epr_AB,g_par1_epr_AB 
      double precision chi_cntr(3,nbf)
      integer nat,slc_spinpolAO
      external get_par1_EPR,
     &         get_H11munu_EPR,get_R01munu_EPR,
     &         get_Fji,zora_getv_EPR,
     &         dft_zoraNMR_write,
     &         dft_zoraNMR_write_AB
      logical done_Fji,do_prntNMRCS
      if(.not.rtdb_get(rtdb,'zora:do_prntNMRCS',        ! FA
     &                 mt_log,1,do_prntNMRCS))          ! FA
     &  do_prntNMRCS= .false.      
c ---- get spin-densty: g_sdens -------- START
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &                      'EPRgpar1: g_sdens',
     $                      0,0,g_sdens))
     $       call errquit('EPR-1: g_epr1', 0,
     &                    GA_ERR)
      if(.not.rtdb_get(rtdb,'zora:slc_spinpolAO',        ! FA
     &                 mt_int,1,slc_spinpolAO))          ! FA
     &  slc_spinpolAO= 0 ! 0=A-B,1=A,2=-B selecting spin-density matrix      

       if (ga_nodeid().eq.0)
     &   write(*,11) slc_spinpolAO
 11      format('In dft_zora_EPR:: slc_spinpolAO=',i3)

            if (slc_spinpolAO .eq. 1) then
       if (ga_nodeid().eq.0) write(*,*) 'WARNING: SLC Alpha contrib'
       call ga_add( 1.0d0,g_densZ4(1),
     &              0.0d0,g_densZ4(2),g_sdens)
       else if (slc_spinpolAO .eq. 2) then
       if (ga_nodeid().eq.0) write(*,*) 'WARNING: SLC Beta-neg contrib'
       call ga_add(  0.0d0,g_densZ4(1),
     &              -1.0d0,g_densZ4(2),g_sdens)
       else if (slc_spinpolAO .eq. 0) then ! normal calc (DEFAULT)
       if (ga_nodeid().eq.0) write(*,*) 'WARNING: SLC A-B contrib'
       call ga_add( 1.0d0,g_densZ4(1),
     &             -1.0d0,g_densZ4(2),g_sdens)
       else
        write(*,12) slc_spinpolAO
 12     format('Error in dft_zora_EPR:: slc_spinpolAO=',i3,
     &         'It should be: 0,1 or 2')
        stop
       endif
c ---- get spin-densty: g_sdens -------- END

       do i=1,3
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &                      'EPRgpar1: g_epr1',
     $                      0,0,g_epr1(i)))
     $       call errquit('EPR-1: g_epr1', 0,
     &                    GA_ERR)
         call ga_zero(g_epr1(i))
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &                      'EPRgpar3: g_epr3',
     $                      0,0,g_epr3(i)))
     $       call errquit('EPR-3: g_epr3', 0,
     &                    GA_ERR)
         call ga_zero(g_epr3(i))
         do j=1,3
          if (.not. ga_create(mt_dbl,nbf,nbf,
     &                      'EPRgpar2: g_epr2',
     $                      0,0,g_epr2(i,j)))
     $       call errquit('EPR-1: g_epr2', 0,
     &                    GA_ERR)
          call ga_zero(g_epr2(i,j))
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &                      'EPRgpar4: g_epr4',
     $                      0,0,g_epr4(i,j)))
     $       call errquit('EPR-1: g_epr4', 0,
     &                    GA_ERR)
          call ga_zero(g_epr4(i,j))
         enddo ! end-loop-j
       enddo ! end-loop-i

c -------- reading basis function centers --------- START
      status=geom_ncent(geom,nat) ! out: nat
      call get_chi_centers(chi_cntr,ao_bas_han,nbf,geom,nat) ! out: g_v
c -------- reading basis function centers --------- END
c ------ Read Knucl   for including ONLY nuclear part in K ZORA ----- START
c Note.- stored in rel_input.F(rel_input(rtdb))
         Knucl=.false.
         status=rtdb_get(rtdb,'zora:Knucl',mt_log,1,Knucl) ! Check if gaussian nucl model requested
         if (ga_nodeid().eq.0)
     &     write(*,*) 'In dft_zora_EPR:: zora:Knucl=',Knucl
c ------ Read Knucl   for including ONLY nuclear part in K ZORA ----- END
      call zora_getv_EPR(rtdb,g_dens_at, 
     &                   g_epr1,g_epr2, ! out: munu matrix
     &                   g_epr3,g_epr4, ! out: munu matrix
     &                   chi_cntr,      ! in : R_nu-basis-centers
     &                   Knucl,         ! in : = T do Kzora with Nucl pot ONLY
     &                   nexc) 

       goto 111

       do i=1,3
        if (ga_nodeid().eq.0) then
         write(*,1) i
  1      format('--- g_epr1(',i3,')------ START')
        endif
        call ga_print(g_epr1(i))
        if (ga_nodeid().eq.0) then
         write(*,2) i
  2      format('--- g_epr1(',i3,')------ END')
         write(*,3) i
  3      format('--- g_epr3(',i3,')------ START')
        endif
        call ga_print(g_epr3(i))
        if (ga_nodeid().eq.0) then
         write(*,4) i
  4      format('--- g_epr3(',i3,')------ END')
        endif
         do j=1,3
          if (ga_nodeid().eq.0) then
           write(*,5) i,j
  5        format('--- g_epr2(',i3,',',i3,')------ START')
          endif
          call ga_print(g_epr2(i,j))
          if (ga_nodeid().eq.0) then
           write(*,6) i,j
  6        format('--- g_epr2(',i3,',',i3,')------ END')
          endif
          if (ga_nodeid().eq.0) then
           write(*,7) i,j
  7        format('--- g_epr4(',i3,',',i3,')------ START')
          endif
          call ga_print(g_epr4(i,j))
          if (ga_nodeid().eq.0) then
           write(*,8) i,j
  8        format('--- g_epr4(',i3,',',i3,')------ END')
          endif
         enddo ! end-loop-j
       enddo ! end-loop-i     

 111   continue
       call get_dia12_par1_EPR(ga_dia_epr,                 ! out
     &                         ga_para1_epr,               ! out
     &                         g_dia_epr_AB,               ! out
     &                         g_par1_epr_AB,              ! out
     &                         ga_h01_epr,                 ! out
     &                         g_sdens,
     &                         g_densZ4,
     &                         noc,                        ! in
     &                         focc,                       ! in
     &                         g_epr1,g_epr2,g_epr3,g_epr4,! in
     &                         geom,ao_bas_han,nbf,ipol,
     &                         lbl_nlmogshift,rtdb)        ! in
c     if (ga_nodeid().eq.0) then
c      write(*,*) '----OUT:ga_para1_epr-ppt----- START'
c     endif
c     call ga_print(ga_para1_epr)
c     if (ga_nodeid().eq.0) then
c      write(*,*) '----OUT:ga_para1_epr-ppt----- END'
c     endif
c      if (ga_nodeid().eq.0) then
c       write(*,*) '----OUT:ga_h01_epr----- START'
c      endif
c      call ga_print(ga_h01_epr)
c      if (ga_nodeid().eq.0) then
c       write(*,*) '----OUT:ga_h01_epr----- END'
c      endif
      if (.not.(done_Fji)) then
      call get_Fji(ga_Fji, ! OUT: munu-mat-Fji
     &             rtdb,g_dens_at,nexc,
     &             geom,ao_bas_han,nbf)
      endif
      if (do_prntNMRCS) then  ! printNMRCS-if-requested-START
         if (ga_nodeid().eq.0) then
           write(*,*) '----OUT:ga_Fji------ START'
         endif
         call ga_print(ga_Fji)
         if (ga_nodeid().eq.0) then
           write(*,*) '----OUT:ga_Fji------ END'
         endif
      endif                               ! printNMRCS-if-requested-END  
c ---Destroying ga arrays ------- START
         if (.not. ga_destroy(g_sdens)) call errquit(
     &     'EPRgsdens: ga_destroy failed ',0, GA_ERR)
        do i=1,3
         if (.not. ga_destroy(g_epr1(i))) call errquit(
     &     'EPRgdia1: ga_destroy failed ',0, GA_ERR)
         if (.not. ga_destroy(g_epr3(i))) call errquit(
     &     'EPRgdia3: ga_destroy failed ',0, GA_ERR)
         do j=1,3
          if (.not. ga_destroy(g_epr2(i,j))) call errquit(
     &      'EPRgdia2: ga_destroy failed ',0, GA_ERR)
          if (.not. ga_destroy(g_epr4(i,j))) call errquit(
     &      'EPRgdia4: ga_destroy failed ',0, GA_ERR)
         enddo
        enddo
c ---Destroying ga arrays ------- END
c       == get filename for the zora data ==
        nat_slc=1 ! ONLY for gshift
        type_nmrdata=3 ! =1,2,3=shieldings,hyperfine,gshift
c       Note.- lbl_nmrgshift defined in zora.fh
        call util_file_name(lbl_nmrgshift,.false.,.false.,
     &                      zorafilename)
        if (.not.dft_zoraNMR_write(zorafilename,
     &     type_nmrdata, ! =1,2,3=shieldings,hyperfine,gshift
     &     nbf,nat_slc,g_AtNr,ga_dia_epr,ga_para1_epr,
     &     ga_h01_epr,ga_Fji))
     &     call errquit('get_EPRg: dft_zoraNMR_write failed',
     &                  0,DISK_ERR)
c ----- writing dia,par AB tensors ------- START
        call util_file_name(lbl_nmrgshift_AB,.false.,.false.,
     &                      zorafilename)
        if (.not. dft_zoraNMR_write_AB(zorafilename,
     &       type_nmrdata, ! =1,2,3=shieldings,hyperfine,gshift
     &       nbf,nat_slc,g_AtNr,
     &       g_dia_epr_AB,g_par1_epr_AB))
     &   call errquit('get_EPRg: dft_zoraNMR_write failed',
     &          0,DISK_ERR)
c ----- writing dia,par AB tensors ------- END
c ---- Destroy stored ga arrays ------ START
           if (.not. ga_destroy(ga_dia_epr)) call errquit(
     &    'get_EPRg: ga_destroy failed ',0, GA_ERR)  
        if (.not. ga_destroy(ga_para1_epr)) call errquit(
     &    'get_EPRg: ga_destroy failed ',0, GA_ERR)  
        if (.not. ga_destroy(ga_h01_epr)) call errquit(
     &    'get_EPRg: ga_destroy failed ',0, GA_ERR)   
           if (.not. ga_destroy(g_dia_epr_AB)) call errquit(
     &    'get_EPRg: ga_destroy failed ',0, GA_ERR)  
        if (.not. ga_destroy(g_par1_epr_AB)) call errquit(
     &    'get_EPRg: ga_destroy failed ',0, GA_ERR)       
c ---- Destroy stored ga arrays ------ END
      return
      end

      subroutine get_dia12_par1_EPR(
     &                           g_dia_epr,     ! out
     &                           g_par1_epr,    ! out
     &                           g_dia_epr_AB , ! out : Alpha and Beta dia tensors
     &                           g_par1_epr_AB, ! out : Alpha and Beta par tensors    
     &                           g_h01_epr,     ! out
     &                           g_sdens,       ! in: spin-density
     &                           g_densZ4,      ! in: density matrices
     &                           noc,           ! in: nr. occupied states
     &                           focc,          ! in: occupation values
     &                           g_epr1,        ! in
     &                           g_epr2,        ! in
     &                           g_epr3,        ! in
     &                           g_epr4,        ! in
     &                           geom,          ! in : geometry       handle
     &                           ao_bas_han,    ! in : basis function handle
     &                           nbf,           ! in : nr. basis functions
     &                           ipol,          ! in : nr. polarizations
     &                           lbl_nlmogshift,! in : for g-shift NLMO analysis
     &                           rtdb)          ! in : rtdb handle
c Calculate 
c 1. \sum_{munu} (P_{munu,AA}-P_{munu,BB})*
c                i/(2c)<chi_mu| (r_nu x R_{numu})_k h_t^{01,EPR} |chi_nu>
c    ga output: 
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh" 
#include "geom.fh" 
      external
     &        fill_munuFCSD,fill_munuPSOSO,
     &        dft_zoraGshift_NLMOAnalysis_write,
     &        util_file_name
      integer g_densZ4(3),vectors(2)
      integer g_sdens    ! in : spin-density
      integer g_dia_epr  ! out
      integer g_par1_epr ! out
      integer g_dia_epr_AB  ! out
      integer g_par1_epr_AB ! out
      integer g_h01_epr  ! out
      integer g_dia1_epr
      integer g_dia2_epr 
      integer geom,ao_bas_han
      integer g_epr1(3),g_v(3)
      integer g_t1,g_t2,g_t3
      integer g_epr2(3,3),g_epr3(3),g_epr4(3,3)
      integer g_Lepr3(3,3),g_M_epr(3,3)
c ----- for g-shift NLMO analysis ----------- START
      character lbl_nlmogshift*(*)
      integer ndata,rtdb  
      character*255 zorafilename
      integer n_munu,ndir,ndir1,gshiftfile     
      integer g_munuEPRdia,g_munuEPRpar1 
      integer g_munuEPRHpar,g_par1A(3,3),g_c1
      logical dft_zoraGshift_NLMOAnalysis_write,
     &        dft_zoraGshift_NLMOAnalysis_read
c ----- for g-shift NLMO analysis ----------- END
      integer nbf,ipol,noc(2)
      double precision focc(nbf*ipol) ! contains occupation values if modified
      integer i,j,a,b,k,t,nat,count,count1
      integer l_buf,k_buf,cbuf,
     &        l_buf1,k_buf1,cbuf1,
     &        l_buf2,k_buf2,cbuf2,
     &        l_buf3,k_buf3,cbuf3,
     &        l_buf4,k_buf4,cbuf4,
     &        l_buf5,k_buf5,cbuf5
      integer chunk(2),dims(2)
      integer alo(3),ahi(3),ld(2)
      integer blo(3),bhi(3)
      integer ind_tmn(2,3)
      logical status
      double precision toac(2),coeffpol,toscale
      double precision ppt,val1,val2,factor
      double precision q1,q2,q3 
      integer ind,ispin
      double precision ac_occ(ipol),val

      data ind_tmn / 2, 3,  ! tmn=123
     &               3, 1,  ! tmn=231
     &               1, 2 / ! tmn=312
      data ppt     /1.0d+03/  
c ---- calculate total occupation alpha and beta --- START
       do i=1,ipol
        val=0.0d0
        ind=nbf*(i-1)
        do j=1,noc(i)         
         val=val+focc(ind+j)
        enddo
        ac_occ(i)=val
       enddo
c      if (ga_nodeid().eq.0) then
c       write(*,2) ac_occ(1),ac_occ(2),noc(1),noc(2)
c2      format('In get_dia12_par1_EPR:ac_occ=(',
c    &         f15.8,',',f15.8,')  nocc=(',
c    &         i4,',',i4,')')
c      endif
c ---- calculate total occupation alpha and beta --- END
c      if (noc(1) .ne. noc(2)) then
      if (ac_occ(1) .ne. ac_occ(2)) then
c      coeffpol=4.0d0/(noc(1)-noc(2))
       coeffpol=4.0d0/(ac_occ(1)-ac_occ(2))
      else
       write(6,1) noc(1),noc(2)
 1     format('Error in get_H11munu_EPR(): ', 
     &        'noc=(',i3,',',i3,') ',
     &        '-> closed shell system not allowed!')
       call errquit('get_dia12_par1_EPR error ',0,0)
      endif 
      factor=(-0.25d0)*coeffpol*ppt
c     if (ga_nodeid().eq.0)
c    &  write(*,*) 'coeffpol=',coeffpol

       dims(1) =nbf
       chunk(1)=nbf 
       do i=1,3
        if (.not. nga_create(mt_dbl,1,dims,"Array v",chunk,g_v(i)))
     $       call errquit('EPRgpar1: g_v', 0,
     &                    GA_ERR)
       enddo
       status=geom_ncent(geom,nat) ! out: nat
       call get_chi_centers_ga(g_v,ao_bas_han,nbf,geom,nat) ! out: g_v
c -------- Creating scratch ga-arrays ------------ START
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &                       'EPRgdia2: g_t1',
     $                       0,0,g_t1))
     $       call errquit('EPRgdia2: g_t1', 0,
     &                    GA_ERR)
         call ga_zero(g_t1)
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &                       'EPRgdia2: g_t2',
     $                       0,0,g_t2))
     $       call errquit('EPRgdia2: g_t2', 0,
     &                    GA_ERR)
         call ga_zero(g_t2)
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &                       'EPRgdia2: g_t3',
     $                       0,0,g_t3))
     $       call errquit('EPRgdia2: g_t3', 0,
     &                    GA_ERR)
         call ga_zero(g_t3)
       do i=1,3
        do j=1,3
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &                      'EPRgdia2: g_Lepr3',
     $                      0,0,g_Lepr3(i,j)))
     $       call errquit('EPRgdia2: g_Lepr3', 0,
     &                    GA_ERR)
         call ga_zero(g_Lepr3(i,j))
        if (.not. ga_create(mt_dbl,nbf,nbf,
     &                      'EPRgdia2: g_M_epr',
     $                      0,0,g_M_epr(i,j)))
     $       call errquit('EPRgdia2: g_M_epr', 0,
     &                    GA_ERR)
         call ga_zero(g_M_epr(i,j))
c---------- For g-shift NLMO analysis -------- START
        gshiftfile=0 ! not doing NLMO analysis by default
        status=rtdb_get(rtdb,'prop:gshiftfile',mt_int,1,gshiftfile) ! for NLMO analysis
        if (gshiftfile.eq.1) then
          if (.not. ga_create(mt_dbl,nbf,nbf,
     &                        'EPRgdia2: g_par1A',
     $                        0,0,g_par1A(i,j)))
     $       call errquit('EPRgdia2: g_par1A', 0,
     &                    GA_ERR)
          call ga_zero(g_par1A(i,j))
        endif
c---------- For g-shift NLMO analysis -------- END
        enddo
       enddo
c -------- Creating scratch ga-arrays ------------ END
c -------- creating g_h01 -------------------- START
      alo(1) = nbf
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbf
      ahi(2) = nbf
      ahi(3) = 3
      if (.not.nga_create(MT_DBL,3,ahi,'H01 matrix',
     &   alo,g_h01_epr)) call errquit(
     &     'get_d2p1: nga_create failed g_h01_epr',0,GA_ERR)
      call ga_zero(g_h01_epr)
      if(.not.ma_alloc_get(mt_dbl,nbf*nbf,'get_d2p1',
     &                    l_buf2,k_buf2))
     &    call errquit('gp1: ma failed',911,MA_ERR)
c ---------- For g-shift-NLMO analysis --------- START
        gshiftfile=0 ! not doing NLMO analysis by default
        status=rtdb_get(rtdb,'prop:gshiftfile',mt_int,1,gshiftfile) ! for NLMO analysis
        if (gshiftfile.eq.1) then
         ndir=3 ! xx,yy,zz,xy,xz,yz
         n_munu=nbf*(nbf+1)/2*ndir
         if (.not. ga_create(mt_dbl,1,n_munu,
     &        'get_dia12_par1_EPR: g_munu',0,0,g_munuEPRHpar))
     $    call errquit('get_dia12_par1_EPR:',0,GA_ERR)
         call ga_zero(g_munuEPRHpar)
        endif
      count1=1 ! counter for NLMO analysis
c ---------- For g-shift-NLMO analysis --------- END
      do t=1,3
       alo(1)=1
       ahi(1)=nbf
       alo(2)=1
       ahi(2)=nbf
       alo(3)=t
       ahi(3)=t
       ld(1)=nbf
       ld(2)=nbf
       call ga_get(g_epr1(t),1,nbf,1,nbf,dbl_mb(k_buf2),nbf)
       call nga_put(g_h01_epr,alo,ahi,dbl_mb(k_buf2),ld) ! store g_h01_epr
c ------------ NLMO analysis para term (PSOSO) ------------- START
        gshiftfile=0 ! not doing NLMO analysis by default
        status=rtdb_get(rtdb,'prop:gshiftfile',mt_int,1,gshiftfile) ! for NLMO analysis
        if (gshiftfile.eq.1) then
         call fill_munuPSOSO(g_munuEPRHpar, !out   : array with matrices (g_epr1--> g_munuEPRHpar)
     &                       count1,        !in/out: counting data stored in g_munuEPRHpar (dummy counter here)
     &                       g_epr1(t),     ! in:  nbf
     &                       nbf)
c WARNING: I need to scale g_munuEPRHpar by factor
        endif
c ------------ NLMO analysis para term (PSOSO) ------------- END
      enddo ! end-loop-t  
c -------- creating g_h01 -------------------- END
c -------- define g_dia2,(l_buf,k_buf) --- START
      alo(1) =  3
      alo(2) = -1
      alo(3) = -1
      ahi(1) =  3
      ahi(2) =  3
      ahi(3) =  1
      if (.not.nga_create(MT_DBL,3,ahi,'g_dia_epr matrix',
     &                    alo,g_dia_epr)) 
     &    call errquit('EPRgdia2: nga_create failed g_dia1_epr',
     &            0,GA_ERR)
      if (.not.nga_create(MT_DBL,3,ahi,'g_dia1_epr matrix',
     &                    alo,g_dia1_epr)) 
     &    call errquit('EPRgdia2: nga_create failed g_dia1_epr',
     &            0,GA_ERR)
      if (.not.nga_create(MT_DBL,3,ahi,'g_dia2_epr matrix',
     &                    alo,g_dia2_epr)) 
     &    call errquit('EPRgdia2: nga_create failed g_dia1_epr',
     &            0,GA_ERR)
      if (.not.nga_create(MT_DBL,3,ahi,'g_par1_epr matrix',
     &                    alo,g_par1_epr)) 
     &    call errquit('EPRpar1: nga_create failed g_par1_epr',
     &            0,GA_ERR)
c ------ Define AB tensors ------------ START
      alo(1) =  3
      alo(2) = -1
      alo(3) = -1
      ahi(1) =  3
      ahi(2) =  3
      ahi(3) =  2
      if (.not.nga_create(MT_DBL,3,ahi,'g_dia_epr matrix',
     &                    alo,g_dia_epr_AB)) 
     &    call errquit('EPRgdia2: nga_create failed g_dia1_eprAB',
     &            0,GA_ERR)
      if (.not.nga_create(MT_DBL,3,ahi,'g_par1_epr matrix',
     &                    alo,g_par1_epr_AB)) 
     &    call errquit('EPRpar1: nga_create failed g_par1_eprAB',
     &            0,GA_ERR)
      if(.not.ma_alloc_get(mt_dbl,2*3*3,'EPRgdia2:buf4',
     &                    l_buf4,k_buf4))
     &    call errquit('EPRgdia1: ma failed',911,MA_ERR)
      if(.not.ma_alloc_get(mt_dbl,2*3*3,'EPRgpar1:buf5',
     &                    l_buf5,k_buf5))
     &    call errquit('EPRgpar1: ma failed',911,MA_ERR)
c ------ Define AB tensors ------------ END
      if(.not.ma_alloc_get(mt_dbl,3*3,'EPRgdia2:buf1',
     &                    l_buf,k_buf))
     &    call errquit('EPRgdia1: ma failed',911,MA_ERR)
      if(.not.ma_alloc_get(mt_dbl,3*3,'EPRgpar1:buf1',
     &                    l_buf1,k_buf1))
     &    call errquit('EPRgpar1: ma failed',911,MA_ERR)
      if(.not.ma_alloc_get(mt_dbl,3*3,'EPRgdia1:buf1',
     &                    l_buf3,k_buf3))
     &    call errquit('EPRgpar1: ma failed',911,MA_ERR)
c -------- define g_dia2,(l_buf,k_buf) --- END
c ====get g_Lepr3(k,t)========================START
c ---- get dia1 ----------------------- START
      do k=1,3
       do t=1,3
        call ga_scale(g_epr4(k,t),-1.0d0)
       enddo ! end-loop-t
      enddo ! end-loop-k
      call ga_zero(g_t1)
      do k=1,3
        call ga_add(1.0d0,g_t1,1.0d0,g_epr4(k,k),g_t1)
      enddo
      do k=1,3
        call ga_add(-1.0d0,g_t1,+1.0d0,g_epr4(k,k),
     &              g_epr4(k,k))
      enddo
c ---- get dia1 ----------------------- END
      do k=1,3
       do t=1,3
        call ga_copy(g_epr3(k),g_t1)
        call get_scld_A(g_t1,g_v(t),g_t2) ! out g_t1
        call ga_scale(g_t1,-1.0d0)
        call ga_copy(g_t1,g_Lepr3(k,t))
       enddo ! end-loop-t
      enddo ! end-loop-k
      call ga_zero(g_t1)
      do k=1,3
       call ga_add(1.0d0,g_t1,1.0d0,g_Lepr3(k,k),g_t1)
      enddo
      do k=1,3
       call ga_add(-1.0d0,g_t1,+1.0d0,g_Lepr3(k,k),
     &             g_Lepr3(k,k))
      enddo
c ====get g_Lepr3(k,t)========================END
c ++++get g_M_epr(k,t)++++++++++++++++++++++++START
      do k=1,3
       a=ind_tmn(1,k)
       b=ind_tmn(2,k)
       do t=1,3   
        call ga_copy(g_epr2(a,t),g_t1)
        call get_scld_A(g_t1,g_v(b),g_t2) ! output : g_t1
        call ga_copy(g_epr2(b,t),g_t2)
        call get_scld_A(g_t2,g_v(a),g_t3) ! output : g_t2
c g_M(k)= R_{numu,b'}[N_{munu,mn,a'}-N_{munu,nm,a'}]-
c         R_{numu,a'}[N_{munu,mn,b'}-N_{munu,nm,b'}]
        call ga_add(1.0d0,g_t1,-1.0d0,g_t2,g_M_epr(k,t))
c out: g_M_epr()
        call get_3rdterm_R(g_epr1(t),g_v,
     &                     a,b,
     &                     g_t2,g_t3) ! out: g_t3  
        q2=ga_ddot(g_sdens,g_t3)*factor
c -- For g-shift NLMO analysis ---- START
        gshiftfile=0 ! not doing NLMO analysis by default
        status=rtdb_get(rtdb,'prop:gshiftfile',mt_int,1,gshiftfile) ! for NLMO analysis
        if (gshiftfile.eq.1) then
         call ga_copy(g_t3,g_par1A(k,t))
        endif
c -- For g-shift NLMO analysis ---- END
        cbuf1=k_buf1+(k-1)*3+t-1
        dbl_mb(cbuf1)=q2
c ------par-AB---- START
        do ispin=1,ipol
        cbuf4=k_buf4+(ispin-1)*9+(k-1)*3+t-1 
        dbl_mb(cbuf4)=ga_ddot(g_densZ4(ispin),g_t3)*factor
        enddo
c ------par-AB---- END
c ++++++checking (g_t3,g_M_epr,g_Lepr3) ++++++++ START
        q1=ga_ddot(g_sdens,g_M_epr(k,t))*factor
        q3=ga_ddot(g_sdens,g_Lepr3(k,t))*factor
c ++++++checking (g_t3,g_M_epr,g_Lepr3) ++++++++ END
        call ga_add(1.0d0,g_M_epr(k,t),
     &             -1.0d0,g_t3,
     &              g_M_epr(k,t))
        call ga_add(1.0d0,g_M_epr(k,t),
     &              1.0d0,g_Lepr3(k,t),
     &              g_M_epr(k,t))  
        val2=ga_ddot(g_sdens,g_M_epr(k,t))*factor
c ------- get-dia1-val1 ------- START
        val1=ga_ddot(g_sdens,g_epr4(k,t))*factor 
        cbuf3=k_buf3+(k-1)*3+t-1
        dbl_mb(cbuf3)=val1     
c ------- get-dia1-val1 ------- END
c ------dia-AB---- START
         call ga_add(1.0d0,g_M_epr(k,t),
     &               1.0d0,g_epr4(k,t),
     &               g_M_epr(k,t))    
        do ispin=1,ipol
         cbuf5=k_buf5+(ispin-1)*9+(k-1)*3+t-1 
         dbl_mb(cbuf5)=ga_ddot(g_densZ4(ispin),g_M_epr(k,t))*factor
        enddo
c ------dia-AB---- END
        if (ga_nodeid().eq.0) then
          write(*,11) k,t,
     &                q1,-q2,
     &                q1-q2, ! = aotmp2-ADF 
     &                q3,
     &                val2,val1,val2+val1
 11       format('(q1,-q2,q1-q2,q3',
     &           ',dia2,dia1,dia2+dia1)(',i3,',',i3,')=(',
     &           f12.8,',',f12.8,',',f12.8,',',
     &           f12.8,',',f12.8,',',
     &           f12.8,',',f12.8,')')
        endif
        cbuf=k_buf+(k-1)*3+t-1
        dbl_mb(cbuf)=val2
       enddo ! t loop
      enddo ! k loop
c -------- g-shift-NLMO analysis: store dia ---------------- START
        gshiftfile=0 ! not doing NLMO analysis by default
        status=rtdb_get(rtdb,'prop:gshiftfile',mt_int,1,gshiftfile) ! for NLMO analysis
        if (gshiftfile.eq.1) then
         ndir=6 ! xx,yy,zz,xy,xz,yz
         n_munu=nbf*(nbf+1)/2*ndir
         if (.not. ga_create(mt_dbl,1,n_munu,
     &        'get_dia12_par1_EPR: g_munu',0,0,g_munuEPRdia))
     $    call errquit('get_dia12_par1_EPR:',0,GA_ERR)
         call ga_zero(g_munuEPRdia)
         if (.not. ga_create(mt_dbl,1,n_munu,
     &        'get_dia12_par1_EPR: g_munu',0,0,g_munuEPRpar1))
     $    call errquit('get_dia12_par1_EPR:',0,GA_ERR)
         call ga_zero(g_munuEPRpar1)
c --- Scaling and Symmetrizing dia,par1 for NLMO analysis--- START
         toscale=factor*0.5d0 ! scaling and averaging factor
         do k=1,3
          do t=k+1,3
           call ga_add(toscale,g_M_epr(k,t),
     &                 toscale,g_M_epr(t,k),g_t1)
           call ga_copy(g_t1,g_M_epr(k,t))
           call ga_copy(g_t1,g_M_epr(t,k))
           call ga_add(toscale,g_par1A(k,t),
     &                 toscale,g_par1A(t,k),g_t1)
           call ga_copy(g_t1,g_par1A(k,t))
           call ga_copy(g_t1,g_par1A(t,k))
          enddo
         enddo
         do k=1,3
          call ga_scale(g_M_epr(k,k),factor)
          call ga_scale(g_par1A(k,k),factor)
         enddo
c --- Scaling and Symmetrizing dia,par1 for NLMO analysis--- END
         call ga_scale(g_munuEPRHpar,factor)
         count=1
         call fill_munuFCSD(g_munuEPRdia, !out   : array with matrices (g_M_epr--> g_munuEPRdia)
     &                      count,        !in/out: counting data stored in g_munuEPRdia
     &                      g_M_epr,      !   in : nbf
     &                      nbf) 
         count=1
         call fill_munuFCSD(g_munuEPRpar1,!out   : array with matrices (g_par1A--> g_munuEPRpar1)
     &                      count,       !in/out: counting data stored in g_munuEPRpar1
     &                      g_par1A,      !   in : nbf
     &                      nbf) 
c -------> write g-shift NLMO analysis data here: 
c          indata: g_munuEPRdia,g_munuEPRpar1  with ndir  = 6 : xx,yy,zz,xy,xz,yz
c                  g_munuEPRHpar               with ndir1 = 3 : x,y,z
         ndir =6  ! xx,yy,zz,xy,xz,yz
         ndir1=3  ! x,y,z
         ndata=1 !  =1 write dia,par1,sdens =2 write g_c1
         call util_file_name(lbl_nlmogshift,.false.,.false.,
     &                       zorafilename)
         if (.not.dft_zoraGshift_NLMOAnalysis_write(
     &       zorafilename, ! in: filename
     &                nbf, ! in: nr basis functions
     &               ndir, ! in: nr of directions: 6 = xx yy zz xy xz yz for g_munuFCSD
     &              ndir1, ! in: nr of directions: 3 = x y z             for g_munuPSOSO
     &                noc, ! in: dummy not used yet here
     &              ndata, ! in: =1 write dia,par1,sdens =2 write g_c1
     &       g_munuEPRdia, ! in: munu          dia 
     &      g_munuEPRpar1, ! in: munu 1st term para 
     &      g_munuEPRHpar, ! in: munu para for perturbed AO operator H10 
     &               ipol, ! in: nr. of polarizations
     &            vectors, ! in: dummy not used yet here
     &               g_c1, ! in: dummy not used yet here
     &            g_sdens))! in: spin density
     &   call errquit('get_Htensor_fast: dft_zoraHYPNLMO_write failed',
     &                0,DISK_ERR)
         do i=1,3
          do j=1,3
            if (.not. ga_destroy(g_par1A(i,j))) call errquit( 
     &          'EPRgdia2: ga_destroy failed ',0, GA_ERR)  
          enddo
         enddo
        if (.not. ga_destroy(g_munuEPRdia)) call errquit(
     &    'EPRgdia2: ga_destroy failed ',0, GA_ERR)
        if (.not. ga_destroy(g_munuEPRpar1)) call errquit(
     &    'EPRgdia2: ga_destroy failed ',0, GA_ERR)
        if (.not. ga_destroy(g_munuEPRHpar)) call errquit(
     &    'EPRgdia2: ga_destroy failed ',0, GA_ERR)
        endif
c -------- g-shift-NLMO analysis: store dia ---------------- END
c ++++get g_M_epr(k,t)++++++++++++++++++++++++END
      alo(1)=1
      ahi(1)=3
      alo(2)=1
      ahi(2)=3
      alo(3)=1
      ahi(3)=1
      ld(1)=3
      ld(2)=3
      call nga_put(g_dia1_epr,alo,ahi,dbl_mb(k_buf3),ld)
      call nga_put(g_dia2_epr,alo,ahi,dbl_mb(k_buf) ,ld)
      call nga_put(g_par1_epr,alo,ahi,dbl_mb(k_buf1),ld)
      call ga_add(1.0d0,g_dia1_epr,1.0d0,g_dia2_epr,g_dia_epr)
c ---- store par,dia AB ---------- START
      alo(1)=1
      ahi(1)=3
      alo(2)=1
      ahi(2)=3
      alo(3)=1
      ahi(3)=2
      ld(1)=3
      ld(2)=3
      call nga_put(g_par1_epr_AB,alo,ahi,dbl_mb(k_buf4),ld)
      call nga_put(g_dia_epr_AB,alo,ahi,dbl_mb(k_buf5),ld)
c ---- store par,dia AB ---------- END
      goto 10
      if (ga_nodeid().eq.0) then
       write(*,*) '----OUT:ga_dia_epr-ppt----- START'
      endif
      call ga_print(g_dia_epr)
      if (ga_nodeid().eq.0) then
       write(*,*) '----OUT:ga_dia_epr-ppt----- END'
      endif
      if (ga_nodeid().eq.0) then
       write(*,*) '----OUT:ga_dia1_epr----- START'
      endif
      call ga_print(g_dia1_epr)
      if (ga_nodeid().eq.0) then
       write(*,*) '----OUT:ga_dia1_epr----- END'
      endif 
      if (ga_nodeid().eq.0) then
       write(*,*) '----OUT:g_dia2_epr-ppt----- START'
      endif
      call ga_print(g_dia2_epr)
      if (ga_nodeid().eq.0) then
       write(*,*) '----OUT:g_dia2_epr-ppt----- END'
      endif
c ------- print dia,epr AB --------- START
      if (ga_nodeid().eq.0) then
       write(*,*) '----OUT:ga_dia_epr_AB-ppt----- START'
      endif
      call ga_print(g_dia_epr_AB)
      if (ga_nodeid().eq.0) then
       write(*,*) '----OUT:ga_dia_epr_AB-ppt----- END'
      endif
      if (ga_nodeid().eq.0) then
       write(*,*) '----OUT:ga_par1_epr_AB-ppt----- START'
      endif
      call ga_print(g_par1_epr_AB)
      if (ga_nodeid().eq.0) then
       write(*,*) '----OUT:ga_par1_epr_AB-ppt----- END'
      endif
 10   continue
c ------- print dia,epr AB --------- END
c----deallocate memory - FA
       if (.not.ma_free_heap(l_buf)) call errquit
     &    ('get_d2p1, ma_free_heap of l_buf1 failed',
     &      911,MA_ERR)
       if (.not.ma_free_heap(l_buf1)) call errquit
     &    ('get_d2p1, ma_free_heap of l_buf1 failed',
     &      911,MA_ERR)
       if (.not.ma_free_heap(l_buf2)) call errquit
     &    ('get_d2p1, ma_free_heap of l_buf1 failed',
     &      911,MA_ERR)
       if (.not.ma_free_heap(l_buf3)) call errquit
     &    ('get_d2p1, ma_free_heap of l_buf1 failed',
     &      911,MA_ERR)
       if (.not.ma_free_heap(l_buf4)) call errquit
     &    ('get_d2p1, ma_free_heap of l_buf4 failed',
     &      911,MA_ERR)
       if (.not.ma_free_heap(l_buf5)) call errquit
     &    ('get_d2p1, ma_free_heap of l_buf5 failed',
     &      911,MA_ERR)
c ---Destroying ga arrays ------- START
        if (.not. ga_destroy(g_dia1_epr)) call errquit(
     &    'EPRgdia2: ga_destroy failed ',0, GA_ERR)
        if (.not. ga_destroy(g_dia2_epr)) call errquit(
     &    'EPRgdia2: ga_destroy failed ',0, GA_ERR)
        if (.not. ga_destroy(g_t1)) call errquit(
     &    'EPRgdia2: ga_destroy failed ',0, GA_ERR)
        if (.not. ga_destroy(g_t2)) call errquit(
     &    'EPRgdia2: ga_destroy failed ',0, GA_ERR)
        if (.not. ga_destroy(g_t3)) call errquit(
     &    'EPRgdia2: ga_destroy failed ',0, GA_ERR)
       do i=1,3
         if (.not. ga_destroy(g_v(i))) call errquit(
     &    'EPRgdia2: ga_destroy failed ',0, GA_ERR)  
        do j=1,3
         if (.not. ga_destroy(g_M_epr(i,j))) call errquit(
     &    'EPRgdia2: ga_destroy failed ',0, GA_ERR)  
         if (.not. ga_destroy(g_Lepr3(i,j))) call errquit(
     &    'EPRgdia2: ga_destroy failed ',0, GA_ERR) 
        enddo ! end-loop-j
       enddo ! end-loop-i 
c ---Destroying ga arrays ------- END
      return
      end
c -------------- for g-shift NLMO analysis ------------- START
      logical function dft_zoraGshift_NLMOAnalysis_write(
     &           filename, ! in: filename
     &                nbf, ! in: nr basis functions
     &               ndir, ! in: nr of directions: 6 = xx yy zz xy xz yz for g_munuFCSD
     &              ndir1, ! in: nr of directions: 3 = x y z for g_munuPSOSO
     &               nocc, ! in: nocc(i) i=1,2 nr. occupations
     &              ndata, ! in: writing order =1,2
     &       g_munuEPRdia, ! in: munu dia
     &      g_munuEPRpar1, ! in: munu 1st term para
     &      g_munuEPRHpar, ! in: munu para for perturbed AO operator H10
     &               npol, ! in: nr. of polarizations
     &            vectors, ! in: MOs
     &               g_c1, ! in: perturbed MO coeffs
     &            g_sdens) ! in: spin density
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
      integer g_munuEPRdia,g_munuEPRpar1,
     &        g_munuEPRHpar,g_sdens,
     &        vectors(npol),g_c1    
      integer ndir,ndir1,nbf,nlst,ndata,ntot,nocc(2)
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
       nmat=nlst*ndir
       if (.not. ma_alloc_get(
     &        mt_dbl,nmat,'dft_zoraNLMO_writegshift',
     &        l_mat,k_mat))
     $  call errquit('dft_zoraNLMO_writegshift: ma failed', 
     &               nmat, MA_ERR)
       nmat2=nlst*ndir1
       if (.not. ma_alloc_get(
     &        mt_dbl,nmat2,'dft_zoraNLMO_writegshift',
     &        l_mat2,k_mat2))
     $  call errquit('dft_zoraNLMO_writegshift: ma failed', 
     &               nmat2, MA_ERR)
       nmat1=nbf*nbf
       if (.not. ma_alloc_get(
     &        mt_dbl,nmat,'dft_zoraNLMO_writegshift',
     &        l_mat1,k_mat1))
     $  call errquit('dft_zoraNLMO_writegshift: ma failed', 
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
        call ga_get(g_munuEPRdia,1,1,1,nmat,
     &              dbl_mb(k_mat),1)
        call swrite(unitno,dbl_mb(k_mat),nmat)
        call ga_get(g_munuEPRpar1,1,1,1,nmat,
     &              dbl_mb(k_mat),1)
        call swrite(unitno,dbl_mb(k_mat),nmat)
        call ga_get(g_munuEPRHpar,1,1,1,nmat2,
     &              dbl_mb(k_mat2),1)
        call swrite(unitno,dbl_mb(k_mat2),nmat2)
        call ga_get(g_sdens,1,nbf,1,nbf,
     &              dbl_mb(k_mat1),nbf)
        call swrite(unitno,dbl_mb(k_mat1),nmat1)
c     Deallocate the temporary buffer
c ----- Using ma_free_heap ------------ START
       if (.not. ma_free_heap(l_mat))
     $  call errquit('dft_zoraNLMO_writegshift: ma free_heap failed', 
     &               911, MA_ERR)
       if (.not. ma_free_heap(l_mat1))
     $  call errquit('dft_zoraNLMO_writegshift: ma free_heap failed', 
     &               911, MA_ERR)
       if (.not. ma_free_heap(l_mat2))
     $  call errquit('dft_zoraNLMO_writegshift: ma free_heap failed', 
     &               911, MA_ERR)
c ----- Using ma_free_heap ------------ END
      else if (ndata.eq.2) then
       if (.not. ma_alloc_get(
     &        mt_dbl,nbf,'dft_zoraNLMO_writegshift',
     &        l_mo,k_mo))
     $  call errquit('dft_zoraNLMO_writegshift: l_mo failed', 
     &               nbf,MA_ERR)
        ndir = 3 ! x,y,z
        ntot = nocc(1)+nocc(2)    
        nmat = nbf*ndir*ntot
        if (.not. ma_alloc_get(
     &        mt_dbl,nmat,'dft_zoraNLMO_writegshift',l_c1,k_c1))
     $   call errquit('dft_zoraNLMO_writegshift: ma failed', 
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
     $   call errquit('dft_zoraNLMO_writegshift: ma free_heap failed', 
     &               911, MA_ERR)
      endif
c     Close the file
       close(unitno,err=1002)
       ok = 1
      end if
c     Broadcast status to other nodes
 10   call ga_brdcst(Msg_Vec_Stat+MSGINT, ok, inntsize, 0) ! Propagate status
      call ga_sync()
      dft_zoraGshift_NLMOAnalysis_write = (ok .eq. 1)
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

      logical function dft_zoraGshift_NLMOAnalysis_read(
     &           filename, ! in: filename
     &                nbf, ! in: nr basis functions
     &               ndir, ! in: nr of directions: 6 = xx yy zz xy xz yz for g_munuFCSD
     &              ndir1, ! in: nr od directions: 3 = x y z for g_munuPSOSO
     &               nocc, ! in: nocc(i) i=1,2
     &               npol, ! in: nr polarizations
     &       g_munuEPRdia, ! out: munu dia
     &      g_munuEPRpar1, ! out: munu 1st term para
     &      g_munuEPRHpar, ! out: munu para for perturbed AO operator H10
     &            vectors, ! out: MOs
     &               g_c1, ! out: perturbed MO coeffs
     &            g_sdens) ! out: spin density
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
      integer g_munuEPRdia,g_munuEPRpar1,
     &        g_munuEPRHpar,g_sdens,
     &        vectors(npol),g_c1
      integer nbf,nbf_read,nlst,nlst_read,
     &        ndir,ndir_read,
     &        ndir1,ndir1_read,   
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
      nmat=nlst*ndir
       if (.not. ga_create(mt_dbl,1,nmat,
     &   'dft_zoraNLMO_read: g_munuEPRdia',0,0,g_munuEPRdia)) 
     $   call errquit('dft_zoraNLMO_read: g_munuEPRdia',0,GA_ERR)
        call ga_zero(g_munuEPRdia)  
      nmat=nlst*ndir
       if (.not. ga_create(mt_dbl,1,nmat,
     &   'dft_zoraNLMO_read: g_munuEPRpar1',0,0,g_munuEPRpar1)) 
     $   call errquit('dft_zoraNLMO_read: g_munuEPRpar1',0,GA_ERR)
        call ga_zero(g_munuEPRpar1)  
      nmat2=nlst*ndir1
       if (.not. ga_create(mt_dbl,1,nmat2,
     &   'dft_zoraNLMO_read: g_munuEPRHpar',0,0,g_munuEPRHpar)) 
     $   call errquit('dft_zoraNLMO_read: g_munuEPRHpar',0,GA_ERR)
        call ga_zero(g_munuEPRHpar)   
       if (.not. ga_create(mt_dbl,nbf,nbf,
     &   'dft_zoraNLMO_read: g_sdens',0,0,g_sdens)) 
     $   call errquit('dft_zoraNLMO_read: g_sdens',0,GA_ERR)
        call ga_zero(g_sdens)  
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
 22    format(/' Read ZORA Gshift NLMO data from ',a/)
       call util_flush(luout)
c      Open the file
       open(unitno, status='old', form='unformatted', file=filename,
     $        err=1000)
c      Read in some basics to check if they are consistent with the calculation
       read(unitno, err=1001, end=1001) nbf_read
       read(unitno, err=1001, end=1001) nlst_read
       read(unitno, err=1001, end=1001) ndir_read
       read(unitno, err=1001, end=1001) ndir1_read
c      Error checks
       if ((nbf_read      .ne. nbf)   .or.
     &     (nlst_read     .ne. nlst)  .or.
     &     (ndir_read     .ne. ndir)  .or.
     &     (ndir1_read    .ne. ndir1)) goto 1003
c ---> ma_alloc_get: to allocate memory
c ---> ma_free_heap: to release allocated memory
       nmat=nlst*ndir
       if (.not. ma_alloc_get(mt_dbl,nmat, ! allocate memory
     &    'dft_zoraNLMO_read',l_mat,k_mat))
     $  call errquit('dft_zoraNLMO_read: ma failed', 
     &               nmat, MA_ERR)
       nmat2=nlst*ndir1
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
       call ga_put(g_munuEPRdia,1,1,1,nmat,dbl_mb(k_mat),1)
       call sread(unitno,dbl_mb(k_mat),nmat)
       call ga_put(g_munuEPRpar1,1,1,1,nmat,dbl_mb(k_mat),1)
       call sread(unitno,dbl_mb(k_mat2),nmat2)
       call ga_put(g_munuEPRHpar,1,1,1,nmat2,dbl_mb(k_mat2),1)
       call sread(unitno,dbl_mb(k_mat1),nmat1)
       call ga_put(g_sdens,1,nbf,1,nbf,dbl_mb(k_mat1),nbf)
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
      dft_zoraGshift_NLMOAnalysis_read = ok .eq. 1
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
c -------------- for g-shift NLMO analysis ------------- END
      subroutine print_EPRg_version()
c
      implicit none
c
#include "stdio.fh"
c
      write(LuOut,*)
      call util_print_centered(LuOut,
     $   'EPR GShift', 23, .true.)
      write(LuOut,*)
c
      return
      end
c
c $Id$
