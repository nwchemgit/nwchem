      subroutine get_EPRg(rtdb,
     &                    g_dens_at,
     &                    nexc,
     &                    geom, 
     &                    ao_bas_han,
     &                    nbf,
     &                    noc,
     &                    ipol)
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
c  ga_dia_epr   ! OUTPUT
c  ga_para1_epr ! OUTPUT
c  ga_h01_epr   ! OUTPUT
c  ga_Fji       ! OUTPUT
c ------ defined in zora.fh ---- END
      integer g_dia1_epr,g_dia2_epr    
      integer g_densZ4(3),g_dens_at(2) 
      integer rtdb    
      integer geom
      integer ao_bas_han    
      integer nbf,noc(2),ipol,nexc
      integer i,j,nat_slc
      integer alo(3),ahi(3),ld(2)
      logical status  
      integer g_epr1(3),g_epr3(3),
     &        g_epr2(3,3),g_epr4(3,3)
      external get_densZ4,get_par1_EPR,
     &         get_H11munu_EPR,get_R01munu_EPR,
     &         get_Fji

      call get_densZ4(rtdb,ao_bas_han,geom,g_densZ4)

c +++++++++ TEST : ev. all EPR AOs +++ START
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
      if (ga_nodeid().eq.0)
     &    write(*,*) 'BEF zora_getv_EPR:'
      call zora_getv_EPR(rtdb,g_dens_at, 
     &                   g_epr1,g_epr2, ! out: munu matrix
     &                   g_epr3,g_epr4, ! out: munu matrix
     &                   nexc) 
c      goto 111

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

c === allocate arrays to store diamagnetic tensor === START
      alo(1) =  3
      alo(2) = -1
      alo(3) = -1
      ahi(1) =  3
      ahi(2) =  3
      ahi(3) =  1
      if (.not.nga_create(MT_DBL,3,ahi,'g_dia_epr matrix',alo,
     &    ga_dia_epr)) 
     &    call errquit('get_EPRg_SRZ4: nga_create failed g_dia',
     &                 0,GA_ERR)
c === allocate arrays to store diamagnetic tensor === END      
       if (ga_nodeid().eq.0) 
     &     write(*,*) "==== get_H11munu_EPR ======= START"
       call get_H11munu_EPR(g_dia1_epr, ! out
     &                      g_epr3,g_epr4,g_densZ4,noc,   
     &                      geom,ao_bas_han,nbf,ipol)   
      if (ga_nodeid().eq.0) then
       write(*,*) '----OUT:ga_dia1_epr----- START'
      endif
      call ga_print(g_dia1_epr)
      if (ga_nodeid().eq.0) then
       write(*,*) '----OUT:ga_dia1_epr----- END'
      endif 
       if (ga_nodeid().eq.0) 
     &     write(*,*) "==== get_H11munu_EPR ======= END"

      if (ga_nodeid().eq.0) 
     &     write(*,*) "==== get_H01_EPR ======= START"
       call get_par1_EPR(ga_para1_epr,             ! out
     &                   ga_h01_epr,               ! out
     &                   g_epr1,
     &                   g_densZ4,noc,             ! in
     &                   geom,ao_bas_han,nbf,ipol) ! in

      if (ga_nodeid().eq.0) then
       write(*,*) '----OUT:ga_para1_epr-ppm----- START'
      endif
      call ga_print(ga_para1_epr)
      if (ga_nodeid().eq.0) then
       write(*,*) '----OUT:ga_para1_epr-ppm----- END'
      endif
      if (ga_nodeid().eq.0) then
       write(*,*) '----OUT:ga_h01_epr----- START'
      endif
      call ga_print(ga_h01_epr)
      if (ga_nodeid().eq.0) then
       write(*,*) '----OUT:ga_h01_epr----- END'
      endif

      if (ga_nodeid().eq.0) 
     &     write(*,*) "==== get_H01_EPR ======= END"

       if (ga_nodeid().eq.0) 
     &     write(*,*) "==== get_R01munu_EPR ======= START"
       call get_R01munu_EPR(g_dia2_epr,               ! out
     &                      g_densZ4,noc,             ! in
     &                      g_epr1,g_epr2,g_epr3, ! in
     &                      geom,ao_bas_han,nbf,ipol) ! in

      call ga_add(1.0d0,g_dia1_epr,1.0d0,g_dia2_epr,ga_dia_epr)
      if (ga_nodeid().eq.0) then
       write(*,*) '----OUT:ga_dia_epr-ppm----- START'
      endif
      call ga_print(ga_dia_epr)
      if (ga_nodeid().eq.0) then
       write(*,*) '----OUT:ga_dia_epr-ppm----- END'
      endif
      if (ga_nodeid().eq.0) then
       write(*,*) '----OUT:g_dia2_epr-ppm----- START'
      endif
      call ga_print(g_dia2_epr)
      if (ga_nodeid().eq.0) then
       write(*,*) '----OUT:g_dia2_epr-ppm----- END'
      endif
      if (ga_nodeid().eq.0) 
     &    write(*,*) "==== get_R01munu_EPR ======= END"

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
c ---Destroying ga arrays ------- START
        do i=1,3
         if (.not. ga_destroy(g_epr1(i))) call errquit(
     &     'EPRgdia1: ga_destroy failed ',0, GA_ERR)
         if (.not. ga_destroy(g_epr3(i))) call errquit(
     &     'EPRgdia1: ga_destroy failed ',0, GA_ERR)
         do j=1,3
          if (.not. ga_destroy(g_epr2(i,j))) call errquit(
     &      'EPRgdia1: ga_destroy failed ',0, GA_ERR)
          if (.not. ga_destroy(g_epr4(i,j))) call errquit(
     &      'EPRgdia1: ga_destroy failed ',0, GA_ERR)
         enddo
        enddo
c ---Destroying ga arrays ------- END
      return
      end

      subroutine get_H11munu_EPR(g_dia1_epr,   ! out
     &                           g_epr3,g_epr4,! in
     &                           g_densZ4,     ! in: e-density
     &                           noc,          ! in: occ nr.
     &                           geom,
     &                           ao_bas_han,     
     &                           nbf,          ! in: nr. basis functions
     &                           ipol)         ! in: nr. polarizations
c Calculate 1st term of diamagnetic EPR tensor, g_{kt}^{d,1}
c g_{kt}^{d,1}= \sum_{mu nu} (P_{mu nu,AA}-P_{mu nu,BB}) <chi_mu|h_{kt,nu}^{11,EPR}|chi_{nu}>
c h_{kt,nu}^{11,EPR}=i \delta_{kt} \sum_n p_n K/(4c) r_n-r_n K/(4c) p_n -
c                                        (p_k K/(4c) r_t-r_t K/(4c) p_k)
c = g_dia1_epr
c FA-11-02-10

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh" 
#include "geom.fh" 
      integer g_dia1_epr ! output : 1st term of diamagnetic g-tensor
      integer g_epr3(3),g_epr4(3,3)
      integer g_densZ4(3)
      integer noc(2)
      integer nbf,ipol,nat
      integer ao_bas_han,geom
      integer g_tmp,g_tmp1,g_v(3)
      integer i,j,k,t
      integer l_buf,k_buf,cbuf
      integer alo(3),ahi(3),ld(2)
      integer ind_tmn(2,3)
      integer dims(3),chunk(3)
      double precision coeffpol
      double precision toac(2),ppm
      logical status
      external get_chi_centers_ga
      data ind_tmn / 2, 3,  ! tmn=123
     &               3, 1,  ! tmn=231
     &               1, 2 / ! tmn=312    
      data ppm     /26.62566914d+00/  
      if (noc(1) .ne. noc(2)) then
      coeffpol=4.0d0/(noc(1)-noc(2))
      else
       write(*,1) noc(1),noc(2)
 1     format('Error in get_H11munu_EPR(): ', 
     &        'noc=(',i3,',',i3,') ',
     &        '-> closed shell system not allowed!')
       stop
      endif 
      write(*,*) 'coeffpol=',coeffpol
c ------------ get g_v --------- START       
       dims(1) =nbf
       chunk(1)=nbf 
       do i=1,3
        if (.not. nga_create(mt_dbl,1,dims,"Array v",chunk,g_v(i)))
     $       call errquit('EPRgpar1: g_v', 0,
     &                    GA_ERR)
       enddo
       status=geom_ncent(geom,nat) ! out: nat
       call get_chi_centers_ga(g_v,ao_bas_han,nbf,geom,nat) ! out: g_v
c ------------ get g_v --------- END
        if (.not. ga_create(mt_dbl,nbf,nbf,
     &                      'EPRgdia1: g_tmp',
     $                      0,0,g_tmp))
     $       call errquit('EPRgdia1: g_tmp', 0,
     &                    GA_ERR)
        call ga_zero(g_tmp)
        if (.not. ga_create(mt_dbl,nbf,nbf,
     &                      'EPRgdia1: g_tmp1',
     $                      0,0,g_tmp1))
     $       call errquit('EPRgdia1: g_tmp1', 0,
     &                    GA_ERR)
        call ga_zero(g_tmp1)
      alo(1) =  3
      alo(2) = -1
      alo(3) = -1
      ahi(1) =  3
      ahi(2) =  3
      ahi(3) =  1
      if (.not.nga_create(MT_DBL,3,ahi,'g_dia1_epr matrix',
     &                    alo,g_dia1_epr)) 
     &    call errquit('EPRgdia1: nga_create failed g_dia1_epr',
     &            0,GA_ERR)
      if(.not.MA_push_get(MT_DBL,3*3,'EPRgdia1:buf1',
     &                    l_buf,k_buf))
     $     call errquit('EPRgdia1: ma failed',
     &                  3*3, MA_ERR)
      do k=1,3
       do t=1,3
        call ga_copy(g_epr3(k),g_tmp)
        call ga_transpose(g_tmp,g_tmp1)
        call ga_add(1.0d0,g_tmp,1.0d0,g_tmp1,
     &              g_tmp)               
        call ga_scale_cols(g_tmp,g_v(t)) ! R_{nu} <chi_mu|K|dchi_nu/dr_n>   
        call ga_transpose(g_epr4(k,t),g_tmp1)
        call ga_add(1.0d0,g_epr4(k,t),1.0d0,g_tmp1,
     &              g_epr4(k,t))    
        call ga_add(-1.0d0,g_epr4(k,t),+1.0d0,g_tmp,
     &              g_epr4(k,t))        
       enddo ! end-loop-t
      enddo ! end-loop-k
      call ga_zero(g_tmp)
      do k=1,3
        call ga_add(1.0d0,g_tmp,1.0d0,g_epr4(k,k),g_tmp)
      enddo
      do k=1,3
        call ga_add(-1.0d0,g_tmp,+1.0d0,g_epr4(k,k),
     &              g_epr4(k,k))
      enddo
      do k=1,3
       do t=1,3   
c       g_epr4(k,t)= B_{mu nu,kt} + Transpose(B_{mu nu,kt})
c       B_{mu nu,kt}=\delta_{kt} \sum_n <chi_{mu}| K (r-R_{nu})_n|dchi_{nu}/dr_n>-
c                    <chi_{mu}| K (r-R_{nu})_t|dchi_{nu}/dr_k>
        do j=1,ipol
         toac(j)=ga_ddot(g_densZ4(j),g_epr4(k,t)) 
        enddo ! end-loop-spin
        cbuf=k_buf+(k-1)*3+t-1
        dbl_mb(cbuf)=(toac(1)-toac(2))*(-0.25d0)*
     &               coeffpol ! *ppt              
           if (ga_nodeid().eq.0) then
            write(*,16) k,t,dbl_mb(cbuf)
 16         format('gdia1-EPR(',i3,',',i3,')=',f15.8)
           endif      
       enddo ! end-loop-t
      enddo ! end-loop-k
        alo(1)=1
        ahi(1)=3
        alo(2)=1
        ahi(2)=3
        alo(3)=1
        ahi(3)=1
        ld(1)=3
        ld(2)=3
        call nga_put(g_dia1_epr,alo,ahi,dbl_mb(k_buf),ld)
        if (ga_nodeid().eq.0) 
     &     write(*,*) '----- g_dia1_epr ---- START'
        call ga_print(g_dia1_epr)
        if (ga_nodeid().eq.0) 
     &     write(*,*) '----- g_dia1_epr ---- END'
c ----- working out off diagonal elements (k ne t) ---- END
c----deallocate memory - FA
      if (.not. MA_pop_stack(l_buf)) call errquit
     $     ('EPRgdia1: pop failed', 0, GA_ERR)
c ---Destroying ga arrays ------- START
        if (.not. ga_destroy(g_tmp)) call errquit(
     &    'EPRgdia1: ga_destroy failed ',0, GA_ERR)
        if (.not. ga_destroy(g_tmp1)) call errquit(
     &    'EPRgdia1: ga_destroy failed ',0, GA_ERR)
c ---Destroying ga arrays ------- END
      return
      end

      subroutine get_par1_EPR(g_par1_epr,  ! out
     &                        g_h01_epr,   ! out
     &                        g_epr1,      ! in
     &                        g_densZ4,    ! in: e-density
     &                        noc,         ! in: occ. nr.
     &                        geom,
     &                        ao_bas_han,
     &                        nbf,
     &                        ipol)
c Calculate 
c 1. (R_mu,a R_nu,b-R_mu,b R_nu,a) <chi_mu| h_t^{01,EPR}|chi_nu>
c    ga output: 
c 2. <chi_mu| h_t^{01,EPR}|chi_nu>
c    ga_output:
c kab=123,231,312
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh" 
#include "geom.fh" 
#include "zora.fh" 
      external get_chi_centers_ga,get_3rdterm_R 
      integer g_par1_epr,g_epr1(3),g_h01_epr
      integer g_v(3),g_t1,g_t2
      integer geom,ao_bas_han,nat
      integer g_densZ4(3)
      integer i,j,nbf,ipol,noc(2)
      integer l_buf,k_buf,cbuf,
     &        l_buf1,k_buf1
      integer alo(3),ahi(3),ld(2)
      integer chunk(2),dims(2)
      double precision coeffpol
      integer a,b,k,t,ind_tmn(2,3)
      double precision toac(2),ppm
      data ppm     /26.62566914d+00/   
      data ind_tmn / 2, 3,  ! tmn=123
     &               3, 1,  ! tmn=231
     &               1, 2 / ! tmn=312
      logical status
      if (noc(1) .ne. noc(2)) then
      coeffpol=4.0d0/(noc(1)-noc(2))
      else
       write(*,1) noc(1),noc(2)
 1     format('Error in get_par1_EPR(): ', 
     &        'noc=(',i3,',',i3,') ',
     &        '-> closed shell system not allowed!')
       stop
      endif 
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
     &                    l_buf1,k_buf1))
     &    call errquit('gp1: ma failed',911,MA_ERR)
      do t=1,3
       alo(1)=1
       ahi(1)=nbf
       alo(2)=1
       ahi(2)=nbf
       alo(3)=t
       ahi(3)=t
       ld(1)=nbf
       ld(2)=nbf
       call ga_get(g_epr1(t),1,nbf,1,nbf,dbl_mb(k_buf1),nbf)
       call nga_put(g_h01_epr,alo,ahi,dbl_mb(k_buf1),ld) ! store g_h01_epr
      enddo
c --------- creating scratch ga-arrays ------- START
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &                       'EPRgpar1: g_t1',
     $                       0,0,g_t1))
     $       call errquit('EPRgpar1: g_t1', 0,
     &                    GA_ERR)
         call ga_zero(g_t1)
         if (.not. ga_create(mt_dbl,nbf,nbf,
     &                       'EPRgpar1: g_t2',
     $                       0,0,g_t2))
     $       call errquit('EPRgpar1: g_t2', 0,
     &                    GA_ERR)
         call ga_zero(g_t2)
c --------- creating scratch ga-arrays ------- END
      alo(1) =  3
      alo(2) = -1
      alo(3) = -1
      ahi(1) =  3
      ahi(2) =  3
      ahi(3) = 1
      if (.not.nga_create(MT_DBL,3,ahi,'g_par1_epr matrix',
     &                    alo,g_par1_epr)) 
     &    call errquit('EPRgpar1: nga_create failed g_par1_epr',
     &            0,GA_ERR)
      if(.not.MA_push_get(MT_DBL,3*3,'EPRgpar1:buf1',
     &                    l_buf,k_buf))
     $     call errquit('EPRgpar1: ma failed',
     &                  3*3, MA_ERR)
       dims(1) =nbf
       chunk(1)=nbf 
       do i=1,3
        if (.not. nga_create(mt_dbl,1,dims,"Array v",chunk,g_v(i)))
     $       call errquit('EPRgpar1: g_v', 0,
     &                    GA_ERR)
       enddo
       status=geom_ncent(geom,nat) ! out: nat
       call get_chi_centers_ga(g_v,ao_bas_han,nbf,geom,nat) ! out: g_v
       do k=1,3                 
        do t=1,3
         a=ind_tmn(1,k)
         b=ind_tmn(2,k)
         do j=1,ipol    
          call get_3rdterm_R(g_epr1(t),g_v,
     &                       a,b,g_t1,g_t2)   ! out: g_t3  
          toac(j)=ga_ddot(g_densZ4(j),g_t2) 
         enddo ! end-loop-j  
         cbuf=k_buf+(k-1)*3+t-1
         dbl_mb(cbuf)=(toac(1)-toac(2))*(-0.25d0)*
     &                coeffpol
         if (ga_nodeid().eq.0) then
          write(*,16) k,t,dbl_mb(cbuf)
 16       format('gpara1-EPR(',i3,',',i3,')=',f15.8)
         endif     
        enddo ! end-loop-t
       enddo ! end-loop-k
       alo(1)=1
       ahi(1)=3
       alo(2)=1
       ahi(2)=3
       alo(3)=1
       ahi(3)=1
       ld(1)=3
       ld(2)=3
       call nga_put(g_par1_epr,alo,ahi,dbl_mb(k_buf),ld)
c----deallocate memory - FA
      if (.not. MA_pop_stack(l_buf)) call errquit
     $     ('EPRgdia1: pop failed', 0, GA_ERR)
       if (.not.ma_free_heap(l_buf1)) call errquit
     &    ('get_d2p1, ma_free_heap of l_buf1 failed',
     &      911,MA_ERR)
c ---Destroying ga arrays ------- START
        if (.not. ga_destroy(g_t1)) call errquit(
     &    'EPRgdia2: ga_destroy failed ',0, GA_ERR)
        if (.not. ga_destroy(g_t2)) call errquit(
     &    'EPRgdia2: ga_destroy failed ',0, GA_ERR)
      return
      end

      subroutine get_R01munu_EPR(g_dia2_epr,    ! out
     &                           g_densZ4,noc,  ! in
     &                           g_h01_epr,     ! in
     &                           g_epr2,g_epr3, ! in
     &                           geom,
     &                           ao_bas_han,
     &                           nbf,
     &                           ipol)
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
      external get_chi_centers_ga
      integer g_dia2_epr ! out
      integer geom,ao_bas_han
      integer g_densZ4(3),noc(2)
      integer g_h01_epr(3)
      integer g_v(3)
      integer g_t1,g_t2,g_t3
      integer g_epr2(3,3),g_epr3(3)
      integer g_Lepr3(3,3)
      integer g_M_epr(3,3)
      integer nbf,ipol
      integer i,j,a,b,k,t,nat
      integer l_buf,k_buf,cbuf
      integer chunk(2),dims(2)
      integer alo(3),ahi(3),ld(2)
      integer ind_tmn(2,3)
      logical status
      double precision toac(2),coeffpol
      double precision ppm
      data ind_tmn / 2, 3,  ! tmn=123
     &               3, 1,  ! tmn=231
     &               1, 2 / ! tmn=312
      data ppm     /26.62566914d+00/   
      if (noc(1) .ne. noc(2)) then
      coeffpol=4.0d0/(noc(1)-noc(2))
      else
       write(*,1) noc(1),noc(2)
 1     format('Error in get_H11munu_EPR(): ', 
     &        'noc=(',i3,',',i3,') ',
     &        '-> closed shell system not allowed!')
       stop
      endif 
      write(*,*) 'coeffpol=',coeffpol
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
        enddo
       enddo
c -------- Creating scratch ga-arrays ------------ END
c -------- define g_dia2,(l_buf,k_buf) --- START
      alo(1) =  3
      alo(2) = -1
      alo(3) = -1
      ahi(1) =  3
      ahi(2) =  3
      ahi(3) =  1
      if (.not.nga_create(MT_DBL,3,ahi,'g_dia1_epr matrix',
     &                    alo,g_dia2_epr)) 
     &    call errquit('EPRgdia2: nga_create failed g_dia1_epr',
     &            0,GA_ERR)
      if(.not.MA_push_get(MT_DBL,3*3,'EPRgdia1:buf1',
     &                    l_buf,k_buf))
c -------- define g_dia2,(l_buf,k_buf) --- END
     $     call errquit('EPRgdia1: ma failed',
     &                  3*3, MA_ERR)
       dims(1) =nbf
       chunk(1)=nbf 
       do i=1,3
        if (.not. nga_create(mt_dbl,1,dims,"Array v",chunk,g_v(i)))
     $       call errquit('EPRgpar1: g_v', 0,
     &                    GA_ERR)
       enddo
       status=geom_ncent(geom,nat) ! out: nat
       call get_chi_centers_ga(g_v,ao_bas_han,nbf,geom,nat) ! out: g_v
c ====get g_Lepr3(k,t)========================START
      do k=1,3
       do t=1,3
        call ga_copy(g_epr3(k),g_t1)
        call get_scld_A(g_t1,g_v(t),g_t2) ! out g_t1
        call ga_scale(g_t1,-1.0d0)
        call ga_copy(g_t1,g_Lepr3(k,t))
       enddo
      enddo
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
        call get_3rdterm_R(g_h01_epr(t),g_v,
     &                     a,b,
     &                     g_t2,g_t3) ! out: g_t3  
        call ga_add(1.0d0,g_M_epr(k,t),1.0,g_t3,
     &              g_M_epr(k,t))
        call ga_add(1.0d0,g_M_epr(k,t),
     &              1.0d0,g_Lepr3(k,t),
     &              g_M_epr(k,t))  
       enddo ! t loop
      enddo ! k loop
c ++++get g_M_epr(k,t)++++++++++++++++++++++++END
      do k=1,3
       do t=1,3
        do j=1,ipol
            toac(j)=ga_ddot(g_densZ4(j),g_M_epr(k,t))
        enddo ! end-loop-j
        cbuf=k_buf+(k-1)*3+t-1
        dbl_mb(cbuf)=(toac(1)-toac(2))*(-0.25d0)*
     &                coeffpol   
        if (ga_nodeid().eq.0) then
         write(*,16) k,t,dbl_mb(cbuf)
 16      format('gdia2-EPR(',i3,',',i3,')=',f15.8)
        endif      
       enddo  ! end-loop-t
      enddo ! end-loop-k
      alo(1)=1
      ahi(1)=3
      alo(2)=1
      ahi(2)=3
      alo(3)=1
      ahi(3)=1
      ld(1)=3
      ld(2)=3
      call nga_put(g_dia2_epr,alo,ahi,dbl_mb(k_buf),ld)
c----deallocate memory - FA
      if (.not. MA_pop_stack(l_buf)) call errquit
     $     ('EPRgdia1: pop failed', 0, GA_ERR)
c ---Destroying ga arrays ------- START
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

       subroutine print_EPRg_version()
       implicit none
       write(*,*) "+++++++++++++++++++++++++++++++++++++"
       write(*,*) "+ EPR-gshift Version 11-10-10 12:00 +"
       write(*,*) "+ Author: Fredy W. Aquino           +"
       write(*,*) "+++++++++++++++++++++++++++++++++++++"
       return
       end
