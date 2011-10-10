       subroutine dft_scaleMO_so(rtdb,           ! IN
     &                           g_moso,         ! IN     : MO vectors
     &                           occ,            ! IN/OUt : occupancies
     &                           g_densso,       ! OUT    : spin-orbit density matrix
     &                           nbf_mo,         ! IN     : nr. basis functions
     &                           nTotOcc,        ! IN     : nr. occ
     &                           switch_sclMO_so)! OUT    : switch 1,0=ON,OFF
c      Purpose: Fractional occupation routine
c               entered from input script:
c               occupations nA nB
c               occA-1 occB-1
c               ...
c               occA-X  occB-nB
c               occA-nA
c      Author : Fredy Aquino
c      Date   : 02-15-11
       implicit none
#include "errquit.fh"
       integer rtdb
#include "mafdecls.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
       double precision occ(nbf_mo) ! occupancies
       integer switch_sclMO_so ! switch 1,0=ON,OFF scaling MOs
                               !                   with occupations
       integer g_densso(2)     ! spin-orbit density matrix
       integer scale_switch    ! input
       integer nbf_mo          ! input
       integer g_moso(2)       ! input/output
       integer g_moso1(2) ! scratch ga array
       integer i,ispin,iorb,iorb1,nTotOcc,nbf_ao
       logical status
       integer l_Iocc,k_Iocc,nfocc_in(2), 
     &         nfocc_tot,ind_min,ind_max               
       double precision ac_occ,             
     &                  diff_occ,error_Iocc 
       external dft_densm_so
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +++++++ test: occupations keyword -- FA-02-10-11 ++ START
c ---- FA-01-31-11 ----- check vars ---- START 
c +++++++ test reading from rtdb occupations settings +++ START
c       if (ga_nodeid().eq.0)
c    &    write(*,*) 'In dft_scaleMO:: BEF rtdb_get'
        status = rtdb_get(rtdb,'focc:occupations',
     &                    mt_int,2,nfocc_in)
        nfocc_tot=nfocc_in(1)+nfocc_in(2)
        if (status) then
c        if (ga_nodeid().eq.0) then  
c         write(*,7) nfocc_in(1),nfocc_in(2),nfocc_tot
c7        format('rtdb-check: nfocc=(',i4,',',i4,',',i4,')')
c        endif
        if (.not. ma_push_get(mt_dbl,nfocc_tot,'occup1',
     &                        l_Iocc,k_Iocc))
     &  call errquit('occ_input ma_push_get fail k_Iocc',
     &               0,MA_ERR)
        endif   
        status = rtdb_get(rtdb,'focc:occup list',
     &                    mt_dbl,nfocc_tot,
     &                    dbl_mb(k_Iocc))    
       if (status) then ! --- START-fracc-occ-routine
         do i=1,2
          if(.not.ga_create(mt_dbl,nbf_mo,nbf_mo,
     &                      'Movecs Re',0,0, 
     &                      g_moso1(i)))     
     &    call errquit('scaleMO_so: error creating Movecs Re',0,
     &                 GA_ERR)
          call ga_copy(g_moso(i),g_moso1(i))
         enddo
c       if (ga_nodeid().eq.0)
c    &    write(*,*) 'In dft_scaleMO:: AFT rtdb_get'
c ------- check nel=sum(n_i) + int_mb(k_nIocc+i-1) =[1,nbf_mo]---- START
      ac_occ=0.0d0
      do i=1,nfocc_tot
c      if (ga_nodeid().eq.0) then
c       write(*,118) i,dbl_mb(k_Iocc+i-1)
c118    format('focc(',i3,')=',f15.8)
c      endif
       ac_occ=ac_occ+dbl_mb(k_Iocc+i-1)
      enddo
      diff_occ=abs(nTotOcc*1.0d0-ac_occ)
c     if (ga_nodeid().eq.0) then
c       write(*,18) nfocc_tot,nTotOcc,ac_occ,diff_occ
c 18    format('(nfocc_tot,nTotOcc,ac_occ,diff_occ)=(',
c    &       i3,',',i3,',',f15.8,',',f15.8,')')
c     endif
      error_Iocc=1e-2
      if (diff_occ .gt. error_Iocc) then ! TRUE IF-commented for the moment
       write(*,*) 'Error in dft_scf_so:',
     &            ' occupations keyword problem: ',
     &            'ac_occ from occupations keyword',
     &            ' is not close to nTotOcc'
       write(*,191) nTotOcc,ac_occ,diff_occ,error_Iocc
  191  format('(nTotOcc,ac_occ,diff_occ,error_Iocc)=(',i4,',',f15.8,
     &        ',',f15.8,',',f15.8)
       stop
      else
c ---- Calculate Sqrt(n_i)
c      if (ga_nodeid().eq.0)
c    &    write(*,*) 'Calculating sqrt(n_i):'
       do i=1,nfocc_tot
c       if (ga_nodeid().eq.0) then
c        write(*,24) i,dbl_mb(k_Iocc+i-1),sqrt(dbl_mb(k_Iocc+i-1))
c24      format('(n_i,sqrt(n_i))(',i3,')=(',f15.8,',',f15.8,')')
c       endif
        dbl_mb(k_Iocc+i-1)=sqrt(dbl_mb(k_Iocc+i-1))
       enddo
      endif
c ------- check nel=sum(n_i) + int_mb(k_nIocc+i-1) =[1,nbf_mo]---- END
c      write(*,*) 'BEF-scaling check g_moso ---------- START'
c      call ga_print(g_moso(1))
c      call ga_print(g_moso(2))
c      write(*,*) 'BEF-scaling check g_moso ---------- END'
c      if (ga_nodeid().eq.0) then
c        write(*,*) 'BEF. updating occupations:'
         do i=1,nbf_mo
c         write(*,1) i,occ(i)
c 1       format('occ(',i3,')=',f15.8)
          occ(i)=0.0d0 ! reset occupations
         enddo
c      endif
        ind_min=2
        ind_max=1
        if (nfocc_in(1).lt.nfocc_in(2)) then
         ind_min=1 
         ind_max=2
        endif
      do ispin=1,2
       iorb=ispin
       do i=1,nfocc_in(ind_min)
c       if (ga_nodeid().eq.0) then
c        write(*,111) i,iorb,dbl_mb(k_Iocc+iorb-1)
c111     format('CHECK-RTDB-occ:occ(',i3,',',i3,')=',f15.8)
c       endif
        call ga_scale_patch(g_moso1(1),1,nbf_mo,iorb,iorb,
     &                      dbl_mb(k_Iocc+iorb-1))
        call ga_scale_patch(g_moso1(2),1,nbf_mo,iorb,iorb,
     &                      dbl_mb(k_Iocc+iorb-1))
        occ(iorb)=dbl_mb(k_Iocc+iorb-1)*dbl_mb(k_Iocc+iorb-1) ! update occ
        iorb=iorb+2
       enddo ! end-loop-i
      end do   ! ispin-loop   
c ---- Storing the unpaired electrons ---- START
        iorb=2*nfocc_in(ind_min)+1 
        iorb1=iorb
        if (nfocc_in(1).lt.nfocc_in(2)) then
         iorb1=iorb+1
        endif
        do i=1,nfocc_in(ind_max)-nfocc_in(ind_min)
c        if (ga_nodeid().eq.0) then
c         write(*,151) i,iorb,iorb1,dbl_mb(k_Iocc+iorb-1)
c151      format('CHECK-RTDB-occ:occ-1(',
c    &           i3,',',i3,',',i3,')=',f15.8)
c        endif
         call ga_scale_patch(g_moso1(1),1,nbf_mo,iorb1,iorb1,
     &                       dbl_mb(k_Iocc+iorb-1))
         call ga_scale_patch(g_moso1(2),1,nbf_mo,iorb1,iorb1,
     &                       dbl_mb(k_Iocc+iorb-1))
         occ(iorb1)=dbl_mb(k_Iocc+iorb-1)*dbl_mb(k_Iocc+iorb-1) ! update occ -unpaired electron
         iorb=iorb+1
         iorb1=iorb1+2
        enddo ! end-loop-i
c ---- Storing the unpaired electrons ---- END
c ---- Calculate spin-orbit density matrix
       nbf_ao=nbf_mo/2
c      if (ga_nodeid().eq.0)
c    &  write(*,*) 'nbf_ao=',nbf_ao
       call dft_densm_so(g_densso,g_moso1,nbf_ao,nTotOcc)    
c----- Update occupations ----------- START
c      if (ga_nodeid().eq.0) then
c        write(*,*) 'AFT. updating occupations:'
c        do i=1,nbf_mo
c         write(*,2) i,occ(i)
c 2       format('occ(',i3,')=',f15.8)
c        enddo
c      endif
c----- Update occupations ----------- END
c ----- Free memory
          if (.not.ma_pop_stack(l_Iocc)) call
     &    errquit('dft_scaleMO_so: ma_pop_stack l_Iocc',0, MA_ERR)
        do i=1,2
          if (.not. ga_destroy(g_moso1(i))) 
     &       call errquit
     &       ('scaleMO_so: could not destroy g_moso_tmp',
     &        0, GA_ERR)
        enddo
        switch_sclMO_so=1 ! set ON
       endif ! --- END-fracc-occ-routine
       return
       end
c $Id$
