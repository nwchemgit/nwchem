       subroutine dft_scaleMO(rtdb,           ! IN
     &                        g_vecs,         ! IN     : MO vectors
     &                        occ,            ! IN/OUt : occupancies
     &                        g_dens,         ! OUT    : density matrix
     &                        noc,            ! IN/OUT : nr occupations
     &                        nbf,            ! IN     : nr. basis functions
     &                        npol,           ! IN     : nr of polarizations
     &                        nTotOcc)        ! IN     : nr total of occupations
c      Purpose: Use fractional occupations input using
c               occupations keyword for dft module 
c               entered from input script:
c               occupations nA nB
c               occA-1 occB-1
c               ...
c               occA-X  occB-nB
c               occA-nA
c      Author : Fredy Aquino
c      Date   : 02-22-11 (noc_A > noc_B)           
c      Date   : 05-11-11 arbirtrary (noc_A,noc_B) 
c      LIMITATIONS : This routine ONLY works for unrestricted calc.
       implicit none
#include "errquit.fh"
       integer rtdb
#include "mafdecls.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
       double precision occ(nbf*npol) ! occupancies
       integer g_dens(2)       ! spin-orbit density matrix
       integer scale_switch    ! input
       integer nbf             ! input
       integer npol            ! nr of polarizations
       integer g_vecs(npol)    ! input/output
       integer g_vecs1(npol)   ! scratch ga array
       integer noc(2)
       integer i,ispin,iorb,iorb1,nTotOcc,nbf_ao,
     &         ind_min,ind_max
       logical status
       integer l_Iocc,k_Iocc,nfocc_in(2), 
     &         nfocc_tot                
       double precision ac_occ,             
     &                  diff_occ,error_Iocc 

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
         do i=1,npol
          if(.not.ga_create(mt_dbl,nbf,nbf,
     &                      'Movecs1',0,0, 
     &                      g_vecs1(i)))     
     &    call errquit('scaleMO: error creating Movecs1',0,
     &                 GA_ERR)
          call ga_copy(g_vecs(i),g_vecs1(i))
         enddo
c       if (ga_nodeid().eq.0)
c    &    write(*,*) 'In dft_scaleMO:: AFT rtdb_get'
c ------- check nel=sum(n_i) + int_mb(k_nIocc+i-1) =[1,nbf]---- START
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
       write(*,*) 'Error in dft_scaleMO:',
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
c ------- check nel=sum(n_i) + int_mb(k_nIocc+i-1) =[1,nbf]---- END
c      if (ga_nodeid().eq.0) then
c        write(*,*) 'BEF. updating occupations:'
c        write(*,4) noc(1),noc(2)
c 4      format('noc=(',i5,',',i5,')') 
         do i=1,nbf*npol
c         write(*,1) i,occ(i)
c 1       format('occ(',i3,')=',f15.8)
          occ(i)=0.0d0 ! reset occupations
         enddo
c      endif
c      do ispin=1,npol
c        if (ga_nodeid().eq.0)
c     &     write(*,*) 'BEF-scaling check g_vecs1(',
c     &                ispin,') ---------- START'
c        call ga_print(g_vecs1(ispin))
c        if (ga_nodeid().eq.0)
c     &     write(*,*) 'BEF-scaling check g_vecs1(',
c     &                ispin,') ---------- END'
c      enddo ! end-loop-ispin
        ind_min=2
        ind_max=1
        if (nfocc_in(1).lt.nfocc_in(2)) then
         ind_min=1 
         ind_max=2
        endif
      do ispin=1,npol
       iorb1=nbf*(ispin-1)+1
       iorb=ispin
       do i=1,nfocc_in(ind_min)
c       if (ga_nodeid().eq.0) then
c        write(*,111) i,iorb,iorb1,dbl_mb(k_Iocc+iorb-1)
c111     format('CHECK-RTDB-occ:occ(',
c    &           i3,',',i3,',',i3,')=',f15.8)
c       endif
        call ga_scale_patch(g_vecs1(ispin),1,nbf,iorb1,iorb1,
     &                      dbl_mb(k_Iocc+iorb-1))
        occ(iorb1)=dbl_mb(k_Iocc+iorb-1)*dbl_mb(k_Iocc+iorb-1) ! update occ
        iorb=iorb+2
        iorb1=iorb1+1
       enddo ! end-loop-i
      end do   ! ispin-loop   
c ---- Storing the unpaired electrons ---- START
        iorb =2*nfocc_in(ind_min)+1              
        iorb1=nfocc_in(ind_min)+1
        if (nfocc_in(1).lt.nfocc_in(2)) then
         iorb1=nbf+nfocc_in(ind_min)+1
        endif
        do i=1,nfocc_in(ind_max)-nfocc_in(ind_min)
c        if (ga_nodeid().eq.0) then
c         write(*,151) i,iorb,iorb1,dbl_mb(k_Iocc+iorb-1)
c151      format('CHECK-RTDB-occ:occ-1(',
c    &           i3,',',i3,',',i3,')=',f15.8)
c        endif
         call ga_scale_patch(g_vecs1(ind_max),1,nbf,iorb1,iorb1,
     &                       dbl_mb(k_Iocc+iorb-1))
         occ(iorb1)=dbl_mb(k_Iocc+iorb-1)*dbl_mb(k_Iocc+iorb-1) ! update occ -unpaired electron
         iorb=iorb+1
         iorb1=iorb1+1
        enddo ! end-loop-i
c ---- Storing the unpaired electrons ---- END
c      do ispin=1,npol
c        if (ga_nodeid().eq.0)
c     &     write(*,*) 'AFT-scaling check g_vecs1(',
c     &                ispin,') ---------- START'
c        call ga_print(g_vecs1(ispin))
c        if (ga_nodeid().eq.0)
c     &     write(*,*) 'AFT-scaling check g_vecs1(',
c     &                ispin,') ---------- END'
c      enddo ! end-loop-ispin
c ---- Calculate density matrix ------------ START
       do ispin=1,npol
         call ga_sync()
         call ga_dgemm('N','T',nbf,nbf,nfocc_in(ispin),
     &                 1.0d0,g_vecs1(ispin),g_vecs(ispin), 
     &                 0.0d0,g_dens(ispin))
         noc(ispin)=nfocc_in(ispin)  
       enddo ! end-loop-ispin
c ---- Calculate density matrix ------------ END
c----- Update occupations ----------- START
c      if (ga_nodeid().eq.0) then
c        write(*,*) 'AFT. updating occupations:'
c        write(*,3) noc(1),noc(2)
c 3      format('noc=(',i5,',',i5,')') 
c        do i=1,nbf*npol
c         write(*,2) i,occ(i)
c 2       format('occ(',i3,')=',f15.8)
c        enddo
c      endif
c----- Update occupations ----------- END
c ----- Free memory
          if (.not.ma_pop_stack(l_Iocc)) call
     &    errquit('dft_scaleMO_so: ma_pop_stack l_Iocc',0, MA_ERR)
        do i=1,npol
          if (.not. ga_destroy(g_vecs1(i))) 
     &       call errquit
     &       ('scaleMO: could not destroy g_vecs1',
     &        0, GA_ERR)
        enddo
       endif ! --- END-fracc-occ-routine
       return
       end
