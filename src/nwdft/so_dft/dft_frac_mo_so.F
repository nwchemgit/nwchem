c
c      Fractional MO occupations handler
c      occup nA nB
c          occA-1 occB-1
c          ...
c          occA-X  occB-nB
c          occA-nA
c      end
c
c
       subroutine dft_frac_mo_so(rtdb,           ! IN
     &                           g_moso,         ! IN     : MO vectors
     &                           focc,           ! IN/OUT : occupancies
     &                           g_densso,       ! OUT    : spin-orbit density matrix
     &                           nbf,            ! IN     : nr. basis functions
     &                           ntot_occ)       ! IN     : nr. occ
c
       implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
c
       integer rtdb
       integer nbf          ! input
       double precision focc(nbf) ! occupancies
       integer g_densso(2)     ! spin-orbit density matrix
       integer g_moso(2)       ! input/output
       integer g_moso_tmp(2) ! scratch ga array
       integer i,ispin,iorb,iorb1,ntot_occ,nbf_ao
       logical status
       integer l_occ,k_occ,nfocc_in(2), 
     &         nfocc,ind_min,ind_max               
       double precision user_occ,diff_occ,error_occ 
       logical loccup,loccup_list
       external dft_densm_so
c
c      get the user defined occupations
       nfocc = 0
       do i =1,2
         nfocc_in(i) = 0
       end do
c
       loccup = rtdb_get(rtdb,'focc:occup',mt_int,2,nfocc_in)
       if (.not. loccup) return
       nfocc=nfocc_in(1)+nfocc_in(2)
       if (nfocc.eq.0) return
c
c      memory for the user occupations
       if (.not. ma_push_get(mt_dbl,nfocc,'occups',l_occ,k_occ))
     &  call errquit('occ_input ma_push_get fail k_occ',0,MA_ERR)
c
c      get user occupation list
       loccup_list = rtdb_get(rtdb,'focc:occup_list',mt_dbl,nfocc,
     &   dbl_mb(k_occ))
       if (.not.loccup_list) return
c 
       do i=1,2
        if(.not.ga_create(mt_dbl,nbf,nbf,'Movecs Re',0,0,g_moso_tmp(i)))
     &  call errquit('scaleMO_so: error creating Movecs Re',0,GA_ERR)
        call ga_copy(g_moso(i),g_moso_tmp(i))
       enddo
c
c      tally up user occupation
       user_occ=0.0d0
       do i=1,nfocc
         user_occ=user_occ+dbl_mb(k_occ+i-1)
       enddo
       diff_occ=abs(ntot_occ*1.0d0-user_occ)
c
       do i=1,nfocc
        dbl_mb(k_occ+i-1)=sqrt(dbl_mb(k_occ+i-1))
       enddo
c
       do i=1,nbf
        focc(i)=0.0d0 ! reset occupations
       enddo
c
       ind_min=2
       ind_max=1
       if (nfocc_in(1).lt.nfocc_in(2)) then
         ind_min=1 
         ind_max=2
       endif
c
       do ispin=1,2
        iorb=ispin
        do i=1,nfocc_in(ind_min)
         call ga_scale_patch(g_moso_tmp(1),1,nbf,iorb,iorb,
     &                      dbl_mb(k_occ+iorb-1))
         call ga_scale_patch(g_moso_tmp(2),1,nbf,iorb,iorb,
     &                      dbl_mb(k_occ+iorb-1))
         focc(iorb)=dbl_mb(k_occ+iorb-1)*dbl_mb(k_occ+iorb-1) ! update occ
         iorb=iorb+2
        enddo ! end-loop-i
      end do   ! ispin-loop   
c
c     unpaired electrons
      iorb=2*nfocc_in(ind_min)+1 
      iorb1=iorb
      if (nfocc_in(1).lt.nfocc_in(2)) iorb1=iorb+1
c
      do i=1,nfocc_in(ind_max)-nfocc_in(ind_min)
         call ga_scale_patch(g_moso_tmp(1),1,nbf,iorb1,iorb1,
     &                       dbl_mb(k_occ+iorb-1))
         call ga_scale_patch(g_moso_tmp(2),1,nbf,iorb1,iorb1,
     &                       dbl_mb(k_occ+iorb-1))
         focc(iorb1)=dbl_mb(k_occ+iorb-1)*dbl_mb(k_occ+iorb-1) ! update occ -unpaired electron
         iorb=iorb+1
         iorb1=iorb1+2
      enddo ! end-loop-i
c
c     calculate spin-orbit density matrix
      nbf_ao=nbf/2
      call dft_densm_so(g_densso,g_moso_tmp,nbf_ao,ntot_occ)    
c
c     clean up
      if (.not.ma_pop_stack(l_occ)) 
     &  call errquit('dft_scaleMO_so: ma_pop_stack l_occ',0, MA_ERR)
       do i=1,2
        if (.not. ga_destroy(g_moso_tmp(i))) 
     &   call errquit('scaleMO_so: could not destroy g_moso_tmp',
     &   0,GA_ERR)
      enddo
c
      return
      end
c $Id$
