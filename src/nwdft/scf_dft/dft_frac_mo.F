c
c      Fractional MO occupations handler
c      occup nA nB
c          occA-1 occB-1
c          ...
c          occA-X  occB-nB
c          occA-nA
c      end
c
       subroutine dft_frac_mo(rtdb,           ! IN     : rtdb handle
     &                        g_vecs,         ! IN     : MO vectors
     &                        focc,           ! OUT    : occupancies
     &                        nbf,            ! IN     : nr. basis functions
     &                        npol,           ! IN     : nr of polarizations
     &                        ntot_occ)       ! IN     : nr total of occupations
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
       integer nbf             ! input
       integer npol            ! nr of polarizations
       double precision focc(nbf*npol) ! fractional occupancies
       integer g_dens(npol)       ! spin-orbit density matrix
       integer g_vecs(npol)    ! input/output
       integer i,ispin,iorb,iorb1,ntot_occ,nbf_ao,ind_min,ind_max
       logical loccup_list,loccup,lerror_check
       integer l_occ,k_occ,nfocc_in(npol),nfocc
       double precision user_occ,diff_occ,error_occ 
c
c      get the user defined occupations
       nfocc = 0
       do i =1,npol
         nfocc_in(i) = 0
       end do
c
       loccup = rtdb_get(rtdb,'focc:occup',mt_int,2,nfocc_in) ! read both nfocc_in(1) and nfocc_in(2)
       if (.not. loccup) return
       if (npol.eq.1) then
          nfocc=nfocc_in(1)
       else
          nfocc=nfocc_in(1)+nfocc_in(2)
       end if
       if (nfocc.eq.0) return
c
c      memory for the user occupations
       if (.not. ma_push_get(mt_dbl,nfocc,'occups',l_occ,k_occ))
     &  call errquit('dft_frac_mo: k_occ error',0,MA_ERR)
c
c      get user occupation list
       loccup_list = rtdb_get(rtdb,'focc:occup_list',mt_dbl,nfocc,
     &   dbl_mb(k_occ))
       if (.not.loccup_list) return
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
c      initialize fractional occupations
       do i=1,nbf*npol
         focc(i)=0.0d0
       enddo
c
c      determine majority and minority component
       ind_min=2
       ind_max=1
       if (nfocc_in(1).lt.nfocc_in(2)) then
         ind_min=1 
         ind_max=2
       endif
c
       do ispin=1,npol
        iorb1=nbf*(ispin-1)+1
        iorb=ispin
        do i=1,nfocc_in(ind_min)
          focc(iorb1)=dbl_mb(k_occ+iorb-1)*dbl_mb(k_occ+iorb-1)
          iorb=iorb+2
          iorb1=iorb1+1
        end do ! i loop
       end do   ! ispin loop   
c
       iorb =2*nfocc_in(ind_min)+1              
       iorb1=nfocc_in(ind_min)+1
       if (nfocc_in(1).lt.nfocc_in(2)) then
         iorb1=nbf+nfocc_in(ind_min)+1
       endif
c
c      unoccupied electron
       do i=1,nfocc_in(ind_max)-nfocc_in(ind_min)
         focc(iorb1)=dbl_mb(k_occ+iorb-1)*dbl_mb(k_occ+iorb-1)
         iorb=iorb+1
         iorb1=iorb1+1
       enddo ! i loop
c
c      clean up
       if (.not.ma_pop_stack(l_occ)) 
     &     call errquit('dft_frac_mo: ma_pop_stack l_occ',0, MA_ERR)
c
       return
       end
c $Id$
