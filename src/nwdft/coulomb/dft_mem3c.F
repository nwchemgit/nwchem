      logical function dft_mem3c(rtdb,
     I     natoms,npol,oprint_parm,oprint_3c2e,
     O     n3c_int,n3c_dbl,
     O     l_3ceri,k_3ceri, l_3cwhat,k_3cwhat)
* $Id$
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "msgids.fh"
#include "cdft.fh"
#include "basP.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "util.fh"
      integer rtdb   ! [in]
      integer natoms ! [in]
      integer npol ! [in]
      integer n3c_int ! [out]
      integer n3c_dbl ! [out]
      logical oprint_parm,oprint_3c2e ![in]
      integer l_3ceri,k_3ceri ! [out]
      integer l_3cwhat,k_3cwhat ! [out]
      integer max_component,max_elem_ang_scr
      integer nscr
      integer avail,availm
      double precision availm_r
      integer dft_n3cint
      external dft_n3cint
      double precision dft_n3cdbl
      integer dft_nao2_max
      external dft_n3cdbl,dft_nao2_max
      double precision n3c_dbl_r,n3c_dbl_max
      integer nao2_max
      integer nshpair_sum
      integer me,nproc,icount
      integer maxg,scrmx
      integer ishc,itype,nprimo,nshbfc,isphere
      integer deficit,atom_c,sh_lo_c,sh_hi_c
      logical spherical_ao
      integer n3c_int_in,n3c_dbl_in
c
      me=ga_nodeid()
      nproc=ga_nnodes()
      dft_mem3c=.true.
c     
c     Determine how big a buffer can be allocated to 3-center
c     2e- integrals.
c     
c     - amount needed for all incore:
c     
      n3c_dbl_r = dft_n3cdbl()
      n3c_int = dft_n3cint()
c     either read no. of batches or make it at least 1000
      if(.not.rtdb_get(rtdb,'dft:n3c_int',mt_int,1,n3c_int_in))
     N     n3c_int_in=1000
      n3c_int=max(n3c_int,n3c_int_in)
c     either read no. of batches or make it at least 1000
      if(.not.rtdb_get(rtdb,'dft:n3c_dbl',mt_int,1,n3c_dbl_in))
     N     n3c_dbl_in=1000
      n3c_dbl_r=max(n3c_dbl_r,1d0*n3c_dbl_in)
      n3c_dbl_max=n3c_dbl_r
c     need to find max among all procs since n3cdbl might be different
      call ga_dgop(19641964, n3c_dbl_max, 1, 'max')
c     
c     find - (minimum)amount local available memory on all nodes
c     
      avail = MA_inquire_avail(mt_dbl)
      call ga_igop(msg_min_stack_avail, avail, 1, 'min')
c     
c     estimate and subtract off amount needed for DIIS
c
      availm = avail - ((nfock+4)*nbf_ao*nbf_ao)/nproc
c     
c        estimate and subtract off amount needed for XC numerical integration
c     in xc_quadv0
c
      availm = availm - (natoms*(natoms+1)/2 + 13*natoms + 
     &     3*nqmax*(7*ipol +  npol + natoms + 
     &     nbf_ao_mxnbf_ce + 4) + 
     &     nbf_ao_mxcont + nbf_ao_mxprim + 
     &     2*nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce)
c     
c     estimate and subtract off amount needed for XC numerical integration
c     in xc_quadv0_a
c     
      max_component = (nbf_ao_mxang+1)*(nbf_ao_mxang+2)/2 
      if (nbf_ao_mxang .eq. 0)then
         max_elem_ang_scr = max_component * 3
      elseif (nbf_ao_mxang .le. 3)then
         max_elem_ang_scr = max_component * 9
      else                      ! general case
         max_elem_ang_scr = max_component * 28
      endif 
c     
      nscr = 3*nqmax*nbf_ao_mxcont + 
     &     max(3*nqmax*nbf_ao_mxprim,nqmax*max_elem_ang_scr) + 1
c     
c     The big chunk is the memory needed for new_eval_gbsets 
c     which is roughly 4*nqmax*nbf_ao.  This is reduced by
c     screening (and chunking up the angular grid) and is 
c     computed at the end of xc_setquad to be 4*max_pr_mbfnq.
c     
      availm = availm - (nqmax*(natoms + 3*nbf_ao_mxnbf_ce + 1) +
     &     4*max_pr_mbfnq + 
     &     nbf_ao + nscr)
c     
c        Subtract off a few extra bits
c     
      availm = availm - 100000
c     
      if(availm.lt.0)then
         availm = 0
      endif
      availm_r = dble(availm)
      availm_r = min(availm_r,n3c_dbl_max)
      deficit=0
      if (availm_r.lt.n3c_dbl_max)then
         deficit=nint(n3c_dbl_max)-availm_r
c
c     get amount of local MA in Mbytes need to get incore done
c
         deficit=max((deficit*8+deficit/2)/1024/1024,1)
c     
c     force direct if  memory is not sufficient
c     
         availm_r = 0
      endif
c
c     check if availm_r is big enough for max batch
c
      call int_mem_2e3c(maxg, scrmx)

      n3c_dbl = min(availm_r,n3c_dbl_r)
      if (me.eq.0 .and. oprint_parm)
     &     write(LuOut,3228)avail, availm, 
     N     n3c_dbl,nint(n3c_dbl_max),n3c_int
      
      if(deficit.ne.0.and.me.eq.0) write(luout,3230) deficit
c     
c     Loops are parallelized over the products of atoms
c     (check for zero ... must be at least 1).
c     

      if (direct.or.(deficit.ne.0))then
         dft_mem3c = .false.
         n3c_dbl = 1
         n3c_int = 1
      endif
      if (.not.MA_Push_Get
     &     (MT_Dbl,n3c_dbl,'3c ERI space',l_3cERI,k_3cERI))
     &     call errquit('dft_scf: push_get failed', 12, MA_ERR)
c     
      if (.not.MA_Push_Get
     &     (MT_int,n3c_int,'3c what space',l_3cwhat,k_3cwhat))
     &     call errquit('dft_scf:push_get failed', 13, MA_ERR)
c     
      if (dft_mem3c)then
         if (me.eq.0 .and. oprint_3c2e)
     &        write(LuOut,3229)n3c_dbl*1.d-6
      endif
c     if(me.eq.0) call MA_summarize_allocated_blocks()
         return
 3228 format(10x,'Minimum dble words available (all nodes) is: ',i15,
     &     /,10x,'         This is reduced (for later use) to: ',i15,
     &     /,10x,'            proc 0 Suggested buffer size is: ',i15,
     &     /,10x,'               Max Suggested buffer size is: ',i15,
     &     /,10x,'                    no. integral batches is: ',i15)

 3230 format(/7x,
     &     'Extra Local Memory (stack+heap) needed for incore:',i6,
     &     ' Mbytes ')
 3229 format(/,6x,f10.3,' MW buffer allocated for incore 3-center '
     &     /,10x,'2e- integral storage on stack. ')
         end
c     hack for spherical since there is always a cartesian internal bit
      integer function dft_nao2_max(basis,natoms)
      implicit none
#include "bas.fh"
#include "errquit.fh"
      integer basis
      integer natoms
c
      integer atom_c,atom_d
      integer ishc,itype,nprimo,nshbfc,isphere
      integer sh_lo_c, sh_hi_c
      integer nshpair_sum
c
      dft_nao2_max=0
      do atom_c=1,natoms
         if (.not. bas_ce2cnr( basis, atom_c, sh_lo_c, sh_hi_c)) 
     &        call errquit('Exiting in dft_fitcd',110, BASIS_ERR)
         nshpair_sum=0
         do ishc=sh_lo_c,sh_hi_c
            if(.not.bas_continfo(basis,ishc,
     &           itype,nprimo,nshbfc,isphere))
     &           call errquit('Exiting in fitcd.',44, CALC_ERR)
            
            nshpair_sum = nshpair_sum +
     P           ((itype+1)*(itype+2))/2*nshbfc
         enddo
         dft_nao2_max=max(dft_nao2_max,nshpair_sum)
      enddo
      dft_nao2_max=dft_nao2_max**2
      return
      end
