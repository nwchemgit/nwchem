      logical function dft_mem3c(
     I     natoms,npol,oprint_parm,oprint_3c2e,
     O     n3c_int,n3c_dbl,n_semi_bufs,
     O     l_3ceri,k_3ceri, l_3cwhat,k_3cwhat)
* $Id: dft_mem3c.F,v 1.8 2004-04-28 05:09:36 edo Exp $
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "msgids.fh"
#include "cdft.fh"
#include "basP.fh"
#include "bas.fh"
      integer natoms ! [in]
      integer npol ! [in]
      integer n3c_int ! [out]
      integer n3c_dbl ! [out]
      integer  n_semi_bufs ! [out]
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
      external dft_n3cdbl
      double precision n3c_dbl_r
      integer nao2_max,batch_max
      integer me,nproc,icount
      integer maxg,scrmx
      integer ishc,itype,nprimo,nshbfc,isphere
      integer deficit
      logical spherical_ao
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
      availm_r = min(availm_r,n3c_dbl_r)
      deficit=0
      if (availm_r.lt.n3c_dbl_r)then
         deficit=n3c_dbl_r-availm_r
c
c     get amount of local MA in Mbytes need to get incore done
c
         deficit=max((deficit*8+deficit/2)/1024/1024,1)
c     
c     cannot all be done incore so let's restrict memory use to a
c     max of 5 MW for semi-direct (I/O)
c     
         availm_r = min(availm_r,5d+6)
      endif
c
c     check if availm_r is big enough for max batch
c
      call int_mem_2e3c(maxg, scrmx)

      n3c_dbl = nint(availm_r)
      if (me.eq.0 .and. oprint_parm)
     &     write(LuOut,3228)avail, availm, n3c_dbl
      
      if(deficit.ne.0.and.me.eq.0) write(luout,3230) deficit
c     
c     Loops are parallelized over the products of AO shells
c     (check for zero ... must be at least 1).
c     
      icount = (nshells_ao + mod(nshells_ao,2))/(2*nproc)
      icount = max(icount,1)
      if (nproc.eq.1)icount = 1
c     
      nao2_max = nbf_ao_mxnbf_cn**2
      spherical_ao=bas_spherical(ao_bas_han)
      if(spherical_ao) then
         nao2_max=0
         do ishc=1,nshells_ao
            if(.not.bas_continfo(ao_bas_han,ishc,
     &           itype,nprimo,nshbfc,isphere))
     &           call errquit('Exiting in fitcd.',44, CALC_ERR)
           nao2_max=max(nao2_max,
     ,           ((itype+1)*(itype+2))/2*nshbfc)
         enddo
         nao2_max=nao2_max**2
      endif
c     
c     compute size of largest block of integrals computed
c     
      batch_max = nao2_max*nbf_cd*icount
c     
      if ((n3c_dbl.lt.batch_max).or.direct)then
         dft_mem3c = .false.
         n3c_dbl = 1
         n3c_int = 1
      endif
      if (availm_r.lt.n3c_dbl_r)then
c     
c     cannot all be done incore, compute number of buffers 
c     required if semi-direct (I/O)
c     
         n_semi_bufs = nint(n3c_dbl_r/availm_r) 
c         n_semi_bufs = nint(n3c_dbl_r/availm_r) + 1
      else
         n_semi_bufs = 0
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
     &     /,10x,'                   Suggested buffer size is: ',i15)
 3230 format(/7x,
     &     'Extra Local Memory (stack+heap) needed for incore:',i6,
     &     ' Mbytes ')
 3229 format(/,6x,f10.3,' MW buffer allocated for incore 3-center '
     &     /,10x,'2e- integral storage on stack. ')
         end
