      subroutine xc_rep_init(rtdb, nbf_ao,ipol,goforit)
c
C$Id: xc_replicated.F,v 1.1 2002-05-15 16:43:03 edo Exp $
c
      implicit none
c      
      integer rtdb
c
#include "mafdecls.fh"
#include "rtdb.fh"
#include "dftpara.fh"
#include "global.fh"
#include "util.fh"
#include "stdio.fh"
      integer nbf_ao,ipol
      logical goforit
c
      integer xcrep_glob
      logical xcreplicated_l,oprint
c
      oprint = util_print('xcreplicated', print_high)
c
c     go replicated for xc ?
c
      if(.not.goforit) then
         xcreplicated=.false.
         return
      endif
      if (.not. rtdb_get(rtdb, 'dft:xcreplicated', mt_log, 1, 
     .     xcreplicated_l))
     &     xcreplicated_l=.false.
      if(xcreplicated_l) then
         nbf_ld=nbf_ao
c     
c     allocate replicated Fock matrix
c
         xcreplicated_l=ma_push_get(MT_Dbl,nbf_ao*nbf_ao*ipol,
     .        'xcmat',l_repxc(1),k_repxc(1))
         if(xcreplicated_l) then
            call dcopy(nbf_ao*nbf_ao*ipol,0d0,0,dbl_mb(k_repxc(1)),1)
            k_repxc(2)=k_repxc(1)+nbf_ao*nbf_ao
         endif
         
      endif
c
c     check if everyone is ok
c
      xcrep_glob=0
      if(xcreplicated_l) xcrep_glob=1
      call ga_igop(375,xcrep_glob,1, '+')
      xcreplicated=xcrep_glob.eq.ga_nnodes()
c
c     deallocate MA if cannot go replicated globally
c
      if(.not.xcreplicated.and.xcreplicated_l) then
         if (.not.ma_chop_stack(l_repxc(1)))
     &        call errquit('xc_rep_init: cannot chop stack',1)
      endif
      if(oprint.and.ga_nodeid().eq.0) then
         if(xcreplicated) then
            write(luout,*) ' going xcreplicated '
         else
            write(luout,*) ' not enough memory for xcreplicated '
         endif
      endif
      return
      end

      logical function xc_rep_close(rtdb, nbf_ao,ipol,g_vxc)
c
C$Id: xc_replicated.F,v 1.1 2002-05-15 16:43:03 edo Exp $
c
      implicit none
c      
      integer rtdb
c
#include "mafdecls.fh"
#include "global.fh"
#include "dftpara.fh"
      integer g_vxc(2),ipol,nbf_ao
c
      integer i
c
      xc_rep_close=.true.
      if(xcreplicated) then
         call ga_dgop(374,dbl_mb(k_repxc(1)),nbf_ao*nbf_ao*ipol, '+')
         if (ga_nodeid() .eq. 0) then
            do i = 1, ipol
               call ga_put(g_vxc(i), 1, nbf_ao, 1, nbf_ao, 
     $              dbl_mb(k_repxc(i)), nbf_ao)
            end do
         end if
!         if (.not.ma_chop_stack(l_repxc(1)))
!     &        call errquit('xc_getv: cannot chop stack',1)
         xc_rep_close=ma_chop_stack(l_repxc(1))
      endif
      return
      end

