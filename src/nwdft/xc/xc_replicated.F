      subroutine xc_rep_init(rtdb, geom, basis, nbf_ao,ipol,goforit)
c
C$Id: xc_replicated.F,v 1.2 2003-06-05 22:06:21 edo Exp $
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
      integer geom,basis
c
      integer xcrep_glob
      logical xcreplicated_l,oprint
      integer plist,dims(2),map(2),nblock(2),l_map,k_map
c
      oprint = util_print('xcreplicated', print_high)
c
c     go replicated for xc ?
c
      if(.not.goforit) then
         xcreplicated=.false.
         return
      endif
      if (.not. rtdb_get(rtdb, 'dft:xcrepoldway', mt_log, 1, 
     .     oldway))
     &     oldway=.true.
      if (.not. rtdb_get(rtdb, 'dft:xcreplicated', mt_log, 1, 
     .     xcreplicated_l))
     &     xcreplicated_l=.false.
      if(xcreplicated_l) then
         if(oldway) then
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
         else
c
c     use of ga_mirror
c
            plist=ga_mirror_config()
            dims(1)=nbf_ao
            dims(2)=nbf_ao
            if(.not.ma_push_get(MT_int,ga_cluster_nnodes()*2,
     .           'ga block maps',l_map,k_map))
     .           call errquit('gacratblock: not enough local mem',0)
            call util_atom_blocked(ga_nnodes(),
     .           geom, basis, int_mb(k_map), nblock)

            xcreplicated_l=nga_create_irreg_config(MT_DBL, 2, dims,
     .           'replicated xc GA',
     .           int_mb(k_map),nblock,
     .           plist, k_repxc(1))
            if(ipol.eq.2) then
               xcreplicated_l=nga_create_irreg_config(MT_DBL, 2, dims,
     .              'replicated xc GA',
     .              map,nblock,
     .              plist, k_repxc(2))
            endif
         if (.not.ma_chop_stack(l_map))
     &        call errquit('xc_rep_init: cannot chop stack',0)

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
         if(oldway) then
         if (.not.ma_chop_stack(l_repxc(1)))
     &        call errquit('xc_rep_init: cannot chop stack',1)
         else
cedo         if (.not.ga_destroy(l_repxc(1)))
cedo     &        call errquit('xc_rep_init: cannot chop stack',1)
         endif
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
C$Id: xc_replicated.F,v 1.2 2003-06-05 22:06:21 edo Exp $
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
         if(oldway) then
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
         else
            do i = 1, ipol
               call ga_merge_mirrored(k_repxc(i))
               call ga_copy(k_repxc(i),g_vxc(i))
               if (.not.ga_destroy(k_repxc(i)))
     &              call errquit('xc_rep_close: cannot gadestr',i)
            enddo
         endif
      endif
      return
      end

