      subroutine xc_rep_init(rtdb, geom, basis, g_dens,
     ,     nbf_ao,ipol,goforit,dograd)
c
C$Id: xc_replicated.F,v 1.10 2003-06-19 19:57:38 edo Exp $
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
      logical dograd ! true if just gradient (not energy)
      integer geom,basis
      integer g_dens(*)
c
      integer xcrep_glob,ii
      logical xcreplicated_l,oprint
      integer plist,dims(2),nblock(2),l_map,k_map,procsmp
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
      if(.not.xcreplicated_l) then
         xcreplicated=.false.
         return
      endif
      if(xcreplicated_l) then
         xcreplicated_l=.false.
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
            procsmp=ga_cluster_nprocs(ga_cluster_nodeid())
            if(.not.ma_push_get(MT_int,procsmp*2,
     .           'ga block maps',l_map,k_map)) then
               call errquit('xcrepl: not enough local mem',
     .              procsmp*2)
               endif
            call util_atom_blocked(procsmp,
     .           geom, basis, int_mb(k_map), nblock)

c
c     compress map array
c
            do ii=0,nblock(2)-1
               int_mb(k_map+nblock(1)+ii)=
     =              int_mb(k_map+procsmp+ii)
            enddo
            do ii=1,ipol
#ifdef GAIRREG
               if(.not.dograd)
     .              if(.not.nga_create_irreg_config(MT_DBL, 2, dims,
     .              'replicated xc GA',
     .              int_mb(k_map),nblock,plist,k_repxc(ii))) goto 1984
               if(.not.nga_create_irreg_config(MT_DBL, 2, dims,
     .              'replicated DM GA',
     .              int_mb(k_map),nblock,plist,g_repdm(ii))) goto 1984
#else
               nblock(1)=0
               nblock(2)=0
               if(.not.dograd) then
                  if(.not.nga_create_config(MT_DBL, 2, dims,
     .                 'replicated xc GA',
     .                 nblock,plist, k_repxc(ii))) goto 1984
                  call ga_zero(k_repxc(ii))
               endif
               if(.not.nga_create_config(MT_DBL, 2, dims,
     .              'replicated DM GA',
     .              nblock,plist, g_repdm(ii))) goto 1984
               call ga_copy(g_dens(ii),g_repdm(ii))
#endif
#ifdef DEBUG
               call ga_print_distribution(k_repxc(ii))
#endif
            enddo
            if (.not.ma_pop_stack(l_map))
     &           call errquit('xc_rep_init: cannot pop stack',11)

         endif
         
      endif
      xcreplicated_l=.true.
c
c     check if everyone is ok
c
 1984 continue
      xcrep_glob=0
      if(xcreplicated_l) xcrep_glob=1
      call ga_igop(375,xcrep_glob,1, '+')
      xcreplicated=xcrep_glob.eq.ga_nnodes()
c
c     deallocate MA if cannot go replicated globally
c
      if(.not.xcreplicated.and.xcreplicated_l) then
         if(oldway) then
         if (.not.ma_pop_stack(l_repxc(1)))
     &        call errquit('xc_rep_init: cannot pop stack',21)
         endif
      endif
      if(oprint.and.ga_nodeid().eq.0) then
         if(xcreplicated) then
            write(luout,*) ' going xcreplicated '
         else
            write(luout,*) ' not enough memory for xcreplicated '
         endif
      endif
#ifdef DEBUG
      call ga_summarize(1)
      call ma_summarize_allocated_blocks()
#endif
      return
      end

      logical function xc_rep_close(rtdb, nbf_ao,ipol,g_vxc,dograd)
c
C$Id: xc_replicated.F,v 1.10 2003-06-19 19:57:38 edo Exp $
c
      implicit none
c      
      integer rtdb
c
#include "mafdecls.fh"
#include "global.fh"
#include "dftpara.fh"
      integer g_vxc(2),ipol,nbf_ao
      logical dograd ! true if just gradient (not energy)
c
      integer i
c
      xc_rep_close=.true.
      if(xcreplicated) then
         if(oldway) then
         if(.not.dograd) then
         call ga_dgop(374,dbl_mb(k_repxc(1)),nbf_ao*nbf_ao*ipol, '+')
         if (ga_nodeid() .eq. 0) then
            do i = 1, ipol
               call ga_put(g_vxc(i), 1, nbf_ao, 1, nbf_ao, 
     $              dbl_mb(k_repxc(i)), nbf_ao)
            end do
         end if
         xc_rep_close=ma_chop_stack(l_repxc(1))
         endif
         else
            do i = 1, ipol
               if(.not.dograd) then
                  call ga_merge_mirrored(k_repxc(i))
                  call ga_copy(k_repxc(i),
     &                 g_vxc(i))
                  if (.not.ga_destroy(k_repxc(i)))
     &                 call errquit('xc_rep_close: cannot gadestr',i)
                  endif
               if (.not.ga_destroy(g_repdm(i)))
     &              call errquit('xc_rep_close: cannot gadestr',i)
            enddo
         endif
c
c     avoid getting things out of control
c
      xcreplicated=.false.
#ifdef DEBUG
      call ga_summarize(1)
      call ma_summarize_allocated_blocks()
#endif
      endif
      return
      end

