C
C undefined the first define line if you get undefined symbols
C
C#define OLD_GAAPI 1
      subroutine xc_rep_init(rtdb, geom, basis, g_dens,g_vxc,
     ,     nbf_ao,ipol,goforit,dograd)
c
C$Id: xc_replicated.F,v 1.29 2006-04-25 01:41:21 edo Exp $
c
      implicit none
#include "errquit.fh"
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
      integer g_dens(*),g_vxc(*)
c
      integer xcrep_glob,ii
      logical xcreplicated_l,oprint,xcrepirrg,xcrep_ga
      logical util_mirrmat
      external util_mirrmat
      integer plist,dims(2),nblock(2)
      integer l_map,k_map,procsmp
      integer gnblock(2), membyte
      integer idum,g_ndim,idumm(10),dm_ndim,xc_ndim
c
      oprint = util_print('xcreplicated', print_high)

      k_repxc(1)=-1
c     nonblocking
c     non blocking
      if (.not. rtdb_get(rtdb, 'dft:nbacc', mt_log, 1,
     N  dftnbacc)) dftnbacc=.false.
      if (.not. rtdb_get(rtdb, 'dft:nbget', mt_log, 1,
     N  dftnbget)) dftnbget=.false.

c     we do not mirror 3D GAs for now

      call nga_inquire(g_dens(1), idum, dm_ndim, idumm)
      call nga_inquire(g_vxc(1), idum, xc_ndim, idumm)

c  Check number of nodes to decide to use Mirror Array
        if(ga_cluster_nnodes().eq.1) then
         xcreplicated=.false.
         return
        endif
        

      if (.not. rtdb_get(rtdb, 'dft:dorepdm', mt_log, 1, 
     .     dorepdm))
     &     dorepdm=.true.
      if (.not. rtdb_get(rtdb, 'dft:dorepxc', mt_log, 1, 
     .     dorepxc))
     &     dorepxc=.true.
      if(dorepxc) dorepxc=xc_ndim.eq.2
      if (.not. rtdb_get(rtdb, 'dft:xcreplicated', mt_log, 1, 
     .     xcreplicated_l))
     &     xcreplicated_l=.false.
      if (.not. rtdb_get(rtdb, 'dft:xcrepirrg', mt_log, 1, 
     .     xcrepirrg)) xcrepirrg=.false.

         xcreplicated_l=.false.
	 xcrep_ga=.false.

c     
c     use of ga_mirror
c
#ifdef OLD_GAAPI
         plist=ga_mirror_config()
#else
         plist=ga_pgroup_get_mirror()
#endif
         dims(1)=nbf_ao
         dims(2)=nbf_ao
         if(xcrepirrg) then
            procsmp=ga_cluster_nprocs(ga_cluster_nodeid())
            if(.not.ma_push_get(MT_int,procsmp*2,
     .           'ga block maps',l_map,k_map)) then
               call errquit('xcrepl: not enough local mem',
     .              procsmp*2, MA_ERR)
            endif

            call util_atom_blocked(procsmp,
     .        geom, basis, int_mb(k_map), nblock)
            
c     
c     compress map array
c     
            do ii=0,nblock(2)-1
               int_mb(k_map+nblock(1)+ii)=
     =              int_mb(k_map+procsmp+ii)
            enddo
         endif

         do ii=1,ipol
            gnblock(1)=0
            gnblock(2)=0

            if(dorepxc) then
               if(.not.dograd) then
                  if(xcrepirrg) then
	         xcrep_ga=nga_create_irreg_config(MT_DBL, 2, dims,
     .                    'replicatedxcGA',int_mb(k_map),nblock,plist,
     .                 k_repxc(ii))
	         if(.not.xcrep_ga) goto 1984
                 call ga_zero(k_repxc(ii))
            else
c check available memory before creat mirrored arroy for Fock matrix
#ifdef DEBUG
            membyte=ga_memory_avail()
            call ga_igop(1,membyte,1,'+')
	    if(ga_nodeid()==0)print*,membyte/1000000, 'MB ENTER K_REP'
#endif
                 xcrep_ga=util_mirrmat(1,g_dens(ii),k_repxc(ii),
     .                                  .false.,.true.)

               if(.not.xcrep_ga) goto 1984
c check available memory after array being created(or not created)
#ifdef DEBUG
            if(ga_nodeid()==0) print*,'xcrep_ga',xcrep_ga
            membyte=ga_memory_avail()
            call ga_igop(1,membyte,1,'+')
            if(ga_nodeid()==0)print*,membyte/1000000, 'MB LEFT K_REP'
#endif
                  endif
               endif
            endif


            if(dorepdm) then
               if(xcrepirrg) then
	        xcreplicated_l=nga_create_irreg_config(MT_DBL, 2, dims,
     .                 'replicated DM GA',int_mb(k_map),nblock,plist,
     .                 g_repdm(ii))
	        if(.not.xcreplicated_l) goto 1984

                  call ga_mask_sync(.true.,.false.)
                  call ga_copy(g_dens(ii),g_repdm(ii))
               else
c check memory availability before creat mirrored array(density M)
#ifdef DEBUG
	       call ga_sync()
               membyte=ga_memory_avail()
               call ga_igop(1,membyte,1,'+')
               if(ga_nodeid()==0)print*,membyte/1000000, 'MB ENTER grep'
#endif
                  xcreplicated_l=
     X                 util_mirrmat(1,g_dens(ii),g_repdm(ii),
     T                 .true.,.false.)

                 if(.not.xcreplicated_l) goto 1984
c check memory after
#ifdef DEBUG
	      call ga_sync()
              if(ga_nodeid()==0) print*, 'xcid', xcreplicated_l
              membyte=ga_memory_avail()
              call ga_igop(1,membyte,1,'+')
              if(ga_nodeid()==0)print*,membyte/1000000, 'MB LEFT grep'
#endif  

            endif
          end if

         enddo
c
c     check if everyone is ok
c
 1984 continue
         if(xcrepirrg) then
            if (.not.ma_pop_stack(l_map))
     &           call errquit('xc_rep_init: cannot pop stack',11,
     &       MA_ERR)
         endif

      xcrep_glob=0
      if(xcreplicated_l) xcrep_glob=1
      call ga_igop(375,xcrep_glob,1, '+')
      xcreplicated=xcrep_glob.eq.ga_nnodes()

      if(xcreplicated) then

	 if(oprint.and.ga_nodeid().eq.0) then
            write(luout,*) ' going xcreplicated '
	end if

      else
        if(ga_nodeid().eq.0)
     $   write(luout,*) ' Not enough memory for xcreplicated, back
     $ to computing without mirrored array '
c Destroy the ones that are allocated
	do ii=1,ipol
	  if(xcreplicated_l) then
	   if(.not.ga_destroy(g_repdm(ii))) 
     $     call errquit('grep quit', 2, 0)
	  end if
	 if(xcrep_ga) then
	   if(.not.ga_destroy(k_repxc(ii)))
     $     call errquit('grep quit', 3,0)
	 end if
	end do
      endif
#ifdef DEBUG1
	if(ga_nodeid()==0) call ga_summarize(0)
        call ma_summarize_allocated_blocks()
      write(0,*) ga_nodeid(),' XC DM tot', 
     A       dorepxc,dorepdm,xcreplicated
#endif
      return
      end

      logical function xc_rep_close(rtdb, nbf_ao,ipol,g_vxc,dograd)
c
C$Id: xc_replicated.F,v 1.29 2006-04-25 01:41:21 edo Exp $
c
      implicit none
#include "errquit.fh"
c      
      integer rtdb
c
#include "mafdecls.fh"
#include "global.fh"
#include "dftpara.fh"
#include "dftps.fh"
      integer g_vxc(2),ipol,nbf_ao
      logical dograd ! true if just gradient (not energy)
c
      integer i
c
      xc_rep_close=.true.
      if(xcreplicated) then
         do i = 1, ipol
            if(.not.dograd) then
               if(dorepxc) then
      if (odftps) call pstat_on(ps_xmrge)
                    call ga_mask_sync(.true.,.false.)
                  call ga_merge_mirrored(k_repxc(i))
      if (odftps) call pstat_off(ps_xmrge)
      if (odftps) call pstat_on(ps_xcopy)
                    call ga_mask_sync(.true.,.false.)
                  call ga_copy(k_repxc(i),g_vxc(i))
      if (odftps) call pstat_off(ps_xcopy)
                  if (.not.ga_destroy(k_repxc(i)))
     &              call errquit('xc_rep_close: cannot gadestr',i,
     &       GA_ERR)
               endif
            endif
            if(dorepdm) then
               if (.not.ga_destroy(g_repdm(i)))
     &              call errquit('xc_rep_close: cannot gadestr',i,
     &       GA_ERR)
            endif
         enddo
c	if(ga_nodeid()==0) print*,'**',dorepxc,dorepdm
	call ga_sync()
c     
c     avoid getting things out of control
c
      xcreplicated=.false.
      endif
      return
      end

