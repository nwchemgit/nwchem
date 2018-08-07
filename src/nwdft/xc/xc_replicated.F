C
C undefined the first define line if you get undefined symbols
C
      subroutine xc_rep_init(rtdb, geom, basis, g_dens,g_vxc,
     ,     nbf_ao,ipol,goforit,dograd)
c
C$Id$
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
#include "dftps.fh"
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
c
      integer plist,dims(2),nblock(2)
      integer l_map,k_map,procsmp
      integer gnblock(2), membyte
      integer idum,g_ndim,idumm(10),dm_ndim,xc_ndim
      integer ntwo12,howmanyrr,memrylft,memigop
      logical dm_ronly
      integer ilo, ihi, jlo, jhi,jj
c
      oprint = util_print('xcreplicated', print_high)
#ifdef HAVE_SET_GA_PROPERTY
      if (.not.rtdb_get(rtdb,'dft:dm_ronly', mt_log, 1, dm_ronly))
     &        dm_ronly=.false.
c
c     try read-only DM or mirroring
c     
      if(dm_ronly) then
         do ii=1,ipol
            write(6,*) ga_nodeid(), ' set dm ',ii,g_dens(ii)
            call ga_distribution(g_dens(ii),
     .     ga_nodeid(), ilo, ihi, jlo, jhi)
c            call ga_print(g_dens(ii))
            write(6,'(i3,a,4i5)') ga_nodeid(), ' old distr ',
     L           ilo, ihi, jlo, jhi
            call ga_set_property(g_dens(ii),'read_only')
            call ga_distribution(g_dens(ii),
     .     ga_nodeid(), ilo, ihi, jlo, jhi)
            write(6,'(i3,a,4i5)') ga_nodeid(), ' new distr ',
     L           ilo, ihi, jlo, jhi
         enddo
         xcreplicated=.false.
         call ga_sync()
         truerep  =.false.
         truerepdm=.false.
         truerepxc=.false.
         dorepdm  =.false.
         dorepxc  =.false.
      
         return
      endif
#endif
      do ii = 1, ipol
        k_repxc(ii)=-1
        g_repdm(ii)=-1
      enddo
c     non blocking
      if (.not. rtdb_get(rtdb, 'dft:truerepdm', mt_log, 1,
     +  truerepdm)) truerepdm=.false.
      if (.not. rtdb_get(rtdb, 'dft:truerepxc', mt_log, 1,
     +  truerepxc)) truerepxc=.false.
      truerep=truerepdm.or.truerepxc
      if (.not. rtdb_get(rtdb, 'dft:nbacc', mt_log, 1,
     +  dftnbacc)) dftnbacc=.false.
      if (.not. rtdb_get(rtdb, 'dft:nbget', mt_log, 1,
     +  dftnbget)) dftnbget=.false.

c     we do not mirror 3D GAs for now

      call nga_inquire(g_dens(1), idum, dm_ndim, idumm)
      call nga_inquire(g_vxc(1), idum, xc_ndim, idumm)

      xcreplicated=.false.
c     Check number of nodes to decide to use Mirror Array
      if(ga_cluster_nnodes().eq.1.and.(.not.truerep)) return

      if (.not. rtdb_get(rtdb, 'dft:dorepdm', mt_log, 1, 
     &     dorepdm))
     &     dorepdm=.true.
      if (.not. rtdb_get(rtdb, 'dft:dorepxc', mt_log, 1, 
     &     dorepxc))
     &     dorepxc=.true.
      if (rtdb_get(rtdb, 'dft:xcreplicated', mt_log, 1, 
     &     xcreplicated_l)) then
         if(.not.xcreplicated_l.and.(.not.truerep)) return
      else
         xcreplicated_l=.false.
      endif
      if (.not. rtdb_get(rtdb, 'dft:xcrepirrg', mt_log, 1, 
     &     xcrepirrg)) xcrepirrg=.false.
      if(dorepxc) dorepxc=xc_ndim.eq.2
      if(dograd) truerepxc=.false.
      if(truerepxc) truerepxc=xc_ndim.eq.2
      xcreplicated_l=.false.
      xcrep_ga=.false.
      if(truerep) then
c
c         true replicated data aka no mirroring
c
          howmanyrr=0
c         no xc bit in the gradients
          if(truerepxc) howmanyrr=howmanyrr+1
          if(truerepdm) howmanyrr=howmanyrr+1
          ntwo12=(nbf_ao*(nbf_ao+1))/2*ipol
          memrylft=nint(0.92d0*MA_inquire_avail(mt_dbl)-
     -                  howmanyrr*ntwo12)
          memigop=0
          if(memrylft.lt.0) memigop=-memrylft
          call ga_igop(374,memigop,1, 'max')
          if(memigop.gt.0) then
              truerep=.false.
              truerepxc=.false.
              truerepdm=.false.
              xcreplicated=.false.
              if(ga_nodeid().eq.0) then
                  write(luout,*) ' not enough MA memory for Repl.Data'
                  write(luout,1235) int((memigop*
     *                 MA_sizeof(MT_DBL,1,MT_BYTE)*17)/10000000)
                  call util_flush(luout)
 1235             format(' add ',I6,' MBytes')
              endif
              return
          endif
          dorepdm=.false.
          dorepxc=.false.
          dftnbacc=.false.
          dftnbget=.false.
          nbf_ld=nbf_ao
          if(truerepdm) then
c     
c            allocate replicated D matrix
c
             truerepdm=ma_push_get(MT_Dbl,ntwo12*ipol,
     +            'dmat',l_repdm(1),k_repdm(1))
             if(truerepdm) then
                k_repdm(2)=k_repdm(1)+ntwo12/ipol
                if (odftps) call pstat_on(ps_xcopy)
                do ii=1,ipol
                   call xc_bcastdm(geom,basis,nbf_ao,
     +                  g_dens(ii),dbl_mb(k_repdm(ii)))
                enddo
                if (odftps) call pstat_off(ps_xcopy)
             endif
          endif
          if(truerepxc) then
c     
c            allocate replicated Fock matrix
c
             truerepxc=ma_push_get(MT_Dbl,ntwo12,
     .            'xcmat',l_repxc(1),k_repxc(1))
             if(truerepxc) then
                call dcopy(ntwo12,0d0,0,
     D               dbl_mb(k_repxc(1)),1)
                k_repxc(2)=k_repxc(1)+ntwo12/ipol
             endif
          endif
          truerep=truerepdm.or.truerepxc
          xcreplicated_l=truerep

      else
c     
c        use of ga_mirror
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
     +           'ga block maps',l_map,k_map)) then
               call errquit('xcrepl: not enough local mem',
     +              procsmp*2, MA_ERR)
            endif

            call util_atom_blocked(procsmp,
     +           geom, basis, int_mb(k_map), nblock)
            
c     
c           compress map array
c     
            do ii=0,nblock(2)-1
               int_mb(k_map+nblock(1)+ii)=
     +              int_mb(k_map+procsmp+ii)
            enddo
         endif

         do ii=1,ipol
            gnblock(1)=0
            gnblock(2)=0

            if(dorepxc) then
               if(.not.dograd) then
                  if(xcrepirrg) then
                     xcrep_ga=nga_create_irreg_config(MT_DBL, 2, dims,
     +                      'replicatedxcGA',int_mb(k_map),nblock,plist,
     +                      k_repxc(ii))
                     if(.not.xcrep_ga) goto 1984
                     call ga_zero(k_repxc(ii))
                  else
                     xcrep_ga=util_mirrmat(1,g_dens(ii),k_repxc(ii),
     +                                      .false.,.true.)

                     if(.not.xcrep_ga) goto 1984
                  endif
               endif
            else
               xcrep_ga = .true.
            endif


            if(dorepdm) then
               if(xcrepirrg) then
                  xcreplicated_l=nga_create_irreg_config(MT_DBL,2,dims,
     +                 'replicated DM GA',int_mb(k_map),nblock,plist,
     +                 g_repdm(ii))
                  if(.not.xcreplicated_l) goto 1984

                  call ga_mask_sync(.true.,.false.)
                  call ga_copy(g_dens(ii),g_repdm(ii))
               else
                  xcreplicated_l=
     +                 util_mirrmat(1,g_dens(ii),g_repdm(ii),
     +                 .true.,.false.)

                  if(.not.xcreplicated_l) goto 1984

               endif
            else
               xcreplicated_l = .true.
            endif

         enddo
c
c        check if everyone is ok
c
 1984    continue
         if(xcrepirrg) then
            if (.not.ma_pop_stack(l_map))
     &           call errquit('xc_rep_init: cannot pop stack',11,
     &       MA_ERR)
         endif
      endif
      xcrep_glob=0
      if(xcreplicated_l.and.(xcrep_ga.or.truerep)) xcrep_glob=1
      call ga_igop(375,xcrep_glob,1, '+')
      xcreplicated=xcrep_glob.eq.ga_nnodes()

      if(xcreplicated) then

        if(oprint.and.ga_nodeid().eq.0) then
            write(luout,*) ' going xcreplicated '
        end if

c       disable non-blocking stuff
        dftnbacc=.false.
        dftnbget=.false.

      else

        if(oprint.and.ga_nodeid().eq.0)
     $    write(luout,*) ' Not enough memory for xcreplicated, back
     $ to computing without mirrored array '

c       Destroy the ones that are allocated
c       Just attempt the destroy operations. If the array is allocated
c       the destroy operation should succeed, if it is not we do not
c       care that the destroy operation fails.
c
        do ii=1,ipol
          if(.not.ga_destroy(g_repdm(ii))) then
c           we do not need to do anything
          end if
          if(k_repxc(ii).ne.-1) then
          if(.not.ga_destroy(k_repxc(ii))) 
     *            call errquit(' xc_rep: gadestr xc failed',0,GA_ERR)
          end if
        end do
      endif
      return
      end

      logical function xc_rep_close(rtdb, nbf_ao,ipol,
     G     g_vxc,g_dens,dograd)
c
C$Id$
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
#include "rtdb.fh"
      integer g_vxc(2),ipol,nbf_ao
      integer g_dens(2)
      logical dograd ! true if just gradient (not energy)
c
      integer i,j,ilo,ihi,jlo,jhi,numi,numj,icount
      integer idim,jdim,istart
      logical dm_ronly
c
      xc_rep_close=.true.
      if (.not.rtdb_get(rtdb,'dft:dm_ronly', mt_log, 1, dm_ronly))
     &        dm_ronly=.false.
      if(dm_ronly) then
#ifdef HAVE_SET_GA_PROPERTY
         
         do i=1,ipol
            write(6,*) ga_nodeid(), ' unset dm ',i
            call ga_unset_property(g_dens(i))
         enddo
#endif
         xc_rep_close=.true.
      elseif(xcreplicated) then
         if(truerep) then
            if(truerepxc) then
               if (odftps) call pstat_on(ps_xmrge)
               call ga_dgop(374,dbl_mb(k_repxc(1)),
     N              (nbf_ao*(nbf_ao+1))/2*ipol,'+')
               do i = 1, ipol
                  call ga_zero(g_vxc(i))
                  call xc_putget12('put',g_vxc(i),
     N                 ga_nodeid(),nbf_ao,dbl_mb(k_repxc(i)))
               end do
            endif
            if(.not.truerepxc) l_repxc(1)=l_repdm(1)
            xc_rep_close=ma_chop_stack(l_repxc(1))
            if(ga_nodeid().eq.0.and.(.not.xc_rep_close)) then
               write(0,*) ' lrepxc handle no. ',l_repxc(1)
               write(0,*) ' lrepdm handle no. ',l_repdm(1)
               call ma_summarize_allocated_blocks()
            endif
            if(truerepxc) then
               call ga_sync()
               if (odftps) call pstat_off(ps_xmrge)
            endif
c
c     cleanup of logical to avoid mayhem in later calls
c
            truerepxc=.false.
            truerepdm=.false.
            truerep=.false.
            xcreplicated=.false.
            dorepdm=.false.
            dorepxc=.false.
         else
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
                     if (.not.ga_destroy(k_repxc(i))) call 
     E                    errquit('xc_rep_close: cannot gadestr',i,
     &                    GA_ERR)
                  endif
               endif
               if(dorepdm) then
                  if (.not.ga_destroy(g_repdm(i)))
     &                 call errquit('xc_rep_close: cannot gadestr',i,
     &                 GA_ERR)
               endif
            enddo
         endif
         
         if(.not.truerepxc) call ga_sync()
c     
c     avoid getting things out of control
c     
         xcreplicated=.false.
      endif
      return
      end
      subroutine xc_bcastdm(geom,basis,nbf,
     A              g_a,a)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer geom,basis,nbf,g_a,mysize,nnn
      double precision a(*)
      integer l_acopy,k_acopy,i,j,idim,icount,ica
      integer proc, myproc, ilo, ihi, jlo, jhi,jj
c
      myproc = ga_nodeid() 
      proc   = ga_nnodes()
      nnn=(nbf*(nbf+1))/2
      call dcopy(nnn,0d0,0,a,1)
      call xc_putget12('get',g_a,
     N     myproc,nbf,a)
      call ga_dgop(374,a, nnn,'+')
      return
      end
      subroutine xc_putget12(what,g_a,me,nbf_ao,a)
      implicit none
      character*3 what
      double precision a(*)
      integer g_a, me,nbf_ao
c
      integer ilo, ihi, jlo, jhi
      integer idim,jdim,j,istart,icount
c      
      call ga_distribution(g_a,
     .     me, ilo, ihi, jlo, jhi)
      if (ilo.gt.0 .and.jlo.gt.0.and.
     A    ilo.le.ihi.and.jlo.le.ihi) then
c     lower triangle
         idim =  ihi-ilo+1
         jdim =  jhi-jlo+1
         if (idim.gt.0 .and. jdim.gt.0) then
            do j=jlo,jhi
               if(j.le.ihi) then
c     lower triangle
coffdiag
                  if(j.lt.ilo) then
                     istart=ilo
                     icount=((j-1)*(2*(nbf_ao+1)-j))/2+
     P                    istart-j+1
                  else
c     diag
                     istart=j
                     icount=((j-1)*(2*(nbf_ao+1)-j))/2+1
                  endif
                  if(what.eq.'get') then
                     call ga_get(g_a,istart,
     N                    ihi,j,j, 
     $                    a(icount),1)
                  else
                     call ga_put(g_a,istart,
     N                    ihi,j,j, 
     $                    a(icount),1)
                  endif
               endif
            enddo
         endif
      endif
      return
      end
