*
* $Id: neb_drv.F,v 1.2 2005-08-08 23:50:37 marat Exp $
*
*   last modified June 2003 by Renee M. Van Ginhoven

      logical function neb(rtdb)
      implicit none
      integer rtdb

#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"

c     
*     **** local variables ****     
      logical value, newchain, svalue,verlet,oprint
      character    mark
      character*50 bead_list
      integer ii
      integer i,it,nbeads,nion,ng,iterNEB,nebsteps
      integer e1(2),g0(2),g1(2),s(2),t1(2),v1(2),c0(2),c1(2)
      integer mass(2),dti(2)
      real*8  path_energy,path_distance,norm,norm0,time_step
      real*8  Gmax,Grms,Xmax,Xrms,dE,path_energy0

*     **** external functions ****
      logical task_gradient
      external task_gradient
      integer  size_bead_list,nion_bead_list
      real*8   energy_bead_list
      external size_bead_list,nion_bead_list
      external energy_bead_list

      oprint = ga_nodeid() .eq. 0

      if (.not.rtdb_get(rtdb,'neb:stepsize',mt_dbl,1,time_step))
     >  time_step = 10.0d0


*  RRR only initialize if this is a new neb chain!
*     **** initialize neb list ****
      bead_list = 'neb_list'
      newchain = .false.
      if (.not.rtdb_get(rtdb,'neb:nebnew',mt_log,1,newchain))
     >  newchain = .true.

      if(newchain) then
        if (ga_nodeid().eq.0)
     >       write(*,*)'NEW NEB CHAIN, INITIALIZING'

        call neb_initialize(rtdb,bead_list)
      else
         if (ga_nodeid().eq.0)
     >       write(*,*)'EXISTING NEB CHAIN? RESTARTING'

      endif
      newchain = .false.


      IF (ga_nodeid().eq.0)
     >       write(*,*) "hera:",bead_list
      nbeads = size_bead_list(bead_list)
      IF (ga_nodeid().eq.0)
     >        write(*,*) "herb:",bead_list
      nion   = nion_bead_list(bead_list,1)
      ng     = 3*nion*nbeads

*     *** is verlet algorithm used ***
      if (.not.rtdb_get(rtdb,'neb:verlet',mt_log,1,verlet))
     >  verlet = .false.


*     **** allocate space for gradients and coordinates ****
      value  = MA_push_get(mt_dbl,nbeads,'e1',e1(2),e1(1))
      value  = value.and.
     >         MA_push_get(mt_dbl,ng,'g0',g0(2),g0(1))
      value  = value.and.
     >         MA_push_get(mt_dbl,ng,'g1',g1(2),g1(1))
      value  = value.and.
     >         MA_push_get(mt_dbl,ng,'s',s(2),s(1))
      value  = value.and.
     >         MA_push_get(mt_dbl,ng,'t1',t1(2),t1(1))
      value  = value.and.
     >         MA_push_get(mt_dbl,ng,'v1',v1(2),v1(1))
      value  = value.and.
     >         MA_push_get(mt_dbl,ng,'c0',c0(2),c0(1))
      value  = value.and.
     >         MA_push_get(mt_dbl,ng,'c1',c1(2),c1(1))
      value  = value.and.
     >         MA_push_get(mt_dbl,ng,'dti',dti(2),dti(1))
      if (.not.value) 
     >  call errquit('neb failed - increase stack memory',1,0)

      IF (ga_nodeid().eq.0)
     >       write(*,*) "herc:",bead_list,nbeads
*     *** set dti ***
      value  = MA_push_get(mt_dbl,nion,'mass',mass(2),mass(1))      
      if (.not.value) 
     >  call errquit('neb failed - increase stack memory',2,0)
      IF (ga_nodeid().eq.0)
     >       write(*,*) "hercc:",bead_list
      call neb_masses_get(rtdb,dbl_mb(mass(1)))
      IF (ga_nodeid().eq.0)
     >      write(*,*) "herccc:",bead_list

      do i=1,nbeads
        do it=1,nion
          dbl_mb(dti(1)+3*nion*(i-1)+3*(it-1)   )
     >    = -(time_step**2)/dbl_mb(mass(1)+it-1)
          dbl_mb(dti(1)+3*nion*(i-1)+3*(it-1) +1)
     >    = -(time_step**2)/dbl_mb(mass(1)+it-1)
          dbl_mb(dti(1)+3*nion*(i-1)+3*(it-1) +2)
     >    = -(time_step**2)/dbl_mb(mass(1)+it-1)
        end do
      end do
      value = MA_pop_stack(mass(2))      
      if (.not.value) call errquit('neb failed',3,0)
  
      IF (ga_nodeid().eq.0)
     >       write(*,*) "herd:",bead_list,nbeads

*     **** initial step ****
      call runall_bead_list(bead_list,task_gradient)
      call neb_energies_get(bead_list,dbl_mb(e1(1)))
      call neb_coords_get(bead_list,dbl_mb(c1(1)))
      call neb_gradient_get(bead_list,
     >                      dbl_mb(c1(1)),
     >                      dbl_mb(e1(1)),
     >                      dbl_mb(t1(1)),
     >                      dbl_mb(g1(1))) 

      call neb_path_energy(bead_list,
     >                     path_distance,
     >                     path_energy) 

      IF (ga_nodeid().eq.0) THEN
      write(*,*)
      write(*,*) "Initial Path Energy         "
      write(*,*) "----------------------------"
      do i=1,nbeads
         write(*,*) i,dbl_mb(e1(1)+i-1)
      end do
      write(*,*)
      END IF
      call create_xyz_file_bead_list(bead_list)

      norm = dsqrt(ddot(ng,dbl_mb(g1(1)),1,dbl_mb(g1(1)),1))
      IF (ga_nodeid().eq.0)
     >       write(*,*) "Path Energy, Path Distance, |G_neb|:",
     >            path_energy,path_distance,norm



      svalue =  rtdb_get(rtdb,'neb:steps',mt_int,1,nebsteps)      
      if(svalue)then
      iterNEB = nebsteps
      IF (ga_nodeid().eq.0)
     >      write(*,*)'NEB iterations  =',iterNEB
      else
         iterNEB = 5
      IF (ga_nodeid().eq.0)
     >      write(*,*)'SHORTRUN, NEB iterations =',iterNEB
      endif

      call dcopy(ng,dbl_mb(g1(1)),1,dbl_mb(s(1)),1)

      do it=1,iterNEB

*        *** save old forces  and coordinates ***
         call dcopy(ng,dbl_mb(c1(1)),1,dbl_mb(c0(1)),1)
         call dcopy(ng,dbl_mb(g1(1)),1,dbl_mb(g0(1)),1)

         norm0=norm

          if(verlet) THEN
            IF (ga_nodeid().eq.0)
     >         write(*,*)'using verlet algroithm'
            call neb_verlet_update(ng,
     >                         dbl_mb(c0(1)),
     >                         dbl_mb(c1(1)),
     >                         dbl_mb(v1(1)),
     >                         dbl_mb(dti(1)),
     >                         dbl_mb(g1(1)))

            call neb_coords_set(bead_list,dbl_mb(c1(1)))
            call runmid_bead_list(bead_list,task_gradient)
            call neb_energies_get(bead_list,dbl_mb(e1(1)))
            call neb_gradient_get(bead_list,
     >                            dbl_mb(c1(1)),
     >                            dbl_mb(e1(1)),
     >                            dbl_mb(t1(1)),
     >                            dbl_mb(g1(1)))

          else 
            if(oprint) write(*,*)'using cg algorithm'
            call neb_cg_direction(ng,
     >                         dbl_mb(g0(1)),
     >                         dbl_mb(g1(1)),
     >                         dbl_mb(s(1)))

            do ii=1,10

            if(oprint)  write(*,*) "refining time step"
              call neb_move(ng,
     >                      time_step,
     >                      dbl_mb(c0(1)),
     >                      dbl_mb(c1(1)),
     >                      dbl_mb(s(1)))

              call neb_coords_set(bead_list,dbl_mb(c1(1)))
              call runmid_bead_list(bead_list,task_gradient)
              call neb_energies_get(bead_list,dbl_mb(e1(1)))
              call neb_gradient_get(bead_list,
     >                             dbl_mb(c1(1)),
     >                             dbl_mb(e1(1)),
     >                             dbl_mb(t1(1)),
     >                             dbl_mb(g1(1)))

              norm = dsqrt(ddot(ng,dbl_mb(g1(1)),1,dbl_mb(g1(1)),1))
              if(oprint)  write(*,*) "norm=",norm
              if(oprint)  write(*,*) "norm0=",norm0
              if(norm.gt.norm0) then
                time_step=time_step/2.0d0
              else
                call dscal(ng,time_step,dbl_mb(s(1)),1)
                go to 19
              end if
           
            end do
         end if

19       continue         


*        *** RRR write out cumulative path energy
         IF (ga_nodeid().eq.0) THEN
         write(*,*)
         write(*,*) "Path Energy #",it
         write(*,*) "----------------------------"
         do i=1,nbeads
            write(*,*) i,dbl_mb(e1(1)+i-1)
         end do
         write(*,*)
         END IF
         call create_xyz_file_bead_list(bead_list)

         path_energy0 = path_energy
         call neb_path_energy(bead_list,
     >                        path_distance,
     >                        path_energy) 
         dE = path_energy - path_energy0
         call neb_calc_convergence(ng,dbl_mb(g1(1)),
     >                                dbl_mb(c0(1)),
     >                                dbl_mb(c1(1)),
     >                                Gmax,Grms,Xmax,Xrms)

!        IF (ga_nodeid().eq.0) THEN
!             write(*,*) "@  Iteration#:",it
!             write(*,*) "@  Path Energy:",path_energy
!             write(*,*) "@  Path Distance:",path_distance
!             write(*,*) "@  |G_neb|:",norm
!        END IF

         if (oprint) then
           mark = '@'
           if (it .gt. 1) mark = ' '
           write(6,1) mark, mark
           mark = '@'
           write(6,2) mark, it-1, path_energy, dE,
     $       Gmax, Grms, Xrms, Xmax, util_wallsec()
 1         format(
     $        /,a1,' Step     Path Energy   Delta E   Gmax',
     $        '     Grms     Xrms     Xmax   Walltime',
     $        /,a1,' ---- ---------------- -------- --------',
     $        ' -------- -------- -------- --------')
 2         format(
     $        a1,i5,f17.8,1p,d9.1,0p,4f9.5,f9.1,/,
     $        1x,5x,17x,9x,4a9,/)
         endif

  
      end do

      value = value.and.MA_pop_stack(dti(2))
      value = value.and.MA_pop_stack(c1(2))
      value = value.and.MA_pop_stack(c0(2))
      value = value.and.MA_pop_stack(v1(2))
      value = value.and.MA_pop_stack(t1(2))
      value = value.and.MA_pop_stack(s(2))
      value = value.and.MA_pop_stack(g1(2))
      value = value.and.MA_pop_stack(g0(2))
      value = value.and.MA_pop_stack(e1(2))      
      if (.not.value) call errquit('neb failed',4,0)

      call ga_sync()
      neb = .true. 
      end

