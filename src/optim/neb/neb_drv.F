*
* $Id: neb_drv.F,v 1.8 2005-10-03 00:38:47 marat Exp $
*

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
      character*4    mark
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

      nbeads = size_bead_list(bead_list)
      nion   = nion_bead_list(bead_list,1)
      ng     = 3*nion*nbeads

*     *** is verlet algorithm used ***
      if (.not.rtdb_get(rtdb,'neb:verlet',mt_log,1,verlet))
     >  verlet = .false.


*     **** allocate space for gradients and coordinates ****
      value  = MA_alloc_get(mt_dbl,nbeads,'e1',e1(2),e1(1))
      value  = value.and.
     >         MA_alloc_get(mt_dbl,ng,'g0',g0(2),g0(1))
      value  = value.and.
     >         MA_alloc_get(mt_dbl,ng,'g1',g1(2),g1(1))
      value  = value.and.
     >         MA_alloc_get(mt_dbl,ng,'s',s(2),s(1))
      value  = value.and.
     >         MA_alloc_get(mt_dbl,ng,'t1',t1(2),t1(1))
      value  = value.and.
     >         MA_alloc_get(mt_dbl,ng,'v1',v1(2),v1(1))
      value  = value.and.
     >         MA_alloc_get(mt_dbl,ng,'c0',c0(2),c0(1))
      value  = value.and.
     >         MA_alloc_get(mt_dbl,ng,'c1',c1(2),c1(1))
      value  = value.and.
     >         MA_alloc_get(mt_dbl,ng,'dti',dti(2),dti(1))
      if (.not.value) 
     >  call errquit('neb failed - increase stack memory',1,0)

*     *** set dti ***
      value  = MA_alloc_get(mt_dbl,nion,'mass',mass(2),mass(1))      
      if (.not.value) 
     >  call errquit('neb failed - increase stack memory',2,0)
      call neb_masses_get(rtdb,dbl_mb(mass(1)))

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
      value = MA_free_heap(mass(2))      
      if (.not.value) call errquit('neb failed',3,0)
  

*     **** initial step ****
      if (ga_nodeid().eq.0) 
     >   write(*,*) "neb: Calculating Initial Path Energy "
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
      write(*,*) "neb: Initial Path Energy    "
      write(*,*) "neb: -----------------------"
      do i=1,nbeads
         write(*,*) "neb: ",i,dbl_mb(e1(1)+i-1)
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

         IF (ga_nodeid().eq.0)
     >         write(*,*)'neb: iteration #',it
*        *** save old forces  and coordinates ***
         call dcopy(ng,dbl_mb(c1(1)),1,dbl_mb(c0(1)),1)
         call dcopy(ng,dbl_mb(g1(1)),1,dbl_mb(g0(1)),1)

         norm0=norm

          if(verlet) THEN
            IF (ga_nodeid().eq.0)
     >         write(*,*)'neb: using verlet algroithm'
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
            if(oprint) write(*,*)'neb: using cg algorithm'
            call neb_cg_direction(ng,
     >                         dbl_mb(g0(1)),
     >                         dbl_mb(g1(1)),
     >                         dbl_mb(s(1)))

            do ii=1,10

            if(oprint)  write(*,*) "neb: refining time step"
              call neb_move(ng,
     >                      time_step,
     >                      dbl_mb(c0(1)),
     >                      dbl_mb(c1(1)),
     >                      dbl_mb(s(1)))

              call neb_coords_set(bead_list,dbl_mb(c1(1)))
              if(oprint)  write(*,*) "neb: running internal beads"
              call runmid_bead_list(bead_list,task_gradient)
              call neb_energies_get(bead_list,dbl_mb(e1(1)))
              call neb_gradient_get(bead_list,
     >                             dbl_mb(c1(1)),
     >                             dbl_mb(e1(1)),
     >                             dbl_mb(t1(1)),
     >                             dbl_mb(g1(1)))

              norm = dsqrt(ddot(ng,dbl_mb(g1(1)),1,dbl_mb(g1(1)),1))
              if(oprint)  write(*,*) "neb: new gnorm=",norm
              if(oprint)  write(*,*) "neb: old gnorm0=",norm0
              if(norm.gt.norm0) then
                time_step=time_step/2.0d0
                if(oprint)  
     >            write(*,*) "neb: reducing time step",time_step 
              else
                call dscal(ng,time_step,dbl_mb(s(1)),1)
                if(oprint)  
     >            write(*,*) "neb: accepting time step"
                go to 19
              end if
           
            end do
         end if

19       continue         


*        *** RRR write out cumulative path energy
         IF (ga_nodeid().eq.0) THEN
         write(*,*)
         write(*,*) "neb: Path Energy #",it
         write(*,*) "----------------------------"
         do i=1,nbeads
            write(*,*) "neb: ",i,dbl_mb(e1(1)+i-1)
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
           mark = '@neb'
           if (it .gt. 1) mark = ' '
           write(6,1) mark, mark
           mark = '@neb'
           write(6,2) mark, it-1, path_energy, dE,
     $       Gmax, Grms, Xrms, Xmax, util_wallsec()
 1         format(
     $        /,a4,' Step     Path Energy   Delta E   Gmax',
     $        '     Grms     Xrms     Xmax   Walltime',
     $        /,a4,' ---- ---------------- -------- --------',
     $        ' -------- -------- -------- --------')
 2         format(
     $        a4,i5,f17.8,1p,d9.1,0p,4f9.5,f9.1,/,
     $        1x,5x,17x,9x,4a9,/)
         endif

  
      end do

      value = value.and.MA_free_heap(dti(2))
      value = value.and.MA_free_heap(c1(2))
      value = value.and.MA_free_heap(c0(2))
      value = value.and.MA_free_heap(v1(2))
      value = value.and.MA_free_heap(t1(2))
      value = value.and.MA_free_heap(s(2))
      value = value.and.MA_free_heap(g1(2))
      value = value.and.MA_free_heap(g0(2))
      value = value.and.MA_free_heap(e1(2))      
      if (.not.value) call errquit('neb failed',4,0)

      call ga_sync()
      neb = .true. 
      end

