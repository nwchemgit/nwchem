*
* $Id$
*

      logical function neb(rtdb)
      implicit none
      integer rtdb

#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"

c     
*     **** local variables ****     
      logical value, newchain, svalue,verlet,oprint,done
      character*4    mark
      character*50 bead_list
      character*80 title,neb_movecs
      integer ii
      integer i,it,nbeads,nion,ng,nebsteps
      integer e1(2),g0(2),g1(2),s(2),t1(2),v1(2),c0(2),c1(2)
      integer mass(2),dti(2)
      real*8  path_energy,path_distance,norm,norm0,time_step,kbeads
      real*8  Gmax,Grms,Xmax,Xrms,dE,path_energy0
      real*8  nebGmax,nebGrms,nebXmax,nebXrms,nebdE

*     **** external functions ****
      logical task_gradient
      external task_gradient
      integer  size_bead_list,nion_bead_list
      real*8   energy_bead_list
      external size_bead_list,nion_bead_list
      external energy_bead_list

      oprint = ga_nodeid() .eq. 0

      if (oprint) then
         write(luout,*)
         write(luout,*)
         call util_print_centered(luout,
     >        'NWChem Minimum Energy Pathway Program (NEB)',
     >        40,.true.)
         write(luout,*)
         write(luout,*)
      endif

      if (rtdb_cget(rtdb,'title',1,title)) then
         if (oprint) then
            write(luout,*)
            write(luout,*)
            call util_print_centered(6, title, 40, .false.)
            write(luout,*)
            write(luout,*)
         endif
      endif
      if (.not. rtdb_get(rtdb,'neb:gmax',mt_dbl,1,nebgmax))
     >   nebgmax = 0.00045d0
      if (.not. rtdb_get(rtdb,'neb:grms',mt_dbl,1,nebgrms))
     >   nebgrms = 0.00030d0
      if (.not. rtdb_get(rtdb,'neb:xmax',mt_dbl,1,nebxmax))
     >   nebxmax = 0.00180d0
      if (.not. rtdb_get(rtdb,'neb:xrms',mt_dbl,1,nebxrms))
     >   nebxrms = 0.00120d0

      if (.not.rtdb_get(rtdb,'neb:stepsize',mt_dbl,1,time_step))
     >  time_step = 10.0d0
      if (.not.rtdb_get(rtdb,'neb:kbeads',mt_dbl,1,kbeads))
     >  kbeads = 0.1d0
      if (.not.rtdb_get(rtdb,'neb:steps',mt_int,1,nebsteps))      
     >   nebsteps = 5

      if (.not.rtdb_cget(rtdb,'neb:movecs',1,neb_movecs)) then
         call util_file_prefix('movecs',neb_movecs)
      end if


*  RRR only initialize if this is a new neb chain!
*     **** initialize neb list ****
      bead_list = 'neb_list'
      newchain = .false.
      if (.not.rtdb_get(rtdb,'neb:nebnew',mt_log,1,newchain))
     >  newchain = .true.

      if(newchain) then
        if (oprint) write(*,*)'NEW NEB CHAIN, INITIALIZING'
        call neb_initialize(rtdb,bead_list)
      else
         if (oprint) write(*,*)'EXISTING NEB CHAIN? RESTARTING'
      endif
      newchain = .false.

      nbeads = size_bead_list(bead_list)
      nion   = nion_bead_list(bead_list,1)
      ng     = 3*nion*nbeads

*     *** is verlet algorithm used ***
      if (.not.rtdb_get(rtdb,'neb:verlet',mt_log,1,verlet))
     >  verlet = .false.


      if (oprint) then
         write(luout,1) nebgmax,nebgrms,nebxmax,nebxrms,
     >                  time_step,nebsteps,nbeads,kbeads,
     >                  neb_movecs
         write(luout,*)
         write(luout,*)
 1       format(
     >        ' maximum gradient threshold         (gmax) = ', f10.6,/,
     >        ' rms gradient threshold             (grms) = ', f10.6,/,
     >        ' maximum cartesian step threshold   (xmax) = ', f10.6,/,
     >        ' rms cartesian step threshold       (xrms) = ', f10.6,/,
     >        0p,/,
     >        ' step size                       (tepsize) = ', f10.6,/,
     >        ' maximum number of steps         (maxiter) = ', i4,/,
     >        ' number of images in path         (nbeads) = ', i4,/,
     >        ' NEB spring constant in a.u.      (kbeads) = ', f10.6,/,
     >        ' NEB movecs filename                       = ', a)
      end if


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
      if (oprint) write(luout,*) "neb: Calculating Initial Path Energy"
      call runall_bead_list(bead_list,task_gradient)
      call neb_energies_get(bead_list,dbl_mb(e1(1)))
      call neb_coords_get(bead_list,dbl_mb(c1(1)))
      call neb_gradient_get(bead_list,kbeads,
     >                      dbl_mb(c1(1)),
     >                      dbl_mb(e1(1)),
     >                      dbl_mb(t1(1)),
     >                      dbl_mb(g1(1))) 

      call neb_path_energy(bead_list,
     >                     path_distance,
     >                     path_energy) 

      if (oprint) then
      write(luout,*)
      write(luout,*) "neb: Initial Path Energy    "
      write(luout,*) "neb: -----------------------"
      do i=1,nbeads
         write(luout,*) "neb: ",i,dbl_mb(e1(1)+i-1)
      end do
      write(luout,*)
      END IF
      call create_xyz_file_bead_list(bead_list)

      norm = dsqrt(ddot(ng,dbl_mb(g1(1)),1,dbl_mb(g1(1)),1))
      if (oprint) write(luout,*) "Path Energy, Path Distance, |G_neb|:",
     >            path_energy,path_distance,norm

      call dcopy(ng,dbl_mb(g1(1)),1,dbl_mb(s(1)),1)

      it = 0
      done = .false.
      do while (.not.done)
         it = it + 1
         if (oprint) write(luout,*) 'neb: iteration #',it

*        *** save old forces  and coordinates ***
         call dcopy(ng,dbl_mb(c1(1)),1,dbl_mb(c0(1)),1)
         call dcopy(ng,dbl_mb(g1(1)),1,dbl_mb(g0(1)),1)

         norm0=norm

          if(verlet) THEN
            IF (oprint) write(luout,*)'neb: using verlet algroithm'
            call neb_verlet_update(ng,
     >                         dbl_mb(c0(1)),
     >                         dbl_mb(c1(1)),
     >                         dbl_mb(v1(1)),
     >                         dbl_mb(dti(1)),
     >                         dbl_mb(g1(1)))

            call neb_coords_set(bead_list,dbl_mb(c1(1)))
            call runmid_bead_list(bead_list,task_gradient)
            call neb_energies_get(bead_list,dbl_mb(e1(1)))
            call neb_gradient_get(bead_list,kbeads,
     >                            dbl_mb(c1(1)),
     >                            dbl_mb(e1(1)),
     >                            dbl_mb(t1(1)),
     >                            dbl_mb(g1(1)))

          else 
            if (oprint) write(luout,*) 'neb: using cg algorithm'
            call neb_cg_direction(ng,
     >                         dbl_mb(g0(1)),
     >                         dbl_mb(g1(1)),
     >                         dbl_mb(s(1)))

            do ii=1,10

            if (oprint) write(luout,*) "neb: refining time step"
              call neb_move(ng,
     >                      time_step,
     >                      dbl_mb(c0(1)),
     >                      dbl_mb(c1(1)),
     >                      dbl_mb(s(1)))

              call neb_coords_set(bead_list,dbl_mb(c1(1)))
              if (oprint) write(luout,*) "neb: running internal beads"
              call runmid_bead_list(bead_list,task_gradient)
              call neb_energies_get(bead_list,dbl_mb(e1(1)))
              call neb_gradient_get(bead_list,kbeads,
     >                             dbl_mb(c1(1)),
     >                             dbl_mb(e1(1)),
     >                             dbl_mb(t1(1)),
     >                             dbl_mb(g1(1)))

              norm = dsqrt(ddot(ng,dbl_mb(g1(1)),1,dbl_mb(g1(1)),1))
              if(oprint)  write(luout,*) "neb: new gnorm=",norm
              if(oprint)  write(luout,*) "neb: old gnorm0=",norm0
              if(norm.gt.norm0) then
                time_step=time_step/2.0d0
                if (oprint) 
     >             write(luout,*) "neb: reducing time step ",time_step
              else
                call dscal(ng,time_step,dbl_mb(s(1)),1)
                if (oprint) 
     >             write(luout,*) "neb: accepting time step ",time_step
                go to 19
              end if
            end do
         end if
19       continue         


*        *** RRR write out cumulative path energy
         if (oprint) then
            write(luout,*)
            write(luout,*) "neb: Path Energy #",it
            write(luout,*) "----------------------------"
            do i=1,nbeads
               write(luout,*) "neb: ",i,dbl_mb(e1(1)+i-1)
            end do
            write(luout,*)
         end if
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

         if (oprint) then
           mark = '@neb'
           if (it .gt. 1) mark = ' '
           write(luout,10) mark, mark
           mark = '@neb'
           write(luout,11) mark, it, path_energy, dE,
     $       Gmax, Grms, Xrms, Xmax, util_wallsec()
 10        format(
     $        /,a4,' Step     Path Energy   Delta E   Gmax',
     $        '     Grms     Xrms     Xmax   Walltime',
     $        /,a4,' ---- ---------------- -------- --------',
     $        ' -------- -------- -------- --------')
 11        format(
     $        a4,i5,f17.8,1p,d9.1,0p,4f9.5,f9.1,/,
     $        1x,5x,17x,9x,4a9,/)
         endif

        done =          (Gmax.le.nebgmax)
        done = done.and.(Grms.le.nebgrms)
        done = done.and.(Xmax.le.nebxmax)
        done = done.and.(Xrms.le.nebxrms)
        done = done.or.(it.ge.nebsteps)
  
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

      if (.not.rtdb_put(rtdb,'neb:nebnew',mt_log,1,.false.))
     > call errquit('setting neb:nebnew failed',4,RTDB_ERR)

      call ga_sync()
      neb = .true. 
      end

