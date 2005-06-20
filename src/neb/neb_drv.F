*
* $Id: neb_drv.F,v 1.1 2005-06-20 06:17:05 marat Exp $
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
      call runall_bead_list(bead_list)
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
*      stop

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
            call runmid_bead_list(bead_list)
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
              call runmid_bead_list(bead_list)
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

*     *** RRR  SUBRT TO GET PATH ENERGY

      subroutine neb_path_energy(bead_list,dpath,epath)
      implicit none
      character*(*) bead_list
      real*8 dpath
      real*8 epath

#include "mafdecls.fh"

*     **** local variables ****
      logical value
      integer i,index,index_p,nbeads,nion,ng
      integer e(2),c(2),r(2)
      real*8  dist
     
*     **** external functions ****
      integer  size_bead_list,nion_bead_list
      real*8   ddot
      external size_bead_list,nion_bead_list
      external ddot

      nbeads = size_bead_list(bead_list)
      nion   = nion_bead_list(bead_list,1)
      ng     = 3*nion*nbeads

      value  = MA_push_get(mt_dbl,nbeads,'e',e(2),e(1))
      value  = value.and.
     >         MA_push_get(mt_dbl,3*nion,'r',r(2),r(1))
      value  = value.and.
     >         MA_push_get(mt_dbl,ng,'c',c(2),c(1))
      if (.not.value) 
     > call errquit('neb_path_energy failed - increase stack',0,0)
    
      call neb_energies_get(bead_list,dbl_mb(e(1)))
      call neb_coords_get(bead_list,dbl_mb(c(1)))
      
*     **** calculate the path length and path energy ****
      dpath = 0.0d0
      epath = 0.0d0
      do i=1,(nbeads-1)
         index   = (i-1)*3*nion + 1
         index_p = (i  )*3*nion + 1
         call dcopy(3*nion,dbl_mb(c(1)+index_p-1),1,
     >                     dbl_mb(r(1)),1)
         call daxpy(3*nion,(-1.0d0),
     >              dbl_mb(c(1)+index-1),1,
     >              dbl_mb(r(1)),   1)
         dist = dsqrt(ddot(3*nion,dbl_mb(r(1)),1,
     >                            dbl_mb(r(1)),1))
         dpath = dpath + dist
         epath = epath + 0.5d0*dbl_mb(e(1)+i-1)*dist
         epath = epath + 0.5d0*dbl_mb(e(1)+i  )*dist
      end do
      value = value.and.MA_pop_stack(c(2))
      value = value.and.MA_pop_stack(r(2))
      value = value.and.MA_pop_stack(e(2))      
      if (.not.value) call errquit('neb_path_energy failed',1,0)

      return
      end



*     *** RRR END OF ENERGY SUBROUTINE  ***

      subroutine neb_gradient_get(bead_list,c,e,t,g)
      implicit none
      character*(*) bead_list
      real*8 c(*)
      real*8 e(*)
      real*8 t(*)
      real*8 g(*)

*     **** local variables ****
      integer i,index,nbeads,nion
      real*8  k,norm
     
*     **** external functions ****
      integer  size_bead_list,nion_bead_list
      real*8   ddot
      external size_bead_list,nion_bead_list
      external ddot


      k = 0.1d0 !*** what's a good value?

      nbeads = size_bead_list(bead_list)
      nion   = nion_bead_list(bead_list,1)

*     *** get the tangent ****
      call neb_tangent(nbeads,nion,c,e,t)

*     **** normal gradient |g> = |g> - |t><t|g> ****
      call dcopy(nbeads*3*nion,0.0d0,0,g,1)
      do i=2,(nbeads-1)
         index = (i-1)*3*nion + 1
         call gradient_get_bead_list(bead_list,i,g(index))
         norm = ddot(3*nion,t(index),1,g(index),1)
         call daxpy(3*nion,(-norm),
     >              t(index),1,
     >              g(index),1)
      end do

*     *** add to normal gradient the tangent spring gradient *** 
      call neb_add_spring_gradient(nbeads,nion,c,t,k,g)
      
      return
      end


      subroutine neb_add_spring_gradient(nbeads,nion,c,t,k,gs)
      implicit none
      integer nbeads,nion
      real*8 c(*)
      real*8 t(*)
      real*8 k
      real*8 gs(*)

*     **** local variables ****
      integer i,index,index_m,index_p
      integer rp,rm
      real*8  normm,normp

*     **** external functions ***
      real*8   ddot
      external ddot


      rm = 1
      rp = (nbeads-1)*3*nion + 1
      
      do i=2,(nbeads-1)
         index   = (i-1)*3*nion + 1
         index_m = (i-2)*3*nion + 1
         index_p = (i  )*3*nion + 1
         call dcopy(3*nion,c(index),1,t(rm),1)
         call daxpy(3*nion,(-1.0d0),
     >               c(index_m),1,
     >               t(rm),     1)
         call dcopy(3*nion,c(index_p),1,t(rp),1)
         call daxpy(3*nion,(-1.0d0),
     >               c(index),1,
     >               t(rp),   1)
         normm = ddot(3*nion,t(rm),1,t(rm),1)
         normp = ddot(3*nion,t(rp),1,t(rp),1)
         normp = -k*(dsqrt(normp) - dsqrt(normm))
         call daxpy(3*nion,normp,
     >              t(index), 1,
     >              gs(index),1)
      end do
      call dcopy(3*nion,0.0d0,0,t(rm),1)
      call dcopy(3*nion,0.0d0,0,t(rp),1)

      return
      end

      subroutine neb_tangent(nbeads,nion,c,e,t)
      implicit none
      integer nbeads,nion
      real*8 c(*)
      real*8 e(*)
      real*8 t(*)

*     **** local variables ****
      integer i,index,index_m,index_p
      integer rp,rm
      real*8  norm,dVmax,dVmin

*     **** external functions ***
      real*8   ddot
      external ddot

      rm = 1
      rp = (nbeads-1)*3*nion + 1
      
      do i=2,(nbeads-1)
         index   = (i-1)*3*nion + 1
         index_m = (i-2)*3*nion + 1
         index_p = (i  )*3*nion + 1
         call dcopy(3*nion,c(index),1,t(rm),1)
         call daxpy(3*nion,(-1.0d0),
     >               c(index_m),1,
     >               t(rm),     1)
         call dcopy(3*nion,c(index_p),1,t(rp),1)
         call daxpy(3*nion,(-1.0d0),
     >               c(index),1,
     >               t(rp),   1)

         if      ((e(i+1).gt.e(i)).and.(e(i).gt.e(i-1))) then
            call dcopy(3*nion,t(rp),1,t(index),1)
         else if ((e(i-1).gt.e(i)).and.(e(i).gt.e(i+1))) then
            call dcopy(3*nion,t(rm),1,t(index),1)
         else
            
             if (dabs(e(i+1)-e(i)).gt.dabs(e(i-1)-e(i))) then
               dVmax = dabs(e(i+1)-e(i))
               dVmin = dabs(e(i-1)-e(i))
             else
               dVmax = dabs(e(i-1)-e(i))
               dVmin = dabs(e(i+1)-e(i))
             end if

             if (e(i+1).gt.e(i-1)) then
               call dscal(3*nion,dVmax,t(rp),1)
               call dscal(3*nion,dVmin,t(rm),1)
             else
               call dscal(3*nion,dVmin,t(rp),1)
               call dscal(3*nion,dVmax,t(rm),1)
             end if
             call dcopy(3*nion,t(rp),1,t(index),1)
             call daxpy(3*nion,(1.0d0),
     >                  t(rm),1,
     >                  t(index),1)
         end if

*        *** normalize tangent ***
         norm = ddot(3*nion,t(index),1,t(index),1)
         norm = 1.0d0/dsqrt(norm)
         call dscal(3*nion,norm,t(index),1)

      end do
      call dcopy(3*nion,0.0d0,0,t(rm),1)
      call dcopy(3*nion,0.0d0,0,t(rp),1)

      return
      end




*     *** RRR SUBRT TO GET AND ADJUST GRADIENT ***
*     *** RRR CHANGE TO INCLUDE CLIBING IMAGE ***
*     *** RRR CHANGE COMBINE PROJECTION and sforce
*     *** RRR addition into one subroutine

      subroutine neb_gradient_get1(bead_list,c,e,t,g)


      implicit none
      character*(*) bead_list
      real*8 c(*)
      real*8 e(*)
      real*8 t(*)
      real*8 g(*)

#include "rtdb.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "global.fh"

*     **** local variables ****
      integer i,index,index_m, index_p
      integer nbeads, nion, rm, rp
      real*8  k,norm, normm, normp
      real*8  sprf, sprnorm
     
*     **** external functions ****
      integer  size_bead_list,nion_bead_list
      real*8   ddot
      external size_bead_list,nion_bead_list
      external ddot

*    *** RRR SET SPRING CONSTANT HERE
*    *** user-supplied k, default 0.1d0
*    *** currently not working
*    *** use default
*      ispring = rtdb_get(rtdb, 'spring', mt_dbl,1,spring)
*      write(*,*)'is spring constant set? =',ispring,spring
*      if(ispring) then
*         k = spring
*         write(*,*)'NEB spring constant set to  =',k
*      else
         k = 0.1d0
*         write(*,*)'NEB default springs =',k
*      endif      
*      stop





      nbeads = size_bead_list(bead_list)
      nion   = nion_bead_list(bead_list,1)

      rm = 1
      rp = (nbeads-1)*3*nion + 1


*     *** get the tangent ****
      call neb_tangent(nbeads,nion,c,e,t)


*    *** RRR project out tangent PES force

*     **** normal gradient |g> = |g> - |t><t|g> ****
      call dcopy(nbeads*3*nion,0.0d0,0,g,1)

*    *** RRR loop over ions to be moved
      do i=2,(nbeads-1)

*    *** RRR set index numbers
         index = (i-1)*3*nion + 1
         index_m = (i-2)*3*nion + 1
         index_p = (i  )*3*nion + 1

*    *** RRR get original gradients
         call gradient_get_bead_list(bead_list,i,g(index))

*    *** RRR dot product of 
*    *** RRR original force and tangent vector
*    *** RRR to get magnitude of gradient
*    *** RRR  parallel to tangent direction
         norm = ddot(3*nion,t(index),1,g(index),1)

*    *** RRR check relative energy of current bead
*    *** RRR if both nneighbors lower, then
*    *** RRR invert parallel force component
*    *** RRR to climb uphill along MEP
         if ((e(i).gt.e(i-1)).and.(e(i).gt.e(i+1))) then
                  call daxpy(3*nion,2*(-norm),
     >              t(index),1,
     >              g(index),1)
         else
*    *** MV for other images use regular projection ***
             call daxpy(3*nion,-norm,
     >       t(index),1,
     >       g(index),1)
          

*    *** RRR calculate spring forces
                  
                 call dcopy(3*nion,c(index),1,t(rm),1)
                 call daxpy(3*nion,(-1.0d0),
     >               c(index_m),1,
     >               t(rm),     1)
                 call dcopy(3*nion,c(index_p),1,t(rp),1)
                 call daxpy(3*nion,(-1.0d0),
     >               c(index),1,
     >               t(rp),   1)
                 normm = ddot(3*nion,t(rm),1,t(rm),1)
                 normp = ddot(3*nion,t(rp),1,t(rp),1)
                 sprf = -k*(dsqrt(normp) - dsqrt(normm))
                 sprnorm = sprf - norm
                 call daxpy(3*nion,sprnorm,
     >              t(index), 1,
     >              g(index),1)
         endif
        
      end do

*   *** RRR END OF FORCE MODIFICATION
      call dcopy(3*nion,0.0d0,0,t(rm),1)
      call dcopy(3*nion,0.0d0,0,t(rp),1)

      return
      end


*    *** RRR GET TANGENT TO BEAD PATH
*    *** RRR  DEFINE TANGENT AS THE VECTOR
*    *** RRR   TO THE NEIGHBOR BEAD THAT IS
*    *** RRR   HIGHER IN ENERGY
*    *** RRR  IE COMPARE e(i-1) to e(e+1)
*    *** RRR  if e(i-1) > e(e+1) then
*    *** RRR  tan goes from e(i+1) to e(e-1)
    

      subroutine neb_tangent1(nbeads,nion,c,e,t)
      implicit none
      integer nbeads,nion
      real*8 c(*)
      real*8 e(*)
      real*8 t(*)

*     **** local variables ****
      integer i,index,index_m,index_p
      integer rp,rm
      real*8  norm,dVmax,dVmin

*     **** external functions ***
      real*8   ddot
      external ddot

      rm = 1
      rp = (nbeads-1)*3*nion + 1
      
      do i=2,(nbeads-1)
         index   = (i-1)*3*nion + 1
         index_m = (i-2)*3*nion + 1
         index_p = (i  )*3*nion + 1
         call dcopy(3*nion,c(index),1,t(rm),1)
         call daxpy(3*nion,(-1.0d0),
     >               c(index_m),1,
     >               t(rm),     1)
         call dcopy(3*nion,c(index_p),1,t(rp),1)
         call daxpy(3*nion,(-1.0d0),
     >               c(index),1,
     >               t(rp),   1)

         if      ((e(i+1).gt.e(i)).and.(e(i).gt.e(i-1))) then
            call dcopy(3*nion,t(rp),1,t(index),1)
         else if ((e(i-1).gt.e(i)).and.(e(i).gt.e(i+1))) then
            call dcopy(3*nion,t(rm),1,t(index),1)
         else
            
             if (dabs(e(i+1)-e(i)).gt.dabs(e(i-1)-e(i))) then
               dVmax = dabs(e(i+1)-e(i))
               dVmin = dabs(e(i-1)-e(i))
             else
               dVmax = dabs(e(i-1)-e(i))
               dVmin = dabs(e(i+1)-e(i))
             end if

             if (e(i+1).gt.e(i-1)) then
               call dscal(3*nion,dVmax,t(rp),1)
               call dscal(3*nion,dVmin,t(rm),1)
             else
               call dscal(3*nion,dVmin,t(rp),1)
               call dscal(3*nion,dVmax,t(rm),1)
             end if
             call dcopy(3*nion,t(rp),1,t(index),1)
             call daxpy(3*nion,(1.0d0),
     >                  t(rm),1,
     >                  t(index),1)
         end if

*        *** normalize tangent ***
         norm = ddot(3*nion,t(index),1,t(index),1)
         norm = 1.0d0/dsqrt(norm)
         call dscal(3*nion,norm,t(index),1)

      end do
      call dcopy(3*nion,0.0d0,0,t(rm),1)
      call dcopy(3*nion,0.0d0,0,t(rp),1)

      return
      end



      subroutine neb_energies_get(bead_list,e)
      implicit none
      character*(*) bead_list
      real*8 e(*)

*     **** local variables ****
      integer i,nbeads

*     **** external functions ****
      integer  size_bead_list
      real*8   energy_bead_list
      external size_bead_list
      external energy_bead_list

      nbeads = size_bead_list(bead_list)
    
      do i=1,nbeads
         e(i) = energy_bead_list(bead_list,i)
      end do

      return
      end


      subroutine neb_coords_get(bead_list,c)
      implicit none
      character*(*) bead_list
      real*8 c(*)

*     **** local variables ****
      integer i,index,nbeads,nion

*     **** external functions ****
      integer  size_bead_list,nion_bead_list
      external size_bead_list,nion_bead_list

      nbeads = size_bead_list(bead_list)
      nion   = nion_bead_list(bead_list,1)

      do i=1,(nbeads)
        index = (i-1)*3*nion+1
        call coords_get_bead_list(bead_list,i,c(index))
      end do

      return
      end

      subroutine neb_coords_set(bead_list,c)
      implicit none
      character*(*) bead_list
      real*8 c(*)

*     **** local variables ****
      integer i,index,nbeads,nion

*     **** external functions ****
      integer  size_bead_list,nion_bead_list
      external size_bead_list,nion_bead_list

      nbeads = size_bead_list(bead_list)
      nion   = nion_bead_list(bead_list,1)

      do i=1,nbeads
        index = (i-1)*3*nion+1
        call coords_set_bead_list(bead_list,i,c(index))
      end do


      return
      end

      subroutine  neb_masses_get(rtdb,m)
      implicit none
      integer rtdb
      real*8 m(*)

#include "geom.fh"

*     **** local variables ****
      logical value
      integer geom,nion

      value = geom_create(geom,'neb_tmp')
      value = value.and.geom_rtdb_load(rtdb,geom,'neb_start')
      value = value.and.geom_ncent(geom,nion)
      value = value.and.geom_masses_get(geom,nion,m)
      value = value.and.geom_destroy(geom)
      if (.not.value) call errquit('neb_masses failed',0,0)

      call dscal(nion,1822.89d0,m,1)

      return
      end


      subroutine neb_initialize(rtdb, bead_list)
      implicit none
      integer rtdb
      character*(*) bead_list 

#include "nwc_const.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "inp.fh"
c     
c     This routine initializes the common /coptopt/ and
c     also creates and returns the geometry handle
c     
c      integer i, j, num, ma_type, nactive_atoms, l_actlist
      integer nbeads
c      logical ignore
      logical oprint
      character*80 title,neb_movecs
c      character*8 source, test
c      character*32 theory
c     
      oprint = util_print('information', print_low)
     $     .and. (ga_nodeid() .eq. 0)
c
      if (oprint) then
         write(6,*)
         write(6,*)
         call util_print_centered(6,
     $        'NWChem Minimum Energy Pathway Program',
     $        40,.true.)
         write(6,*)
         write(6,*)
      endif
c
      if (rtdb_cget(rtdb,'title',1,title)) then
         if (oprint) then
            write(6,*)
            write(6,*)
            call util_print_centered(6, title, 40, .false.)
            write(6,*)
            write(6,*)
         endif
      endif

      if (.not.rtdb_cget(rtdb,'neb:movecs',1,neb_movecs)) then
         call util_file_prefix('movecs',neb_movecs)
      end if
 
      if (.not.rtdb_get(rtdb,'neb:nbeads',mt_int,1,nbeads)) then
         nbeads = 5
      end if
      if (oprint) then
         write(6,1) 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     >              0,nbeads,
     >              neb_movecs
 1       format(
     $        ' maximum gradient threshold         (gmax) = ', f10.6,/,
     $        ' rms gradient threshold             (grms) = ', f10.6,/,
     $        ' maximum cartesian step threshold   (xmax) = ', f10.6,/,
     $        ' rms cartesian step threshold       (xrms) = ', f10.6,/,
     $        ' fixed trust radius                (trust) = ', f10.6,/,
     $        ' energy precision                  (eprec) = ', 1p,d9.1,
     $        0p,/,
     $        ' maximum number of steps          (nptopt) = ', i4,/,
     $        ' number of images in path         (nbeads) = ', i4,/,
     $        ' NEB movecs filename                       = ', a)
      end if


*     **** create bead_list *** 
      call init_bead_list(rtdb,bead_list,neb_movecs)
      call neb_initial_path(rtdb,bead_list,nbeads)

      return
      end


      subroutine neb_initial_path(rtdb,bead_list,nbeads)
      implicit none
      integer rtdb
      character*(*) bead_list
      integer nbeads

#include "mafdecls.fh"
#include "geom.fh"

*     **** local variables ****
      logical value
      integer i,geom,geomlen,movecslen,nion
      integer r1(2),r2(2),r3(2)
      real*8  t
      character*255 geom_name,movecs_name

*     **** external functions ****
      integer     inp_strlen
      character*7 bead_index_name
      external    inp_strlen
      external    bead_index_name

      value = geom_create(geom,'neb_tmp')
      value = value.and.geom_rtdb_load(rtdb,geom,'neb_start')
      value = value.and.geom_ncent(geom,nion)
      value = value.and.geom_destroy(geom)
      if (.not.value) call errquit('neb_initial_path failed',0,0)

      value = value.and.MA_push_get(mt_dbl, (3*nion), 'r1',r1(2),r1(1))
      value = value.and.MA_push_get(mt_dbl, (3*nion), 'r2',r2(2),r2(1))
      value = value.and.MA_push_get(mt_dbl, (3*nion), 'r3',r3(2),r3(1))
      if (.not.value) call errquit('neb_initial_path failed',1,0)

      value = value.and.geom_create(geom,'neb_tmp')
      value = value.and.geom_rtdb_load(rtdb,geom,'neb_end')
      value = value.and.geom_cart_coords_get(geom,dbl_mb(r2(1)))
      value = value.and.geom_rtdb_load(rtdb,geom,'neb_start')
      value = value.and.geom_cart_coords_get(geom,dbl_mb(r1(1)))
      if (.not.value) call errquit('neb_initial_path failed',2,0)

      do i=1,nbeads
        t = (i-1)/dble(nbeads-1)

        call dcopy(3*nion,dbl_mb(r1(1)),1,dbl_mb(r3(1)),1)
        call dscal(3*nion,(1.0d0-t),dbl_mb(r3(1)),1)
        call daxpy(3*nion,t,dbl_mb(r2(1)),1,dbl_mb(r3(1)),1)

        geom_name   = 'neb_bead'//bead_index_name(i)//':geom'
        movecs_name = 'neb_bead'//bead_index_name(i)//'.movecs'
        geomlen     = inp_strlen(geom_name)
        movecslen   = inp_strlen(movecs_name)
        value = value.and.geom_cart_coords_set(geom,dbl_mb(r3(1)))
        value = value.and.geom_rtdb_store(rtdb,geom,
     >                                    geom_name(1:geomlen))

        call add_bead_list(bead_list,
     >                     movecs_name(1:movecslen),
     >                     geom_name(1:geomlen))

      end do
      value = value.and.geom_destroy(geom)
      value = value.and.MA_pop_stack(r3(2))
      value = value.and.MA_pop_stack(r2(2))
      value = value.and.MA_pop_stack(r1(2))
      if (.not.value) call errquit('neb_initial_path failed',3,0)


      return
      end

      subroutine neb_verlet_update(ng,c0,c1,v1,dti,g1)
        integer ng
        double precision c0(*)
        double precision c1(*)
        double precision v1(*)
        double precision dti(*)
        double precision g1(*)

        integer i

*        *** c1 <- 2*c1 - c0 + dti*g          ***
*        ***    <- c1 + ((c1-c0)/t)*t + dti*g ***

*        *** RRR loop over ions *** 
*        *** NEED: CONSTRAINED DYNAMICS (FROZEN ATOMS)
         do i=1,ng
           v1(i) = c1(i)-c0(i)
           if (v1(i)*g1(i).lt.0.0d0) v1(i) = 0.0d0
         end do
         call dcopy(ng,c1(1),1,c0(1),1)
         do i=1,ng
            c1(i) = 
     >              c1(i)
     >            + v1(i)  
     >            + dti(i)*g1(i)
         end do

      end
      subroutine neb_cg_direction(ng,g0,g1,s)
        integer ng
        double precision g0(*)
        double precision g1(*)
        double precision s(*)

        integer i
        double precision gamma1
        double precision gamma2
        double precision sn

c       *** choosing Polac-Ribiere coeff ***        
        gamma1=0.0d0
        gamma2=0.0d0
        do i=1,ng
           gamma1 = gamma1 + (g1(i)-g0(i))*g1(i)
           gamma2 = gamma2 + g0(i)*g0(i)
        end do

        do i=1,ng
           s(i) = -g1(i) + s(i)*gamma1/gamma2
        end do
        
      end

      subroutine neb_move(ng,dt,c0,c1,s)
        integer ng
        double precision dt
        double precision c0(*)
        double precision c1(*)
        double precision s(*)

        integer i


        do i=1,ng
            c1(i) = c0(i)+dt*s(i)
        end do

      end

      subroutine neb_calc_convergence(ng,g1,c0,c1,Gmax,Grms,Xmax,Xrms)
      implicit none
      integer ng
      double precision g1(*)
      double precision c0(*)
      double precision c1(*)
      double precision Gmax,Grms,Xmax,Xrms

      !***** local variables ****
      integer ii
      double precision dx

      Gmax = 0.0
      do ii=1,ng
        if (dabs(g1(ii)).gt.Gmax) Gmax = dabs(g1(ii))
      end do

      Grms = 0.0
      do ii=1,ng
        Grms = Grms + g1(ii)*g1(ii)
      end do
      Grms = dsqrt(Grms)/dble(ng)

      Xmax = 0.0
      do ii=1,ng
        dx = dabs(c1(ii)-c0(ii))
        if (dabs(g1(ii)).gt.Xmax) Xmax = dx
      end do

      Xrms = 0.0
      do ii=1,ng
        dx = dabs(c1(ii)-c0(ii))
        Xrms = Xrms + dx*dx
      end do
      Xrms = dsqrt(Xrms)/dble(ng)


      return
      end
