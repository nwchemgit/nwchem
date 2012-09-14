C> \ingroup task
C> @{
!
! Zero temperature string method (ZTS) to compute minimum energy paths (MEPs)
! The method works faster if the initial end points are both minima
! This code allows one to input a guess at the "middle" point also
! See:
! http://cims.nyu.edu/~eve2/string.htm (existed in 2011)
! W. E, W. Ren, and E. Vanden-Eijnden, Phys. Rev. B: 66, 052301 (2002).
! W. E, W. Ren, and E. Vanden-Eijnden, J. Chem. Phys.: 126, 164103 (2007).
! http://www.math.princeton.edu/string/index.html (existed in 2011)
!
!   It should be noted that the end points will go downhill because of 
!   this evolution in "time".  The non-end points also slip downhill, but 
!   the reparameterization forces them to still form a curve between
!   the two end points.  With that in mind, they should eventually
!   form a minimum energy path.
!

! The RTDB parameters that control the code are:
! string:tol   : The convergence criteria (based on RMS coordinate change)
! string:stepsize  : The "Time step" small=slow, big=unstable
! string:nbeads  : The number of points on the string more=>slow & stable
! string:maxit     : The number of time steps
! string:interpol  : Linear or spline or other.  Not very important.
! string:freeze1  : freeze endpoint 1
! string:freezen  : freeze endpoint N
! string:linopt   : Turn on/off linear optimization
! string:projection1  : Turn off/on projection of gradient
! string:pathguess  : which method for guessing path from 1 to TS to N


! This code requires at two or three geometry blocks to be defined
! "geometry"  Endpoint #1
! "endgeom"   Endpoint #2
! "midgeom"   Optional transition state
! At this point the code does not verify that they geometries are the same atoms
! in the same order, but they need to be for the code to work right.
! The main routine can handle having all geometries specified, but the 
! driver does support reading them in yet.

      logical function task_string(rtdb)
      implicit none

! External stuff
#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "geom.fh"

! Input parameters that do not get changed
      integer rtdb   ! The runtime database

! Local variables for memory blocks
      integer l_dg, k_dg,l_energy, k_energy, l_dcoords, k_dcoords
      integer l_dVdx, k_dVdx,l_coords_old, k_coords_old
      integer l_wrk, k_wrk,l_coords, k_coords

! Local variables
      integer maxit,nbeads,natoms,middle,geom,geom2,size,i,interpol
      real*8 tol,stepsize
      logical converged, freeze1, freezeN

! Set convergence criteria
      if (.not. rtdb_get(rtdb,'string:tol',mt_dbl,1,tol))
     $      tol = 0.0001d0

! Set the time step size
      if (.not. rtdb_get(rtdb,'string:stepsize',mt_dbl,1,stepsize))
     $      stepsize = 0.1d0

! Number of points on the path
      if (.not. rtdb_get(rtdb,'string:nbeads',mt_int,1,nbeads))
     $      nbeads=10

! Number of steps
      if (.not. rtdb_get(rtdb,'string:maxit',mt_int,1,maxit))
     $      maxit = 25

! They type of interpolation
      if (.not. rtdb_get(rtdb,'string:interpol',mt_int,1,interpol))
     $      interpol = 1

! The general geometry structure
      if (.not. geom_create(geom, 'geometry'))
     &   call errquit('task_string: geom_create?', 911, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &   call errquit('task_string: no geometry ', geom, RTDB_ERR)

! Get number of atoms
      if (.not. geom_ncent(geom,natoms))
     $     call errquit('task_string: natoms?',geom, GEOM_ERR)

! Get memory
      if(.not.ma_push_get(mt_dbl,nbeads,'dg',l_dg, k_dg))
     $     call errquit('task_string: memory', nbeads, MA_ERR)
      if(.not.ma_push_get(mt_dbl,nbeads,'energy',l_energy, k_energy))
     $     call errquit('task_string: memory', nbeads, MA_ERR)
      if(.not.ma_push_get(mt_dbl,nbeads,'dcoord',l_dcoords, k_dcoords))
     $     call errquit('task_string: memory', nbeads, MA_ERR)
      size = 3*natoms*nbeads
      if(.not.ma_push_get(mt_dbl,size,'dVdx',l_dVdx, k_dVdx))
     $     call errquit('task_string: memory', size, MA_ERR)
      if(.not.ma_push_get(mt_dbl,size,'coords_old',l_coords_old,
     $     k_coords_old))
     $     call errquit('task_string: memory', size, MA_ERR)
      if(.not.ma_push_get(mt_dbl,size,'wrk',l_wrk, k_wrk))
     $     call errquit('task_string: memory', size, MA_ERR)
      if(.not.ma_push_get(mt_dbl,size,'coords',l_coords, k_coords))
     $     call errquit('task_string: memory', size, MA_ERR)

! Set the coordinates of geometry 1
      if (.not. geom_cart_coords_get(geom, dbl_mb(k_coords)))
     $        call errquit('task_string: geom?',0, GEOM_ERR)

! Set the coordinates of last geometry
      size = k_coords + (nbeads-1)*3*natoms
      if (.not. geom_create(geom2, 'geom2'))
     &   call errquit('task_string: geom_create?', 921, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom2, 'endgeom'))
     &   call errquit('task_string: no endgeom ', geom2, RTDB_ERR)
      if (.not. geom_cart_coords_get(geom2, dbl_mb(size)))
     $        call errquit('task_string: geom?',0, GEOM_ERR)

! Middle geometry -- if it exists
      if (.not. geom_rtdb_load(rtdb, geom2, 'midgeom')) then
         middle = 0
      else
        middle = (nbeads/2)+1
        size = k_coords + (middle-1)*3*natoms
        if (.not. geom_cart_coords_get(geom2, dbl_mb(size)))
     $        call errquit('task_string: geom?',0, GEOM_ERR)
      endif

      if (.not. geom_destroy(geom2))
     &   call errquit('task_string: geom_destroy?', geom2, GEOM_ERR)

! Do we freeze bead 1 or N?
      if (.not. rtdb_get(rtdb, 'string:freeze1', mt_log,1,freeze1))
     &        freeze1 = .false.
      if (.not. rtdb_get(rtdb, 'string:freezen', mt_log,1,freezeN))
     &        freezeN = .false.

! Warning: No way to read in all nbead geomtries

! Call the code
      if (.not. rtdb_put(rtdb,'opt:string',mt_log,1,.true.))
     $     call errquit('task_string: rtdb_put failed',0,
     &       RTDB_ERR)
      call zts_meps(maxit,nbeads,tol,stepsize,natoms,middle,
     &    dbl_mb(k_coords),geom,rtdb,dbl_mb(k_dg),dbl_mb(k_energy),
     &    dbl_mb(k_dVdx),dbl_mb(k_dcoords),dbl_mb(k_coords_old),
     &    dbl_mb(k_wrk),interpol,freeze1,freezeN,converged)
      if (.not. rtdb_put(rtdb,'opt:string',mt_log,1,.false.))
     $     call errquit('task_string: rtdb_put failed',0,
     &       RTDB_ERR)
! Return memory - Being lazy, just pop is all in one fell swoop
      if(.not. ma_chop_stack(l_dg))
     &   call errquit('task_string: ma?',99,MA_ERR)

      if (.not. geom_destroy(geom))
     &   call errquit('task_string: geom_destroy?', geom, GEOM_ERR)

! Done
      task_string = converged
      return
      end
C> @}

      subroutine zts_meps(maxit,nbeads,tol,stepsize,natoms,middle,
     &    coords,geom,rtdb,dg,energy,dVdx,dcoords,coords_old,wrk,
     &    interpol,freeze1,freezeN,converged)
      implicit none
      
! External stuff
      logical task_gradient, task_energy
      external task_gradient, task_energy
      real*8 zts_distance
      external zts_distance
      integer ga_nodeid
      external ga_nodeid
#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "geom.fh"
#include "util.fh"
#include "inp.fh"

! Input parameters that do not get changed
      integer maxit  ! Maximum number of iterations (200?)
      integer nbeads  ! Number of images or beads along the path (3 or more)
      real*8 tol ! Tolerance for convergence  (1.0d-7)
      real*8 stepsize ! The "time" step size (0.0001)
      integer natoms ! The number of atoms
      integer middle ! The "middle" coordinates bead number - 0 if none, -1 all
                     !      of the coordinates are filled
      integer rtdb   ! The runtime database
      integer interpol ! The style of interpolation
      logical freeze1, freezeN ! Freeze ends of string

! Input/output
      real*8  coords(3,natoms,nbeads)
      integer geom   ! A geometry record, should match molecule of interest
      logical converged  ! Did we converge
      
! Variables used locally, but not passed in or out
      real*8 dg(nbeads) ! 0,.5,1 for three, 0,1/3,2/3,1 for four
      real*8 energy(nbeads)
      real*8 dVdx(3,natoms,nbeads)
      real*8 dcoords(nbeads)
      real*8 coords_old(3,natoms,nbeads)
      real*8 wrk(3,natoms,nbeads)

! Local variables
      integer i, j, k, nstep
      real*8 tmp, tmp2,t
      logical badgeom
      real*8 g1(3),g2(3),g3(3)
      integer iistart, iiend
      real*8 hess, e0, e1, a2, alpha, dsgrad
      logical linopt
      logical projection1,found,status
      integer pathguess

      character     ch_tmp
      character*255 filename,full_filename,perm_name,movecs_name
      integer movecslen,permlen

      character*7 bead_index_name
      external    bead_index_name

! Setup the problem
! "Secret" options with defaults that generally do not need changed
      if (.not. rtdb_get(rtdb, 'string:linopt', mt_log,1,linopt))
     &        linopt = .true.
      if (.not.rtdb_get(rtdb,'string:projection1',mt_log,1,projection1))
     &        projection1 = .false.
      if (.not. rtdb_get(rtdb,'string:pathguess',mt_int,1,pathguess))
     $      pathguess = 2

      do i = 1, nbeads
        dg(i) = dble(i-1)/dble(nbeads-1)
      enddo
      
      if (middle .eq. 0) then  ! Flexible linear interpolation to get TS
        call zts_guessall(natoms,nbeads,coords,geom)
      else if (middle .gt. 1 .and. middle .lt. nbeads) then ! TS guess
        call zts_guess(natoms,nbeads,middle,coords,pathguess)
      else if (middle .eq. -1) then  ! We know exactly what we are doing
        ! Do nothing.  All coordinates were setup on entry.
      else  ! This makes no sense
        call errquit('zts_meps: fatal error', 0, 0)
      endif

*     **** set .string_epath filename ****s
      call util_file_prefix('string_epath',filename)
      call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                             full_filename)

! Now we do the work
      if (ga_nodeid() .eq. 0) then
        write(*,'(a)')'@ String method.'
        write(*,'(a,f9.5)')'@ Temperature          = ', 0.0d0
        write(*,'(a,f9.5)')'@ Covergence Tolerance = ', tol
        write(*,'(a,f9.5)')'@ Step Size            = ', stepsize
        write(*,'(a,i9)')  '@ Maximum Time Steps   = ', maxit
        write(*,'(a,i9)')  '@ Number of replicas   = ', nbeads
        write(*,'(a,i9)')  '@ String Interpolator  = ', interpol
        if (freeze1) then
          write(*,'(a)')   '@ First Replica        = frozen'
        else
          write(*,'(a)')   '@ First Replica        = moves'
        endif
        if (freezeN) then
          write(*,'(a)')   '@ Last Replica         = frozen'
        else
          write(*,'(a)')   '@ Last Replica         = moves'
        endif
        write(*,'(a)')
     &    '@  Step   xrms          E end         '
     &     //'E middle          E end          E max'
      endif
      if (maxit .lt. 1) then ! Just printing out initial geometries
         do i = 1, nbeads
           energy(i) = 0
         enddo
      endif
      do nstep = 1, maxit
!        Look for atoms that are too close
         do i = 1, nbeads
          call zts_Robinson(natoms,geom,coords(1,1,i),badgeom)
         enddo
!        Project out net motion
         do i = 2, nbeads
          call zts_min_motion(natoms,coords(1,1,i),coords(1,1,i-1),geom)
         enddo
         ! Save old coordinates to check convergence
         do i = 1, nbeads
           do j = 1, natoms
             do k = 1, 3
               coords_old(k,j,i) = coords(k,j,i)
             enddo
           enddo
         enddo
         ! Calculated the energy and forces for each bead
         do i = 1, nbeads
           if ( .not. geom_cart_coords_set( geom, coords(1,1,i) ) )
     &        call errquit('zts_meps: set geometry',0,GEOM_ERR)
           call sym_geom_project(geom, 1.0d-6)
           if (.not. geom_rtdb_store(rtdb, geom, 'geometry'))
     &       call errquit('zts_meps: store geometry',0,GEOM_ERR)
           if ((i .eq. 1 .and. freeze1 .and. nstep .ne. 1) .or. 
     &         (i .eq. nbeads .and. freezeN .and. nstep .ne. 1)) then
              ! Skip the gradient for frozen stuff, unless first step - need E
           else

*            *** copy bead movecs to taskmovecs ***
             call util_file_prefix('movecs',perm_name)
             movecs_name = "w"//bead_index_name(i)//perm_name
             call util_file_name_resolve(perm_name, .false.)
             call util_file_name_resolve(movecs_name, .false.)
             permlen = inp_strlen(perm_name)
             movecslen = inp_strlen(movecs_name)
             if (ga_nodeid().eq.0) then
               inquire(file=movecs_name,exist=status)
               if (status) then
                  call util_file_copy(movecs_name(1:movecslen),
     >                                perm_name(1:permlen))
               end if
             end if

             call ga_sync()
             if ( .not. task_gradient(rtdb))
     &          call errquit('zts_meps: gradient failed',0, 0)

*            *** copy taskmovecs to bead movecs ****
             if (ga_nodeid().eq.0) then
                inquire(file=perm_name,exist=status)
                if (status) then
                   call util_file_copy(perm_name(1:permlen),
     >                                 movecs_name(1:movecslen))
                end if
             end if

             if (.not. rtdb_get(rtdb,'task:energy',mt_dbl,1,energy(i)))
     &          call errquit('zts_meps: failed getting energy',0,
     &          RTDB_ERR)
             if (.not. rtdb_get(rtdb, 'task:gradient', mt_dbl, 3*natoms,
     &           dVdx(1,1,i))) call errquit(
     &           'zts_meps: failed getting gradient',0,RTDB_ERR)
           endif
         enddo
         if (freeze1) then
           do j = 1, natoms
             do k = 1, 3
               dVdx(k,j,1) = 0.0d0
             enddo
           enddo
         endif
         if (freezeN) then 
           do j = 1, natoms
             do k = 1, 3
               dVdx(k,j,nbeads) = 0.0d0 
             enddo
           enddo
         endif
         ! Projection of gradient onto perpedicular to path?
         if (projection1) then
           do i = 2, nbeads - 1
             do j = 1, natoms
               do k = 1, 3
                 g1(k) = coords(k,j,i) - coords(k,j,i+1)
                 g2(k) = coords(k,j,i) - coords(k,j,i-1)
               enddo
               call cross_product(g1,g2,g3) ! Vector perpendicular to plane
               tmp = 0.0d0
               do k = 1, 3
                 tmp = tmp + g3(k)**2
               enddo
               if (tmp .gt. 1.0d-16) then
                 tmp = 1.0d0/sqrt(tmp)
                 do k = 1, 3
                   g3(k) = g3(k)*tmp ! Normalize
                 enddo
                 ! Projection
                 call cross_product(dVdx(1,j,i),g3,g2)
                 call cross_product(g3,g2,dVdx(1,j,i))
               endif
             enddo
           enddo
         endif
         ! Linear extrapolation - based upon linopt=1 in driver
         if (linopt) then
           if (freeze1) then
              iistart = 2
           else
              iistart = 1
           endif
           if (freezeN) then
              iiend   = nbeads-1
           else
              iiend   = nbeads
           endif
           do i = iistart,iiend
             alpha = 1.0d0
             dsgrad = 0.0d0
             do j = 1, natoms
               do k = 1, 3
                  dsgrad = dVdx(k,j,i)**2
               enddo
             enddo
             dsgrad=-stepsize*dsgrad ! ds*g (with ds=-stepsize*g)
             e0 = energy(i)
c
c     Compute energy at alpha*ds
c
             do j = 1, natoms
               do k = 1, 3
                 wrk(k,j,1)=coords(k,j,i)-alpha*stepsize*dVdx(k,j,i) ! Force=-Gradient
               enddo
             enddo
             if ( .not. geom_cart_coords_set( geom, wrk) )
     &          call errquit('zts_meps: set geometry',0,GEOM_ERR)
             call sym_geom_project(geom, 1.0d-6)
             if (.not. geom_rtdb_store(rtdb, geom, 'geometry'))
     &         call errquit('zts_meps: store geometry',0,GEOM_ERR)

*            *** copy bead movecs to taskmovecs ***
             call util_file_prefix('movecs',perm_name)
             movecs_name = "w"//bead_index_name(i)//perm_name
             call util_file_name_resolve(perm_name, .false.)
             call util_file_name_resolve(movecs_name, .false.)
             permlen = inp_strlen(perm_name)
             movecslen = inp_strlen(movecs_name)
             if (ga_nodeid().eq.0) then
               inquire(file=movecs_name,exist=status)
               if (status) then
                  call util_file_copy(movecs_name(1:movecslen),
     >                                perm_name(1:permlen))
               end if
             end if

             call ga_sync()
             if ( .not. task_energy(rtdb))
     &          call errquit('zts_meps: energy failed',0, 0)

*            *** copy taskmovecs to bead movecs ****
             if (ga_nodeid().eq.0) then
                inquire(file=perm_name,exist=status)
                if (status) then
                   call util_file_copy(perm_name(1:permlen),
     >                                 movecs_name(1:movecslen))
                end if
             end if

             if (.not. rtdb_get(rtdb,'task:energy',mt_dbl,1,e1))
     &          call errquit('zts_meps: failed get energy',0,RTDB_ERR)
c
c            Fit to a parabola as a function of alpha
c
             hess = (e1-e0-alpha*dsgrad)/alpha**2
c
             if (abs(e1-e0) .lt. 1.0d-7) then
                a2 = alpha ! If no big deal, be conservative
             else if (hess .gt. 0d0) then
c
c            Positive curvature.
c
                a2 = -dsgrad/(2d0*hess)  ! a2 will be positive
                if (a2 .gt. 4.0d0*alpha) then
                   a2 = 4.0d0*alpha! Limit
                endif
             else
c
c            Negative curvature and downhill
c
                 a2 = 4.0d0*alpha
             endif
c            scale gradient
             alpha = a2
             do j = 1, natoms
               do k = 1, 3
                 dVdx(k,j,i) = dVdx(k,j,i) * alpha
               enddo
             enddo 
           enddo
         endif
         ! Take a step (evolve it in "time").  Warning: Very simple minded
         do i = 1, nbeads
           do j = 1, natoms
             do k = 1, 3
               coords(k,j,i) = coords(k,j,i) - stepsize*dVdx(k,j,i) ! Force=-Gradient
             enddo
           enddo
         enddo
!        Project out net motion - Can rotate/translater frozen N bead
         do i = 2, nbeads
          call zts_min_motion(natoms,coords(1,1,i),coords(1,1,i-1),geom)
         enddo
         !reparameterize - This is where the magic happens
         dcoords(1) = 0.0d0
         do i = 2, nbeads
           dcoords(i) = 
     &      zts_distance(natoms,coords(1,1,i),coords(1,1,i-1))
         enddo
         do i = 2, nbeads
           dcoords(i) = dcoords(i) + dcoords(i-1) ! Now relative to bead 1
         enddo
         do i = 1, nbeads
           dcoords(i) = dcoords(i)/dcoords(nbeads) ! Runs from zero to one
         enddo
         do i = 2, nbeads-1
           if (dcoords(i) .le. dcoords(i-1)) then  ! paranoid
             dcoords(i) = (99.0d0*dcoords(i-1)+dcoords(i+1))/100.0d0
           endif
         enddo
         if (dcoords(nbeads) .le. dcoords(nbeads-1)) then  ! paranoid
          dcoords(nbeads-1) = (dcoords(nbeads-2)+dcoords(nbeads))/2.0d0
         endif
         ! We now interpolate new points along the string
         ! This adds more points near regions of big dcoords
         ! This really a normalized arc length on a spline
         if (interpol .eq. 1) then
           call zts_linear_fit(natoms*3,nbeads,coords,dcoords,wrk)
         else if (interpol .eq. 2) then
           call zts_spline_fit(natoms*3,nbeads,coords,dcoords,dg)
         else if (interpol .eq. 3) then
           call zts_akima_fit(natoms*3,nbeads,coords,dcoords,dg)
         else
           call errquit('zts_meps: fatal error', 0, 0)
         end if
!        Project out net motion - Can rotate/translater frozen N bead
         do i = 2, nbeads
          call zts_min_motion(natoms,coords(1,1,i),coords(1,1,i-1),geom)
         enddo
         ! Check convergence
         tmp = 0.0d0
         do i = 1, nbeads
           tmp = tmp +
     &      zts_distance(natoms,coords(1,1,i),coords_old(1,1,i))
         enddo
         tmp = tmp/dble(nbeads*natoms)
         if (ga_nodeid() .eq. 0) then
           tmp2 = energy(1)
           do i = 2, nbeads
             if (energy(i) .gt. tmp2) tmp2=energy(i)
           enddo
           write(*,'(a2,i5,1x,f8.6,1x,f14.7,1x,f14.7,1x,f14.7,1x,f14.7)'
     &            )'@ ',nstep,tmp,energy(1),
     &                            energy(1+nbeads/2),
     &                            energy(nbeads),tmp2

*           ***************************************
*           ***** print out .string_epath file ****
*           ***************************************
            inquire(file=full_filename,exist=found)
*           **** FILE already exists - parse to EOF ****
            if (found) then
              open(unit=19,file=full_filename,form='formatted',
     >             status='old')
              do while(.true.)
                read(19,*,ERR=30,END=30) ch_tmp
              end do
 30           continue
#if defined(FUJITSU_SOLARIS) || defined(PSCALE) || defined(SOLARIS) || defined(__crayx1) || defined(GCC46)
              backspace 19
#endif
            write(19,*) " "
*           **** .string_epath FILE does not exist ****
            else
              open(unit=19,file=full_filename,form='formatted')
            end if
            write(19,*)
     > "#-------------------------------------------------------"
            write(19,*)  "# String Path iteration  = ",nstep
            write(19,*)  "# String Interpolator    = ",interpol
            write(19,*)  "# nbeads                 = ",nbeads
            write(19,*)  "# step Size              = ",stepsize
            write(19,*)  "# convergence            = ",tmp
            write(19,*)  "# convergence tol (Xrms) = ",tol
            write(19,*)
     > "#-------------------------------------------------------"
            do i=1,nbeads
               t = (i-1)/dble(nbeads-1)
               write(19,*) t,energy(i)
            end do
            close(19)

         endif
         if ( tmp .lt. tol) goto 747
      enddo
 747  continue
      if (nstep .ge. maxit) then
         if (ga_nodeid() .eq. 0)
     &   write(*,'(a)') '@ The string calculation failed to converge'
         converged = .false.
      else
         if(ga_nodeid() .eq. 0)
     &   write(*,'(a)') '@ The string calculation converged'
         converged = .true.
      endif


*     **** write out xyz trajectory ****
      call util_file_prefix('string_final.xyz',filename)
      call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                             full_filename)
      if (ga_nodeid() .eq. 0) then
      open(unit=19,file=full_filename,form='formatted')

       write(*,*)
       do i = 1, nbeads
         write(*,'(a16,i5,a21,f20.12)') '@ Bead number = ', i ,
     &       ' Potential Energy = ', energy(i)
         if ( .not. geom_cart_coords_set( geom, coords(1,1,i) ) )
     &        call errquit('zts_meps: set geometry',i,GEOM_ERR)
         if ( .not. geom_print_xyz(geom,6))
     &        call errquit('zts_meps: print geometry',i,GEOM_ERR)
         if ( .not. geom_print_xyz(geom,19))
     &        call errquit('zts_meps: print geometry',i,GEOM_ERR)
       enddo
       close(19)
       end if


*     **** write out final path energies ****
      call util_file_prefix('string_final_epath',filename)
      call util_file_name_noprefix(filename,.false.,
     >                             .false.,
     >                             full_filename)
      if (ga_nodeid() .eq. 0) then
         open(unit=19,file=full_filename,form='formatted')
         write(19,*)
     > "#-------------------------------------------------------"
         write(19,*)  "# String Path iteration  = ",nstep
         write(19,*)  "# String Interpolator    = ",interpol
         write(19,*)  "# nbeads                 = ",nbeads
         write(19,*)  "# step Size              = ",stepsize
         write(19,*)  "# convergence            = ",tmp
         write(19,*)  "# convergence tol (Xrms) = ",tol
         write(19,*)
     > "#-------------------------------------------------------"
         do i=1,nbeads
            t = (i-1)/dble(nbeads-1)
            write(19,*) t,energy(i)
         end do
         close(19)
      endif

      if (ga_nodeid() .eq. 0) write(*,*)
      
      end
      
      subroutine zts_linear_fit(nvar,nbeads,Y,Xin,wrk)
      implicit none
      ! This will return Y evenly spaced in X, instead of spaced Xin
      ! This uses linear interpolation.

      ! input
      integer nvar
      integer nbeads
      real*8 Xin(nbeads)
      ! input/output
      real*8 Y(nvar,nbeads)
      ! Work storage
      real*8 wrk(nvar,nbeads)
      ! local variables
      integer ivar, ibead, point
      real*8 X, XinRight1, XinLeft1, frac

      point = 1
      do ibead = 2, nbeads-1
        ! Find where this point is in the data
        X = dble(ibead-1)/dble(nbeads-1)  ! runs from zero to one
 11     if (Xin(point+1) .le. X) then
          point = point + 1
          goto 11
        endif
        if (point .eq. nbeads) point = nbeads-1 ! paranoid
        ! Solve linear interpolation function for data
        XinLeft1 = Xin(point)
        XinRight1 = Xin(point+1)
        frac = (X-XinLeft1)/(XinRight1-XinLeft1)
        do ivar = 1, nvar
          wrk(ivar,ibead) = (1.0d0-frac)*Y(ivar,point)+
     &                      frac*Y(ivar,point+1)
        enddo 
      enddo
      do ibead = 2, nbeads-1 ! end points do not change
        do ivar = 1, nvar
          Y(ivar,ibead) = wrk(ivar,ibead)
        enddo
      enddo
      
      return
      end
      
      subroutine zts_spline_fit(nvar,nbeads,Y,Xin,dg)
      implicit none
      ! This will return Y evenly spaced in X, instead of spaced Xin
      ! This uses spline interpolation.

      ! input
      integer nvar
      integer nbeads
      real*8 Xin(nbeads)
      real*8 dg(nbeads)
      ! input/output
      real*8 Y(nvar,nbeads)
      ! local variables
      integer ivar, ibead
      real*8 X
      
      integer MAX_N
      PARAMETER (MAX_N=750)
      real*8 Y_tmp(MAX_N),fit(MAX_N)

      if (nbeads .gt. MAX_N)
     &     call errquit('zts_meps: MAX_N exceeded', 0, 0)

      do ivar = 1, nvar  ! loop over all the coordinates
        do ibead = 1, nbeads
          Y_tmp(ibead) = Y(ivar,ibead)
        enddo
        call zts_spline_Setup(Xin,Y_tmp,nbeads,fit) ! (x,y,n,fit)

        do ibead = 2, nbeads-1
          X = dg(ibead)
          call zts_spline(Xin,Y_tmp,nbeads,fit,X,Y(ivar,ibead)) ! (x,y,n,fit,xin,yout)
        enddo
      enddo
      return
      end
      
      subroutine zts_spline_Setup(x,y,n,fit)
      ! Input
      integer n
      real*8 x(n),y(n)
      ! Output
      real*8 fit(n)
      ! Local
      integer MAX_N
      PARAMETER (MAX_N=750)
      integer i
      real*8 tmp1,tmp2,u(MAX_N)
! This routine solves the tri-diagonal problem of a natural spline
      fit(1)=0.0d0
      u(1)=0.0d0
      
      if (n .gt. MAX_N)
     &     call errquit('zts_meps: MAX_N exceeded', 0, 0)

      do i=2,n-1
        tmp1=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        tmp2=tmp1*fit(i-1)+2
        fit(i)=(tmp1-1.0d0)/tmp2
        u(i)=(6.0d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/
     &        (x(i)-x(i-1)))/(x(i+1)-x(i-1))-tmp1*u(i-1))/tmp2
      enddo

      fit(n)=0.0d0
      do i=n-1,1,-1
        fit(i)=fit(i)*fit(i+1)+u(i)
      enddo
      return
      end

      subroutine zts_spline(x,y,n,fit,xin,yout)
      integer n
      real*8 xin,yout,x(n),fit(n),y(n)
      integer try,high,low
      real*8 a,b,h
      low=1
      high=n
      do while (high-low .gt. 1)
        try=(high+low)/2
        if(x(try).gt.xin)then
          high=try
        else
          low=try
        endif
      enddo
      h=x(high)-x(low)
      if (h.eq.0.0d0) then
        yout=(y(low)+y(high))/2.0d0 ! That is not good
      else
        a=(x(high)-xin)/h
        b=(xin-x(low))/h
        yout=a*y(low)+b*y(high)+((a**3-a)*fit(low)+(b**3-b)*fit(high))*
     &           (h**2)/6.0d0
      endif
      return
      end

      subroutine zts_min_motion(natoms,new ,old, geom)
      implicit none
#include "errquit.fh"
      integer natoms, geom
      real*8 new(3,natoms), old(3,natoms)

      real*8 improve
      integer isystype
      logical geom_systype_get
      external geom_systype_get

      if (.not. geom_systype_get(geom, isystype))
     $     call errquit('zts_min_motion: isystype?',0, GEOM_ERR)
      
      improve = 1.0d10
      do while (improve .gt. 0.001d0)
         improve = 0.0d0
         call zts_xyz(1,natoms,new,old,improve)
         call zts_xyz(2,natoms,new,old,improve)
         call zts_xyz(3,natoms,new,old,improve)
         if (isystype .eq. 0) then
           call zts_xyz(4,natoms,new,old,improve)
           call zts_xyz(5,natoms,new,old,improve)
           call zts_xyz(6,natoms,new,old,improve)
         else if (isystype .eq. 1) then
           call zts_xyz(1,natoms,new,old,improve) ! Does not change X
         endif
      enddo

      return
      end

      subroutine zts_xyz(which,natoms,new,old,improve)
      implicit none
      real*8 zts_distance
      external zts_distance
      integer which
      integer natoms
      real*8 new(3,natoms), old(3,natoms)
      real*8 improve
      
      real*8 tmp(3,natoms)
      real*8  left,  right,  try, middle
      real*8 xleft, xright, xtry
      real*8 initial, final

      initial = zts_distance(natoms,old,new)
      xleft   = -0.1d0
      xright  =  0.1d0
      middle  = initial

 99   continue  ! Walk unto we bracket the minimum
      call zts_copy_coords(natoms,new,tmp,which,xleft)
      left = zts_distance(natoms,old,tmp)
      if (left .lt. middle) then
         middle = left
         call zts_copy_coords(natoms,tmp,new,0,0.0d0)
         goto 99
      endif
      call zts_copy_coords(natoms,new,tmp,which,xright)
      right = zts_distance(natoms,old,tmp)
      if (right .lt. middle) then
         middle = right
         call zts_copy_coords(natoms,tmp,new,0,0.0d0)
         goto 99
      endif

! Now we have a middle that is better than left and right - go hunting
! xmiddle is always zero, xleft is negative, xright is positive
! We don't try to be smart with this linear search, since
! it is fast and non-linear

      do while ( (xright - xleft) .gt. 0.001d0)
        xtry = xright*0.5d0
        call zts_copy_coords(natoms,new,tmp,which,xtry)
        try = zts_distance(natoms,old,tmp)
        if (try .ge. middle) then
          right  = try
          xright = xtry
        else
          xleft  = -xtry
          left   = middle
          xright = xright-xtry
          right  = right
          middle = try
          call zts_copy_coords(natoms,tmp,new,0,0.0d0)
        endif
        xtry = xleft*0.5d0
        call zts_copy_coords(natoms,new,tmp,which,xtry)
        try = zts_distance(natoms,old,tmp)
        if (try .ge. middle) then
          left  = try
          xleft = xtry
        else
          xright  = -xtry
          right   = middle
          xleft   = xleft-xtry
          left    = left
          middle  = try
          call zts_copy_coords(natoms,tmp,new,0,0.0d0)
        endif
      enddo

      final = zts_distance(natoms,old,new)
      improve = improve + (initial - final)
      return
      end

! Make this a function so that we can easily mass weight it
! or otherwise play with it
      real*8 function zts_distance(n,c1,c2)
      implicit none
#include "errquit.fh"
      integer n
      real*8 c1(3,n),c2(3,n)
      integer j
      real*8 tmp

! Warning:  Does not deal with periodic wrapping yet

      tmp = 0.0d0
      do j = 1, n
           tmp = tmp + sqrt(
     &          (c1(1,j) - c2(1,j))**2 +
     &          (c1(2,j) - c2(2,j))**2 +
     &          (c1(3,j) - c2(3,j))**2)
      enddo
      zts_distance = tmp
      return
      end
      
! By abstracting this out, we can treat angles and distances the same above

      subroutine zts_copy_coords(n,in,out,which,shift)
      implicit none
      integer n
      real*8 in(3,n),out(3,n)
      integer which
      real*8 shift
      
      integer i,j,k
      real*8 tmp1, tmp2
      
      
      do i = 1, n
          out(1,i) = in(1,i)
          out(2,i) = in(2,i)
          out(3,i) = in(3,i)
      enddo
      if (which .eq. 0 .or. shift .eq. 0.0d0) then
!       Nothing to do
      else if (which .ge. 1 .and. which .le. 3) then
        do i = 1, n
          out(which,i) = out(which,i) + shift
        enddo
      else if (which .ge. 4 .and. which .le. 6) then ! Angles
        j = which - 3
        k = which - 2
        if (k .eq. 4) k = 1
        do i = 1, n
           tmp1 = out(j,i)*cos(shift) - out(k,i)*sin(shift)
           tmp2 = out(j,i)*sin(shift) + out(k,i)*cos(shift)
           out(j,i) = tmp1
           out(k,i) = tmp2
        enddo
      else
        call errquit('zts_copy_coords: bad which',0,0)
      endif
      return
      end

! If atoms are too close, then fix it
! Subroutine is named after the phrase "Danger, Will Robinson"
! Warning, no support for periodic images yet
       subroutine zts_Robinson(natoms,geom,coords,stillbad)
       implicit none

! Global includes
#include "nwc_const.fh"
       logical geom_cart_get_charges
       external geom_cart_get_charges
! Input only
       integer natoms
       integer geom
! Output only
       logical stillbad
! Input/Output
       real*8 coords(3,natoms)
! Local variables
       integer max_bad
       parameter (max_bad = 250)
       real*8 bad(4,max_bad)
       real*8 typei, typej
       integer is_h
       real*8 min_h_h, min_x_h, min_x_x, min_dist
       integer trys
       integer num_bad
       real*8 siner, coser, cosed, sined 
       real*8 r
       integer i,j,try,ij
       real*8 charges(nw_max_atom)

! Warning: These assume Angstom units and non-fractional
! Otherwise they are way too hard a constraint
       min_h_h = 0.7d0
       min_x_h = 0.9d0
       min_x_x = 1.0d0

       trys = 10

       if (.not. geom_cart_get_charges(geom,natoms,charges))
     1        call errquit('zts_Robinson:charges',i, 0)

       do try = 1, trys
        stillbad=.false.
       
        num_bad=0

!       Find the bad ones

        do i=1,natoms
         typei = charges(i)
         do j=1,i-1
           typej = charges(j)
           r = sqrt((coords(1,j) - coords(1,i))**2
     1          +   (coords(2,j) - coords(2,i))**2
     2          +   (coords(3,j) - coords(3,i))**2)
           is_h = 0
           if (typei .le. 2) is_h = is_h + 1
           if (typej .le. 2) is_h = is_h + 1

           if (((r .lt. min_h_h) .and. (is_h .eq. 2)) .or.
     1         ((r .lt. min_x_h) .and. (is_h .eq. 1)) .or.
     2         ((r .lt. min_x_x) .and. (is_h .eq. 0))) then
               num_bad=num_bad+1
               if(num_bad.gt.max_bad) then
                  call errquit('zts_Robinson:too many bad',num_bad, 0)
               endif
               bad(1,num_bad)=i
               bad(2,num_bad)=j
               bad(3,num_bad)=typei
               bad(4,num_bad)=typej
           endif
         enddo
        enddo

!       Fix the bad ones

        do ij=1,num_bad
          stillbad = .true.
          typei = bad(3,ij)
          typej = bad(4,ij)
          if (typej .lt. typei) then  ! Move lighter atom
            i = nint(bad(2,ij))
            j = nint(bad(1,ij))
          else
            i = nint(bad(1,ij))
            j = nint(bad(2,ij))
          endif
          is_h = 0
          if (typei .le. 2) is_h = is_h + 1
          if (typej .le. 2) is_h = is_h + 1

          if (is_h .eq. 2) then
            min_dist = min_h_h
          else if (is_h .eq. 1) then
            min_dist = min_x_h
          else
            min_dist = min_x_x
          endif
          if      (try .le. 2) then
             min_dist = min_dist / 3.0d0
          else if (try. le. 4) then
             min_dist = min_dist / 2.0d0
          else if (try .le. 5) then
             min_dist = min_dist / 1.5d0
          endif
          ! What we have now
          r = sqrt((coords(1,i)-coords(1,j))**2 + 
     1             (coords(2,i)-coords(2,j))**2 +
     2             (coords(3,i)-coords(3,j))**2)
          if (r .lt. 0.0001d0) then
           ! different amounts in each direction to make sure symmetry is broken
           if (coords(3,i) .gt. coords(3,j)) then
             coords(3,i) = coords(3,i) + 0.001d0
           else
             coords(3,i) = coords(3,i) - 0.001d0
           endif
           if (coords(2,i) .gt. coords(2,j)) then
             coords(2,i) = coords(2,i) + 0.0001d0
           else
             coords(2,i) = coords(2,i) - 0.0001d0
           endif
           if (coords(1,i) .gt. coords(1,j)) then
             coords(1,i) = coords(1,i) + 0.00001d0
           else
             coords(1,i) = coords(1,i) - 0.00001d0
           endif
          else if (r .lt. min_dist) then
           siner = (coords(3,i)-coords(3,j))/r
           coser = sqrt(1.0d0-siner*siner)
           if(coser .lt. 1.0d-10) then ! Force 0*0/0 to be zero
             cosed = 0.0d0
             sined = 0.0d0
           else
             cosed = (coords(1,i)-coords(1,j))/(r*coser)
             sined = (coords(2,i)-coords(2,j))/(r*coser)
           endif
           ! What we get instead
           coords(3,i) = min_dist*siner + coords(3,j)
           coords(2,i) = (min_dist * coser) * sined + coords(2,j)
           coords(1,i) = (min_dist * coser) * cosed + coords(1,j)
          endif
        enddo
        if (.not. stillbad) goto 395
       enddo ! "try" loop

 395   continue

       return
       end


      subroutine zts_akima_fit(nvar,nbeads,Y,Xin,dg)
      implicit none
      ! This will return Y evenly spaced in X, instead of spaced Xin
      ! This uses Akima spline interpolation.

      ! input
      integer nvar
      integer nbeads
      real*8 Xin(nbeads)
      real*8 dg(nbeads)
      ! input/output
      real*8 Y(nvar,nbeads)
      ! local variables
      integer ivar, ibead
      real*8 X
      
      integer MAX_N
      PARAMETER (MAX_N=750)
      real*8 Y_tmp(MAX_N)

      if (nbeads .gt. MAX_N)
     &     call errquit('zts_meps: MAX_N exceeded', 0, 0)
! WARNING: Would be faster, if we did a "setup" like natural spline code
      do ivar = 1, nvar  ! loop over all the coordinates
        do ibead = 1, nbeads
          Y_tmp(ibead) = Y(ivar,ibead)
        enddo
        do ibead = 2, nbeads-1
          X = dg(ibead)
          call zts_akima_spline(nbeads,Xin,Y_tmp,X,Y(ivar,ibead))
        enddo
      enddo
      return
      end
      

!
! Akima cubic spline interpolation - This style avoids oscillations
! See: Akima, Hiroshi
! A New Method of Interpolation and Smooth Curve Fitting Based on Local Procedures
! J. ACM 17(4)
! October 1970
! 589-602
! doi = 10.1145/321607.321609
!

      subroutine zts_akima_spline(ndim,x,y,xin,yout)
      implicit none
! Input
      integer ndim
      real*8 x(ndim)
      real*8 y(ndim)
      real*8 xin
! Output
      real*8 yout
! Local
      integer low, high, try, i
      real*8 dydx(5),ddydx(4),t(2)
      real*8 tmp, h, dx

! First we find where we are
      low=1
      high=ndim
      do while (high-low .gt. 1)
        try=(high+low)/2
        if(x(try).gt.xin)then
          high=try
        else
          low=try
        endif
      enddo
      try = high ! Value to the right, not left

      if(try.eq.2) then
          do i=3,5
            dydx(i) = (y(try-3+i)-y(try-4+i)) / (x(try-3+i)-x(try-4+i))
          end do
          dydx(2) = 2.0d0*dydx(3) - dydx(4)
          dydx(1) = 2.0d0*dydx(2) - dydx(3)
      else if(try.eq.3) then
          do i=2,5
            dydx(i) = (y(try-3+i)-y(try-4+i)) / (x(try-3+i)-x(try-4+i))
          end do
          dydx(1) = 2.0d0*dydx(2) - dydx(3)
      else if(try.eq.ndim) then
          do i=1,3
            dydx(i) = (y(try-3+i)-y(try-4+i)) / (x(try-3+i)-x(try-4+i))
          end do
          dydx(4) = 2.0d0*dydx(3) - dydx(2)
          dydx(5) = 2.0d0*dydx(4) - dydx(3)
      else if(try.eq.ndim-1) then
          do i=1,4
            dydx(i) = (y(try-3+i)-y(try-4+i)) / (x(try-3+i)-x(try-4+i))
          end do
          dydx(5) = 2.0d0*dydx(4) - dydx(3)
      else
          do i=1,5
            dydx(i) = (y(try-3+i)-y(try-4+i)) / (x(try-3+i)-x(try-4+i))
          end do
      end if
      do i=1,4
          ddydx(i) = abs(dydx(i+1) - dydx(i))
      end do
      do i=1,2
          tmp = ddydx(i+2) + ddydx(i)
          if(tmp .eq. 0) then
            t(i) = 0.0d0
          else
            t(i) = (ddydx(i+2)*dydx(i+1)+ddydx(i)*dydx(i+2))/tmp
          end if
      end do
      h = x(try)-x(try-1)
      dx = xin - x(try-1)
      yout = y(try-1)
     1      + (dx)   *(t(1))
     2      + (dx**2)*(3.0d0*dydx(3)-2.0d0*t(1)-t(2))/h
     3      + (dx**3)*(t(1)+t(2)-2.0d0*dydx(3))      /(h**2)
    
      return
      end

      subroutine zts_guess(natoms,nbeads,middle,coords,guess)
      implicit none
! Input
      integer natoms, nbeads, middle,guess
! Input/Output
      real*8 coords(3,natoms,nbeads)
! Local variables
      integer i,j,k
      real*8 dgg
      real*8 a,b,c
      real*8 pi
      parameter (pi=3.141592653589793238462643383279d0)

!     Simple linear interpolation

      if (guess .eq. 1) then
        do i = 2, middle-1
         dgg = dble(i-1)/dble(middle-1)
         do j = 1, natoms
           do k = 1, 3
            coords(k,j,i) = coords(k,j,1)*(1-dgg) +
     &                      coords(k,j,middle)*dgg
           enddo
         enddo
        enddo
        do i = middle+1, nbeads-1
         dgg = dble(i-middle)/dble(nbeads-middle)
         do j = 1, natoms
           do k = 1, 3
            coords(k,j,i) = coords(k,j,middle)*(1-dgg) +
     &                      coords(k,j,nbeads)*dgg
           enddo
         enddo
        enddo
      else if (guess .eq. 2) then ! Quadratic in x,y,z
        do j = 1, natoms
          do k = 1, 3
            dgg = dble(middle-1)/dble(nbeads-1)
            a=coords(k,j,1)
            b=((coords(k,j,middle)-coords(k,j,1))-
     &         (coords(k,j,nbeads)-coords(k,j,1))*dgg**2) /
     &        (dgg-dgg**2) 
            c=(coords(k,j,nbeads)-coords(k,j,1)) - b
            do i = 2, nbeads-1
              dgg = dble(i-1)/dble(nbeads-1)
              coords(k,j,i) = a+b*dgg+c*dgg**2
            enddo
          enddo
        enddo
      else if (guess .eq. 3) then ! Quadratic both parts with dx=0 at ends
        ! this just turns out to be linear with more points need ends
        do i = 2, middle-1
         dgg = dble(i-1)/dble(middle-1)
         do j = 1, natoms
           do k = 1, 3
            a = coords(k,j,1)
            b = coords(k,j,middle) - coords(k,j,1)
            coords(k,j,i) = a + b*(3*dgg**2 + 2*dgg**3)
           enddo
         enddo
        enddo
        do i = middle+1, nbeads-1
         dgg = dble(i-middle)/dble(nbeads-middle)
         do j = 1, natoms
           do k = 1, 3
            a = coords(k,j,middle)
            b = coords(k,j,nbeads) - coords(k,j,middle)
            coords(k,j,i) = a + b*(3*dgg**2 + 2*dgg**3)
           enddo
         enddo
        enddo
      else if (guess .eq. 4) then ! f=a+b*x+c*sin(pi*x)
        do j = 1, natoms
          do k = 1, 3
            dgg = dble(middle-1)/dble(nbeads-1)
            a=coords(k,j,1)
            b=coords(k,j,nbeads) - coords(k,j,1)
            c=((coords(k,j,middle)-coords(k,j,1)) - b*dgg)/sin(pi*dgg)
            do i = 2, nbeads-1
              dgg = dble(i-1)/dble(nbeads-1)
              coords(k,j,i) = a+b*dgg+c*sin(pi*dgg)
            enddo
          enddo
        enddo
      else
        call errquit('zts_guess: Bad Guess type', guess, 0)
      endif

      return
      end

      subroutine zts_guessall(natoms,nbeads,coords,geom)
      implicit none
!     Input
      integer natoms, nbeads, geom
!     Input/Output
      real*8 coords(3,natoms,nbeads)
!     Local variabls
      integer i, ix, iix, k, j
      integer MAX_N
      PARAMETER (MAX_N=750)
      logical done(MAX_N)
      logical badgeom
      real*8 dgg

      if (nbeads .gt. MAX_N)
     &     call errquit('zts_meps: MAX_N exceeded', 0, 0)
      done(1) = .true.
      do i = 2, nbeads - 1
        done(i) = .false.
      enddo
      done(nbeads) = .true.
 625  continue ! Loop start
        ix = -1
        iix= -1
        do i = 1, nbeads ! Find un-done region
          if (ix .eq. -1 .and. .not. done(i)) then
            ix = i-1
          endif
          if (ix .ne. -1 .and. done(i)) then
            iix = i
            goto 626
          endif
        enddo
        goto 627 ! Nothing found, we are done
 626    continue
        i = (iix+ix)/2
        dgg = dble(i-ix)/dble(iix-ix)
        do j = 1, natoms
          do k = 1, 3
            coords(k,j,i) = coords(k,j,ix )*(1-dgg) +
     &                      coords(k,j,iix)*dgg
          enddo
        enddo
        call zts_Robinson(natoms,geom,coords(1,1,i),badgeom)
        done(i) = .true.
        goto 625
      ! End of loop
 627  continue

        return
        end
c $Id$
