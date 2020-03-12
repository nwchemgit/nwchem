      subroutine ga_lkain_2cpl3_ext(rtdb, 
     &                          g_x, g_b, 
     &                          g_x_im, g_b_im,
     &                          product, precond, 
     $                          tol, mmaxsub, maxiter, 
     &                          odiff, oprint, omega, limag,
     &                          lifetime, gamwidth, ncomp)

c     $Id$
c  Written by J. Autschbach, SUNY Buffalo
c  Improvements made
c          by F. Aquino,     Northwestern University 
c          03-15-12
c  Note.- Modifying/Improving ga_lkain_2cpl3()

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
      integer ncomp             ! [input] no. of components to treat
      integer rtdb              ! [input] database handle
      integer g_x(ncomp)        ! [input/output] Initial guess/solution
      integer g_x_im(ncomp)     ! not used
      integer g_b(ncomp)        ! [input] Right-hand side vectors 
      integer g_b_im(ncomp)     ! not used
      double precision omega    ! [input] coupling parameter
      logical limag             ! [input] imaginary perturbation?
      logical lifetime          ! [input] consider damping or not?
      double precision gamwidth ! [input] damping parameter

      external product          ! [input] product routine
      external precond          ! [input] preconditioner routine
      double precision tol      ! [input] convergence threshold
      integer mmaxsub           ! [input] maximum subspace dimension
      integer maxiter           ! [input] maximum no. of iterations
      logical odiff             ! [input] use differences in product
      logical oprint            ! [input] print flag

      integer ipm
c
c     Solves the linear equations A(X)=0 for multiple vectors.
c
c ... jochen:
c     This is a modified version of ga_lkain from file ga_it2.F
c     This version allows to solve a coupled set of equations, i.e.
c     there are two right-hand vectors and two initial guesses and two
c     solutions which are coupled. The coupling is mediated by a
c     parameter omega in the call to the preconditioner 
c     (elsewhere, omega is simply called "frequency")
c
c     note: when called from cphf_solve3, odiff = .false. on input
c
c     call product(acc,g_x, g_Ax)
c     . acc is the accuracy trequired for each element of the product
c     . g_x contains the vectors and g_Ax should be filled
c     .     with the product vectors.  The no. of vectors (columns) in
c     . g_x might differ from the no. of vectors input to ga_lkain().
c
c     call precond(g_x,shift)
c     . apply preconditioning directly to the vectors in g_x with the
c     . coupling parameter omega
c
c     On input g_x should contain an initial guess.  It returns the
c     solution.
c
c     maxsub should be at least 3*nvec and can be beneficially increased
c     to about 10*nvec.
c
c     Needs to be extended to store the sub-space vectors out-of-core
c     at least while the product() routine is being executed.

      integer iter, n, n2, nvec, nsub, isub, type, maxsub, 
     &        ntmp1, ntmp2

c ... jochen: for convenience, now most arrays have two components.
c     that might be changed later if memory becomes an issue
      integer g_y,g_Ay,g_r2,
     &        g_Ax(ncomp),g_r(ncomp),  
     &        g_xold(ncomp),g_Axold(ncomp),g_Ax_im(ncomp)
      double precision rmax,acc
      logical converged
      logical odebug,debug,converge_precond
      logical debug1

      character*255 filestub,filesoln
      character*4 digit4
      logical  file_write_ga, file_read_ga
      external file_write_ga, file_read_ga
c
      logical solver_restart
      external solver_restart
c 
      logical do_restart
c
c     =================================================================

      debug = (.false. .and. ga_nodeid().eq.0) ! for code development
      debug1=.false.

c     check input key if we should check for convergence
c     after the preconditioner has been applied to the residual
      if (.not. rtdb_get(rtdb, 'aoresponse:precond',    mt_log, 1,
     &                            converge_precond))
     &  converge_precond = .false.
      
      if (debug) write (6,*) 'ga_lkain_2cpl3 omega =',omega
      if (debug) write (6,*) 'ga_lkain_2cpl3 limag =',limag
      if (debug) write (6,*) 'ga_lkain_2cpl3 lifetime,gamwidth',
     &   lifetime,gamwidth
      if (debug) write (6,*) 'ga_lkain_2cpl3 converge_precond',
     &   converge_precond

      if (lifetime) call errquit('ga_lkain_2cpl3 called with damping',
     &   0,UNKNOWN_ERR)
c     
      odebug = util_print('debug lsolve', print_never) .and. 
     $   ga_nodeid().eq.0
      if (.not. rtdb_get(rtdb, 'cphf:acc',    mt_dbl, 1,
     &                            acc)) acc = 1d-4*tol
c     
      call ga_inquire(g_x(1), type, n, nvec)

      if (ncomp.gt.1) then
        call ga_inquire(g_x(2), type, ntmp1, ntmp2)      
c       ... jochen: do a sanity check on the array dimensions
        if (ntmp1.ne.n .or. ntmp2.ne.nvec) call errquit
     &     ('ga_lkain_2cpl:inconsistent dimensions of g_x components',
     &     nvec,CALC_ERR)
      endif
      
c     later we combine the two components to vectors of double
c     length if we have two components, otherwise not:
      n2 = n
      if (ncomp.gt.1) n2 = n+n                  

      maxsub = mmaxsub          ! So don't modify input scalar arg
      if (maxsub .lt. 3*nvec) maxsub = 3*nvec
      maxsub = (maxsub/nvec)*nvec
c     
      if (oprint .and. ga_nodeid().eq.0) then
        write(6,1) n2, nvec, maxsub, tol, util_wallsec()
    1   format(//,'Iterative solution of linear equations',/,
     $     '  No. of variables', i9,/,
     $     '  No. of equations', i9,/,
     $     '  Maximum subspace', i9,/,
     $     '       Convergence', 1p,d9.1,/,
     $     '        Start time', 0p,f9.1,/)
        call util_flush(6)
      end if
c     
      do ipm = 1,ncomp
        if (.not. ga_create(MT_DBL, n, nvec, 'lkain_2cpl: Ax',
     $     0, 0, g_Ax(ipm)))
     $     call errquit('lkain: failed allocating Ax', nvec,
     &     GA_ERR)
        if (.not. ga_create(MT_DBL, n, nvec, 'lkain_2cpl: r',
     $     0, 0, g_r(ipm)))
     $     call errquit('lkain_2cpl: failed allocating r', nvec,
     &     GA_ERR)
        if (odiff) then
          if (.not. ga_create(MT_DBL, n, nvec, 'lkain_2cpl: xold',
     $       0, 0, g_xold(ipm)))
     $       call errquit('lkain: failed allocating xold', nvec,
     &       GA_ERR)
          if (.not. ga_create(MT_DBL, n, nvec, 'lkain_2cpl: xold',
     $       0, 0, g_Axold(ipm)))
     $       call errquit('lkain: failed allocating Axold', nvec,
     &       GA_ERR)
          call ga_zero(g_xold(ipm))
          call ga_zero(g_Axold(ipm))
        end if                  ! odiff
        call ga_zero(g_Ax(ipm))
        call ga_zero(g_r(ipm))
c       
      enddo                     ! ipm = 1,ncomp
      
c     allocate g_y, g_Ay, and g_r2 with dimension n2 to hold
c     the number of components 
      if (.not. ga_create(MT_DBL, n2, maxsub, 'lkain_2cpl: Y', 
     $   0, 0, g_y))
     $   call errquit('lkain: failed allocating subspace', maxsub,
     &   GA_ERR)
      if (.not. ga_create(MT_DBL, n2, maxsub, 'lkain_2cpl: Ay', 
     $   0, 0, g_Ay))
     $   call errquit('lkain: failed allocating subspace2', maxsub,
     &   GA_ERR)
      if (.not. ga_create(MT_DBL, n2, nvec, 'lkain_2cpl: r2',
     $   0, 0, g_r2))
     $   call errquit('lkain_2cpl: failed allocating r2', nvec,
     &   GA_ERR)
      
      call ga_zero(g_y)
      call ga_zero(g_Ay)
      call ga_zero(g_r2)
      call ga_sync()
c
c     Solution file
c
      if (.not. rtdb_cget(rtdb, 'solver:filestub', 1, filestub))
     &       filestub = 'lkain_soln'
      if (.not. rtdb_cget(rtdb, 'solver:filesoln', 1, filesoln))
     &       filesoln = 'lkain_soln'
#if 0
      call util_file_name(filestub,.false.,.false.,filesoln)
#else
      call cphf_fname(filestub,filesoln)
#endif
      if (ga_nodeid().eq.0) write(luout,*) "ga_lkain filestub:",filestub
      if (ga_nodeid().eq.0) write(luout,*) "ga_lkain filesoln:",filesoln
c
c     Check if this is a restart
c
      if (solver_restart(rtdb)) then
         do_restart = .true.
c         write(6,*) ' attempt reading restart '
         do ipm = 1,ncomp   
            write(digit4,'(".",i3.3)') ipm
            if(.not.file_read_ga(
     P               filesoln(1:inp_strlen(filesoln))//digit4,g_x(ipm)
     C               )) call errquit
     $     ('ga_lkain:could not read solution',1, DISK_ERR)
           enddo
        if (do_restart) then
          if (ga_nodeid().eq.0) 
     &     write(luout,*) "Restarting solution from: ", 
     P               filesoln(1:inp_strlen(filesoln))//digit4
        else
          if (ga_nodeid().eq.0)
     &     write(luout,*) "Error in restart solution: ", filesoln
        end if  ! do_restart
      end if  ! solver_restart
c     
      if (oprint .and. ga_nodeid().eq.0) then
        write(6,2)
        call util_flush(6)
    2   format(/
     $     '   iter   nsub   residual    time ',/,
     $     '   ----  ------  --------  --------- ')
      end if
c     
      nsub = 0
      converged = .false.
c     
c     ---------------------
c     start interation loop
c     ---------------------
c     
      do iter = 1, maxiter       
c       
c ... jochen: here in the iteration loops we keep track
c       of two components of the solution vector, ipm = 1 and 2
c       (ipm stands for + (plus) and - (minus) components)
c       
        if (odiff) then
          do ipm = 1,ncomp   
            call ga_add(1.0d0,g_x(ipm), 
     &                 -1.0d0,g_xold(ipm),
     &                        g_x(ipm))
            call ga_sync()
          enddo ! end-loop-ncomp
        endif
c       
c ... jochen: call product routine with initial or intermediate
c       solution vector: g_x and g_Ax MUST have two components here
        
        if (debug) write (6,*) 'calling product from ga_lkain_2cpl'

        call product(acc, 
     &               g_x   , g_Ax, 
     &               g_x_im, g_Ax_im, 
     &               omega, limag,
     &               lifetime, gamwidth, ncomp)

        if (debug) write (6,*) 'returning product from ga_lkain_2cpl'

c       g_r is zeroed below so we should make sure to do the same
c       with g_r2 here
        call ga_zero(g_r2)
        
        do ipm = 1,ncomp         
          if (odiff) then
            call ga_add(1.0d0, g_Ax(ipm), 
     &                  1.0d0, g_Axold(ipm), 
     &                         g_Ax(ipm))
            call ga_add(1.0d0, g_x(ipm),  
     &                  1.0d0, g_xold(ipm),  
     &                         g_x(ipm))
            call ga_sync()
            call ga_copy(g_x(ipm), g_xold(ipm))
            call ga_copy(g_Ax(ipm), g_Axold(ipm))
          end if
          call ga_zero(g_r(ipm))
c
c         g_Ax = g_b if the system is solved. During the first cycle,
c         g_Ax is calculated from the initial guess
          call ga_add(1.0d0, g_b(ipm),
     &               -1.0d0, g_Ax(ipm), 
     &                       g_r(ipm)) ! The residual        
        enddo                   ! ipm = 1,ncomp
 
c        if (ga_nodeid().eq.0)
c     &   write(*,*) 'FA BEF get_precond_rmax'

        call get_precond_rmax_re(
     &                  rmax,    ! out: max(g_r,g_r_im)   
     &                  g_r,     ! in : real part of g_zr
     &                  g_Ax,    ! in : real part of g_Az
     &                  precond, ! in : name of preconditioner routine
     &        converge_precond,  ! in : =.true. prec->max
     &                  omega,   ! in : omega
     &                  ncomp,   ! in : nr. components
     &                  iter,    ! in : nr. iteration
     &                  debug1)  ! in : =.true. -> allow debug printouts

c        if (ga_nodeid().eq.0)
c     &   write(*,*) 'FA AFT get_precond_rmax'

c       JEM: Putting rmax into rtdb
        if (.not. rtdb_put(rtdb, 'lkain:rmax', mt_dbl, 1, rmax))
     $    call errquit('ga_lkain_2cpl3_ext: rmax put failed', 1,
     $                 RTDB_ERR) 

        if (oprint .and. ga_nodeid().eq.0) then
          write(6,3) iter, nsub+nvec, rmax, util_wallsec()
          call util_flush(6)
    3     format(' ', i5, i7, 3x,1p,d9.2,0p,f10.1,5x,i3)
        end if

c       stop iterations if residual is smaller than criterion
        do ipm = 1,ncomp
           write(digit4,'(".",i3.3)') ipm
           if(.not.file_write_ga(
     P              filesoln(1:inp_strlen(filesoln))//digit4,g_x(ipm)
     C          )) call errquit
     $          ('ga_lkain:could not write solution',1, DISK_ERR)
        enddo
        if (rmax .lt. tol) then
          converged = .true.
          goto 100
        end if
        
c       Copy the vectors to the subspace work area      
        call updating_Az1_z1_zr1(
     &                      g_Ay,  ! in/ou:
     &                      g_y,   ! in/ou:
     &                      g_r2,  ! in/ou:
     &                      g_Ax,  ! in   :
     &                      g_x,   ! in   :
     &                      g_r,   ! in   :
     &                      nvec,  ! in   :
     &                      ncomp, ! in   :
     &                      nsub,  ! in   :
     &                      n)     ! in   :
        
       nsub = nsub + nvec

       call solve_xlineq(
     &              g_x,   ! in/out: updated solution
     &              g_Ay,  ! in    : history of g_Az
     &              g_y,   ! in    : history of g_z
     &              g_r2,  ! in    : history of g_zr
     &              nsub,  ! in    : subspace length
     &              nvec,  ! in    : increment of subspace
     &              ncomp, ! in    : nr. components
     &              n,     ! in    : nr. elements per comp.
     &              iter,  ! in    : iteration nr.
     &              debug1)! in   : =.true. show debug printouts

        if (nsub .eq. maxsub) then
c       
c       Reduce the subspace as necessary
c       
c ====== FA: left-shifting patch ==== START
c Note.- matrices Ay,y shift to left nvec positions
c        removing leftmost patch of dimension: n4 x nvec
c         if (ga_nodeid().eq.0)
c     &    write(*,*) 'FA-matrix-nvec-left-shifting:'
         do isub = nvec+1, maxsub, nvec
          call ga_copy_patch('n',g_Ay,1,n2,isub,isub+nvec-1, 
     $                           g_Ay,1,n2,isub-nvec,isub-1)
          call ga_copy_patch('n',g_y ,1,n2,isub,isub+nvec-1, 
     $                           g_y ,1,n2,isub-nvec,isub-1)
         enddo ! end-loop-isub
c ====== FA: left-shifting patch ==== END
         nsub = nsub - nvec
        end if                  ! (nsub .eq. maxsub)       
      end do                    ! iter = 1,maxiter
  100 continue                  ! jump here if converged   
c     deallocate workspace:     
c
c        Save intermediate solution
c
      do ipm = 1,ncomp
         write(digit4,'(".",i3.3)') ipm
         if(.not.file_write_ga(
     P               filesoln(1:inp_strlen(filesoln))//digit4,g_x(ipm)
     C               )) call errquit
     $     ('ga_lkain:could not write solution',1, DISK_ERR)
      enddo
      do ipm = 1,ncomp
        if (odiff) then
          if (.not. ga_destroy(g_xold(ipm))) call errquit
     &       ('lkain_2cpl: destroy',1, GA_ERR)
          if (.not. ga_destroy(g_Axold(ipm))) call errquit
     &       ('lkain_2cpl: destroy',2,GA_ERR)
        end if
        if (.not. ga_destroy(g_Ax(ipm))) call errquit
     &     ('lkain_2cpl: destroy',20, GA_ERR)
        if (.not. ga_destroy(g_r(ipm))) call errquit
     &     ('lkain_2cpl: destroy',5, GA_ERR)
c       
      enddo                     ! ipm = 1,2
      
      if (.not. ga_destroy(g_Ay)) call errquit
     &   ('lkain_2cpl: destroy Ay',3, GA_ERR)
      if (.not. ga_destroy(g_y)) call errquit
     &   ('lkain_2cpl: destroy r',4, GA_ERR)
      if (.not. ga_destroy(g_r2)) call errquit
     &   ('lkain_2cpl: destroy r2',6, GA_ERR)

      if (.not. converged) then
        if (ga_nodeid().eq.0) then
          write (luout,*) 'WARNING: CPKS procedure is NOT converged'
          write (luout,*) '  I will proceed, but check your results!'
        endif
      endif    
      end

      subroutine copy_r2tor(g_r2,
     &                      g_r,
     &                      g_r_im,
     &                      ncomp,
     &                      nvec,
     &                      n,
     &                      lifetime)
c
c Author: Fredy W. Aquino, Northwestern University
c Date  : 03-15-12

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
       integer ncomp
       integer g_r2,g_r(ncomp),g_r_im(ncomp)
       integer nvec,n,m1,m2,ipm
       logical lifetime
        m1=1
        m2=n
        do ipm=1,ncomp
         call ga_copy_patch('n',g_r2    ,m1,m2,1,nvec, 
     $                          g_r(ipm),1 ,n ,1,nvec)
         m1=m1+n
         m2=m2+n
        enddo !end-loop-ipm
        if (lifetime) then
        do ipm=1,ncomp
         call ga_copy_patch('n',g_r2       ,m1,m2,1,nvec, 
     $                          g_r_im(ipm),1 ,n ,1,nvec)
         m1=m1+n
         m2=m2+n
        enddo !end-loop-ipm
        endif ! end-if-lifetime
      return
      end

      subroutine copy_rtor2(g_r2,
     &                      g_r,
     &                      g_r_im,
     &                      ncomp,
     &                      nvec,
     &                      n,
     &                      lifetime)
c
c Author: Fredy W. Aquino, Northwestern University
c Date  : 03-15-12

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
       integer ncomp
       integer g_r2,g_r(ncomp),g_r_im(ncomp)
       integer nvec,n,m1,m2,ipm
       logical lifetime
        m1=1
        m2=n
        do ipm=1,ncomp
         call ga_copy_patch('n',g_r(ipm),1 ,n ,1,nvec,
     &                          g_r2    ,m1,m2,1,nvec)
         m1=m1+n
         m2=m2+n
        enddo !end-loop-ipm
        if (lifetime) then
        do ipm=1,ncomp
         call ga_copy_patch('n',g_r_im(ipm),1 ,n ,1,nvec,
     &                          g_r2       ,m1,m2,1,nvec)
         m1=m1+n
         m2=m2+n
        enddo !end-loop-ipm
        endif ! end-if-lifetime
      return
      end

      subroutine copy_AxxtoAyy(g_Ax,g_Ax_im,
     &                         g_x,g_x_im,
     &                         g_Ay,g_y,
     &                         nvec,
     &                         ncomp,
     &                         nsub,
     &                         n,
     &                         lifetime)
c
c Author: Fredy W. Aquino, Northwestern University
c Date  : 03-15-12

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
       integer ncomp
       integer g_Ax(ncomp),g_Ax_im(ncomp), 
     &         g_x(ncomp),g_x_im(ncomp),
     &         g_Ay,g_y
       integer nvec,n,p1,p2,m1,m2,nsub,ipm
       logical lifetime
        p1=nsub+1
        p2=nsub+nvec
        m1=1
        m2=n
        do ipm=1,ncomp
         call ga_copy_patch('n',g_Ax(ipm),1 ,n ,1 ,nvec, 
     $                          g_Ay     ,m1,m2,p1,p2)  
         call ga_copy_patch('n',g_x(ipm) ,1 ,n ,1 ,nvec, 
     $                          g_y      ,m1,m2,p1,p2)   
         m1=m1+n
         m2=m2+n    
        enddo ! end-loop-ipm
        if (lifetime) then
        do ipm=1,ncomp
         call ga_copy_patch('n',g_Ax_im(ipm),1 ,n ,1 ,nvec, 
     $                          g_Ay        ,m1,m2,p1,p2)  
         call ga_copy_patch('n',g_x_im(ipm) ,1 ,n ,1 ,nvec, 
     $                          g_y         ,m1,m2,p1,p2)   
         m1=m1+n
         m2=m2+n    
        enddo ! end-loop-ipm
        endif ! end-if-lifetime
      return
      end

      subroutine update_g_x1(g_r2,
     &                      g_x,
     &                      g_x_im,
     &                      ncomp,
     &                      nvec,
     &                      n,
     &                      lifetime)
c
c Author: Fredy W. Aquino, Northwestern University
c Date  : 03-15-12

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
      integer ipm,ncomp,nvec,n,m1,m2
      integer g_x(ncomp),g_x_im(ncomp),
     &        g_r2
      logical lifetime
         m1=1
         m2=n
        do ipm=1,ncomp
         call ga_add_patch(1.0d0,g_r2    ,m1,m2,1,nvec, 
     $                     1.0d0,g_x(ipm),1 ,n ,1,nvec,
     $                           g_x(ipm),1 ,n ,1,nvec)
         m1=m1+n
         m2=m2+n
        enddo !end-loop-ipm
        if (lifetime) then
        do ipm=1,ncomp
         call ga_add_patch(1.0d0,g_r2       ,m1,m2,1,nvec, 
     $                     1.0d0,g_x_im(ipm),1 ,n ,1,nvec,
     $                           g_x_im(ipm),1 ,n ,1,nvec)
         m1=m1+n
         m2=m2+n
        enddo !end-loop-ipm
        endif ! end-if-lifetime
      return
      end

      subroutine ga_lkain_2cpl3_damp_cmplx(
     &                   rtdb, 
     &                   g_x, 
     &                   g_b, 
     &                   g_x_im, 
     &                   g_b_im,
     &                   product,precond, ! in: routines
     $                   tol,
     &                   mmaxsub,
     &                   maxiter, 
     &                   odiff,oprint, 
     &                   omega, 
     &                   limag,
     &                   lifetime,  ! damp means complex, it is redundant
     &                   gamwidth, 
     &                   ncomp)    ! ncomp=2 (+/-)
c
c Purpose: Getting (g_x,g_x_im) by solving recursively a complex linear equation
c          and reducing more memory cost.
c          --> Modified from ga_lkain_2cpl3()
c Author : Fredy W. Aquino, Northwestern University
c Date   : 03-15-12

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "stdio.fh"
#include "rtdb.fh"
      integer ncomp             ! [input] no. of components to treat
      integer rtdb              ! [input] database handle
      integer g_x(ncomp)        ! [input/output] Initial guess/solution Re
      integer g_x_im(ncomp)     ! [input/output] Initial guess/solution Im
      integer g_b(ncomp)        ! [input] Right-hand side vectors Re
      integer g_b_im(ncomp)     ! [input] Right-hand side vectors Im
      double precision omega    ! [input] coupling parameter
      logical limag             ! [input] imaginary perturbation?
      logical lifetime          ! [input] consider damping or not?
      double precision gamwidth ! [input] damping parameter
      external product          ! [input] product routine
      external precond          ! [input] preconditioner routine
      double precision tol      ! [input] convergence threshold
      integer mmaxsub           ! [input] maximum subspace dimension
      integer maxiter           ! [input] maximum no. of iterations
      logical odiff             ! [input] use differences in product
      logical oprint            ! [input] print flag
c
c     Solves the linear equations A(X)=0 for multiple vectors.
c
c ... jochen:
c     This is a modified version of ga_lkain from file ga_it2.F
c     This version allows to solve a coupled set of equations, i.e.
c     there are two right-hand vectors and two initial guesses and two
c     solutions which are coupled. The coupling is mediated by a
c     parameter omega in the call to the preconditioner 
c     (elsewhere, omega is simply called "frequency")
c
c ... jochen: the above comment is from ga_lkain_2cpl3. This here is
c     a modified version of that routine and takes care of a real and an
c     imaginary part for each frequency component. I.e. now arrays
c     have four components ... 
c
c     note: when called from cphf_solve3, odiff = .false. on input
c
c     call product(acc,g_x, g_Ax)
c     . acc is the accuracy trequired for each element of the product
c     . g_x contains the vectors and g_Ax should be filled
c     .     with the product vectors.  The no. of vectors (columns) in
c     . g_x might differ from the no. of vectors input to ga_lkain().
c
c     call precond(g_x,shift)
c     . apply preconditioning directly to the vectors in g_x with the
c     . coupling parameter omega
c
c     On input g_x should contain an initial guess.  It returns the
c     solution.
c
c     maxsub should be at least 3*nvec and can be beneficially increased
c     to about 10*nvec.
c
c     Needs to be extended to store the sub-space vectors out-of-core
c     at least while the product() routine is being executed.
c       
c ... jochen: here in the iteration loops we keep track
c       of two components of the solution vector, ipm = 1 and 2
c       (ipm stands for + (plus) and - (minus) components)
      integer iter,n,n1,
     &        nvec, nsub, isub, type, maxsub, ipm,
     &        ntmp1, ntmp2

c ... jochen: for convenience, now most arrays have two components.
c     that might be changed later if memory becomes an issue
      integer g_xold(ncomp), g_Axold(ncomp)
      integer g_r(ncomp) ,g_r_im(ncomp), 
     &        g_Ax(ncomp),g_Ax_im(ncomp)
      double precision rmax, rmax1, rmax2, acc
      logical converged
      logical odebug, debug, converge_precond
      double complex val_cmplx
      logical debug1
      integer p1,p2,m1,m2,stat_solve
      
      integer g_z(ncomp),g_Az(ncomp),g_zr(ncomp),
     &        g_z1,g_Az1,g_zr1
      double precision omg(2)
      external conv2complex,
     &         updating_Az1_z1_zr1,
     &         solve_zlineq,
     &         solve_zlineq_sep,
     &         get_precond_rmax,
     &         conv2reim

c     later we combine the two components to vecors of double
c     length and combine again Re and Im, i.e. 
c     the dimension is up to 4*n

      call ga_inquire(g_x(1), type, n, nvec) ! get (n,nvec)

      if (ga_nodeid().eq.0) then
       write(*,14) n,ncomp,maxiter,mmaxsub,lifetime
   14  format('(n,ncomp,maxiter,mmaxsub,lifetime)=(',
     &          i3,',',i3,',',i5,',',i5,',',L1,')')
      endif
      n1=ncomp*n

      maxsub = mmaxsub          ! So don't modify input scalar arg
      if (maxsub .lt. 3*nvec) maxsub = 3*nvec
      maxsub = (maxsub/nvec)*nvec 

      if (ga_nodeid().eq.0) then
       write(*,1023) n1,mmaxsub,maxsub,n
 1023  format('(n1,mmaxsub,maxsub,n)=(',
     &        i15,',',i5,',',i5,',',i15,')')
      endif

      debug1 = .false. ! no printouts
c      debug1 = .true. ! allow debugging printouts
c
c     =================================================================

      debug = (.false. .and. ga_nodeid().eq.0) ! for code development

c     check input key if we should check for convergence
c     after the preconditioner has been applied to the residual
      if (.not. rtdb_get(rtdb, 'aoresponse:precond',    mt_log, 1,
     &                            converge_precond))
     &  converge_precond = .false.
      
      if (debug) write (6,*) 'ga_lkain_2cpl_damp omega =',omega
      if (debug) write (6,*) 'ga_lkain_2cpl_damp limag =',limag
      if (debug) write (6,*) 'ga_lkain_2cpl_damp lifetime =',lifetime
      if (debug) write (6,*) 'ga_lkain_2cpl_damp gamwidth =',gamwidth
      if (debug) write (6,*) 'ga_lkain_2cpl_damp ncomp =', ncomp
      if (debug) write (6,*) 'ga_lkain_2cpl3 converge_precond',
     &   converge_precond
c
c     exit if this is the wrong routine to call (lifetime switch
c     must be set)
      if (.not.lifetime) call errquit
     &   ('ga_lkain_2cpl_damp but lifetime=.F.',0,UNKNOWN_ERR)

c     make sure odiff is false (never tested for odiff = .true.)
      if (odiff) call errquit
     &   ('ga_lkain_2cpl_damp odiff=.T.',0,UNKNOWN_ERR)
c     
      odebug = util_print('debug lsolve', print_never) .and. 
     $   ga_nodeid().eq.0
c     
      if (.not. rtdb_get(rtdb, 'cphf:acc',    mt_dbl, 1,
     &                            acc)) acc = 0.01d0*tol
c ------- create (zre,zim) ---------- START
        do ipm=1,ncomp
         if (.not. ga_create(MT_DCPL,n,nvec, 'lkain_2cpl: z',
     $     0, 0, g_z(ipm)))
     $     call errquit('lkain: failed allocating z', nvec,
     &     GA_ERR)
           call ga_zero(g_z(ipm))
         if (.not. ga_create(MT_DCPL,n,nvec, 'lkain_2cpl: Az',
     $     0, 0, g_Az(ipm)))
     $     call errquit('lkain: failed allocating Az', nvec,
     &     GA_ERR)
           call ga_zero(g_Az(ipm))
         if (.not. ga_create(MT_DCPL,n,nvec, 'lkain_2cpl: zr',
     $     0, 0, g_zr(ipm)))
     $     call errquit('lkain: failed allocating zr', nvec,
     &     GA_ERR)
           call ga_zero(g_zr(ipm))
        enddo ! end-loop-ipm
       if (debug1) then
        do ipm=1,ncomp
        if (ga_nodeid().eq.0) 
     &  write(*,*) '--------g_b-re(',ipm,') -------- START'
        call ga_print(g_b(ipm))
        if (ga_nodeid().eq.0) 
     &  write(*,*) '--------g_b-re(',ipm,') -------- END'
        enddo ! end-loop-ipm
        do ipm=1,ncomp
        if (ga_nodeid().eq.0) 
     &  write(*,*) '--------g_b-im(',ipm,') -------- START'
        call ga_print(g_b_im(ipm))
        if (ga_nodeid().eq.0) 
     &  write(*,*) '--------g_b-im(',ipm,') -------- END'
        enddo ! end-loop-ipm
       endif ! end-if-debug1

      if (.not. ga_create(MT_DCPL,n1,maxsub, 'lkain_2cpl: z1', 
     $   0, 0, g_z1))
     $   call errquit('lkain: failed alloc subspace-z1',maxsub,
     &   GA_ERR)
      if (.not. ga_create(MT_DCPL,n1,maxsub, 'lkain_2cpl: Az1', 
     $   0, 0, g_Az1))
     $   call errquit('lkain: failed alloc subspace Az1',maxsub,
     &   GA_ERR)
      if (.not. ga_create(MT_DCPL,n1,nvec, 'lkain_2cpl: zr2',
     $   0, 0, g_zr1))
     $   call errquit('lkain_2cpl: failed allocating zr1', nvec,
     &   GA_ERR)  
      call ga_zero(g_z1)
      call ga_zero(g_Az1)
      call ga_zero(g_zr1)
c ------- create (zre,zim) ---------- END
      if (ncomp.gt.1) then
        call ga_inquire(g_x(2), type, ntmp1, ntmp2)       
c       ... jochen: do a sanity check on the array dimensions
        if (ntmp1.ne.n .or. ntmp2.ne.nvec) call errquit
     &     ('ga_lkain_2cpl:inconsistent dimensions of g_x components',
     &     nvec,CALC_ERR)
      endif     
    
      if (oprint .and. ga_nodeid().eq.0) then
        write(6,1) n1, nvec, maxsub, tol, util_wallsec()
    1   format(//,'Iterative solution of linear equations',/,
     $     '  No. of variables', i9,/,
     $     '  No. of equations', i9,/,
     $     '  Maximum subspace', i9,/,
     $     '       Convergence', 1p,d9.1,/,
     $     '        Start time', 0p,f9.1,/)
        call util_flush(6)
      end if
c     
      do ipm = 1,ncomp
        if (.not. ga_create(MT_DBL, n, nvec, 'lkain_2cpl: Ax',
     $     0, 0, g_Ax(ipm)))
     $     call errquit('lkain: failed allocating Ax', nvec,
     &     GA_ERR)
        if (.not. ga_create(MT_DBL, n, nvec, 'lkain_2cpl: r',
     $     0, 0, g_r(ipm)))
     $     call errquit('lkain_2cpl: failed allocating r', nvec,
     &     GA_ERR)
        if (lifetime) then
          if (.not. ga_create(MT_DBL, n, nvec, 'lkain_2cpl: Ax_im',
     $       0, 0, g_Ax_im(ipm)))
     $       call errquit('lkain: failed allocating Ax_im', nvec,
     &       GA_ERR)
          if (.not. ga_create(MT_DBL, n, nvec, 'lkain_2cpl: r_im',
     $       0, 0, g_r_im(ipm)))
     $       call errquit('lkain_2cpl: failed allocating r_im', nvec,
     &       GA_ERR)
        endif                   ! lifetime
        
        if (odiff) then
c         jochen: this part and all subsequent "odiff" parts were
c         never adapted for the imaginary components
          if (.not. ga_create(MT_DBL, n, nvec, 'lkain_2cpl: xold',
     $       0, 0, g_xold(ipm)))
     $       call errquit('lkain: failed allocating xold', nvec,
     &       GA_ERR)
          if (.not. ga_create(MT_DBL, n, nvec, 'lkain_2cpl: xold',
     $       0, 0, g_Axold(ipm)))
     $       call errquit('lkain: failed allocating Axold', nvec,
     &       GA_ERR)
          call ga_zero(g_xold(ipm))
          call ga_zero(g_Axold(ipm))
        end if                  ! odiff

        call ga_zero(g_Ax(ipm))
        call ga_zero(g_r(ipm))
        if (lifetime) then
          call ga_zero(g_Ax_im(ipm))
          call ga_zero(g_r_im(ipm))
        endif
c       
      enddo                     ! ipm = 1,ncomp
     
      if (oprint .and. ga_nodeid().eq.0) then
        write(6,2)
        call util_flush(6)
    2   format(/
     $     '   iter   nsub   residual    time ',/,
     $     '   ----  ------  --------  --------- ')
      end if
c     
      nsub = 0
      converged = .false.
c     
c     ---------------------
c     start interation loop
c     ---------------------
c     
      do iter = 1, maxiter
        if (odiff) then
          do ipm = 1,ncomp   
            call ga_add( 1.0d0,g_x(ipm), 
     &                  -1.0d0,g_xold(ipm),  
     &                         g_x(ipm))
            call ga_sync()
          enddo
        endif
c       
c ... jochen: call product routine with initial or intermediate
c       solution vector: g_x and g_Ax MUST have dimension two here
c       even if only one of them is used
          if (debug) then   
            do ipm=1,ncomp
             if (ga_nodeid().eq.0) then
              write(*,112) iter,ipm
  112         format('------ prod-g_x-1(',i3,',',i3,')------ START')
             endif
             call ga_print(g_x(ipm))
             if (ga_nodeid().eq.0) then
              write(*,113) iter,ipm
  113         format('------ prod-g_x-1(',i3,',',i3,')------ END')
             endif
            enddo ! end-loop-ipm
         endif ! end-if-debug
        
        if (debug) write (6,*)
     &     'calling product from ga_lkain_2cpl_damp'
        call product(acc, 
     &               g_x,      ! in  : x
     &               g_Ax,     ! out : product A x  
     &               g_x_im,   ! in  : x_im
     &               g_Ax_im,  ! out : product A x_im
     &               omega,    ! in  :
     &               limag,    ! in  :
     &               lifetime, ! in  : =.true. -> x is complex
     &               gamwidth, ! in  :
     &               ncomp)    ! in  : nr. components
         if (debug) then
            do ipm=1,ncomp
             if (ga_nodeid().eq.0) then
              write(*,116) iter,ipm
  116         format('------ prod-g_x-2(',i3,',',i3,')------ START')
             endif
             call ga_print(g_x(ipm))
             if (ga_nodeid().eq.0) then
              write(*,117) iter,ipm
  117         format('------ prod-g_x-2(',i3,',',i3,')------ END')
             endif
             if (ga_nodeid().eq.0) then
              write(*,118) iter,ipm
  118         format('------ prod-g_Ax-2(',i3,',',i3,')------ START')
             endif
             call ga_print(g_Ax(ipm))
             if (ga_nodeid().eq.0) then
              write(*,119) iter,ipm
  119         format('------ prod-g_Ax-2(',i3,',',i3,')------ END')
             endif
             enddo ! end-loop-ipm
           endif ! end-if-debug  

        if (debug) write (6,*)
     &     'returning product from ga_lkain_2cpl_damp'
        
        do ipm = 1,ncomp
          
          if (odiff) then
c           jochen: odiff stuff presently ignored
            call ga_add(1.0d0,g_Ax(ipm), 
     &                  1.0d0,g_Axold(ipm), 
     &                        g_Ax(ipm))
            call ga_add(1.0d0,g_x(ipm),  
     &                  1.0d0,g_xold(ipm),
     &                        g_x(ipm))
            call ga_copy(g_x(ipm), g_xold(ipm))
            call ga_copy(g_Ax(ipm), g_Axold(ipm))
          end if                ! odiff

          call ga_zero(g_r(ipm))
          if (lifetime) call ga_zero(g_r_im(ipm))

c         g_r will be the quantity -Ax + b, i.e. if the equation system
c         Ax = b is solved then this vector will be zero
c
c         During the first cycle,
c         g_Ax is calculated from the initial guess for which the
c         preconditioner has already been applied (to be more clear:
c         we have divided the perturbation matrix elements by orbital
c         energy denominators, including the frequency term, 
c         and assigned real and imaginary parts accordingly)
          call ga_add( 1.0d0,g_b(ipm),
     &                -1.0d0,g_Ax(ipm), 
     &                       g_r(ipm))    ! The residual, Real part
          call ga_add( 1.0d0,g_b_im(ipm),
     &                -1.0d0,g_Ax_im(ipm), 
     &                       g_r_im(ipm)) ! The residual, Im part

          if (debug) then
             if (ga_nodeid().eq.0) then
              write(*,120) iter,ipm
  120         format('------ prod-g_b(',i3,',',i3,')------ START')
             endif
             call ga_print(g_b(ipm))
             if (ga_nodeid().eq.0) then
              write(*,121) iter,ipm
  121         format('------ prod-g_b(',i3,',',i3,')------ END')
             endif
             if (ga_nodeid().eq.0) then
              write(*,122) iter,ipm
  122         format('------ prod-g_r(',i3,',',i3,')------ START')
             endif
             call ga_print(g_r(ipm))
             if (ga_nodeid().eq.0) then
              write(*,123) iter,ipm
  123         format('------ prod-g_r(',i3,',',i3,')------ END')
             endif
          endif ! end-if-debug      
        enddo                   ! ipm = 1,ncomp

c        if (ga_nodeid().eq.0)
c     &   write(*,*) 'FA BEF get_precond_rmax'

        call get_precond_rmax(
     &                  rmax,    ! out: max(g_r,g_r_im)   
     &                  g_r,     ! in : real part of g_zr
     &                  g_r_im,  ! in : imag part of g_zr
     &                  g_Ax,    ! in : real part of g_Az
     &                  g_Ax_im, ! in : imag part of g_Az 
     &                  precond, ! in : name of preconditioner routine
     &        converge_precond,  ! in : =.true. prec->max
     &                  omega,   ! in : omega
     &                  gamwidth,! in : gamwidth
     &                  ncomp,   ! in : nr. components
     &                  iter,    ! in : nr. iteration
     &                  debug1)  ! in : =.true. -> allow debug printouts

c        if (ga_nodeid().eq.0)
c     &   write(*,*) 'FA AFT get_precond_rmax'

c -------- printout per iteration -------------- START
        if (oprint .and. ga_nodeid().eq.0) then
          write(6,3) iter, nsub+nvec, rmax, util_wallsec()
          call util_flush(6)
    3     format(' ', i5, i7, 3x,1p,d9.2,0p,f10.1,5x,i3)
        end if
c -------- printout per iteration -------------- END
c        if (ga_nodeid().eq.0) then
c         write(*,10) iter,nsub+nvec,rmax
c 10      format('(iter,nsub+nvec,rmax)=(',i5,',',i7,',',f15.8,')')       
c        endif
c       stop iterations if residual is smaller than criterion
        if (rmax .lt. tol) then
          converged = .true.
          goto 100
        end if

c ======== get complex: g_z,g_Az,g_zr ====== START
c Do: (x,x_im) -> z
c (x,x_im) -> (k_zre,k_zim) -> z
       call conv2complex(g_z,    ! out: = complx(g_x,g_x_im)
     &                   g_x,    ! in : real      arr
     &                   g_x_im, ! in : imaginary arr
     &                   n,      ! in : n    rows
     &                   nvec,   ! in : nvec columns
     &                   ncomp)  ! in : nr. components
       call conv2complex(g_Az,   ! out: = complx(g_Ax,g_Ax_im)
     &                   g_Ax,   ! in : real      arr
     &                   g_Ax_im,! in : imaginary arr
     &                   n,      ! in : n    rows
     &                   nvec,   ! in : nvec columns
     &                   ncomp)  ! in : nr. components
       call conv2complex(g_zr,   ! out: = complx(g_Ax,g_Ax_im)
     &                   g_r,    ! in : real      arr
     &                   g_r_im, ! in : imaginary arr
     &                   n,      ! in : n    rows
     &                   nvec,   ! in : nvec columns
     &                   ncomp)  ! in : nr. components
c ======== get complex: g_z,g_Az,g_zr ====== END

        if (debug1) then
         do ipm=1,ncomp
          if (ga_nodeid().eq.0) then
           write(*,3001) ipm,iter
 3001      format('---------g_z-0(',i3,',',i3,')-----START')
          endif
          call ga_print(g_z(ipm))
          if (ga_nodeid().eq.0) then
           write(*,3002) ipm,iter
 3002      format('---------g_z-0(',i3,',',i3,')-----END')
          endif
          if (ga_nodeid().eq.0) then
           write(*,3003) ipm,iter
 3003      format('---------g_Az-0(',i3,',',i3,')-----START') 
          endif
          call ga_print(g_Az(ipm))
          if (ga_nodeid().eq.0) then
           write(*,3004) ipm,iter
 3004      format('---------g_Az-0(',i3,',',i3,')-----END')
          endif
          if (ga_nodeid().eq.0) then
           write(*,2800) ipm,iter
 2800      format('---------g_zr-0(',i3,',',i3,')-----START') 
          endif
          call ga_print(g_zr(ipm))
          if (ga_nodeid().eq.0) then
           write(*,2801) ipm,iter
 2801      format('---------g_zr-0(',i3,',',i3,')-----END')
          endif
         enddo ! end-loop-ipm
        endif ! end-if-debug1  
        call updating_Az1_z1_zr1(
     &                      g_Az1, ! in/ou:
     &                      g_z1,  ! in/ou:
     &                      g_zr1, ! in/ou:
     &                      g_Az,  ! in   :
     &                      g_z,   ! in   :
     &                      g_zr,  ! in   :
     &                      nvec,  ! in   :
     &                      ncomp, ! in   :
     &                      nsub,  ! in   :
     &                      n)     ! in   :
          if (debug1) then
          if (ga_nodeid().eq.0) then
           write(*,2700) iter
 2700      format('---------g_z1-0(',i3,')-----START')
          endif
          call ga_print(g_z1)
          if (ga_nodeid().eq.0) then
           write(*,2701) iter
 2701      format('---------g_z1-0(',i3,')-----END')
          endif
          if (ga_nodeid().eq.0) then
           write(*,2702) iter
 2702      format('---------g_Az1-0(',i3,')-----START')
          endif
          call ga_print(g_Az1)
          if (ga_nodeid().eq.0) then
           write(*,2703) iter
 2703      format('---------g_Az1-0(',i3,')-----END')
          endif
          if (ga_nodeid().eq.0) then
           write(*,2704) iter
 2704      format('---------g_zr1-0(',i3,')-----START')
          endif
          call ga_print(g_zr1)
          if (ga_nodeid().eq.0) then
           write(*,2705) iter
 2705      format('---------g_zr1-0(',i3,')-----END')
          endif
          endif ! end-if-debug1

        nsub = nsub + nvec
c ============ doing complex-linear solving ======= START
c        if (ga_nodeid().eq.0)
c     &   write(*,*) 'BEFORE solve_zlineq ...'
         call solve_zlineq(
     &              g_z,   ! in/out: updated solution
     &              g_Az1, ! in    : history of g_Az
     &              g_z1,  ! in    : history of g_z
     &              g_zr1, ! in    : history of g_zr
     &              nsub,  ! in    : subspace length
     &              nvec,  ! in    : increment of subspace
     &              ncomp, ! in    : nr. components
     &              n,     ! in    : nr. elements per comp.
     &              iter,  ! in    : iteration nr.
     &              debug1)! in   : =.true. show debug printouts
c        if (ga_nodeid().eq.0)
c     &   write(*,*) 'AFTER solve_zlineq ...'
c ============ doing complex-linear solving ======= END
c ======= g_z --> (g_x,g_x_im) ========= START
c       if (ga_nodeid().eq.0)
c     &  write(*,*) 'BEF conv2reim-x'
       call conv2reim(g_x,   ! out : real      arr
     &                g_x_im,! out : imaginary arr
     &                g_z,   ! in  : = complx(g_xre,g_xim)
     &                n,     ! in  : n    rows
     &                nvec,  ! in  : nvec columns
     &                ncomp) ! in  : nr. components   
c       if (ga_nodeid().eq.0)
c     &  write(*,*) 'AFT conv2reim-x'
c ======= g_z --> (g_x,g_x_im) ========= END
        if (debug1) then
         do ipm=1,ncomp
          if (ga_nodeid().eq.0) then
           write(*,4000) ipm,iter
 4000      format('---------g_xre-AFT-0(',i3,',',i3,')-----START')
          endif
          call ga_print(g_x(ipm))
          if (ga_nodeid().eq.0) then
           write(*,4001) ipm,iter
 4001      format('---------g_xre-AFT-0(',i3,',',i3,')-----END')
          endif
          if (ga_nodeid().eq.0) then
           write(*,4002) ipm,iter
 4002      format('---------g_xim-AFT-0(',i3,',',i3,')-----START')
          endif
          call ga_print(g_x_im(ipm))
          if (ga_nodeid().eq.0) then
           write(*,4003) ipm,iter
 4003      format('---------g_xim-AFT-0(',i3,',',i3,')-----END')
          endif
         enddo ! end-loop-ipm
        endif ! end-if-debug1

c        if (iter.eq.2) then
c         if (ga_nodeid().eq.0) 
c     &      write(*,*) 'FA-check-zlinear-solver'
c          stop
c        endif

        if (nsub .eq. maxsub) then
c       
c       Reduce the subspace as necessary
c       
c ====== FA: left-shifting patch ==== START
c Note.- matrices Ay,y shift to left nvec positions
c        removing leftmost patch of dimension: n4 x nvec
         if (ga_nodeid().eq.0)
     &    write(*,*) 'FA-matrix-nvec-left-shifting:'
         do isub = nvec+1, maxsub, nvec
          call ga_copy_patch('n',g_Az1,1,n1,isub,isub+nvec-1, 
     $                           g_Az1,1,n1,isub-nvec,isub-1)
          call ga_copy_patch('n',g_z1 ,1,n1,isub,isub+nvec-1, 
     $                           g_z1 ,1,n1,isub-nvec,isub-1)
         enddo ! end-loop-isub
c ====== FA: left-shifting patch ==== END
         nsub = nsub - nvec
        end if                  ! (nsub .eq. maxsub) 
      enddo                     ! iter = 1,maxiter
  100 continue                  ! jump here if converged  
c     deallocate workspace:
c     
      do ipm = 1,ncomp
        if (odiff) then
          if (.not. ga_destroy(g_xold(ipm))) call errquit
     &       ('lkain_2cpl: destroy',1, GA_ERR)
          if (.not. ga_destroy(g_Axold(ipm))) call errquit
     &       ('lkain_2cpl: destroy',2,GA_ERR)
        end if
        if (.not. ga_destroy(g_Ax(ipm))) call errquit
     &     ('lkain_2cpl: destroy',20, GA_ERR)
        if (.not. ga_destroy(g_r(ipm))) call errquit
     &     ('lkain_2cpl: destroy',5, GA_ERR)
        if (lifetime) then
          if (.not. ga_destroy(g_Ax_im(ipm))) call errquit
     &       ('lkain_2cpl: destroy',201, GA_ERR)
          if (.not. ga_destroy(g_r_im(ipm))) call errquit
     &       ('lkain_2cpl: destroy',51, GA_ERR)
        endif      
      enddo                     ! ipm = 1,2

       do ipm=1,ncomp
        if (.not. ga_destroy(g_Az(ipm))) call errquit
     &   ('lkain_2cpl3-cmplx: destroy Az',3, GA_ERR)  
        if (.not. ga_destroy(g_z(ipm))) call errquit
     &   ('lkain_2cpl3-cmplx: destroy z',3, GA_ERR)
        if (.not. ga_destroy(g_zr(ipm))) call errquit
     &   ('lkain_2cpl3-cmplx: destroy zr',3, GA_ERR)        
       enddo ! end-loop-ipm     
      if (.not. ga_destroy(g_Az1)) call errquit
     &   ('lkain_2cpl3-cmplx: destroy Az1',3, GA_ERR)
      if (.not. ga_destroy(g_z1)) call errquit
     &   ('lkain_2cpl3-cmplx: destroy z1',4, GA_ERR)
      if (.not. ga_destroy(g_zr1)) call errquit
     &   ('lkain_2cpl3-cmplx: destroy zr1',6, GA_ERR)    

      if (.not. converged) then
        if (ga_nodeid().eq.0) then
          write (luout,*) 'WARNING: CPKS procedure is NOT converged'
          write (luout,*) '  I will proceed, but check your results!'
        endif
      endif
c     
      end

c ========================================================
c ========= Reduce memory consumption ============== START
      subroutine ga_lkain_2cpl3_damp_cmplx_redmem(
     &                   rtdb,
     &                   g_z, ! in/out: solution
     &                   g_zb,! in    : b (of Ax=b)  
     &                   product, ! in: routine to compute Az
     &                   precond, ! in: routine to do energy scaling Az,r
     $                   tol,
     &                   mmaxsub,
     &                   maxiter, 
     &                   odiff,oprint, 
     &                   omega, 
     &                   limag,
     &                   lifetime,  ! damp means complex, it is redundant
     &                   gamwidth, 
     &                   ncomp,    ! ncomp=2 (+/-)
     &                   npol,
     &                   nvir,
     &                   nocc)
c
c Author : Fredy W. Aquino, Northwestern University
c Purpose: Getting g_z by solving recursively a complex linear equation
c          and reducing memory cost.
c          --> Modified from ga_lkain_2cpl3()
c Date   : 03-15-12

      implicit none
#include "errquit.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
#include "util.fh"
#include "stdio.fh"
#include "inp.fh"
#include "rtdb.fh"
      integer ncomp       ! [input] no. of components to treat
      integer g_z(ncomp)  ! [input/output] Initial guess/solution (Re,Im)      
      integer g_zb(ncomp),! [input] b of Ax=b
     &        g_z1,       ! Scratch GA contains history of z  in (n1,mmaxsub) 
     &        g_Az1,      ! Scratch GA contains history of Az in (n1,mmaxsub) 
     &        g_zr1       ! Scratch GA (r= b-Ax) error of size (n1,nvec)
      integer g_zb1
      integer npol,
     &        nvir(npol),nocc(npol)
c Note.- In g_z1,g_Az1 a (n1,nvec) block is added per iteration.
      integer rtdb              ! [input] database handle
      double precision omega    ! [input] coupling parameter
      logical limag             ! [input] imaginary perturbation?
      logical lifetime          ! [input] consider damping or not?
      double precision gamwidth ! [input] damping parameter
      external product          ! [input] product routine
      external precond          ! [input] preconditioner routine
      double precision tol      ! [input] convergence threshold
      integer mmaxsub           ! [input] maximum subspace dimension
      integer maxiter           ! [input] maximum no. of iterations
      logical odiff             ! [input] use differences in product
      logical oprint            ! [input] print flag
c
c     Solves the linear equations A(X)=0 for multiple vectors.
c
c ... jochen:
c     This is a modified version of ga_lkain from file ga_it2.F
c     This version allows to solve a coupled set of equations, i.e.
c     there are two right-hand vectors and two initial guesses and two
c     solutions which are coupled. The coupling is mediated by a
c     parameter omega in the call to the preconditioner 
c     (elsewhere, omega is simply called "frequency")
c
c ... jochen: the above comment is from ga_lkain_2cpl3. This here is
c     a modified version of that routine and takes care of a real and an
c     imaginary part for each frequency component. I.e. now arrays
c     have four components ... 
c
c     note: when called from cphf_solve3, odiff = .false. on input
c
c     call product(acc,g_x, g_Ax)
c     . acc is the accuracy trequired for each element of the product
c     . g_x contains the vectors and g_Ax should be filled
c     .     with the product vectors.  The no. of vectors (columns) in
c     . g_x might differ from the no. of vectors input to ga_lkain().
c
c     call precond(g_x,shift)
c     . apply preconditioning directly to the vectors in g_x with the
c     . coupling parameter omega
c
c     On input g_x should contain an initial guess.  It returns the
c     solution.
c
c     maxsub should be at least 3*nvec and can be beneficially increased
c     to about 10*nvec.
c
c     Needs to be extended to store the sub-space vectors out-of-core
c     at least while the product() routine is being executed.
c       
c ... jochen: here in the iteration loops we keep track
c       of two components of the solution vector, ipm = 1 and 2
c       (ipm stands for + (plus) and - (minus) components)
      integer iter,n,n1,
     &        nvec, nsub, isub, 
     &        type, maxsub, ipm,
     &        nsub_file

c ... jochen: for convenience, now most arrays have two components.
c     that might be changed later if memory becomes an issue
      double precision rmax, acc
      logical converged, odebug, debug, 
     &        converge_precond, debug1
      double complex val_cmplx
      integer p1,p2,m1,m2,dim1,dim2,nblock
      double complex one_cmplx,mone_cmplx,zero_cmplx
      logical dft_CPHF2_read,
     &        dft_CPHF2_write,
     &        dft_CPHF2_read2fix
      real ran1
      integer status_gasvd,idum
      double precision factor_x
      external solve_zlineq1,
     &         get_precond_rmax_zin,
     &         dft_CPHF2_read,
     &         dft_CPHF2_write
      external copy_complx2real_4redmem,
     &         solve_xlineq_4redmem,
     &         solve_zlineq_KAIN1,
     &         solve_zlineq_KAIN3
      integer g_Ax1,g_x1,g_xr1,slcKAIN

      logical status,flag2readfile
      integer index4cphf,checkorth,cphf3write
      character*255 aorespfilename
      character*(*) lbl_cphfaoresp
      character*255 lbl_cphfaoresp1
      integer iimoderaman,iiistepraman 
      logical ramanspc                
      parameter(lbl_cphfaoresp='aoresp_fiao_f')

      slcKAIN=2 ! Choose KAIN linear solver (1,2,3)

      one_cmplx =dcmplx( 1.0d0,0.0d0)
      mone_cmplx=dcmplx(-1.0d0,0.0d0)
      zero_cmplx=dcmplx( 0.0d0,0.0d0)
      if (.not. rtdb_get(rtdb, 'cphf:cphf3write',    mt_int, 1,
     &                   cphf3write))
     &  cphf3write = 0 ! assigns 0 if unsuccessfull read from rtdb
      if (.not. rtdb_get(rtdb, 'cphf:checkorth',    mt_int, 1,
     &                   checkorth))
     &  checkorth = 0 ! assigns 0 if unsuccessfull read from rtdb
      call ga_inquire(g_z(1),type,n,nvec) ! get (n,nvec) n=sum(nocc*nvirt(i) i=1,npol)
      n1=ncomp*n
      maxsub = mmaxsub          ! So don't modify input scalar arg
      if (maxsub .lt. 3*nvec) maxsub = 3*nvec
      maxsub = (maxsub/nvec)*nvec 
      debug1 = .false. ! no printouts
c
c     =================================================================

      debug = (.false. .and. ga_nodeid().eq.0) ! for code development

c     check input key if we should check for convergence
c     after the preconditioner has been applied to the residual
      if (.not. rtdb_get(rtdb, 'aoresponse:precond',    mt_log, 1,
     &                            converge_precond))
     &  converge_precond = .false.
      
      if (debug) write (6,*) 'ga_lkain_2cpl_damp omega =',omega
      if (debug) write (6,*) 'ga_lkain_2cpl_damp limag =',limag
      if (debug) write (6,*) 'ga_lkain_2cpl_damp lifetime =',lifetime
      if (debug) write (6,*) 'ga_lkain_2cpl_damp gamwidth =',gamwidth
      if (debug) write (6,*) 'ga_lkain_2cpl_damp ncomp =', ncomp
      if (debug) write (6,*) 'ga_lkain_2cpl3 converge_precond',
     &   converge_precond
c
c     exit if this is the wrong routine to call (lifetime switch
c     must be set)
      if (.not.lifetime) call errquit
     &   ('ga_lkain_2cpl_damp but lifetime=.F.',0,UNKNOWN_ERR)

c     make sure odiff is false (never tested for odiff = .true.)
      if (odiff) call errquit
     &   ('ga_lkain_2cpl_damp odiff=.T.',0,UNKNOWN_ERR)
c     
      odebug = util_print('debug lsolve', print_never) .and. 
     $   ga_nodeid().eq.0
c     
      if (.not. rtdb_get(rtdb, 'cphf:acc',    mt_dbl, 1,
     &                            acc)) acc = 0.0001d0*tol
c ------- create (zre,zim) ---------- START

c ++++++ added for solve_zlineq_KAIN1 +++ START
      if (slcKAIN.eq.2) then
       if (.not. ga_create(MT_DBL,2*n1,maxsub, 'lkain_dbl: x1', 
     $   0, 0, g_x1))
     $   call errquit('lkain: failed alloc subspace-x1',maxsub,
     &   GA_ERR)
       if (.not. ga_create(MT_DBL,2*n1,maxsub, 'lkain_dbl: Ax1', 
     $   0, 0, g_Ax1))
     $   call errquit('lkain: failed alloc subspace Az1',maxsub,
     &   GA_ERR)
       if (.not. ga_create(MT_DBL,2*n1,nvec, 'lkain_dbl: xr1',
     $   0, 0, g_xr1))
     $   call errquit('lkain_2cpl: failed allocating zr1', nvec,
     &   GA_ERR)  
      endif
c ++++++ added for solve_zlineq_KAIN1 +++ END

      if (.not. ga_create(MT_DCPL,n1,maxsub, 'lkain_2cpl: z1', 
     $   0, 0, g_z1))
     $   call errquit('lkain: failed alloc subspace-z1',maxsub,
     &   GA_ERR)
      if (.not. ga_create(MT_DCPL,n1,maxsub, 'lkain_2cpl: Az1', 
     $   0, 0, g_Az1))
     $   call errquit('lkain: failed alloc subspace Az1',maxsub,
     &   GA_ERR)
      if (.not. ga_create(MT_DCPL,n1,nvec, 'lkain_2cpl: zr2',
     $   0, 0, g_zr1))
     $   call errquit('lkain_2cpl: failed allocating zr1', nvec,
     &   GA_ERR)  
      if (.not. ga_create(MT_DCPL,n1,nvec, 'lkain_2cpl: zr2',
     $   0, 0, g_zb1))
     $   call errquit('lkain_2cpl: failed allocating zb1', nvec,
     &   GA_ERR)  
      call ga_zero(g_zb1)
      m1=1
      m2=n
      do ipm=1,ncomp
         call ga_copy_patch('n',g_zb(ipm),1 ,n ,1,nvec,
     &                          g_zb1    ,m1,m2,1,nvec)    
       m1=m1+n
       m2=m2+n
      enddo ! end-loop-ipm
      call ga_zero(g_z1)
      call ga_zero(g_Az1)
      call ga_zero(g_zr1)
      call ga_sync()
c ------- create (zre,zim) ---------- END
    
      if (oprint .and. ga_nodeid().eq.0) then
        write(6,1) n1, nvec, maxsub, tol, util_wallsec()
    1   format(//,'Iterative solution of linear equations',/,
     $     '  No. of variables', i9,/,
     $     '  No. of equations', i9,/,
     $     '  Maximum subspace', i9,/,
     $     '       Convergence', 1p,d9.1,/,
     $     '        Start time', 0p,f9.1,/)
        call util_flush(6)
      end if
     
      if (oprint .and. ga_nodeid().eq.0) then
        write(6,2)
        call util_flush(6)
    2   format(/
     $     '   iter   nsub   residual    time ',/,
     $     '   ----  ------  --------  --------- ')
      end if
c     
      nsub = 0
      converged = .false.
c     
c     ---------------------
c     start interation loop
c     ---------------------
c  
c 000000000000000 getting cphf filename to store 00000000 START
      if (.not. rtdb_get(rtdb,'cphf3-aores:guess1',
     & mt_int,1,index4cphf))  index4cphf = 0

      ramanspc=.false.
      status=rtdb_get(rtdb,'raman:aores0',mt_log,1,ramanspc)
      if (ramanspc) then
       if (.not. rtdb_get(rtdb,'raman:aores1', 
     &          mt_int, 1,iimoderaman)) call
     $     errquit('ga_lkain_2cpl3_redmem1: failed to read iimoderaman', 
     &             0, RTDB_ERR)      
       if (.not. rtdb_get(rtdb,'raman:aores2', 
     &          mt_int, 1,iiistepraman)) call
     $   errquit('ga_lkain_2cpl3_redmem1: failed to read iiistepraman', 
     &             0, RTDB_ERR)    
       write(lbl_cphfaoresp1,'(a13,i1,"_",i4.4,"-",i1)') 
     &                      lbl_cphfaoresp,index4cphf,
     &                      iimoderaman,iiistepraman                  
      else
       write(lbl_cphfaoresp1,'(a13,i1)') lbl_cphfaoresp,index4cphf    
      endif
       call util_file_name(lbl_cphfaoresp1,
     &                     .false.,.false.,aorespfilename)

      nsub_file=0 ! reset value in all nodes
      nsub=0      ! reset value in all nodes
      flag2readfile=.false.
      if (.not. dft_CPHF2_read(
     &           aorespfilename, ! in: filename
     &           n,        ! in: sum_{i=1,npol} nocc(i)*nvirt(i)
     &           ncomp,    ! in: nr. components
     &           nvec,     ! in: nr. of directions = 3
     &           n1,       ! in: =n*ncomp
     &           nsub,     ! ou: last subspace index (nsub+1)= nr of subspaces stored
     &           nsub_file,! ou: last subspace read from file
     &           maxsub,   ! in: maximum subspace
     &           g_z1,     ! ou: history matrix z
     &           g_Az1))   ! ou: history matrix Az
     &  then
c       if (ga_nodeid().eq.0)
c     &  write(*,1999) aorespfilename(1:inp_strlen(aorespfilename))
c 1999   format('File ',a,
c     &          ' does not exist, proceed to generate (z1,Az1)')
c ------ g_z0 --> g_z1 ----- START
c Copying initial guess
         nsub_file=0
         nsub=0
         m1=1
         m2=n
         p1=nsub+1
         p2=nsub+nvec
         do ipm=1,ncomp

          if (debug1) then
          if (ga_nodeid().eq.0) then
           write(*,2770) ipm
 2770      format('---------g_z-guess(',i3,')-----START')
          endif
          call ga_print(g_z(ipm))
          if (ga_nodeid().eq.0) then
           write(*,2701) ipm
 2771      format('---------g_z-guess(',i3,')-----END')
          endif
          endif ! end-if-debug1

          call ga_copy_patch('n',g_z(ipm),1 ,n ,1 ,nvec, 
     $                           g_z1    ,m1,m2,p1,p2)
          m1=m1+n
          m2=m2+n          
         enddo ! end-loop-ipm
c ------ g_z0 --> g_z1 ----- END
      else
         call ga_sync()
c Note.- I need to propagate nsub
         call ga_igop(6,nsub_file,1,'+') ! node0 nsub ne 0, eq 0 every other node    
         call ga_igop(6,nsub,1,'+')      ! node0 nsub ne 0, eq 0 every other node    
         flag2readfile=.true.
         nblock=nsub/3+1
         m1=1
         m2=n
         p1=nsub+1
         p2=nsub+nvec
         do ipm=1,ncomp
          call ga_copy_patch('n',g_z1    ,m1,m2,p1,p2,
     &                           g_z(ipm),1 ,n ,1 ,nvec)
     
          if (debug1) then
           if (ga_nodeid().eq.0) then
            write(*,2790) ipm
 2790       format('---------g_z-guess(',i3,')-----START')
           endif
           call ga_print(g_z(ipm))
           if (ga_nodeid().eq.0) then
            write(*,2791) ipm
 2791       format('---------g_z-guess(',i3,')-----END')
            endif
          endif ! end-if-debug1

          m1=m1+n
          m2=m2+n          
         enddo ! end-loop-ipm
         if (debug1) then
          if (ga_nodeid().eq.0)
     &      write(*,*) '---------g_z1-read-from-file-----START'
          call ga_print(g_z1)
          if (ga_nodeid().eq.0)
     &      write(*,*) '---------g_z1-read-from-file-----END'
          if (ga_nodeid().eq.0)
     &      write(*,*) '---------g_Az1-read-from-file-----START'
          call ga_print(g_Az1)
          if (ga_nodeid().eq.0)
     &      write(*,*) '---------g_Az1-read-from-file-----END'
         endif ! end-if-debug1
      endif  
c 000000000000000 getting cphf filename to store 00000000 END
c +++++++++++++++++++++++++++++++++++++++++++++++++ START
c ========== complex linear solver iteration =========
c +++++++++++++++++++++++++++++++++++++++++++++++++ START
      do iter = 1, maxiter       
        if (debug) write (6,*)
     &     'calling product from ga_lkain_2cpl_damp'
c Note.- product=rohf_hessv3_cmplx,uhf_hessv3_cmplx
          if (debug1) then
           do ipm=1,ncomp
            if (ga_nodeid().eq.0) then
             write(*,2775) ipm,iter
 2775        format('----g_z-toprod(',i3,',',i3,')-----START')
            endif
            call ga_print(g_z(ipm))
            if (ga_nodeid().eq.0) then
             write(*,2776) ipm,iter
 2776        format('----g_z-toprod(',i3,',',i3,')-----END')
            endif
           enddo ! end-loop-ipm
          endif ! end-if-debug1

        call product(acc, 
     &               g_z,      ! in  : x
     &               g_Az1,    ! out : product A x  
     &               nsub,
     &               omega,    ! in  :
     &               limag,    ! in  :
     &               lifetime, ! in  : =.true. -> x is complex
     &               gamwidth, ! in  :
     &               ncomp,    ! in  : nr. components
     &               iter)

        if (debug) write (6,*)
     &     'returning product from ga_lkain_2cpl_damp'

          p1=nsub+1
          p2=nsub+nvec
          m1=1
          m2=n      
          do ipm = 1,ncomp         
c         g_r will be the quantity -Ax + b, i.e. if the equation system
c         Ax = b is solved then this vector will be zero
c
c         During the first cycle,
c         g_Ax is calculated from the initial guess for which the
c         preconditioner has already been applied (to be more clear:
c         we have divided the perturbation matrix elements by orbital
c         energy denominators, including the frequency term, 
c         and assigned real and imaginary parts accordingly)
c ========= get new (r,r_im) ======== START
           call ga_add_patch( one_cmplx,g_zb(ipm),1 ,n ,1 ,nvec, 
     $                       mone_cmplx,g_Az1    ,m1,m2,p1,p2,
     $                                  g_zr1    ,m1,m2,1 ,nvec)
           m1=m1+n
           m2=m2+n          
c ========= get new (r,r_im) ======== END
          enddo                   ! ipm = 1,ncomp

          if (debug1) then
          if (ga_nodeid().eq.0) then
           write(*,2773) iter
 2773      format('---------g_Az1-aft-prod(',i3,')-----START')
          endif
          call ga_print(g_Az1)
          if (ga_nodeid().eq.0) then
           write(*,2778) iter
 2778      format('---------g_Az1-aft-prod(',i3,')-----END')
          endif
           do ipm=1,ncomp
            if (ga_nodeid().eq.0) then
             write(*,2779) ipm,iter
 2779        format('----g_zb(',i3,',',i3,')-----START')
            endif
            call ga_print(g_zb(ipm))
            if (ga_nodeid().eq.0) then
             write(*,2880) ipm,iter
 2880        format('----g_zb(',i3,',',i3,')-----END')
            endif
           enddo ! end-loop-ipm
            if (ga_nodeid().eq.0) then
             write(*,2782) iter
 2782        format('----g_zr1-chk(',i3,')-----START')
            endif
            call ga_print(g_zr1)
            if (ga_nodeid().eq.0) then
             write(*,2783) iter
 2783        format('----g_zr1-chk(',i3,')-----END')
            endif
          endif ! end-if-debug1        

        call get_precond_rmax_zin(
     &                  rmax,    ! out: max(g_r,g_r_im)   
     &                  g_zr1,   ! in : complex+accumulated g_zr
     &                  g_Az1,   ! in : complex+accumulated g_Az
     &                  nsub,    ! in : pointer to current (g_zr,g_Az)
     &                  precond, ! in : name of preconditioner routine
     &        converge_precond,  ! in : =.true. prec->max
     &                  omega,   ! in : omega
     &                  gamwidth,! in : gamwidth
     &                  ncomp,   ! in : nr. components
     &                  npol,    ! in : nr. polarizations (1 or 2)
     &                  nvir,    ! in : nr. virtual  MOs
     &                  nocc,    ! in : nr. occupied MOs
     &                  n,       ! in : =sum_i (nocc * nvir)(i) i=1,npol
     &                  nvec,    ! in : =3 (x,y,z)
     &                  iter,    ! in : nr. iteration
     &                  debug)   ! in : =.true. -> allow debug printouts

         if (debug1) then
          if (ga_nodeid().eq.0) then
           write(*,1775) iter
 1775      format('---------g_z1-to-file(',i3,')-----START')
          endif
          call ga_print(g_z1)
          if (ga_nodeid().eq.0) then
           write(*,1776) iter
 1776      format('---------g_z1-to-file(',i3,')-----END')
          endif
          if (ga_nodeid().eq.0) then
           write(*,1773) iter
 1773      format('---------g_Az1-to-file(',i3,')-----START')
          endif
          call ga_print(g_Az1)
          if (ga_nodeid().eq.0) then
           write(*,1774) iter
 1774      format('---------g_Az1-to-file(',i3,')-----END')
          endif
         endif ! end-if-debug1

        if ((.not.(flag2readfile .and. iter.eq.1).and.
     &            cphf3write.eq.1) .or.
     &      (cphf3write.eq.2 .and. rmax.lt.tol)) ! store only last (g_z1,g_Az1) block
     &  then
          status=dft_CPHF2_write(
     &           aorespfilename, ! in: filename
     &           n,        ! in: sum_{i=1,npol} nocc(i)*nvirt(i)
     &           ncomp,    ! in: nr. components
     &           nvec,     ! in: nr. of directions = 3
     &           n1,       ! in: =n*ncomp
     &           nsub,     ! in: nsub
     &           nsub_file,! in: last subspace index (nsub+1)= nr of subspaces stored
     &           g_z1,     ! in: history matrix z
     &           g_Az1)    ! in: history matrix Az
          nsub_file=nsub_file+nvec
        endif ! end-if-write-block

c       JEM: Putting rmax into rtdb
        if (.not. rtdb_put(rtdb, 'lkain:rmax', mt_dbl, 1, rmax))
     $    call errquit(
     $           'ga_lkain_2cpl3_damp_cmplx_redmem: rmax put failed',
     $                 1, RTDB_ERR) 


c -------- printout per iteration -------------- START
        if (oprint .and. ga_nodeid().eq.0) then

          if (debug1) then
           write(6,4) iter, nsub+nvec, rmax, util_wallsec()
           call util_flush(6)
    4      format('FA-chk: ', i5, i7, 3x,1p,d9.2,0p,f10.1,5x,i3)
          endif     

          write(6,3) iter, nsub+nvec, rmax, util_wallsec()
          call util_flush(6)
    3     format(' ', i5, i7, 3x,1p,d9.2,0p,f10.1,5x,i3)
        end if
c -------- printout per iteration -------------- END

c       stop iterations if residual is smaller than criterion
        if (rmax .lt. tol) then
          converged = .true.
c ========== g_z1 --> g_z ======== START
c Note.- Extract last (n1,nvec) block from g_z1 and
c        put it in g_z (this is the solution to Ax=b)
         m1=1
         m2=n
         p1=nsub+1
         p2=nsub+nvec
         do ipm=1,ncomp
          call ga_copy_patch('n',g_z1    ,m1,m2,p1,p2, 
     $                           g_z(ipm),1 ,n ,1,nvec)
          m1=m1+n
          m2=m2+n          
         enddo ! end-loop-ipm
          if (debug1) then
           do ipm=1,ncomp
            if (ga_nodeid().eq.0) then
             write(*,2785) ipm,iter
 2785        format('----g_z-SOLUTION(',i3,',',i3,')-----START')
            endif
            call ga_print(g_z(ipm))
            if (ga_nodeid().eq.0) then
             write(*,2786) ipm,iter
 2786        format('----g_z-SOLUTION(',i3,',',i3,')-----END')
            endif
           enddo ! end-loop-ipm
          endif ! end-if-debug1
c ========== g_z1 --> g_z ======== END
          goto 100
        end if

          if (debug1) then
          if (ga_nodeid().eq.0) then
           write(*,2700) iter
 2700      format('---------g_z1-0(',i3,')-----START')
          endif
          call ga_print(g_z1)
          if (ga_nodeid().eq.0) then
           write(*,2701) iter
 2701      format('---------g_z1-0(',i3,')-----END')
          endif
          if (ga_nodeid().eq.0) then
           write(*,2702) iter
 2702      format('---------g_Az1-0(',i3,')-----START')
          endif
          call ga_print(g_Az1)
          if (ga_nodeid().eq.0) then
           write(*,2703) iter
 2703      format('---------g_Az1-0(',i3,')-----END')
          endif
          if (ga_nodeid().eq.0) then
           write(*,2704) iter
 2704      format('---------g_zr1-0(',i3,')-----START')
          endif
          call ga_print(g_zr1)
          if (ga_nodeid().eq.0) then
           write(*,2705) iter
 2705      format('---------g_zr1-0(',i3,')-----END')
          endif
          endif ! end-if-debug1

         nsub = nsub + nvec
         select case(slcKAIN)
         case(1) ! Real solver
          call solve_zlineq_KAIN1(
     &              g_Ax1,    ! in/out: history of Ax1  (real)
     &              g_x1,     ! in/out: history of  x1  (real)
     &              g_xr1,    ! in/out: history of  xr1 (real) residual
     &              g_Az1,    ! in    : history of g_Az
     &              g_z1,     ! in/out: history of g_z adding a block g_z1
     &              g_zr1,    ! in    : g_zr
     &              nsub,     ! in    : subspace length
     &              nvec,     ! in    : increment of subspace
     &              ncomp,    ! in    : nr. components
     &              n,        ! in    : nr. elements per comp.
     &              iter,     ! in    : iteration nr.
     &              checkorth,! in    : =1 display: z1^t*(z1c) (goes to zero is it converges)
     &              debug1)   ! in    : =.true. show debug printouts
         case(2) ! complex-solver
          call solve_zlineq1(
     &              g_Az1,    ! in    : history of g_Az
     &              g_z1,     ! in/out: history of g_z adding a block g_z1
     &              g_zr1,    ! in    : g_zr
     &              nsub,     ! in    : subspace length
     &              nvec,     ! in    : increment of subspace
     &              ncomp,    ! in    : nr. components
     &              n,        ! in    : nr. elements per comp.
     &              iter,     ! in    : iteration nr.
     &              checkorth,! in    : =1 display: z1^t*(z1c) (goes to zero is it converges)
     &              debug1)   ! in    : =.true. show debug printouts
         case(3) ! complex-solver-fixed (fixed KAIN with differences)
          call solve_zlineq_KAIN3(
     &              g_Az1,    ! in    : history of g_Az
     &              g_z1,     ! in/out: history of g_z adding a block g_z1
     &              g_zr1,    ! in    : g_zr
     &              nsub,     ! in    : subspace length
     &              nvec,     ! in    : increment of subspace
     &              ncomp,    ! in    : nr. components
     &              n,        ! in    : nr. elements per comp.
     &              iter,     ! in    : iteration nr.
     &              checkorth,! in    : =1 display: z1^t*(z1c) (goes to zero is it converges)
     &              debug1)   ! in    : =.true. show debug printouts
         case (:0) ! All other values
         case (4:)
         call errquit(
     &       'ga_lkain_2cpl3_damp_cmplx_redmem: slcKAIN not 1,2 or 3',
     &        555, RTDB_ERR)
      end select

c -------- get g_z from g_z1 --- START
      m1=1
      m2=n
      p1=nsub+1
      p2=nsub+nvec
      call ga_inquire(g_z1,type,dim1,dim2)
      do ipm=1,ncomp
       call ga_copy_patch('n',g_z1    ,m1,m2,p1,p2,
     &                        g_z(ipm),1 ,n ,1 ,nvec)
       m1=m1+n
       m2=m2+n          
      enddo ! end-loop-ipm
c -------- get g_z from g_z1 --- START

        if (nsub .eq. maxsub-nvec) then
c       
c       Reduce the subspace as necessary
c       
c ====== left-shifting patch ==== START
c Note.- matrices Ay,y shift to left nvec positions
c        removing leftmost patch of dimension: n4 x nvec
         do isub = nvec+1, maxsub, nvec
          call ga_copy_patch('n',g_Az1,1,n1,isub,isub+nvec-1, 
     $                           g_Az1,1,n1,isub-nvec,isub-1)
          call ga_copy_patch('n',g_z1 ,1,n1,isub,isub+nvec-1, 
     $                           g_z1 ,1,n1,isub-nvec,isub-1)
         enddo ! end-loop-isub
c ====== left-shifting patch ==== END
         nsub = nsub - nvec
        end if                  ! (nsub .eq. maxsub) 

      enddo                     ! iter = 1,maxiter
c +++++++++++++++++++++++++++++++++++++++++++++++++ END
c ========== complex linear solver iteration =========
c +++++++++++++++++++++++++++++++++++++++++++++++++ END

  100 continue                  ! jump here if converged 
  
      if (.not. converged) then
        if (ga_nodeid().eq.0) then
          write (luout,*) 'WARNING: CPKS procedure is NOT converged'
          write (luout,*) '  I will proceed, but check your results!'
        endif
c ======= Still write non-converge (g_z1,g_Az1) ==== START
         if (cphf3write.eq.2) then
          nsub=nsub-nvec ! point to previous sub-space
          status=dft_CPHF2_write(
     &           aorespfilename, ! in: filename
     &           n,              ! in: sum_{i=1,npol} nocc(i)*nvirt(i)
     &           ncomp,          ! in: nr. components
     &           nvec,           ! in: nr. of directions = 3
     &           n1,             ! in: =n*ncomp
     &           nsub,           ! in: nsub
     &           nsub_file,      ! in: last subspace index (nsub+1)= nr of subspaces stored
     &           g_z1,           ! in: history matrix z
     &           g_Az1)          ! in: history matrix Az
         endif
c ======= Still write non-converge (g_z1,g_Az1) ==== END
      endif   
      if (.not. ga_destroy(g_zb1)) call errquit
     &   ('lkain_2cpl3-cmplx: destroy zb1',3, GA_ERR)
      if (.not. ga_destroy(g_Az1)) call errquit
     &   ('lkain_2cpl3-cmplx: destroy Az1',3, GA_ERR)
      if (.not. ga_destroy(g_z1)) call errquit
     &   ('lkain_2cpl3-cmplx: destroy z1',4, GA_ERR)
      if (.not. ga_destroy(g_zr1)) call errquit
     &   ('lkain_2cpl3-cmplx: destroy zr1',6, GA_ERR)
c ++++++ added for solve_zlineq_KAIN1 +++ START
      if (slcKAIN.eq.2) then
       if (.not. ga_destroy(g_Ax1)) call errquit
     &   ('lkain_2cpl3-real: destroy Ax1',3, GA_ERR)
       if (.not. ga_destroy(g_x1)) call errquit
     &   ('lkain_2cpl3-real: destroy x1',4, GA_ERR)
       if (.not. ga_destroy(g_xr1)) call errquit
     &   ('lkain_2cpl3-real: destroy xr1',6, GA_ERR)  
      endif 
c ++++++ added for solve_zlineq_KAIN1 +++ END   
      end
c Auxiliar routine for redmem()
      subroutine solve_zlineq_KAIN1(
     &              g_Ax1,    ! in/out: history of Ax1  (real)
     &              g_x1,     ! in/out: history of  x1  (real)
     &              g_xr1,    ! in/out: history of  xr1 (real) residual
     &              g_Az1,    ! in    : history of g_Az
     &              g_z1,     ! in/out: history of g_z adding a block g_z1
     &              g_zr1,    ! in    : g_zr
     &              nsub,     ! in    : subspace length
     &              nvec,     ! in    : increment of subspace
     &              ncomp,    ! in    : nr. components
     &              n,        ! in    : nr. elements per comp.
     &              iter,     ! in    : iteration nr.
     &              checkorth,! in    : =1 display: z1^t*(z1c) (goes to zero is it converges)
     &              debug1)   ! in    : =.true. show debug printouts
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
        integer g_Az1,g_z1,g_zr1,
     &          g_Ax1,g_x1,g_xr1
        integer nsub,nvec,ncomp,n,
     &          iter,checkorth
        logical debug1
        external solve_xlineq_4redmem,
     &           copy_complx2real_4redmem

      call copy_complx2real_4redmem(
     &                         g_Ax1,    ! ou : history of Az1  (real)
     &                         g_x1,     ! ou : history of  z1  (real)
     &                         g_xr1,    ! ou : history of  zr1 (real) residual
     &                         g_Az1,    ! in : history of Az1  (complex)
     &                         g_z1,     ! in : history of  z1  (complex)
     &                         g_zr1,    ! in : history of  zr1 (complex) residual
     &                         nsub,
     &                         nvec,
     &                         ncomp,
     &                         n)

        if (debug1) then
          if (ga_nodeid().eq.0) then
           write(*,2706) iter
 2706      format('---------g_x1-0(',i3,')-----START')
          endif
          call ga_print(g_x1)
          if (ga_nodeid().eq.0) then
           write(*,2707) iter
 2707      format('---------g_x1-0(',i3,')-----END')
          endif
          if (ga_nodeid().eq.0) then
           write(*,2708) iter
 2708      format('---------g_Ax1-0(',i3,')-----START')
          endif
          call ga_print(g_Ax1)
          if (ga_nodeid().eq.0) then
           write(*,2709) iter
 2709      format('---------g_Ax1-0(',i3,')-----END')
          endif
          if (ga_nodeid().eq.0) then
           write(*,2710) iter
 2710      format('---------g_xr1-0(',i3,')-----START')
          endif
          call ga_print(g_xr1)
          if (ga_nodeid().eq.0) then
           write(*,2711) iter
 2711      format('---------g_xr1-0(',i3,')-----END')
          endif
        endif ! end-if-debug1

c ++++++++++++++++++++++++++++++++++++++
c ++++ Solve real linear system +++START
c ++++++++++++++++++++++++++++++++++++++
         call solve_xlineq_4redmem(
     &              g_Ax1,    ! in    : history of products Ax
     &              g_x1,     ! in    : history of solution  x
     &              g_xr1,    ! in    : current residual  r=Ax-b
     &              g_z1,     ! in/out: history of g_z adding a block g_z1
     &              nsub,     ! in    : subspace length
     &              nvec,     ! in    : increment of subspace
     &              ncomp,    ! in    : nr. components
     &              n,        ! in    : nr. elements per comp.
     &              iter,     ! in    : iteration nr.
     &              debug1)   ! in    : =.true. show debug printouts
c ++++++++++++++++++++++++++++++++++++++
c ++++ Solve real linear system +++END
c ++++++++++++++++++++++++++++++++++++++
      return
      end
      subroutine solve_zlineq_KAIN3(
     &              g_Az1,    ! in    : history of g_Az
     &              g_z1,     ! in/out: history of g_z adding a block g_z1
     &              g_zr1,    ! in    : g_zr
     &              nsub,     ! in    : subspace length
     &              nvec,     ! in    : increment of subspace
     &              ncomp,    ! in    : nr. components
     &              n,        ! in    : nr. elements per comp.
     &              iter,     ! in    : iteration nr.
     &              checkorth,! in    : =1 display: z1^t*(z1c) (goes to zero is it converges)
     &              debug1)   ! in    : =.true. show debug printouts

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
        integer g_Az1,g_z1,g_zr1
        integer nsub,nvec,ncomp,n,
     &          iter,checkorth
        logical debug1
        external solve_zlineq1,
     &           solve_zlineq1_fixed

        if (iter.eq.1) then
          call solve_zlineq1(
     &              g_Az1,    ! in    : history of g_Az
     &              g_z1,     ! in/out: history of g_z adding a block g_z1
     &              g_zr1,    ! in    : g_zr
     &              nsub,     ! in    : subspace length
     &              nvec,     ! in    : increment of subspace
     &              ncomp,    ! in    : nr. components
     &              n,        ! in    : nr. elements per comp.
     &              iter,     ! in    : iteration nr.
     &              checkorth,! in    : =1 display: z1^t*(z1c) (goes to zero is it converges)
     &              debug1)   ! in    : =.true. show debug printouts
        else ! else-if-iter 
          call solve_zlineq1_fixed(
     &              g_Az1,    ! in    : history of g_Az
     &              g_z1,     ! in/out: history of g_z adding a block g_z1
     &              g_zr1,    ! in    : g_zr
     &              nsub,     ! in    : subspace length
     &              nvec,     ! in    : increment of subspace
     &              ncomp,    ! in    : nr. components
     &              n,        ! in    : nr. elements per comp.
     &              iter,     ! in    : iteration nr.
     &              checkorth,! in    : =1 display: z1^t*(z1c) (goes to zero is it converges)
     &              debug1)   ! in    : =.true. show debug printouts
        endif ! end-if-iter

      return 
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++
c ------- copy_complx2real_4redmem ----------- START 
c ++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine solve_xlineq_4redmem(
     &              g_Ax1,    ! in    : history of products Ax
     &              g_x1,     ! in    : history of solution  x
     &              g_xr1,    ! in    : current residual  r=Ax-b
     &              g_z1,     ! in/out: history of g_z adding a block g_z1
     &              nsub,     ! in    : subspace length
     &              nvec,     ! in    : increment of subspace
     &              ncomp,    ! in    : nr. components
     &              n,        ! in    : nr. elements per comp.
     &              iter,     ! in    : iteration nr.
     &              debug1)   ! in    : =.true. show debug printouts
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
      external update_g_z1_4redmem1,
     &         update_g_z1_4redmem2,
     &         ga_svd_solve_seq
      integer iter,nsub,nvec,ncomp,n,n1,n4
      integer g_Ax1,g_x1,g_xr1,g_z1,
     &        g_aa,g_bb,g_cc
      logical debug1

        n1=n*ncomp
        n4=2*n1 ! for re+im with two components
        if (.not. ga_create(MT_DBL, nsub, nsub,
     &     'lkain_2cpl3_damp: A', 0, 0, g_aa))
     $     call errquit('lkain: allocating g_a?', nsub, GA_ERR)
        if (.not. ga_create(MT_DBL, nsub, nvec,
     &     'lkain_2cpl3_damp: B', 0, 0,g_bb))
     $     call errquit('lkain: allocating g_bb?', nsub, GA_ERR)
        if (.not. ga_create(MT_DBL, nsub, nvec,
     &     'lkain_2cpl3_damp: C', 0, 0, g_cc))
     $     call errquit('lkain: allocating g_c?', nsub, GA_ERR)
        call ga_zero(g_aa)
        call ga_zero(g_bb)
        call ga_zero(g_cc)
        call ga_dgemm('t','n',nsub,nsub,n4,1.0d0,
     &                g_x1,g_Ax1,0.0d0,g_aa)
        call ga_dgemm('t','n',nsub,nvec,n4,1.0d0,
     &                g_x1,g_xr1,0.0d0,g_bb)

        if (debug1) then
         if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_a(',iter,')-----START' 
         call ga_print(g_aa)
         if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_a(',iter,')-----END'
         if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_b(',iter,')-----START' 
         call ga_print(g_bb)
         if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_b(',iter,')-----END'
        endif ! end-if-debug1

         call ga_svd_solve_seq(g_aa,g_bb,g_cc,1d-14)

        if (debug1) then
         if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_c-old(',iter,')-----START' 
         call ga_print(g_cc)
         if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_c-old(',iter,')-----END'
        endif ! end-if-debug1
c       
c       Form and add the correction, in parts, onto the solution      
c FA: Step 5:
        if (debug1) then
         if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_r2-BEF(',iter,')-----START' 
         call ga_print(g_xr1)
         if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_r2-BEF(',iter,')-----END'
        endif ! end-if-debug1

        call ga_dgemm('n','n',n4,nvec,nsub,-1.0d0,
     &                g_Ax1,g_cc,1.0d0,g_xr1)

        if (debug1) then
         if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_r2-AFT(',iter,')-----START' 
         call ga_print(g_xr1)
         if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_r2-AFT(',iter,')-----END'
        endif ! end-if-debug1
c       
c       copy components of g_r2 into g_r before adding g_r to  g_x
      call update_g_z1_4redmem1(
     &                         g_z1,    ! ou : old solution to update
     &                         g_xr1,   ! in : update
     &                         nsub,
     &                         nvec,
     &                         ncomp,
     &                         n)
         if (debug1) then
           if (ga_nodeid().eq.0) then
           write(*,4700) iter
 4700      format('---------g_z1-1(',i3,')-----START')
           endif
           call ga_print(g_z1)
           if (ga_nodeid().eq.0) then
           write(*,4701) iter
 4701      format('---------g_z1-1(',i3,')-----END')
           endif
         endif ! end-if-debug1
c FA: Step 8:
        call ga_dgemm('n','n',n4,nvec,nsub,1.0d0,
     &                g_x1,g_cc,0.0d0,g_xr1)

        if (debug1) then
         if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_y c(',iter,')-----START' 
         call ga_print(g_xr1)
         if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_y c(',iter,')-----END'
        endif ! end-if-debug1

c       copy components of g_r2 into g_r before adding g_r to  g_x
      call update_g_z1_4redmem2(
     &                         g_z1,    ! ou : old solution to update
     &                         g_xr1,   ! in : update
     &                         nsub,
     &                         nvec,
     &                         ncomp,
     &                         n)

         if (debug1) then
          if (ga_nodeid().eq.0) then
           write(*,4702) iter
 4702      format('---------g_z1-2(',i3,')-----START')
          endif
          call ga_print(g_z1)
          if (ga_nodeid().eq.0) then
           write(*,4703) iter
 4703      format('---------g_z1-2(',i3,')-----END')
          endif 
         endif ! end-if-debug1
     
        if (.not. ga_destroy(g_aa)) call errquit
     &     ('lkain_2cpl: a',0, GA_ERR)
        if (.not. ga_destroy(g_bb)) call errquit
     &     ('lkain_2cpl: b',0, GA_ERR)
        if (.not. ga_destroy(g_cc)) call errquit
     &     ('lkain_2cpl: c',0, GA_ERR)
      return
      end
      subroutine copy_complx2real_4redmem(
     &                         g_Ax1,    ! ou : history of Az1  (real)
     &                         g_x1,     ! ou : history of  z1  (real)
     &                         g_xr1,    ! ou : history of  zr1 (real) residual
     &                         g_Az1,    ! in : history of Az1  (complex)
     &                         g_z1,     ! in : history of  z1  (complex)
     &                         g_zr1,    ! in : history of  zr1 (complex) residual
     &                         nsub,
     &                         nvec,
     &                         ncomp,
     &                         n)
c
c Author : Fredy W. Aquino, Northwestern University
c Purpose: Translate history matrices from complex to real
c          so that I can test old solver for c-KAIN coefficients
c          using real variables
c          g_Az1 --> g_Ax1 (history of products  Az)
c          g_z1  --> g_x1  (history of solutions z )
c          g_zr1 --> g_xr1 (residual for current iteration)
c          dimension(g_Ax1) = 2 dimension(g_Az1)
c          because structure of g_Ax1 = [ re im ...]
c          similarly for g_x1,g_xr1
c Date   : 03-10-14

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
       integer n1,n,ncomp,nvec,nsub,
     &         idat,idat1,
     &         ivec,ivec1,shift
       integer l_z,k_z
       integer g_Az1,g_z1,g_zr1,
     &         g_Ax1,g_x1,g_xr1
       double precision  val_re,val_im

       shift=nsub-nvec
       n1=n*ncomp

      if (.not.MA_Push_Get(mt_dcpl,n1,'copy_complx2real_4redmem: l_z',
     &                     l_z,k_z))
     &  call errquit('copy_complx2real_4redmem: cannot allocate l_z',
     &               n1, MA_ERR)

         do ivec=1,nvec
          ivec1=shift+ivec
          call ga_get(g_Az1,1,n1,ivec1,ivec1,dcpl_mb(k_z),1)
          do idat=1,n1
           idat1=n1+idat
           val_re=dreal(dcpl_mb(k_z+idat-1))
           val_im=dimag(dcpl_mb(k_z+idat-1))
           call ga_put(g_Ax1,idat ,idat ,ivec1,ivec1,val_re,1)
           call ga_put(g_Ax1,idat1,idat1,ivec1,ivec1,val_im,1)
          enddo ! end-loop-idat
         enddo ! end-loop-ivec
         do ivec=1,nvec
          ivec1=shift+ivec
          call ga_get(g_z1,1,n1,ivec1,ivec1,dcpl_mb(k_z),1)
          do idat=1,n1
           idat1=n1+idat
           val_re=dreal(dcpl_mb(k_z+idat-1))
           val_im=dimag(dcpl_mb(k_z+idat-1))
           call ga_put(g_x1,idat ,idat ,ivec1,ivec1,val_re,1)
           call ga_put(g_x1,idat1,idat1,ivec1,ivec1,val_im,1)
          enddo ! end-loop-idat
         enddo ! end-loop-ivec
         call ga_zero(g_xr1)
         do ivec=1,nvec
           call ga_get(g_zr1,1,n1,ivec,ivec,dcpl_mb(k_z),1)
          do idat=1,n1
           idat1=n1+idat
           val_re=dreal(dcpl_mb(k_z+idat-1))
           val_im=dimag(dcpl_mb(k_z+idat-1))
           call ga_put(g_xr1,idat ,idat ,ivec,ivec,val_re,1)
           call ga_put(g_xr1,idat1,idat1,ivec,ivec,val_im,1)
          enddo ! end-loop-idat
         enddo ! end-loop-ivec
      if (.not.ma_pop_stack(l_z))
     $  call errquit('copy_complx2real_4redmem: pop problem with l_z',
     &               555,MA_ERR)
      return
      end
      subroutine update_g_z1_4redmem1(
     &                         g_z1,    ! ou : old solution to update
     &                         g_xr1,   ! in : update
     &                         nsub,
     &                         nvec,
     &                         ncomp,
     &                         n)
c
c Author : Fredy W. Aquino, Northwestern University
c Purpose: Update g_z1 (history vector of solutions by adding a 
c          complex block n1 x nvec
c Date   : 03-10-14

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
       integer n2,n1,n,ncomp,nvec,nsub,
     &         idat,
     &         ivec,shift
       integer l_x,k_x,p1,p2,q1,q2
       integer g_z1,g_xr1,g_z2
       double complex  val_cmplx,one_cmplx

       one_cmplx =dcmplx( 1.0d0,0.0d0)
       shift=nsub-nvec
       n2=2*n*ncomp
       n1=n*ncomp
      if (.not. ga_create(MT_DCPL,n1,nvec, 'update_g_z1_4redmem1: z2', 
     $   0, 0, g_z2))
     $   call errquit('update_g_z1_4redmem1: failed alloc subspace-z1',
     &               nvec,GA_ERR)
      if (.not.MA_Push_Get(mt_dbl,n2,'update_g_z1_4redmem1: l_x',
     &                     l_x,k_x))
     &  call errquit('update_g_z1_4redmem1: cannot allocate l_x',
     &               n2, MA_ERR)
         call ga_zero(g_z2)
         do ivec=1,nvec
           call ga_get(g_xr1,1,n2,ivec,ivec,dbl_mb(k_x),1)
          do idat=1,n1
           val_cmplx=dcmplx(dbl_mb(k_x+idat-1),
     &                      dbl_mb(k_x+idat-1+n1))
           call ga_put(g_z2,idat ,idat ,ivec,ivec,val_cmplx,1)
          enddo ! end-loop-idat
         enddo ! end-loop-ivec
         p1=nsub-nvec+1
         p2=nsub-nvec+nvec
         q1=p1+nvec
         q2=p2+nvec
         call ga_add_patch(one_cmplx,g_z2,1,n1,1,nvec, 
     $                     one_cmplx,g_z1,1,n1,p1,p2,
     $                               g_z1,1,n1,q1,q2)
      if (.not.ma_pop_stack(l_x))
     $  call errquit('update_g_z1_4redmem1: pop problem with l_x',
     &               555,MA_ERR)
      if (.not. ga_destroy(g_z2)) call errquit
     &   ('update_g_z1_4redmem1: destroy z2',3, GA_ERR) 
      return
      end
      subroutine update_g_z1_4redmem2(
     &                         g_z1,    ! ou : old solution to update
     &                         g_xr1,   ! in : update
     &                         nsub,
     &                         nvec,
     &                         ncomp,
     &                         n)
c
c Author : Fredy W. Aquino, Northwestern University
c Purpose: Update g_z1 (history vector of solutions by adding a 
c          complex block n1 x nvec
c Date   : 03-10-14

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
       integer n2,n1,n,ncomp,nvec,nsub,
     &         idat,
     &         ivec,shift
       integer l_x,k_x,p1,p2,q1,q2
       integer g_z1,g_xr1,g_z2
       double complex  val_cmplx,one_cmplx

       one_cmplx =dcmplx( 1.0d0,0.0d0)
       shift=nsub-nvec
       n2=2*n*ncomp
       n1=n*ncomp
      if (.not. ga_create(MT_DCPL,n1,nvec, 'update_g_z1_4redmem2: z2', 
     $   0, 0, g_z2))
     $   call errquit('update_g_z1_4redmem2: failed alloc subspace-z1',
     &               nvec,GA_ERR)
      if (.not.MA_Push_Get(mt_dbl,n2,'update_g_z1_4redmem2: l_x',
     &                     l_x,k_x))
     &  call errquit('update_g_z1_4redmem2: cannot allocate l_x',
     &               n2, MA_ERR)
         call ga_zero(g_z2)
         do ivec=1,nvec
           call ga_get(g_xr1,1,n2,ivec,ivec,dbl_mb(k_x),1)
          do idat=1,n1
           val_cmplx=dcmplx(dbl_mb(k_x+idat-1),
     &                      dbl_mb(k_x+idat-1+n1))
           call ga_put(g_z2,idat ,idat ,ivec,ivec,val_cmplx,1)
          enddo ! end-loop-idat
         enddo ! end-loop-ivec
         p1=nsub-nvec+1
         p2=nsub-nvec+nvec
         q1=p1+nvec
         q2=p2+nvec
         call ga_add_patch(one_cmplx,g_z2,1,n1,1,nvec, 
     $                     one_cmplx,g_z1,1,n1,q1,q2,
     $                               g_z1,1,n1,q1,q2)
      if (.not.ma_pop_stack(l_x))
     $  call errquit('update_g_z1_4redmem2: pop problem with l_x',
     &               555,MA_ERR)
      if (.not. ga_destroy(g_z2)) call errquit
     &   ('update_g_z1_4redmem2: destroy z2',3, GA_ERR) 
      return
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++
c ------- copy_complx2real_4redmem ----------- START 
c ++++++++++++++++++++++++++++++++++++++++++++++++++
c Note.- Differences bet XXredmem and XXredmem1 routines:
c        In ga_lkain_2cpl3_damp_cmplx_redmem 
c        the product routine is: uhf_hessv3_cmplx
c        which uses as output g_Az1 (history matrix of Az products)
c        In ga_lkain_2cpl3_damp_cmplx_redmem1
c        the product routine is: uhf_hessv3_cmplx1
c        which uses as output g_Az (Az product from ith iteration)

      subroutine ga_lkain_2cpl3_damp_cmplx_redmem1(
     &                   rtdb,
     &                   g_z,     ! in/out: solution
     &                   g_zb,    ! in    : b (of Ax=b)  
     &                   product, ! in    : routine to compute Az
     &                   precond, ! in    : routine to do energy scaling Az,r
     $                   tol,
     &                   mmaxsub,
     &                   maxiter, 
     &                   odiff,oprint, 
     &                   omega, 
     &                   limag,
     &                   lifetime,  ! damp means complex, it is redundant
     &                   gamwidth, 
     &                   ncomp,    ! ncomp=2 (+/-)
     &                   npol,
     &                   nvir,
     &                   nocc)
c
c Purpose: Getting g_z by solving recursively a complex linear equation
c          and reducing more memory cost.
c          --> Modified from ga_lkain_2cpl3()
c Author : Fredy W. Aquino, Northwestern University
c Date   : 03-15-12

      implicit none
#include "errquit.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
#include "util.fh"
#include "stdio.fh"
#include "inp.fh"
#include "rtdb.fh"
      integer ncomp       ! [input] no. of components to treat
      integer g_z(ncomp)  ! [input/output] Initial guess/solution (Re,Im)      
      integer g_zb(ncomp),! [input] b of Ax=b
     &        g_z1,       ! Scratch GA contains history of z  in (n1,mmaxsub) 
     &        g_Az1,      ! Scratch GA contains history of Az in (n1,mmaxsub) 
     &        g_zr1,      ! Scratch GA (r= b-Ax) error of size (n1,nvec)
     &        g_Az        ! Scratch GA (store Az-ith product)
      integer g_zb1
      integer npol,
     &        nvir(npol),nocc(npol)
c Note.- In g_z1,g_Az1 a (n1,nvec) block is added per iteration.
      integer rtdb              ! [input] database handle
      double precision omega    ! [input] coupling parameter
      logical limag             ! [input] imaginary perturbation?
      logical lifetime          ! [input] consider damping or not?
      double precision gamwidth ! [input] damping parameter
      external product          ! [input] product routine
      external precond          ! [input] preconditioner routine
      double precision tol      ! [input] convergence threshold
      integer mmaxsub           ! [input] maximum subspace dimension
      integer maxiter           ! [input] maximum no. of iterations
      logical odiff             ! [input] use differences in product
      logical oprint            ! [input] print flag
c
c     Solves the linear equations A(X)=0 for multiple vectors.
c
c ... jochen:
c     This is a modified version of ga_lkain from file ga_it2.F
c     This version allows to solve a coupled set of equations, i.e.
c     there are two right-hand vectors and two initial guesses and two
c     solutions which are coupled. The coupling is mediated by a
c     parameter omega in the call to the preconditioner 
c     (elsewhere, omega is simply called "frequency")
c
c ... jochen: the above comment is from ga_lkain_2cpl3. This here is
c     a modified version of that routine and takes care of a real and an
c     imaginary part for each frequency component. I.e. now arrays
c     have four components ... 
c
c     note: when called from cphf_solve3, odiff = .false. on input
c
c     call product(acc,g_x, g_Ax)
c     . acc is the accuracy trequired for each element of the product
c     . g_x contains the vectors and g_Ax should be filled
c     .     with the product vectors.  The no. of vectors (columns) in
c     . g_x might differ from the no. of vectors input to ga_lkain().
c
c     call precond(g_x,shift)
c     . apply preconditioning directly to the vectors in g_x with the
c     . coupling parameter omega
c
c     On input g_x should contain an initial guess.  It returns the
c     solution.
c
c     maxsub should be at least 3*nvec and can be beneficially increased
c     to about 10*nvec.
c
c     Needs to be extended to store the sub-space vectors out-of-core
c     at least while the product() routine is being executed.
c       
c ... jochen: here in the iteration loops we keep track
c       of two components of the solution vector, ipm = 1 and 2
c       (ipm stands for + (plus) and - (minus) components)
      integer iter,n,n1,
     &        nvec, nsub, isub, 
     &        type, maxsub, ipm,
     &        nsub_file

c ... jochen: for convenience, now most arrays have two components.
c     that might be changed later if memory becomes an issue
      double precision rmax, acc
      logical converged, odebug, debug, 
     &        converge_precond, debug1
      double complex val_cmplx,num
      double precision ac
      integer p1,p2,m1,m2,dim1,dim2,nblock
      double complex one_cmplx,mone_cmplx,zero_cmplx
      logical dft_CPHF2_read,
     &        dft_CPHF2_write,
     &        dft_CPHF2_read2fix
      real ran1
      integer status_gasvd,idum
      double precision factor_x
      external solve_zlineq1,conv2reim_rhs,
     &         get_precond_rmax_zin,
     &         dft_CPHF2_read,
     &         dft_CPHF2_write,
     &         dft_CPHF2_read2fix
      logical status,flag2readfile
      integer index4cphf,checkorth,cphf3write,
     &        csub
      character*255 aorespfilename
      character*(*) lbl_cphfaoresp
      character*255 lbl_cphfaoresp1
      integer iimoderaman,iiistepraman 
      logical ramanspc                 
      parameter(lbl_cphfaoresp='aoresp_fiao_f') 
      one_cmplx =dcmplx( 1.0d0,0.0d0)
      mone_cmplx=dcmplx(-1.0d0,0.0d0)
      zero_cmplx=dcmplx( 0.0d0,0.0d0)
      if (.not. rtdb_get(rtdb, 'cphf:cphf3write',mt_int,1,
     &                   cphf3write))
     &  cphf3write = 0 ! assigns 0 if unsuccessfull read from rtdb
      if (.not. rtdb_get(rtdb, 'cphf:checkorth',mt_int,1,
     &                   checkorth))
     &  checkorth = 0  ! assigns 0 if unsuccessfull read from rtdb
      call ga_inquire(g_z(1),type,n,nvec) ! get (n,nvec) n=sum(nocc*nvirt(i) i=1,npol)
      n1=ncomp*n
      maxsub = mmaxsub          ! So don't modify input scalar arg
      if (maxsub .lt. 3*nvec) maxsub = 3*nvec
      maxsub = (maxsub/nvec)*nvec 
      debug1=.false.
      debug = (.false. .and. ga_nodeid().eq.0) ! for code development

c     check input key if we should check for convergence
c     after the preconditioner has been applied to the residual
      if (.not. rtdb_get(rtdb, 'aoresponse:precond',    mt_log, 1,
     &                            converge_precond))
     &  converge_precond = .false.
      
      if (debug) write (6,*) 'ga_lkain_2cpl_damp omega =',omega
      if (debug) write (6,*) 'ga_lkain_2cpl_damp limag =',limag
      if (debug) write (6,*) 'ga_lkain_2cpl_damp lifetime =',lifetime
      if (debug) write (6,*) 'ga_lkain_2cpl_damp gamwidth =',gamwidth
      if (debug) write (6,*) 'ga_lkain_2cpl_damp ncomp =', ncomp
      if (debug) write (6,*) 'ga_lkain_2cpl3 converge_precond',
     &   converge_precond
c
c     exit if this is the wrong routine to call (lifetime switch
c     must be set)
      if (.not.lifetime) call errquit
     &   ('ga_lkain_2cpl_damp but lifetime=.F.',0,UNKNOWN_ERR)

c     make sure odiff is false (never tested for odiff = .true.)
      if (odiff) call errquit
     &   ('ga_lkain_2cpl_damp odiff=.T.',0,UNKNOWN_ERR)
c     
      odebug = util_print('debug lsolve', print_never) .and. 
     $   ga_nodeid().eq.0
c     
      if (.not. rtdb_get(rtdb, 'cphf:acc',    mt_dbl, 1,
     &                            acc)) acc = 0.0001d0*tol
c ------- create (zre,zim) ---------- START
      if (.not. ga_create(MT_DCPL,n1,nvec, 'lkain_2cpl: Az', 
     $   0, 0, g_Az))
     $   call errquit('lkain: failed alloc subspace Az',nvec,
     &   GA_ERR)
      if (.not. ga_create(MT_DCPL,n1,maxsub, 'lkain_2cpl: z1', 
     $   0, 0, g_z1))
     $   call errquit('lkain: failed alloc subspace-z1',maxsub,
     &   GA_ERR)
      if (.not. ga_create(MT_DCPL,n1,maxsub, 'lkain_2cpl: Az1', 
     $   0, 0, g_Az1))
     $   call errquit('lkain: failed alloc subspace Az1',maxsub,
     &   GA_ERR)
      if (.not. ga_create(MT_DCPL,n1,nvec, 'lkain_2cpl: zr2',
     $   0, 0, g_zr1))
     $   call errquit('lkain_2cpl: failed allocating zr1', nvec,
     &   GA_ERR)  
      if (.not. ga_create(MT_DCPL,n1,nvec, 'lkain_2cpl: zr2',
     $   0, 0, g_zb1))
     $   call errquit('lkain_2cpl: failed allocating zb1', nvec,
     &   GA_ERR)  
      call ga_zero(g_zb1)
      m1=1
      m2=n
      do ipm=1,ncomp
       call ga_copy_patch('n',g_zb(ipm),1 ,n ,1,nvec,
     &                        g_zb1    ,m1,m2,1,nvec)    
       m1=m1+n
       m2=m2+n
      enddo ! end-loop-ipm
      call ga_zero(g_z1)
      call ga_zero(g_Az1)
      call ga_zero(g_zr1)
      call ga_sync()
c ------- create (zre,zim) ---------- END
    
      if (oprint .and. ga_nodeid().eq.0) then
        write(6,1) n1, nvec, maxsub, tol, util_wallsec()
    1   format(//,'Iterative solution of linear equations',/,
     $     '  No. of variables', i9,/,
     $     '  No. of equations', i9,/,
     $     '  Maximum subspace', i9,/,
     $     '       Convergence', 1p,d9.1,/,
     $     '        Start time', 0p,f9.1,/)
        call util_flush(6)
      end if
     
      if (oprint .and. ga_nodeid().eq.0) then
        write(6,2)
        call util_flush(6)
    2   format(/
     $     '   iter   nsub   residual    time ',/,
     $     '   ----  ------  --------  --------- ')
      end if
c     
      nsub = 0
      converged = .false.
c     
c     ---------------------
c     start interation loop
c     ---------------------
c  
c 000000000000000 getting cphf filename to store 00000000 START
      if (.not. rtdb_get(rtdb,'cphf3-aores:guess1',
     & mt_int,1,index4cphf))  index4cphf = 0

      ramanspc=.false.
      status=rtdb_get(rtdb,'raman:aores0',mt_log,1,ramanspc)
      if (ramanspc) then
       if (.not. rtdb_get(rtdb,'raman:aores1', 
     &          mt_int, 1,iimoderaman)) call
     $     errquit('ga_lkain_2cpl3_redmem1: failed to read iimoderaman', 
     &             0, RTDB_ERR)      
       if (.not. rtdb_get(rtdb,'raman:aores2', 
     &          mt_int, 1,iiistepraman)) call
     $   errquit('ga_lkain_2cpl3_redmem1: failed to read iiistepraman', 
     &             0, RTDB_ERR)    
       write(lbl_cphfaoresp1,'(a13,i1,"_",i4.4,"-",i1)') 
     &                      lbl_cphfaoresp,index4cphf,
     &                      iimoderaman,iiistepraman                    
      else
       write(lbl_cphfaoresp1,'(a13,i1)') lbl_cphfaoresp,index4cphf    
      endif
       call util_file_name(lbl_cphfaoresp1,
     &                     .false.,.false.,aorespfilename)
      nsub_file=0 ! reset value in all nodes
      nsub=0      ! reset value in all nodes
      flag2readfile=.false.
      if (.not. dft_CPHF2_read(
     &           aorespfilename, ! in: filename
     &           n,        ! in: sum_{i=1,npol} nocc(i)*nvirt(i)
     &           ncomp,    ! in: nr. components
     &           nvec,     ! in: nr. of directions = 3
     &           n1,       ! in: =n*ncomp
     &           nsub,     ! ou: last subspace index (nsub+1)= nr of subspaces stored
     &           nsub_file,! ou: last subspace read from file
     &           maxsub,   ! in: maximum subspace
     &           g_z1,     ! ou: history matrix z
     &           g_Az1))   ! ou: history matrix Az
     &  then
c          if (ga_nodeid().eq.0)
c     &     write(*,1999) aorespfilename(1:inp_strlen(aorespfilename))
c 1999      format('File ',a,
c     &          ' does not exist, proceed to generate (z1,Az1)')
c ------ g_z0 --> g_z1 ----- START
c Copying initial guess
         nsub_file=0
         nsub=0
         m1=1
         m2=n
         p1=nsub+1
         p2=nsub+nvec
         do ipm=1,ncomp

          if (debug1) then
           if (ga_nodeid().eq.0) then
           write(*,2770) ipm
 2770      format('---------g_z-guess(',i3,')-----START')
           endif
           call ga_print(g_z(ipm))
           if (ga_nodeid().eq.0) then
           write(*,2771) ipm
 2771      format('---------g_z-guess(',i3,')-----END')
           endif
          endif ! end-if-debug1

          call ga_copy_patch('n',g_z(ipm),1 ,n ,1 ,nvec, 
     $                           g_z1    ,m1,m2,p1,p2)
          m1=m1+n
          m2=m2+n          
         enddo ! end-loop-ipm
c ------ g_z0 --> g_z1 ----- END
      else
c Note.- After reading (g_z1,g_Az1) I need to use precond routine
c        which will do energy scaling.
         call ga_sync()
c Note.- I need to propagate nsub
         call ga_igop(6,nsub_file,1,'+') ! node0 nsub ne 0, eq 0 every other node    
         call ga_igop(6,nsub,1,'+')      ! node0 nsub ne 0, eq 0 every other node    
         flag2readfile=.true.
         nblock=nsub/3+1
         m1=1
         m2=n
         p1=nsub+1
         p2=nsub+nvec
         do ipm=1,ncomp
          call ga_copy_patch('n',g_z1    ,m1,m2,p1,p2,
     &                           g_z(ipm),1 ,n ,1 ,nvec)
     
          if (debug1) then
           if (ga_nodeid().eq.0) then
            write(*,2790) ipm
 2790       format('---------g_z-guess(',i3,')-----START')
           endif
           call ga_print(g_z(ipm))
           if (ga_nodeid().eq.0) then
            write(*,2791) ipm
 2791       format('---------g_z-guess(',i3,')-----END')
            endif
          endif ! end-if-debug1
          m1=m1+n
          m2=m2+n          
         enddo ! end-loop-ipm
         if (debug1) then
          if (ga_nodeid().eq.0)
     &      write(*,*) '---------g_z1-read-from-file-----START'
          call ga_print(g_z1)
          if (ga_nodeid().eq.0)
     &      write(*,*) '---------g_z1-read-from-file-----END'
          if (ga_nodeid().eq.0)
     &      write(*,*) '---------g_Az1-read-from-file-----START'
          call ga_print(g_Az1)
          if (ga_nodeid().eq.0)
     &      write(*,*) '---------g_Az1-read-from-file-----END'
         endif ! end-if-debug1
      endif  
c 000000000000000 getting cphf filename to store 00000000 END
c +++++++++++++++++++++++++++++++++++++++++++++++++ START
c ========== complex linear solver iteration =========
c +++++++++++++++++++++++++++++++++++++++++++++++++ START

      do iter = 1, maxiter
        if (debug) write (6,*)
     &     'calling product from ga_lkain_2cpl_damp'
c Note.- product=rohf_hessv3_cmplx,uhf_hessv3_cmplx

          if (debug1) then
           do ipm=1,ncomp
            if (ga_nodeid().eq.0) then
             write(*,2775) ipm,iter
 2775        format('----g_z-toprod(',i3,',',i3,')-----START')
            endif
            call ga_print(g_z(ipm))
            if (ga_nodeid().eq.0) then
             write(*,2776) ipm,iter
 2776        format('----g_z-toprod(',i3,',',i3,')-----END')
            endif
           enddo ! end-loop-ipm
          endif ! end-if-debug1

        call ga_zero(g_Az) ! reset g_Az

        call product(acc, 
     &               g_z,      ! in  : x
     &               g_Az,     ! out : product A x  
     &               omega,    ! in  :
     &               limag,    ! in  :
     &               lifetime, ! in  : =.true. -> x is complex
     &               gamwidth, ! in  :
     &               ncomp)    ! in  : nr. components

        if (debug) write (6,*)
     &     'returning product from ga_lkain_2cpl_damp'

          p1=nsub+1
          p2=nsub+nvec
          m1=1
          m2=n      
c --------- copy g_Az --> g_Az1 ------- START
          call ga_copy_patch('n',g_Az ,1,n1,1,nvec, 
     &                           g_Az1,1,n1,p1,p2)   
c --------- copy g_Az --> g_Az1 ------- END   
          do ipm = 1,ncomp         
c         g_r will be the quantity -Ax + b, i.e. if the equation system
c         Ax = b is solved then this vector will be zero
c
c         During the first cycle,
c         g_Ax is calculated from the initial guess for which the
c         preconditioner has already been applied (to be more clear:
c         we have divided the perturbation matrix elements by orbital
c         energy denominators, including the frequency term, 
c         and assigned real and imaginary parts accordingly)
c ========= get new (r,r_im) ======== START
           call ga_add_patch( one_cmplx,g_zb(ipm),1 ,n ,1 ,nvec, 
     $                       mone_cmplx,g_Az1    ,m1,m2,p1,p2,
     $                                  g_zr1    ,m1,m2,1 ,nvec)
           m1=m1+n
           m2=m2+n          
c ========= get new (r,r_im) ======== END
          enddo                   ! ipm = 1,ncomp

         if (debug1) then
          if (ga_nodeid().eq.0) then
           write(*,2773) iter
 2773      format('---------g_Az1-aft-prod(',i3,')-----START')
          endif
          call ga_print(g_Az1)
          if (ga_nodeid().eq.0) then
           write(*,2778) iter
 2778      format('---------g_Az1-aft-prod(',i3,')-----END')
          endif
           do ipm=1,ncomp
            if (ga_nodeid().eq.0) then
             write(*,2779) ipm,iter
 2779        format('----g_zb(',i3,',',i3,')-----START')
            endif
            call ga_print(g_zb(ipm))
            if (ga_nodeid().eq.0) then
             write(*,2880) ipm,iter
 2880        format('----g_zb(',i3,',',i3,')-----END')
            endif
           enddo ! end-loop-ipm
            if (ga_nodeid().eq.0) then
             write(*,2782) iter
 2782        format('----g_zr1-chk(',i3,')-----START')
            endif
            call ga_print(g_zr1)
            if (ga_nodeid().eq.0) then
             write(*,2783) iter
 2783        format('----g_zr1-chk(',i3,')-----END')
            endif
         endif ! end-if-debug1        

        call get_precond_rmax_zin(
     &                  rmax,    ! out: max(g_r,g_r_im)   
     &                  g_zr1,   ! in : complex+accumulated g_zr
     &                  g_Az1,   ! in : complex+accumulated g_Az
     &                  nsub,    ! in : pointer to current (g_zr,g_Az)
     &                  precond, ! in : name of preconditioner routine
     &        converge_precond,  ! in : =.true. prec->max
     &                  omega,   ! in : omega
     &                  gamwidth,! in : gamwidth
     &                  ncomp,   ! in : nr. components
     &                  npol,    ! in : nr. polarizations (1 or 2)
     &                  nvir,    ! in : nr. virtual  MOs
     &                  nocc,    ! in : nr. occupied MOs
     &                  n,       ! in : =sum_i (nocc * nvir)(i) i=1,npol
     &                  nvec,    ! in : =3 (x,y,z)
     &                  iter,    ! in : nr. iteration
     &                  debug)   ! in : =.true. -> allow debug printouts

         if (debug1) then
          if (ga_nodeid().eq.0) then
           write(*,1775) iter
 1775      format('---------g_z1-to-file(',i3,')-----START')
          endif
          call ga_print(g_z1)
          if (ga_nodeid().eq.0) then
           write(*,1776) iter
 1776      format('---------g_z1-to-file(',i3,')-----END')
          endif
          if (ga_nodeid().eq.0) then
           write(*,1773) iter
 1773      format('---------g_Az1-to-file(',i3,')-----START')
          endif
          call ga_print(g_Az1)
          if (ga_nodeid().eq.0) then
           write(*,1774) iter
 1774      format('---------g_Az1-to-file(',i3,')-----END')
          endif
         endif ! end-if-debug1

c Note.- ".not.(flag2readfile .and. iter.eq.1) .and. cphf3write.eq.1" means
c        if successfully read data (g_z1,g_Az1) from file
c        then skip iter=1 to avoid storing repeteadly the last block.
        if ((.not.(flag2readfile .and. iter.eq.1).and.
     &            cphf3write.eq.1) .or.
     &      (cphf3write.eq.2 .and. rmax.lt.tol)) ! store only last (g_z1,g_Az1) block
     &  then
          status=dft_CPHF2_write(
     &           aorespfilename, ! in: filename
     &           n,              ! in: sum_{i=1,npol} nocc(i)*nvirt(i)
     &           ncomp,          ! in: nr. components
     &           nvec,           ! in: nr. of directions = 3
     &           n1,             ! in: =n*ncomp
     &           nsub,           ! in: nsub
     &           nsub_file,      ! in: last subspace index (nsub+1)= nr of subspaces stored
     &           g_z1,           ! in: history matrix z
     &           g_Az1)          ! in: history matrix Az     
          nsub_file=nsub_file+nvec
        endif ! end-if-write-block
c -------- printout per iteration -------------- START
        if (oprint .and. ga_nodeid().eq.0) then
          write(6,3) iter, nsub+nvec, rmax, util_wallsec()
          call util_flush(6)
    3     format(' ', i5, i7, 3x,1p,d9.2,0p,f10.1,5x,i3)
        end if
c -------- printout per iteration -------------- END
c       stop iterations if residual is smaller than criterion 
        if (rmax .lt. tol) then
          converged = .true.
c ========== g_z1 --> g_z ======== START
c Note.- Extract last (n1,nvec) block from g_z1 and
c        put it in g_z (this is the solution to Ax=b)
         m1=1
         m2=n
         p1=nsub+1
         p2=nsub+nvec
         do ipm=1,ncomp
          call ga_copy_patch('n',g_z1    ,m1,m2,p1,p2, 
     $                           g_z(ipm),1 ,n ,1,nvec)
          m1=m1+n
          m2=m2+n          
         enddo ! end-loop-ipm
          if (debug1) then
           do ipm=1,ncomp
            if (ga_nodeid().eq.0) then
             write(*,2785) ipm,iter
 2785        format('----g_z-SOLUTION(',i3,',',i3,')-----START')
            endif
            call ga_print(g_z(ipm))
            if (ga_nodeid().eq.0) then
             write(*,2786) ipm,iter
 2786        format('----g_z-SOLUTION(',i3,',',i3,')-----END')
            endif
           enddo ! end-loop-ipm
          endif ! end-if-debug1
c ========== g_z1 --> g_z ======== END
          goto 100
        end if

        if (debug1) then
          if (ga_nodeid().eq.0) then
           write(*,2700) iter
 2700      format('---------g_z1-0(',i3,')-----START')
          endif
          call ga_print(g_z1)
          if (ga_nodeid().eq.0) then
           write(*,2701) iter
 2701      format('---------g_z1-0(',i3,')-----END')
          endif
          if (ga_nodeid().eq.0) then
           write(*,2702) iter
 2702      format('---------g_Az1-0(',i3,')-----START')
          endif
          call ga_print(g_Az1)
          if (ga_nodeid().eq.0) then
           write(*,2703) iter
 2703      format('---------g_Az1-0(',i3,')-----END')
          endif
          if (ga_nodeid().eq.0) then
           write(*,2704) iter
 2704      format('---------g_zr1-0(',i3,')-----START')
          endif
          call ga_print(g_zr1)
          if (ga_nodeid().eq.0) then
           write(*,2705) iter
 2705      format('---------g_zr1-0(',i3,')-----END')
          endif
        endif ! end-if-debug1

          nsub = nsub + nvec

          call solve_zlineq1( ! Using complex linear solver
     &              g_Az1,    ! in    : history of g_Az
     &              g_z1,     ! in/out: history of g_z adding a block g_z1
     &              g_zr1,    ! in    : g_zr
     &              nsub,     ! in    : subspace length
     &              nvec,     ! in    : increment of subspace
     &              ncomp,    ! in    : nr. components
     &              n,        ! in    : nr. elements per comp.
     &              iter,     ! in    : iteration nr.
     &              checkorth,! in    : =1 display: z1^t*(z1c) (goes to zero as it converges)
     &              debug1)   ! in    : =.true. show debug printouts
c -------- get g_z from g_z1 --- START
        m1=1
        m2=n
        p1=nsub+1
        p2=nsub+nvec
        call ga_inquire(g_z1,type,dim1,dim2)
        do ipm=1,ncomp
         call ga_copy_patch('n',g_z1    ,m1,m2,p1,p2,
     &                          g_z(ipm),1 ,n ,1 ,nvec)
         m1=m1+n
         m2=m2+n          
        enddo ! end-loop-ipm

        if (nsub .eq. maxsub-nvec) then
c       
c       Reduce the subspace as necessary
c       
c ====== FA: left-shifting patch ==== START
c Note.- matrices Ay,y shift to left nvec positions
c        removing leftmost patch of dimension: n4 x nvec
         do isub = nvec+1, maxsub, nvec
          call ga_copy_patch('n',g_Az1,1,n1,isub,isub+nvec-1, 
     $                           g_Az1,1,n1,isub-nvec,isub-1)
          call ga_copy_patch('n',g_z1 ,1,n1,isub,isub+nvec-1, 
     $                           g_z1 ,1,n1,isub-nvec,isub-1)
         enddo ! end-loop-isub
c ====== FA: left-shifting patch ==== END
         nsub = nsub - nvec
        end if                  ! (nsub .eq. maxsub) 

      enddo                     ! iter = 1,maxiter
c +++++++++++++++++++++++++++++++++++++++++++++++++ END
c ========== complex linear solver iteration =========
c +++++++++++++++++++++++++++++++++++++++++++++++++ END

  100 continue                  ! jump here if converged 

      if (.not. converged) then
        if (ga_nodeid().eq.0) then
          write (luout,*) 'WARNING: CPKS procedure is NOT converged'
          write (luout,*) '  I will proceed, but check your results!'
        endif
c ======= Still write non-converge (g_z1,g_Az1) ==== START
         if (cphf3write.eq.2) then
          nsub=nsub-nvec ! point to previous sub-space
          status=dft_CPHF2_write(
     &           aorespfilename, ! in: filename
     &           n,              ! in: sum_{i=1,npol} nocc(i)*nvirt(i)
     &           ncomp,          ! in: nr. components
     &           nvec,           ! in: nr. of directions = 3
     &           n1,             ! in: =n*ncomp
     &           nsub,           ! in: nsub
     &           nsub_file,      ! in: last subspace index (nsub+1)= nr of subspaces stored
     &           g_z1,           ! in: history matrix z
     &           g_Az1)          ! in: history matrix Az
         endif
c ======= Still write non-converge (g_z1,g_Az1) ==== END
      endif 
      if (.not. ga_destroy(g_zb1)) call errquit
     &   ('lkain_2cpl3-cmplx: destroy zb1',3, GA_ERR)
      if (.not. ga_destroy(g_Az1)) call errquit
     &   ('lkain_2cpl3-cmplx: destroy Az1',3, GA_ERR)
      if (.not. ga_destroy(g_z1)) call errquit
     &   ('lkain_2cpl3-cmplx: destroy z1',4, GA_ERR)
      if (.not. ga_destroy(g_zr1)) call errquit
     &   ('lkain_2cpl3-cmplx: destroy zr1',6, GA_ERR)   
      if (.not. ga_destroy(g_Az)) call errquit
     &   ('lkain_2cpl3-cmplx: destroy Az',3, GA_ERR) 
      end
c ========= Reduce memory consumption ============== END
c ========================================================
c ++++++++++++++++++++clean routine++++++++++++++++++++ END

      subroutine conv2complex(g_z,  ! out: = complx(g_xre,g_xim)
     &                        g_xre,! in : real      arr
     &                        g_xim,! in : imaginary arr
     &                        n,    ! in : n    rows
     &                        nvec, ! in : nvec columns
     &                        ncomp)! in : nr. components
c
c    Purpose: Convert into complex array
c             (g_xre,g_xim) --> g_z
c             structure of g_xre(ipm),g_xim(ipm):  (n,nvec)
c             ipm=1,ncomp
c             Usual values: nvec=3 (x,y,z) ncomp=2
c
c    Author: Fredy W. Aquino, Northwestern University
c    Date  : 04-08-12

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
      integer ipm,ivec,idat,
     &        n,nvec,ncomp,
     &        l_zre,k_zre,
     &        l_zim,k_zim   
      integer g_xre(ncomp),
     &        g_xim(ncomp),g_z(ncomp)   
      double complex val_cmplx

      if (.not.MA_Push_Get(mt_dbl,n,'hessv jfacs',l_zre,k_zre))
     &     call errquit('conv2complex: cannot allocate zre',
     &                  n, MA_ERR)
      if (.not.MA_Push_Get(mt_dbl,n,'hessv kfacs',l_zim,k_zim))
     &     call errquit('conv2complex: cannot allocate zim',
     &                  n, MA_ERR)
        do ipm=1,ncomp
         call ga_zero(g_z(ipm))
         do ivec=1,nvec
          call ga_get(g_xre(ipm),1,n,ivec,ivec,dbl_mb(k_zre),n)
          call ga_get(g_xim(ipm),1,n,ivec,ivec,dbl_mb(k_zim),n)
          do idat=1,n
           val_cmplx=dcmplx(dbl_mb(k_zre+idat-1),
     &                      dbl_mb(k_zim+idat-1))
           call ga_put(g_z(ipm),idat,idat,ivec,ivec,val_cmplx,1)
          enddo ! end-loop-idat
         enddo ! end-loop-ivec
        enddo ! end-loop-ipm
      if (.not.ma_pop_stack(l_zim))
     $  call errquit('conv2complex: pop problem with l_zim',
     &               555,MA_ERR)
      if (.not.ma_pop_stack(l_zre))
     $  call errquit('conv2complex: pop problem with l_zre',
     &               555,MA_ERR)
      return
      end

      subroutine conv2complex1(g_z,  ! out: = complx(g_xre,g_xim)
     &                         g_xre,! in : real      arr
     &                         g_xim,! in : imaginary arr
     &                         nsub, ! in  : pointer to block
     &                         nvir, ! in  : nr. virtual  MOs
     &                         nocc, ! in  : nr. occupied MOs
     &                         ipm,  ! in  : =1,2 components indices
     &                         n,    ! in : n    rows
     &                         nvec) ! in : nvec columns
c
c    Purpose: Convert into complex array
c             (g_xre,g_xim) --> g_z
c             structure of g_xre(ipm),g_xim(ipm):  (n,nvec)
c             ipm=1,ncomp
c             Usual values: nvec=3 (x,y,z) ncomp=2
c
c    Author: Fredy W. Aquino, Northwestern University
c    Date  : 04-08-12

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
      integer ipm,ivec,ivec1,idat,i,j,
     &        idat1,b1,b2,j1,nocc,nvir,ioff,ioff1,
     &        n,nvec,ncomp,nsub,
     &        l_zre,k_zre,
     &        l_zim,k_zim   
      integer g_xre,g_xim,g_z   
      double complex val_cmplx

      if (.not.MA_Push_Get(mt_dbl,nvir,'hessv jfacs',l_zre,k_zre))
     &     call errquit('conv2complex: cannot allocate zre',
     &                  nvir, MA_ERR)
      if (.not.MA_Push_Get(mt_dbl,nvir,'hessv kfacs',l_zim,k_zim))
     &     call errquit('conv2complex: cannot allocate zim',
     &                  nvir, MA_ERR)
         b1=nsub+1
         b2=nsub+nvec
         ivec1=1
         ioff1=(ipm-1)*n
         do ivec=b1,b2
          do i = ga_nodeid()+1,nocc,ga_nnodes()
           ioff = (i-1)*nvir + 1     
           call ga_get(g_xre,ioff,ioff+nvir-1,ivec1,ivec1,
     $                 dbl_mb(k_zre),nvir)
           call ga_get(g_xim,ioff,ioff+nvir-1,ivec1,ivec1,
     $                 dbl_mb(k_zim),nvir)
           do j=1,nvir
            j1=ioff1+ioff+j-1
            val_cmplx=dcmplx(dbl_mb(k_zre+j-1),
     &                       dbl_mb(k_zim+j-1))
            call ga_put(g_z,j1,j1,ivec,ivec,val_cmplx,1)
           enddo ! end-loop-j
          enddo ! end-loop-i
          ivec1=ivec1+1
         enddo ! end-loop-ivec
      if (.not.ma_pop_stack(l_zim))
     $  call errquit('conv2complex: pop problem with l_zim',
     &               555,MA_ERR)
      if (.not.ma_pop_stack(l_zre))
     $  call errquit('conv2complex: pop problem with l_zre',
     &               555,MA_ERR)
      return
      end

      subroutine conv2complex1_u(
     &                         g_z,  ! out: = complx(g_xre,g_xim)
     &                         g_xre,! in : real      arr
     &                         g_xim,! in : imaginary arr
     &                         nsub, ! in : pointer to block
     &                         npol, ! in : nr. polarizations
     &                         nvir, ! in : nr. virtual  MOs
     &                         nocc, ! in : nr. occupied MOs
     &                         ipm,  ! in : =1,2 components indices
     &                         n,    ! in : n    rows
     &                         nvec) ! in : nvec columns
c
c    Purpose: Convert into complex array
c             (g_xre,g_xim) --> g_z
c             structure of g_xre(ipm),g_xim(ipm):  (n,nvec)
c             ipm=1,ncomp
c             Usual values: nvec=3 (x,y,z) ncomp=2
c
c    Author: Fredy W. Aquino, Northwestern University
c    Date  : 04-08-12

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
      integer ipm,ivec,ivec1,idat,i,j,
     &        idat1,b1,b2,j1,
     &        ipol,npol,nocc(npol),nvir(npol),
     &        ioff,ioff1,shift,
     &        n,nvec,ncomp,nsub,
     &        l_zre,k_zre,
     &        l_zim,k_zim   
      integer g_xre,g_xim,g_z   
      double complex val_cmplx
         b1=nsub+1
         b2=nsub+nvec
         ivec1=1
         ioff1=(ipm-1)*n ! n=sum_{i=1,npol} (nocc*nvir)(i)
         do ivec=b1,b2
          do ipol=1,npol
           if (.not.MA_Push_Get(mt_dbl,nvir(ipol),
     &              'hessv jfacs',l_zre,k_zre))
     &     call errquit('conv2complex: cannot allocate zre',
     &                  nvir(ipol), MA_ERR)
           if (.not.MA_Push_Get(mt_dbl,nvir(ipol),
     &              'hessv kfacs',l_zim,k_zim))
     &     call errquit('conv2complex: cannot allocate zim',
     &                  nvir(ipol), MA_ERR)
           shift=nocc(1)*nvir(1)*(ipol-1)
           do i = ga_nodeid()+1,nocc(ipol),ga_nnodes()
            ioff = shift+(i-1)*nvir(ipol) + 1     
            call ga_get(g_xre,ioff,ioff+nvir(ipol)-1,ivec1,ivec1,
     $                  dbl_mb(k_zre),nvir(ipol))
            call ga_get(g_xim,ioff,ioff+nvir(ipol)-1,ivec1,ivec1,
     $                  dbl_mb(k_zim),nvir(ipol))
            do j=1,nvir(ipol)
             j1=ioff1+ioff+j-1
             val_cmplx=dcmplx(dbl_mb(k_zre+j-1),
     &                        dbl_mb(k_zim+j-1))
             call ga_put(g_z,j1,j1,ivec,ivec,val_cmplx,1)
            enddo ! end-loop-j
           enddo ! end-loop-i
          if (.not.ma_pop_stack(l_zim))
     $     call errquit('conv2complex: pop problem with l_zim',
     &               555,MA_ERR)
          if (.not.ma_pop_stack(l_zre))
     $     call errquit('conv2complex: pop problem with l_zre',
     &               555,MA_ERR)
          enddo ! end-loop-ipol
          ivec1=ivec1+1
         enddo ! end-loop-ivec

      return
      end

      subroutine conv2complex2(g_z,    ! out: = complx(g_xre,g_xim)
     &                         g_xreim,! in : real      arr
     &                         indrm)  ! in  : =1 -> re =2 -> im
c
c    Purpose: Convert into complex array
c             (g_xre,g_xim) --> g_z
c             structure of g_xre(ipm),g_xim(ipm):  (n,nvec)
c             ipm=1,ncomp
c             Usual values: nvec=3 (x,y,z) ncomp=2
c    This routine does: g_z= g_z + g_xreim
c
c    Author: Fredy W. Aquino, Northwestern University
c    Date  : 04-08-12

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
      integer ipm,indrm,
     &        ivec,idat,
     &        n,nvec,
     &        l_xreim,k_xreim,
     &        l_z,k_z,
     &        nocc,nvir,i,j,ioff   
      integer g_xreim,
     &        g_z,g_a,type
      double complex val_cmplx,one_cmplx
      one_cmplx =dcmplx( 1.0d0,0.0d0)
      if (indrm.ne.1 .and.
     &    indrm.ne.2) then
       call errquit('conv2complex2: indrm ne 1 or 2',
     &                  0,MA_ERR) 
      endif
      if (.not. ga_create(MT_DCPL,1,1,
     &     'conv2complex2: A',0,0,g_a))
     $     call errquit('solve_rlineq: failed allocating g_a', 
     &                  1,GA_ERR)
      call ga_inquire(g_z,type,n,nvec) ! get (n,nvec)
      if (.not.MA_Push_Get(mt_dbl,n,'hessv jfacs',l_xreim,k_xreim))
     &     call errquit('conv2complex2: cannot allocate xreim',
     &                  n, MA_ERR)
      if (.not.MA_Push_Get(mt_dcpl,n,'hessv kfacs',l_z,k_z))
     &     call errquit('conv2complex2: cannot allocate z',
     &                  n, MA_ERR)
      if   (indrm.eq.1) then ! updating only REAL part
       do ivec=1,nvec
        call ga_get(g_xreim,1,n,ivec,ivec,dbl_mb(k_xreim),n)
        call ga_get(g_z    ,1,n,ivec,ivec,dcpl_mb(k_z),n)
        do idat=1,n
         val_cmplx=dcmplx(dbl_mb(k_xreim+idat-1),0.0d0)
         call ga_put(g_a,1,1,1,1,val_cmplx,1)
         call ga_add_patch(one_cmplx,g_z,idat,idat,ivec,ivec,
     &                     one_cmplx,g_a,1   ,1   ,1   ,1   ,
     &                               g_z,idat,idat,ivec,ivec)
        enddo ! end-loop-idat
       enddo ! end-loop-ivec
      else if (indrm.eq.2) then ! updating only IMAG part
       do ivec=1,nvec
        call ga_get(g_xreim,1,n,ivec,ivec,dbl_mb(k_xreim),n)
        call ga_get(g_z    ,1,n,ivec,ivec,dcpl_mb(k_z),n)
          do idat=1,n
           val_cmplx=dcmplx(0.0d0,dbl_mb(k_xreim+idat-1))
           call ga_put(g_a,1,1,1,1,val_cmplx,1)
           call ga_add_patch(one_cmplx,g_z,idat,idat,ivec,ivec,
     &                       one_cmplx,g_a,1   ,1   ,1   ,1   ,
     &                                 g_z,idat,idat,ivec,ivec)
          enddo ! end-loop-idat
       enddo ! end-loop-ivec
      endif ! end-if-indrm
      if (.not.ma_pop_stack(l_z))
     $  call errquit('conv2complex2: pop problem with l_z',
     &               555,MA_ERR)
      if (.not.ma_pop_stack(l_xreim))
     $  call errquit('conv2complex2: pop problem with l_xreim',
     &               555,MA_ERR)
      return
      end

      subroutine conv2complex3(g_z,    ! out: = complx(g_xre,g_xim)
     &                         g_xreim,! in : real      arr
     &                         nvir,   ! in  : nr. virtual  MOs
     &                         nocc,   ! in  : nr. occupied MOs
     &                         indrm)  ! in  : =1 -> re =2 -> im
c
c    Purpose: Convert into complex array
c             (g_xre,g_xim) --> g_z
c             structure of g_xre(ipm),g_xim(ipm):  (n,nvec)
c             ipm=1,ncomp
c             Usual values: nvec=3 (x,y,z) ncomp=2
c    This routine does: g_z= g_xreim (copies either RE or IM)
c
c    Author: Fredy W. Aquino, Northwestern University
c    Date  : 04-08-12

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
      integer ipm,indrm,
     &        ivec,
     &        n,nvec,
     &        l_xreim,k_xreim,
     &        l_z,k_z,
     &        nocc,nvir,i,j,j1,ioff   
      integer g_xreim,
     &        g_z,type
      double precision val_real,val_imag
      double complex val_cmplx,one_cmplx
      one_cmplx =dcmplx( 1.0d0,0.0d0)
      if (indrm.ne.1 .and.
     &    indrm.ne.2) then
       call errquit('conv2complex2: indrm ne 1 or 2',
     &                  0,MA_ERR) 
      endif
      call ga_inquire(g_z,type,n,nvec) ! get (n,nvec)
      if (.not.MA_Push_Get(mt_dbl,nvir,'hessv jfacs',l_xreim,k_xreim))
     &     call errquit('conv2complex3: cannot allocate xreim',
     &                  nvir, MA_ERR)
      if (.not.MA_Push_Get(mt_dcpl,nvir,'hessv kfacs',l_z,k_z))
     &     call errquit('conv2complex3: cannot allocate z',
     &                  nvir, MA_ERR)
      if   (indrm.eq.1) then ! updating only REAL part
       do ivec=1,nvec
        do i = ga_nodeid()+1,nocc,ga_nnodes()
         ioff = (i-1)*nvir + 1     
         call ga_get(g_xreim,ioff,ioff+nvir-1,ivec,ivec,
     $               dbl_mb(k_xreim),nvir)
         call ga_get(g_z,ioff,ioff+nvir-1,ivec,ivec,
     $               dcpl_mb(k_z),nvir)
         do j=1,nvir
          j1=ioff+j-1
          val_imag =dimag(dcpl_mb(k_z+j-1))
          val_cmplx=dcmplx(dbl_mb(k_xreim+j-1),val_imag)
          call ga_put(g_z,j1,j1,ivec,ivec,val_cmplx,1)
         enddo ! end-loop-j
        enddo ! end-loop-i
       enddo ! end-loop-ivec
      else if (indrm.eq.2) then ! updating only IMAG part
       do ivec=1,nvec
        do i = ga_nodeid()+1,nocc,ga_nnodes()
         ioff = (i-1)*nvir + 1     
         call ga_get(g_xreim,ioff,ioff+nvir-1,ivec,ivec,
     $               dbl_mb(k_xreim),nvir)
         call ga_get(g_z,ioff,ioff+nvir-1,ivec,ivec,
     $               dcpl_mb(k_z),nvir)
         do j=1,nvir
          j1=ioff+j-1
          val_real =dreal(dcpl_mb(k_z+j-1))
          val_cmplx=dcmplx(val_real,dbl_mb(k_xreim+j-1))
          call ga_put(g_z,j1,j1,ivec,ivec,val_cmplx,1)
         enddo ! end-loop-j
        enddo ! end-loop-i
       enddo ! end-loop-ivec
      endif ! end-if-indrm
      if (.not.ma_pop_stack(l_z))
     $  call errquit('conv2complex3: pop problem with l_z',
     &               555,MA_ERR)
      if (.not.ma_pop_stack(l_xreim))
     $  call errquit('conv2complex3: pop problem with l_xreim',
     &               555,MA_ERR)
      return
      end

      subroutine conv2complex4(g_z,    ! out: = history matrix complex
     &                         g_xreim,! in : real      arr
     &                         nsub,   ! in  : subblock index
     &                         ipm,    ! in  : = 1,2 to access slctd component
     &                         nvir,   ! in  : nr. virtual  MOs
     &                         nocc,   ! in  : nr. occupied MOs
     &                         indrm)  ! in  : =1 -> re =2 -> im
c
c    Purpose: Convert into complex array
c             (g_xre,g_xim) --> g_z
c             structure of g_xre(ipm),g_xim(ipm):  (n,nvec)
c             ipm=1,ncomp
c             Usual values: nvec=3 (x,y,z) ncomp=2
c    This routine does: g_z= g_xreim (copies either RE or IM)
c
c    Author: Fredy W. Aquino, Northwestern University
c    Date  : 04-08-12
c
c    dim(g_z) = (n1,maxsub)  n1=n*ncomp  maxsub=maxiter*nvec
c               n=nocc*nvirt maxiter=10 (usually) nvec=3 (x,y,z)
c               ncomp=2

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
      integer ipm,indrm,
     &        ivec,ivec1,
     &        n,nvec,
     &        l_xreim,k_xreim,
     &        l_z,k_z,nsub,
     &        nocc,nvir,i,j,j1,
     &        ioff,ioff1,ioff2  
      integer g_xreim,
     &        g_z,type
      double precision val_real,val_imag
      double complex val_cmplx,one_cmplx
      one_cmplx =dcmplx( 1.0d0,0.0d0)
      if (indrm.ne.1 .and.
     &    indrm.ne.2) then
       call errquit('conv2complex2: indrm ne 1 or 2',
     &                  0,MA_ERR) 
      endif
      call ga_inquire(g_xreim,type,n,nvec) ! get (n,nvec)
      if (.not.MA_Push_Get(mt_dbl,nvir,'hessv jfacs',l_xreim,k_xreim))
     &     call errquit('conv2complex3: cannot allocate xreim',
     &                  nvir, MA_ERR)
      if (.not.MA_Push_Get(mt_dcpl,nvir,'hessv kfacs',l_z,k_z))
     &     call errquit('conv2complex3: cannot allocate z',
     &                  nvir, MA_ERR)
      ioff1=(ipm-1)*n
      if   (indrm.eq.1) then ! updating only REAL part
       ivec1=1
       do ivec=nsub+1,nsub+nvec
        do i = ga_nodeid()+1,nocc,ga_nnodes()
         ioff =(i-1)*nvir + 1   
         ioff2=ioff1+ioff  
         call ga_get(g_xreim,ioff,ioff+nvir-1,ivec1,ivec1,
     $               dbl_mb(k_xreim),nvir)
         call ga_get(g_z,ioff2,ioff2+nvir-1,ivec,ivec,
     $               dcpl_mb(k_z),nvir)
         do j=1,nvir
          j1=ioff1+ioff+j-1
          val_imag =dimag(dcpl_mb(k_z+j-1))
          val_cmplx=dcmplx(dbl_mb(k_xreim+j-1),val_imag)
          call ga_put(g_z,j1,j1,ivec,ivec,val_cmplx,1)
         enddo ! end-loop-j
        enddo ! end-loop-i
        ivec1=ivec1+1
       enddo ! end-loop-ivec
      else if (indrm.eq.2) then ! updating only IMAG part
       ivec1=1
       do ivec=nsub+1,nsub+nvec
        do i = ga_nodeid()+1,nocc,ga_nnodes()
         ioff =(i-1)*nvir + 1   
         ioff2=ioff1+ioff  
         call ga_get(g_xreim,ioff,ioff+nvir-1,ivec1,ivec1,
     $               dbl_mb(k_xreim),nvir)
         call ga_get(g_z,ioff2,ioff2+nvir-1,ivec,ivec,
     $               dcpl_mb(k_z),nvir)
         do j=1,nvir
          j1=ioff1+ioff+j-1
          val_real =dreal(dcpl_mb(k_z+j-1))
          val_cmplx=dcmplx(val_real,dbl_mb(k_xreim+j-1))
          call ga_put(g_z,j1,j1,ivec,ivec,val_cmplx,1)
         enddo ! end-loop-j
        enddo ! end-loop-i
        ivec1=ivec1+1
       enddo ! end-loop-ivec
      endif ! end-if-indrm
      if (.not.ma_pop_stack(l_z))
     $  call errquit('conv2complex3: pop problem with l_z',
     &               555,MA_ERR)
      if (.not.ma_pop_stack(l_xreim))
     $  call errquit('conv2complex3: pop problem with l_xreim',
     &               555,MA_ERR)
      return
      end

      subroutine conv2complex4_u(
     &                         g_z,    ! out: = history matrix complex
     &                         g_xreim,! in : real      arr
     &                         shift,  ! in : = 0 spin 1, (nocc*nvir)(1) spin 2
     &                         nsub,   ! in : subblock index
     &                         ipm,    ! in : = 1,2 to access slctd component
     &                         nvir,   ! in : nr. virtual  MOs
     &                         nocc,   ! in : nr. occupied MOs
     &                         indrm)  ! in : =1 -> re =2 -> im
c
c    Purpose: Convert into complex array
c             (g_xre,g_xim) --> g_z
c             structure of g_xre(ipm),g_xim(ipm):  (n,nvec)
c             ipm=1,ncomp
c             Usual values: nvec=3 (x,y,z) ncomp=2
c    This routine does: g_z= g_xreim (copies either RE or IM)
c
c    Author: Fredy W. Aquino, Northwestern University
c    Date  : 04-08-12
c
c    dim(g_z) = (n1,maxsub)  n1=n*ncomp  maxsub=maxiter*nvec
c               n=nocc*nvirt maxiter=10 (usually) nvec=3 (x,y,z)
c               ncomp=2

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
      integer ipm,indrm,
     &        ivec,ivec1,
     &        n,nvec,shift,
     &        l_xreim,k_xreim,
     &        l_z,k_z,nsub,
     &        nocc,nvir,i,j,j1,
     &        ioff,ioff1,ioff2  
      integer g_xreim,
     &        g_z,type
      double precision val_real,val_imag
      double complex val_cmplx,one_cmplx
      one_cmplx =dcmplx( 1.0d0,0.0d0)
      if (indrm.ne.1 .and.
     &    indrm.ne.2) then
       call errquit('conv2complex2: indrm ne 1 or 2',
     &                  0,MA_ERR) 
      endif
      call ga_inquire(g_xreim,type,n,nvec) ! get (n,nvec)
      if (.not.MA_Push_Get(mt_dbl,nvir,'hessv jfacs',l_xreim,k_xreim))
     &     call errquit('conv2complex3: cannot allocate xreim',
     &                  nvir, MA_ERR)
      if (.not.MA_Push_Get(mt_dcpl,nvir,'hessv kfacs',l_z,k_z))
     &     call errquit('conv2complex3: cannot allocate z',
     &                  nvir, MA_ERR)
      ioff1=shift+(ipm-1)*n
      if   (indrm.eq.1) then ! updating only REAL part
       ivec1=1
       do ivec=nsub+1,nsub+nvec
        do i = ga_nodeid()+1,nocc,ga_nnodes()
         ioff =(i-1)*nvir + 1   
         ioff2=ioff1+ioff  
         call ga_get(g_xreim,ioff,ioff+nvir-1,ivec1,ivec1,
     $               dbl_mb(k_xreim),nvir)
         call ga_get(g_z,ioff2,ioff2+nvir-1,ivec,ivec,
     $               dcpl_mb(k_z),nvir)
         do j=1,nvir
          j1=ioff1+ioff+j-1
          val_imag =dimag(dcpl_mb(k_z+j-1))
          val_cmplx=dcmplx(dbl_mb(k_xreim+j-1),val_imag)
          call ga_put(g_z,j1,j1,ivec,ivec,val_cmplx,1)
         enddo ! end-loop-j
        enddo ! end-loop-i
        ivec1=ivec1+1
       enddo ! end-loop-ivec
      else if (indrm.eq.2) then ! updating only IMAG part
       ivec1=1
       do ivec=nsub+1,nsub+nvec
        do i = ga_nodeid()+1,nocc,ga_nnodes()
         ioff =(i-1)*nvir + 1   
         ioff2=ioff1+ioff  
         call ga_get(g_xreim,ioff,ioff+nvir-1,ivec1,ivec1,
     $               dbl_mb(k_xreim),nvir)
         call ga_get(g_z,ioff2,ioff2+nvir-1,ivec,ivec,
     $               dcpl_mb(k_z),nvir)
         do j=1,nvir
          j1=ioff1+ioff+j-1
          val_real =dreal(dcpl_mb(k_z+j-1))
          val_cmplx=dcmplx(val_real,dbl_mb(k_xreim+j-1))
          call ga_put(g_z,j1,j1,ivec,ivec,val_cmplx,1)
         enddo ! end-loop-j
        enddo ! end-loop-i
        ivec1=ivec1+1
       enddo ! end-loop-ivec
      endif ! end-if-indrm
      if (.not.ma_pop_stack(l_z))
     $  call errquit('conv2complex3: pop problem with l_z',
     &               555,MA_ERR)
      if (.not.ma_pop_stack(l_xreim))
     $  call errquit('conv2complex3: pop problem with l_xreim',
     &               555,MA_ERR)
      return
      end

      subroutine conv2complex4_u1(
     &                         g_z,    ! out: = history matrix complex
     &                         g_xreim,! in : real      arr
     &                         nsub,   ! in : subblock index
     &                         ipm,    ! in : = 1,2 to access slctd component
     &                         npol,   ! in : nr. polarizations
     &                         nvir,   ! in : nr. virtual  MOs
     &                         nocc,   ! in : nr. occupied MOs
     &                         indrm)  ! in : =1 -> re =2 -> im
c
c    Purpose: Convert into complex array
c             (g_xre,g_xim) --> g_z
c             structure of g_xre(ipm),g_xim(ipm):  (n,nvec)
c             ipm=1,ncomp
c             Usual values: nvec=3 (x,y,z) ncomp=2
c    This routine does: g_z= g_xreim (copies either RE or IM)
c
c    Author: Fredy W. Aquino, Northwestern University
c    Date  : 04-08-12
c
c    dim(g_z) = (n1,maxsub)  n1=n*ncomp  maxsub=maxiter*nvec
c               n=nocc*nvirt maxiter=10 (usually) nvec=3 (x,y,z)
c               ncomp=2

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
      integer ipm,indrm,
     &        ivec,ivec1,
     &        n,nvec,shift,
     &        l_xreim,k_xreim,
     &        l_z,k_z,nsub,
     &        ipol,npol,nocc(npol),nvir(npol),
     &        i,j,j1,
     &        ioff,ioff1,ioff2  
      integer g_xreim,
     &        g_z,type
      double precision val_real,val_imag
      double complex val_cmplx,one_cmplx
      one_cmplx =dcmplx( 1.0d0,0.0d0)
      if (indrm.ne.1 .and.
     &    indrm.ne.2) then
       call errquit('conv2complex2: indrm ne 1 or 2',
     &                  0,MA_ERR) 
      endif
      call ga_inquire(g_xreim,type,n,nvec) ! get (n,nvec)
 
      ioff1=(ipm-1)*n
      if   (indrm.eq.1) then ! updating only REAL part
       ivec1=1
       do ivec=nsub+1,nsub+nvec
        do ipol=1,npol
         if (.not.MA_Push_Get(mt_dbl,nvir(ipol),'hessv jfacs',
     &       l_xreim,k_xreim))
     &     call errquit('conv2complex3: cannot allocate xreim',
     &                  nvir(ipol), MA_ERR)
         if (.not.MA_Push_Get(mt_dcpl,nvir(ipol),'hessv kfacs',
     &        l_z,k_z))
     &     call errquit('conv2complex3: cannot allocate z',
     &                  nvir(ipol), MA_ERR)
         shift=nocc(1)*nvir(1)*(ipol-1)    
         do i = ga_nodeid()+1,nocc(ipol),ga_nnodes()
          ioff =shift+(i-1)*nvir(ipol) + 1   
          ioff2=ioff1+ioff  
          call ga_get(g_xreim,ioff,ioff+nvir(ipol)-1,ivec1,ivec1,
     $                dbl_mb(k_xreim),nvir(ipol))
          call ga_get(g_z,ioff2,ioff2+nvir(ipol)-1,ivec,ivec,
     $                dcpl_mb(k_z),nvir(ipol))
          do j=1,nvir(ipol)
           j1=ioff2+j-1
           val_imag =dimag(dcpl_mb(k_z+j-1))
           val_cmplx=dcmplx(dbl_mb(k_xreim+j-1),val_imag)
           call ga_put(g_z,j1,j1,ivec,ivec,val_cmplx,1)
          enddo ! end-loop-j
         enddo ! end-loop-i
         if (.not.ma_pop_stack(l_z))
     $   call errquit('conv2complex3: pop problem with l_z',
     &               555,MA_ERR)
         if (.not.ma_pop_stack(l_xreim))
     $   call errquit('conv2complex3: pop problem with l_xreim',
     &               555,MA_ERR)
        enddo ! end-loop-ipol
        ivec1=ivec1+1
       enddo ! end-loop-ivec
      else if (indrm.eq.2) then ! updating only IMAG part
       ivec1=1
       do ivec=nsub+1,nsub+nvec
        do ipol=1,npol
         if (.not.MA_Push_Get(mt_dbl,nvir(ipol),'hessv jfacs',
     &       l_xreim,k_xreim))
     &     call errquit('conv2complex3: cannot allocate xreim',
     &                  nvir(ipol), MA_ERR)
         if (.not.MA_Push_Get(mt_dcpl,nvir(ipol),'hessv kfacs',
     &        l_z,k_z))
     &     call errquit('conv2complex3: cannot allocate z',
     &                  nvir(ipol), MA_ERR)
         shift=nocc(1)*nvir(1)*(ipol-1)    
         do i = ga_nodeid()+1,nocc(ipol),ga_nnodes()
          ioff =shift+(i-1)*nvir(ipol) + 1   
          ioff2=ioff1+ioff  
          call ga_get(g_xreim,ioff,ioff+nvir(ipol)-1,ivec1,ivec1,
     $                dbl_mb(k_xreim),nvir(ipol))
          call ga_get(g_z,ioff2,ioff2+nvir(ipol)-1,ivec,ivec,
     $                dcpl_mb(k_z),nvir(ipol))
          do j=1,nvir(ipol)
           j1=ioff2+j-1
           val_real =dreal(dcpl_mb(k_z+j-1))
           val_cmplx=dcmplx(val_real,dbl_mb(k_xreim+j-1))
           call ga_put(g_z,j1,j1,ivec,ivec,val_cmplx,1)
          enddo ! end-loop-j
         enddo ! end-loop-i
         if (.not.ma_pop_stack(l_z))
     $   call errquit('conv2complex3: pop problem with l_z',
     &               555,MA_ERR)
         if (.not.ma_pop_stack(l_xreim))
     $   call errquit('conv2complex3: pop problem with l_xreim',
     &               555,MA_ERR)
        enddo ! end-loop-ipol
        ivec1=ivec1+1
       enddo ! end-loop-ivec
      endif ! end-if-indrm

      return
      end

      subroutine update_gz_reorim(g_z,    ! out: = complx(g_xre,g_xim)
     &                            g_xreim,! in : real      arr 
     &                            indrm,  ! in : =1 -> re =2 -> im
     &                            scl,    ! in : scaling factor 
     &                            nvir,
     &                            nocc,
     &                            ivec)
c
c    Purpose: Convert into complex array
c             (g_xre,g_xim) --> g_z
c             structure of g_xre,g_xim:  (n,nvec)  n=nvir*nocc
c             ipm=1,ncomp
c             Usual values: nvec=3 (x,y,z) ncomp=2
c
c    Author: Fredy W. Aquino, Northwestern University
c    Date  : 04-08-12
c
c Note.- To be used in rohf_hessv_2e3_opt_cmplx() 
c        located in ddscf/rohf_hessv3.F
c        To mimic,
c          call ga_mat_to_vec(g_tmp1,1,nvir,1,nclosed,
c     $                       g_ax_re(ipm),1,ivec,four,'+')  [ scl=four RDFT]
c        for complex g_z (instead of g_ax_re)
c        g_xreim=g_tmp1

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
      integer ipm,indrm,
     &        ivec,idat,nvec,
     &        n,n1,nvir,nocc,
     &        l_xreim,k_xreim,
     &        l_z,k_z,i,j,j1,ioff  
      integer g_xreim,
     &        g_z,g_a,type
      double precision scl
      double complex val_cmplx,one_cmplx,scl_cmplx
      one_cmplx =dcmplx(1.0d0,0.0d0)
      scl_cmplx =dcmplx(scl,0.0d0)
      if (indrm.ne.1 .and.
     &    indrm.ne.2) then
       call errquit('conv2complex2: indrm ne 1 or 2',
     &                  0,MA_ERR) 
      endif
      if (.not. ga_create(MT_DCPL,1,1,
     &     'conv2complex2: A',0,0,g_a))
     $     call errquit('solve_rlineq: failed allocating g_a', 
     &                  1,GA_ERR)
      call ga_inquire(g_z,type,n,nvec) ! get (n,nvec)
      if (.not.MA_Push_Get(mt_dbl,nvir,'hessv jfacs',l_xreim,k_xreim))
     &     call errquit('conv2complex2: cannot allocate xreim',
     &                  nvir, MA_ERR)
      if (.not.MA_Push_Get(mt_dcpl,nvir,'hessv kfacs',l_z,k_z))
     &     call errquit('conv2complex2: cannot allocate z',
     &                  nvir, MA_ERR)
      if   (indrm.eq.1) then ! updating only REAL part
c ++++++++++++++++++
c NOTE.- Assumming the ordering in dbl_mb(k_xreim) is same as in dcpl_mb(k_z)
c ++++++++++++++++++
       do i = ga_nodeid()+1,nocc,ga_nnodes()
        ioff = (i-1)*nvir + 1   
        call ga_get(g_xreim,1,nvir,i,i,dbl_mb(k_xreim),nvir)  
        call ga_get(g_z,ioff,ioff+nvir-1,ivec,ivec,
     &              dcpl_mb(k_z),nvir)
        do j=1,nvir
          val_cmplx=dcmplx(dbl_mb(k_xreim+j-1),0.0d0)
          call ga_put(g_a,1,1,1,1,val_cmplx,1)
          j1=ioff+j-1
          call ga_add_patch(one_cmplx,g_z,j1,j1,ivec,ivec,
     &                      scl_cmplx,g_a,1 ,1 ,1   ,1   ,
     &                                g_z,j1,j1,ivec,ivec)
        enddo ! end-loop-j
       enddo ! end-loop-i
      else if (indrm.eq.2) then ! updating only IMAG part
       do i = ga_nodeid()+1,nocc,ga_nnodes()
        ioff = (i-1)*nvir + 1   
        call ga_get(g_xreim,1,nvir,i,i,dbl_mb(k_xreim),nvir)  
        call ga_get(g_z,ioff,ioff+nvir-1,ivec,ivec,
     &              dcpl_mb(k_z),nvir)
        do j=1,nvir
          val_cmplx=dcmplx(0.0d0,dbl_mb(k_xreim+j-1))
          call ga_put(g_a,1,1,1,1,val_cmplx,1)
          j1=ioff+j-1
          call ga_add_patch(one_cmplx,g_z,j1,j1,ivec,ivec,
     &                      scl_cmplx,g_a,1 ,1 ,1   ,1   ,
     &                                g_z,j1,j1,ivec,ivec)
        enddo ! end-loop-j
       enddo ! end-loop-i
      endif ! end-if-indrm
      if (.not.ma_pop_stack(l_z))
     $  call errquit('conv2complex2: pop problem with l_z',
     &               555,MA_ERR)
      if (.not.ma_pop_stack(l_xreim))
     $  call errquit('conv2complex2: pop problem with l_xreim',
     &               555,MA_ERR)
      return
      end

      subroutine update_gz_reorim1(g_z,    ! out: = complx(g_xre,g_xim)
     &                             g_xreim,! in : real      arr 
     &                             indrm,  ! in : =1 -> re =2 -> im
     &                             nsub,   ! in : index to sub-block in g_z
     &                             ipm,    ! in : = 1 or 2 index for component
     &                             n,      ! in : = nocc*nvir
     &                             scl,    ! in : scaling factor 
     &                             nvir,
     &                             nocc,
     &                             ivec)
c
c    Purpose: Convert into complex array
c             (g_xre,g_xim) --> g_z
c             structure of g_xre,g_xim:  (n,nvec)  n=nvir*nocc
c             ipm=1,ncomp
c             Usual values: nvec=3 (x,y,z) ncomp=2
c
c    Author: Fredy W. Aquino, Northwestern University
c    Date  : 04-08-12
c Note.- To be used in rohf_hessv_2e3_opt_cmplx() 
c        located in ddscf/rohf_hessv3.F
c        To mimic,
c          call ga_mat_to_vec(g_tmp1,1,nvir,1,nclosed,
c     $                       g_ax_re(ipm),1,ivec,four,'+')  [ scl=four RDFT]
c        for complex g_z (instead of g_ax_re)
c        g_xreim=g_tmp1
c   g_z is history matrix of dim(n1,maxsub)
c       n1=n*ncomp maxsub=maxiter*nvec    
c       n=nocc*nvir maxiter=10 (usually) nvec=3 (x,y,z)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
      integer ipm,indrm,
     &        ivec,ivec1,nvec,
     &        n,n1,nvir,nocc,nsub,
     &        l_xreim,k_xreim,
     &        l_z,k_z,i,j,j1,ioff,ioff1  
      integer g_xreim,
     &        g_z,type
      double precision scl,val_zre,val_zim
      double complex val_cmplx,one_cmplx,scl_cmplx
      one_cmplx =dcmplx(1.0d0,0.0d0)
      scl_cmplx =dcmplx(scl,0.0d0)
      if (indrm.ne.1 .and.
     &    indrm.ne.2) then
       call errquit('update_gz_reorim1: indrm ne 1 or 2',
     &                  0,MA_ERR) 
      endif
      if (.not.MA_Push_Get(mt_dbl,nvir,'hessv jfacs',l_xreim,k_xreim))
     &     call errquit('conv2complex2: cannot allocate xreim',
     &                  nvir, MA_ERR)
      if (.not.MA_Push_Get(mt_dcpl,nvir,'hessv kfacs',l_z,k_z))
     &     call errquit('conv2complex2: cannot allocate z',
     &                  nvir, MA_ERR)
       ioff1=(ipm-1)*n
       ivec1=nsub+ivec
      if   (indrm.eq.1) then ! updating only REAL part
c ++++++++++++++++++
c NOTE.- Assumming the ordering in dbl_mb(k_xreim) is same as in dcpl_mb(k_z)
c ++++++++++++++++++
       do i = ga_nodeid()+1,nocc,ga_nnodes()
        ioff1=(ipm-1)*n
        ivec1=nsub+ivec
        ioff = ioff1+(i-1)*nvir + 1   
        call ga_get(g_xreim,1,nvir,i,i,dbl_mb(k_xreim),nvir)  
        call ga_get(g_z,ioff,ioff+nvir-1,ivec1,ivec1,
     &              dcpl_mb(k_z),nvir)
        do j=1,nvir
          val_zre=dreal(dcpl_mb(k_z+j-1))
          val_zim=dimag(dcpl_mb(k_z+j-1))
          val_cmplx=dcmplx(val_zre+scl*dbl_mb(k_xreim+j-1),val_zim)
          j1=ioff+j-1
          call ga_put(g_z,j1,j1,ivec1,ivec1,val_cmplx,1)
        enddo ! end-loop-j
       enddo ! end-loop-i
      else if (indrm.eq.2) then ! updating only IMAG part
       do i = ga_nodeid()+1,nocc,ga_nnodes()
        ioff = ioff1+(i-1)*nvir + 1   
        call ga_get(g_xreim,1,nvir,i,i,dbl_mb(k_xreim),nvir)  
        call ga_get(g_z,ioff,ioff+nvir-1,ivec1,ivec1,
     &              dcpl_mb(k_z),nvir)
        do j=1,nvir
          val_zre=dreal(dcpl_mb(k_z+j-1))
          val_zim=dimag(dcpl_mb(k_z+j-1))
          val_cmplx=dcmplx(val_zre,val_zim+scl*dbl_mb(k_xreim+j-1))
          j1=ioff+j-1
          call ga_put(g_z,j1,j1,ivec1,ivec1,val_cmplx,1)
        enddo ! end-loop-j
       enddo ! end-loop-i
      endif ! end-if-indrm
      if (.not.ma_pop_stack(l_z))
     $  call errquit('conv2complex2: pop problem with l_z',
     &               555,MA_ERR)
      if (.not.ma_pop_stack(l_xreim))
     $  call errquit('conv2complex2: pop problem with l_xreim',
     &               555,MA_ERR)
      return
      end

      subroutine conv2reim(g_xre,! out : real      arr
     &                     g_xim,! out : imaginary arr
     &                     g_z,  ! in  : = complx(g_xre,g_xim)
     &                     n,    ! in  : n    rows
     &                     nvec, ! in  : nvec columns
     &                     ncomp)! in  : nr. components
c
c    Purpose: Convert into (g_xre,g_xim)
c             g_z -> (g_xre,g_xim) 
c             structure of g_xre(ipm),g_xim(ipm):  (n,nvec)
c             ipm=1,ncomp
c             Usual values: nvec=3 (x,y,z) ncomp=2
c
c    Author: Fredy W. Aquino, Northwestern University
c    Date  : 04-08-12

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
      integer ipm,ivec,idat,
     &        n,nvec,ncomp,
     &        l_z,k_z 
      integer g_xre(ncomp),
     &        g_xim(ncomp),g_z(ncomp)   
      double precision  val_re,val_im

      if (.not.MA_Push_Get(mt_dcpl,n,'conv2reim l_z',
     &                     l_z,k_z))
     &     call errquit('conv2complex: cannot allocate zre',
     &                  n, MA_ERR)
        do ipm=1,ncomp
         call ga_zero(g_xre(ipm))
         call ga_zero(g_xim(ipm))
         do ivec=1,nvec
          call ga_get(g_z(ipm),1,n,ivec,ivec,dcpl_mb(k_z),n)
          do idat=1,n
           val_re=dreal(dcpl_mb(k_z+idat-1))
           val_im=dimag(dcpl_mb(k_z+idat-1))
           call ga_put(g_xre(ipm),idat,idat,ivec,ivec,val_re,1)
           call ga_put(g_xim(ipm),idat,idat,ivec,ivec,val_im,1)
          enddo ! end-loop-idat
         enddo ! end-loop-ivec
        enddo ! end-loop-ipm
      if (.not.ma_pop_stack(l_z))
     $  call errquit('conv2complex: pop problem with l_zim',
     &               555,MA_ERR)
      return
      end

      subroutine conv2reim_rhs(
     &                     g_xre,! out : real      arr
     &                     g_xim,! out : imaginary arr
     &                     g_z,  ! in  : = complx(g_xre,g_xim)
     &                     n,    ! in  : n    rows
     &                     nvec, ! in  : nvec columns
     &                     ncomp,! in  : nr. components
     &                     nsub) ! in  : =1,2=g_b,g_z index to subspace
c
c    Purpose: Convert into (g_xre,g_xim)
c             g_z -> (g_xre,g_xim) 
c             structure of g_xre(ipm),g_xim(ipm):  (n,nvec)
c             ipm=1,ncomp
c             Usual values: nvec=3 (x,y,z) ncomp=2
c
c    Author: Fredy W. Aquino, Northwestern University
c    Date  : 05-07-12

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
      integer ipm,ivec,idat,nsub,shift,ivec1,
     &        n,nvec,ncomp,
     &        l_z,k_z 
      integer g_xre(ncomp),
     &        g_xim(ncomp),g_z(ncomp)   
      double precision  val_re,val_im

      if (.not.MA_Push_Get(mt_dcpl,n,'conv2reim l_z',
     &                     l_z,k_z))
     &     call errquit('conv2complex: cannot allocate zre',
     &                  n, MA_ERR)
        shift=nvec*(nsub-1)
        do ipm=1,ncomp
         call ga_zero(g_xre(ipm))
         call ga_zero(g_xim(ipm))
         do ivec=1,nvec
          ivec1=shift+ivec
          call ga_get(g_z(ipm),1,n,ivec,ivec,dcpl_mb(k_z),1)
          do idat=1,n
           val_re=dreal(dcpl_mb(k_z+idat-1))
           val_im=dimag(dcpl_mb(k_z+idat-1))
           call ga_put(g_xre(ipm),idat,idat,ivec1,ivec1,val_re,1)
           call ga_put(g_xim(ipm),idat,idat,ivec1,ivec1,val_im,1)
          enddo ! end-loop-idat
         enddo ! end-loop-ivec
        enddo ! end-loop-ipm
      if (.not.ma_pop_stack(l_z))
     $  call errquit('conv2complex: pop problem with l_zim',
     &               555,MA_ERR)
      return
      end

      subroutine conv2reim1(g_xre,! out : real      arr
     &                      g_xim,! out : imaginary arr
     &                      g_z,  ! in  : = complx(g_xre,g_xim)
     &                      nsub, ! in  : pointer to block
     &                      nvir, ! in  : nr. virtual  MOs
     &                      nocc, ! in  : nr. occupied MOs
     &                      ipm,  ! in  : =1,2 components indices
     &                      n,    ! in  : n    rows
     &                      nvec) ! in  : nvec columns
c
c     Purpose: Extract into (g_xre,g_xim)
c             g_z -> (g_xre,g_xim) 
c             g_z = g_zr1 or g_Az1
c             dim(g_zr1)=(n1,nvec)    nvec=3 (x,y and z)
c             dim(g_Az1)=(n1,maxsub)  maxsub=10*nvec
c             n1=n*ncomp ncomp=2 (usually it stands for +/- solutions)
c             structure of g_xre(ipm),g_xim(ipm):  (n,nvec)
c             ipm=1,ncomp
c             For case g_z=g_zr1: nsub should be 0
c             For case g_z=g_Az1: nsub should be real value
c             Usual values: nvec=3 (x,y,z) ncomp=2
c
c    Author: Fredy W. Aquino, Northwestern University
c    Date  : 04-21-12

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
      integer ipm,ivec,ivec1,
     &        a1,a2,b1,b2,
     &        n,nvec,ncomp,nsub,
     &        nocc,nvir,i,j,j1,
     &        ioff,ioff1,ioff2,
     &        l_z,k_z 
      integer g_xre,g_xim,g_z  
      double precision val_re,val_im

      if (.not.MA_Push_Get(mt_dcpl,nvir,'conv2reim l_z',
     &                     l_z,k_z))
     &     call errquit('conv2complex: cannot allocate zre',
     &                  nvir, MA_ERR)
         call ga_zero(g_xre)
         call ga_zero(g_xim)
         a1=(ipm-1)*n+1
         a2=a1+n-1
         b1=nsub+1
         b2=nsub+nvec
         ivec1=1
         ioff1=(ipm-1)*n
         do ivec=b1,b2
          do i = ga_nodeid()+1,nocc,ga_nnodes()
           ioff = (i-1)*nvir + 1 
           ioff2=ioff1+ioff    
           call ga_get(g_z,ioff2,ioff2+nvir-1,ivec,ivec,
     &                 dcpl_mb(k_z),nvir)
           do j=1,nvir
            j1=ioff+j-1
            val_re=dreal(dcpl_mb(k_z+j-1))
            val_im=dimag(dcpl_mb(k_z+j-1))
            call ga_put(g_xre,j1,j1,ivec1,ivec1,val_re,1)    
            call ga_put(g_xim,j1,j1,ivec1,ivec1,val_im,1)       
           enddo ! end-loop-j
          enddo ! end-loop-i
          ivec1=ivec1+1
         enddo ! end-loop-ivec
      if (.not.ma_pop_stack(l_z))
     $  call errquit('conv2complex: pop problem with l_zim',
     &               555,MA_ERR)
      return
      end

      subroutine conv2reim1_u(
     &                  g_xre,! out : real      arr
     &                  g_xim,! out : imaginary arr
     &                  g_z,  ! in  : = complx(g_xre,g_xim)
     &                  nsub, ! in  : pointer to block
     &                  npol, ! in  : nr. polarizations
     &                  nvir, ! in  : nr. virtual  MOs
     &                  nocc, ! in  : nr. occupied MOs
     &                  ipm,  ! in  : =1,2 components indices
     &                  n,    ! in  : n    rows
     &                  nvec) ! in  : nvec columns
c
c     Purpose: Extract into (g_xre,g_xim)
c             g_z -> (g_xre,g_xim) 
c             g_z = g_zr1 or g_Az1
c             dim(g_zr1)=(n1,nvec)    nvec=3 (x,y and z)
c             dim(g_Az1)=(n1,maxsub)  maxsub=10*nvec
c             n1=n*ncomp ncomp=2 (usually it stands for +/- solutions)
c             structure of g_xre(ipm),g_xim(ipm):  (n,nvec)
c             ipm=1,ncomp
c             For case g_z=g_zr1: nsub should be 0
c             For case g_z=g_Az1: nsub should be real value
c             Usual values: nvec=3 (x,y,z) ncomp=2
c
c    Author: Fredy W. Aquino, Northwestern University
c    Date  : 04-21-12

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
      integer ipm,ivec,ivec1,
     &        a1,a2,b1,b2,
     &        n,nvec,ncomp,nsub,
     &        ipol,npol,
     &        nocc(npol),nvir(npol),
     &        i,j,j1,shift,
     &        ioff,ioff1,ioff2,
     &        l_z,k_z 
      integer g_xre,g_xim,g_z  
      double precision val_re,val_im
         call ga_zero(g_xre)
         call ga_zero(g_xim)
         b1=nsub+1
         b2=nsub+nvec
         ivec1=1
         ioff1=(ipm-1)*n ! n=sum_{i=1,npol} (nocc*nvir)_i
         do ivec=b1,b2
          do ipol=1,npol
            if (.not.MA_Push_Get(mt_dcpl,nvir(ipol),
     &          'conv2reim l_z',l_z,k_z))
     &      call errquit('conv2complex: cannot allocate zre',
     &                  nvir(ipol), MA_ERR)
           shift=nocc(1)*nvir(1)*(ipol-1)
           do i = ga_nodeid()+1,nocc(ipol),ga_nnodes()
            ioff = shift+(i-1)*nvir(ipol) + 1 
            ioff2=ioff1+ioff    
            call ga_get(g_z,ioff2,ioff2+nvir(ipol)-1,ivec,ivec,
     &                  dcpl_mb(k_z),nvir(ipol))
            do j=1,nvir(ipol)
             j1=ioff+j-1
             val_re=dreal(dcpl_mb(k_z+j-1))
             val_im=dimag(dcpl_mb(k_z+j-1))
             call ga_put(g_xre,j1,j1,ivec1,ivec1,val_re,1)    
             call ga_put(g_xim,j1,j1,ivec1,ivec1,val_im,1)       
            enddo ! end-loop-j
           enddo ! end-loop-i
           if (.not.ma_pop_stack(l_z))
     $      call errquit('conv2complex: pop problem with l_zim',
     &               555,MA_ERR)
          enddo ! end-loop-ipol
          ivec1=ivec1+1
         enddo ! end-loop-ivec

      return
      end

      subroutine getreorim(g_xreim,! out : real or im arr
     &                     g_z,    ! in  : = complx(g_xre,g_xim)
     &                     nvir,   ! in  : nr. virtual  MOs
     &                     nocc,   ! in  : nr. occupied MOs
     &                     indrm)  ! in  : =1 -> re =2 -> im
c
c    Purpose: Convert into (g_xre,g_xim)
c             g_z -> (g_xre,g_xim) 
c             structure of g_xre(ipm),g_xim(ipm):  (n,nvec)
c             ipm=1,ncomp
c             Usual values: nvec=3 (x,y,z) ncomp=2
c
c    Author: Fredy W. Aquino, Northwestern University
c    Date  : 04-08-12
c --> This only works for closed shell where we have one single
c     set of (nvir,nocc)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
      integer ipm,ivec,indrm,
     &        n,nvec, ! obtained from ga_inquire(g_z)
     &        ncomp,
     &        l_z,k_z 
      integer g_xreim,g_z,type,
     &        i,j,j1,ioff,nocc,nvir  
      double precision  val_re,val_im
      if (indrm.ne.1 .and.
     &    indrm.ne.2) then
       call errquit('getreorim: indrm ne 1 or 2',
     &                  0,MA_ERR) 
      endif
      call ga_inquire(g_z,type,n,nvec) ! get (n,nvec)
      if (.not.MA_Push_Get(mt_dcpl,nvir,'conv2reim l_z',
     &                     l_z,k_z))
     &     call errquit('getreorim: cannot allocate k_z',
     &                  nvir, MA_ERR)
      call ga_zero(g_xreim)
      if   (indrm.eq.1) then ! copying only REAL part
       do ivec=1,nvec
        do i = ga_nodeid()+1,nocc,ga_nnodes()
         ioff = (i-1)*nvir + 1    
         call ga_get(g_z,ioff,ioff+nvir-1,ivec,ivec,
     &               dcpl_mb(k_z),nvir)
         do j=1,nvir
          val_re=dreal(dcpl_mb(k_z+j-1))
          j1=ioff+j-1
          call ga_put(g_xreim,j1,j1,ivec,ivec,val_re,1)         
         enddo ! end-loop-j
        enddo ! end-loop-i
       enddo ! end-loop-ivec
      else if (indrm.eq.2) then ! copying only IMAG part
       do ivec=1,nvec
        do i = ga_nodeid()+1,nocc,ga_nnodes()
         ioff = (i-1)*nvir + 1     
         call ga_get(g_z,ioff,ioff+nvir-1,ivec,ivec,
     $               dcpl_mb(k_z),nvir)
         do j=1,nvir
          val_im=dimag(dcpl_mb(k_z+j-1))
          j1=ioff+j-1
          call ga_put(g_xreim,j1,j1,ivec,ivec,val_im,1)         
         enddo ! end-loop-j
        enddo ! end-loop-i
       enddo ! end-loop-ivec  
      endif ! end-if-indrm
      if (.not.ma_pop_stack(l_z))
     $  call errquit('getreorim: pop problem with l_zim',
     &               555,MA_ERR)
      return
      end

      subroutine getreorim_u(g_xreim,! out : real or im arr
     &                       g_z,    ! in  : = complx(g_xre,g_xim)
     &                       shift,  ! in  : = 0 -> ipol=1, nocc*nvirt(1) -> ipol=2
     &                       nvir,   ! in  : nr. virtual  MOs
     &                       nocc,   ! in  : nr. occupied MOs
     &                       indrm)  ! in  : =1 -> re =2 -> im
c
c    Purpose: Convert into (g_xre,g_xim)
c             g_z -> (g_xre,g_xim) 
c             structure of g_xre(ipm),g_xim(ipm):  (n,nvec)
c             ipm=1,ncomp
c             Usual values: nvec=3 (x,y,z) ncomp=2
c
c    Author: Fredy W. Aquino, Northwestern University
c    Date  : 04-08-12
c --> This only works for closed shell where we have one single
c     set of (nvir,nocc)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
      integer ipm,ivec,indrm,
     &        n,nvec, ! obtained from ga_inquire(g_z)
     &        ncomp,
     &        l_z,k_z 
      integer g_xreim,g_z,type,
     &        i,j,j1,ioff,ioff2,nocc,nvir,
     &        shift  
      double precision  val_re,val_im
      if (indrm.ne.1 .and.
     &    indrm.ne.2) then
       call errquit('getreorim: indrm ne 1 or 2',
     &                  0,MA_ERR) 
      endif
      call ga_inquire(g_z,type,n,nvec) ! get (n,nvec)
      if (.not.MA_Push_Get(mt_dcpl,nvir,'conv2reim l_z',
     &                     l_z,k_z))
     &     call errquit('getreorim: cannot allocate k_z',
     &                  nvir, MA_ERR)
      call ga_zero(g_xreim)
      if   (indrm.eq.1) then ! copying only REAL part
       do ivec=1,nvec
        do i = ga_nodeid()+1,nocc,ga_nnodes()
         ioff  = (i-1)*nvir + 1 
         ioff2 = shift+ioff  
         call ga_get(g_z,ioff2,ioff2+nvir-1,ivec,ivec,
     &               dcpl_mb(k_z),nvir)
         do j=1,nvir
          val_re=dreal(dcpl_mb(k_z+j-1))
          j1=ioff+j-1
          call ga_put(g_xreim,j1,j1,ivec,ivec,val_re,1)         
         enddo ! end-loop-j
        enddo ! end-loop-i
       enddo ! end-loop-ivec
      else if (indrm.eq.2) then ! copying only IMAG part
       do ivec=1,nvec
        do i = ga_nodeid()+1,nocc,ga_nnodes()
         ioff = (i-1)*nvir + 1 
         ioff2 = shift+ioff     
         call ga_get(g_z,ioff2,ioff2+nvir-1,ivec,ivec,
     $               dcpl_mb(k_z),nvir)
         do j=1,nvir
          val_im=dimag(dcpl_mb(k_z+j-1))
          j1=ioff+j-1
          call ga_put(g_xreim,j1,j1,ivec,ivec,val_im,1)         
         enddo ! end-loop-j
        enddo ! end-loop-i
       enddo ! end-loop-ivec  
      endif ! end-if-indrm
      if (.not.ma_pop_stack(l_z))
     $  call errquit('getreorim: pop problem with l_zim',
     &               555,MA_ERR)
      return
      end

      subroutine getreorim_u1(
     &                       g_xreim,! out : real or im arr
     &                       g_z,    ! in  : = complx(g_xre,g_xim)
     &                       npol,   ! in  : nr. polarizations
     &                       nvir,   ! in  : nr. virtual  MOs
     &                       nocc,   ! in  : nr. occupied MOs
     &                       indrm)  ! in  : =1 -> re =2 -> im
c
c    Purpose: Convert into (g_xre,g_xim)
c             g_z -> (g_xre,g_xim) 
c             structure of g_xre(ipm),g_xim(ipm):  (n,nvec)
c             ipm=1,ncomp
c             Usual values: nvec=3 (x,y,z) ncomp=2
c
c    Author: Fredy W. Aquino, Northwestern University
c    Date  : 04-08-12
c --> This only works for closed shell where we have one single
c     set of (nvir,nocc)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
      integer ipm,ivec,indrm,
     &        n,nvec, ! obtained from ga_inquire(g_z)
     &        ncomp,
     &        l_z,k_z 
      integer g_xreim,g_z,type,
     &        i,j,j1,ioff,
     &        ipol,npol,nocc(npol),nvir(npol),
     &        shift  
      double precision  val_re,val_im
      if (indrm.ne.1 .and.
     &    indrm.ne.2) then
       call errquit('getreorim: indrm ne 1 or 2',
     &                  0,MA_ERR) 
      endif
      call ga_inquire(g_z,type,n,nvec) ! get (n,nvec)
      call ga_zero(g_xreim)
      if   (indrm.eq.1) then ! copying only REAL part
       do ivec=1,nvec
        do ipol=1,npol
          if (.not.MA_Push_Get(mt_dcpl,nvir(ipol),
     &        'conv2reim l_z',l_z,k_z))
     &     call errquit('getreorim: cannot allocate k_z',
     &                  nvir(ipol), MA_ERR)
         shift=nocc(1)*nvir(1)*(ipol-1)
         do i = ga_nodeid()+1,nocc(ipol),ga_nnodes()
          ioff = shift+(i-1)*nvir(ipol) + 1    
          call ga_get(g_z,ioff,ioff+nvir(ipol)-1,ivec,ivec,
     &                dcpl_mb(k_z),nvir(ipol))
          do j=1,nvir(ipol)
           val_re=dreal(dcpl_mb(k_z+j-1))
           j1=ioff+j-1
           call ga_put(g_xreim,j1,j1,ivec,ivec,val_re,1)         
          enddo ! end-loop-j
         enddo ! end-loop-i
         if (.not.ma_pop_stack(l_z))
     $    call errquit('getreorim: pop problem with l_zim',
     &               555,MA_ERR)
        enddo ! end-loop-ipol
       enddo ! end-loop-ivec
      else if (indrm.eq.2) then ! copying only IMAG part
       do ivec=1,nvec
        do ipol=1,npol
          if (.not.MA_Push_Get(mt_dcpl,nvir(ipol),
     &        'conv2reim l_z',l_z,k_z))
     &     call errquit('getreorim: cannot allocate k_z',
     &                  nvir(ipol), MA_ERR)
         shift=nocc(1)*nvir(1)*(ipol-1)
         do i = ga_nodeid()+1,nocc(ipol),ga_nnodes()
          ioff = shift+(i-1)*nvir(ipol) + 1     
          call ga_get(g_z,ioff,ioff+nvir(ipol)-1,ivec,ivec,
     $                dcpl_mb(k_z),nvir(ipol))
          do j=1,nvir(ipol)
           val_im=dimag(dcpl_mb(k_z+j-1))
           j1=ioff+j-1
           call ga_put(g_xreim,j1,j1,ivec,ivec,val_im,1)         
          enddo ! end-loop-j
         enddo ! end-loop-i
         if (.not.ma_pop_stack(l_z))
     $    call errquit('getreorim: pop problem with l_zim',
     &               555,MA_ERR)
        enddo ! end-loop-ipol
       enddo ! end-loop-ivec  
      endif ! end-if-indrm

      return
      end

      subroutine getreorim1(g_xreim,! out : real or im arr
     &                      g_z,    ! in  : = complx(g_xre,g_xim)
     &                      nsub,   ! in  : subblock index
     &                      ipm,    ! in  : = 1,2 to access slctd component
     &                      nvir,   ! in  : nr. virtual  MOs
     &                      nocc,   ! in  : nr. occupied MOs
     &                      indrm)  ! in  : =1 -> re =2 -> im
c
c    Purpose: Convert into (g_xre,g_xim)
c             g_z -> (g_xre,g_xim) 
c             structure of g_xre(ipm),g_xim(ipm):  (n,nvec)
c             ipm=1,ncomp
c             Usual values: nvec=3 (x,y,z) ncomp=2
c
c    Author: Fredy W. Aquino, Northwestern University
c    Date  : 04-24-12
c    g_z : history matrix (g_Az1 or g_z1)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
      integer ipm,ivec,ivec1,indrm,
     &        n,nvec, ! obtained from ga_inquire(g_z)
     &        ncomp,nsub,
     &        l_z,k_z 
      integer g_xreim,g_z,type,
     &        i,j,j1,ioff,ioff1,ioff2,nocc,nvir  
      double precision  val_re,val_im
      if (indrm.ne.1 .and.
     &    indrm.ne.2) then
       call errquit('getreorim1: indrm ne 1 or 2',
     &                  0,MA_ERR) 
      endif
      call ga_inquire(g_xreim,type,n,nvec) ! get (n,nvec)

      if (.not.MA_Push_Get(mt_dcpl,nvir,'conv2reim l_z',
     &                     l_z,k_z))
     &     call errquit('getreorim: cannot allocate k_z',
     &                  nvir, MA_ERR)
      call ga_zero(g_xreim)
      if   (indrm.eq.1) then ! copying only REAL part
       ivec1=1
       ioff1=(ipm-1)*n
       do ivec=nsub+1,nsub+nvec
        do i = ga_nodeid()+1,nocc,ga_nnodes()
         ioff = (i-1)*nvir + 1   
         ioff2=ioff1+ioff 
         call ga_get(g_z,ioff2,ioff2+nvir-1,ivec,ivec,
     &               dcpl_mb(k_z),nvir)
         do j=1,nvir
          val_re=dreal(dcpl_mb(k_z+j-1))
          j1=ioff+j-1
          call ga_put(g_xreim,j1,j1,ivec1,ivec1,val_re,1)         
         enddo ! end-loop-j
        enddo ! end-loop-i
        ivec1=ivec1+1
       enddo ! end-loop-ivec
      else if (indrm.eq.2) then ! copying only IMAG part
       ivec1=1
       ioff1=(ipm-1)*n
       do ivec=nsub+1,nsub+nvec
        do i = ga_nodeid()+1,nocc,ga_nnodes()
         ioff = (i-1)*nvir + 1   
         ioff2=ioff1+ioff   
         call ga_get(g_z,ioff2,ioff2+nvir-1,ivec,ivec,
     $               dcpl_mb(k_z),nvir)
         do j=1,nvir
          val_im=dimag(dcpl_mb(k_z+j-1))
          j1=ioff+j-1
          call ga_put(g_xreim,j1,j1,ivec1,ivec1,val_im,1)         
         enddo ! end-loop-j
        enddo ! end-loop-i
        ivec1=ivec1+1
       enddo ! end-loop-ivec  
      endif ! end-if-indrm
      if (.not.ma_pop_stack(l_z))
     $  call errquit('getreorim: pop problem with l_zim',
     &               555,MA_ERR)
      return
      end

      subroutine getreorim1_u1(
     &                      g_xreim,! out : real or im arr
     &                      g_z,    ! in  : = complx(g_xre,g_xim)
     &                      nsub,   ! in  : subblock index
     &                      ipm,    ! in  : = 1,2 to access slctd component
     &                      npol,
     &                      nvir,   ! in  : nr. virtual  MOs
     &                      nocc,   ! in  : nr. occupied MOs
     &                      indrm)  ! in  : =1 -> re =2 -> im
c
c    Purpose: Convert into (g_xre,g_xim)
c             g_z -> (g_xre,g_xim) 
c             structure of g_xre(ipm),g_xim(ipm):  (n,nvec)
c             ipm=1,ncomp
c             Usual values: nvec=3 (x,y,z) ncomp=2
c
c    Author: Fredy W. Aquino, Northwestern University
c    Date  : 04-24-12
c    g_z : history matrix (g_Az1 or g_z1)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
      integer ipm,ivec,ivec1,indrm,
     &        n,nvec, ! obtained from ga_inquire(g_z)
     &        ncomp,nsub,
     &        l_z,k_z 
      integer g_xreim,g_z,type,
     &        i,j,j1,ioff,ioff1,ioff2,
     &        ipol,npol,nocc(npol),nvir(npol),
     &        shift  
      double precision  val_re,val_im
      if (indrm.ne.1 .and.
     &    indrm.ne.2) then
       call errquit('getreorim1_u1: indrm ne 1 or 2',
     &                  0,MA_ERR) 
      endif
      call ga_inquire(g_xreim,type,n,nvec) ! get (n,nvec)
      call ga_zero(g_xreim)
      if   (indrm.eq.1) then ! copying only REAL part
       ivec1=1
       ioff1=(ipm-1)*n ! n=sum_{i=1,npol} (nocc*nvir)(i)
       do ivec=nsub+1,nsub+nvec
        do ipol=1,npol
         if (.not.MA_Push_Get(mt_dcpl,nvir(ipol),'conv2reim l_z',
     &                     l_z,k_z))
     &     call errquit('getreorim: cannot allocate k_z',
     &                  nvir(ipol), MA_ERR)
         shift=nocc(1)*nvir(1)*(ipol-1)
         do i = ga_nodeid()+1,nocc(ipol),ga_nnodes()
          ioff = shift+(i-1)*nvir(ipol) + 1   
          ioff2=ioff1+ioff 
          call ga_get(g_z,ioff2,ioff2+nvir(ipol)-1,ivec,ivec,
     &               dcpl_mb(k_z),1)
          do j=1,nvir(ipol)
           val_re=dreal(dcpl_mb(k_z+j-1))
           j1=ioff+j-1
           call ga_put(g_xreim,j1,j1,ivec1,ivec1,val_re,1)         
          enddo ! end-loop-j
         enddo ! end-loop-i
         if (.not.ma_pop_stack(l_z))
     $   call errquit('getreorim: pop problem with l_zim',
     &               555,MA_ERR)
        enddo ! end-loop-ipol
        ivec1=ivec1+1
       enddo ! end-loop-ivec
      else if (indrm.eq.2) then ! copying only IMAG part
       ivec1=1
       ioff1=(ipm-1)*n
       do ivec=nsub+1,nsub+nvec
        do ipol=1,npol
         if (.not.MA_Push_Get(mt_dcpl,nvir(ipol),'conv2reim l_z',
     &                     l_z,k_z))
     &     call errquit('getreorim: cannot allocate k_z',
     &                  nvir(ipol), MA_ERR)
         shift=nocc(1)*nvir(1)*(ipol-1)
         do i = ga_nodeid()+1,nocc(ipol),ga_nnodes()
          ioff = shift+(i-1)*nvir(ipol) + 1   
          ioff2=ioff1+ioff   
          call ga_get(g_z,ioff2,ioff2+nvir(ipol)-1,ivec,ivec,
     $                dcpl_mb(k_z),1)
          do j=1,nvir(ipol)
           val_im=dimag(dcpl_mb(k_z+j-1))
           j1=ioff+j-1
           call ga_put(g_xreim,j1,j1,ivec1,ivec1,val_im,1)         
          enddo ! end-loop-j
         enddo ! end-loop-i
         if (.not.ma_pop_stack(l_z))
     $   call errquit('getreorim: pop problem with l_zim',
     &               555,MA_ERR)
        enddo ! end-loop-ipol
        ivec1=ivec1+1
       enddo ! end-loop-ivec  
      endif ! end-if-indrm
      return
      end

      subroutine getreorim1_u(
     &                      g_xreim,! out : real or im arr
     &                      g_z,    ! in  : = complx(g_xre,g_xim)
     &                      nsub,   ! in  : subblock index
     &                      shift,  ! in  : = 0 for spin 1 ncomp*(nocc*nvir)(1) for spin 2
     &                      ipm,    ! in  : = 1,2 to access slctd component
     &                      nvir,   ! in  : nr. virtual  MOs
     &                      nocc,   ! in  : nr. occupied MOs
     &                      indrm)  ! in  : =1 -> re =2 -> im
c
c    Purpose: Convert into (g_xre,g_xim)
c             g_z -> (g_xre,g_xim) 
c             structure of g_xre(ipm),g_xim(ipm):  (n,nvec)
c             ipm=1,ncomp
c             Usual values: nvec=3 (x,y,z) ncomp=2
c
c    Author: Fredy W. Aquino, Northwestern University
c    Date  : 04-24-12
c    g_z : history matrix (g_Az1 or g_z1)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
      integer ipm,ivec,ivec1,indrm,
     &        n,nvec, ! obtained from ga_inquire(g_z)
     &        ncomp,nsub,shift,
     &        l_z,k_z 
      integer g_xreim,g_z,type,
     &        i,j,j1,ioff,ioff1,ioff2,nocc,nvir  
      double precision  val_re,val_im
      if (indrm.ne.1 .and.
     &    indrm.ne.2) then
       call errquit('getreorim1: indrm ne 1 or 2',
     &                  0,MA_ERR) 
      endif
      call ga_inquire(g_xreim,type,n,nvec) ! get (n,nvec)

      if (.not.MA_Push_Get(mt_dcpl,nvir,'conv2reim l_z',
     &                     l_z,k_z))
     &     call errquit('getreorim: cannot allocate k_z',
     &                  nvir, MA_ERR)
      call ga_zero(g_xreim)
      if   (indrm.eq.1) then ! copying only REAL part
       ivec1=1
       ioff1=(ipm-1)*n+shift
       do ivec=nsub+1,nsub+nvec
        do i = ga_nodeid()+1,nocc,ga_nnodes()
         ioff = (i-1)*nvir + 1   
         ioff2=ioff1+ioff 
         call ga_get(g_z,ioff2,ioff2+nvir-1,ivec,ivec,
     &               dcpl_mb(k_z),1)
         do j=1,nvir
          val_re=dreal(dcpl_mb(k_z+j-1))
          j1=ioff+j-1
          call ga_put(g_xreim,j1,j1,ivec1,ivec1,val_re,1)         
         enddo ! end-loop-j
        enddo ! end-loop-i
        ivec1=ivec1+1
       enddo ! end-loop-ivec
      else if (indrm.eq.2) then ! copying only IMAG part
       ivec1=1
       ioff1=shift+(ipm-1)*n
       do ivec=nsub+1,nsub+nvec
        do i = ga_nodeid()+1,nocc,ga_nnodes()
         ioff = (i-1)*nvir + 1   
         ioff2=ioff1+ioff   
         call ga_get(g_z,ioff2,ioff2+nvir-1,ivec,ivec,
     $               dcpl_mb(k_z),1)
         do j=1,nvir
          val_im=dimag(dcpl_mb(k_z+j-1))
          j1=ioff+j-1
          call ga_put(g_xreim,j1,j1,ivec1,ivec1,val_im,1)         
         enddo ! end-loop-j
        enddo ! end-loop-i
        ivec1=ivec1+1
       enddo ! end-loop-ivec  
      endif ! end-if-indrm
      if (.not.ma_pop_stack(l_z))
     $  call errquit('getreorim: pop problem with l_zim',
     &               555,MA_ERR)
      return
      end

      subroutine updating_Az1_z1_zr1(
     &                      g_Az1, ! in/ou:
     &                      g_z1,  ! in/ou:
     &                      g_zr1, ! in/ou:
     &                      g_Az,  ! in   :
     &                      g_z,   ! in   :
     &                      g_zr,  ! in   :
     &                      nvec,  ! in   :
     &                      ncomp, ! in   :
     &                      nsub,  ! in   :
     &                      n)     ! in   :
c
c Author : Fredy W. Aquino, Northwestern University
c Purpose: Update (g_Az1,g_z1,g_zr1)
c Date   : 03-15-12

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
       integer ncomp,nvec,nsub,n
       integer g_Az1,g_z1,g_zr1,
     &         g_Az(ncomp),
     &         g_z(ncomp),
     &         g_zr(ncomp)
       integer p1,p2,m1,m2,ipm
       p1=nsub+1
       p2=nsub+nvec
       m1=1
       m2=n
       do ipm=1,ncomp
         call ga_copy_patch('n',g_Az(ipm),1 ,n ,1 ,nvec, 
     $                          g_Az1    ,m1,m2,p1,p2)  
         call ga_copy_patch('n',g_z(ipm) ,1 ,n ,1 ,nvec, 
     $                          g_z1     ,m1,m2,p1,p2)   
         call ga_copy_patch('n',g_zr(ipm),1 ,n ,1 ,nvec, 
     $                          g_zr1    ,m1,m2,1 ,nvec)   
         m1=m1+n
         m2=m2+n    
       enddo ! end-loop-ipm
       return
       end
c -------------- solve_rlineq ---------------- START
       subroutine solve_xlineq(
     &              g_x,   ! in/out: updated solution
     &              g_Ax1, ! in    : history of g_Az
     &              g_x1,  ! in    : history of g_z
     &              g_xr1, ! in    : history of g_zr
     &              nsub,  ! in    : subspace length
     &              nvec,  ! in    : increment of subspace
     &              ncomp, ! in    : nr. components
     &              n,     ! in    : nr. elements per comp.
     &              iter,  ! in    : iteration nr.
     &              debug1)! in   : =.true. show debug printouts
c
c Author : Fredy W. Aquino, Northwestern University
c Purpose: Solve 'complex' linear equation using real 'history' GA arrays
c          
c Date   : 03-15-12
c
c    dim(g_Az1)=dim(z1)=(ncomp*n,maxsub)
c    dim(g_zr1)=(ncomp*n,nvec)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
        integer nsub,nvec,ncomp,n,n1,iter,ipm
        integer g_a,g_b,g_c,
     &          g_x1,g_Ax1,g_xr1,
     &          g_x(ncomp)
        logical debug1
        double precision one,mone,zero
        parameter (one=1.0d0,mone=-1.0d0,zero=0.0d0)
        external ga_svd_solve_seq,update_g_x
c       Form and solve the subspace equations using SVD in order
c       to manage near linear dependence in the subspace.
        n1=ncomp*n     
        if (.not. ga_create(MT_DBL, nsub, nsub,
     &     'solve_rlineq: A',0,0,g_a))
     $     call errquit('solve_rlineq: failed allocating g_a', 
     &                  nsub,GA_ERR)
        if (.not. ga_create(MT_DBL, nsub, nvec,
     &     'solve_rlineq: B',0,0,g_b))
     $     call errquit('solve_rlineq: failed allocating g_b', 
     &                  nsub,GA_ERR)
        if (.not. ga_create(MT_DBL, nsub, nvec,
     &     'solve_rlineq: C',0,0,g_c))
     $     call errquit('solve_rlineq: failed allocating g_c', 
     &                  nsub,GA_ERR)
        call ga_zero(g_a)
        call ga_zero(g_b)
        call ga_zero(g_c)

        if (debug1) then
         if (ga_nodeid().eq.0) 
     &   write(*,*) '-------BEF:g_x1(',iter,')-------START'
         call ga_print(g_x1)
         if (ga_nodeid().eq.0) 
     &   write(*,*) '-------BEF:g_x1(',iter,')-------END'
        endif ! end-if-debug1
        call ga_dgemm('t','n',nsub,nsub,n1,one,
     &                g_x1,g_Ax1,zero,g_a)
        call ga_dgemm('t','n',nsub,nvec,n1,one,
     &                g_x1,g_xr1,zero,g_b)

        if (debug1) then
          if (ga_nodeid().eq.0)
     &    write(*,*) '-------g_a--------START'
          call ga_print(g_a)
          if (ga_nodeid().eq.0)
     &    write(*,*) '-------g_a--------END'
          if (ga_nodeid().eq.0)
     &    write(*,*) '-------g_b--------START'
          call ga_print(g_b)
          if (ga_nodeid().eq.0)
     &    write(*,*) '-------g_b--------END'
        endif ! end-if-debug

c       The threshold used here should reflect the accuracy in the
c       products.  If very accurate products are used,
c       then there is big
c       advantage for small cases (maxsub close to n) in using a very
c       small threshold in the SVD solve (e.g., 1e-14), but for more
c       realistic examples (maxsub << n) there is only a little
c       advantage and in the precence of real noise in the products
c       screening with a realistic threshold is important.

        call ga_svd_solve_seq(g_a,g_b,g_c,1d-14)

        if (debug1) then
        if (ga_nodeid().eq.0)
     &   write(*,*) '-------g_c(',iter,')--------START'
        call ga_print(g_c)
        if (ga_nodeid().eq.0)
     &   write(*,*) '-------g_c(',iter,')--------END'
        endif ! end-if-debug1

        call ga_dgemm('n','n',n1,nvec,nsub,mone,
     &                g_Ax1,g_c,one,g_xr1)

        if (debug1) then
          if (ga_nodeid().eq.0) then
           write(*,10) iter
   10       format('---------g_xr1-1(',i3,')-----START')
          endif
          call ga_print(g_xr1)
          if (ga_nodeid().eq.0) then
           write(*,6) iter
    6      format('---------g_xr1-1(',i3,')-----END')
          endif
        endif ! end-if-debug1

        call update_g_x(g_x,  ! in/ou: solution updated
     &                  g_xr1,! in   : added to g_z
     &                  ncomp,! in   : nr. components
     &                  nvec, ! in   : (x,y,z)
     &                  n)    ! in   : vector length
        if (debug1) then
         do ipm=1,ncomp
          if (ga_nodeid().eq.0) then
           write(*,2) ipm,iter
    2      format('---------g_x-1(',i3,',',i3,')-----START')
          endif
          call ga_print(g_x(ipm))
          if (ga_nodeid().eq.0) then
           write(*,3) ipm,iter
    3      format('---------g_x-1(',i3,',',i3,')-----END')
          endif
         enddo ! end-loop-ipm
        endif ! end-if-debug1
        call ga_zero(g_xr1)
        call ga_dgemm('n','n',n1,nvec,nsub,one,
     &                g_x1,g_c,zero,g_xr1)

        if (debug1) then
          if (ga_nodeid().eq.0) then
           write(*,7) iter
    7      format('---------g_xr1-2(',i3,')-----START')
          endif
          call ga_print(g_xr1)
          if (ga_nodeid().eq.0) then
           write(*,9) iter
    9      format('---------g_xr1-2(',i3,')-----END')
          endif
        endif ! end-if-debug1

        call update_g_x(g_x,  ! in/ou: solution updated
     &                  g_xr1,! in   : added to g_z
     &                  ncomp,! in   : nr. components
     &                  nvec, ! in   : (x,y,z)
     &                  n)    ! in   : vector length

        if (debug1) then
        do ipm=1,ncomp
          if (ga_nodeid().eq.0) then
           write(*,4) ipm,iter
    4      format('---------g_x-2(',i3,',',i3,')-----START')
          endif
          call ga_print(g_x(ipm))
          if (ga_nodeid().eq.0) then
           write(*,5) ipm,iter
    5      format('---------g_x-2(',i3,',',i3,')-----END')
          endif
        enddo ! end-loop-ipm
        endif ! end-if-debug1

        if (.not. ga_destroy(g_a)) call errquit
     &     ('solve_zlineq: a',0, GA_ERR)
        if (.not. ga_destroy(g_b)) call errquit
     &     ('solve_zlineq: b',0, GA_ERR)
        if (.not. ga_destroy(g_c)) call errquit
     &     ('solve_zlineq: c',0, GA_ERR)
       return
       end
c -------------- solve_rlineq ---------------- END

       subroutine solve_zlineq(
     &              g_z,   ! in/out: updated solution
     &              g_Az1, ! in    : history of g_Az
     &              g_z1,  ! in    : history of g_z
     &              g_zr1, ! in    : history of g_zr
     &              nsub,  ! in    : subspace length
     &              nvec,  ! in    : increment of subspace
     &              ncomp, ! in    : nr. components
     &              n,     ! in    : nr. elements per comp.
     &              iter,  ! in    : iteration nr.
     &              debug1)! in   : =.true. show debug printouts
c
c Author : Fredy W. Aquino, Northwestern University
c Purpose: Solve complex linear equation using 'history' GA arrays
c Date   : 03-15-12
c
c    dim(g_Az1)=dim(z1)=(ncomp*n,maxsub)
c    dim(g_zr1)=(ncomp*n,nvec)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
        integer nsub,nvec,ncomp,n,n1,iter,ipm
        integer g_a,g_b,g_c,
     &          g_z1,g_Az1,g_zr1,
     &          g_z(ncomp)
        logical debug1
        double complex one_cmplx,mone_cmplx,zero_cmplx
        external ga_svd_solve_seq_cmplx,update_g_z
c       Form and solve the subspace equations using SVD in order
c       to manage near linear dependence in the subspace.
c  
      one_cmplx =dcmplx( 1.0d0,0.0d0)
      mone_cmplx=dcmplx(-1.0d0,0.0d0)
      zero_cmplx=dcmplx( 0.0d0,0.0d0)
        n1=ncomp*n     
        if (.not. ga_create(MT_DCPL, nsub, nsub,
     &     'solve_zlineq: A',0,0,g_a))
     $     call errquit('solve_zlineq: failed allocating g_a', 
     &                  nsub,GA_ERR)
        if (.not. ga_create(MT_DCPL, nsub, nvec,
     &     'solve_zlineq: B',0,0,g_b))
     $     call errquit('solve_zlineq: failed allocating g_b', 
     &                  nsub,GA_ERR)
        if (.not. ga_create(MT_DCPL, nsub, nvec,
     &     'solve_zlineq: C',0,0,g_c))
     $     call errquit('solve_zlineq: failed allocating g_c', 
     &                  nsub,GA_ERR)
        call ga_zero(g_a)
        call ga_zero(g_b)
        call ga_zero(g_c)
        if (debug1) then
        if (ga_nodeid().eq.0) 
     &   write(*,*) '-------BEF:g_z1(',iter,')-------START'
        call ga_print(g_z1)
        if (ga_nodeid().eq.0) 
     &   write(*,*) '-------BEF:g_z1(',iter,')-------END'
        endif ! end-if-debug1
        call get_cconjugate(g_z1) ! out: complex-conjugate of g_z1
        call ga_zgemm('t','n',nsub,nsub,n1,one_cmplx,
     &                g_z1,g_Az1,zero_cmplx,g_a)
        call ga_zgemm('t','n',nsub,nvec,n1,one_cmplx,
     &                g_z1,g_zr1,zero_cmplx,g_b)
        if (debug1) then
        if (ga_nodeid().eq.0)
     &   write(*,*) '-------g_a--------START'
        call ga_print(g_a)
        if (ga_nodeid().eq.0)
     &   write(*,*) '-------g_a--------END'
        if (ga_nodeid().eq.0)
     &   write(*,*) '-------g_b--------START'
        call ga_print(g_b)
        if (ga_nodeid().eq.0)
     &   write(*,*) '-------g_b--------END'
        endif ! end-if-debug

        call ga_svd_solve_seq_cmplx(g_a,g_b,g_c,1d-14) 

        if (debug1) then
        if (ga_nodeid().eq.0)
     &   write(*,*) '-------g_c(',iter,')--------START'
        call ga_print(g_c)
        if (ga_nodeid().eq.0)
     &   write(*,*) '-------g_c(',iter,')--------END'
        endif ! end-if-debug1

        call ga_zgemm('n','n',n1,nvec,nsub,mone_cmplx,
     &                g_Az1,g_c,one_cmplx,g_zr1)

        if (debug1) then
          if (ga_nodeid().eq.0) then
           write(*,10) iter
   10       format('---------g_zr1-1(',i3,')-----START')
          endif
          call ga_print(g_zr1)
          if (ga_nodeid().eq.0) then
           write(*,6) iter
    6      format('---------g_zr1-1(',i3,')-----END')
          endif
        endif ! end-if-debug1

        call update_g_z(g_z,  ! in/ou: solution updated
     &                  g_zr1,! in   : added to g_z
     &                  ncomp,! in   : nr. components
     &                  nvec, ! in   : (x,y,z)
     &                  n)    ! in   : vector length
        if (debug1) then
        do ipm=1,ncomp
          if (ga_nodeid().eq.0) then
           write(*,2) ipm,iter
    2      format('---------g_z-1(',i3,',',i3,')-----START')
          endif
          call ga_print(g_z(ipm))
          if (ga_nodeid().eq.0) then
           write(*,3) ipm,iter
    3      format('---------g_z-1(',i3,',',i3,')-----END')
          endif
        enddo ! end-loop-ipm
        endif ! end-if-debug1
        call ga_zero(g_zr1)
        call get_cconjugate(g_z1) ! put back g_z1
        call ga_zgemm('n','n',n1,nvec,nsub,one_cmplx,
     &                g_z1,g_c,zero_cmplx,g_zr1)

        if (debug1) then
          if (ga_nodeid().eq.0) then
           write(*,7) iter
    7      format('---------g_zr1-2(',i3,')-----START')
          endif
          call ga_print(g_zr1)
          if (ga_nodeid().eq.0) then
           write(*,9) iter
    9      format('---------g_zr1-2(',i3,')-----END')
          endif
        endif ! end-if-debug1

        call update_g_z(g_z,  ! in/ou: solution updated
     &                  g_zr1,! in   : added to g_z
     &                  ncomp,! in   : nr. components
     &                  nvec, ! in   : (x,y,z)
     &                  n)    ! in   : vector length

        if (debug1) then
        do ipm=1,ncomp
          if (ga_nodeid().eq.0) then
           write(*,4) ipm,iter
    4      format('---------g_z-2(',i3,',',i3,')-----START')
          endif
          call ga_print(g_z(ipm))
          if (ga_nodeid().eq.0) then
           write(*,5) ipm,iter
    5      format('---------g_z-2(',i3,',',i3,')-----END')
          endif
        enddo ! end-loop-ipm
        endif ! end-if-debug1

        if (.not. ga_destroy(g_a)) call errquit
     &     ('solve_zlineq: a',0, GA_ERR)
        if (.not. ga_destroy(g_b)) call errquit
     &     ('solve_zlineq: b',0, GA_ERR)
        if (.not. ga_destroy(g_c)) call errquit
     &     ('solve_zlineq: c',0, GA_ERR)
       return
       end
c +++++++++++++++++++++++++++++++++++++++++++++++++++
c +++++++++ FA-12-06-13: fix-KAIN +++++++++++++ START
c +++++++++++++++++++++++++++++++++++++++++++++++++++
       subroutine getdiffs_Az1z1(
     &              g_Az1,    ! in/out: history of g_Az
     &              g_z1,     ! in/out: history of g_z adding a block g_z1
     &              nsub,     ! in    : subspace length
     &              nvec,     ! in    : increment of subspace
     &              ncomp,    ! in    : nr. components
     &              n3,       ! in    : nr. elements per comp.
     &              op)       ! in    : = -1   DO differences
                              !         = +1 UNDO differences
c Purpose: Compute KAIN differences in (g_Az1,g_z1)_k k=0,...,n-1
c          using (g_z1)_n
c Note.- n1=nsub/3  should be pointing to n+1 
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
        integer g_z1,g_Az1,op,
     &          nsub,nvec,ncomp,
     &          n1,n2,n3,i
        integer p1,p2,q1,q2
        double complex op_cmplx,one_cmplx

        if (op.ne.-1.and.op.ne.1) then
         call errquit(
     &     'getdiffs_Az1z1: failed op ne +1 or -1')
        endif
        if      (op.eq.-1) then
         op_cmplx=dcmplx(-1.0d0,0.0d0)
        else if (op.eq.1) then
         op_cmplx=dcmplx(+1.0d0,0.0d0)
        endif
        if (nsub .lt. 6) then
         call errquit('getdiffs_Az1z1: failed nsub lt 6') 
        endif
        one_cmplx =dcmplx( 1.0d0,0.0d0)
        n1=ncomp*n3    
        n2=nsub/nvec+1  ! = n+1  CONDITION: nsub>=6
c --- Compute differences up to n1-2=n-1 -- START
       do i=1,n2-2
         p1=(i-1)*nvec+1
         p2=p1+nvec-1
         q1=nsub-nvec+1
         q2=q1+nvec-1
         call ga_add_patch( op_cmplx, g_z1,1,n1,q1,q2, 
     $                     one_cmplx, g_z1,1,n1,p1,p2,
     $                                g_z1,1,n1,p1,p2)
         call ga_add_patch( op_cmplx,g_Az1,1,n1,q1,q2, 
     $                     one_cmplx,g_Az1,1,n1,p1,p2,
     $                               g_Az1,1,n1,p1,p2)
       enddo ! end-loop-i
c --- Compute differences up to n1-2=n-1 -- END       
       return
       end

       subroutine solve_zlineq1_fixed(
     &              g_Az1,    ! in    : history of g_Az
     &              g_z1,     ! in/out: history of g_z adding a block g_z1
     &              g_zr1,    ! in    : g_zr
     &              nsub,     ! in    : subspace length
     &              nvec,     ! in    : increment of subspace
     &              ncomp,    ! in    : nr. components
     &              n,        ! in    : nr. elements per comp.
     &              iter,     ! in    : iteration nr.
     &              checkorth,! in    : =1 display: z1^t*(z1c) (goes to zero is it converges)
     &              debug1)   ! in    : =.true. show debug printouts
c
c Author : Fredy W. Aquino, Northwestern University
c Purpose: Solve complex linear equation using history GA arrays 
c          and reducing memory cost
c Date   : 03-15-12
c
c    dim(g_Az1)=dim(z1)=(ncomp*n,maxsub)
c    dim(g_zr1)=(ncomp*n,nvec)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
        integer nsub,nsub1,nvec,ncomp,
     &          p1,p2,q1,q2,
     &          n,n1,iter,ipm,checkorth
        integer g_a,g_b,g_c,
     &          g_z1,g_Az1,g_zr1
        logical debug1
        double complex one_cmplx,mone_cmplx,zero_cmplx
        external ga_svd_solve_seq_cmplx,
     &           update_g_z,
     &           toview_orthz1c,
     &           getdiffs_Az1z1
c       Form and solve the subspace equations using SVD in order
c       to manage near linear dependence in the subspace.
c  
        one_cmplx =dcmplx( 1.0d0,0.0d0)
        mone_cmplx=dcmplx(-1.0d0,0.0d0)
        zero_cmplx=dcmplx( 0.0d0,0.0d0)
        n1=ncomp*n    
        nsub1=nsub-nvec ! nvec=3
        if (.not. ga_create(MT_DCPL, nsub1, nsub1,
     &     'solve_zlineq: A',0,0,g_a))
     $     call errquit('solve_zlineq: failed allocating g_a', 
     &                  nsub1,GA_ERR)
        if (.not. ga_create(MT_DCPL, nsub1, nvec,
     &     'solve_zlineq: B',0,0,g_b))
     $     call errquit('solve_zlineq: failed allocating g_b', 
     &                  nsub1,GA_ERR)
        if (.not. ga_create(MT_DCPL, nsub1, nvec,
     &     'solve_zlineq: C',0,0,g_c))
     $     call errquit('solve_zlineq: failed allocating g_c', 
     &                  nsub1,GA_ERR)
        call ga_zero(g_a)
        call ga_zero(g_b)
        call ga_zero(g_c)
        if (debug1) then
         if (ga_nodeid().eq.0) 
     &   write(*,*) '-------BEF:g_z1(',iter,')-------START'
         call ga_print(g_z1)
         if (ga_nodeid().eq.0) 
     &   write(*,*) '-------BEF:g_z1(',iter,')-------END'
         if (ga_nodeid().eq.0) 
     &    write(*,*) '-------BEF:g_Az1(',iter,')-------START'
         call ga_print(g_Az1)
         if (ga_nodeid().eq.0) 
     &    write(*,*) '-------BEF:g_Az1(',iter,')-------END'
         if (ga_nodeid().eq.0) 
     &    write(*,*) '-------BEF:g_zr1(',iter,')-------START'
         call ga_print(g_zr1)
         if (ga_nodeid().eq.0) 
     &    write(*,*) '-------BEF:g_zr1(',iter,')-------END'
        endif ! end-if-debug1
        call getdiffs_Az1z1(
     &              g_Az1,    ! in/out: history of g_Az
     &              g_z1,     ! in/out: history of g_z adding a block g_z1
     &              nsub,     ! in    : subspace length
     &              nvec,     ! in    : increment of subspace
     &              ncomp,    ! in    : nr. components
     &              n,        ! in    : nr. elements per comp.
     &              -1)       ! in    : = -1   DO differences
                              !         = +1 UNDO differences

        call get_cconjugate(g_z1) ! out: complex-conjugate of g_z1
        call ga_zgemm('t','n',nsub1,nsub1,n1,one_cmplx,
     &                g_z1,g_Az1,zero_cmplx,g_a)
        call ga_zgemm('t','n',nsub1,nvec,n1,one_cmplx,
     &                g_z1,g_zr1,zero_cmplx,g_b)

        if (debug1) then
        if (ga_nodeid().eq.0)
     &   write(*,10) iter
   10    format('-------g_a(',i4,')--------START')
        call ga_print(g_a)
        if (ga_nodeid().eq.0)
     &   write(*,*) iter
   11    format('-------g_a(',i4,')--------END')
        if (ga_nodeid().eq.0)
     &   write(*,12) iter
   12    format('-------g_b(',i4,')--------START')
        call ga_print(g_b)
        if (ga_nodeid().eq.0)
     &   write(*,13) iter
   13    format('-------g_b(',i4,')--------END')
        endif ! end-if-debug

        call ga_svd_solve_seq_cmplx(g_a,g_b,g_c,1d-14)

        if (debug1) then
        if (ga_nodeid().eq.0)
     &   write(*,14) iter
   14    format('-------g_c(',i4,')--------START')
        call ga_print(g_c)
        if (ga_nodeid().eq.0)
     &   write(*,15) iter
   15    format('-------g_c(',i4,')--------END')
        endif ! end-if-debug1

        call ga_zgemm('n','n',n1,nvec,nsub1,
     &                mone_cmplx,g_Az1,g_c,
     &                one_cmplx ,g_zr1)
c 000000000000 check-orthonogality-1 000000000000 START
c Compute: g_z1^t . (g_zr1-Az1c)=0 by construction
c 000000000000 check-orthonogality-1 000000000000 END
        call get_cconjugate(g_z1) ! put back g_z1 as it was

        if (debug1) then
          if (ga_nodeid().eq.0) then
           write(*,20) iter
   20      format('---------g_zr1-1(',i3,')-----START')
          endif
          call ga_print(g_zr1)
          if (ga_nodeid().eq.0) then
           write(*,6) iter
    6      format('---------g_zr1-1(',i3,')-----END')
          endif
        endif ! end-if-debug1

c --- construct new (n1,nvec) block-in g_z1------ START
c Note.- It uses previous (n1,nvec) block in g_z1 and
c        g_zr1(=g_Az1 * g_c)
         p1=nsub-nvec+1
         p2=p1+nvec-1
         q1=p1+nvec
         q2=p2+nvec
         call ga_add_patch(one_cmplx,g_zr1,1,n1,1,nvec, 
     $                     one_cmplx,g_z1 ,1,n1,p1,p2,
     $                               g_z1 ,1,n1,q1,q2)
c --- construct new (n1,nvec) block-in g_z1------ END
        if (debug1) then
          if (ga_nodeid().eq.0) then
           write(*,2) iter
    2      format('---------g_z1-1(',i3,')-----START')
          endif
          call ga_print(g_z1)
          if (ga_nodeid().eq.0) then
           write(*,3) iter
    3      format('---------g_z1-1(',i3,')-----END')
          endif
        endif ! end-if-debug1
        call ga_zero(g_zr1)
        call ga_zgemm('n','n',n1,nvec,nsub1,one_cmplx,
     &                g_z1,g_c,zero_cmplx,g_zr1)
c 000000000000 check-orthonogality-2 000000000000 START
c Compute: g_z1^t . z1c:
        if (debug1) then
         if (ga_nodeid().eq.0)
     &    write(*,*) '-------z1c(',iter,')--------START'
         call ga_print(g_zr1)
         if (ga_nodeid().eq.0)
     &    write(*,*) '-------z1c(',iter,')--------END'
        endif ! end-if-debug 

        call get_cconjugate(g_z1) ! conjugate for next op
        call ga_zero(g_b)
        call ga_zgemm('t','n',nsub1,nvec,n1,one_cmplx,
     &                g_z1,g_zr1,zero_cmplx,g_b)   
        call get_cconjugate(g_z1) ! put back as it was
        if (checkorth.eq.1) then ! display z1^t*(z1c)
          call toview_orthz1c(
     &                    g_b,
     &                  nsub1,
     &                   nvec,
     &                   iter)
        endif
c        endif ! end-if-debug1
c 000000000000 check-orthonogality-2 000000000000 END
        if (debug1) then
          if (ga_nodeid().eq.0) then
           write(*,7) iter
    7      format('---------g_zr1-2(',i3,')-----START')
          endif
          call ga_print(g_zr1)
          if (ga_nodeid().eq.0) then
           write(*,9) iter
    9      format('---------g_zr1-2(',i3,')-----END')
          endif
        endif ! end-if-debug1

c --- update new (n1,nvec) block-in g_z1------ START
         p1=nsub-nvec+1
         p2=p1+nvec-1
         q1=p1+nvec
         q2=p2+nvec
         call ga_add_patch(one_cmplx,g_zr1,1,n1,1,nvec, 
     $                     one_cmplx,g_z1 ,1,n1,q1,q2,
     $                               g_z1 ,1,n1,q1,q2)
c --- update new (n1,nvec) block-in g_z1------ END

        call getdiffs_Az1z1(
     &              g_Az1,    ! in/out: history of g_Az
     &              g_z1,     ! in/out: history of g_z adding a block g_z1
     &              nsub,     ! in    : subspace length
     &              nvec,     ! in    : increment of subspace
     &              ncomp,    ! in    : nr. components
     &              n,        ! in    : nr. elements per comp.
     &              1)        ! in    : = -1   DO differences
                              !         = +1 UNDO differences
        if (debug1) then
          if (ga_nodeid().eq.0) then
           write(*,4) iter
    4      format('---------g_z1-2(',i3,')-----START')
          endif
          call ga_print(g_z1)
          if (ga_nodeid().eq.0) then
           write(*,5) iter
    5      format('---------g_z1-2(',i3,')-----END')
          endif
        endif ! end-if-debug1

        if (.not. ga_destroy(g_a)) call errquit
     &     ('solve_zlineq: a',0, GA_ERR)
        if (.not. ga_destroy(g_b)) call errquit
     &     ('solve_zlineq: b',0, GA_ERR)
        if (.not. ga_destroy(g_c)) call errquit
     &     ('solve_zlineq: c',0, GA_ERR)
       return
       end

       subroutine solve_zlineq1(
     &              g_Az1,    ! in    : history of g_Az
     &              g_z1,     ! in/out: history of g_z adding a block g_z1
     &              g_zr1,    ! in    : g_zr
     &              nsub,     ! in    : subspace length
     &              nvec,     ! in    : increment of subspace
     &              ncomp,    ! in    : nr. components
     &              n,        ! in    : nr. elements per comp.
     &              iter,     ! in    : iteration nr.
     &              checkorth,! in    : =1 display: z1^t*(z1c) (goes to zero is it converges)
     &              debug1)   ! in    : =.true. show debug printouts
c
c Author : Fredy W. Aquino, Northwestern University
c Purpose: Solve complex linear equation using history GA arrays 
c          and reducing memory cost
c Date   : 03-15-12
c
c    dim(g_Az1)=dim(z1)=(ncomp*n,maxsub)
c    dim(g_zr1)=(ncomp*n,nvec)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
        integer nsub,nvec,ncomp,
     &          p1,p2,q1,q2,
     &          n,n1,iter,ipm,checkorth
        integer g_a,g_b,g_c,
     &          g_z1,g_Az1,g_zr1
        logical debug1
        double complex one_cmplx,mone_cmplx,zero_cmplx
        external ga_svd_solve_seq_cmplx,
     &           update_g_z,
     &           toview_orthz1c
c       Form and solve the subspace equations using SVD in order
c       to manage near linear dependence in the subspace.
c  
        one_cmplx =dcmplx( 1.0d0,0.0d0)
        mone_cmplx=dcmplx(-1.0d0,0.0d0)
        zero_cmplx=dcmplx( 0.0d0,0.0d0)
        n1=ncomp*n    
        if (.not. ga_create(MT_DCPL, nsub, nsub,
     &     'solve_zlineq: A',0,0,g_a))
     $     call errquit('solve_zlineq: failed allocating g_a', 
     &                  nsub,GA_ERR)
        if (.not. ga_create(MT_DCPL, nsub, nvec,
     &     'solve_zlineq: B',0,0,g_b))
     $     call errquit('solve_zlineq: failed allocating g_b', 
     &                  nsub,GA_ERR)
        if (.not. ga_create(MT_DCPL, nsub, nvec,
     &     'solve_zlineq: C',0,0,g_c))
     $     call errquit('solve_zlineq: failed allocating g_c', 
     &                  nsub,GA_ERR)
        call ga_zero(g_a)
        call ga_zero(g_b)
        call ga_zero(g_c)
        if (debug1) then
         if (ga_nodeid().eq.0) 
     &   write(*,*) '-------BEF:g_z1(',iter,')-------START'
         call ga_print(g_z1)
         if (ga_nodeid().eq.0) 
     &   write(*,*) '-------BEF:g_z1(',iter,')-------END'
         if (ga_nodeid().eq.0) 
     &    write(*,*) '-------BEF:g_Az1(',iter,')-------START'
         call ga_print(g_Az1)
         if (ga_nodeid().eq.0) 
     &    write(*,*) '-------BEF:g_Az1(',iter,')-------END'
         if (ga_nodeid().eq.0) 
     &    write(*,*) '-------BEF:g_zr1(',iter,')-------START'
         call ga_print(g_zr1)
         if (ga_nodeid().eq.0) 
     &    write(*,*) '-------BEF:g_zr1(',iter,')-------END'
        endif ! end-if-debug1

        call get_cconjugate(g_z1) ! out: complex-conjugate of g_z1
        call ga_zgemm('t','n',nsub,nsub,n1,one_cmplx,
     &                g_z1,g_Az1,zero_cmplx,g_a)
        call ga_zgemm('t','n',nsub,nvec,n1,one_cmplx,
     &                g_z1,g_zr1,zero_cmplx,g_b)
        if (debug1) then
        if (ga_nodeid().eq.0)
     &   write(*,10) iter
   10    format('-------g_a(',i4,')--------START')
        call ga_print(g_a)
        if (ga_nodeid().eq.0)
     &   write(*,*) iter
   11    format('-------g_a(',i4,')--------END')
        if (ga_nodeid().eq.0)
     &   write(*,12) iter
   12    format('-------g_b(',i4,')--------START')
        call ga_print(g_b)
        if (ga_nodeid().eq.0)
     &   write(*,13) iter
   13    format('-------g_b(',i4,')--------END')
        endif ! end-if-debug

        call ga_svd_solve_seq_cmplx(g_a,g_b,g_c,1d-14)

        if (debug1) then
        if (ga_nodeid().eq.0)
     &   write(*,14) iter
   14    format('-------g_c(',i4,')--------START')
        call ga_print(g_c)
        if (ga_nodeid().eq.0)
     &   write(*,15) iter
   15    format('-------g_c(',i4,')--------END')
        endif ! end-if-debug1

        call ga_zgemm('n','n',n1,nvec,nsub,mone_cmplx,
     &                g_Az1,g_c,one_cmplx,g_zr1)
c 000000000000 check-orthonogality-1 000000000000 START
c Compute: g_z1^t . (g_zr1-Az1c)=0 by construction
c 000000000000 check-orthonogality-1 000000000000 END
        call get_cconjugate(g_z1) ! put back g_z1 as it was

        if (debug1) then
          if (ga_nodeid().eq.0) then
           write(*,20) iter
   20       format('---------g_zr1-1(',i3,')-----START')
          endif
          call ga_print(g_zr1)
          if (ga_nodeid().eq.0) then
           write(*,6) iter
    6      format('---------g_zr1-1(',i3,')-----END')
          endif
        endif ! end-if-debug1

c --- construct new (n1,nvec) block-in g_z1------ START
c Note.- It uses previous (n1,nvec) block in g_z1 and
c        g_zr1(=g_Az1 * g_c)
         p1=nsub-nvec+1
         p2=nsub-nvec+nvec
         q1=p1+nvec
         q2=p2+nvec
         call ga_add_patch(one_cmplx,g_zr1,1,n1,1,nvec, 
     $                     one_cmplx,g_z1 ,1,n1,p1,p2,
     $                               g_z1 ,1,n1,q1,q2)
c --- construct new (n1,nvec) block-in g_z1------ END
        if (debug1) then
          if (ga_nodeid().eq.0) then
           write(*,2) iter
    2      format('---------g_z1-1(',i3,')-----START')
          endif
          call ga_print(g_z1)
          if (ga_nodeid().eq.0) then
           write(*,3) iter
    3      format('---------g_z1-1(',i3,')-----END')
          endif
        endif ! end-if-debug1
        call ga_zero(g_zr1)
        call ga_zgemm('n','n',n1,nvec,nsub,one_cmplx,
     &                g_z1,g_c,zero_cmplx,g_zr1)
c 000000000000 check-orthonogality-2 000000000000 START
c Compute: g_z1^t . z1c:
        if (debug1) then
         if (ga_nodeid().eq.0)
     &    write(*,*) '-------z1c(',iter,')--------START'
         call ga_print(g_zr1)
         if (ga_nodeid().eq.0)
     &    write(*,*) '-------z1c(',iter,')--------END'
        endif ! end-if-debug 

        call get_cconjugate(g_z1) ! conjugate for next op
        call ga_zero(g_b)
        call ga_zgemm('t','n',nsub,nvec,n1,one_cmplx,
     &                g_z1,g_zr1,zero_cmplx,g_b)   
        call get_cconjugate(g_z1) ! put back as it was
        if (checkorth.eq.1) then ! display z1^t*(z1c)
          call toview_orthz1c(
     &                    g_b,
     &                   nsub,
     &                   nvec,
     &                   iter)
        endif
c        endif ! end-if-debug1
c 000000000000 check-orthonogality-2 000000000000 END
        if (debug1) then
          if (ga_nodeid().eq.0) then
           write(*,7) iter
    7      format('---------g_zr1-2(',i3,')-----START')
          endif
          call ga_print(g_zr1)
          if (ga_nodeid().eq.0) then
           write(*,9) iter
    9      format('---------g_zr1-2(',i3,')-----END')
          endif
        endif ! end-if-debug1

c --- update new (n1,nvec) block-in g_z1------ START
         p1=nsub-nvec+1
         p2=nsub-nvec+nvec
         q1=p1+nvec
         q2=p2+nvec
         call ga_add_patch(one_cmplx,g_zr1,1,n1,1,nvec, 
     $                     one_cmplx,g_z1 ,1,n1,q1,q2,
     $                               g_z1 ,1,n1,q1,q2)
c --- update new (n1,nvec) block-in g_z1------ END
        if (debug1) then
          if (ga_nodeid().eq.0) then
           write(*,4) iter
    4      format('---------g_z1-2(',i3,')-----START')
          endif
          call ga_print(g_z1)
          if (ga_nodeid().eq.0) then
           write(*,5) iter
    5      format('---------g_z1-2(',i3,')-----END')
          endif
        endif ! end-if-debug1

        if (.not. ga_destroy(g_a)) call errquit
     &     ('solve_zlineq: a',0, GA_ERR)
        if (.not. ga_destroy(g_b)) call errquit
     &     ('solve_zlineq: b',0, GA_ERR)
        if (.not. ga_destroy(g_c)) call errquit
     &     ('solve_zlineq: c',0, GA_ERR)
       return
       end

       subroutine toview_orthz1c(
     &                 g_b,
     &                 nsub,
     &                 nvec,
     &                 iter)
c
c Author : Fredy W. Aquino, Northwestern University
c Purpose: Visualize g_b with 8 decimals.
c Date   : 03-15-12

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
c       Note.- g_z= z1^t . (z1c) it is becoming zero
c              as the iteration increases
        integer nsub,nvec,g_b,
     &          l_z,k_z,
     &          i,j,iter
        double precision valre,valim
      if (.not.MA_Push_Get(mt_dcpl,nvec,'g_bre',l_z,k_z))
     &     call errquit('toview_orthz1c: cannot allocate zre',
     &                  nvec, MA_ERR)    
        if (ga_nodeid().eq.0)
     &   write(*,1) iter
    1    format('-------z1^t*(z1c)(',i4,')--------START')
        do i=1,nsub
         call ga_get(g_b,i,i,1,nvec,dcpl_mb(k_z),1)
         if (ga_nodeid().eq.0) then
          write(*,'(i8,":",3(" (",f14.8,",",f14.8,") "))') 
     &      i,(dreal(dcpl_mb(k_z+j-1)),dimag(dcpl_mb(k_z+j-1)),j=1,nvec)
         endif
        enddo ! end-loop-idata  
        if (ga_nodeid().eq.0)
     &   write(*,2) iter
    2    format('-------z1^t*(z1c)(',i4,')--------END')
      if (.not.ma_pop_stack(l_z))
     $  call errquit('toview_orthz1c: pop problem with l_zim',
     &               555,MA_ERR)
       return
       end

       subroutine toview_orthz1c_short(
     &                 g_b,
     &                 nblocks,
     &                 nvec,
     &                 iter)
c
c Author : Fredy W. Aquino, Northwestern University
c Purpose: Visualize g_b with 8 decimals.
c Date   : 03-15-12

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
c       Note.- g_z= z1^t . (z1c) it is becoming zero
c              as the iteration increases
        integer nblocks,nvec,g_b,
     &          l_z,k_z,
     &          i,j,iter
        double precision valre,valim
      if (.not.MA_Push_Get(mt_dcpl,nvec,'g_bre',l_z,k_z))
     &     call errquit('toview_orthz1c: cannot allocate zre',
     &                  nvec, MA_ERR)    
          if (ga_nodeid().eq.0) then
           write(*,16) iter
   16      format('---------z1^t*(z1c)(',i3,')-----START')
          endif
        do i=1,nblocks
         call ga_get(g_b,i,i,1,nvec,dcpl_mb(k_z),1)
         if (ga_nodeid().eq.0) then
          write(*,'(i8,":",3(" (",f14.8,",",f14.8,") "))') 
     &      i,(dreal(dcpl_mb(k_z+j-1)),dimag(dcpl_mb(k_z+j-1)),j=1,nvec)
         endif
        enddo ! end-loop-idata  
          if (ga_nodeid().eq.0) then
           write(*,17) iter
   17      format('---------z1^t*(z1c)(',i3,')-----END')
          endif
      if (.not.ma_pop_stack(l_z))
     $  call errquit('toview_orthz1c: pop problem with l_zim',
     &               555,MA_ERR)
       return
       end

      subroutine ga_svd_solve_seq_cmplx(
     &                          g_a, ! in : a of a x = b
     &                          g_b, ! in : b of a x = b
     &                          g_x, ! out: x of a x = b
     &                          tol) ! in : tolerance
c
c Author: Fredy W. Aquino, Northwestern University
c Date  : 03-15-12
c
c    Note.- Adapted to do complex calc. from ga_svd_solve_seq()
c           located in ga_it2.F

      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
      integer g_a, g_b, g_x
      double precision tol
c
c     Solve for X from the linear equations
c
c     A*X = B
c
c     A(m,n)*X(n,nvec) = B(m,nvec)
c
c     Where A is a general real matrix (not necessarily square, or
c     symmetric, or full rank) and X and B are matrices with one or more
c     columns representing the solutions and right hand sides.  Singular
c     values of A less than tol are neglected.  X is returned.
c
c     If the SVD of A is U*values*VT, then the solution
c     is of the form
c
c     V*(1/values)*UT*B
c
c     where the reciprocal of values less than tol are neglected.
c     
      integer m,n,nn,type,nvec,nsing,i,
     &        l_val,k_val,
     &        g_u,g_vt,g_tmp
      logical oprint
      double complex one_cmplx,zero_cmplx
      external ga_svd_seq_cmplx,ga_scale_lh_cmplx,
     &         get_cconjugate

      oprint = util_print('debug svdsolve', print_high) .and.
     $     ga_nodeid().eq.0

      call ga_inquire(g_a, type, m, n)
      call ga_inquire(g_b, type, nn, nvec)

      if (nn .ne. n) call errquit('gasvdsol: b does not conform',nn,
     &       GA_ERR)
      nsing = min(m,n)
      if (.not. ma_push_get(MT_DBL, nsing, 'gasvdsol', 
     &                     l_val, k_val))
     $     call errquit('gasvdsol: val',nsing, MA_ERR)
      if (.not. ga_create(MT_DCPL,m,nsing,'gasvd',0,0,g_u))
     $     call errquit('gasvdsol: u',m*nsing, GA_ERR)
      if (.not. ga_create(MT_DCPL,nsing,n,'gasvd',0,0,g_vt))
     $     call errquit('gasvdsol: u',nsing*n, GA_ERR)
      if (.not. ga_create(MT_DCPL,nsing,nvec,'gasvd',0,0,g_tmp))
     $     call errquit('gasvdsol: tmp',nsing*nvec, GA_ERR)
      call ga_zero(g_tmp)

      call ga_svd_seq_cmplx(g_a,g_u,g_vt,dbl_mb(k_val))

      do i = 0, nsing-1
         if (dbl_mb(k_val+i) .lt. tol) then
            if (ga_nodeid() .eq. 0 .and. oprint) then
              write(6,*) ' neglecting ', i+1, dbl_mb(k_val+i)
            endif
            dbl_mb(k_val+i) = 0.0d0
         else
            dbl_mb(k_val+i) = 1.0d0/dbl_mb(k_val+i)
         end if
      end do
c Ax=b from SVD: A= U w V^t -> x=A^{-1}b  A^{-1}=V w^{-1} U^t 
c Note.- Using property: U^{-1}=U^t V^{-1}=V^t
c        Using property: (AB)^{-1}=B^{-1}A^{-1} (FA-04-06-12)
      one_cmplx =dcmplx(1.0d0,0.0d0)
      zero_cmplx=dcmplx(0.0d0,0.0d0)
c Note.- U^t --> U^H (complex conjugate + transposed)
      call get_cconjugate(g_u) 
      call ga_zgemm('t','n',nsing,nvec,m,one_cmplx,g_u,g_b,
     &              zero_cmplx,g_tmp)                 ! U^t b          -> g_tmp
      call ga_scale_lh_cmplx(g_tmp,dbl_mb(k_val))     ! w^{-1} U^t b   -> g_tmp
      call ga_zero(g_x)
c Note.-since g_vt = V^H to get V from V^H : V=complex conjugate + transpose of V^H
      call get_cconjugate(g_vt) 
      call ga_zgemm('t','n',n,nvec,nsing,one_cmplx,g_vt,g_tmp,
     &              zero_cmplx,g_x)                   ! V w^{-1} U^t b -> g_x
      if (.not. ga_destroy(g_tmp)) call errquit('gasvdsol: des',1,
     &       GA_ERR)
      if (.not. ga_destroy(g_u)) call errquit('gasvdsol: des',2,
     &       GA_ERR)
      if (.not. ga_destroy(g_vt)) call errquit('gasvdsol: des',3,
     &       GA_ERR)
      if (.not. ma_pop_stack(l_val)) call errquit('gasvdsol: pop',4,
     &       GA_ERR)

      end

      subroutine ga_svd_seq_cmplx(g_a, g_u, g_vt, values)
c Author: Fredy W. Aquino, Northwestern University
c Date  : 03-15-12
c
c Note.- Adapted from ga_svd_seq
c        to handle complex g_a FA-04-08-12

      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer g_a, g_u, g_vt
      double precision values(*)
      external zgesvd
c     
c     Perform SVD on rectangular matrix
c
c     nsing = min(n,m)
c     g_a(m,n)      --- input matrix
c     g_u(m,nsing)  --- left singular vectors (output)
c     g_vt(nsing,n) --- right singular vectors transposed (output)
c     values(nsing) --- singular values (output)
c
c     A = U*values*VT
c
c     A possible parallel algorithm is to diagonalize ATA to get
c     V and AAT to get U --- both have values**2 as eigenvalues.

      integer n, m, type,n_rwork,
     &        l_a,k_a,
     &        l_u,k_u,
     &        l_vt,k_vt, 
     $        l_work,k_work,
     $        l_rwork,k_rwork,
     &        lwork, info, nsing

      integer i ! for debugging purpose
      call ga_inquire(g_a, type, m, n)
      nsing = min(m,n)
      call ga_sync() ! FA-added 07-01-12
c 000000000000000000000000000000000000000000000000000000000
      if (ga_nodeid() .eq. 0) then ! 000000 node-0 00 START
         lwork = 10*max(m,n)

         if (.not. ma_push_get(MT_DCPL, m*n, 'gasvd1',l_a,k_a))
     $        call errquit('gasvd: a',m*n, MA_ERR)
         if (.not. ma_push_get(MT_DCPL, m*nsing, 'gasvd2',l_u,k_u))
     $        call errquit('gasvd: u',m*nsing, MA_ERR)
         if (.not. ma_push_get(MT_DCPL, nsing*n, 'gasvd3',l_vt,k_vt))
     $        call errquit('gasvd: vt',nsing*n, MA_ERR)
         if (.not. ma_push_get(MT_DCPL, lwork, 'gasvd4',l_work,k_work))
     $        call errquit('gasvd: work',lwork, MA_ERR)
         n_rwork=5*min(m,n)
         if (.not. ma_push_get(MT_DBL, n_rwork, 'gasvd5', 
     &                         l_rwork, k_rwork))
     $        call errquit('gasvd: rwork',n_rwork, MA_ERR)

         call ga_get(g_a, 1, m, 1, n, dcpl_mb(k_a), m)
         call zgesvd('s','s',m,n,dcpl_mb(k_a),m,values,
     $               dcpl_mb(k_u),m,dcpl_mb(k_vt),nsing,
     $               dcpl_mb(k_work),lwork,dbl_mb(k_rwork),info)
         if (info .ne. 0) then
          call errquit('gasvd: failed', info, MEM_ERR)
         endif
         call ga_put(g_u,  1, n,     1, nsing, dcpl_mb(k_u),  n)
         call ga_put(g_vt, 1, nsing, 1, m,     dcpl_mb(k_vt), n)
c Deallocating l_a and (l_u,l_vt,l_work,l_rwork)
         if (.not. ma_chop_stack(l_a)) call errquit('gasvd a',0,
     &       MA_ERR)
      end if ! 000000000000000000000000000000 node-0 00 END
c 000000000000000000000000000000000000000000000000000000000
      call ga_sync()
      call ga_brdcst(1,values,n*8,0)
      call ga_sync()    
      end

      subroutine ga_scale_lh_cmplx(g_a,vector)
c
c Author: Fredy W. Aquino, Northwestern University
c Date  : 03-15-12
c
c Note.- Adapted from ga_scale_lh in ga_extra.F
c        to handle complex g_a FA-04-08-12

      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
      integer g_a
      double precision vector(*)
c
c     Do full matrix times diagonal matrix with the diagonal
c     matrix on the left hand side stored as a vector
c     This boils down to a row wise scaling of the g_a
c    
c     do i=1,nbf
c        do j=1,nbf
c           g_a(i,j)=g_a(i,j)*vector(i)
c        enddo
c     enddo

      integer ma_type, dim1, dim2, n
      integer i, j, ilo, ihi, jlo, jhi
      double complex number

      call ga_sync()
      call ga_inquire(g_a, ma_type, dim1, dim2)
      if (ma_type.ne.mt_dcpl) then
         write(0,*) ' ma_type ',ma_type,' mt_dcpl ',mt_dcpl
        call errquit
     $     ('ga_scale_lh_cmplx: array is not complex', g_a,0)
      endif
      n = dim1
c
c     Extract and sum the diags local to each process

      call ga_distribution(g_a, ga_nodeid(), ilo, ihi, jlo, jhi)
      if (ilo.gt.0 .and. jlo.gt.0) then
         do i = ilo,ihi
            do j = jlo,jhi
               call ga_get(g_a, i, i, j, j, number, 1)
               number = dcmplx(dreal(number)*vector(i),
     &                         dimag(number)*vector(i))
               call ga_put(g_a,i, i, j, j, number, 1)
            enddo
         enddo
      endif
      call ga_sync()
      end

      subroutine update_g_x(g_x,  ! in/ou: solution updated
     &                      g_xr1,! in   : added to g_z
     &                      ncomp,! in   : nr. components
     &                      nvec, ! in   : (x,y,z)
     &                      n)    ! in   : vector length
c
c Author: Fredy W. Aquino, Northwestern University
c Date  : 03-15-12

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
      integer ipm,ncomp,nvec,n,m1,m2
      integer g_x(ncomp),g_xr1
      double precision one
      one=1.0d0
         m1=1
         m2=n
        do ipm=1,ncomp
         call ga_add_patch(one,g_xr1   ,m1,m2,1,nvec, 
     $                     one,g_x(ipm),1 ,n ,1,nvec,
     $                         g_x(ipm),1 ,n ,1,nvec)
         m1=m1+n
         m2=m2+n
        enddo !end-loop-ipm
      return
      end

      subroutine update_g_z(g_z,  ! in/ou: solution updated
     &                      g_zr1,! in   : added to g_z
     &                      ncomp,! in   : nr. components
     &                      nvec, ! in   : (x,y,z)
     &                      n)    ! in   : vector length
c
c Author: Fredy W. Aquino, Northwestern University
c Date  : 03-15-12

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
      integer ipm,ncomp,nvec,n,m1,m2
      integer g_z(ncomp),g_zr1
      double complex one_cmplx
      one_cmplx=dcmplx(1.0d0,0.0d0)
         m1=1
         m2=n
        do ipm=1,ncomp
         call ga_add_patch(one_cmplx,g_zr1   ,m1,m2,1,nvec, 
     $                     one_cmplx,g_z(ipm),1 ,n ,1,nvec,
     $                               g_z(ipm),1 ,n ,1,nvec)
         m1=m1+n
         m2=m2+n
        enddo !end-loop-ipm
      return
      end

      subroutine get_cconjugate(g_a)
c
c Author: Fredy W. Aquino, Northwestern University
c Date  : 03-15-12

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
      integer g_a,i,j,ilo,ihi,jlo,jhi
      integer dim1,dim2,ma_type
      double complex number
c Purpose: Compute complex conjugate of g_a
      call ga_sync()
      call ga_inquire(g_a, ma_type, dim1, dim2)
      if (ma_type.ne.mt_dcpl) then
         write(0,*) ' ma_type ',ma_type,' mt_dcpl ',mt_dcpl
        call errquit
     $     ('get_cconjugate: array is not complex', g_a,0)
      endif
      call ga_distribution(g_a,ga_nodeid(),ilo,ihi,jlo,jhi)
      if (ilo.gt.0 .and. jlo.gt.0) then
         do i = ilo,ihi
            do j = jlo,jhi
               call ga_get(g_a,i,i,j,j,number,1)
               number = dcmplx( dreal(number),
     &                         -dimag(number))
               call ga_put(g_a,i,i,j,j,number,1)
            enddo
         enddo
      endif
      call ga_sync()
      return
      end

      subroutine get_modulus(g_a)
c
c Author: Fredy W. Aquino, Northwestern University
c Date  : 03-15-12

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
      integer g_a,i,j,ilo,ihi,jlo,jhi
      integer dim1,dim2,ma_type
      double complex number
      double precision mod
c Purpose: Compute complex conjugate of g_a
      call ga_sync()
      call ga_inquire(g_a, ma_type, dim1, dim2)
      if (ma_type.ne.mt_dcpl) then
         write(0,*) ' ma_type ',ma_type,' mt_dcpl ',mt_dcpl
        call errquit
     $     ('get_cconjugate: array is not complex', g_a,0)
      endif
      call ga_distribution(g_a,ga_nodeid(),ilo,ihi,jlo,jhi)
      if (ilo.gt.0 .and. jlo.gt.0) then
         do i = ilo,ihi
            do j = jlo,jhi
               call ga_get(g_a,i,i,j,j,number,1)
               mod=dreal(number)**2+
     &             dimag(number)**2
               if (mod .gt. 0.0d0) then
                 mod=dsqrt(mod)
               else
                 mod=0.0d0
               endif
               number = dcmplx(mod,0.0d0)
               call ga_put(g_a,i,i,j,j,number,1)
            enddo
         enddo
      endif
      call ga_sync()
      return
      end

      subroutine getrmax_z(rmax, ! out: max(all elements g_zr)
     &                     g_zr, ! in : 
     &                     n,    ! in : nr. rows       in g_zr
     &                     nvec, ! in : nr. cols       in g_zr
     &                     ncomp)! in P nr. components in g_zr   
c
c Author: Fredy W. Aquino, Northwestern University
c Date  : 03-15-12
                
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"
      double precision dat,rmax
      integer ncomp,n,nvec,nreim,
     &        i,cc,ipm,g_zr(ncomp),
     &        g_arr ! scratch arr 
      integer l_a,k_a
      external ga_maxelt
      nreim=2
         if (.not. ma_push_get(MT_DCPL, n*nvec,'gasvd',l_a,k_a))
     $        call errquit('getrmax_z: a',n*nvec, MA_ERR)
         if (.not. ga_create(MT_DBL,nreim*n*nvec*ncomp,1,
     &             'getrmax_z: arr',0, 0, g_arr))
     $     call errquit('lkain: failed allocating zb', 1,
     &     GA_ERR)
           call ga_zero(g_arr)
      cc=0
      do ipm=1,ncomp
       call ga_get(g_zr(ipm),1,n,1,nvec,dcpl_mb(k_a),nvec)
       do i=1,n*nvec
        dat=dreal(dcpl_mb(k_a+i-1))
        call ga_put(g_arr,cc,cc,1,1,dat,1)
        cc=cc+1
       enddo ! end-loop-i
       do i=1,n*nvec
        dat=dimag(dcpl_mb(k_a+i-1))
        call ga_put(g_arr,cc,cc,1,1,dat,1)
        cc=cc+1
       enddo ! end-loop-i
      enddo ! end-loop-ipm
      call ga_maxelt(g_arr,rmax)
      if (ga_nodeid().eq.0)
     & write(*,*) 'In getrmax_z: rmax=',rmax

        if (.not. ga_destroy(g_arr)) call errquit
     &     ('getrmax_z: arr',0, GA_ERR)
         if (.not. ma_chop_stack(l_a)) call errquit('getrmax_z ma',0,
     &       MA_ERR)
      return
      end
c -------------- get_precond_rmax_re -------------- START
      subroutine get_precond_rmax_re(
     &                  rmax,    ! out: max(g_r,g_r_im)   
     &                  g_r,     ! in : real part of g_zr
     &                  g_Ax,    ! in : real part of g_Az
     &                  precond, ! in : name of preconditioner routine
     &        converge_precond,  ! in : =.true. prec->max
     &                  omega,   ! in : omega
     &                  ncomp,   ! in : nr. components
     &                  iter,    ! in : nr. iteration
     &                  debug)   ! in : =.true. -> allow debug printouts
c
c Author: Fredy W. Aquino, Northwestern University
c Date  : 03-15-12

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
      integer iter,ipm,cc,ncomp
      double precision rmax,omg(2),
     &                 omega,rmx(2)
      integer g_r(ncomp),g_Ax(ncomp)
      logical debug,converge_precond
      external precond ! preconditioner routine
c       convergence checking:
c       find the largest element of the residual either 
c       before or after the call to the preconditioner
         omg(1)=-omega
         omg(2)= omega
        if (converge_precond) then
         do ipm=1,ncomp
          call precond(g_r(ipm),omg(ipm))
         enddo ! end-loop-ipm
        endif ! end-if-conver_precond
c ----- find Absolute maximum -------- START
         rmx(1)=0.0d0
         rmx(2)=0.0d0
         do ipm=1,ncomp
          call ga_maxelt(g_r(ipm),rmx(ipm))  
         enddo ! end-loop-ipm
         rmax = max(rmx(1),rmx(2))  
c ----- find Absolute maximum -------- END
        if (.not.converge_precond) then
         do ipm=1,ncomp
          call precond(g_r(ipm),omg(ipm))
         enddo ! end-loop-ipm
        endif ! end-if-conver_precond
        do ipm=1,ncomp
         call precond(g_Ax(ipm),omg(ipm))
        enddo ! end-loop-ipm
      return
      end
c -------------- get_precond_rmax_re -------------- END
      subroutine get_precond_rmax(
     &                  rmax,    ! out: max(g_r,g_r_im)   
     &                  g_r,     ! in : real part of g_zr
     &                  g_r_im,  ! in : imag part of g_zr 
     &                  g_Ax,    ! in : real part of g_Az
     &                  g_Ax_im, ! in : imag part of g_Az 
     &                  precond, ! in : name of preconditioner routine
     &        converge_precond,  ! in : =.true. prec->max
     &                  omega,   ! in : omega
     &                  gamwidth,! in : gamwidth
     &                  ncomp,   ! in : nr. components
     &                  iter,    ! in : nr. iteration
     &                  debug)   ! in : =.true. -> allow debug printouts
c
c Author: Fredy W. Aquino, Northwestern University
c Date  : 03-15-12

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
      integer iter,ipm,cc,ncomp
      double precision rmax,omg(2),
     &                 omega,gamwidth
      integer g_r(ncomp),g_r_im(ncomp),
     &        g_Ax(ncomp),g_Ax_im(ncomp)
      logical debug,converge_precond
      external get_maxzarr,
     &         precond ! preconditioner routine
        omg(1)=-omega
        omg(2)= omega
c       convergence checking:
c       find the largest element of the residual either 
c       before or after the call to the preconditioner
        if (converge_precond) then
         do ipm=1,ncomp
          call precond(g_r(ipm),g_r_im(ipm),omg(ipm),gamwidth)
         enddo ! end-loop-ipm
        endif ! end-if-conver_precond
          call get_maxzarr(
     &               rmax,  ! ou: max(g_re,g_im)
     &               g_r,   ! in: real      part
     &               g_r_im,! in: imaginary part
     &               ncomp, ! in: nr. components
     &               iter,  ! in: iteration nr.
     &               debug) ! in: =.true. -> allow debugging printouts
        if (.not.converge_precond) then
         do ipm=1,ncomp
          call precond(g_r(ipm),g_r_im(ipm) ,omg(ipm),gamwidth)
         enddo ! end-loop-ipm
        endif ! end-if-conver_precond
        do ipm=1,ncomp
         call precond(g_Ax(ipm),g_Ax_im(ipm),omg(ipm),gamwidth)
        enddo ! end-loop-ipm
      return
      end

      subroutine get_maxzarr(
     &               rmax,  ! ou: max(g_re,g_im)
     &               g_re,  ! in: real      part
     &               g_im,  ! in: imaginary part
     &               ncomp, ! in: nr. components
     &               iter,  ! in: iteration nr.
     &               debug) ! in: =.true. -> allow debugging printouts
c Note.- ga_maxelt ->  max(abs(a(i,j)),value)
c
c Author: Fredy W. Aquino, Northwestern University
c Date  : 03-15-12

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
      logical debug
      integer i,iter,ipm,cc,n,ncomp,
     &        l_dat,k_dat
      double precision rmax
      integer g_re(ncomp),g_im(ncomp)
      external ga_maxelt
      n=2*ncomp
      if (.not.MA_Push_Get(mt_dbl,n,'hessv jfacs',l_dat,k_dat))
     &     call errquit('get_maxzarr: cannot allocate dat',
     &                  n, MA_ERR)
        cc=0
        do ipm=1,ncomp
         call ga_maxelt(g_re(ipm),dbl_mb(k_dat+cc  ))
         call ga_maxelt(g_im(ipm),dbl_mb(k_dat+cc+1))
         cc=cc+2
        enddo ! end-loop-ncomp
        rmax=-1.0d0 ! any negative number so that it pick rmax_arr(1)
        do i=1,2*ncomp
         if (dbl_mb(k_dat+i-1).gt.rmax) rmax=dbl_mb(k_dat+i-1)
        enddo !  end-loop-i
        if (debug) then
         if (ga_nodeid().eq.0) then
          write(*,1) iter,rmax,
     &               dbl_mb(k_dat  ),dbl_mb(k_dat+1),
     &               dbl_mb(k_dat+2),dbl_mb(k_dat+3)
    1     format('(iter,rmax,rmax_arr)=(',
     &           i3,',',f15.8,',[',f15.8,',',f15.8,',',
     &           f15.8,',',f15.8,'])')
         endif
        endif ! end-if-debug
      if (.not.ma_pop_stack(l_dat))
     $  call errquit('get_maxzarr: pop problem with l_dat',555,
     &       MA_ERR)
      return
      end

      subroutine get_precond_rmax_zin(
     &                  rmax,    ! out: max(g_r,g_r_im)   
     &                  g_zr1,   ! in : complex+accumulated g_zr
     &                  g_Az1,   ! in : complex+accumulated g_Az
     &                  nsub,    ! in : pointer to current (g_zr,g_Az)
     &                  precond, ! in : name of preconditioner routine
     &        converge_precond,  ! in : =.true. prec->max
     &                  omega,   ! in : omega
     &                  gamwidth,! in : gamwidth
     &                  ncomp,   ! in : nr. components
     &                  npol,    ! in : nr. polarizations
     &                  nvir,    ! in : nr. virtual  MOs
     &                  nocc,    ! in : nr. occupied MOs
     &                  n,       ! in : =nocc * nvir
     &                  nvec,    ! in : =3 (x,y,z)
     &                  iter,    ! in : nr. iteration
     &                  debug)   ! in : =.true. -> allow debug printouts
c
c Author: Fredy W. Aquino, Northwestern University
c Date  : 03-15-12

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
      integer iter,ipm,cc,ncomp,nmx,i,
     &        npol,nvir(npol),nocc(npol)
      double precision rmax,
     &                 omg(2),
     &                 gam(2),
     &                 omega,gamwidth
      integer nsub,n,nvec,
     &        l_max,k_max,
     &        g_zr1,g_Az1,
     &        g_dre,g_dim ! scratch GA used for (g_r or g_Ax)
      logical debug,converge_precond
      external conv2reim1_u,conv2complex1_u,
     &         precond ! preconditioner routine

c --> Create (g_dre,g_dim) scratch GA arrays
        if (.not. ga_create(MT_DBL,n,nvec,
     &   'get_precond_rmax_zin: g_dre',0,0,g_dre))
     $    call errquit('get_precond_rmax_zin: failed alloc g_dre', 
     &                  nvec,GA_ERR)
      if (.not. ga_create(MT_DBL,n,nvec,
     &   'get_precond_rmax_zin: g_dre',0,0,g_dim))
     $    call errquit('get_precond_rmax_zin: failed alloc g_dim', 
     &                  nvec,GA_ERR)
      nmx=2*ncomp  
      if (.not.MA_Push_Get(mt_dbl,nmx,'hessv jfacs',l_max,k_max))
     &     call errquit('get_precond_rmax_zin: cannot allocate max',
     &                  nmx, MA_ERR)
        omg(1)=-omega
        omg(2)= omega
        gam(1)=-gamwidth
        gam(2)= gamwidth
c       convergence checking:
c       find the largest element of the residual either 
c       before or after the call to the preconditioner
        if (converge_precond) then
         cc=0
         do ipm=1,ncomp
c ------- extract g_zr1 --> (g_dre,g_dim) ------- START
            call conv2reim1_u(
     &                    g_dre,  ! out : real      arr
     &                    g_dim,  ! out : imaginary arr
     &                    g_zr1,  ! in  : = complx(g_xre,g_xim)
     &                    0,      ! in  : pointer to block
     &                    npol,   ! in  : nr. polarizations
     &                    nvir,   ! in  : nr. virtual  MOs
     &                    nocc,   ! in  : nr. occupied MOs
     &                    ipm,    ! in  : =1,2 components indices
     &                    n,      ! in  : n    rows
     &                    nvec)   ! in  : nvec columns
c ------- extract g_zr1 --> (g_dre,g_dim) ------- END
c             call precond(g_dre,g_dim,omg(ipm),gamwidth)
             call precond(g_dre,g_dim,omg(ipm),gam(ipm)) ! FA-03-12-14
c            endif
c ------- Collect max values -------------------- START
            call ga_maxelt(g_dre,dbl_mb(k_max+cc  ))
            call ga_maxelt(g_dim,dbl_mb(k_max+cc+1))
            cc=cc+2
c ------- Collect max values -------------------- END
c ------- update  (g_dre,g_dim) --> g_zr1 ------- START
            call conv2complex1_u(
     &                    g_zr1,  ! out: = complx(g_xre,g_xim)
     &                    g_dre,  ! in : real      arr
     &                    g_dim,  ! in : imaginary arr
     &                    0,      ! in  : pointer to block
     &                    npol,   ! in : nr. polarizations
     &                    nvir,   ! in  : nr. virtual  MOs
     &                    nocc,   ! in  : nr. occupied MOs
     &                    ipm,    ! in  : =1,2 components indices
     &                    n,      ! in : n    rows
     &                    nvec)   ! in : nvec columns
c ------- update  (g_dre,g_dim) --> g_zr1 ------- END
         enddo ! end-loop-ipm
        endif ! end-if-conver_precond
        if (.not.converge_precond) then
         cc=0
         do ipm=1,ncomp
           call conv2reim1_u(
     &                    g_dre,  ! out : real      arr
     &                    g_dim,  ! out : imaginary arr
     &                    g_zr1,  ! in  : = complx(g_xre,g_xim)
     &                    0,      ! in  : pointer to block
     &                    npol,   ! in  : nr. polarizations
     &                    nvir,   ! in  : nr. virtual  MOs
     &                    nocc,   ! in  : nr. occupied MOs
     &                    ipm,    ! in  : =1,2 components indices
     &                    n,      ! in  : n    rows
     &                    nvec)   ! in  : nvec columns
c ------- Collect max values -------------------- START
           call ga_maxelt(g_dre,dbl_mb(k_max+cc  ))
           call ga_maxelt(g_dim,dbl_mb(k_max+cc+1))
           cc=cc+2
c ------- Collect max values -------------------- END
c           call precond(g_dre,g_dim,omg(ipm),gamwidth)
           call precond(g_dre,g_dim,omg(ipm),gam(ipm)) ! FA-03-12-14
           call conv2complex1_u(
     &                    g_zr1,! out: = complx(g_xre,g_xim)
     &                    g_dre,! in : real      arr
     &                    g_dim,! in : imaginary arr
     &                    0,    ! in  : pointer to block
     &                    npol, ! in  : nr. polarizations
     &                    nvir, ! in  : nr. virtual  MOs
     &                    nocc, ! in  : nr. occupied MOs
     &                    ipm,  ! in  : =1,2 components indices
     &                    n,    ! in : n    rows
     &                    nvec) ! in : nvec columns
         enddo ! end-loop-ipm
        endif ! end-if-conver_precond
c ----- obtain rmax ------- START
        rmax=-1.0d0 ! any negative number so that it pick rmax_arr(1)
        do i=1,2*ncomp
         if (dbl_mb(k_max+i-1).gt.rmax) rmax=dbl_mb(k_max+i-1)
        enddo !  end-loop-i
        if (debug) then
         if (ga_nodeid().eq.0) then
          write(*,1) iter,rmax,
     &               dbl_mb(k_max  ),dbl_mb(k_max+1),
     &               dbl_mb(k_max+2),dbl_mb(k_max+3)
    1     format('(iter,rmax,rmax_arr)=(',
     &           i3,',',f15.8,',[',f15.8,',',f15.8,',',
     &           f15.8,',',f15.8,'])')
         endif
        endif ! end-if-debug
c ----- obtain rmax ------- END
        do ipm=1,ncomp
          call conv2reim1_u(
     &                    g_dre,  ! out : real      arr
     &                    g_dim,  ! out : imaginary arr
     &                    g_Az1,  ! in  : = complx(g_xre,g_xim)
     &                    nsub,   ! in  : pointer to block
     &                    npol,   ! in  : nr. polarizations
     &                    nvir,   ! in  : nr. virtual  MOs
     &                    nocc,   ! in  : nr. occupied MOs
     &                    ipm,    ! in  : =1,2 components indices
     &                    n,      ! in  : n    rows
     &                    nvec)   ! in  : nvec columns
c             call precond(g_dre,g_dim,omg(ipm),gamwidth)
             call precond(g_dre,g_dim,omg(ipm),gam(ipm)) ! FA-03-12-14
          call conv2complex1_u(
     &                    g_Az1,! out: = complx(g_xre,g_xim)
     &                    g_dre,! in : real      arr
     &                    g_dim,! in : imaginary arr
     &                    nsub, ! in : pointer to block
     &                    npol, ! in : nr. polarizations
     &                    nvir, ! in : nr. virtual  MOs
     &                    nocc, ! in : nr. occupied MOs
     &                    ipm,  ! in : =1,2 components indices
     &                    n,    ! in : n    rows
     &                    nvec) ! in : nvec columns
        enddo ! end-loop-ipm
        if (.not. ga_destroy(g_dre)) call errquit
     &     ('get_precond_rmax_zin: g_dre',0, GA_ERR)
        if (.not. ga_destroy(g_dim)) call errquit
     &     ('get_precond_rmax_zin: g_dim',0, GA_ERR)
      if (.not.ma_pop_stack(l_max))
     $  call errquit('get_precond_rmax_zin: pop problem with l_max',
     &               555,MA_ERR)
      return
      end
c =============== FA:complex solver ==================== END
