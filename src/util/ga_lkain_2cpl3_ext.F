      subroutine ga_lkain_2cpl3_ext(rtdb, 
     &                          g_x, g_b, 
     &                          g_x_im, g_b_im,
     &                          product, precond, 
     $                          tol, mmaxsub, maxiter, 
     &                          odiff, oprint, omega, limag,
     &                          lifetime, gamwidth, ncomp)

c     $Id: ga_lkain_2cpl3.F 19707 2010-10-29 17:59:36Z d3y133 $
c  Written by J. Autschbach, SUNY Buffalo
c  Improvements made
c          by F. Aquino,     Northwestern University 
c          03-15-12
c  Note.- Modifying/Improving ga_lkain_2cpl3()
c --> Experimental (not published yet)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
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
      external solve_xlineq,updating_Az1_z1_zr1,
     &         get_precond_rmax_re
c
c     =================================================================

      debug = (.false. .and. ga_nodeid().eq.0) ! for code development

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
     &                            acc)) acc = 0.01d0*tol
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

c           if (ga_nodeid().eq.0) then
c             write(*,10) iter
c 10          format('------ g_x-BEF-product(',i3,')---START')
c           endif
c            do ipm=1,ncomp
c             call ga_print(g_x(ipm))
c             call ga_print(g_Ax(ipm))
c             if (lifetime) then
c              call ga_print(g_x_im(ipm))
c              call ga_print(g_Ax_im(ipm))
c             endif
c            enddo ! end-loop-ipm
c           if (ga_nodeid().eq.0) then
c             write(*,11) iter
c 11          format('------ g_x-BEF-product(',i3,')---END')
c           endif

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
    
        if (oprint .and. ga_nodeid().eq.0) then
          write(6,3) iter, nsub+nvec, rmax, util_wallsec()
          call util_flush(6)
    3     format(' ', i5, i7, 3x,1p,d9.2,0p,f10.1,5x,i3)
        end if

c       stop iterations if residual is smaller than criterion
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
          write (luout,*) '  I will proceed, but check your results '//
     &                       'carefully!!!'
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
c --> Experimental (not published yet)

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
c --> Experimental (not published yet)

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
c --> Experimental (not published yet)

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
c --> Experimental (not published yet)

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
c ============================================================
c =============== FA:complex solver ==================== START
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +++++++++++++++ ga_lkain_2cpl3_damp_cmplx_debug +++++++++++ START
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c Note.- This routine is kept only for debugging purpose 
c        It works fine but it stops in 5th iteration.
c        At the end of 5th iteration it will call
c        rohf_hessv_precond3_damp_anti() which does the inverse
c        transformation done in rohf_hessv_precond3_damp()
c        The xxx_anti() routine works fine, I tested it doing
c        rohf_hessv_precond3_damp() --> rohf_hessv_precond3_damp_anti()
c        to see if it undo properly the first transformation and it
c        does.  This routine was written to see if I can see some
c        relationship between g_x(1) and g_x(2) after undoing energy
c        transformation.
c        The routine is set to use the complex linear solver.
c        If we want to see back the real linear solver to work, 
c        do the following:
c        a) Comment  : goto 778 ! avoid-real-linear-solver-1
c        b) Uncomment: c       goto 9111 ! skip so that it does not use complex solution
c        c) Comment  : goto 777 ! avoid-real-linear-solver-2
c        For the case of water (h2o) with (xc,bs)=(b3lyp,6-31G*)
c        I get the following convergence behaviour:
c         A) Using real linear solver:
c           (iter,rmax,tol)=(   1,     0.29449271,     0.00000100)
c           (iter,rmax,tol)=(   2,     0.21804001,     0.00000100)
c           (iter,rmax,tol)=(   3,     0.07201125,     0.00000100)
c           (iter,rmax,tol)=(   4,     0.00492173,     0.00000100)
c           (iter,rmax,tol)=(   5,     0.00157834,     0.00000100)
c           (iter,rmax,tol)=(   6,     0.00076566,     0.00000100)
c           (iter,rmax,tol)=(   7,     0.00000797,     0.00000100)
c         B) Using complex linear solver:
c           (iter,rmax,tol)=(   1,     0.29449271,     0.00000100)
c           (iter,rmax,tol)=(   2,     0.21802670,     0.00000100)
c           (iter,rmax,tol)=(   3,     0.07199298,     0.00000100)
c           (iter,rmax,tol)=(   4,     0.00491850,     0.00000100)
c           (iter,rmax,tol)=(   5,     0.00044391,     0.00000100)
c           (iter,rmax,tol)=(   6,     0.00001302,     0.00000100)
c           (iter,rmax,tol)=(   7,     0.00000090,     0.00000100)
c        From those results faster convergence is observed if using
c        complex linear solver.
c  FA-04-13-12
      subroutine ga_lkain_2cpl3_damp_cmplx_debug(
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
     &                   ncomp)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "stdio.fh"
#include "rtdb.fh"
c
      integer rtdb              ! [input] database handle
      integer g_x(2)            ! [input/output] Initial guess/solution Re
      integer g_x_im(2)         ! [input/output] Initial guess/solution Im
      integer g_b(2)            ! [input] Right-hand side vectors Re
      integer g_b_im(2)         ! [input] Right-hand side vectors Im
      double precision omega    ! [input] coupling parameter
      logical limag             ! [input] imaginary perturbation?
      logical lifetime          ! [input] consider damping or not?
      double precision gamwidth ! [input] damping parameter
      integer ncomp             ! [input] no. of components to treat
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

      integer iter, n, n1,n2, nvec, nsub, isub, type, maxsub, ipm,
     &   ntmp1, ntmp2, n3, n4

c ... jochen: for convenience, now most arrays have two components.
c     that might be changed later if memory becomes an issue
      integer g_a, g_bb,g_c, 
     &        g_xold(2), g_Axold(2)
      integer g_r(2) ,g_r_im(2), 
     &        g_Ax(2),g_Ax_im(2),
     &        g_y, g_Ay, g_r2
      double precision rmax, rmax1, rmax2, acc
      logical converged
      logical odebug, debug, converge_precond
      double complex val_cmplx
      logical debug1
      integer p1,p2,m1,m2,stat_solve
      
      integer g_z(ncomp),g_Az(ncomp),g_zr(ncomp),g_zb(ncomp),
     &        g_z1,g_Az1,g_zr1
      double precision omg(ncomp)
      external copy_rtor2,copy_r2tor,copy_AxxtoAyy,
     &         update_g_x1,conv2complex,
     &         updating_Az1_z1_zr1,
     &         solve_zlineq,get_precond_rmax,
     &         conv2reim,
     &         rohf_hessv_precond3_damp_anti

c     later we combine the two components to vecors of double
c     length and combine again Re and Im, i.e. 
c     the dimension is up to 4*n

      call ga_inquire(g_x(1), type, n, nvec) ! get (n,nvec)

      if (ga_nodeid().eq.0) then
       write(*,14) n,ncomp,maxiter,mmaxsub,lifetime
 14    format('(n,ncomp,maxiter,mmaxsub,lifetime)=(',
     &          i3,',',i3,',',i5,',',i5,',',L1,')')
      endif
c ========= to be removed ============ start
      n2 = n
      if (ncomp.gt.1) n2 = 2 * n   
      n3 = n
      if (lifetime .or. ncomp.gt.1) n3 = 2 * n
      if (lifetime .and. ncomp.gt.1) n3 = 3 * n
      n4 = n
      if (ncomp.gt.1 .or. lifetime) n4 = 2* n
      if (lifetime .and. ncomp.gt.1) n4 = 4 * n
      if (debug) write (6,*) 'n1n2n3n4',n,n2,n3,n4
c ========= to be removed ============ end
      maxsub = mmaxsub          ! So don't modify input scalar arg
      if (maxsub .lt. 3*nvec) maxsub = 3*nvec
      maxsub = (maxsub/nvec)*nvec 

      if (ga_nodeid().eq.0) then
       write(*,1023) n4,mmaxsub,maxsub,n
 1023  format('(n4,mmaxsub,maxsub,n)=(',
     &        i15,',',i5,',',i5,',',i15,')')
      endif

      debug1 = .true. ! no printouts
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
         if (.not. ga_create(MT_DCPL,n,nvec, 'lkain_2cpl: zb',
     $     0, 0, g_zb(ipm)))
     $     call errquit('lkain: failed allocating zb', nvec,
     &     GA_ERR)
           call ga_zero(g_zb(ipm))
        enddo ! end-loop-ipm
       call conv2complex(g_zb,   ! out: = complx(g_x,g_x_im)
     &                   g_b,    ! in : real      arr
     &                   g_b_im, ! in : imaginary arr
     &                   n,      ! in : n    rows
     &                   nvec,   ! in : nvec columns
     &                   ncomp)  ! in : nr. components

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
       do ipm=1,ncomp
       if (ga_nodeid().eq.0) 
     &  write(*,*) '--------g_zb(',ipm,') -------- START'
       call ga_print(g_zb(ipm))
       if (ga_nodeid().eq.0) 
     &  write(*,*) '--------g_zb(',ipm,') -------- END'
       enddo ! end-loop-ipm
       endif ! end-if-debug1
      n1=ncomp*n
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
        write(6,1) n4, nvec, maxsub, tol, util_wallsec()
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
      
c     allocate g_y, g_Ay, and g_r2 with dimension n4 to hold
c     all necessary components simultaneously
      if (.not. ga_create(MT_DBL, n4, maxsub, 'lkain_2cpl: Y', 
     $   0, 0, g_y))
     $   call errquit('lkain: failed allocating subspace', maxsub,
     &   GA_ERR)
      if (.not. ga_create(MT_DBL, n4, maxsub, 'lkain_2cpl: Ay', 
     $   0, 0, g_Ay))
     $   call errquit('lkain: failed allocating subspace2', maxsub,
     &   GA_ERR)
      if (.not. ga_create(MT_DBL, n4, nvec, 'lkain_2cpl: r2',
     $   0, 0, g_r2))
     $   call errquit('lkain_2cpl: failed allocating r2', nvec,
     &   GA_ERR)
      
      call ga_zero(g_y)
      call ga_zero(g_Ay)
      call ga_zero(g_r2)
      call ga_sync()
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
        if (ga_nodeid().eq.0)
     &   write(*,*) 'FA-chk-odiff=',odiff

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
 112          format('------ prod-g_x-1(',i3,',',i3,')------ START')
             endif
             call ga_print(g_x(ipm))
             if (ga_nodeid().eq.0) then
              write(*,113) iter,ipm
 113          format('------ prod-g_x-1(',i3,',',i3,')------ END')
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
 116          format('------ prod-g_x-2(',i3,',',i3,')------ START')
             endif
             call ga_print(g_x(ipm))
             if (ga_nodeid().eq.0) then
              write(*,117) iter,ipm
 117          format('------ prod-g_x-2(',i3,',',i3,')------ END')
             endif
             if (ga_nodeid().eq.0) then
              write(*,118) iter,ipm
 118          format('------ prod-g_Ax-2(',i3,',',i3,')------ START')
             endif
             call ga_print(g_Ax(ipm))
             if (ga_nodeid().eq.0) then
              write(*,119) iter,ipm
 119          format('------ prod-g_Ax-2(',i3,',',i3,')------ END')
             endif
             enddo ! end-loop-ipm
           endif ! end-if-debug  

        if (debug) write (6,*)
     &     'returning product from ga_lkain_2cpl_damp'

c       g_r is zeroed below so we should make sure to do the same
c       with g_r2 here
        call ga_zero(g_r2)
        call ga_sync()
        
        do ipm = 1,ncomp
          
          if (odiff) then
c           jochen: odiff stuff presently ignored
            call ga_add(1.0d0,g_Ax(ipm), 
     &                  1.0d0,g_Axold(ipm), 
     &                        g_Ax(ipm))
            call ga_add(1.0d0,g_x(ipm),  
     &                  1.0d0,g_xold(ipm),
     &                        g_x(ipm))
            call ga_sync()
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

          call ga_sync()
c FA: Step 1:
          call ga_add( 1.0d0,g_b(ipm),
     &                -1.0d0,g_Ax(ipm), 
     &                       g_r(ipm))    ! The residual, Real part
          call ga_add( 1.0d0,g_b_im(ipm),
     &                -1.0d0,g_Ax_im(ipm), 
     &                       g_r_im(ipm)) ! The residual, Im part

          if (debug) then
             if (ga_nodeid().eq.0) then
              write(*,120) iter,ipm
 120          format('------ prod-g_b(',i3,',',i3,')------ START')
             endif
             call ga_print(g_b(ipm))
             if (ga_nodeid().eq.0) then
              write(*,121) iter,ipm
 121          format('------ prod-g_b(',i3,',',i3,')------ END')
             endif
             if (ga_nodeid().eq.0) then
              write(*,122) iter,ipm
 122          format('------ prod-g_r(',i3,',',i3,')------ START')
             endif
             call ga_print(g_r(ipm))
             if (ga_nodeid().eq.0) then
              write(*,123) iter,ipm
 123          format('------ prod-g_r(',i3,',',i3,')------ END')
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
c       stop iterations if residual is smaller than criterion
        if (rmax .lt. tol) then
          converged = .true.
          goto 100
        end if
               
c       Copy the vectors to the subspace work area
c   
c ---- FA-copy ((Ax,Ax_im),(x,x_im)) -> (Ay,y) ------ START
c FA: Step 2:
        if (debug1) then
         do ipm=1,ncomp
          if (ga_nodeid().eq.0) then
           write(*,3000) ipm,iter
 3000      format('---------g_xre-0(',i3,',',i3,')-----START')
          endif
          call ga_print(g_x(ipm))
          if (ga_nodeid().eq.0) then
           write(*,2999) ipm,iter
 2999      format('---------g_xre-0(',i3,',',i3,')-----END')
          endif
          if (ga_nodeid().eq.0) then
           write(*,2996) ipm,iter
 2996      format('---------g_xim-0(',i3,',',i3,')-----START')
          endif
          call ga_print(g_x_im(ipm))
          if (ga_nodeid().eq.0) then
           write(*,2997) ipm,iter
 2997      format('---------g_xim-0(',i3,',',i3,')-----END')
          endif
          if (ga_nodeid().eq.0) then
           write(*,2995) ipm,iter
 2995      format('---------g_Axre-0(',i3,',',i3,')-----START') 
          endif
          call ga_print(g_Ax(ipm))
          if (ga_nodeid().eq.0) then
           write(*,2994) ipm,iter
 2994      format('---------g_Axre-0(',i3,',',i3,')-----END')
          endif
          if (ga_nodeid().eq.0) then
           write(*,2993) ipm,iter
 2993      format('---------g_Axim-0(',i3,',',i3,')-----START') 
          endif
          call ga_print(g_Ax_im(ipm))
          if (ga_nodeid().eq.0) then
           write(*,2992) ipm,iter
 2992      format('---------g_Axim-0(',i3,',',i3,')-----END')
          endif
          if (ga_nodeid().eq.0) then
           write(*,2805) ipm,iter
 2805      format('---------g_r-re-0(',i3,',',i3,')-----START') 
          endif
          call ga_print(g_r(ipm))
          if (ga_nodeid().eq.0) then
           write(*,2806) ipm,iter
 2806      format('---------g_r-re-0(',i3,',',i3,')-----END')
          endif
          if (ga_nodeid().eq.0) then
           write(*,2807) ipm,iter
 2807      format('---------g_r-im-0(',i3,',',i3,')-----START') 
          endif
          call ga_print(g_r_im(ipm))
          if (ga_nodeid().eq.0) then
           write(*,2808) ipm,iter
 2808      format('---------g_r-im-0(',i3,',',i3,')-----END')
          endif
         enddo ! end-loop-ipm
        endif ! end-if-debug1  
c +++++++++++++ test: build a complex x +++++++++++++ START
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
c +++++++++++++ test: build a complex x +++++++++++++ END

c ================= to be removed =========== START
        goto 778 ! avoid-real-linear-solver-1
        call copy_AxxtoAyy(g_Ax,g_Ax_im,
     &                     g_x,g_x_im,
     &                     g_Ay,g_y,
     &                     nvec,
     &                     ncomp,
     &                     nsub,
     &                     n,
     &                     lifetime)
c ---- FA-copy ((Ax,Ax_im),(x,x_im)) -> (Ay,y) ------ END
c ----- FA--- (g_r,g_r_im) -> g_r2 ---------- START
c FA: Step 3:
         call copy_rtor2(g_r2,
     &                   g_r,
     &                   g_r_im,
     &                   ncomp,
     &                   nvec,
     &                   n,
     &                   lifetime)
c ----- FA--- (g_r,g_r_im) -> g_r2 ---------- END 
        if (debug1) then
        if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_y-0(',iter,')-----START' 
        call ga_print(g_y)
        if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_y-0(',iter,')-----END'          
        if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_Ay-0(',iter,')-----START' 
        call ga_print(g_Ay)
        if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_Ay-0(',iter,')-----END'          
        if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_r2-0(',iter,')-----START' 
        call ga_print(g_r2)
        if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_r2-0(',iter,')-----END'    
        endif ! end-if-debug1
 778    continue
c ================= to be removed =========== END


        nsub = nsub + nvec
c ============ doing complex-linear solving ======= START
        if (ga_nodeid().eq.0)
     &   write(*,*) 'BEFORE solve_zlineq ...'
c +++++++++++++ keep this as reminder +++++++++ START
c        if (ga_nodeid().eq.0)
c     &   write(*,*) 'FA-test ga_lu_solve ...'
c        call ga_lu_solve('n',g_a,g_bb)

c        stat_solve=ga_solve(g_a,g_bb)
c        if (ga_nodeid().eq.0) then
c         write(*,14) stat_solve
c 14      format('FA-test ga_lu_solve(Choleski->LU) stat=',i3)
c        endif
c        call ga_copy(g_bb,g_c)
c +++++++++++++ keep this as reminder +++++++++ END
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

c       goto 9111 ! skip so that it does not use complex solution
       if (ga_nodeid().eq.0)
     &  write(*,*) 'BEF conv2reim-x'
       call conv2reim(g_x,   ! out : real      arr
     &                g_x_im,! out : imaginary arr
     &                g_z,   ! in  : = complx(g_xre,g_xim)
     &                n,     ! in  : n    rows
     &                nvec,  ! in  : nvec columns
     &                ncomp) ! in  : nr. components   
       if (ga_nodeid().eq.0)
     &  write(*,*) 'AFT conv2reim-x'
9111    continue

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
        if (ga_nodeid().eq.0)
     &   write(*,*) 'AFTER solve_zlineq ...'
c ============ doing complex-linear solving ======= END   

c =============================================================
c ===================== to be removed =================== START
        goto 777 ! avoid-real-linear-solver-2
c       Form and solve the subspace equations using SVD in order
c       to manage near linear dependence in the subspace.
c       
        if (.not. ga_create(MT_DBL, nsub, nsub,
     &     'lkain_2cpl3_damp: A', 0, 0, g_a))
     $     call errquit('lkain: allocating g_a?', nsub, GA_ERR)
        if (.not. ga_create(MT_DBL, nsub, nvec,
     &     'lkain_2cpl3_damp: B', 0, 0,g_bb))
     $     call errquit('lkain: allocating g_bb?', nsub, GA_ERR)
        if (.not. ga_create(MT_DBL, nsub, nvec,
     &     'lkain_2cpl3_damp: C', 0, 0, g_c))
     $     call errquit('lkain: allocating g_c?', nsub, GA_ERR)
        call ga_zero(g_a)
        call ga_zero(g_bb)
        call ga_zero(g_c)
        call ga_dgemm('t','n',nsub,nsub,n4,1.0d0,
     &                g_y,g_Ay,0.0d0,g_a)
        call ga_dgemm('t','n',nsub,nvec,n4,1.0d0,
     &                g_y,g_r2,0.0d0,g_bb)

        if (odebug) call ga_print(g_a)
        if (odebug) call ga_print(g_c)
c       
c       The threshold used here should reflect the accuracy in the
c       products.  If very accurate products are used,
c       then there is big
c       advantage for small cases (maxsub close to n) in using a very
c       small threshold in the SVD solve (e.g., 1e-14), but for more
c       realistic examples (maxsub << n) there is only a little
c       advantage and in the precence of real noise in the products
c       screening with a realistic threshold is important.
        if (debug1) then
        if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_a(',iter,')-----START' 
        call ga_print(g_a)
        if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_a(',iter,')-----END'
        if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_bb(',iter,')-----START' 
        call ga_print(g_bb)
        if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_bb(',iter,')-----END'
        endif ! end-if-debug1
c ++++++++++++++++++++++++++++++++++++++++++++++++++++    
c ++++++++++++++ real linear solver ++++++++++++++++++ 
         call ga_svd_solve_seq(g_a,g_bb,g_c,1d-14)

c        if (ga_nodeid().eq.0)
c     &   write(*,*) 'FA-test ga_lu_solve ...'
c        call ga_lu_solve('n',g_a,g_bb)

c        stat_solve=ga_solve(g_a,g_bb)
c        if (ga_nodeid().eq.0) then
c         write(*,14) stat_solve
c 14      format('FA-test ga_lu_solve(Choleski->LU) stat=',i3)
c        endif
c        call ga_copy(g_bb,g_c)
c ++++++++++++++++++++++++++++++++++++++++++++++++++++
c ++++++++++++++++++++++++++++++++++++++++++++++++++++
        if (debug1) then
        if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_c-old(',iter,')-----START' 
        call ga_print(g_c)
        if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_c-old(',iter,')-----END'
        endif ! end-if-debug1
        
        if (odebug) call ga_print(g_c)
c       
c       Form and add the correction, in parts, onto the solution      
c FA: Step 5:
        if (debug1) then
        if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_r2-BEF(',iter,')-----START' 
        call ga_print(g_r2)
        if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_r2-BEF(',iter,')-----END'
        endif ! end-if-debug1

        call ga_dgemm('n','n',n4,nvec,nsub,-1.0d0,
     &                g_Ay,g_c,1.0d0,g_r2)

        if (debug1) then
        if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_r2-AFT(',iter,')-----START' 
        call ga_print(g_r2)
        if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_r2-AFT(',iter,')-----END'
        if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_x-BEF-1(',iter,')-----START' 
        call ga_print(g_x)
        if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_x-BEF-1(',iter,')-----END'
        endif ! end-if-debug1
c       
c       copy components of g_r2 into g_r before adding g_r to  g_x
       call update_g_x1(g_r2,    ! in   :
     &                 g_x,     ! in/ou: updated
     &                 g_x_im,  ! in/ou: updated
     &                 ncomp,   ! in   : nr. components
     &                 nvec,    ! in   : (x,y,z) = 3
     &                 n,       ! in   : size of block (e.g. RE1,RE2,IM1 or IM2)
     &                 lifetime)! in   : = .true. if complex
        if (debug1) then
        do ipm=1,ncomp
          if (ga_nodeid().eq.0) then
           write(*,332) ipm,iter
 332       format('---------g_x-1(',i3,',',i3,')-----START')
          endif
          call ga_print(g_x(ipm))
          if (ga_nodeid().eq.0) then
           write(*,333) ipm,iter
 333       format('---------g_x-1(',i3,',',i3,')-----END')
          endif
        enddo ! end-loop-ipm
        do ipm=1,ncomp
          if (ga_nodeid().eq.0) then
           write(*,334) ipm,iter
 334       format('---------g_x-im-1(',i3,',',i3,')-----START')
          endif
          call ga_print(g_x_im(ipm))
          if (ga_nodeid().eq.0) then
           write(*,335) ipm,iter
 335       format('---------g_x-im-1(',i3,',',i3,')-----END')
          endif
        enddo ! end-loop-ipm
        endif ! end-if-debug1
c FA: Step 8:
        call ga_dgemm('n','n',n4,nvec,nsub,1.0d0,
     &                g_y,g_c,0.0d0,g_r2)
        if (debug1) then
        if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_y c(',iter,')-----START' 
        call ga_print(g_r2)
        if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_y c(',iter,')-----END'
        endif ! end-if-debug1
c       copy components of g_r2 into g_r before adding g_r to  g_x
       call update_g_x1(g_r2,    ! in   :
     &                 g_x,     ! in/ou: updated
     &                 g_x_im,  ! in/ou: updated
     &                 ncomp,   ! in   : nr. components
     &                 nvec,    ! in   : (x,y,z) = 3
     &                 n,       ! in   : size of block (e.g. RE1,RE2,IM1 or IM2)
     &                 lifetime)! in   : = .true. if complex

        if (debug1) then
        do ipm=1,ncomp
          if (ga_nodeid().eq.0) then
           write(*,336) ipm,iter
 336       format('---------g_x-2(',i3,',',i3,')-----START')
          endif
          call ga_print(g_x(ipm))
          if (ga_nodeid().eq.0) then
           write(*,337) ipm,iter
 337       format('---------g_x-2(',i3,',',i3,')-----END')
          endif
        enddo ! end-loop-ipm
        do ipm=1,ncomp
          if (ga_nodeid().eq.0) then
           write(*,338) ipm,iter
 338       format('---------g_x-im-2(',i3,',',i3,')-----START')
          endif
          call ga_print(g_x_im(ipm))
          if (ga_nodeid().eq.0) then
           write(*,339) ipm,iter
 339       format('---------g_x-im-2(',i3,',',i3,')-----END')
          endif
        enddo ! end-loop-ipm
        endif ! end-if-debug1
       
        if (.not. ga_destroy(g_a)) call errquit
     &     ('lkain_2cpl: a',0, GA_ERR)
        if (.not. ga_destroy(g_bb)) call errquit
     &     ('lkain_2cpl: b',0, GA_ERR)
        if (.not. ga_destroy(g_c)) call errquit
     &     ('lkain_2cpl: c',0, GA_ERR)
c ===================== to be removed =================== END
c =============================================================
 777    continue

        if (ga_nodeid().eq.0) then
         write(*,17) iter,nsub,maxsub,nvec,n4
 17      format('(iter,nsub,maxsub,nvec,n4)=(',
     &          i3,',',i4,',',i4,',',i3,',',i4,')')
        endif


        if (debug1) then
         if (iter.eq.2) then
          write(*,*) 'STOP-test-linsolver'
          stop
         endif
        endif ! end-if-debug1

        goto 1998

        if (debug1) then
         if (iter.eq.5) then

          if (ga_nodeid().eq.0)
     &     write(*,*) 'FA-BEF rohf_hessv_precond3_damp_anti'
          omg(1)=-omega
          omg(2)= omega
          do ipm=1,ncomp
          call rohf_hessv_precond3_damp_anti(
     &                        g_x(ipm), 
     &                        g_x_im(ipm), 
     &                        omg(ipm), 
     &                        gamwidth)
          enddo ! end-loop-ipm
          if (ga_nodeid().eq.0)
     &     write(*,*) 'FA-AFT rohf_hessv_precond3_damp_anti'
        if (debug1) then
        do ipm=1,ncomp
          if (ga_nodeid().eq.0) then
           write(*,436) ipm,iter
 436       format('---------g_x-anti(',i3,',',i3,')-----START')
          endif
          call ga_print(g_x(ipm))
          if (ga_nodeid().eq.0) then
           write(*,437) ipm,iter
 437       format('---------g_x-anti(',i3,',',i3,')-----END')
          endif
        enddo ! end-loop-ipm
        do ipm=1,ncomp
          if (ga_nodeid().eq.0) then
           write(*,438) ipm,iter
 438       format('----g_x-im-anti(',i3,',',i3,')-----START')
          endif
          call ga_print(g_x_im(ipm))
          if (ga_nodeid().eq.0) then
           write(*,439) ipm,iter
 439       format('----g_x-im-anti(',i3,',',i3,')-----END')
          endif
        enddo ! end-loop-ipm
        endif ! end-if-debug1
          write(*,*) 'STOP-test-linsolver'
          stop
         endif
        endif ! end-if-debug1 

1998    continue

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

         do isub = nvec+1, maxsub, nvec
          call ga_copy_patch('n',g_Ay,1,n4,isub,isub+nvec-1, 
     $                           g_Ay,1,n4,isub-nvec,isub-1)
          call ga_copy_patch('n',g_y ,1,n4,isub,isub+nvec-1, 
     $                           g_y ,1,n4,isub-nvec,isub-1)
         enddo ! end-loop-isub
c ====== FA: left-shifting patch ==== END
         nsub = nsub - nvec
        end if                  ! (nsub .eq. maxsub) 
      enddo                     ! iter = 1,maxiter
  100 continue                  ! jump here if converged
c     
c     
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

      if (.not. ga_destroy(g_Ay)) call errquit
     &   ('lkain_2cpl: destroy Ay',3, GA_ERR)
      if (.not. ga_destroy(g_y)) call errquit
     &   ('lkain_2cpl: destroy r',4, GA_ERR)
      if (.not. ga_destroy(g_r2)) call errquit
     &   ('lkain_2cpl: destroy r2',6, GA_ERR)
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

      if (debug) then
        write (6,*) 'ga_lkain_2cpl3_damp: solution vectors:'
        call ga_print(g_x(1))
        if (ncomp.gt.1) call ga_print(g_x(2))
        if (lifetime) then
          call ga_print(g_x_im(1))
          if (ncomp.gt.1) call ga_print(g_x_im(2))
        endif
      endif
      
      
c ... jochen: disable this error exit during debuging phase
c     but print a warning instead

c     error exit if this hasn't converged :-(
c     if (.not. converged) call errquit('lkain_2cpl: not converged',0,
c     &    CALC_ERR)

      if (.not. converged) then
        if (ga_nodeid().eq.0) then
          write (luout,*) 'WARNING: CPKS procedure is NOT converged'
          write (luout,*) '  I will proceed, but check your results '//
     &                       'carefully!!!'
        endif
      endif
c     
      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +++++++++++++++ ga_lkain_2cpl3_damp_cmplx_debug +++++++++++ END
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c ++++++++++++++++++++clean routine++++++++++++++++++++ START
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
c --> Experimental (not published yet)

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
      double precision omg(ncomp)
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
 14    format('(n,ncomp,maxiter,mmaxsub,lifetime)=(',
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
 112          format('------ prod-g_x-1(',i3,',',i3,')------ START')
             endif
             call ga_print(g_x(ipm))
             if (ga_nodeid().eq.0) then
              write(*,113) iter,ipm
 113          format('------ prod-g_x-1(',i3,',',i3,')------ END')
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
 116          format('------ prod-g_x-2(',i3,',',i3,')------ START')
             endif
             call ga_print(g_x(ipm))
             if (ga_nodeid().eq.0) then
              write(*,117) iter,ipm
 117          format('------ prod-g_x-2(',i3,',',i3,')------ END')
             endif
             if (ga_nodeid().eq.0) then
              write(*,118) iter,ipm
 118          format('------ prod-g_Ax-2(',i3,',',i3,')------ START')
             endif
             call ga_print(g_Ax(ipm))
             if (ga_nodeid().eq.0) then
              write(*,119) iter,ipm
 119          format('------ prod-g_Ax-2(',i3,',',i3,')------ END')
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
 120          format('------ prod-g_b(',i3,',',i3,')------ START')
             endif
             call ga_print(g_b(ipm))
             if (ga_nodeid().eq.0) then
              write(*,121) iter,ipm
 121          format('------ prod-g_b(',i3,',',i3,')------ END')
             endif
             if (ga_nodeid().eq.0) then
              write(*,122) iter,ipm
 122          format('------ prod-g_r(',i3,',',i3,')------ START')
             endif
             call ga_print(g_r(ipm))
             if (ga_nodeid().eq.0) then
              write(*,123) iter,ipm
 123          format('------ prod-g_r(',i3,',',i3,')------ END')
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
          write (luout,*) '  I will proceed, but check your results '//
     &                       'carefully!!!'
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
c --> Experimental (not published yet)

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
      external solve_zlineq1,conv2reim_rhs,
     &         get_precond_rmax_zin,
     &         dft_CPHF2_read,
     &         dft_CPHF2_write,
     &         solve_zlineq1_sep_spins,
     &         dft_CPHF2_read2fix,
     &         solve_zlineq_PopleAlg,
     &         ran1
      logical status,flag2readfile
      integer index4cphf,checkorth,cphf3write
      character*255 aorespfilename
      character*(*) lbl_cphfaoresp
      character*255 lbl_cphfaoresp1
      integer iimoderaman,iiistepraman ! FA-06-18-12-raman defs
      logical ramanspc                 ! FA-06-18-12-raman defs
      parameter(lbl_cphfaoresp='aoresp_fiao_f')
      one_cmplx =dcmplx( 1.0d0,0.0d0)
      mone_cmplx=dcmplx(-1.0d0,0.0d0)
      zero_cmplx=dcmplx( 0.0d0,0.0d0)
      if (.not. rtdb_get(rtdb, 'cphf:cphf3write',    mt_int, 1,
     &                   cphf3write))
     &  cphf3write = 0 ! assigns 0 if unsuccessfull read from rtdb
c      if (ga_nodeid().eq.0)
c     & write(*,*) 'In ga_lkain_2cpl3_damp_cmplx_redmem: cphf3write=',
c     &           cphf3write
      if (.not. rtdb_get(rtdb, 'cphf:checkorth',    mt_int, 1,
     &                   checkorth))
     &  checkorth = 0 ! assigns 0 if unsuccessfull read from rtdb
c      if (ga_nodeid().eq.0)
c     & write(*,*) 'In ga_lkain_2cpl3_damp_cmplx_redmem: checkorth=',
c     &           checkorth

      call ga_inquire(g_z(1),type,n,nvec) ! get (n,nvec) n=sum(nocc*nvirt(i) i=1,npol)

c      if (ga_nodeid().eq.0) then
c       write(*,14) n,ncomp,maxiter,mmaxsub,lifetime
c 14    format('(n,ncomp,maxiter,mmaxsub,lifetime)=(',
c     &          i3,',',i3,',',i5,',',i5,',',L1,')')
c      endif
      n1=ncomp*n
      maxsub = mmaxsub          ! So don't modify input scalar arg
      if (maxsub .lt. 3*nvec) maxsub = 3*nvec
      maxsub = (maxsub/nvec)*nvec 

c      if (ga_nodeid().eq.0) then
c       write(*,1023) n1,mmaxsub,maxsub,n
c 1023  format('(n1,mmaxsub,maxsub,n)=(',
c     &        i15,',',i5,',',i5,',',i15,')')
c      endif

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
c      if (.not. rtdb_get(rtdb,'cphf3-aores:guess1', 
c     &          mt_int, 1,index4cphf)) call
c     $     errquit('ga_lkain_2cpl3_redmem: failed to read index4cphf ', 
c     &             0, RTDB_ERR)

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
c       write(*,112) iimoderaman,iiistepraman,lbl_cphfaoresp1
c 112   format('(iimoderaman,iiistepraman,lbl_cphfaoresp1)=(',
c     &         i4,',',i4,',',a,')')                  
      else
       write(lbl_cphfaoresp1,'(a13,i1)') lbl_cphfaoresp,index4cphf    
      endif
       call util_file_name(lbl_cphfaoresp1,
     &                     .false.,.false.,aorespfilename)
c      if (ga_nodeid().eq.0) then
c        write(*,*) 'ga_lkain_2cpl3_redmem: lbl_cphfaoresp1=',
c     &            lbl_cphfaoresp1,' aorespfilename=',
c     &            aorespfilename
c      endif  

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
       if (ga_nodeid().eq.0)
     &  write(*,1999) aorespfilename(1:inp_strlen(aorespfilename))
 1999   format('File ',a,
     &          ' does not exist, proceed to generate (z1,Az1)')
c ------ g_z0 --> g_z1 ----- START
c Copying initial guess
         nsub_file=0
         nsub=0
         m1=1
         m2=n
         p1=nsub+1
         p2=nsub+nvec
         do ipm=1,ncomp
c      write(*,799) m1,m2,p1,p2,ncomp,nvec,n,ipm
c 799  format('(m1,m2,p1,p2,ncomp,nvec,n,ipm)=(',
c     &         i12,',',i12,',',i12,',',i12,',',
c     &         i12,',',i12,',',i12,',',i12,')')
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
c         if (ga_nodeid().eq.0) then
c          write(*,2780) nblock,nsub
c 2780      format('Copying nblock=',i4,' g_z1--> g_z nsub=',i4)
c          endif
         m1=1
         m2=n
         p1=nsub+1
         p2=nsub+nvec
         do ipm=1,ncomp
c          if (ga_nodeid().eq.0) then
c           write(*,199) m1,m2,p1,p2,ncomp,nvec,n,ipm
c 199       format('copying: (m1,m2,p1,p2,ncomp,nvec,n,ipm)=(',
c     &         i12,',',i12,',',i12,',',i12,',',
c     &         i12,',',i12,',',i12,',',i12,')')
c          endif
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
c        if (ga_nodeid().eq.0)
c     &   write(*,*) 'BEF-prod'
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
c        if (ga_nodeid().eq.0)
c     &   write(*,*) 'AFT-prod'
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
c        if (ga_nodeid().eq.0)
c     &   write(*,*) 'FA BEF get_precond_rmax_zin'
c        if (ga_nodeid().eq.0) then
c        write(*,2777) iter,converge_precond
c 2777   format('ga_lkain_2cpl3-redmem:(iter,converge_precond)=(',
c     &         i3,',',L,')')
c        endif

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

c        if (ga_nodeid().eq.0) then
c         write(*,1782) flag2readfile,iter,cphf3write,rmax,tol,
c     &                 .not.(flag2readfile .and. iter.eq.1)
c 1782    format('(flag2readfile,iter,cphf3write,rmax,tol)=(',
c     &           L,',',i3,',',i2,',',f15.8,',',f15.8,
c     &           ',',L,')')
c        endif

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

c        if (ga_nodeid().eq.0)
c     &   write(*,*) 'FA AFT get_precond_rmax_zin'
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
c ========== g_z1 --> g_z ======== START
c Note.- Extract last (n1,nvec) block from g_z1 and
c        put it in g_z (this is the solution to Ax=b)
c         if (ga_nodeid().eq.0)
c     &    write(*,*) 'CONVERGENCE ACHIEVED!'
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
c ============ doing complex-linear solving ======= START
c         goto 1235
c        if (ga_nodeid().eq.0)
c     &   write(*,77) iter,nsub
c 77      format('BEFORE solve_zlineq1 iter=',i4,' nsub=',i4,')')
c         debug1=.true.

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

c        if (ga_nodeid().eq.0)
c     &   write(*,78) iter,nsub
c 78      format('AFTER solve_zlineq1 iter=',i4,' nsub=',i4,')')
c -------- get g_z from g_z1 --- START
      m1=1
      m2=n
      p1=nsub+1
      p2=nsub+nvec
      call ga_inquire(g_z1,type,dim1,dim2)
c      if (ga_nodeid().eq.0) then
c       write(*,79) m1,m2,p1,p2,nsub,dim1,dim2
c 79    format('(m1,m2,p1,p2,nsub,dim1,dim2=(',
c     &        i3,',',i3,',',i3,',',
c     &        i3,',',i3,',',i3,',',i3,')')
c      endif
      do ipm=1,ncomp
       call ga_copy_patch('n',g_z1    ,m1,m2,p1,p2,
     &                        g_z(ipm),1 ,n ,1 ,nvec)
       m1=m1+n
       m2=m2+n          
      enddo ! end-loop-ipm

c -------- get g_z from g_z1 --- START
c        if (ga_nodeid().eq.0)
c     &   write(*,*) 'AFTER solve_zlineq1 ...'
c ============ doing complex-linear solving ======= END
c        if (iter.eq.8) then
c         if (ga_nodeid().eq.0) 
c     &      write(*,*) 'FA-check-zlinear-solver-redmem'
c          stop
c        endif

        if (nsub .eq. maxsub-nvec) then
c       
c       Reduce the subspace as necessary
c       
c ====== FA: left-shifting patch ==== START
c Note.- matrices Ay,y shift to left nvec positions
c        removing leftmost patch of dimension: n4 x nvec
c         if (ga_nodeid().eq.0)
c     &    write(*,*) 'left-shifting, removing oldest subspace'
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
          write (luout,*) '  I will proceed, but check your results '//
     &                       'carefully!!!'
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
      end
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
c --> Experimental (not published yet)

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
     &         solve_zlineq1_sep_spins,
     &         dft_CPHF2_read2fix,
     &         precond_z1Az1,
     &         uhf_hessv_precond3_damp_anti,
     &         uhf_hessv_precond3_damp_check,
     &         ran1
      logical status,flag2readfile
      integer index4cphf,checkorth,cphf3write,
     &        csub
      character*255 aorespfilename
      character*(*) lbl_cphfaoresp
      character*255 lbl_cphfaoresp1
      integer iimoderaman,iiistepraman ! FA-06-18-12-raman defs
      logical ramanspc                 ! FA-06-18-12-raman defs
      parameter(lbl_cphfaoresp='aoresp_fiao_f') 
      one_cmplx =dcmplx( 1.0d0,0.0d0)
      mone_cmplx=dcmplx(-1.0d0,0.0d0)
      zero_cmplx=dcmplx( 0.0d0,0.0d0)
      if (.not. rtdb_get(rtdb, 'cphf:cphf3write',mt_int,1,
     &                   cphf3write))
     &  cphf3write = 0 ! assigns 0 if unsuccessfull read from rtdb
c      if (ga_nodeid().eq.0)
c     & write(*,*) 'In ga_lkain_2cpl3_damp_cmplx_redmem1: cphf3write=',
c     &           cphf3write

      if (.not. rtdb_get(rtdb, 'cphf:checkorth',mt_int,1,
     &                   checkorth))
     &  checkorth = 0 ! assigns 0 if unsuccessfull read from rtdb
c      if (ga_nodeid().eq.0)
c     & write(*,*) 'In ga_lkain_2cpl3_damp_cmplx_redmem1: checkorth=',
c     &           checkorth

      call ga_inquire(g_z(1),type,n,nvec) ! get (n,nvec) n=sum(nocc*nvirt(i) i=1,npol)

c      if (ga_nodeid().eq.0) then
c       write(*,14) n,ncomp,maxiter,mmaxsub,lifetime
c 14    format('(n,ncomp,maxiter,mmaxsub,lifetime)=(',
c     &          i3,',',i3,',',i5,',',i5,',',L1,')')
c      endif
      n1=ncomp*n
      maxsub = mmaxsub          ! So don't modify input scalar arg
      if (maxsub .lt. 3*nvec) maxsub = 3*nvec
      maxsub = (maxsub/nvec)*nvec 

c      if (ga_nodeid().eq.0) then
c       write(*,1023) n1,mmaxsub,maxsub,n
c 1023  format('(n1,mmaxsub,maxsub,n)=(',
c     &        i15,',',i5,',',i5,',',i15,')')
c      endif

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
c      if (.not. rtdb_get(rtdb,'cphf3-aores:guess1',
c     &          mt_int, 1,index4cphf)) call
c     $     errquit('ga_lkain_2cpl3_redmem: failed to read index4cphf ',
c     &             0, RTDB_ERR)

      if (.not. rtdb_get(rtdb,'cphf3-aores:guess1',
     & mt_int,1,index4cphf))  index4cphf = 0
c
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
c       write(*,112) iimoderaman,iiistepraman,lbl_cphfaoresp1
c 112   format('(iimoderaman,iiistepraman,lbl_cphfaoresp1)=(',
c     &         i4,',',i4,',',a,')')                  
      else
       write(lbl_cphfaoresp1,'(a13,i1)') lbl_cphfaoresp,index4cphf    
      endif
       call util_file_name(lbl_cphfaoresp1,
     &                     .false.,.false.,aorespfilename)

c      if (ga_nodeid().eq.0) then
c        write(*,*) 'ga_lkain_2cpl3_redmem: lbl_cphfaoresp1=',
c     &            lbl_cphfaoresp1,' aorespfilename=',
c     &            aorespfilename
c      endif  

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
          if (ga_nodeid().eq.0)
     &     write(*,1999) aorespfilename(1:inp_strlen(aorespfilename))
 1999      format('File ',a,
     &          ' does not exist, proceed to generate (z1,Az1)')
c ------ g_z0 --> g_z1 ----- START
c Copying initial guess
         nsub_file=0
         nsub=0
         m1=1
         m2=n
         p1=nsub+1
         p2=nsub+nvec
         do ipm=1,ncomp
c      write(*,799) m1,m2,p1,p2,ncomp,nvec,n,ipm
c 799  format('(m1,m2,p1,p2,ncomp,nvec,n,ipm)=(',
c     &         i12,',',i12,',',i12,',',i12,',',
c     &         i12,',',i12,',',i12,',',i12,')')
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
c Note.- After reading (g_z1,g_Az1) I need to use precond routine
c        which will do energy scaling.
c -------- Energy scaling --------- START
        goto 3703
        if (ga_nodeid().eq.0) then
         write(*,24) nsub,nsub_file
  24     format('(nsub,nsub_file)=(',i4,',',i4,')')    
        endif
        do csub=0,nsub,nvec
         if (ga_nodeid().eq.0) then
           write(*,23) csub
  23       format('E-scaling csub=',i4)
         endif

         call precond_z1Az1(
     &                  g_z1,    ! in : complex+accumulated g_zr
     &                  g_Az1,   ! in : complex+accumulated g_Az
     &                  csub,    ! in : pointer to current (g_z1,g_Az1)
     &                  precond, ! in : preconditioning    
     &                  omega,   ! in : omega
     &                  gamwidth,! in : gamwidth
     &                  ncomp,   ! in : nr. components
     &                  npol,    ! in : nr. polarizations
     &                  nvir,    ! in : nr. virtual  MOs
     &                  nocc,    ! in : nr. occupied MOs
     &                  n,       ! in : =nocc * nvir
     &                  nvec,    ! in : =3 (x,y,z)
     &                  debug)   ! in : =.true. -> allow debug printouts

        enddo ! end-loop-csub
 3703   continue
c -------- Energy scaling --------- END
         call ga_sync()
c Note.- I need to propagate nsub
         call ga_igop(6,nsub_file,1,'+') ! node0 nsub ne 0, eq 0 every other node    
         call ga_igop(6,nsub,1,'+')      ! node0 nsub ne 0, eq 0 every other node    
         flag2readfile=.true.
         nblock=nsub/3+1
c         if (ga_nodeid().eq.0) then
c          write(*,2780) nblock,nsub
c 2780      format('Copying nblock=',i4,' g_z1--> g_z nsub=',i4)
c         endif
         m1=1
         m2=n
         p1=nsub+1
         p2=nsub+nvec
         do ipm=1,ncomp
c          if (ga_nodeid().eq.0) then
c           write(*,199) m1,m2,p1,p2,ncomp,nvec,n,ipm
c 199       format('copying: (m1,m2,p1,p2,ncomp,nvec,n,ipm)=(',
c     &         i12,',',i12,',',i12,',',i12,',',
c     &         i12,',',i12,',',i12,',',i12,')')
c          endif
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
c        if (ga_nodeid().eq.0)
c     &   write(*,*) 'FA BEF get_precond_rmax_zin'
c        if (ga_nodeid().eq.0) then
c        write(*,2777) iter,converge_precond
c 2777   format('ga_lkain_2cpl3-redmem:(iter,converge_precond)=(',
c     &         i3,',',L,')')
c        endif

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
c        if (ga_nodeid().eq.0)
c     &   write(*,*) 'FA AFT get_precond_rmax_zin'
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
c ========== g_z1 --> g_z ======== START
c Note.- Extract last (n1,nvec) block from g_z1 and
c        put it in g_z (this is the solution to Ax=b)
c         if (ga_nodeid().eq.0)
c     &    write(*,*) 'CONVERGENCE ACHIEVED!'
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
c ============ doing complex-linear solving ======= START
c        if (ga_nodeid().eq.0)
c     &   write(*,77) iter,nsub
c 77      format('BEFORE solve_zlineq1 iter=',i4,' nsub=',i4,')')

          call solve_zlineq1(
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

c        if (ga_nodeid().eq.0)
c     &   write(*,78) iter,nsub
c 78      format('AFTER solve_zlineq1 iter=',i4,' nsub=',i4,')')
c -------- get g_z from g_z1 --- START
      m1=1
      m2=n
      p1=nsub+1
      p2=nsub+nvec
      call ga_inquire(g_z1,type,dim1,dim2)
c      if (ga_nodeid().eq.0) then
c       write(*,79) m1,m2,p1,p2,nsub,dim1,dim2
c 79    format('(m1,m2,p1,p2,nsub,dim1,dim2=(',
c     &        i3,',',i3,',',i3,',',
c     &        i3,',',i3,',',i3,',',i3,')')
c      endif
      do ipm=1,ncomp
       call ga_copy_patch('n',g_z1    ,m1,m2,p1,p2,
     &                        g_z(ipm),1 ,n ,1 ,nvec)
       m1=m1+n
       m2=m2+n          
      enddo ! end-loop-ipm

c -------- get g_z from g_z1 --- START
c        if (ga_nodeid().eq.0)
c     &   write(*,*) 'AFTER solve_zlineq1 ...'
c ============ doing complex-linear solving ======= END
c        if (iter.eq.18) then
c         if (ga_nodeid().eq.0) 
c     &      write(*,*) 'FA-check-zlinear-solver-redmem'
c          stop
c        endif

        if (nsub .eq. maxsub-nvec) then
c       
c       Reduce the subspace as necessary
c       
c ====== FA: left-shifting patch ==== START
c Note.- matrices Ay,y shift to left nvec positions
c        removing leftmost patch of dimension: n4 x nvec
c         if (ga_nodeid().eq.0)
c     &    write(*,*) 'left-shifting, removing oldest subspace'
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
          write (luout,*) '  I will proceed, but check your results '//
     &                       'carefully!!!'
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

      subroutine precond_z1Az1(
     &                  g_z1,    ! in : complex+accumulated g_zr
     &                  g_Az1,   ! in : complex+accumulated g_Az
     &                  nsub,    ! in : pointer to current (g_zr,g_Az)
     &                  precond, ! in : name of preconditioner or anti-preconditioner routine        
     &                  omega,   ! in : omega
     &                  gamwidth,! in : gamwidth
     &                  ncomp,   ! in : nr. components
     &                  npol,    ! in : nr. polarizations
     &                  nvir,    ! in : nr. virtual  MOs
     &                  nocc,    ! in : nr. occupied MOs
     &                  n,       ! in : =nocc * nvir
     &                  nvec,    ! in : =3 (x,y,z)
     &                  debug)   ! in : =.true. -> allow debug printouts
c
c Purpose:      Apply energy scaling  using      precond routine
c          Undo Apply energy scaling  using anti-precond routine
c
c Author : Fredy W. Aquino, Northwestern University
c Date   : 03-15-12
c --> Experimental (not published yet)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
      integer ipm,ncomp,i,
     &        npol,nvir(npol),nocc(npol)
      double precision omg(ncomp),
     &                 omega,gamwidth
      integer nsub,n,nvec,
     &        g_z1,g_Az1,
     &        g_dre,g_dim ! scratch GA used for (g_r or g_Ax)
      logical debug
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
        omg(1)=-omega
        omg(2)= omega
        do ipm=1,ncomp
          call conv2reim1_u(
     &                    g_dre,  ! out : real      arr
     &                    g_dim,  ! out : imaginary arr
     &                    g_z1,   ! in  : = complx(g_xre,g_xim)
     &                    nsub,   ! in  : pointer to block
     &                    npol,   ! in  : nr. polarizations
     &                    nvir,   ! in  : nr. virtual  MOs
     &                    nocc,   ! in  : nr. occupied MOs
     &                    ipm,    ! in  : =1,2 components indices
     &                    n,      ! in  : n    rows
     &                    nvec)   ! in  : nvec columns
c         if (ga_nodeid().eq.0)
c     &    write(*,*) '-------BEF-precod(',ipm,') ----- START'
c          call ga_print(g_dre)
c          call ga_print(g_dim)
c         if (ga_nodeid().eq.0)
c     &    write(*,*) '-------BEF-precod(',ipm,') ----- END'

          call precond(g_dre,g_dim,omg(ipm),gamwidth)

c         if (ga_nodeid().eq.0)
c     &    write(*,*) '-------AFT-precod(',ipm,') ----- START'
c          call ga_print(g_dre)
c          call ga_print(g_dim)
c         if (ga_nodeid().eq.0)
c     &    write(*,*) '-------AFT-precod(',ipm,') ----- END'

          call conv2complex1_u(
     &                    g_z1, ! out: = complx(g_xre,g_xim)
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
          call precond(g_dre,g_dim,omg(ipm),gamwidth)
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
      return
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
c --> Experimental (not published yet)

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
          call ga_get(g_xre(ipm),1,n,ivec,ivec,dbl_mb(k_zre),1)
          call ga_get(g_xim(ipm),1,n,ivec,ivec,dbl_mb(k_zim),1)
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
c --> Experimental (not published yet)

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
     $                 dbl_mb(k_zre),1)
           call ga_get(g_xim,ioff,ioff+nvir-1,ivec1,ivec1,
     $                 dbl_mb(k_zim),1)
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
c --> Experimental (not published yet)

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
     $                  dbl_mb(k_zre),1)
            call ga_get(g_xim,ioff,ioff+nvir(ipol)-1,ivec1,ivec1,
     $                  dbl_mb(k_zim),1)
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
c --> Experimental (not published yet)

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
c ---- Alg-1-START
        call ga_get(g_xreim,1,n,ivec,ivec,dbl_mb(k_xreim),1)
        call ga_get(g_z    ,1,n,ivec,ivec,dcpl_mb(k_z),1)
        do idat=1,n
         val_cmplx=dcmplx(dbl_mb(k_xreim+idat-1),0.0d0)
         call ga_put(g_a,1,1,1,1,val_cmplx,1)
         call ga_add_patch(one_cmplx,g_z,idat,idat,ivec,ivec,
     &                     one_cmplx,g_a,1   ,1   ,1   ,1   ,
     &                               g_z,idat,idat,ivec,ivec)
        enddo ! end-loop-idat
c ---- Alg-1-END
c ---- Alg-2-START
c --> It needs (nocc,nvir) info
c       do i = ga_nodeid()+1,nocc,ga_nnodes()
c        ioff = (i-1)*nvir + 1     
c        call ga_get(g_xreim,ioff,ioff+nvir-1,ivec,ivec,
c     $              dbl_mb(k_xreim),1)
c        call ga_get(g_z,ioff,ioff+nvir-1,ivec,ivec,
c     $              dbl_mb(k_z),1)
c        do j=1,nvir
c         val_cmplx=dcmplx(dbl_mb(k_xreim+j-1),0.0d0)
c          call ga_put(g_a,1,1,1,1,val_cmplx,1)
c          call ga_add_patch(one_cmplx,g_z,j,j,ivec,ivec,
c     &                      one_cmplx,g_a,1,1,1   ,1   ,
c     &                                g_z,j,j,ivec,ivec)   
c        enddo ! end-loop-j
c       enddo ! end-loop-i
c ---- Alg-2-END
       enddo ! end-loop-ivec
      else if (indrm.eq.2) then ! updating only IMAG part
       do ivec=1,nvec
        call ga_get(g_xreim,1,n,ivec,ivec,dbl_mb(k_xreim),1)
        call ga_get(g_z    ,1,n,ivec,ivec,dcpl_mb(k_z),1)
          do idat=1,n
           val_cmplx=dcmplx(0.0d0,dbl_mb(k_xreim+idat-1))
           call ga_put(g_a,1,1,1,1,val_cmplx,1)
           call ga_add_patch(one_cmplx,g_z,idat,idat,ivec,ivec,
     &                       one_cmplx,g_a,1   ,1   ,1   ,1   ,
     &                                 g_z,idat,idat,ivec,ivec)
          enddo ! end-loop-idat
c ---- Alg-2-START
c --> It needs (nocc,nvir) info
c       do i = ga_nodeid()+1,nocc,ga_nnodes()
c        ioff = (i-1)*nvir + 1     
c        call ga_get(g_xreim,ioff,ioff+nvir-1,ivec,ivec,
c     $              dbl_mb(k_xreim),1)
c        call ga_get(g_z,ioff,ioff+nvir-1,ivec,ivec,
c     $              dbl_mb(k_z),1)
c        do j=1,nvir
c          val_cmplx=dcmplx(0.0d0,dbl_mb(k_xreim+j-1))
c          call ga_put(g_a,1,1,1,1,val_cmplx,1)
c          call ga_add_patch(one_cmplx,g_z,j,j,ivec,ivec,
c     &                      one_cmplx,g_a,1,1,1   ,1   ,
c     &                                g_z,j,j,ivec,ivec)   
c        enddo ! end-loop-j
c       enddo ! end-loop-i
c ---- Alg-2-END
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
c --> Experimental (not published yet)

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
     $               dbl_mb(k_xreim),1)
         call ga_get(g_z,ioff,ioff+nvir-1,ivec,ivec,
     $               dcpl_mb(k_z),1)
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
     $               dbl_mb(k_xreim),1)
         call ga_get(g_z,ioff,ioff+nvir-1,ivec,ivec,
     $               dcpl_mb(k_z),1)
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
c --> Experimental (not published yet)

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
     $               dbl_mb(k_xreim),1)
         call ga_get(g_z,ioff2,ioff2+nvir-1,ivec,ivec,
     $               dcpl_mb(k_z),1)
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
     $               dbl_mb(k_xreim),1)
         call ga_get(g_z,ioff2,ioff2+nvir-1,ivec,ivec,
     $               dcpl_mb(k_z),1)
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
c --> Experimental (not published yet)

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
     $               dbl_mb(k_xreim),1)
         call ga_get(g_z,ioff2,ioff2+nvir-1,ivec,ivec,
     $               dcpl_mb(k_z),1)
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
     $               dbl_mb(k_xreim),1)
         call ga_get(g_z,ioff2,ioff2+nvir-1,ivec,ivec,
     $               dcpl_mb(k_z),1)
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
c --> Experimental (not published yet)

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
     $                dbl_mb(k_xreim),1)
          call ga_get(g_z,ioff2,ioff2+nvir(ipol)-1,ivec,ivec,
     $                dcpl_mb(k_z),1)
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
     $                dbl_mb(k_xreim),1)
          call ga_get(g_z,ioff2,ioff2+nvir(ipol)-1,ivec,ivec,
     $                dcpl_mb(k_z),1)
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
c --> Experimental (not published yet)

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
        call ga_get(g_xreim,1,nvir,i,i,dbl_mb(k_xreim),1)  
        call ga_get(g_z,ioff,ioff+nvir-1,ivec,ivec,
     &              dcpl_mb(k_z),1)

c        write(*,2) (dbl_mb(k_xreim+j-1),j=1,nvir)
c 2      format(1x,f15.8)
c        write(*,3) (dcpl_mb(k_z+j-1),j=1,nvir)
c 3      format(1x,'(',f15.8,',',f15.8,')')

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
        call ga_get(g_xreim,1,nvir,i,i,dbl_mb(k_xreim),1)  
        call ga_get(g_z,ioff,ioff+nvir-1,ivec,ivec,
     &              dcpl_mb(k_z),1)

c        write(*,*) 'IM iocc=',i
c        write(*,4) (dbl_mb(k_xreim+j-1),j=1,nvir)
c 4      format(1x,f15.8)
c        write(*,5) (dcpl_mb(k_z+j-1),j=1,nvir)
c 5      format(1x,'(',f15.8,',',f15.8,')')

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
c --> Experimental (not published yet)

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
        call ga_get(g_xreim,1,nvir,i,i,dbl_mb(k_xreim),1)  
        call ga_get(g_z,ioff,ioff+nvir-1,ivec1,ivec1,
     &              dcpl_mb(k_z),1)
        do j=1,nvir
          val_zre=dreal(dcpl_mb(k_z+j-1))
          val_zim=dimag(dcpl_mb(k_z+j-1))
          val_cmplx=dcmplx(val_zre+scl*dbl_mb(k_xreim+j-1),val_zim)
          j1=ioff+j-1
          call ga_put(g_z,j1,j1,ivec1,ivec1,val_cmplx)
        enddo ! end-loop-j
       enddo ! end-loop-i
      else if (indrm.eq.2) then ! updating only IMAG part
       do i = ga_nodeid()+1,nocc,ga_nnodes()
        ioff = ioff1+(i-1)*nvir + 1   
        call ga_get(g_xreim,1,nvir,i,i,dbl_mb(k_xreim),1)  
        call ga_get(g_z,ioff,ioff+nvir-1,ivec1,ivec1,
     &              dcpl_mb(k_z),1)
        do j=1,nvir
          val_zre=dreal(dcpl_mb(k_z+j-1))
          val_zim=dimag(dcpl_mb(k_z+j-1))
          val_cmplx=dcmplx(val_zre,val_zim+scl*dbl_mb(k_xreim+j-1))
          j1=ioff+j-1
          call ga_put(g_z,j1,j1,ivec1,ivec1,val_cmplx)
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
c --> Experimental (not published yet)

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
          call ga_get(g_z(ipm),1,n,ivec,ivec,dcpl_mb(k_z),1)
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
c --> Experimental (not published yet)

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
c --> Experimental (not published yet)

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
     &                 dcpl_mb(k_z),1)
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
c --> Experimental (not published yet)

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
     &                  dcpl_mb(k_z),1)
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
c --> Experimental (not published yet)

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
     &               dcpl_mb(k_z),1)
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
     $               dcpl_mb(k_z),1)
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
c --> Experimental (not published yet)

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
     &               dcpl_mb(k_z),1)
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
     $               dcpl_mb(k_z),1)
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
c --> Experimental (not published yet)

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
     &                dcpl_mb(k_z),1)
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
     $                dcpl_mb(k_z),1)
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
c --> Experimental (not published yet)

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
       ioff1=(ipm-1)*n
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
c --> Experimental (not published yet)

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
c --> Experimental (not published yet)

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
c --> Experimental (not published yet)

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
c --> Experimental (not published yet)

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
 10         format('---------g_xr1-1(',i3,')-----START')
          endif
          call ga_print(g_xr1)
          if (ga_nodeid().eq.0) then
           write(*,6) iter
 6         format('---------g_xr1-1(',i3,')-----END')
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
 2         format('---------g_x-1(',i3,',',i3,')-----START')
          endif
          call ga_print(g_x(ipm))
          if (ga_nodeid().eq.0) then
           write(*,3) ipm,iter
 3         format('---------g_x-1(',i3,',',i3,')-----END')
          endif
         enddo ! end-loop-ipm
        endif ! end-if-debug1
        call ga_zero(g_xr1)
        call ga_dgemm('n','n',n1,nvec,nsub,one,
     &                g_x1,g_c,zero,g_xr1)

        if (debug1) then
          if (ga_nodeid().eq.0) then
           write(*,7) iter
 7         format('---------g_xr1-2(',i3,')-----START')
          endif
          call ga_print(g_xr1)
          if (ga_nodeid().eq.0) then
           write(*,9) iter
 9         format('---------g_xr1-2(',i3,')-----END')
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
 4         format('---------g_x-2(',i3,',',i3,')-----START')
          endif
          call ga_print(g_x(ipm))
          if (ga_nodeid().eq.0) then
           write(*,5) ipm,iter
 5         format('---------g_x-2(',i3,',',i3,')-----END')
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
c --> Experimental (not published yet)

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
 10         format('---------g_zr1-1(',i3,')-----START')
          endif
          call ga_print(g_zr1)
          if (ga_nodeid().eq.0) then
           write(*,6) iter
 6         format('---------g_zr1-1(',i3,')-----END')
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
 2         format('---------g_z-1(',i3,',',i3,')-----START')
          endif
          call ga_print(g_z(ipm))
          if (ga_nodeid().eq.0) then
           write(*,3) ipm,iter
 3         format('---------g_z-1(',i3,',',i3,')-----END')
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
 7         format('---------g_zr1-2(',i3,')-----START')
          endif
          call ga_print(g_zr1)
          if (ga_nodeid().eq.0) then
           write(*,9) iter
 9         format('---------g_zr1-2(',i3,')-----END')
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
 4         format('---------g_z-2(',i3,',',i3,')-----START')
          endif
          call ga_print(g_z(ipm))
          if (ga_nodeid().eq.0) then
           write(*,5) ipm,iter
 5         format('---------g_z-2(',i3,',',i3,')-----END')
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
c -------- FA-04-21-12: complex solver reducing memory cost --- START
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
c --> Experimental (not published yet)

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
 10      format('-------g_a(',i4,')--------START')
        call ga_print(g_a)
        if (ga_nodeid().eq.0)
     &   write(*,*) iter
 11      format('-------g_a(',i4,')--------END')
        if (ga_nodeid().eq.0)
     &   write(*,12) iter
 12      format('-------g_b(',i4,')--------START')
        call ga_print(g_b)
        if (ga_nodeid().eq.0)
     &   write(*,13) iter
 13      format('-------g_b(',i4,')--------END')
        endif ! end-if-debug
      
c        if (ga_nodeid().eq.0)
c     &   write(*,*) 'BEF ga_svd_solve_seq_cmplx'

        call ga_svd_solve_seq_cmplx(g_a,g_b,g_c,1d-14)

c        if (ga_nodeid().eq.0)
c     &   write(*,*) 'AFT ga_svd_solve_seq_cmplx'

        if (debug1) then
        if (ga_nodeid().eq.0)
     &   write(*,14) iter
 14      format('-------g_c(',i4,')--------START')
        call ga_print(g_c)
        if (ga_nodeid().eq.0)
     &   write(*,15) iter
 15      format('-------g_c(',i4,')--------END')
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
 20         format('---------g_zr1-1(',i3,')-----START')
          endif
          call ga_print(g_zr1)
          if (ga_nodeid().eq.0) then
           write(*,6) iter
 6         format('---------g_zr1-1(',i3,')-----END')
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
 2         format('---------g_z1-1(',i3,')-----START')
          endif
          call ga_print(g_z1)
          if (ga_nodeid().eq.0) then
           write(*,3) iter
 3         format('---------g_z1-1(',i3,')-----END')
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
 7         format('---------g_zr1-2(',i3,')-----START')
          endif
          call ga_print(g_zr1)
          if (ga_nodeid().eq.0) then
           write(*,9) iter
 9         format('---------g_zr1-2(',i3,')-----END')
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
 4         format('---------g_z1-2(',i3,')-----START')
          endif
          call ga_print(g_z1)
          if (ga_nodeid().eq.0) then
           write(*,5) iter
 5         format('---------g_z1-2(',i3,')-----END')
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
c --> Experimental (not published yet)

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
 1       format('-------z1^t*(z1c)(',i4,')--------START')
        do i=1,nsub
         call ga_get(g_b,i,i,1,nvec,dcpl_mb(k_z),1)
         if (ga_nodeid().eq.0) then
          write(*,'(i8,":",3(" (",f14.8,",",f14.8,") "))') 
     &      i,(dreal(dcpl_mb(k_z+j-1)),dimag(dcpl_mb(k_z+j-1)),j=1,nvec)
         endif
        enddo ! end-loop-idata  
        if (ga_nodeid().eq.0)
     &   write(*,2) iter
 2       format('-------z1^t*(z1c)(',i4,')--------END')
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
c --> Experimental (not published yet)

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
 16        format('---------z1^t*(z1c)(',i3,')-----START')
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
 17        format('---------z1^t*(z1c)(',i3,')-----END')
          endif
      if (.not.ma_pop_stack(l_z))
     $  call errquit('toview_orthz1c: pop problem with l_zim',
     &               555,MA_ERR)
       return
       end
c -------- FA-04-21-12: complex solver reducing memory cost --- END

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
c --> Experimental (not published yet)

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

c ---------- Examine output --------- START
c       if (ga_nodeid().eq.0)
c     &  write(*,*) '-------g_u --------- START'
c       call ga_print(g_u)
c       if (ga_nodeid().eq.0)
c     &  write(*,*) '-------g_u --------- END'
c       if (ga_nodeid().eq.0)
c     &  write(*,*) '-------g_vt --------- START'
c       call ga_print(g_vt)
c       if (ga_nodeid().eq.0)
c     &  write(*,*) '-------g_vt --------- END'
c      do i = 0, nsing-1
c       write(*,1) i+1,dbl_mb(k_val+i)
c  1    format('(i+1,w)=(',i3,',',f15.8,')')
c      enddo ! end-loop-i
c      stop
c ---------- Examine output --------- END
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
c --> Experimental (not published yet)

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

c       write(*,1) m,n,lwork,nsing
c 1     format('ga_svd_seq_cmplx:(m,n,lwork,nsing)=(',
c     &        i3,',',i3,',',i3,',',i3,')')

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
c --> Experimental (not published yet)

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
c --> Experimental (not published yet)

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
c --> Experimental (not published yet)

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
c --> Experimental (not published yet)

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
c --> Experimental (not published yet)

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
c --> Experimental (not published yet)
                
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
c --> Experimental (not published yet)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
      integer iter,ipm,cc,ncomp
      double precision rmax,omg(ncomp),
     &                 omega,rmx(ncomp)
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
c --> Experimental (not published yet)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
      integer iter,ipm,cc,ncomp
      double precision rmax,omg(ncomp),
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
c --> Experimental (not published yet)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
      logical debug
      integer i,iter,ipm,cc,n,ncomp,
     &        l_dat,k_dat
      double precision g_re(ncomp),g_im(ncomp),
     &                 rmax
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
  1       format('(iter,rmax,rmax_arr)=(',
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
c --> Experimental (not published yet)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "util.fh"
      integer iter,ipm,cc,ncomp,nmx,i,
     &        npol,nvir(npol),nocc(npol)
      double precision rmax,omg(ncomp),
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
            call precond(g_dre,g_dim,omg(ipm),gamwidth)
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
           call precond(g_dre,g_dim,omg(ipm),gamwidth)
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
  1       format('(iter,rmax,rmax_arr)=(',
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
          call precond(g_dre,g_dim,omg(ipm),gamwidth)
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
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +++++++++++++++ device anti-'preconditioning +++++++++++ START
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c Note 1.- To undo the effect of 'preconditioning' and to
c        see if z(1) and z(2) are different only because of
c        that transformation [z(ipm) ipm=1,ncomp ncomp,nr of components]
c        We want to solve only for one of the components and obtain
c        the other from a closed formula because we think that
c        z(1) and z(2) are dependent each other. FA-04-12-12
c Note 2.- 'preconditioning' refers to energy weighting z vector
c        If z=z_r + i z_i
c        z_r1 = a_r z_r + a_i z_i
c        z_i1 = a_i z_r + a_r z_i
c        with a_r = d_aj/D  a_i=-g/D   D=d_aj**2 + g**2
c       d_aj=4.d0*(fcv(j+vbase) - fcv(i+cbase)) + lshift
c            --> energy difference bet. occupied and virtual KS states
c                it includes a shifting too (lshift)
c       g   = gamwidth
c       By anti-preconditioning we want:
c       (z_r,z_i) from (z_r1,z_i1) upon an inverse transformation:
c       If z1 = A z --> z= A^{-1} z1
c       A^{-1}= [ a_r -a_i; -a_i a_r ]/ D1   D1=a_r**2 - a_i**2
c       This procedure relies on: D1 ne 0 otherwise means
c       linear dependency:  a_r= +/- a_i
      subroutine rohf_hessv_precond3_damp_anti(
     &                        g_x, 
     &                        g_x_im, 
     &                        omega, 
     &                        gamwidth)
c
c ... jochen: modified version of rohf_hessv_precond to go along with
c     cphf_solve3. This preconditioner PRESUMABLY sets up the rhs of
c     the CPKS problem which is where we should add and subtract
c     the frequency from the energy denominator.
c     NOTE: this is a modified version of rohf_hessv_precond3 for
c     calculations with finite-lifetime damping
c
c Author: Fredy W. Aquino, Northwestern University
c Date  : 03-15-12
c Note  : modifying rohf_hessv_precond3()
c --> Experimental (not published yet)

      implicit none
#include "errquit.fh"
#include "crohf.fh"
#include "cscf.fh"
c     
c     $Id: rohf_wrap.F 21915 2012-02-01 22:09:20Z d3y133 $

      integer g_x               ! [in/ou] Re rhs vector of eq. system
      integer g_x_im            ! [in/ou] Im rhs vector of eq. system
      double precision omega    ! [input] frequency
      double precision gamwidth ! [input] damping constant
c     
      integer gtype, grow, gcol
      double precision crohf_lshift, crohf_gamwidth
c     
c     Check
c     
      if (crohf_init_flag.ne.1)
     $     call errquit(
     $     'rohf_hessv_precond: ROHF internal block invalid',0,
     &       UNKNOWN_ERR)
      call ga_inquire(g_x,gtype,grow,gcol)
      if (grow.ne.crohf_vlen)
     $     call errquit('rohf_hessv_precond: invalid vector length',0,
     &       UNKNOWN_ERR)
c     
c ... jochen: in contrast to rohf_hessv_precond we add
c     omega here and use the result in the call to
c     rohf_hdiag_scale. Note: in rohf_hdiag_scale there is a division
c     by [ 4 (epsilon_a - epsilon_i) + omega ] so we should use
c     four times omega here to be consistent with the TDDFT formalism.
c
c ... jochen: Note:    
c     level shifting might have some strange effects here 
c     
      crohf_lshift = lshift + 4d0 * omega
c
c     use 4 * lifetime parameter to be consistent with
c     using 4 * omega
      crohf_gamwidth = 4d0 * gamwidth

      call rohf_hdiag_scale_damp_anti(
     &          nmo, nclosed,nopen,
     $          crohf_lshift,crohf_gamwidth, 
     $          crohf_g_fcv,crohf_g_fpv,crohf_g_fcp,
     $          g_x,g_x_im)
      end

       subroutine rohf_hdiag_scale_damp_anti(nmo, 
     &                                       nclosed, 
     &                                       nopen, 
     &                                       lshift,
     &                                       gamwidth,
     $                                       g_fcv, 
     &                                       g_fpv, 
     &                                       g_fcp, 
     &                                       g_x, 
     &                                       g_x_im)
C     $Id: rohf_hscale.F 19707 2010-10-29 17:59:36Z d3y133 $
c
c ... jochen: modified version of rohf_hdiag_scale to allow for
c     a situation where we have damping and a coupled set of equations.
c 
c     further, the original routine contained the check: 
c
c     denominator = 4.d0*(fcv(j+vbase) - fcv(i+cbase)) +
c    $                 lshift
c     if (denominator .lt. 0.1d0) denominator = 0.1d0
c     dbl_mb(k_x+j-1) = dbl_mb(k_x+j-1) / denominator
c
c     where I changed the number to 0.00001 instead, and to print
c     a warning! (ONLY C-V PIECE)
c
c     NOTE: THIS ROUTINE PRESENTLY WORKS ONLY FOR A CLOSED SHELL SYSTEM
c
c Author: Fredy W. Aquino, Northwestern University
c Date  : 03-15-12
c Note  : modifying rohf_hdiag_scale()
c --> Experimental (not published yet)

      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
c     
      integer nmo, nclosed, nopen
      double precision lshift, gamwidth, small, sign
      parameter (small = 1d-5)
      integer g_fcv, g_fpv, g_fcp
      integer g_x, g_x_im
      logical warning
c     
c     
c     Scale the vector x by the ROHF 1e-Hessian diagonal
c     
c               x                    x
c                pq                   ia
c     x   =  --------           = ---------    [ closed-virtual ]
c      pq      (1e)                 cv    cv
c             H                   4F  - 4F
c              pq,pq                aa    ii
c     
c               x
c                pa
c         = ---------    [ partial-virtual ]
c             pv    pv
c           2F  - 2F
c             aa   pp
c     
c               x
c                ip
c         = ---------    [ closed-partial ]
c             cp    cp
c           2F  - 2F
c             pp    ii
c     
      integer nvir
      integer i, j, ioff, ioff1,cbase, obase, vbase
      integer l_fcv, k_fcv
      integer l_fcp, k_fcp
      integer l_fpv, k_fpv
      integer l_x, k_x
c ... jochen: we need some memory for the imaginary part, too
      integer l_x_im, k_x_im
      double precision fcv, fcp, fpv, denominator, 
     &                 factor_re,factor_im, 
     &                 tmpre,tmpim,val,val1,rval
      integer ivec, nvec, gtype, vlen
      fcv(i) = dbl_mb(k_fcv + i - 1) ! Statement function definitions
      fpv(i) = dbl_mb(k_fpv + i - 1)
      fcp(i) = dbl_mb(k_fcp + i - 1)
c
c     =================================================================
      warning = .false.
      call ga_inquire(g_x,gtype,vlen,nvec)    
      nvir = nmo - nclosed - nopen
      cbase = 0                 ! Offsets into diagonals
      obase = nclosed
      vbase = nclosed + nopen
c     
c     Get diagonals of each matrix

      if (.not.ma_push_get(MT_DBL,nmo,'rohf: fcv',l_fcv,k_fcv))
     $     call errquit('hdiag_scale_3im: cannot allocate',0, MA_ERR)
      if (.not.ma_push_get(MT_DBL,nmo,'rohf: x',l_x,k_x))
     $     call errquit('hdiag_scale_3im: cannot allocate',0, MA_ERR)
      if (.not.ma_push_get(MT_DBL,nmo,'rohf: x_im',l_x_im,k_x_im))
     $     call errquit('hdiag_scale_3im: cannot allocate',0, MA_ERR)
     
      call ga_get_diagonal(g_fcv, dbl_mb(k_fcv))
c
c     Scale closed-virtual piece
      if (nvir*nclosed .gt. 0) then
       do ivec = 1, nvec
         do i = ga_nodeid()+1, nclosed, ga_nnodes()
           ioff = (i-1)*nvir + 1
           ioff1= ioff+nvir-1
c          put elements of g_x, g_x_im to local memory:
           call ga_get(g_x   ,ioff,ioff1,ivec,ivec,dbl_mb(k_x)   ,1)
           call ga_get(g_x_im,ioff,ioff1,ivec,ivec,dbl_mb(k_x_im),1)
           do j=1,nvir
c            orbital energy differences:
             denominator = 4.d0*(fcv(j+vbase)-fcv(i+cbase))+lshift
c            check for small HOMO-LUMO gap:
             if (abs(denominator).lt.small) then
               sign = 1d0
               if (denominator.lt.0d0) sign = -1d0
               denominator = small * sign
               if (ga_nodeid().eq.0 .and. .not.warning)
     &            write (LuOut,9000) 
               warning = .true.
             endif             
c            take care of the damping here: set up scale_factor
c            and multiply instead of dividing by denominator
             val =(denominator**2 + gamwidth**2)
             val1=(denominator**2 - gamwidth**2)
             rval=val/val1
             factor_re = denominator*rval
             factor_im = gamwidth   *rval
c            scale real part, mix with imag.; and vice versa
             tmpre = dbl_mb(k_x   +j-1)
             tmpim = dbl_mb(k_x_im+j-1)
             dbl_mb(k_x   +j-1)= tmpre*factor_re + tmpim*factor_im
             dbl_mb(k_x_im+j-1)= tmpre*factor_im + tmpim*factor_re 
           enddo ! end-loop-j (un-occupied MOs)
           call ga_put(g_x   ,ioff,ioff1,ivec,ivec,dbl_mb(k_x),1)
           call ga_put(g_x_im,ioff,ioff1,ivec,ivec,dbl_mb(k_x_im),1)
         enddo ! end-loop-i (occupied MOs)
       enddo ! end-loop-ivec
      endif  ! end-if-nvir-nclosed 
      if (.not.ma_pop_stack(l_x_im))
     $     call errquit('rohf_hdiag_scale_im3: cannot pop stack xim',
     &   0, MA_ERR)
      if (.not.ma_pop_stack(l_x))
     $     call errquit('rohf_hdiag_scale_im3: cannot pop stack lx',
     &   0, MA_ERR)
      if (.not.ma_pop_stack(l_fcv))
     $     call errquit('rohf_hdiag_scale_im3: cannot pop stack fcv',
     &   0, MA_ERR)

      return
c
 9000 format (/' *** WARNING: SMALL OCC-VIRT ORBITAL ENERGY DIFFERENCE',
     &   ' DETECTED'/' (Including frequency factor for response)'/
     &   ' This is going to lead to diverging respose properties'/
     &   ' I have adjusted the denominator to a value of 1E-5')

      end
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +++++++++++++++ device anti-'preconditioning +++++++++++ END
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
