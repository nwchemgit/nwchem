      subroutine ga_lkain_2cpl3(rtdb, 
     &                          g_x, g_b, 
     &                          g_x_im, g_b_im,
     &                          product, precond, 
     $                          tol, mmaxsub, maxiter, 
     &                          odiff, oprint, omega, limag,
     &                          lifetime, gamwidth, ncomp)

c     $Id$

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "stdio.fh"
#include "rtdb.fh"
c
      integer rtdb              ! [input] database handle
      integer g_x(2)            ! [input/output] Initial guess/solution
      integer g_x_im(2)         ! not used
      integer g_b(2)            ! [input] Right-hand side vectors 
      integer g_b_im(2)         ! not used
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
      integer g_y, g_Ay, 
     &        g_Ax(2), g_r(2), g_r2, 
     &        g_a, g_bb,g_c, 
     &        g_xold(2), g_Axold(2), g_Ax_im(2)
      double precision rmax, rmax1, rmax2, acc,rmx(2)
      logical converged
      logical odebug, debug, converge_precond

      double precision omg(2)
      integer dsp
c
c     =================================================================
c
      debug = (.false. .and. ga_nodeid().eq.0) ! for code development

c      debug=.true.

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
c
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

        call ga_sync()
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
          call ga_sync()
c
c         g_Ax = g_b if the system is solved. During the first cycle,
c         g_Ax is calculated from the initial guess
          call ga_add(1.0d0, g_b(ipm),
     &               -1.0d0, g_Ax(ipm), 
     &                       g_r(ipm)) ! The residual
          call ga_sync()        
        enddo                   ! ipm = 1,ncomp
        call ga_sync()

c       convergence checking:
c       find the largest element of the residual either 
c       before or after the call to the preconditioner
         omg(1)=-omega
         omg(2)= omega

         if (converge_precond) then
          do ipm=1,ncomp
           call precond(g_r(ipm),omg(ipm))
          enddo ! end-loop-ipm
         endif ! end-if-converge_precond
         rmx(1)=0.0d0
         rmx(2)=0.0d0
          do ipm=1,ncomp
           call ga_maxelt(g_r(ipm),rmx(ipm))   
          enddo ! end-loop-ipm

         rmax1=rmx(1)
         rmax2=rmx(2)

         rmax = max(rmax1, rmax2)  
      
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
c       
c ... jochen: changed 0d0 to omega in the calls below
c       (there were only two in ga_lkain, one with g_aX, one with g_r)
c       for array g_r the preconditioner call is only necessary
c       in case converge_precond is .false.
        omg(1)=-omega
        omg(2)= omega
        do ipm=1,ncomp
         call precond(g_Ax(ipm),omg(ipm))
         if (.not.converge_precond) 
     &   call precond(g_r(ipm) ,omg(ipm))
        enddo ! end-loop-ipm
  
        call ga_sync()
        
c       Copy the vectors to the subspace work area      
        do ipm=1,ncomp
         dsp=n*(ipm-1)
         call ga_copy_patch('n', 
     $                      g_Ax(ipm),1    ,n    ,1     ,nvec, 
     $                      g_Ay     ,1+dsp,n+dsp,nsub+1,nsub+nvec)
         call ga_copy_patch('n', 
     $                      g_x(ipm),1    ,n    ,1     ,nvec, 
     $                      g_y     ,1+dsp,n+dsp,nsub+1,nsub+nvec)
         call ga_copy_patch('n', 
     $                      g_r(ipm),1    ,n    ,1,nvec, 
     $                      g_r2    ,1+dsp,n+dsp,1,nvec)
        enddo ! end-loop-ipm
        
        nsub = nsub + nvec
c       
c       Form and solve the subspace equations using SVD in order
c       to manage near linear dependence in the subspace.
c       
        if (.not. ga_create(MT_DBL, nsub, nsub,
     &     'lkain_2cpl: A', 0, 0, g_a))
     $     call errquit('lkain: allocating g_a?', nsub, GA_ERR)
        if (.not. ga_create(MT_DBL, nsub, nvec,
     &     'lkain_2cpl: B', 0, 0,g_bb))
     $     call errquit('lkain: allocating g_bb?', nsub, GA_ERR)
        if (.not. ga_create(MT_DBL, nsub, nvec,
     &     'lkain_2cpl: C', 0, 0, g_c))
     $     call errquit('lkain: allocating g_c?', nsub, GA_ERR)
        call ga_zero(g_a)
        call ga_zero(g_bb)
        call ga_zero(g_c)
        call ga_sync()
        call ga_dgemm('t','n',nsub,nsub,n2,1.0d0,
     &                g_y,g_Ay,0.0d0,g_a)
        call ga_dgemm('t','n',nsub,nvec,n2,1.0d0,
     &                g_y,g_r2,0.0d0,g_bb)
        call ga_sync()
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
c       
        call ga_svd_solve_seq(g_a,g_bb,g_c,1d-14)
        if (odebug) call ga_print(g_c)
c       
c       Form and add the correction, in parts, onto the solution
c       
        call ga_sync()
        call ga_dgemm('n','n',n2,nvec,nsub,-1.0d0,
     &                g_Ay,g_c,1.0d0,g_r2)
        if (odebug) then
          write(6,*) ' The update in the complement '
          call ga_print(g_r2)
        end if
c       
c       copy components of g_r2 into g_r before adding g_r to  g_x
        do ipm=1,ncomp
         dsp=n*(ipm-1)
         call ga_copy_patch('n', 
     $                     g_r2    ,1+dsp,n+dsp,1,nvec, 
     $                     g_r(ipm),1    ,n    ,1,nvec)      
        enddo ! end-loop-ipm
        call ga_sync()       
        do ipm = 1,ncomp
          call ga_add(1.0d0, g_r(ipm), 
     &                1.0d0, g_x(ipm), 
     &                       g_x(ipm))
        enddo ! end-loop-ipm
        
        call ga_sync()
        call ga_dgemm('n','n',n2,nvec,nsub,1.0d0,
     &                g_y,g_c,0.0d0,g_r2)

        call ga_sync()
        if (odebug) then
          write(6,*) ' The update in the subspace '
          call ga_print(g_r2)
        end if
c       
c       copy components of g_r2 into g_r before adding g_r to  g_x    
        do ipm=1,ncomp
         dsp=n*(ipm-1)
         call ga_copy_patch('n', 
     $                     g_r2    ,1+dsp,n+dsp,1,nvec, 
     $                     g_r(ipm),1    ,n    ,1,nvec)      
        enddo ! end-loop-ipm
        call ga_sync()
        do ipm = 1,ncomp
          call ga_add(1.0d0, g_r(ipm), 
     &                1.0d0, g_x(ipm), 
     &                       g_x(ipm))
        enddo ! end-loop-ipm
        call ga_sync()
c        
        if (.not. ga_destroy(g_a)) call errquit
     &     ('lkain_2cpl: a',0, GA_ERR)
        if (.not. ga_destroy(g_bb)) call errquit
     &     ('lkain_2cpl: b',0, GA_ERR)
        if (.not. ga_destroy(g_c)) call errquit
     &     ('lkain_2cpl: c',0, GA_ERR)
c       
c       Reduce the subspace as necessary
c       
        if (nsub .eq. maxsub) then
          do isub = nvec+1, maxsub, nvec
c           component 1: 
            call ga_copy_patch('n', 
     $         g_Ay,    1, n, isub, isub+nvec-1, 
     $         g_Ax(1), 1, n, 1, nvec)
            call ga_sync()
            call ga_copy_patch('n', 
     $         g_Ax(1), 1, n, 1, nvec,
     $         g_Ay,    1, n, isub-nvec, isub-1)
            call ga_sync()
c           component 2: 
            if (ncomp.gt.1) then
              call ga_copy_patch('n', 
     $           g_Ay,    n+1, n2, isub, isub+nvec-1, 
     $           g_Ax(2), 1,   n,  1, nvec)
              call ga_sync()
              call ga_copy_patch('n', 
     $           g_Ax(2), 1,   n,  1, nvec,
     $           g_Ay,    n+1, n2, isub-nvec, isub-1)            
              call ga_sync()
            endif
c           
c           component 1:
            call ga_copy_patch('n', 
     $         g_y,     1, n, isub, isub+nvec-1, 
     $         g_Ax(1), 1, n, 1, nvec)
            call ga_sync()
            call ga_copy_patch('n', 
     $         g_Ax(1), 1, n, 1, nvec,
     $         g_y,     1, n, isub-nvec, isub-1)            
            call ga_sync()
c           component 2:
            if (ncomp.gt.1) then
              call ga_copy_patch('n', 
     $           g_y,     n+1, n2, isub, isub+nvec-1, 
     $           g_Ax(2), 1,   n,  1, nvec)
              call ga_sync()
              call ga_copy_patch('n', 
     $           g_Ax(2), 1,   n,  1, nvec,
     $           g_y,     n+1, n2, isub-nvec, isub-1)
              call ga_sync()
            endif
c           
          end do                ! isub = nvec+1, maxsub, nvec
          nsub = nsub - nvec
        end if                  ! (nsub .eq. maxsub) 
c       
      end do                    ! iter = 1,maxiter
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
c       
      enddo                     ! ipm = 1,2
      
      if (.not. ga_destroy(g_Ay)) call errquit
     &   ('lkain_2cpl: destroy Ay',3, GA_ERR)
      if (.not. ga_destroy(g_y)) call errquit
     &   ('lkain_2cpl: destroy r',4, GA_ERR)
      if (.not. ga_destroy(g_r2)) call errquit
     &   ('lkain_2cpl: destroy r2',6, GA_ERR)
      
      
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

      subroutine ga_lkain_2cpl3_damp(rtdb, 
     &                               g_x, g_b, 
     &                               g_x_im, g_b_im,
     &                               product, precond, 
     $                               tol, mmaxsub, maxiter, 
     &                               odiff, oprint, omega, limag,
     &                               lifetime, gamwidth, ncomp)
c  Written by J. Autschbach, SUNY Buffalo
c  Improvements made
c          by F. Aquino,     Northwestern University 
c          03-15-12
c  Note.- Modifying/Improving ga_lkain_2cpl3_damp() as of old src code
c         appeared on download src 11-20-12
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
      double precision omg(ncomp),
     &                 gam(ncomp)
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
      integer iter, n, n2, nvec, nsub, isub, type, maxsub, ipm,
     &   ntmp1, ntmp2, n3, n4

c ... jochen: for convenience, now most arrays have two components.
c     that might be changed later if memory becomes an issue
      integer g_y, g_Ay, g_Ax(2), g_r(2), g_r2, g_a, g_bb,
     &   g_c, g_xold(2), g_Axold(2)
      integer g_r_im(2), g_Ax_im(2)
      double precision rmax, rmax1, rmax2, acc
      logical converged
      logical odebug, debug, converge_precond
      double complex val_cmplx
      logical debug1
      integer p1,p2,m1,m2,stat_solve
      
c      integer l_zre,k_zre,g_z,
c     &        l_zim,k_zim
      external copy_rtor2,copy_r2tor,copy_AxxtoAyy,
     &         update_g_x1

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

      call ga_inquire(g_x(1), type, n, nvec)

c ------- create (zre,zim) ---------- START
c      if (.not.MA_Push_Get(mt_dbl,n,'hessv jfacs',l_zre,k_zre))
c     &     call errquit('sh-fockbld2: cannot allocate jfac',
c     &                  n, MA_ERR)
c      if (.not.MA_Push_Get(mt_dbl,n,'hessv kfacs',l_zim,k_zim))
c     &     call errquit('sh-fockbld2: cannot allocate kfac',
c     &                  n, MA_ERR)
c        if (.not. ga_create(MT_DCPL,n,nvec, 'lkain_2cpl: z',
c     $     0, 0, g_z))
c     $     call errquit('lkain: failed allocating z', nvec,
c     &     GA_ERR)
c           call ga_zero(g_z)
c ------- create (zre,zim) ---------- END
      if (ncomp.gt.1) then
        call ga_inquire(g_x(2), type, ntmp1, ntmp2)
        
c       ... jochen: do a sanity check on the array dimensions
        if (ntmp1.ne.n .or. ntmp2.ne.nvec) call errquit
     &     ('ga_lkain_2cpl:inconsistent dimensions of g_x components',
     &     nvec,CALC_ERR)
      endif
      
c     later we combine the two components to vecors of double
c     length and combine again Re and Im, i.e. 
c     the dimension is up to 4*n

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
c     
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
 112          format('------ prod-g_x-1-re(',i3,',',i3,')------ START')
             endif
             call ga_print(g_x(ipm))
             if (ga_nodeid().eq.0) then
              write(*,113) iter,ipm
 113          format('------ prod-g_x-1-re(',i3,',',i3,')------ END')
             endif
             if (ga_nodeid().eq.0) then
              write(*,114) iter,ipm
 114          format('------ prod-g_Ax-1-re(',i3,',',i3,')------ START')
             endif
             call ga_print(g_Ax(ipm))
             if (ga_nodeid().eq.0) then
              write(*,115) iter,ipm
 115          format('------ prod-g_Ax-1-re(',i3,',',i3,')------ END')
             endif
             if (lifetime) then

             if (ga_nodeid().eq.0) then
              write(*,312) iter,ipm
 312          format('------ prod-g_x-1-im(',i3,',',i3,')------ START')
             endif
             call ga_print(g_x_im(ipm))
             if (ga_nodeid().eq.0) then
              write(*,313) iter,ipm
 313          format('------ prod-g_x-1-im(',i3,',',i3,')------ END')
             endif
             if (ga_nodeid().eq.0) then
              write(*,314) iter,ipm
 314          format('------ prod-g_Ax-1-im(',i3,',',i3,')------ START')
             endif
             call ga_print(g_Ax_im(ipm))
             if (ga_nodeid().eq.0) then
              write(*,315) iter,ipm
 315          format('------ prod-g_Ax-1-im(',i3,',',i3,')------ END')
             endif

             endif ! end-if-lifetime
            enddo ! end-loop-ipm
         endif ! end-if-debug
        
        if (debug) write (6,*)
     &     'calling product from ga_lkain_2cpl_damp'
        call product(acc, g_x, g_Ax,  g_x_im, g_Ax_im, omega, limag,
     &     lifetime, gamwidth, ncomp)
         if (debug) then
            do ipm=1,ncomp

             if (ga_nodeid().eq.0) then
              write(*,116) iter,ipm
 116          format('------ prod-g_x-2-re(',i3,',',i3,')------ START')
             endif
             call ga_print(g_x(ipm))
             if (ga_nodeid().eq.0) then
              write(*,117) iter,ipm
 117          format('------ prod-g_x-2-re(',i3,',',i3,')------ END')
             endif
             if (ga_nodeid().eq.0) then
              write(*,118) iter,ipm
 118          format('------ prod-g_Ax-2-re(',i3,',',i3,')------ START')
             endif
             call ga_print(g_Ax(ipm))
             if (ga_nodeid().eq.0) then
              write(*,119) iter,ipm
 119          format('------ prod-g_Ax-2-re(',i3,',',i3,')------ END')
             endif
             if (ga_nodeid().eq.0) then
              write(*,316) iter,ipm
 316          format('------ prod-g_x-2-im(',i3,',',i3,')------ START')
             endif
             call ga_print(g_x_im(ipm))
             if (ga_nodeid().eq.0) then
              write(*,317) iter,ipm
 317          format('------ prod-g_x-2-im(',i3,',',i3,')------ END')
             endif
             if (ga_nodeid().eq.0) then
              write(*,318) iter,ipm
 318          format('------ prod-g_Ax-2-im(',i3,',',i3,')------ START')
             endif
             call ga_print(g_Ax_im(ipm))
             if (ga_nodeid().eq.0) then
              write(*,319) iter,ipm
 319          format('------ prod-g_Ax-2-im(',i3,',',i3,')------ END')
             endif
             enddo ! end-loop-ipm

           endif ! end-if-debug
c             if (iter.eq.3) then
c               if (ga_nodeid().eq.0)
c     &          write(*,*) 'FA-STOP-convergence at 3rd iter'
c               stop
c             endif              

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
     &                       g_r(ipm)) ! The residual, Real part

          if (debug) then
             if (ga_nodeid().eq.0) then
              write(*,120) iter,ipm
 120          format('------ prod-g_b-re(',i3,',',i3,')------ START')
             endif
             call ga_print(g_b(ipm))
             if (ga_nodeid().eq.0) then
              write(*,121) iter,ipm
 121          format('------ prod-g_b-re(',i3,',',i3,')------ END')
             endif
             if (ga_nodeid().eq.0) then
              write(*,122) iter,ipm
 122          format('------ prod-g_r-re(',i3,',',i3,')------ START')
             endif
             call ga_print(g_r(ipm))
             if (ga_nodeid().eq.0) then
              write(*,123) iter,ipm
 123          format('------ prod-g_r-re(',i3,',',i3,')------ END')
             endif
          endif ! end-if-debug

          if (lifetime) then
            call ga_add( 1.0d0,g_b_im(ipm),
     &                  -1.0d0,g_Ax_im(ipm), 
     &                         g_r_im(ipm)) ! The residual, Im part
          endif

          if (debug) then
             if (ga_nodeid().eq.0) then
              write(*,320) iter,ipm
 320          format('------ prod-g_b-im(',i3,',',i3,')------ START')
             endif
             call ga_print(g_b_im(ipm))
             if (ga_nodeid().eq.0) then
              write(*,321) iter,ipm
 321          format('------ prod-g_b-im(',i3,',',i3,')------ END')
             endif
             if (ga_nodeid().eq.0) then
              write(*,322) iter,ipm
 322          format('------ prod-g_r-im(',i3,',',i3,')------ START')
             endif
             call ga_print(g_r_im(ipm))
             if (ga_nodeid().eq.0) then
              write(*,323) iter,ipm
 323          format('------ prod-g_r-im(',i3,',',i3,')------ END')
             endif
          endif ! end-if-debug
c         
        enddo                   ! ipm = 1,ncomp

c       convergence checking:
c       find the largest element of the residual either 
c       before or after the call to the preconditioner
        omg(1)=-omega
        omg(2)= omega
        gam(1)=-gamwidth
        gam(2)= gamwidth
        if (converge_precond) then
c         do ipm=1,ncomp
c          call precond(g_r(ipm),g_r_im(ipm),-omega,gamwidth)
c         enddo ! end-loop-ipm
c ---- FA-FIX-12-01-13 --- START
c            call precond(g_r(1),  g_r_im(1),  -omega, gamwidth)
c          if (ncomp.gt.1) then
c            call precond(g_r(2),  g_r_im(2),   omega, gamwidth)
c          endif  
c ---- FA-FIX-12-01-13 --- END
c ---- FA-03-12-14 --- START
         do ipm=1,ncomp
          call precond(g_r(ipm),g_r_im(ipm),omg(ipm),gam(ipm))
         enddo
c ---- FA-03-12-14 --- END
        endif

        call ga_maxelt(g_r(1), rmax1)
        if (ncomp.gt.1) then
          call ga_maxelt(g_r(2), rmax2)
        else
          rmax2 = 0d0
        endif
        rmax = max(rmax1, rmax2)
        if (lifetime) then
          call ga_maxelt(g_r_im(1), rmax1)
          if (ncomp.gt.1) call ga_maxelt(g_r_im(2), rmax2)
          rmax = max(rmax, rmax1)
          rmax = max(rmax, rmax2)
        endif
c ---------- FA-12-03-13 ----- START
          if (debug1) then
           write(6,4) iter, nsub+nvec, rmax, util_wallsec()
           call util_flush(6)
    4      format('FA-chk: ', i5, i7, 3x,1p,d9.2,0p,f10.1,5x,i3)
          endif
c ---------- FA-12-03-13 ----- END         
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
c       
c       Call the preconditioner with the residual as well as with
c       the quantity Ax
c ... jochen: changed 0d0 to omega in the calls below compared to
c       the original routine ga_lkain
c       (there were only two in ga_lkain, one with g_aX, one with g_r)
c       for array g_r the preconditioner call is only necessary
c       in case converge_precond is .false.


c        call precond(g_Ax(1), g_Ax_im(1), -omega, gamwidth)
c        if (.not.converge_precond) 
c     &     call precond(g_r(1),  g_r_im(1),  -omega, gamwidth)
c        if (ncomp.gt.1) then
c          call precond(g_Ax(2), g_Ax_im(2),  omega, gamwidth)
c          if (.not.converge_precond) 
c     &       call precond(g_r(2),  g_r_im(2),   omega, gamwidth)
c        endif
c FA-03-12-14 ---- START  
        do ipm=1,ncomp
         call precond(g_Ax(ipm), g_Ax_im(ipm),omg(ipm), gam(ipm))
        enddo
        if (.not.converge_precond) then
         do ipm=1,ncomp
         call precond(g_r(ipm) , g_r_im(ipm),omg(ipm), gam(ipm))
         enddo
        endif
c FA-03-12-14 ---- END  

        if (debug) write (6,*) 'lkain3_damp: back from precond'
                
c       Copy the vectors to the subspace work area
c   
c ---- FA-copy ((Ax,Ax_im),(x,x_im)) -> (Ay,y) ------ START
c FA: Step 2:
        if (debug1) then
        do ipm = 1,ncomp  

         if (ga_nodeid().eq.0) then
          write(*,324) iter,ipm
 324      format('------ g_x-re(',i3,',',i3,')------ START')
         endif
         call ga_print(g_x(ipm))
         if (ga_nodeid().eq.0) then
          write(*,325) iter,ipm
 325      format('------ g_x-re(',i3,',',i3,')------ START')
         endif 
         if (ga_nodeid().eq.0) then
          write(*,326) iter,ipm
 326      format('------ g_x-im(',i3,',',i3,')------ START')
         endif
         call ga_print(g_x_im(ipm))
         if (ga_nodeid().eq.0) then
          write(*,327) iter,ipm
 327      format('------ g_x-im(',i3,',',i3,')------ START')
         endif 

         if (ga_nodeid().eq.0) then
          write(*,328) iter,ipm
 328      format('------ g_Ax-re(',i3,',',i3,')------ START')
         endif
         call ga_print(g_Ax(ipm))
         if (ga_nodeid().eq.0) then
          write(*,329) iter,ipm
 329      format('------ g_Ax-re(',i3,',',i3,')------ START')
         endif 
         if (ga_nodeid().eq.0) then
          write(*,330) iter,ipm
 330      format('------ g_Ax-im(',i3,',',i3,')------ START')
         endif
         call ga_print(g_Ax_im(ipm))
         if (ga_nodeid().eq.0) then
          write(*,331) iter,ipm
 331     format('------ g_Ax-im(',i3,',',i3,')------ START')
         endif 

         if (ga_nodeid().eq.0) then
          write(*,332) iter,ipm
 332      format('------ g_r-re(',i3,',',i3,')------ START')
         endif
         call ga_print(g_r(ipm))
         if (ga_nodeid().eq.0) then
          write(*,333) iter,ipm
 333      format('------ g_r-re(',i3,',',i3,')------ START')
         endif 
         if (ga_nodeid().eq.0) then
          write(*,334) iter,ipm
 334      format('------ g_r-im(',i3,',',i3,')------ START')
         endif
         call ga_print(g_r_im(ipm))
         if (ga_nodeid().eq.0) then
          write(*,335) iter,ipm
 335     format('------ g_r-im(',i3,',',i3,')------ START')
         endif 

        enddo ! end-loop-ipm
        endif ! end-if-debug1  
     
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

        nsub = nsub + nvec
c       
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

c ++++++++++++++ real linear solver ++++++++++++++++++ 
        call ga_svd_solve_seq(g_a,g_bb,g_c,1d-14)
c ++++++++++++++++++++++++++++++++++++++++++++++++++++

        if (debug1) then
        if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_c(',iter,')-----START' 
        call ga_print(g_c)
        if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_c(',iter,')-----END'
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
        if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_x-AFT-1(',iter,')-----START' 
        call ga_print(g_x)
        if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_x-AFT-1(',iter,')-----END'
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
        if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_x-AFT-2(',iter,')-----START' 
        call ga_print(g_x)
        if (ga_nodeid().eq.0)
     &   write(*,*) '---------g_x-AFT-2(',iter,')-----END'  
        endif ! end-if-debug1
       
        if (.not. ga_destroy(g_a)) call errquit
     &     ('lkain_2cpl: a',0, GA_ERR)
        if (.not. ga_destroy(g_bb)) call errquit
     &     ('lkain_2cpl: b',0, GA_ERR)
        if (.not. ga_destroy(g_c)) call errquit
     &     ('lkain_2cpl: c',0, GA_ERR)
c       
c       Reduce the subspace as necessary
c       
        if (ga_nodeid().eq.0) then
         write(*,17) iter,nsub,maxsub,nvec,n4
 17      format('(iter,nsub,maxsub,nvec,n4)=(',
     &          i3,',',i4,',',i4,',',i3,',',i4,')')
        endif
        if (debug1) then
         if (iter.eq.6) then
          write(*,*) 'STOP-test-linsolver'
          stop
         endif
        endif ! end-if-debug1 

        if (nsub .eq. maxsub) then
c ====== FA: left-shifting patch ==== START
c Note.- matrices Ay,y shift to left nvec positions
c        removing leftmost patch of dimension: n4 x nvec
         if (ga_nodeid().eq.0)
     &    write(*,*) 'FA-matrix-nvec-left-shifting:'
         do isub = nvec+1, maxsub, nvec
          call ga_copy_patch('n',g_Ay,1,n4,isub,isub+nvec-1, 
     $                           g_Ay,1,n4,isub-nvec,isub-1)
          call ga_copy_patch('n',g_y ,1,n4,isub,isub+nvec-1, 
     $                           g_y ,1,n4,isub-nvec,isub-1)
c ====== FA: left-shifting patch ==== END
         enddo ! end-loop-isub
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

#if 1
c     error exit if this hasn't converged :-(
      if (.not. converged) call errquit('lkain_2cpl: not converged',0,
     &    CALC_ERR)
#else
      if (.not. converged) then
        if (ga_nodeid().eq.0) then
          write (luout,*) 'WARNING: CPKS procedure is NOT converged'
          write (luout,*) '  I will proceed, but check your results '//
     &                       'carefully!!!'
        endif
      endif
#endif
c     
      end
