c +++++++++++   ga_lkain_ext() ++++++++++++ START
c Note.- Useful for AOresponse calculations (FA-12-10-12)
      subroutine ga_lkain_ext(
     &               rtdb,
     &               g_x, 
     &               g_b, 
     &               product, 
     &               precond, 
     $               tol, 
     &               mmaxsub, 
     &               maxiter, 
     &               odiff, 
     &               oprint)
c
c Note.- Modified from original aoresponse source code
c        for extension to spin-unrestricted case
c        original aoresponse source code was written by 
c        J. Autschbach et.al. and appears on nwchem-devtrunk (date:03-02-12)
c Author of changes: Fredy W. Aquino, Northwestern University
c Date             : 03-15-12

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "rtdb.fh"

      integer g_x               ! [input/output] Initial guess/solution
      integer g_b               ! [input] Right-hand side vectors
      external product          ! [input] product routine
      external precond          ! [input] preconditioner routine
      double precision tol      ! [input] convergence threshold
      integer mmaxsub           ! [input] maximum subspace dimension
      integer maxiter           ! [input] maximum no. of iterations
      logical odiff             ! [input] use differences in product
      logical oprint            ! [input] print flag
      integer rtdb
      logical debug
      external ga_svd_solve_seq
c
c     Solves the linear equations A(X)=0 for multiple vectors.
c
c     call product(acc,g_x, g_Ax)
c     . acc is the accuracy trequired for each element of the product
c     . g_x contains the vectors and g_Ax should be filled
c     .     with the product vectors.  The no. of vectors (columns) in
c     . g_x might differ from the no. of vectors input to ga_lkain().
c
c     call precond(g_x,shift)
c     . apply preconditioning directly to the vectors in g_x with the
c     . optional shift (not used here but used by the diagonalizer)
c
c     On input g_x should contain an initial guess.  It returns the
c     solution.
c
c     maxsub should be at least 3*nvec and can be beneficially increased
c     to about 10*nvec.
c
c     Needs to be extended to store the sub-space vectors out-of-core
c     at least while the product() routine is being executed.

      integer iter, n, nvec, nsub, isub, type, maxsub
      integer g_y, g_Ay, g_Ax, g_r, g_a, g_bb, g_c, g_xold, g_Axold
      double precision factor_x
      double precision rmax, acc, ga_svd_tol
      logical converged
      logical odebug
      debug=.false. ! FA-added4debugging

c      debug=.true.

      odebug = util_print('debug lsolve', print_never) .and. 
     $     ga_nodeid().eq.0
cold      if (.not.rtdb_get(rtdb,'cphf:acc',mt_dbl,1,acc)) acc=0.01d0*tol
      if (.not.rtdb_get(rtdb,'cphf:acc',mt_dbl,1,acc)) acc=1d-4*tol
      call ga_inquire(g_x, type, n, nvec)
      maxsub = mmaxsub          ! So don't modify input scalar arg
      if (maxsub .lt. 3*nvec) maxsub = 3*nvec
      maxsub = (maxsub/nvec)*nvec

      if (.not.rtdb_get(rtdb,'cphf:ga_svd_tol',mt_dbl,1,
     &                  ga_svd_tol)) then
c       See comment just before the ga_svd_solve_seq call to
c       understand these choices.
        if ((100*maxsub).lt.n) then
          ga_svd_tol = 1d-7
        else
          ga_svd_tol = 1d-14
        endif 
      endif

      if (oprint .and. ga_nodeid().eq.0) then
         write(6,1) n, nvec, maxsub, maxiter, tol, util_wallsec()
 1       format(//,'Iterative solution of linear equations',/,
     $        '  No. of variables', i9,/,
     $        '  No. of equations', i9,/,
     $        '  Maximum subspace', i9,/,
     $        '        Iterations', i9,/,
     $        '       Convergence', 1p,d9.1,/,
     $        '        Start time', 0p,f9.1,/)
         call util_flush(6)
      end if

      if (.not. ga_create(MT_DBL, n, maxsub, 'lkain: Y', 
     $     0, 0, g_y))
     $     call errquit('lkain: failed allocating subspace', maxsub,
     &       GA_ERR)
      if (.not. ga_create(MT_DBL, n, maxsub, 'lkain: Ay', 
     $     0, 0, g_Ay))
     $     call errquit('lkain: failed allocating subspace2', maxsub,
     &       GA_ERR)
      if (.not. ga_create(MT_DBL, n, nvec, 'lkain: Ax',
     $     0, 0, g_Ax))
     $     call errquit('lkain: failed allocating subspace3', nvec,
     &       GA_ERR)
      if (.not. ga_create(MT_DBL, n, nvec, 'lkain: r',
     $     0, 0, g_r))
     $     call errquit('lkain: failed allocating subspace4', nvec,
     &       GA_ERR)
      if (odiff) then
         if (.not. ga_create(MT_DBL, n, nvec, 'lkain: xold',
     $        0, 0, g_xold))
     $        call errquit('lkain: failed allocating subspace5', nvec,
     &       GA_ERR)
         if (.not. ga_create(MT_DBL, n, nvec, 'lkain: xold',
     $        0, 0, g_Axold))
     $        call errquit('lkain: failed allocating subspace6', nvec,
     &       GA_ERR)
         call ga_zero(g_xold)
         call ga_zero(g_Axold)
      end if
      call ga_zero(g_y)
      call ga_zero(g_Ay)
      call ga_zero(g_Ax)
      call ga_zero(g_r)

      if (oprint .and. ga_nodeid().eq.0) then
         write(6,2)
         call util_flush(6)
 2       format(/
     $        '   iter   nsub   residual    time',/,
     $        '   ----  ------  --------  ---------')
      end if
      nsub = 0
      converged = .false.
          if (debug) then  
           if (ga_nodeid().eq.0)
     &      write(*,*) '---- g_x guess-------- START'
            call ga_print(g_x)
           if (ga_nodeid().eq.0)
     &      write(*,*) '---- g_x guess--------  END'
          endif ! end-if-debug

      do iter = 1, maxiter
         if (odiff) then
            call ga_add(1.0d0, g_x, -1.0d0, g_xold,  g_x)
         end if
c ======== visualize input to product ============ START
          if (debug) then
             if (ga_nodeid().eq.0) then
              write(*,112) iter
 112          format('------ prod-g_x-1(',i3,')------ START')
             endif
             call ga_print(g_x)
             if (ga_nodeid().eq.0) then
              write(*,113) iter
 113          format('------ prod-g_x-1(',i3,')------ END')
             endif
             if (ga_nodeid().eq.0) then
              write(*,114) iter
 114          format('------ prod-g_Ax-1(',i3,')------ START')
             endif
             call ga_print(g_Ax)
             if (ga_nodeid().eq.0) then
              write(*,115) iter
 115          format('------ prod-g_Ax-1(',i3,')------ END')
             endif
          endif ! end-if-debug
c ======== visualize input to product ============ END

         call product(acc,g_x, g_Ax)
          
c ======== visualize output to product ============ START
          if (debug) then
             if (ga_nodeid().eq.0) then
              write(*,116) iter
 116          format('------ prod-g_x-2(',i3,')------ START')
             endif
             call ga_print(g_x)
             if (ga_nodeid().eq.0) then
              write(*,117) iter
 117          format('------ prod-g_x-2(',i3,')------ END')
             endif
             if (ga_nodeid().eq.0) then
              write(*,118) iter
 118          format('------ prod-g_Ax-2(',i3,')------ START')
             endif
             call ga_print(g_Ax)
             if (ga_nodeid().eq.0) then
              write(*,119) iter
 119          format('------ prod-g_Ax-2(',i3,')------ END')
             endif
          endif ! end-if-debug
c ======== visualize output to product ============ END

         if (odiff) then
            call ga_add(1.0d0, g_Ax, 1.0d0, g_Axold, g_Ax)
            call ga_add(1.0d0, g_x,  1.0d0, g_xold,  g_x)
            call ga_copy(g_x, g_xold)
            call ga_copy(g_Ax, g_Axold)
         end if
          call ga_zero(g_r)
          call ga_add_patch( 1.0d0,g_b ,1,n,1,nvec, 
     $                      -1.0d0,g_Ax,1,n,1,nvec,
     $                             g_r ,1,n,1,nvec) ! The residual

         if (debug) then  
             if (ga_nodeid().eq.0) then
              write(*,120) iter
 120          format('------ prod-g_b(',i3,')------ START')
             endif
             call ga_print(g_b)
             if (ga_nodeid().eq.0) then
              write(*,121) iter
 121          format('------ prod-g_b(',i3,')------ END')
             endif
             if (ga_nodeid().eq.0) then
              write(*,122) iter
 122          format('------ prod-g_r(',i3,')------ START')
             endif
             call ga_print(g_r)
             if (ga_nodeid().eq.0) then
              write(*,123) iter
 123          format('------ prod-g_r(',i3,')------ END')
             endif
         endif ! end-if-debug

         call ga_maxelt(g_r, rmax)
         if (oprint .and. ga_nodeid().eq.0) then
            write(6,3) iter, nsub+nvec, rmax, util_wallsec()
            call util_flush(6)
 3          format(' ', i5, i7, 3x,1p,d9.2,0p,f10.1)
         end if
         if (rmax .lt. tol) then
            converged = .true.
           if (debug) then
             if (ga_nodeid().eq.0)
     &       write(*,*) 'Convergence achieved !'           
             if (ga_nodeid().eq.0)
     &       write(*,*) '---- g_xsol-------- START'
             call ga_print(g_x)
             if (ga_nodeid().eq.0)
     &       write(*,*) '---- g_xsol--------  END'
           endif ! end-if-debug
            goto 100
         end if
         call precond(g_Ax,0.0d0)
         call precond(g_r,0.0d0)
c
c     Copy the vectors to the subspace work area

         call ga_copy_patch('n', 
     $        g_Ax, 1, n, 1, nvec, 
     $        g_Ay, 1, n, nsub+1, nsub+nvec)
         call ga_copy_patch('n', 
     $        g_x, 1, n, 1, nvec, 
     $        g_y, 1, n, nsub+1, nsub+nvec)
         nsub = nsub + nvec
c
c     Form and solve the subspace equations using SVD in order
c     to manage near linear dependence in the subspace.
c     
         if (.not. ga_create(MT_DBL, nsub, nsub, 'lkain: A', 0, 0, g_a))
     $        call errquit('lkain: allocating g_a?', nsub, GA_ERR)
         if (.not. ga_create(MT_DBL, nsub, nvec, 'lkain: B', 0, 0,g_bb))
     $        call errquit('lkain: allocating g_bb?', nsub, GA_ERR)
         if (.not. ga_create(MT_DBL, nsub, nvec, 'lkain: C', 0, 0, g_c))
     $        call errquit('lkain: allocating g_c?', nsub, GA_ERR)
         call ga_zero(g_a)
         call ga_zero(g_bb)
         call ga_zero(g_c)
         call ga_dgemm('t','n',nsub,nsub,n,1.0d0,g_y,g_Ay,0.0d0,g_a)
         call ga_dgemm('t','n',nsub,nvec,n,1.0d0,g_y,g_r ,0.0d0,g_bb)
         if (odebug) call ga_print(g_a)
         if (odebug) call ga_print(g_c)
c
c     The threshold used here should reflect the accuracy in the
c     products.  If very accurate products are used, then there is big
c     advantage for small cases (maxsub close to n) in using a very
c     small threshold in the SVD solve (e.g., 1e-14), but for more
c     realistic examcall ga_sync() ! FA-added 08-13-12ples (maxsub << n) there is only a little
c     advantage and in the precence of real noise in the products
c     screening with a realistic threshold is important.

         call ga_svd_solve_seq(g_a,g_bb,g_c,ga_svd_tol)

         if (odebug) call ga_print(g_c)
c
c     Form and add the correction, in parts, onto the solution
         call ga_dgemm('n','n',n,nvec,nsub,-1.0d0,g_Ay,g_c,1.0d0,g_r)
         if (odebug) then
            write(6,*) ' The update in the complement '
            call ga_print(g_r)
         end if
         call ga_add(1.0d0, g_r, 1.0d0, g_x, g_x)
         call ga_dgemm('n','n',n,nvec,nsub,1.0d0,g_y,g_c,0.0d0,g_r)
         if (odebug) then
            write(6,*) ' The update in the subspace '
            call ga_print(g_r)
         end if
         call ga_add(1.0d0, g_r, 1.0d0, g_x, g_x)
         if (.not. ga_destroy(g_a)) call errquit('lkain: a',0, GA_ERR)
         if (.not. ga_destroy(g_bb))call errquit('lkain: b',0, GA_ERR)
         if (.not. ga_destroy(g_c)) call errquit('lkain: c',0, GA_ERR)

c     Reduce the subspace as necessary
         if (nsub .eq. maxsub) then
          do isub = nvec+1, maxsub, nvec
            call ga_copy_patch('n',g_Ay,1,n,isub,isub+nvec-1, 
     $                             g_Ay,1,n,isub-nvec,isub-1)
            call ga_copy_patch('n',g_y ,1,n,isub,isub+nvec-1, 
     $                             g_y ,1,n,isub-nvec,isub-1)
          end do
          nsub = nsub - nvec
         end if
      enddo                    ! iter = 1,maxiter
 100  continue

      if (odiff) then
         if (.not. ga_destroy(g_xold)) call errquit('lkain: destroy',1,
     &       GA_ERR)
         if (.not. ga_destroy(g_Axold)) call errquit('lkain: destroy',2,
     &       GA_ERR)
      end if
      if (.not. ga_destroy(g_Ax)) call errquit('lkain: destroy',20,
     &       GA_ERR)
      if (.not. ga_destroy(g_Ay)) call errquit('lkain: destroy',3,
     &       GA_ERR)
      if (.not. ga_destroy(g_y)) call errquit('lkain: destroy',4,
     &       GA_ERR)
      if (.not. ga_destroy(g_r)) call errquit('lkain: destroy',5,
     &       GA_ERR)

      if (.not. converged) call errquit('lkain: not converged',0,
     &       CALC_ERR)

      end
c +++++++++++   ga_lkain_ext() ++++++++++++ END
c $Id$
