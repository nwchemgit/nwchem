*
* $Id$
*
c        Modified driver for use in NWChem (MV)
c
c        
c        Wed Jul 27 11:18:47 PDT 2005
c        *****************************************
c        L-BFGS-B is a code for solving large nonlinear optimization
c             problems with simple bounds on the variables.
c
c        The code can also be used for unconstrained problems and is
c        as efficient for these problems as the earlier limited memory
c                          code L-BFGS.
c
c        This is the simplest driver in the package. It uses all the
c                    default settings of the code.
c
c
c     References:
c
c        [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
c        memory algorithm for bound constrained optimization'',
c        SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
c
c        [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
c        Subroutines for Large Scale Bound Constrained Optimization''
c        Tech. Report, NAM-11, EECS Department, Northwestern University,
c        1994.
c
c                              *  *  *
c
c        NEOS, November 1994. (Latest revision June 1996.)
c        Optimization Technology Center.
c        Argonne National Laboratory and Northwestern University.
c        Written by
c                           Ciyou Zhu
c        in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c     **************
      logical function lbfgsb_driver(irtdb)
      implicit none
#include "geom.fh" 
#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "global.fh"
#include "sym.fh" 
 
      integer          irtdb
 
c     Declare the variables needed by the code.
c       A description of all these variables is given at the end of 
c       the driver.
 
      character*60     task, csave
      character*60     pname
      character*255     xyz
      character*32 theory
      logical          lsave(4)
      integer          m, iprint
      integer isave(44)
      double precision  dsave(29)
      double precision f, factr, pgtol,trust
      integer nwa
      integer i_nbd,h_nbd
      integer i_iwa,h_iwa
      integer i_x,h_x,i_xx
      integer i_index_unique,h_index_unique
      integer i_x0,h_x0
      integer i_l,h_l
      integer i_u,h_u
      integer i_g,h_g,i_gg
      integer i_wa,h_wa

      logical ostress,ostress2,onlystress,usesymmetry
      logical status
      logical task_gradient
      external task_gradient
      logical master
      external lbfgs_converged
      logical  lbfgs_converged
      logical geom_amatrix_set
      external geom_amatrix_set
      logical  geom_grad_cart_to_frac
      external geom_grad_cart_to_frac

c

      integer          i,j
c
      integer igeom,isystype
      integer natom,ncent,nc,ncx,ncent_unique
      logical rtdb_mode

      double precision xrms,xmax
      double precision grms,gmax
      double precision f0


      double precision tol(4)
      integer maxiter
      integer nfmax
      double precision matrix(3,4)

      pname = "lbfgsb_driver"

      master = ga_nodeid().eq.0


      if (.not. rtdb_get(irtdb,'onlystress',mt_log,1,onlystress)) 
     $   onlystress = .false.

      if (.not. rtdb_get(irtdb,'includestress',mt_log,1,ostress)) 
     $   ostress = .false.

      if (.not. rtdb_get(irtdb,'includelattice',mt_log,1,ostress2)) then
         ostress2 = .false.
      end if
      if ((.not.ostress).and.(ostress2)) ostress2 = .false.
      if ((ostress)     .and.(ostress2)) ostress  = .false.

       
      if ((.not.ostress).and.(.not.ostress2)) onlystress = .false.


      if (.not. rtdb_get(irtdb,'driver:gmax_tol',mt_dbl,1,tol(1)))
     $     tol(1) = 0.00045d0
      if (.not. rtdb_get(irtdb,'driver:grms_tol',mt_dbl,1,tol(2)))
     $     tol(2) = 0.0003d0
      if (.not. rtdb_get(irtdb,'driver:xmax_tol',mt_dbl,1,tol(3)))
     $     tol(3) = 0.0018d0
      if (.not. rtdb_get(irtdb,'driver:xrms_tol',mt_dbl,1,tol(4)))
     $     tol(4) = 0.0012d0
c
      if (.not. rtdb_get(irtdb,'driver:trust',mt_dbl,1,trust)) 
     $   trust = 0.5d0
C      trust = 1.0d3

c
      if (.not. rtdb_get(irtdb,'driver:nptopt',mt_int,1,maxiter))
     $     maxiter=20
c
      if (.not. rtdb_cget(irtdb,'driver:xyz',1,xyz))
     $     xyz = ' '
c
      if(.not.geom_create(igeom,'geometry'))
     + call errquit(pname//'Failed to create geometry',0, GEOM_ERR)

      if(.not.geom_rtdb_load(irtdb,igeom,'geometry'))
     + call errquit(pname//'Failed to load geometry',0, GEOM_ERR)

      if(.not. geom_ncent(igeom,ncent) ) 
     >    call errquit(pname//'l-bfgs:geom_ncent',0,GEOM_ERR)


      if(.not. geom_ncent_unique(igeom,ncent_unique) ) 
     >    call errquit(pname//'l-bfgs:geom_ncent_unique',0,GEOM_ERR)

      if (.not. geom_systype_get(igeom, isystype))
     $     call errquit(pname//'systype?',0, GEOM_ERR)

      usesymmetry = (sym_number_ops(igeom).gt.0)
      natom = ncent
      ncx   = 3*ncent

      if (ostress) then
         if (master) write(6,*) 'INCLUDING STRESS !!!!!!!!!!!!!!!!'
         if (master) write(6,*) 
         if (isystype.eq.0) call errquit('NOT A PERIODIC SYSTEM',0,
     >       GEOM_ERR)
         ncx = ncx + 9
         if (onlystress) then
            ncent = 3
         else
            ncent = ncent + 3
         end if
      endif

      if (ostress2) then
         if (master) write(6,*) 'INCLUDING LATTICE GRADIENTS !!!!!!!!'
         if (master) write(6,*)
         if (isystype.eq.0) call errquit('NOT A PERIODIC SYSTEM',0,
     >       GEOM_ERR)
         ncx = ncx + 6
         if (onlystress) then
            ncent = 2
         else
            ncent = ncent + 2
         end if
      endif

      nc    = 3*ncent


      if (.not. rtdb_get(irtdb,'driver:mh',mt_int,1,m))
     $        m = 5

      if (.not. rtdb_get(irtdb,'driver:nfmax',mt_int,1,nfmax))
     $        nfmax = 20

      nwa =  2*m*nc+4*nc+12*m*m+12*m


      if(.not.ma_alloc_get(mt_int,ncent_unique,'i_index_unique',
     >                     h_index_unique,i_index_unique))
     > call errquit(
     >      pname//'Failed to allocate memory for i_index_unique',
     >           ncent_unique, MA_ERR)

      if (.not.geom_uniquecent_get(igeom,ncent_unique,
     >                              int_mb(i_index_unique)))
     > call errquit(pname//'Failed to get unique centers',0,GEOM_ERR)


      if(.not.ma_alloc_get(mt_int,nc,'i_nbd',h_nbd,i_nbd))
     + call errquit( pname//'Failed to allocate memory for nbd',
     + nc, MA_ERR)

      if(.not.ma_alloc_get(mt_int,3*nc,'i_iwa',h_iwa,i_iwa))
     + call errquit( pname//'Failed to allocate memory for iwa',
     + nc, MA_ERR)

      if(.not.ma_alloc_get(mt_dbl,ncx,'i_x',h_x,i_x))
     + call errquit( pname//'Failed to allocate memory for x',
     + nc, MA_ERR)

      if (onlystress) then
         i_xx = i_x + 3*natom
      else
         i_xx = i_x
      end if

      if(.not.ma_alloc_get(mt_dbl,nc,'i_x0',h_x0,i_x0))
     + call errquit( pname//'Failed to allocate memory for x0',
     + nc, MA_ERR)

      if(.not.ma_alloc_get(mt_dbl,nc,'i_l',h_l,i_l))
     + call errquit( pname//'Failed to allocate memory for l',
     + nc, MA_ERR)

      if(.not.ma_alloc_get(mt_dbl,nc,'i_u',h_u,i_u))
     + call errquit( pname//'Failed to allocate memory for u',
     + nc, MA_ERR)

      if(.not.ma_alloc_get(mt_dbl,ncx,'i_g',h_g,i_g))
     + call errquit( pname//'Failed to allocate memory for g',
     + nc, MA_ERR)

      if (onlystress) then
         i_gg = i_g + 3*natom
      else 
         i_gg = i_g 
      end if

      if(.not.ma_alloc_get(mt_dbl,nwa,'i_wa',h_wa,i_wa))
     + call errquit( pname//'Failed to allocate memory for wa',
     + nc, MA_ERR)


c     debuging symmetry stuff 
      if (master) then
c         call sym_print_all(igeom,.true.,.true.,.true.,.false.,.false.)
         if (.not.geom_print(igeom))
     >    call errquit(pname//'Failed to print geommetry',0,GEOM_ERR)
      end if

c      write(*,*) "number of unique atoms=",ncent_unique

c      call sym_group_name(igeom,task)
c      write(*,*) "sym_number_ops=",sym_number_ops(igeom),task
c      do i = 1,natom
c         write(*,*) i," atom =",sym_center_map(igeom,i,1)
c      end do
c      call sym_ops_get(igeom,1,matrix)
c      write(*,*)
c      do i=1,3
c        write(*,*) (matrix(i,j),j=1,4)
c      end do
c      write(*,*)
c
c      call sym_print_ops(igeom)

c     debuging symmetry stuff 

 
c     We wish to have output at every iteration.

      iprint = 1

c     We specify the tolerances in the stopping criteria.

c      factr=1.0d+7
c      pgtol=1.0d-5

      factr=0.0d0
      pgtol=0.0d0

c     We specify the dimension n of the sample problem and the number
c        m of limited memory corrections stored.  (n and m should not
c        exceed the limits nmax and mmax respectively.)
 
 
c     We now provide nbd which defines the bounds on the variables:
c                    l   specifies the lower bounds,
c                    u   specifies the upper bounds. 
 
c     Set bounds if any
      do i=1,nc
         int_mb(i_nbd+i-1)=0
c         dbl_mb(i_l+i-1)=dbl_mb(i_xx+i-1) - 0.005
c         dbl_mb(i_u+i-1)=dbl_mb(i_xx+i-1) + 0.005
      end do

c     We now define the starting point.
      if(.not.geom_cart_coords_get(igeom,dbl_mb(i_x)))
     >   call errquit(pname//'Failed to get coord ',0, GEOM_ERR)

       if (ostress) then
          if (.not. geom_amatrix_get(igeom, dbl_mb(i_x+3*natom)))
     >        call errquit(pname//'failed to get amatrix',0,0)

          if (.not. geom_cart_to_frac(igeom, dbl_mb(i_x)))
     >        call errquit(pname//'frac_to_cart?',0, GEOM_ERR)
       end if

       if (ostress2) then
          if (.not. geom_lattice_get(igeom, dbl_mb(i_x+3*natom)))
     >        call errquit(pname//'failed to get lattice',0,0)

          if (.not. geom_cart_to_frac(igeom, dbl_mb(i_x)))
     >        call errquit(pname//'frac_to_cart?',0, GEOM_ERR)
       end if
c 
       call dcopy(nc,dbl_mb(i_xx),1,dbl_mb(i_x0),1)
       f0 = 0.0d0

c     We start the iteration by initializing task.
c 
      task = 'START'

c        ------- the beginning of the loop ----------
 
      lbfgsb_driver = .false.
 111  continue
      
c     This is the call to the L-BFGS-B code.
 

      if(master) then
        call setulb(nc,m,nfmax,
     +              dbl_mb(i_xx),
     +              dbl_mb(i_l),
     +              dbl_mb(i_u),
     +              int_mb(i_nbd),
     +              f,
     +              dbl_mb(i_gg),
     +              factr,pgtol,
     +              dbl_mb(i_wa),
     +              int_mb(i_iwa),
     +              task,iprint,
     +              csave,lsave,isave,dsave)

        rtdb_mode = rtdb_parallel(.false.)
        if (.not.rtdb_cput(irtdb,"lbfgsb:task",1,task)) 
     >          call errquit(pname//'failed RTDB',0,
     >           RTDB_ERR)
c
        if (task(1:2) .eq. 'FG') then
          if (ostress) then
             if (.not. geom_amatrix_set(igeom, dbl_mb(i_x+3*natom)))
     $        call errquit(pname//'failed to set amatrix',0,0)
          
             if (.not. geom_frac_to_cart(igeom, dbl_mb(i_x)))
     $         call errquit(pname//'frac_to_cart?',0,0)
         
          end if

          if (ostress2) then
             if (.not. geom_lattice_set(igeom, dbl_mb(i_x+3*natom)))
     $        call errquit(pname//'failed to set lattice',0,0)
          
             if (.not. geom_frac_to_cart(igeom, dbl_mb(i_x)))
     $         call errquit(pname//'frac_to_cart?',0,0)
          end if
         


          if(.not.geom_cart_coords_set(igeom,dbl_mb(i_x)))
     +       call errquit(pname//'Failed to set coord ',0, GEOM_ERR)

          if (usesymmetry) then
             call sym_geom_project(igeom,trust)
             if(.not.geom_cart_coords_get(igeom,dbl_mb(i_x)))
     >       call errquit(pname//'Failed to get coord ',0,GEOM_ERR)
          end if



          if(.not.geom_rtdb_store(irtdb,igeom,'geometry'))
     +     call errquit(pname//'Failed to store geometry',0, GEOM_ERR)

          if (ostress.or.ostress2) then
             if (.not. geom_cart_to_frac(igeom, dbl_mb(i_x)))
     $           call errquit(pname//'frac_to_cart?',0, GEOM_ERR)
          end if
        end if
      end if

      call ga_sync()

      rtdb_mode = rtdb_parallel(.true.)
      if (.not.rtdb_cget(irtdb,"lbfgsb:task",1,task)) 
     >        call errquit(pname//'failed rtdb_cget',0,
     >         RTDB_ERR)

c        the minimization routine has returned to request the
c        function f and gradient g values at the current x.
      if (task(1:2) .eq. 'FG') then

         if (.not. task_gradient(irtdb))
     $      call errquit(pname//'task_gradient failed',0, GEOM_ERR)

         if (.not.rtdb_get(irtdb,'task:energy', mt_dbl, 1, f))
     $      call errquit(pname//'could not get energy',0, RTDB_ERR)

c        **** don't need to get gradient if only stress ****
         if (.not.rtdb_get(irtdb, 'task:gradient', mt_dbl, 3*natom,
     $                     dbl_mb(i_g))) 
     $        call errquit(pname//'could not get gradient',0,0)

         if (usesymmetry) call sym_grad_symmetrize(igeom,dbl_mb(i_g))


         if (ostress) then
            if (.not. geom_grad_cart_to_frac(igeom, dbl_mb(i_g)))
     $         call errquit(pname//'frac_to_cart?',0,0)

            if (.not. rtdb_cget(irtdb, 'task:theory', 1, theory))
     $       call errquit(pname//'stress theory not specified',0,
     $                    RTDB_ERR)

            if (theory.eq.'pspw') then
             if (.not. rtdb_get(irtdb, 'pspw:stress', mt_dbl, 9,
     $           dbl_mb(i_g+3*natom))) call errquit
     $           (pname//'could not get pspw stress',0,0)

            else if (theory.eq.'band') then
             if (.not. rtdb_get(irtdb, 'band:stress', mt_dbl, 9,
     $           dbl_mb(i_g+3*natom))) call errquit
     $           (pname//'could not get band stress',0,0)

            else if (theory.eq.'paw') then
             if (.not. rtdb_get(irtdb, 'paw:stress', mt_dbl, 9,
     $           dbl_mb(i_g+3*natom))) call errquit
     $           (pname//'could not get paw stress',0,0)
            else
              call errquit(pname//'no stress in theory'//theory,
     $                     0,RTDB_ERR)
            end if
         end if


         if (ostress2) then
            if (.not. geom_grad_cart_to_frac(igeom, dbl_mb(i_g)))
     $         call errquit(pname//'frac_to_cart?',0,0)

            if (.not. rtdb_cget(irtdb, 'task:theory', 1, theory))
     $       call errquit(pname//'stress theory not specified',0,
     $                    RTDB_ERR)

            if (theory.eq.'pspw') then
             if (.not. rtdb_get(irtdb, 'pspw:lstress', mt_dbl, 6,
     $           dbl_mb(i_g+3*natom))) call errquit
     $           (pname//'could not get pspw lstress',0,0)

            else if (theory.eq.'band') then
             if (.not. rtdb_get(irtdb, 'band:lstress', mt_dbl, 6,
     $           dbl_mb(i_g+3*natom))) call errquit
     $           (pname//'could not get band lstress',0,0)

            else if (theory.eq.'paw') then
             if (.not. rtdb_get(irtdb, 'paw:lstress', mt_dbl, 6,
     $           dbl_mb(i_g+3*natom))) call errquit
     $           (pname//'could not get paw lstress',0,0)
            else
              call errquit(pname//'no lstress in theory'//theory,
     $                     0,RTDB_ERR)
            end if
         end if


c        go back to the minimization routine.
         goto 111
      endif

c     the minimization routine has returned with a new iterate,
c     so we have to save it into rtdb
      if (task(1:5) .eq. 'NEW_X')  then
         rtdb_mode = rtdb_parallel(.false.)
         if(master) then
            call lbfgs_compute_info(nc,
     +                                dbl_mb(i_xx),
     +                                dbl_mb(i_x0),
     +                                dbl_mb(i_gg),
     +                                xmax,xrms,gmax,grms)

            call lbfgs_print_info(isave(30)+1,
     +                               f,f-f0,xmax,xrms,
     +                               gmax,grms,
     +                               tol)

            call lbfgs_print_xyz(isave(30)+1,
     +                               igeom,
     +                               xyz)


            if(lbfgs_converged(xmax,xrms,gmax,grms,tol)) then
               task='CONVERGED: MET CONVERGENCE CRITERIA'
               if (.not.rtdb_cput(irtdb,"lbfgsb:task",1,task)) 
     >            call errquit(pname//'failed RTDB',0,RTDB_ERR)

            else if(isave(30)+1.gt.maxiter) then
               task='STOP: REACHED MAXIMUM NUMBER OF ITERATIONS'
               if (.not.rtdb_cput(irtdb,"lbfgsb:task",1,task))
     >            call errquit(pname//'failed RTDB',0,RTDB_ERR)
            end if

         end if !*master*

         rtdb_mode = rtdb_parallel(.true.)
         call ga_sync()
         if (.not.rtdb_cget(irtdb,"lbfgsb:task",1,task)) 
     >    call errquit(pname//'failed rtdb_cget',0,
     >       RTDB_ERR)

         if(task(1:4).eq."STOP") then
            lbfgsb_driver = .false.
            goto 112

         else if (task(1:4).eq."CONV") then
            lbfgsb_driver = .true.
            goto 112

         else
            call dcopy(nc,dbl_mb(i_xx),1,dbl_mb(i_x0),1)
            f0 = f
            goto 111
         end if

      end if 

c           ---------- the end of the loop -------------
 
c     If task is neither FG nor NEW_X we terminate execution.
c    
112   continue


      if (master) then
        if (.not.geom_print(igeom))
     >  call errquit(pname//'Failed to print final geometry',0,GEOM_ERR)
      end if


      status = .true.
      status = status .and. geom_destroy(igeom)
      status = status .and. ma_free_heap(h_wa)
      status = status .and. ma_free_heap(h_g)
      status = status .and. ma_free_heap(h_u)
      status = status .and. ma_free_heap(h_l)
      status = status .and. ma_free_heap(h_x0)
      status = status .and. ma_free_heap(h_x)
      status = status .and. ma_free_heap(h_iwa)
      status = status .and. ma_free_heap(h_nbd)
      status = status .and. ma_free_heap(h_index_unique)
       
      if(.not. status)   
     + call errquit( pname//'Failed to deallocate memory')

      return 

      end
c======================= The end of driver1 ============================

c     --------------------------------------------------------------
c             DESCRIPTION OF THE VARIABLES IN L-BFGS-B
c     --------------------------------------------------------------
c
c     n is an INTEGER variable that must be set by the user to the
c       number of variables.  It is not altered by the routine.
c
c     m is an INTEGER variable that must be set by the user to the
c       number of corrections used in the limited memory matrix.
c       It is not altered by the routine.  Values of m < 3  are
c       not recommended, and large values of m can result in excessive
c       computing time. The range  3 <= m <= 20 is recommended. 
c
c     x is a DOUBLE PRECISION array of length n.  On initial entry
c       it must be set by the user to the values of the initial
c       estimate of the solution vector.  Upon successful exit, it
c       contains the values of the variables at the best point
c       found (usually an approximate solution).
c
c     l is a DOUBLE PRECISION array of length n that must be set by
c       the user to the values of the lower bounds on the variables. If
c       the i-th variable has no lower bound, l(i) need not be defined.
c
c     u is a DOUBLE PRECISION array of length n that must be set by
c       the user to the values of the upper bounds on the variables. If
c       the i-th variable has no upper bound, u(i) need not be defined.
c
c     nbd is an INTEGER array of dimension n that must be set by the
c       user to the type of bounds imposed on the variables:
c       nbd(i)=0 if x(i) is unbounded,
c              1 if x(i) has only a lower bound,
c              2 if x(i) has both lower and upper bounds, 
c              3 if x(i) has only an upper bound.
c
c     f is a DOUBLE PRECISION variable.  If the routine setulb returns
c       with task(1:2)= 'FG', then f must be set by the user to
c       contain the value of the function at the point x.
c
c     g is a DOUBLE PRECISION array of length n.  If the routine setulb
c       returns with taskb(1:2)= 'FG', then g must be set by the user to
c       contain the components of the gradient at the point x.
c
c     factr is a DOUBLE PRECISION variable that must be set by the user.
c       It is a tolerance in the termination test for the algorithm.
c       The iteration will stop when
c
c        (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
c
c       where epsmch is the machine precision which is automatically
c       generated by the code. Typical values for factr on a computer
c       with 15 digits of accuracy in double precision are:
c       factr=1.d+12 for low accuracy;
c             1.d+7  for moderate accuracy; 
c             1.d+1  for extremely high accuracy.
c       The user can suppress this termination test by setting factr=0.
c
c     pgtol is a double precision variable.
c       On entry pgtol >= 0 is specified by the user.  The iteration
c         will stop when
c
c                 max{|proj g_i | i = 1, ..., n} <= pgtol
c
c         where pg_i is the ith component of the projected gradient.
c       The user can suppress this termination test by setting pgtol=0.
c
c     wa is a DOUBLE PRECISION  array of length 
c       (2mmax + 4)nmax + 12mmax^2 + 12mmax used as workspace.
c       This array must not be altered by the user.
c
c     iwa is an INTEGER  array of length 3nmax used as
c       workspace. This array must not be altered by the user.
c
c     task is a CHARACTER string of length 60.
c       On first entry, it must be set to 'START'.
c       On a return with task(1:2)='FG', the user must evaluate the
c         function f and gradient g at the returned value of x.
c       On a return with task(1:5)='NEW_X', an iteration of the
c         algorithm has concluded, and f and g contain f(x) and g(x)
c         respectively.  The user can decide whether to continue or stop
c         the iteration. 
c       When
c         task(1:4)='CONV', the termination test in L-BFGS-B has been 
c           satisfied;
c         task(1:4)='ABNO', the routine has terminated abnormally
c           without being able to satisfy the termination conditions,
c           x contains the best approximation found,
c           f and g contain f(x) and g(x) respectively;
c         task(1:5)='ERROR', the routine has detected an error in the
c           input parameters;
c       On exit with task = 'CONV', 'ABNO' or 'ERROR', the variable task
c         contains additional information that the user can print.
c       This array should not be altered unless the user wants to
c          stop the run for some reason.  See driver2 or driver3
c          for a detailed explanation on how to stop the run 
c          by assigning task(1:4)='STOP' in the driver.
c
c     iprint is an INTEGER variable that must be set by the user.
c       It controls the frequency and type of output generated:
c        iprint<0    no output is generated;
c        iprint=0    print only one line at the last iteration;
c        0<iprint<99 print also f and |proj g| every iprint iterations;
c        iprint=99   print details of every iteration except n-vectors;
c        iprint=100  print also the changes of active set and final x;
c        iprint>100  print details of every iteration including x and g;
c       When iprint > 0, the file iterate.dat will be created to
c                        summarize the iteration.
c
c     csave  is a CHARACTER working array of length 60.
c
c     lsave is a LOGICAL working array of dimension 4.
c       On exit with task = 'NEW_X', the following information is
c         available:
c       lsave(1) = .true.  the initial x did not satisfy the bounds;
c       lsave(2) = .true.  the problem contains bounds;
c       lsave(3) = .true.  each variable has upper and lower bounds.
c
c     isave is an INTEGER working array of dimension 44.
c       On exit with task = 'NEW_X', it contains information that
c       the user may want to access:
c         isave(30) = the current iteration number;
c         isave(34) = the total number of function and gradient
c                         evaluations;
c         isave(36) = the number of function value or gradient
c                                  evaluations in the current iteration;
c         isave(38) = the number of free variables in the current
c                         iteration;
c         isave(39) = the number of active constraints at the current
c                         iteration;
c
c         see the subroutine setulb.f for a description of other 
c         information contained in isave
c
c     dsave is a DOUBLE PRECISION working array of dimension 29.
c       On exit with task = 'NEW_X', it contains information that
c         the user may want to access:
c         dsave(2) = the value of f at the previous iteration;
c         dsave(5) = the machine precision epsmch generated by the code;
c         dsave(13) = the infinity norm of the projected gradient;
c
c         see the subroutine setulb.f for a description of other 
c         information contained in dsave
c
c     --------------------------------------------------------------
c           END OF THE DESCRIPTION OF THE VARIABLES IN L-BFGS-B
c     --------------------------------------------------------------
c
c     << An example of subroutine 'timer' for AIX Version 3.2 >>
c
c     subroutine timer(ttime)
c     double precision ttime
c     integer itemp, integer mclock
c     itemp = mclock()
c     ttime = dble(itemp)*1.0d-2
c     return
c     end
c-----------------------------------------------------------------------
