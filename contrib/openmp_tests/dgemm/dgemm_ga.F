C>
C> \brief Compute C=A*B in parallel using GA
C>
C>
      program dgemm_ga_prog
      implicit none
C
C     Find MPI definitions
C
#include "mpif.h"
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
C
C     Declare interfaces
C
      interface
         subroutine dgemm_ga(ndim,ldim,g_A,g_B,g_C,a,b,c)
           integer, parameter        :: w = selected_real_kind(12)
           integer, intent(in)       :: ndim ! total dimension
           integer, intent(in)       :: ldim ! local dimension
           integer, intent(in)       :: g_A  ! global array handle for A
           integer, intent(in)       :: g_B  ! global array handle for B
           integer, intent(in)       :: g_C  ! global array handle for C
           real(kind=w), intent(in)  :: a(ldim,ldim)
           real(kind=w), intent(in)  :: b(ldim,ldim)
           real(kind=w), intent(out) :: c(ldim,ldim)
         end subroutine dgemm_ga
      end interface
C
C     Data type
C
      integer, parameter :: w = 8
C
C     Local copies of matrices A, B, and C so we can check the results
C
      integer l_A1, k_A1
      integer l_B1, k_B1
      integer l_C1, k_C1
C
C     Local blocks of matrices A, B, and C
C
      integer l_A, k_A
      integer l_B, k_B
      integer l_C, k_C
C
C     Handles of A, B, and C to hold the distributed version of the
C     matrices
C
      integer :: g_A
      integer :: g_B
      integer :: g_C
C
C     Parameters of the job
C
      integer ndim  !< The dimension of the matrices
      integer ldim  !< The dimension of the local matrices
      integer ncore !< The number of cores
      integer icore !< My core rank
C
C     The core grid is defined as (icore_row,icore_col) and consists
C     of icore_dim x icore_dim (= ncore) cores.
C
      integer icore_dim
      integer icore_row
      integer icore_col
C
C     Miscelaneous
C
      integer ierr
      integer world !< For the MPI_COMM_WORLD communicator
      real(kind=w) :: tst
      integer itst
      integer itag
      integer narg !< The number of command line arguments
      integer istatus !< The status of the command reader
      integer ii, jj !< Counters
      integer global !< Size of shared memory in GA
      integer heap   !< Size of local heap memory in MA
      integer stack  !< Size of local stack memory in MA
      character*(80) carg !< A particular command line argument
      integer ilo, jlo, ihi, jhi !< Lower and upper limits of matrix
                                 !< blocks
      real(kind=w) :: time !< Wallclock time
C
      integer, external :: strlen
C
C     Code
C
      dbl_mb(1) = 0.0d0
      int_mb(1) = 0
      log_mb(1) = .false.
      call pbeginf()
      ndim = -1
      ncore = nnodes()
      icore = nodeid()
c
c     Read the command line options (we need to know how much 
c     memory we need)
c
      if (icore.eq.0) then ! TCGMSG as GA is not initialized yet
        narg = command_argument_count() !< Fortran 2003
        ii = 1
        do while (ii .le. narg)
          call get_command_argument(ii,carg,status=istatus) !< Fortran 2003
          if (carg(1:strlen(carg)).eq."--ndim") then
            call get_command_argument(ii+1,carg,status=istatus)
            read(carg,'(i)',iostat=istatus) ndim
            if (istatus.ne.0) then
              write(*,*)'*** ERROR: read failed, iostat = ',istatus,
     &                  ' carg = ',carg(1:strlen(carg))
            endif
            ii = ii + 2 !< Have read "--ndim <int>" which is two tokens
          else
            write(*,*)'*** ERROR: unknown command line option: ',
     &                carg(1:strlen(carg))
          endif
        enddo
      endif
      call brdcst(100,ndim,mitob(1),0)
      global = (4*ndim*ndim+ncore-1)/ncore ! shared memory per core
      heap   = 1000 ! local heap memory per core, just needs to be > 0
      if (icore.eq.0) then
        ! root needs memory for the full matrices
        stack = 3*ndim*ndim+3*(ndim*ndim+ncore-1)/ncore 
      else
        ! non-root needs memory for just matrix blocks
        stack = 4*(ndim*ndim+ncore-1)/ncore 
      endif
c
c     We need to initialize GA first before we can bomb out if the
c     problem dimension is invalid
c
      call ga_initialize_ltd(ma_sizeof(MT_DBL,global,MT_BYTE))
      if (ndim.le.0) call ga_error("illegal dimension",ndim)
      if ( ga_uses_ma() ) then
        if (.not.ma_init(mt_dbl, stack, heap+global))
     &      call ga_error('ma_init failed (ga_uses_ma=T)',
     &                    stack+heap+global)
      else
        if (.not.ma_init(mt_dbl,stack,heap))
     &      call ga_error('ma_init failed (ga_uses_ma=F)',
     &      stack+heap)
      endif
C
C     Check whether the number of cores generates the required
C     core grid
C
      tst = sqrt(1.0_w*ncore)
      itst = tst
      if (abs(itst*1.0_w-tst).gt.1.0e-10_w) then
        call ga_error('*** ERROR: core grid must be square',10)
      endif
C
C     Workout the dimensions of the core grid and my position in the
C     grid. Note that core grid coordinates run from (0,0) to 
C     (icore_dim-1,icore_dim-1).
C
      icore_dim = itst
      icore_row = mod(icore,icore_dim)
      icore_col = icore/icore_dim
C
C     Set up the source arrays
C
      if (icore.eq.0) then
C
C       Allocate and initialize the arrays A1, B1, and C1
C
        if (.not.ma_push_get(MT_DBL,ndim*ndim,"A1",l_A1,k_A1))
     &    call ga_error("could not allocate A1",
     &                  ma_sizeof(MT_DBL,ndim*ndim,MT_BYTE))
        if (.not.ma_push_get(MT_DBL,ndim*ndim,"B1",l_B1,k_B1))
     &    call ga_error("could not allocate B1",
     &                  ma_sizeof(MT_DBL,ndim*ndim,MT_BYTE))
        if (.not.ma_push_get(MT_DBL,ndim*ndim,"C1",l_C1,k_C1))
     &    call ga_error("could not allocate C1",
     &                  ma_sizeof(MT_DBL,ndim*ndim,MT_BYTE))
        call random_seed()
c       call random_number(reshape(dbl_mb(k_A1:k_A1-1+ndim*ndim),
c    &                     (/ndim,ndim/)))
c       call random_number(reshape(dbl_mb(k_B1:k_B1-1+ndim*ndim),
c    &                     (/ndim,ndim/)))
        call random_number(dbl_mb(k_A1:k_A1-1+ndim*ndim))
        call random_number(dbl_mb(k_B1:k_B1-1+ndim*ndim))
C
      endif
C
C     Now rank 0 has the full A1 and B1 matrices, these matrices will
C     have to be distributed across all ranks. With GA this is really
C     easy: we just need to create distributed arrays and pump the
C     data in.
C
      ldim = ndim/icore_dim
      if (.not.ga_create(MT_DBL,ndim,ndim,"g_A",ldim,ldim,g_A))
     &  call ga_error("could not create g_A",
     &                ma_sizeof(MT_DBL,ldim*ldim,MT_BYTE))
      if (.not.ga_create(MT_DBL,ndim,ndim,"g_B",ldim,ldim,g_B))
     &  call ga_error("could not create g_B",
     &                ma_sizeof(MT_DBL,ldim*ldim,MT_BYTE))
      if (.not.ga_create(MT_DBL,ndim,ndim,"g_C",ldim,ldim,g_C))
     &  call ga_error("could not create g_C",
     &                ma_sizeof(MT_DBL,ldim*ldim,MT_BYTE))
c
c     Now put the data in
c
      if (icore.eq.0) then
        call ga_put(g_A,1,ndim,1,ndim,dbl_mb(k_A1),ndim)
        call ga_put(g_B,1,ndim,1,ndim,dbl_mb(k_B1),ndim)
      endif
      call ga_zero(g_C)
c
c     Create local matrix blocks
c
      if (.not.ma_push_get(MT_DBL,ldim*ldim,"A",l_A,k_A))
     &  call ga_error("could not allocate A",
     &                ma_sizeof(MT_DBL,ldim*ldim,MT_BYTE))
      if (.not.ma_push_get(MT_DBL,ldim*ldim,"B",l_B,k_B))
     &  call ga_error("could not allocate B",
     &                ma_sizeof(MT_DBL,ldim*ldim,MT_BYTE))
      if (.not.ma_push_get(MT_DBL,ldim*ldim,"C",l_C,k_C))
     &  call ga_error("could not allocate C",
     &                ma_sizeof(MT_DBL,ldim*ldim,MT_BYTE))
    
      if (ldim*icore_dim.ne.ndim) then
        write(*,*)'*** ERROR: NDIM is not an integer multiple of LDIM: '
     &           ,ndim,icore_dim,ldim
        call ga_error("*** ERROR: NDIM not an integer multiple of LDIM",
     &                20)
      endif
C     
C     Get everyone on the same page, then do the work
C
      call ga_sync()
      time = -MPI_Wtime()
      call dgemm_ga(ndim,ldim,g_A,g_B,g_C,
     &              dbl_mb(k_A),dbl_mb(k_B),dbl_mb(k_C))
      time = time + MPI_Wtime()
C
C     Collect the resulting matrix C, also easy with GA
C
      if (icore.eq.0) then
C
C       - Rank 0 collects the data from all other ranks
C
        call ga_get(g_C,1,ndim,1,ndim,dbl_mb(k_C1),ndim)
        
      endif
C
      if (icore.eq.0) then
        call dgemm_i8('N','N',ndim,ndim,ndim,-1.0_w,dbl_mb(k_A1),ndim,
     &                dbl_mb(k_B1),ndim,1.0_w,dbl_mb(k_C1),ndim)
        tst = 0.0_w
        do ii = 0, ndim*ndim-1
          tst = tst + dbl_mb(k_C1+ii)**2
        enddo
        tst = sqrt(tst)
        write(*,*)"Norm = ",tst
        write(*,*)"Time = ",time
      endif
C
      if (.not.ma_pop_stack(l_C))
     &  call ga_error("could not release C",100)
      if (.not.ma_pop_stack(l_B))
     &  call ga_error("could not release B",110)
      if (.not.ma_pop_stack(l_A))
     &  call ga_error("could not release A",120)
      if (icore.eq.0) then
        if (.not.ma_pop_stack(l_C1))
     &    call ga_error("could not release C1",130)
        if (.not.ma_pop_stack(l_B1))
     &    call ga_error("could not release B1",140)
        if (.not.ma_pop_stack(l_A1))
     &    call ga_error("could not release A1",150)
      endif
      if (.not.ga_destroy(g_A))
     &  call ga_error("could not release g_A",160)
      if (.not.ga_destroy(g_B))
     &  call ga_error("could not release g_B",170)
      if (.not.ga_destroy(g_C))
     &  call ga_error("could not release g_C",180)
c
      call ga_terminate()
      call pend()
C
      end
C>
C> \brief Return the length of a string
C> 
C> Returns the position of the last non-white-space character
C>
      integer function strlen(string)
      implicit none
      character*(*) string !< [Input] The character string
      integer ii !< Counter
      ii = len(string)
      do while(string(ii:ii).eq." ".and.ii.gt.0)
        ii = ii - 1
      enddo
      strlen = ii
      return
      end
C>
C> \brief The actual parallel matrix-matrix multiply routine
C>
C> This routine calculates \f$C = A*B\f$ where A, B, and C are
C> distributed matrices. The calculation is organized so as to compute
C> the local block of matrix C, and the block of matrices A and B
C> are communicated as needed.
C>
      subroutine dgemm_ga(ndim,ldim,g_A,g_B,g_C,a,b,c)
      use mpi
      implicit none
C
#include "global.fh"
C
      integer ndim !< [Input] The dimension of the distributed matrices
      integer ldim !< [Input] The dimension of the local blocks
      integer g_A  !< [Input] Handle of GA A
      integer g_B  !< [Input] Handle of GA B
      integer g_C  !< [Input] Handle of GA C
C
C     Sort out data types
C
      integer, parameter :: w = 8
C
C     Local blocks of matrices A, B, and C
C
      real(kind=w), intent(inout) :: A(ldim,ldim)
      real(kind=w), intent(inout) :: B(ldim,ldim)
      real(kind=w), intent(out)   :: C(ldim,ldim)
C
      integer col_comm
      integer row_comm
      integer col_comm_irank
      integer row_comm_irank
      integer world_irank
C
      integer ncore !< The number of cores
      integer icore !< My core rank
C
C     The core grid is defined as (icore_row,icore_col) and consists
C     of icore_dim x icore_dim (= ncore) cores.
C
      integer icore_dim
      integer icore_row
      integer icore_col
C
      integer ilo, ihi, jlo, jhi
      integer world
      integer key
      integer itst
      integer ierr
      integer ii, jj !< Counters
C
      ncore = ga_nnodes()
      icore = ga_nodeid()
C
C     Check whether the number of cores generates the required
C     core grid
C
      itst = sqrt(1.0_w*ncore)
C
C     Workout the dimensions of the core grid and my position in the
C     grid. Note that core grid coordinates run from (0,0) to 
C     (icore_dim-1,icore_dim-1).
C
      icore_dim = itst
      icore_row = mod(icore,icore_dim)
      icore_col = icore/icore_dim
C
C
C     A bit confusing:
C     - The rank in the row communicator corresponds to the column
C     - The rank in the column communicator corresponds to the row
C
      C = 0.0_w
C
      do ii = 0, icore_dim-1
C
C       Fetch a block of A
C
        ilo = icore_row*ldim+1
        ihi = icore_row*ldim+ldim
        jlo = ii       *ldim+1
        jhi = ii       *ldim+ldim
        call ga_get(g_A,ilo,ihi,jlo,jhi,a,ldim)
C
C       Fetch a block of B
C
        ilo = ii       *ldim+1
        ihi = ii       *ldim+ldim
        jlo = icore_col*ldim+1
        jhi = icore_col*ldim+ldim
        call ga_get(g_B,ilo,ihi,jlo,jhi,b,ldim)
C
        call dgemm_i8('N','N',ldim,ldim,ldim,1.0_w,A,ldim,B,ldim,
     +                1.0_w,C,ldim)
      enddo
C
C     Store the result
C
      ilo = icore_row*ldim+1
      ihi = icore_row*ldim+ldim
      jlo = icore_col*ldim+1
      jhi = icore_col*ldim+ldim
      call ga_put(g_C,ilo,ihi,jlo,jhi,c,ldim)
C
C     Synchronize
C
      call ga_sync()
C
      end
