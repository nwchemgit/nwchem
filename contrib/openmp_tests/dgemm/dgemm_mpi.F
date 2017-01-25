C>
C> \brief Compute C=A*B in parallel using MPI
C>
C>
      program dgemm_mpi_prog
      use mpi
      implicit none
C
C     Find MPI definitions
C
cinclude "mpif.h"
C
C     Declare interfaces
C
      interface
         subroutine dgemm_mpi(ndim,ldim,A,B,C)
           integer, parameter :: w = selected_real_kind(12)
           integer,      intent(in)    :: ndim
           integer,      intent(in)    :: ldim
           real(kind=w), intent(in)    :: A(:,:)
           real(kind=w), intent(in)    :: B(:,:)
           real(kind=w), intent(out)   :: C(:,:)
         end subroutine dgemm_mpi
      end interface
C
C     Data type
C
      integer, parameter :: w = 8
C
C     Local copies of matrices A, B, and C so we can check the results
C
      real(kind=w), allocatable :: A1(:,:)
      real(kind=w), allocatable :: B1(:,:)
      real(kind=w), allocatable :: C1(:,:)
C
C     Blocks of A, B, and C to hold the distributed version of the
C     matrices
C
      real(kind=w), allocatable :: A(:,:)
      real(kind=w), allocatable :: B(:,:)
      real(kind=w), allocatable :: C(:,:)
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
      integer mpistat(MPI_STATUS_SIZE) !< The MPI status
      integer ii, jj !< Counters
      character*(80) carg !< A particular command line argument
      integer ilo, jlo, ihi, jhi !< Lower and upper limits of matrix
                                 !< blocks
      real(kind=w) :: time !< Wallclock time
C
      integer, external :: strlen
C
C     Code
C
      world = MPI_COMM_WORLD
      call MPI_Init(ierr)
      call MPI_Comm_size(world,ncore,ierr)
      call MPI_Comm_rank(world,icore,ierr)
C
C     Check whether the number of cores generates the required
C     core grid
C
      tst = sqrt(1.0_w*ncore)
      itst = tst
      if (abs(itst*1.0_w-tst).gt.1.0e-10_w) then
        write(*,*)'*** ERROR: core grid must be square'
        call MPI_Abort(world,10,ierr)
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
C     Process the command line options
C
      if (icore.eq.0) then
        narg = command_argument_count() !< Fortran 2003
        ii = 1
        ndim = icore_dim
        do while (ii .le. narg)
          call get_command_argument(ii,carg,status=istatus) !< Fortran 2003
          if (carg(1:strlen(carg)).eq."--ndim") then
            call get_command_argument(ii+1,carg,status=istatus)
            read(carg,'(i10)',iostat=istatus) ndim
            if (istatus.ne.0) then
              write(*,*)'*** ERROR: read failed, iostat = ',istatus,
     &                  ' carg = ',carg(1:strlen(carg))
              call MPI_Abort(world,20,ierr)
            endif
            ii = ii + 2 !< Have read "--ndim <int>" which is two tokens
          else
            write(*,*)'*** ERROR: unknown command line option: ',
     &                carg(1:strlen(carg))
            call MPI_Abort(world,30,ierr)
          endif
        enddo
C
C       Allocate and initialize the arrays A1, B1, and C1
C
        allocate(a1(1:ndim,1:ndim),b1(1:ndim,1:ndim),c1(1:ndim,1:ndim))
        call random_seed()
        call random_number(a1)
        call random_number(b1)
C
      endif
C
C     Now rank 0 has the full A1 and B1 matrices, these matrices will
C     have to be distributed across all ranks.
C
C     - first tell everyone how big the problem is
C
      call MPI_Bcast(ndim,1,MPI_INTEGER,0,world,ierr)
      ldim = ndim/icore_dim
      if (ldim*icore_dim.ne.ndim) then
        write(*,*)'*** ERROR: NDIM is not an integer multiple of LDIM: '
     &           ,ndim,icore_dim,ldim
        call MPI_Abort(world,40,ierr)
      endif
      allocate(a(1:ldim,1:ldim),b(1:ldim,1:ldim),c(1:ldim,1:ldim))
C     
      if (icore.eq.0) then
C
C       - Rank 0 sends the blocks the everyone
C
        do jj = icore_dim-1,0,-1
          do ii = icore_dim-1,0,-1
            ilo = 1+(ii  )*ldim
            ihi =   (ii+1)*ldim
            jlo = 1+(jj  )*ldim
            jhi =   (jj+1)*ldim
            a(1:ldim,1:ldim) = a1(ilo:ihi,jlo:jhi)
            b(1:ldim,1:ldim) = b1(ilo:ihi,jlo:jhi)
            itst = jj*icore_dim+ii
            if (ii.ne.0.or.jj.ne.0) then
              itag = 888
              call MPI_Send(a(:,:),ldim*ldim,
     &                      MPI_DOUBLE_PRECISION,itst,itag,world,ierr)
              itag = 889
              call MPI_Send(b(:,:),ldim*ldim,
     &                      MPI_DOUBLE_PRECISION,itst,itag,world,ierr)
            endif
          enddo
        enddo
      else
C
C       - All ranks .ne. 0 receive their blocks from rank 0
C
        itst = 0
        itag = 888
        call MPI_Recv(a(:,:),ldim*ldim,
     &                MPI_DOUBLE_PRECISION,itst,itag,world,mpistat,ierr)
        itst = 0
        itag = 889
        call MPI_Recv(b(:,:),ldim*ldim,
     &                MPI_DOUBLE_PRECISION,itst,itag,world,mpistat,ierr)
      endif
C
C     Get everyone on the same page, then do the work
C
      call MPI_Barrier(world,ierr)
      time = -MPI_Wtime()
      call dgemm_mpi(ndim,ldim,a,b,c)
      time = time + MPI_Wtime()
C
C     Collect the resulting matrix C
C
      if (icore.eq.0) then
C
C       - Rank 0 collects the data from all other ranks
C
        c1(1:ldim,1:ldim) = c(1:ldim,1:ldim)
        do jj = icore_dim-1,0,-1
          do ii = icore_dim-1,0,-1
            ilo = 1+(ii  )*ldim
            ihi =   (ii+1)*ldim
            jlo = 1+(jj  )*ldim
            jhi =   (jj+1)*ldim
            itst = jj*icore_dim+ii
            if (ii.ne.0.or.jj.ne.0) then
              itag = 890
              call MPI_Recv(c(:,:),ldim*ldim,MPI_DOUBLE_PRECISION,
     &                      itst,itag,world,mpistat,ierr)
              c1(ilo:ihi,jlo:jhi) = c(1:ldim,1:ldim)
            endif
          enddo
        enddo
      else
        itst = 0
        itag = 890
        call MPI_Send(c(:,:),ldim*ldim,
     &                MPI_DOUBLE_PRECISION,itst,itag,world,ierr)
      endif
C
      if (icore.eq.0) then
        call dgemm('N','N',ndim,ndim,ndim,-1.0_w,A1,ndim,B1,ndim,1.0_w,
     &             C1,ndim)
        tst = 0.0_w
        do jj = 1, ndim
          do ii = 1, ndim
            tst = tst + C1(ii,jj)**2
          enddo
        enddo
        tst = sqrt(tst)
        write(*,*)"Norm = ",tst
        write(*,*)"Time = ",time
      endif
C
      call MPI_Finalize(ierr)
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
      subroutine dgemm_mpi(ndim,ldim,a,b,c)
      use mpi
      implicit none
C
cinclude "mpif.h"
C
      integer, intent(in) :: ndim !< [Input] The dimension of the distributed matrices
      integer, intent(in) :: ldim !< [Input] The dimension of the local blocks
C
C     Sort out data types
C
      integer, parameter :: w = selected_real_kind(12)
C
C     Local blocks of matrices A, B, and C
C
      real(kind=w), intent(in)  :: A(:,:)
      real(kind=w), intent(in)  :: B(:,:)
      real(kind=w), intent(out) :: C(:,:)
c
      real(kind=w) :: lA(ldim,ldim) !< local copy of a block of A
      real(kind=w) :: lB(ldim,ldim) !< local copy of a block of B
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
      integer world
      integer key
      integer itst
      integer ierr
      integer ii, jj !< Counters
C
      world = MPI_COMM_WORLD
      key   = MPI_UNDEFINED
      call MPI_Comm_size(world,ncore,ierr)
      call MPI_Comm_rank(world,icore,ierr)
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
      call MPI_Comm_split(world,icore_row,key,row_comm,ierr)
      call MPI_Comm_split(world,icore_col,key,col_comm,ierr)
C
      call MPI_Comm_rank(world,world_irank,ierr)
C
C     A bit confusing:
C     - The rank in the row communicator corresponds to the column
C     - The rank in the column communicator corresponds to the row
C
      call MPI_Comm_rank(row_comm,row_comm_irank,ierr)
      call MPI_Comm_rank(col_comm,col_comm_irank,ierr)
C
      C = 0.0_w
C
      do ii = 0, icore_dim-1
C
C       Communicate a block of A down the row
C
        if(ii.eq.row_comm_irank) then
C
C         I send my block of A
C
          lA(1:ldim,1:ldim) = A(1:ldim,1:ldim)
        else
C
C         I receive a block of A
C
        endif
        call MPI_Bcast(lA(:,:),ldim*ldim,MPI_DOUBLE_PRECISION,
     +                 ii,row_comm,ierr)
C
C       Communicate a block of B down the column
C
        if(ii.eq.col_comm_irank) then
C
C         I send my block of B
C
          lB(1:ldim,1:ldim) = B(1:ldim,1:ldim)
        else
C
C         I receive a block of B
C
        endif
        call MPI_Bcast(lB(:,:),ldim*ldim,MPI_DOUBLE_PRECISION,
     +                 ii,col_comm,ierr)
C
        call dgemm('N','N',ldim,ldim,ldim,1.0_w,lA,ldim,lB,ldim,
     +             1.0_w,C,ldim)
      enddo
C
C     Clean up communicators
C
      call MPI_Comm_free(row_comm,ierr)
      call MPI_Comm_free(col_comm,ierr)
      call MPI_Barrier(world,ierr)
C
      end
