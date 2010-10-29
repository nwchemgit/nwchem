      subroutine mpiInit
*
* $Id$
*
      implicit none
#include "myprocs.fh"
      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, me, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nproc, ierr )
      RETURN
      END

      subroutine mpiStop
      implicit none
#include "myprocs.fh"
      call MPI_FINALIZE(istatus)
      RETURN
      END

      subroutine mysync()
      implicit none
#include "myprocs.fh"
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
      RETURN
      END

      double precision FUNCTION myclock()
      implicit none
      double precision mpi_wtime
      myclock = mpi_wtime()
      RETURN
      END

#ifdef MPIPROCS
      integer function nnodes()
      implicit none
#include "myprocs.fh"
      nnodes = nproc
      return
      end 
#endif

#ifdef MPIPROCS
      integer function nodeid()
      implicit none
#include "myprocs.fh"
      nodeid = me 
      return
      end 
#endif


