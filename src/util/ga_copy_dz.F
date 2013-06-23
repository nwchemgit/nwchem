C
C $Id$
C
C> \brief Copy data from a double precision to a double complex array
C>
C> When calculating the Kohn-Sham wavefunction we need the density
C> in a double complex array. However, typically the density matrix
C> will be stored in a double precision array. This function copies
C> the data from the original double precision array into the 
C> double complex replica.
C>
C> \return Return .true. if successfull, and .false. otherwise
c
      logical function ga_copy_dz(g_in,g_out)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
#ifndef MAPOINTER
#define MAPOINTER integer*8
#endif
      integer g_in  !< [Input] Double precision global array
      integer g_out !< [Output] Double complex global array
c
c     Local
c
      integer iproc !< My processor rank
      integer nproc !< The total number processors
      integer num   !< The number of elements
      integer ndim  !< The number of dimensions
      integer i     !< Counters
c
      integer ilo(ga_max_dim) !< The block lower limits
      integer ihi(ga_max_dim) !< The block higher limits
      integer ild(ga_max_dim) !< The block leading dimensions
c
      integer   l_t !< Memory handle for buffer
      MAPOINTER k_t !< Index of buffer memory
      MAPOINTER k_m !< Index of GA memory
      logical status
      integer type_in, type_out
c
      ga_copy_dz = .false.
c
      iproc = ga_nodeid()
      nproc = ga_nnodes()
      ndim  = ga_ndim(g_in)
c
      call nga_inquire(g_in,type_in,ndim,ild)
      if (type_in.ne.MT_DBL) 
     $  call errquit("ga_copy_dz: input not real",type_in,UERR)
      call nga_inquire(g_out,type_out,ndim,ild)
      if (type_out.ne.MT_DCPL) 
     $  call errquit("ga_copy_dz: output not complex",type_in,UERR)
      call nga_distribution(g_in,iproc,ilo,ihi)
      num = 1
      do i = 1, ndim
         ild(i) = max(0,(ihi(i)-ilo(i)+1))
         num = num*ild(i)
      enddo
      status = ma_push_get(MT_DBL,num,"ks:in:buf",l_t,k_t)
      if (.not.status) return
      call nga_access(g_out,ilo,ihi,k_m,ild)
      call nga_get(g_in,ilo,ihi,dbl_mb(k_t),ild)
      do i = 0, num-1
         dcpl_mb(k_m+i) = cmplx(dbl_mb(k_t+i),0.0d0)
      enddo
      call nga_release(g_out,ilo,ihi)
      status = ma_pop_stack(l_t)
c
      ga_copy_dz = status
c
      return
      end
