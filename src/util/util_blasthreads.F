c     
      subroutine util_blas_set_num_threads(nt)
      implicit none
      integer  nt ! [in]
#if defined(OPENBLAS)
      call openblas_set_num_threads(nt)
#elif defined(CRAYBLAS)
      call crayblas_set_num_threads(nt)
#elif defined(MKL)
      call mkl_set_num_threads(nt)
#elif defined(BLIS)
      call bli_thread_set_num_threads(nt)
#elif defined(FLEXIBLAS)
      call flexiblas_set_num_threads(nt)
#elif defined(ACCELERATE)
      call veclib_set_num_threads(nt)
#elif defined(ARMPL)
      call armpl_set_num_threads(nt)
#elif defined(BLAS_OPENMP)
      call omp_set_num_threads(nt)
#elif defined(BLAS_NOTHREADS)
c do nothing      
#else
      call errquit(' util_blas not ready ',0,0)
#endif
      return
      end
      subroutine util_blas_get_num_threads(nt)
      implicit none
      integer  nt ! [in]
#if defined(OPENBLAS)
      external openblas_get_num_threads
      integer openblas_get_num_threads
      nt=openblas_get_num_threads()
#elif defined(CRAYBLAS)
      external crayblas_get_max_threads
      integer crayblas_get_max_threads
      nt=crayblas_get_max_threads()
#elif defined(MKL) 
      external mkl_get_max_threads
      integer mkl_get_max_threads
      nt=mkl_get_max_threads()
#elif defined(BLIS)
      interface
         function bli_thread_get_num_threads() 
     c        result(res) bind(C, name="bli_thread_get_num_threads")
         use iso_c_binding
         implicit none
         integer(C_INT) :: res
         end function bli_thread_get_num_threads
      end interface
      nt=bli_thread_get_num_threads()
c     return value=-1 -> threading not enabled
#elif defined(FLEXIBLAS)
      external flexiblas_get_num_threads
      integer flexiblas_get_num_threads
      nt=flexiblas_get_num_threads()
#elif defined(ACCELERATE)
      external veclib_get_num_threads
      integer veclib_get_num_threads
      nt=veclib_get_num_threads()
#elif defined(ARMPL)
      external armpl_get_num_threads
      integer armpl_get_num_threads
      nt=armpl_get_num_threads()
#elif defined(BLAS_OPENMP)
      integer omp_get_max_threads
      external omp_get_max_threads
      nt=omp_get_max_threads()
#elif defined(BLAS_NOTHREADS)
      nt=1
#else
      call errquit(' util_blas not ready ',0,0)
#endif
      return
      end
      subroutine util_setthreads_fromenv()
      implicit none
c
      integer nt
      character(64) env_string
      integer stat
      nt = -1
c
      call util_getenv('NWC_NUM_THREADS', env_string)
      if ((len(env_string).ne.0).or.(env_string.ne.' ')) then
        read (env_string,*,iostat=stat) nt
      endif
      if (stat.ne.0) nt = 1
      if (nt.lt.1) nt = 1
      write(6,*) ' NWC_NUM_THREADS is ',nt
      call util_blas_set_num_threads(nt)
      return
      end
