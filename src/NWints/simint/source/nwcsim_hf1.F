C#define DEBUG 1
#ifdef USE_SIMINT
      subroutine nwcsim_hf1(
     n     ncenters,
     c     charges,
     x     x, y, z,
     I     i_basis,i_sh,
     J     j_basis,j_sh,
     S     s, t, v, leri,
     L     dos, dok, dov)
c      use SimintFortran
      use nwcsim_shell
      implicit none
      integer ncenters
      double precision charges(*),x(*),y(*),z(*)
      integer i_basis        ! basis for bra [in]
      integer j_basis       ! basis for ket [in]
      integer i_sh,j_sh
      logical dos, dok, dov
      integer ncomputed
      double precision s(*)
      double precision t(*)
      double precision v(*)
      integer leri
#include "errquit.fh"
#include "mafdecls.fh"
      integer i
      character*10 pname
      parameter(pname="nwcsim_hf1")
#ifdef USE_OPENMP
      integer omp_get_max_threads
      external omp_get_max_threads
#ifdef DEBUG      
      i=omp_get_max_threads()
      write(6,*) 'nwcsim_hf1: omp_get_max_thread=',i
#endif      
#endif      

      if( dos) then
         ncomputed = simint_compute_overlap(
     I        smnt_sh(i_sh,i_basis),
     J        smnt_sh(j_sh,j_basis),
     C        s)
      endif
      
      if( dok) then
         ncomputed = simint_compute_ke(
     I        smnt_sh(i_sh,i_basis),
     J        smnt_sh(j_sh,j_basis),
     C        t)
      endif
      
      if( dov ) then
         ncomputed = simint_compute_potential(
     N        ncenters,
     C        charges,
     X        x, y, z,
     I        smnt_sh(i_sh,i_basis),
     J        smnt_sh(j_sh,j_basis),
     C        v)
      endif
      
#ifdef DEBUG
      if(dos) then
         write(6,*) ncomputed,' seri ', s(1),s(ncomputed)
         endif
      if(dok) then
         write(6,*) ncomputed,' teri ', t(1),t(ncomputed)
         endif
      if(dov) then
         write(6,*) ncomputed,' veri ', v(1),v(ncomputed)
         endif
#endif
      if(ncomputed.gt.leri)  then
         write(6,*) ' simint: neri=',ncomputed, ' leri=',leri
          call errquit(
     P        pname//'wrong simint buffer size ',ncomputed, BASIS_ERR)
          endif
      if(ncomputed.le.0) call errquit(
     P     pname//'return value from simint ',ncomputed, BASIS_ERR)
cwrong      if(ncomputed.gt.leri) call errquit(
c     P     pname//'ncomputed gt leri ',ncomputed*10000+leri, BASIS_ERR)
 
      return
      end
#else
      subroutine nwcsim_hf1()
      implicit none
#include "errquit.fh"
      call errquit(
     $   'simint is not in this build of NWChem',
     $    0, CAPMIS_ERR)
      return
      end
#endif
