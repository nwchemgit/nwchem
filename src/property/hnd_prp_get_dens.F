      subroutine hnd_prp_get_dens(rtdb,geom,basis,g_dens,ndens,scftyp,
     2                            nclosed,nopen,nvirt)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "util.fh"
c
c     Assumes energy has been completed, MO vectors stored
c     and all information is still in the RTDB
c
      integer     rtdb          ! [input] database handle
      integer     geom          ! [input] geometry handle
      integer     basis         ! [input] handles to basis
      integer     g_dens(3)     ! [output] handles to density
      integer     ndens         ! [output] number of active density handles (RHF=1, UHF=3)
      character*3 scftyp        ! [output] type of wave function
      integer nclosed(2),nopen(2),nvirt(2) ! [output] occupation info
c
      integer nbf, nmo, k_occ, l_occ, k_eval, l_eval, vectors(2)
c
c     Get vectors and other information
c     Arrays occ(nbf*2) and evals(nbf*2) are needed
c
      if (.not. bas_numbf(basis,nbf)) call 
     &    errquit('hnd_prp_get_dens: could not get nbf',0, BASIS_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*2,'MO eval',l_eval,k_eval)) call
     &    errquit('hnd_prp_get_dens: ma_push_get l_eval',0, MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*2,'MO occ',l_occ,k_occ)) call
     &    errquit('hnd_prp_get_dens: ma_push_get l_occ',0, MA_ERR)
c
      call hnd_prp_vec_read(rtdb,geom,basis,nbf,nclosed,nopen,nvirt,
     &                      scftyp,vectors,dbl_mb(k_occ),
     &                      dbl_mb(k_eval),nmo)
c
c     Make the density matrix
c
      call hnd_prop_dens_make(rtdb,geom,basis,nbf,nmo,nclosed,nopen,
     &                       nvirt,scftyp,vectors,dbl_mb(k_occ),
     &                       g_dens,ndens)
c
c     Cleanup of MA arrays that are not needed
c
      if (.not.ma_pop_stack(l_occ)) call
     &    errquit('hnd_prp_get_dens: ma_pop_stack l_occ',0, MA_ERR)
      if (.not.ma_pop_stack(l_eval)) call
     &    errquit('hnd_prp_get_dens: ma_pop_stack l_eval',0, MA_ERR)
c
      return
      end
