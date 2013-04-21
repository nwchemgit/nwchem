C> \ingroup nwint
C> @{
C>
C> \brief Normalize the contraction coefficients in a basis set
C>
C> The contraction coefficients in basis sets express a basis function
C> in terms of normalized Gaussian functions. However, basis functions
C> as specified with the exponent are unnormalized. A convenient way
C> to address this is to calculate and absorp the normalization factor
C> in the contraction coefficients. 
C>
C> This function calculates the normalization constants and a set of
C> contraction coefficients that express basis functions in terms of
C> unnormalized Gaussians. The unnormalized contraction coefficients
C> are stored as well because they are used when the basis is printed.
C> The unnormalized contraction coefficients correspond to the data
C> a user entered and therefore yield the basis set representation 
C> that can be recognized as matching the input data. This function
C> has to be called to ensure that the integral codes calculate
C> integrals of normalized Gaussian basis functions. Typically this
C> function is called immediately after loading a basis set from 
C> the runtime database:
C>
C> \code
C>       if (.not. bas_create(basis,'ao basis'))
C>      &    call errquit("ddscf: bas_create failed",0,BASIS_ERR)
C>       if (.not. bas_rtdb_load(rtdb,geom,basis,'ao basis'))
C>      &    call errquit("ddscf: basis load failed",0,RTDB_ERR)
C>       if (.not. int_normalize(rtdb,basis))
C>      &    call errquit("ddscf: normalization failed",0,INT_ERR)
C> \endcode
C>
C> \return Returns .true. if the basis set was normalized successfully,
C> returns .false. otherwise.
C>
      logical function int_normalize(rtdb,basisin)
c $Id$
      implicit none
#include "apiP.fh"
#include "errquit.fh"
#include "bas.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "stdio.fh"
      integer basisin  !< [Input] The basis set handle
      integer rtdb     !< [Input] The runtime database handle
      integer basis, nucont, iucont, type, nprim, ngen, ept, cpt, igen
      integer norm_typ
      logical norm_ok
c
#include "bas_exndcf_sfn.fh"
c
c     Normalize the contraction coefficients in this basis.
c     Eventually we need to store both the normalized and unnormalized
c     coefficients so that the user can print out something recognizable
c
      int_normalize = bas_check_handle(basisin, 'int_normalize')
      if (.not. int_normalize) return
c
      basis = basisin + BASIS_HANDLE_OFFSET
c
c if this is a relativistic basis set, we need the primitive normalization
c but not the general contractions.
c     
      if (.not. rtdb_get(rtdb,'dyall_mod_dir',mt_log,1,dyall_mod_dir))
     &    dyall_mod_dir = .false.
      if (.not. rtdb_get(rtdb,'nesc_1e_approx',mt_log,1,nesc_1e_approx))
     &    nesc_1e_approx = .false.
      if (dyall_mod_dir .and.
     &    ((basisin. eq. ao_bsh) .or. (basisin .eq. sc_bsh) .or.
     &     (basisin .eq. lc_bsh .and. nesc_1e_approx))) then
        norm_typ = BasNorm_rel
      else
        norm_typ = BasNorm_STD
      end if
c
c check and set basis normalization
c
      norm_ok = bas_norm_id(basis).eq.norm_typ
      if (norm_ok) return ! basis has been normalized properly
      norm_ok = bas_norm_id(basis).eq.BasNorm_UN
      if (.not.norm_ok) then
        if (.not. bas_norm_print(basisin)) call errquit
     &      ('int_normalize: bas_norm_print failed',911, BASIS_ERR)
        write(luout,*)' basis has already been normalized somehow'
        call errquit('int_norm: basis normalization error ',911,
     &               BASIS_ERR)
      endif
      bas_norm_id(basis) = norm_typ
c
      nucont = infbs_head(HEAD_NCONT,basis)
      do iucont = 1, nucont
         type  = infbs_cont(CONT_TYPE, iucont,basis)
         nprim = infbs_cont(CONT_NPRIM,iucont,basis) 
         ngen  = infbs_cont(CONT_NGEN, iucont,basis)
         ept   = infbs_cont(CONT_IEXP, iucont,basis)
         cpt   = infbs_cont(CONT_ICFP, iucont,basis)
         if (type .ge. 0) then
            do igen = 0,ngen-1
               call nmcoef(dbl_mb(mb_exndcf(ept, basis)),
     &              dbl_mb(mb_exndcf(cpt+igen*nprim,basis)),
     &              type, nprim, norm_typ)
            enddo
c....... sp
         else if (type .eq. -1) then
           call nmcoef(dbl_mb(mb_exndcf(ept, basis)),
     &           dbl_mb(mb_exndcf(cpt, basis)),0, nprim, norm_typ) 
           call nmcoef(dbl_mb(mb_exndcf(ept, basis)),
     &           dbl_mb(mb_exndcf(cpt+nprim,basis)), 1, nprim, norm_typ)
         else if (type .eq. -2) then
c........ spd
           call nmcoef(dbl_mb(mb_exndcf(ept, basis)),
     &           dbl_mb(mb_exndcf(cpt, basis)), 0, nprim, norm_typ) 
            call nmcoef(dbl_mb(mb_exndcf(ept, basis)),
     &           dbl_mb(mb_exndcf(cpt+nprim,basis)),
     &           1, nprim, norm_typ)
            call nmcoef(dbl_mb(mb_exndcf(ept, basis)),
     &            dbl_mb(mb_exndcf(cpt+2*nprim,basis)),
     &           2, nprim, norm_typ)
         else
            call errquit('int_normalize: invalid contraction type', 0,
     &       BASIS_ERR)
         endif
      enddo
c
      end
C>
C> @}
