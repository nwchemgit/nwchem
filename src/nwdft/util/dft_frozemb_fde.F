cc AJL/Begin/FDE
      subroutine dft_frozemb_fde(rtdb, g_frozemb, g_dens_fde) 

      implicit none

#include "rtdb.fh"
#include "bas.fh"
#include "cdft.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "msgtypesf.h"
#include "msgids.fh"
#include "stdio.fh"
#include "tcgmsg.fh"
c
c #include "dft_conv.fh"

      integer rtdb
      integer g_frozemb

      character*32 pname
      character*255 basis_trans
      character*255 basis_name
      character*255 frozemb_movecs
      character*255 frozemb_movecs_out
      character*7 vecs_or_dens
      character*3 scftype

      integer nExc_fde
      integer g_vxc_fde(4)
      integer itmat
      integer natoms

      integer nelectrons_fde
      integer nclosed_fde
      integer nopen_fde
      integer mult_fde
      logical even_fde

      integer g_dens_fde(2)

      double precision nucl_charge_fde
c      double precision Ecoul_fde
      double precision Exc_fde(2)
      double precision rho_n

      double precision tol2e_sleazy
      integer nmo_fde(2)
      integer noc_fde(2)
      integer k_eval(2)
      integer l_eval
      integer k_ir
      integer l_ir
      integer k_occ
      integer l_occ
      integer k_svals
      integer l_svals
      integer i
      integer ncanorg

      double precision toll_s
      integer n_dep
   
      integer g_movecs_fde(2)
      integer g_gmovecs_fde(2)
  
      logical ldmix

      integer inp_strlen
      external inp_strlen

      Exc_fde(1) = 0d0
      Exc_fde(2) = 0d0
      rho_n = 0d0
      nExc_fde = idecomp + 1
      itmat = 0

      write(luout,*) ''
      write(luout,*) '== Preparing Frozen Embedding Matrix =='

      pname = 'dft_frozemb_fde:'
c
c Nuclear component; easy part
c
      call int_1e_fde(AO_bas_han, AO_bas_han, g_frozemb,
     &               'fde_pot', oskel, FD_bas_han)

      write(luout,*) '-> Nuclear contributions complete'
c
c Setup up some empty arrays for the xc_getv.
c Whilst this functionality is not strictly functional right not,
c we might want the features in the future.
c  
      if (.not.MA_Push_Get(MT_Dbl, nbf_fde, 'FDE ovl eig vals',
     &                     l_svals, k_svals))
     &   call errquit(pname//'Cannot allocate ovl eig vals', 0,
     &       MA_ERR)
      call dfill(nbf_fde, 0.0d0, dbl_mb(k_svals), 1)
c
c First up, we need to construct the density matrix for the FDE.
c Create the GA to hold the FDE density matrix
      if (.not. ga_create(mt_dbl, nbf_fde, nbf_fde, 'FDE Alpha Dens',
     &        0, 0, g_dens_fde(1)))
     &        call errquit(pname//'Error creating ga',0,GA_ERR)
      call ga_zero(g_dens_fde(1))
      if (.not. ga_duplicate(g_dens_fde(1),
     &                       g_movecs_fde(1), 'FDE Alpha Movecs'))
     &     call errquit(pname//'Error creating ga',0,GA_ERR)
      call ga_zero(g_movecs_fde(1))
      if (ipol.gt.1) then
      if(.not.ga_duplicate(g_dens_fde(1),g_dens_fde(2),'FDE Beta Dens'))
     &     call errquit(pname//'ga_duplicate failed',1, GA_ERR)
        call ga_zero(g_dens_fde(2))
      if (.not. ga_duplicate(g_dens_fde(2),
     &                       g_movecs_fde(2), 'FDE Beta Movecs'))
     &     call errquit(pname//'Error creating ga',0,GA_ERR)
      call ga_zero(g_movecs_fde(2))
      endif

c Lifted from dft_scf.F
c
c     get/set orbital overlap tolerance
      call util_tolls(rtdb,.false.,toll_s,n_dep,
     &     dbl_mb(k_svals),nbf_fde)

c     == get info for int2e_ and set sleazy tolerance ==
      if (.not.rtdb_get(rtdb,'dft:tol2e_sleazy',
     &     mt_dbl,1,tol2e_sleazy)) tol2e_sleazy = 1.d-4

      if (.not. bas_name(fd_bas_han, basis_name, basis_trans))
     $     call errquit(pname//' bas_name?', 0, BASIS_ERR)

c     == allocate eigenvalue array, including second pointer to beta block ==
      if (.not.MA_Push_Get(MT_Dbl,nbf_fde*ipol,'eval',l_eval,k_eval(1)))
     &     call errquit(pname//' cannot allocate eval',0, MA_ERR)
      if (ipol .eq. 2)then
         k_eval(2) = k_eval(1) + nbf_fde
      endif

c     == allocate array for irreps ==
      if (.not.MA_Push_Get(mt_int,nbf_fde*ipol,'dft:irreps',l_ir,k_ir))
     &   call errquit(pname//' cannot allocate irreps',0, MA_ERR)
      call ifill(nbf_fde*ipol,1,int_mb(k_ir),1)

c     == alocate array for MO occupancy
      if (.not. MA_Push_Get(MT_Dbl,nbf_fde*ipol,'mo occ',l_occ,k_occ))
     &   call errquit(pname//' failed to alloc',999, MA_ERR)

c
c Lets get the number of electrons now in the FDE, for scaling
c
      if (.not. geom_nuc_charge(geom_fde, nucl_charge_fde))
     &     call errquit(pname//'geom_nuc_charge failed', 0,
     &       GEOM_ERR)

      nelectrons_fde = nint(nucl_charge_fde)
      even_fde=mod(nelectrons_fde,2).eq.0

c We should soft code this at some point.
      if (even_fde) then
        mult_fde = 1
      else
        mult_fde = 2
      endif
      nopen_fde = mult_fde - 1
      nclosed_fde = (nelectrons_fde-nopen_fde)/2

      nmo_fde(1) = nbf_fde
      noc_fde(1) = nclosed_fde
      if (ipol.eq.1) then
        scftype = 'RHF'
      else
        scftype = 'UHF'
c
        noc_fde(1) = nclosed_fde + nopen_fde
        noc_fde(2) = nclosed_fde       
        nmo_fde(2) = nbf_fde
      endif
c 
      if (rtdb_cget(rtdb, 'dft:frozemb:movecs', 1,
     &      frozemb_movecs)) then
        call util_file_name_resolve(frozemb_movecs, .false.)
        frozemb_movecs_out = frozemb_movecs
      else
        frozemb_movecs = 'atomic'
        call util_file_name('movecs', .false.,.false.,
     &                       frozemb_movecs_out)
        call util_file_name_resolve(frozemb_movecs_out, .false.)
      endif

c     == determine guess ==
      toll_s = 0.0d0

c      call dft_guessin(frozemb_movecs(1:inp_strlen(frozemb_movecs)),
      call dft_guessin('atomic',
     ,     ldmix,ncanorg,.false.,vecs_or_dens,
     ,     ipol,nbf_fde,g_movecs_fde,g_gmovecs_fde,
     ,     toll_s,dbl_mb(k_svals))
c
c Lets check if we have an input file, else we'll create from atoms.
c
c Construct density matrix

c Lifted from scf_dft/dft_scf.F
      call scf_vectors_guess(rtdb, tol2e_sleazy, geom_fde, fd_bas_han,
     &  basis_trans, frozemb_movecs(1:inp_strlen(frozemb_movecs)), 
     &  frozemb_movecs_out(1:inp_strlen(frozemb_movecs_out)),
     &  movecs_guess, scftype, nclosed_fde, nopen_fde,
     &  nbf_fde, nmo_fde, noc_fde(1), noc_fde(2), 
     &  k_eval, k_occ, k_ir, g_gmovecs_fde, g_dens_fde,
     &  vecs_or_dens, 'dft', title, oskel, oadapt,
     &  .true.)

C Is everything scaled automatically inside scf_vectors_guess for UHF?
      if (ipol.gt.1) then
c           call ga_copy(g_dens_fde(1),g_dens_fde(2))
c Load the extra electron in the alpha channel
c           call ga_dscal(g_dens_fde(1),
c     &                   dble(nelectrons_fde-
c     &                   nclosed_fde)/dble(nelectrons_fde))
c           call ga_dscal(g_dens_fde(2),dble(nclosed_fde)/
c     &                   dble(nelectrons_fde))
c           if(oskel) call ga_symmetrize(g_dens_fde(2))
c Recombine for now, as that is how we will use them in the ESP
           call ga_dadd(1d0, g_dens_fde(1), 1d0, g_dens_fde(2), 
     &                  g_dens_fde(1))
        call ga_sync()
      endif

      call dft_guessout(nmo_fde, nbf_fde, 
     &                  g_gmovecs_fde, g_movecs_fde, ipol)
c
c g_vxc_fde -- to be added to g_frozemb once calculated
      if (.not. ga_duplicate(g_frozemb, g_vxc_fde(1), 'V_FDE'))
     &  call errquit(pname//'Error creating g_vxc_fde1',0,GA_ERR)
      call ga_zero(g_vxc_fde(1))
c second spin channels are not used as electrostatic embedding is
c independent. This needs to be tidied up once checked through.
c       if (ipol.gt.1) then
c         if (.not. ga_duplicate(g_frozemb, g_vxc_fde(2), 'Vcoul,FDE(B)'))
c     &     call errquit(pname//'Error creating g_vxc_fde2',0,GA_ERR)
c         call ga_zero(g_vxc_fde(2))     
c       endif

      call grid_quadv0_fde(rtdb, g_dens_fde, g_vxc_fde,
     &                  nExc_fde, rho_n, Exc_fde,
     &                  dbl_mb(itmat),1,g_dens_fde)

      call ga_dadd(1d0, g_vxc_fde(1), 1d0, g_frozemb, g_frozemb)

      write(luout,*) '-> Electronic contributions complete'
c
c Deallocate arrays we aren't using
c
c      if (.not. ga_destroy(g_xcinv_fde)) call errquit
c     &   (pname//' could not destroy g_xcinv_fde', 0, GA_ERR)
      if (.not. ga_destroy(g_vxc_fde(1))) call errquit
     &   (pname//' could not destroy g_vxc_fde1', 0, GA_ERR)
      if (.not. ga_destroy(g_movecs_fde(1))) call errquit
     &   (pname//' could not destroy g_movecs_fde1', 0, GA_ERR)
c      if (.not. ga_destroy(g_dens_fde(1))) call errquit
c     &   (pname//' could not destroy g_dens_fde1', 0, GA_ERR)
      if (ipol.gt.1) then
c        if (.not. ga_destroy(g_vxc_fde(2))) call errquit
c     &   (pname//' could not destroy g_vxc_fde2', 0, GA_ERR)
        if (.not. ga_destroy(g_movecs_fde(2))) call errquit
     &   (pname//' could not destroy g_movecs_fde2', 0, GA_ERR)
c        if (.not. ga_destroy(g_dens_fde(2))) call errquit
c     &   (pname//' could not destroy g_dens_fde2', 0, GA_ERR)
      endif
      if (.not.ma_pop_stack(l_occ))
     & call errquit(pname//'Cannot pop stack l_occ', 12, MA_ERR)
      if (.not.ma_pop_stack(l_ir))
     & call errquit(pname//'Cannot pop stack l_ir', 12, MA_ERR)
      if (.not.ma_pop_stack(l_eval))
     & call errquit(pname//'Cannot pop stack l_eval', 12, MA_ERR)
      if (.not.ma_pop_stack(l_svals))
     & call errquit(pname//'Cannot pop stack l_svals', 12, MA_ERR)

C Just to test this doesn't alter the code at all in terms of parameters      
c      call ga_zero(g_frozemb)

      write(luout,*) '== Created Frozen Embedding Matrix =='
      write(luout,*) ''
      call util_flush(luout)

      return
      end
c Lifted and edited from grid_quadv0
      subroutine grid_quadv0_fde(rtdb, g_dens, g_vxc, nexc, rho_n,
     &                           Exc,tmat,fde_option,g_dens_fde)
c
      integer rtdb
      integer g_dens(2),nexc,g_vxc(*)
      integer g_dens_fde(2)
      double precision rho_n
      double precision Exc(2)
      double precision tmat(*)  ! XC fit
      integer fde_option
      integer calc_type
c fde_option: 1 = ESP
c             2 = QM with Exc FDE
c             3 = FDE with Exc FDE
c             4 = QM/MM with Exc FDE
c             5 = QM with Ts FDE
c             6 = FDE with Ts FDE
c             7 = QM/MM with Ts TDE
      calc_type = fde_option + 6
c Due to offset by 6 in subsequent subroutines

      if (fde_option.eq.1.or.fde_option.gt.4) then
c Set nExc to 1 as we are only collecting Ts/ESP
        call grid_quadv0_gen(rtdb, g_dens, g_vxc, 1, rho_n, Exc,
     &       1,calc_type,tmat,.false.,.false.,g_dens_fde)
      else
        call grid_quadv0_gen(rtdb, g_dens, g_vxc, nexc, rho_n, Exc,
     &       1,calc_type,tmat,.false.,.false.,g_dens_fde)
      endif
      return
      end
cc AJL/End
