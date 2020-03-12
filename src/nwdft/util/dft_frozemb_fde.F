cc AJL/Begin/FDE

      subroutine dft_frozemb_fde_init(rtdb,test_sic,g_vxc,g_vxc_fde,
     &                                g_xcinv_fde)
c First subroutine to do all the checks and make sure we can run the
c calculation. In an ideal world we'd slowly reduce this list of
c limitations, but it is necessary for now to ensure the results are
c physically meaningful and the calculations will run efficiently.
      implicit none

#include "rtdb.fh"
#include "cdft.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"

      integer rtdb
      integer test_sic
      integer g_vxc(*)
      integer g_vxc_fde(*)
      integer g_xcinv_fde
c local
      integer i
      logical xcreplicated_l
      logical cam_exch
      character*32 pname
      pname = 'dft_frozemb_fde_init: '

c Incompatible errors. With testing and/or editing these could be removed
      if (lxdm.eq.1)
     & call errquit(pname//'XDM is not compatible with FDE',1,INPUT_ERR)
      if (test_sic.ne.0)
     & call errquit(pname//'SIC is not compatible with FDE',1,INPUT_ERR)
      if (cs00.or.lb94)
     & call errquit(pname//'asymptotic corrections are not
     &      compatible with FDE', 10, INPUT_ERR)
      if (rtdb_get(rtdb, 'dft:cam_exch', mt_log, 1, cam_exch)) then
        if (cam_exch)
     &     call errquit(pname//'range-separated functionals are not
     &     compatible with FDE', 10, INPUT_ERR)
        endif

c Disable XC replica so parallelisation works
c Firstly we'll check if XC replication is specifically turned on and
c make a note of this under dft:frozemb for our further use later
      if (rtdb_get(rtdb, 'dft:xcreplicated', mt_log, 1, 
     &  xcreplicated_l)) then
c Store old value for xcreplicated
        if (.not.rtdb_put(rtdb, 'dft:frozemb:xcreplicated', mt_log, 1, 
     &    xcreplicated_l)) 
     &    call errquit(pname//'rtdb_put failed', 0, RTDB_ERR) 
      endif
c And now actually disable XC replication
      if (.not.rtdb_put(rtdb, 'dft:xcreplicated', mt_log, 1, .false.))
     &  call errquit(pname//'rtdb_put failed', 0, RTDB_ERR)

c Create the FDE Vxc matrix, as needed for non-additive terms
      do i=1,ipol 
        if (.not. ga_duplicate(g_vxc(i), g_vxc_fde(i), 'FDE VXC'))
     &    call errquit(pname//'failed duplicate',g_vxc_fde(i),GA_ERR)
        call ga_zero(g_vxc_fde(i))
      enddo 

c Is this necessary? It might make more sense just to throw an error if
c XCFIT is being used, as we do not know what the outcome is.
      if (XCFIT)
     &   call errquit(pname//'XCFIT is not compatible with FDE',1,
     &   INPUT_ERR)
c Duplicate error for Charge Density Fitting. To be tested out.
      if (CDFIT)
     &   call errquit(pname//'CDFIT is not compatible with FDE',1,
     &   INPUT_ERR)
c We'll include the functionality for XCFIT just in case it works in
c the future, but for now I'd rather prevent it running in case there
c are problems rather than getting it to work without testing it.
c NOTE: If re-enabled, remember to also edit fde_end where gxcinv_fde is
c destroyed at the end of the calculation
c
c      if (XCFIT) then
c        if (.not. ga_duplicate(g_xcinv, g_xcinv_fde, 'VXCINV FDE'))
c     &    call errquit(pname//'failed duplicate',g_xcinv_fde, GA_ERR)
c        call ga_zero(g_xcinv_fde)
c      endif

      return

      end

      subroutine dft_frozemb_fde(rtdb, g_frozemb, g_dens, g_dens_fde) 
c Calcule the electrostatic frozen density embedding matrix, due to both
c the nuclear and electronic contributions. This is done sequentially as
c -> Work out nuclear contribution from FDE geometry
c -> Work out electronic from previous wavefunction or linear
c    combination of atomic orbitals - former is preferable
c These are then summed up and added in, unchanged to the Fock matrix 
c when constructed at each iteration of the self-consistent cycle therein
c
c Other things to note, and for the to do list:
c - currently this only deals with singlet or doublet electronic
c systems. For singlet the unpaired electron should be in the alpha
c channel. It would be nice to make this more dynamic but I'm not sure how
c this can be achieved as I'm not very familiar with the movecs structure
c - it also will not work with charged frozen density environments
      implicit none

#include "rtdb.fh"
#include "bas.fh"
#include "cdft.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "stdio.fh"

      integer rtdb
      integer g_frozemb

      character*32 pname
      character*255 basis_trans
      character*255 basis_name
      character*255 frozemb_movecs
      character*255 frozemb_movecs_out
      character*7 vecs_or_dens
      character*4 scftype

      integer nExc_fde
c      integer g_electronic_fde(4)
      integer g_electronic_fde
      integer itmat
      integer natoms

      integer nelectrons_fde
      integer nclosed_fde
      integer nopen_fde
      integer mult_fde
      logical even_fde

      integer g_dens(2)
      integer g_dens_fde(2)

      double precision nucl_charge_fde
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

      integer me
      me = ga_nodeid()

      Exc_fde(1) = 0d0
      Exc_fde(2) = 0d0
      rho_n = 0d0
      nExc_fde = idecomp + 1
      itmat = 0

      if (me.eq.0) then
        write(luout,*) ''
        write(luout,*) '== Preparing Frozen Density Embedding Matrix =='
      endif

      pname = 'dft_frozemb_fde: '
c
c Nuclear component; easy part
c
      call int_1e_fde(AO_bas_han, AO_bas_han, g_frozemb,
     &               'fde_pot', oskel, FD_bas_han)

      if (me.eq.0)
     &  write(luout,*) '(*) Nuclear contributions complete'
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
      if (.not. ga_duplicate(g_dens_fde(1),
     &                       g_movecs_fde(1), 'FDE Alpha Movecs'))
     &     call errquit(pname//'Error creating ga',0,GA_ERR)
c
      call ga_zero(g_dens_fde(1))
      call ga_zero(g_movecs_fde(1))
c
      if (ipol.gt.1) then
c Beta Density
        if(.not.ga_duplicate(g_dens_fde(1),
     &                       g_dens_fde(2),'FDE Beta Dens'))
     &     call errquit(pname//'ga_duplicate failed',1, GA_ERR)
        if (.not. ga_duplicate(g_dens_fde(2),
     &                       g_movecs_fde(2), 'FDE Beta Movecs'))
     &     call errquit(pname//'Error creating ga',0,GA_ERR)
c
        call ga_zero(g_dens_fde(2))
        call ga_zero(g_movecs_fde(2))
c
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

c We should soft code this so we can deal with different multiplicity
c and charge. Currently we are restricted to neutral systems with
c multiplicity of 1 or 2
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
C Doesn't look so, thus we do it manually here.
      if (ipol.gt.1) then
           call ga_copy(g_dens_fde(1),g_dens_fde(2))
c Load the extra electron in the alpha channel
#if 0
c rescaling disabled after fixing the call to scf_vec_guess
           call ga_dscal(g_dens_fde(1),
     &                   dble(nelectrons_fde-
     &                   nclosed_fde)/dble(nelectrons_fde))
           call ga_dscal(g_dens_fde(2),dble(nclosed_fde)/
     &                   dble(nelectrons_fde))
           if(oskel) call ga_symmetrize(g_dens_fde(2))
#endif           
c Recombine for now, as that is how we will use them in the ESP
           call ga_dadd(1d0, g_dens_fde(1), 1d0, g_dens_fde(2), 
     &                   g_dens_fde(1))
      endif

      call dft_guessout(nmo_fde, nbf_fde, 
     &                  g_gmovecs_fde, g_movecs_fde, ipol)
c
c g_vxc_fde -- to be added to g_frozemb once calculated
      if (.not. ga_duplicate(g_frozemb, g_electronic_fde, 'V_FDE'))
     &  call errquit(pname//'Error creating g_electronic_fde',0,GA_ERR)
      call ga_zero(g_electronic_fde)
c second spin channels are not used as electrostatic embedding is
c independent. This needs to be tidied up once checked through.
c       if (ipol.gt.1) then
c         if (.not. ga_duplicate(g_frozemb, g_electronic_fde(2), 'V_FDE(B)'))
c     &     call errquit(pname//'Error creating g_electronic_fde2',0,GA_ERR)
c         call ga_zero(g_electronic_fde(2))     
c       endif

      call grid_quadv0_fde(rtdb, g_dens, g_electronic_fde,
     &                  nExc_fde, rho_n, Exc_fde,
     &                  dbl_mb(itmat), 1, g_dens_fde)

      call ga_dadd(1d0, g_electronic_fde, 1d0, g_frozemb, g_frozemb)

      if (me.eq.0) 
     &  write(luout,*) '(*) Electronic contributions complete'
c
c Deallocate arrays 
c
      if (.not. ga_destroy(g_electronic_fde)) call errquit
     &   (pname//' could not destroy g_electronic_fde', 0, GA_ERR)
      if (.not. ga_destroy(g_movecs_fde(1))) call errquit
     &   (pname//' could not destroy g_movecs_fde1', 0, GA_ERR)
c
      if (ipol.gt.1) then
c        if (.not. ga_destroy(g_electronic_fde(2))) call errquit
c     &   (pname//' could not destroy g_electronic_fde2', 0, GA_ERR)
        if (.not. ga_destroy(g_movecs_fde(2))) call errquit
     &   (pname//' could not destroy g_movecs_fde2', 0, GA_ERR)
      endif
c
c Pop off arrays used for allocated the density of the system
c
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
c      call ga_print(g_frozemb)

      if (me.eq.0) then
        write(luout,*) '== Created Frozen Density Embedding Matrix =='
        write(luout,*) ''
        call util_flush(luout)
      endif

      return
      end

      subroutine grid_quadv0_fde(rtdb, g_dens, g_vxc, nexc, rho_n,
     &                           Exc,tmat,fde_option,g_dens_fde)

      implicit none
c Lifted and edited from grid_quadv0
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
c
      subroutine dft_frozemb_fde_esp(nqpts,qxyz,qwght,
     &           mbf_ao,chi_ao,g_dens_fde,g_vxc)
c Works out the potential for the FDE over the QM density grid
      implicit none

c Called from grid_quadv0b to:
c  -> Calculate potential at each point, embedding potential, 
c     from FDE_bas_han (pass the xyz to hnd_elf_con w/ nder =0)
c  -> Use acc_fock to accumulate this per basis (acc_fock, dft_frozemb)
c  -> Then put this into a ga (ga_put/ga_acc, dft_frozemb)
c In the long term it would be good if we could use the grid
c pre-screening for points with zero potential, and to use mbf != nbf,
c but for now I have used the simplest implementation possible

#include "cdft.fh"
#include "errquit.fh"
c#include "grid_quad.fh"
c#include "cgridfile.fh"
#include "stdio.fh"
c#include "bas.fh"
#include "mafdecls.fh"
c#include "global.fh"
c#include "util.fh"
c#include "grid_sic.fh"
c#include "dftps.fh"
c#include "geom.fh"

      integer nqpts
      double precision qxyz(3,nqpts), qwght(nqpts) ! Copied from grid_quad.fh
      integer mbf_ao
      double precision chi_ao(nqpts,mbf_ao)
      integer g_dens_fde(2),g_vxc(*) ! Copied from grid_quad.fh
c local
      integer l_epot, k_epot
      integer l_fock_fde, k_fock_fde
      integer l_ovlp, k_ovlp
      character*32 pname
      
      pname = 'dft_frozemb_fde_esp'

      if (.not. ma_push_get(mt_dbl,nqpts,'epot pnt',l_epot,k_epot))
     &     call errquit(pname//'ma failed',911,MA_ERR)
c      call dfill(nqpts,0,dbl_mb(k_epot),1)

      call hnd_elfcon_fde(FD_bas_han, geom_fde, g_dens_fde(1), qxyz,
     &                    nqpts, dbl_mb(k_epot), 0)
c Copied from dft_frozemb.F
      if (.not.MA_Push_Get(mt_dbl, mbf_ao*mbf_ao,
     &     'fock_fde', l_fock_fde, k_fock_fde))
     &     call errquit(pname//'fock_fde ma failed',911,MA_ERR)
      call dfill(mbf_ao*mbf_ao,0,dbl_mb(k_fock_fde),1)
      if (.not.MA_Push_Get(mt_dbl, mbf_ao*mbf_ao,
     &     'ovlp', l_ovlp, k_ovlp))
     &     call errquit(pname//'ovlp ma failed',911,MA_ERR)
      call dfill(mbf_ao*mbf_ao,0,dbl_mb(k_ovlp),1)
c       calculate fock and overlap
      call acc_fock(chi_ao, qwght, mbf_ao, nqpts,
     &              dbl_mb(k_epot), dbl_mb(k_fock_fde), dbl_mb(k_ovlp))
c       this is sequential, so lets accummumlate - not just put.
      call ga_acc(g_vxc(1),1,mbf_ao,1,mbf_ao,dbl_mb(k_fock_fde),
     &            mbf_ao, 1d0)

      if (.not.ma_chop_stack(l_epot))
     &     call errquit(pname//'cannot chop stack', 1, MA_ERR)


      return

      end
c
      subroutine dft_frozemb_fde_getvxc(rtdb, Exc_fde, nExc, iVxc_opt, 
     &                   g_xcinv_fde, g_dens, g_vxc_fde, IOLGC, 
     &                   fde_option, g_dens_fde)
c Designed to make the content in dft_fockbld a bit clearer
c This will manage all the pre- and post-processing when calculating
c the non-additive terms in xc energy calculations

      implicit none
#include "global.fh"
#include "cdft.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "util.fh"
      integer rtdb
      double precision Exc_fde(2)
      integer nExc
      integer iVxc_opt
      integer g_xcinv_fde
      integer g_dens
      integer g_vxc_fde(*)
      logical IOLGC
      integer fde_option
c fde_option: 1 = ESP
c             2 = QM with Exc FDE
c             3 = FDE with Exc FDE
c             4 = QM/MM with Exc FDE
c             5 = QM with Ts FDE
c             6 = FDE with Ts FDE
c             7 = QM/MM with Ts TDE
      integer g_dens_fde(*)
c local
      double precision rho_n_fde
      character*32 output_string

      integer me
      me = ga_nodeid()

      if (ipol.gt.1)
c     Restore alpha and beta densities.
     &  call ga_dadd(1d0, g_dens_fde(1), -1d0,
     &                    g_dens_fde(2), g_dens_fde(1))

      call ga_zero(g_vxc_fde(1))
      if (ipol.gt.1) call ga_zero(g_vxc_fde(2))

      rho_n_fde = 0.0d0
      Exc_fde(1) = 0.0d0
      Exc_fde(2) = 0.0d0

      call xc_getvxc(rtdb, Exc_fde, nExc, iVxc_opt, g_xcinv_fde,
     &               g_dens, g_vxc_fde, IOLGC, rho_n_fde,
     &               fde_option, g_dens_fde)

      if(util_print('frozemb nadd', print_debug)) then
        if (me.eq.0) then
          write(luout,*)
          if (fde_option.eq.2) then
            output_string = 'Exc(QM)'
          elseif (fde_option.eq.3) then
            output_string = 'Exc(FDE)'
          elseif (fde_option.eq.4) then
            output_string = 'Exc(Non-add)'
          elseif (fde_option.eq.5) then
            output_string = 'Ts(QM)'
          elseif (fde_option.eq.6) then
            output_string = 'Ts(FDE)'
          elseif (fde_option.eq.7) then
            output_string = 'Ts(Non-add)'
          else 
            call errquit('Invalid value for FDE_option?',1,INPUT_ERR)
          endif

          write(luout,*) output_string//' = ', Exc_fde
          write(luout,*)
        endif
c
        call ga_print(g_vxc_fde(1))
        if (ipol.gt.1)  call ga_print(g_vxc_fde(2))
c
        if (me.eq.0) write(luout,*)
        call util_flush(luout)
      endif

      if (ipol.gt.1)
c     Recombine the FDE densities
     &  call ga_dadd(1d0, g_dens_fde(1), 1d0,
     &                    g_dens_fde(2), g_dens_fde(1))

c Make sure we are all together.
c      call ga_sync()
      return

      end

      subroutine dft_frozemb_fde_end(rtdb,g_vxc_fde,g_dens_fde,
     &                               g_xcinv_fde)
c Subroutine to just deallocate all the FDE arrays when we are done.
c We'll also use this as an opportunity to remove the value of
c xcreplicated from the rtdb, which we put there to get parallelisation
c working in fde_init
      implicit none

#include "cdft.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"

      integer rtdb
      integer g_vxc_fde(*)
      integer g_dens_fde(*)
      integer g_xcinv_fde
c local
      integer i
      logical xcreplicated_l
      character*32 pname
      pname = 'dft_frozemb_fde_end: '

c We'll check if XC replication is stored and restore this
      if (rtdb_get(rtdb, 'dft:frozemb:xcreplicated', mt_log, 1, 
     &  xcreplicated_l)) then
c Restore old value for xcreplicated
        if (.not.rtdb_put(rtdb, 'dft:xcreplicated', mt_log, 1, 
     &    xcreplicated_l)) 
     &    call errquit(pname//'rtdb_put failed', 0, RTDB_ERR)
c Remove the stored value from the frozemb database
        if (.not.rtdb_delete(rtdb, 'dft:frozemb:xcreplicated')) 
     &     call errquit(pname//'rtdb_delete failed', 0, RTDB_ERR)
      else
c Remove any trace of xcreplication if we manually altered the setting
        if (.not.rtdb_delete(rtdb, 'dft:xcreplicated')) 
     &     call errquit(pname//'rtdb_delete failed', 0, RTDB_ERR)
      endif

      do i=1,ipol
        if (.not. ga_destroy(g_vxc_fde(i))) call errquit
     &     (pname//'could not destroy g_vxc_fde', 0, GA_ERR)
        if (.not. ga_destroy(g_dens_fde(i))) call errquit
     &     (pname//'could not destroy g_dens_fde', 0, GA_ERR)
      enddo

c This is currently disabled in fde_init. To be reinstated if XCFIT
c works correctly with FDE
c      if (XCFIT) then
c        if (.not. ga_destroy(g_xcinv_fde)) call errquit
c     &     (pname//'could not destroy g_xcinv_fde', 0, GA_ERR)
c      endif !xcfit
 
      return

      end

cc AJL/End
