c
c     Initialize and Finalize integrals
c     HvD, NG: 2012
c
      subroutine tddft_grad_init_ints(rtdb,ihdl_bfao,ihdl_geom,
     +                                      tol2e)
c
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "stdio.fh"
c
c     Input
c
      integer rtdb      ! runtime database handle
      integer ihdl_bfao ! AO basis set handle
      integer ihdl_geom ! geometry handle
      double precision tol2e ! 2-electron integral tolerance
c
      integer calc_type ! calculation type for Fock builder
c
      character*32 pname
      pname = "tddft_grad_init_ints: "
c
c     Code
      if (.not.rtdb_get(rtdb,'fock_xc:calc_type',mt_int,1,calc_type)) 
     +  calc_type=0
      if (.not.rtdb_put(rtdb,'fock_xc:calc_type',mt_int,1,2))
     +  call errquit(pname//'failed to set calc_type 2',0, RTDB_ERR)
      if (.not.rtdb_put(rtdb,'fock_xc:calc_type_save',mt_int,1,
     +    calc_type)) 
     +  call errquit(pname//'failed to set calc_type_save',0,RTDB_ERR)
c
      call int_init(rtdb,1,ihdl_bfao)
      call schwarz_init(ihdl_geom,ihdl_bfao)
      call scf_get_fock_param(rtdb,tol2e)
c
      end
c
      subroutine tddft_grad_finalize_ints(rtdb)
c
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"
c
c     Input
c
      integer rtdb      ! runtime database handle
      integer calc_type ! calculation type for Fock builder
c
      character*32 pname
      pname = "tddft_grad_finalize_ints: "
c
c     Code
c
      if (.not.rtdb_get(rtdb,'fock_xc:calc_type_save',mt_int,1,
     +  calc_type)) 
     + call errquit(pname//'failed to get calc_type_save',0,RTDB_ERR)
      if (.not.rtdb_put(rtdb,'fock_xc:calc_type',mt_int,1,calc_type))
     + call errquit(pname//'failed to reset calc_type',0,RTDB_ERR)
      if (.not.rtdb_delete(rtdb,'fock_xc:calc_type_save')) 
     + call errquit(pname//'failed to delete calc_type_save',0,RTDB_ERR)
c
      call fock_2e_tidy(rtdb)
      call int_terminate()
      call schwarz_tidy()
c
      end
c $Id$
