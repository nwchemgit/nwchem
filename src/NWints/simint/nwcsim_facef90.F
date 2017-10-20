#ifdef USE_SIMINT
      module nwcsim_shell
      use SimintFortran
#include "nwc_const.fh"
      type(c_simint_shell), target :: smnt_sh(nw_max_shells)
      integer :: nwcsim_noshell
      end module  nwcsim_shell
      subroutine nwcsim_init(rtdb,nbas,bases)
      use nwcsim_shell
      implicit none
#include "bas.fh"
#include "geom.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
#include "apiP.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
c      use iso_c_binding
      integer rtdb,nbas,bases(3)
c      type(c_simint_multi_shellpair), target :: msh
      double precision :: alpha(3), coef(3)
      integer :: ibasis,basis,bas
      logical cart_2_sphe
      integer :: geom, natoms
      integer :: num_der, max_ang
      integer :: mxmx
      integer :: iat,icset,ic1,ic2, nprim
      integer :: l, ncontr, isphere
      character*16 tag 
      integer l_coef,k_coef,l_exp,k_exp
      double precision coord(3)
      double precision q
      character*14 pname
      pname="nwcsim_f90init"
c     assuming cartesian for now
c     stick to nbas=1 for now
         if(nbas.gt.1) call errquit(
     C        ' simint interface not ready for multiple basis',0,0)
c     init
      call simint_init()
c
c     create shells simint structure
c
      do ibasis=1,nbas
         basis = bases(ibasis) 
         bas = basis + BASIS_HANDLE_OFFSET
         cart_2_sphe=bas_spherical(bas)
#if 0
         if(cart_2_sphe) call errquit(
     C        ' simint interface not ready for spherical ',0,0)
#endif
         if (.not.bas_nprim_cn_max(basis, mxmx))
     C        call errquit(' bas_nprim_cn_max failed ',0,BASIS_ERR)
         if (.not. MA_Push_Get(MT_Dbl,mxmx,'simexp',l_exp,k_exp))
     C        call errquit(' simint wrk ',mxmx,MA_ERR)
         if (.not. MA_Push_Get(MT_Dbl,mxmx,'simcoef',l_coef,k_coef))
     C        call errquit(' simint wrk ',mxmx, MA_ERR)
         geom  = ibs_geom(bas)
c
         natoms   =  ncenter(geom)
         nwcsim_noshell=0
         max_ang=0
         do iat=1,natoms
            if (.not.bas_ce2cnr(basis,iat,ic1,ic2))
     &           call errquit(pname//'Exiting ',11, BASIS_ERR)
            do icset = ic1,ic2
               if (.not.bas_continfo(basis, icset,
     &              l, nprim, ncontr, isphere))
     &              call errquit(pname//'Exiting ',3, BASIS_ERR)
               max_ang=max(l,max_ang)
               if (.not.bas_get_exponent(basis, icset, dbl_mb(k_exp)))
     &              call errquit(pname//'Exiting ',7, BASIS_ERR)
               if(.not. bas_get_coeff(basis,icset,dbl_mb(k_coef)))
     &              call errquit(pname//'Exiting ',8, BASIS_ERR)
               if (.not. geom_cent_get(geom, iat, tag,
     &              coord, q))call errquit
               nwcsim_noshell=nwcsim_noshell+1
               call simint_initialize_shell(smnt_sh(nwcsim_noshell))
               
               call simint_create_shell(nprim, l , 
     C              coord(1), coord(2), coord(3),
     &              dbl_mb(k_exp), dbl_mb(k_coef), 
     L              smnt_sh(nwcsim_noshell))
            enddo
         enddo
         if(.not.ma_chop_stack(l_exp)) call errquit(
     E        ' pop stack failed ',l_exp,0)
      enddo
c
c     normalization not needed since nwchem and simint use the same
cedo      call simint_normalize_shells(nwcsim_noshell,smnt_sh)
c
c     memory allocation
c
c     0 deriv
      num_der=0

      mem_2e4c = max(mem_2e4c,
     S     simint_eri_worksize(num_der, max_ang))
      memb_2e4c = mem_2e4c + mem_2e4c/10  ! +10% to be safe
      if(ga_nodeid().eq.0) then
         write(luout,*)
         write(luout,*) ' Using Simint Integral package'
         write(luout,*) ' simint: mem_2e4c   ',mem_2e4c
         write(luout,*) ' simint: memb_2e4c  ',memb_2e4c
         write(luout,*)
         call util_flush(luout)
      endif

      return
      end
      subroutine nwcsim_terminate()
      use nwcsim_shell
      implicit none
c
      integer icsh
c
      do icsh=1,nwcsim_noshell
         call simint_free_shell(smnt_sh(icsh))
      enddo

      call simint_finalize()
      return
      end
#else 
c
c     stubs when simint is not used
c
      subroutine nwcsim_init(rtdb,nbas,bases)
      implicit none
#include "errquit.fh"
      integer rtdb,nbas,bases
      return
      end
      subroutine nwcsim_terminate()
      implicit none
#include "errquit.fh"
      return
      end
#endif
      
