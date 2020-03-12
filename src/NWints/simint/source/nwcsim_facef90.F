#ifdef USE_SIMINT
      subroutine nwcsim_init(rtdb,nbas,bases,num_der)
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
#include "rtdb.fh"
#include "stdio.fh"
#include "nwcsim.fh"
c      use iso_c_binding
      integer rtdb,nbas,bases(3)
      integer :: num_der
c      type(c_simint_multi_shellpair), target :: msh
      double precision :: alpha(3), coef(3)
      integer :: ibasis,basis,bas
      logical cart_2_sphe
      integer :: geom, natoms
      integer :: max_der, max_ang
      integer :: mxmx
      integer :: iat,icset,ic1,ic2, nprim
      integer :: l, ncontr, isphere
      character*16 tag 
      integer l_coef,k_coef,l_exp,k_exp
      double precision coord(3)
      double precision q
      integer nwcsim_maxam
      external nwcsim_maxam
      character*14 pname
      pname="nwcsim_f90init"
c     fetch stuff from rtdb
      if (.not.rtdb_get(rtdb,'int:smnt_screen_tol',MT_DBL,1,
     s     smnt_screen_tol))  smnt_screen_tol=1d-22
#define SIMINT_SCREEN_NONE         0
#define SIMINT_SCREEN_SCHWARZ      1
#define SIMINT_SCREEN_FASTSCHWARZ  2
      if (.not.rtdb_get(rtdb,'int:smnt_screen_method',MT_INT,1,
     s     smnt_screen_method))  smnt_screen_method=
     S     SIMINT_SCREEN_FASTSCHWARZ
c reset screen_tol when screening is off
      if(smnt_screen_method.eq.SIMINT_SCREEN_NONE)smnt_screen_tol=0d0

c     stick to nbas>=2 for now
         if(nbas.gt.2) call errquit(
     C        ' simint interface not ready for no basis gt 2',0,0)
c     init
      max_der=0
#ifdef SIMINT_GRADIENT
      max_der=1
#endif
      if(num_der.gt.max_der) call errquit(
     C        ' simint interface not ready for derivative ',
     D num_der,0)
      if(.not.nwcsim_initialized) then
      call simint_init()
c     needed for 1-e integrals
      call igamma_init()

      nwcsim_initialized=.true.
c
c     create shells simint structure
c
      nwcsim_nbas=nbas
      do ibasis=1,nbas
         basis = bases(ibasis) 
         bas = basis + BASIS_HANDLE_OFFSET
         nwcsim_bas(ibasis)=bas
c     offset=-565
         cart_2_sphe=bas_spherical(bas)
         if (.not.bas_nprim_cn_max(basis, mxmx))
     C        call errquit(' bas_nprim_cn_max failed ',0,BASIS_ERR)
         if (.not. MA_Push_Get(MT_Dbl,mxmx,'simexp',l_exp,k_exp))
     C        call errquit(' simint wrk ',mxmx,MA_ERR)
         if (.not. MA_Push_Get(MT_Dbl,mxmx,'simcoef',l_coef,k_coef))
     C        call errquit(' simint wrk ',mxmx, MA_ERR)
         geom  = ibs_geom(bas)
c
         natoms   =  ncenter(geom)
         nwcsim_noshell(bas)=0
         max_ang=0
         do iat=1,natoms
            if (.not.bas_ce2cnr(basis,iat,ic1,ic2))
     &           call errquit(pname//'Exiting ',11, BASIS_ERR)
            do icset = ic1,ic2
               if (.not.bas_continfo(basis, icset,
     &              l, nprim, ncontr, isphere))
     &              call errquit(pname//'Exiting ',3, BASIS_ERR)
               if(l.gt.nwcsim_maxam()) call errquit(
     p              pname//' simint library maxam too small ',
     p              nwcsim_maxam(),BASIS_ERR)
               max_ang=max(l,max_ang)
               if (.not.bas_get_exponent(basis, icset, dbl_mb(k_exp)))
     &              call errquit(pname//'Exiting ',7, BASIS_ERR)
               if(.not. bas_get_coeff(basis,icset,dbl_mb(k_coef)))
     &              call errquit(pname//'Exiting ',8, BASIS_ERR)
               if (.not. geom_cent_get(geom, iat, tag,
     &              coord, q))call errquit
               nwcsim_noshell(bas)=nwcsim_noshell(bas)+1
               call simint_initialize_shell(
     S              smnt_sh(nwcsim_noshell(bas),bas))
               
               call simint_create_shell(nprim, l , 
     C              coord(1), coord(2), coord(3),
     &              dbl_mb(k_exp), dbl_mb(k_coef), 
     L              smnt_sh(nwcsim_noshell(bas),bas))
c     dummy shell for 3c-
               dbl_mb(k_exp)=0d0
               dbl_mb(k_coef)=1d0
               call simint_create_shell(1, 0 , 
     C              coord(1), coord(2), coord(3),
     &              dbl_mb(k_exp), dbl_mb(k_coef), 
     L              zero_sh(nwcsim_noshell(bas),bas))
            enddo
         enddo
         if(.not.ma_chop_stack(l_exp)) call errquit(
     E        ' pop stack failed ',l_exp,0)
c
c     memory allocation
c
         isz_2e4c = max(isz_2e4c,
     S        simint_eri_worksize(num_der, max_ang))
         mem_2e4c = max(mem_2e4c,
     S        simint_eri_workmem(num_der, max_ang))
      enddo ! basis loop
      endif
c
c     normalization not needed since nwchem and simint use the same
cedo      call simint_normalize_shells(nwcsim_noshell,smnt_sh)
c     SIMINT_PRIM_SCREEN_STAT needs 4 more doubles
c      isz_2e4c        = isz_2e4c + 4
c      
      call util_align(isz_2e4c,SIMINT_SIMD_LEN)
      call util_align(mem_2e4c,SIMINT_SIMD_LEN)
c
      iszb_2e4c=isz_2e4c
      
      if(num_der.eq.1) then
         memb_2e4c = mem_2e4c + mem_2e4c/5
      else
         memb_2e4c = mem_2e4c + mem_2e4c/10 ! +10% to be safe
      endif
      call util_align(memb_2e4c,SIMINT_SIMD_LEN)
      mem_2e3c = mem_2e4c
      mem_2e2c = mem_2e4c
      if(ga_nodeid().eq.0) then
         write(luout,*)
         write(luout,*) ' Using Simint Integral package'
         write(luout,*) ' num_der',num_der
         write(luout,*) ' simint: mem_2e4c    ',mem_2e4c
         write(luout,*) ' simint: memb_2e4c   ',memb_2e4c
         write(luout,*) ' simint: isz_2e4c   ',isz_2e4c
         write(luout,*) ' simint: iszb_2e4c  ',iszb_2e4c
         write(luout,*) ' screen_method ',smnt_screen_method
         write(luout,*) ' screen_tol ',smnt_screen_tol
         write(luout,*)
         call util_flush(luout)
      endif

      return
      end
      subroutine nwcsim_terminate()
      use nwcsim_shell
      implicit none
#include "errquit.fh"
      integer nbas     ! [in]
      integer bases(3) ! [in]
c
      integer icsh,ibasis,bas
      character*14 pname
      pname="nwcsim_termina"
c
      if(nwcsim_initialized) then
      do ibasis=1,nwcsim_nbas
         bas=nwcsim_bas(ibasis)
            do icsh=1,nwcsim_noshell(bas)
               call simint_free_shell(smnt_sh(icsh,bas))
               call simint_free_shell(zero_sh(icsh,bas))
            enddo
      enddo
      endif
      call simint_finalize()
      nwcsim_initialized=.false.
      return
      end
      integer function nwcsim_maxam()
      implicit none
c
      nwcsim_maxam=SIMINT_MAXAM
      return
      end
      logical function nwcsim_initdone()
      use nwcsim_shell
      implicit none
c
      nwcsim_initdone=nwcsim_initialized
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
      subroutine nwcsim_terminate(nbas,bases)
      implicit none
#include "errquit.fh"
      integer nbas,bases
      return
      end
      logical function  nwcsim_initdone()
      implicit none
c
      nwcsim_initdone=.false.
      return
      end
      integer function nwcsim_maxam()
      implicit none
c
      nwcsim_maxam=-1
      return
      end
#endif
      logical function int_forcenwints(rtdb,rtdb_errmsg)
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
      integer rtdb
      character*(*) rtdb_errmsg
c
      logical out1
c
      int_forcenwints=.false.
c     Force texas integrals to false for CAM calculations
      if (.not.rtdb_put(rtdb,'int:cando_txs',mt_log,1,.false.))
     &     then
         rtdb_errmsg='cannot set int:cando_txs'
         return
      endif
c     unset cando_nw if defined to avoid Simint
      if (rtdb_get(rtdb,'int:cando_nw',mt_log,1,out1)) then
         if (.not.rtdb_delete(rtdb,'int:cando_nw')) then
            rtdb_errmsg='cannot del int:cando_nw'
            return
         endif
         if(ga_nodeid().eq.0) write(luout,*) ' cando_nw unset'
      endif
      int_forcenwints=.true.
      call ga_sync()
      return
      end
      logical function int_disablesim(rtdb,rtdb_errmsg)
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
      integer rtdb
      character*(*) rtdb_errmsg
c
      logical cando_txs,cando_nw
c
      int_disablesim=.false.
c
c     check if simint is enabled:
c     cando_nw=f & cando_txs=f
c     not needed if Simint is not enabled
c      
#ifdef USE_SIMINT
      if (.not.rtdb_get(rtdb,'int:cando_txs',mt_log,1,cando_txs))
     &  cando_txs=.true.
      if (.not.rtdb_get(rtdb,'int:cando_nw',mt_log,1,cando_txs))
     &  cando_nw=.true.
c     unset cando_nw if defined to avoid Simint
      if((.not.cando_txs).and.(.not.cando_nw)) then
         if (.not.rtdb_delete(rtdb,'int:cando_nw')) then
            rtdb_errmsg='cannot del int:cando_nw'
            return
         endif
         if(ga_nodeid().eq.0) write(luout,*) ' cando_nw unset'
         if (.not.rtdb_delete(rtdb,'int:cando_txs')) then
            rtdb_errmsg='cannot del int:cando_txs'
            return
         endif   
         if(ga_nodeid().eq.0) write(luout,*) ' cando_txs unset'
      endif
#endif      
      int_disablesim=.true.
      call ga_sync()
      return
      end
      
