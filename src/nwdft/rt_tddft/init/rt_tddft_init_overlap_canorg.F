C
C     rt_tddft_init_overlap_canorg.F
C
C     Initialize overlap matrix and find eigenvalues, then initialize
C     canonical orthogonalization.
C
C
      subroutine rt_tddft_init_overlap_canorg (params)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "cdft.fh"
#include "rt_tddft.fh"


C     == In/out ==
      type(rt_params_t), intent(inout) ::  params  !geom params stored in here

      
C     == Parameters ==
      character(*), parameter :: pname =
     $     "rt_tddft_init_overlap_canorg: "

      
C     == External ==
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked

      
C     == Variables ==
      integer me
      integer g_s, g_svecs
      integer isvals, lsvals
      logical atmblk

C     (tmp)
      double precision toll_s
      integer n_dep
C      integer g_scr, g_s12, g_sn12



      me = ga_nodeid ()


C
C     Build overlap matrix (ripped from dft_main0d)
C     
      if (.not. rtdb_get(params%rtdb, 'dft:atomblock',mt_log,1,atmblk))
     &     atmblk=.true.

      if(atmblk) then
         g_s = ga_create_atom_blocked(geom, AO_bas_han, 'AO ovl')
      else
         if (.not. ga_create(mt_dbl, nbf_ao, nbf_ao, 'AO ovl',
     &        0, 0, g_s))
     &        call errquit(pname//'Error creating ga',0,GA_ERR)
      endif
      
      call ga_zero(g_s)
      
      if (.not.MA_Push_Get(MT_Dbl, nbf_ao, 'ovl eig vals', lsvals, 
     &   isvals))
     &   call errquit(pname//'Cannot allocate ovl eig vals', 0,
     &       MA_ERR)
      
      call dfill(nbf_ao, 0.0d0, dbl_mb(isvals), 1)
      
      call int_1e_ga(AO_bas_han, AO_bas_han, g_s, 'overlap', oskel)
      if (oskel) call sym_symmetrize(geom, AO_bas_han, .false., g_s)

      
C     
C     Diagonalize overlap (partially ripped from dft_main0d).  We just
C     use the stock diagonalizer.
C     
      if (.not. ga_duplicate(g_s, g_svecs, 'AO ovl eig vecs'))
     &     call errquit(pname//'Error creating ga',0,GA_ERR)
      call ga_zero(g_svecs)

CXXX  [KAL]: valrgrind picking up unintialized values when using parallel diag routine??      
C      call ga_diag_std(g_s, g_svecs, Dbl_MB(isvals))
      call ga_diag_std_seq (g_s, g_svecs, Dbl_MB(isvals))


C
C     Now that we have the overlap eigenvalues/vectors, initialize
C     canonical orthogonalization.  This will set the value of
C     params%nbf_mo.
C      
      call canorg_init (params, dbl_mb(isvals), g_svecs)


C
C     Store overlap for future use and destroy eigenvals/vecs.
C
      params%g_s = g_s


C     XXX TMP XXX
c$$$      if (.not. ga_duplicate (g_s, g_scr, "scr"))
c$$$     $     call errquit (pname//"failed to create scr",0,0)
c$$$
c$$$      if (.not. ga_duplicate (g_s, g_s12, "s12"))
c$$$     $     call errquit (pname//"failed to create s12",0,0)
c$$$
c$$$      if (.not. ga_duplicate (g_s, g_sn12, "sn12"))
c$$$     $     call errquit (pname//"failed to create sn12",0,0)

      
      call util_tolls (params%rtdb, .false., toll_s,
     $     n_dep, dbl_mb(isvals), params%nbf_ao)


C     S^{-1/2}
c$$$      call ga_zero (g_scr)
c$$$      call diis_bld12_so(toll_s, dbl_mb(isvals), g_svecs, g_sn12, 
c$$$     &     g_scr, nbf_ao, 2)

C     S^{1/2}
c$$$      call ga_zero (g_scr)
c$$$      call diis_bld12_so(toll_s, dbl_mb(isvals), g_svecs, g_s12, 
c$$$     &     g_scr, nbf_ao, 3)


C      params%g_s12 = g_s12
C      params%g_sn12 = g_sn12


      if (.not. ga_destroy (g_svecs))
     $     call errquit (pname//"Failed to destroy Svecs", 0, 0)

c$$$      if (.not. ga_destroy (g_scr))
c$$$     $     call errquit (pname//"Failed to destroy scr", 0, 0)

      if (.not.ma_pop_stack(lsvals))
     $     call errquit(pname//'cannot pop stack',0, MA_ERR)

      end subroutine








CXXX  [KAL]: EXPERIMENTAL SO ROUTINE

      subroutine rt_tddft_init_overlap_canorg_new (params)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "cdft.fh"
#include "rt_tddft.fh"


C     == In/out ==
      type(rt_params_t), intent(inout) ::  params  !geom params stored in here

      
C     == Parameters ==
      character(*), parameter :: pname =
     $     "rt_tddft_init_overlap_canorg: "

      
C     == External ==
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked

      
C     == Variables ==
      integer me
      integer g_s, g_svecs, g_s_super
      integer isvals, lsvals
      logical atmblk
      integer i


      me = ga_nodeid ()


C
C     Build overlap matrix (ripped from dft_main0d)
C     note "geom" (active geometry) is in common block
C     
      if (.not. rtdb_get(params%rtdb, 'dft:atomblock',mt_log,1,atmblk))
     &     atmblk=.true.

      if(atmblk) then
         g_s = ga_create_atom_blocked(geom, AO_bas_han, 'AO ovl')
      else
         if (.not. ga_create(mt_dbl, nbf_ao, nbf_ao, 'AO ovl',
     &        0, 0, g_s))
     &        call errquit(pname//'Error creating ga',0,GA_ERR)
      endif
      
      call ga_zero(g_s)
      
      call int_1e_ga(AO_bas_han, AO_bas_han, g_s, 'overlap', oskel)
      if (oskel) call sym_symmetrize(geom, AO_bas_han, .false., g_s)


C     (note 2*nbf_ao)
      if (.not.MA_Push_Get(MT_Dbl, 2*nbf_ao, 'ovl eig vals', lsvals, 
     &   isvals))
     &   call errquit(pname//'Cannot allocate ovl eig vals', 0,
     &       MA_ERR)
      
      call dfill(2*nbf_ao, 0.0d0, dbl_mb(isvals), 1)



C
C     Upconvert overlap matrix to supermatrix form:
C
C     [ S 0 ]
C     [ 0 S ]
C
C      call rt_tddft_so_upconvert (params, g_s)

      if (.not. ga_create(mt_dbl, 2*nbf_ao, 2*nbf_ao, 'AO ovl (super)',
     &     0, 0, g_s_super))
     &     call errquit(pname//'Error creating g_s_super',0,GA_ERR)
      
      call ga_zero (g_s_super)

      
C     (top left)
      call ga_dadd_patch (
     $     1d0, g_s, 1, nbf_ao, 1, nbf_ao,
     $     1d0, g_s_super, 1, nbf_ao, 1, nbf_ao,
     $     g_s_super, 1, nbf_ao, 1, nbf_ao)


C     (bottom right)
      call ga_dadd_patch (
     $     1d0, g_s, 1, nbf_ao, 1, nbf_ao,
     $     1d0, g_s_super, nbf_ao+1, 2*nbf_ao, nbf_ao+1, 2*nbf_ao,
     $     g_s_super, nbf_ao+1, 2*nbf_ao, nbf_ao+1, 2*nbf_ao)

      
      
C     
C     Diagonalize the super overlap (partially ripped from dft_main0d).
C     We just use the stock diagonalizer.
C     
      if (.not. ga_duplicate(g_s_super, g_svecs, 'AO ovl eig vecs'))
     &     call errquit(pname//'Error creating ga',0,GA_ERR)
      call ga_zero(g_svecs)

CXXX  [KAL]: valrgrind picking up unintialized values when using parallel diag routine??      
C      call ga_diag_std(g_s, g_svecs, Dbl_MB(isvals))
      call ga_diag_std_seq (g_s_super, g_svecs, Dbl_MB(isvals))


C
C     Now that we have the overlap eigenvalues/vectors, initialize
C     canonical orthogonalization.  This will set the value of
C     params%nbf_mo.
C      
C      call canorg_init (params, dbl_mb(isvals), g_svecs)
C      XXX MAKE NEW CANORG WHICH TAKES SUPER MATS AND DOESNT UPCONVERT

      call ga_print (g_svecs)

      do i = 1, 2*nbf_ao
         if (ga_nodeid().eq.0) then
            write (6,*) i, dbl_mb(isvals+i-1)
         endif
      enddo
      
      call rt_tddft_halt('rt_tddft_init_overlap_canorg.F')

C      call canorg_init_new (params, dbl_mb(isvals), g_svecs)


      
C
C     Store overlap for future use and destroy eigenvals/vecs.
C
      params%g_s = g_s_super

      if (.not. ga_destroy (g_svecs))
     $     call errquit (pname//"Failed to destroy Svecs", 0, 0)

      if (.not.ma_pop_stack(lsvals))
     &     call errquit(pname//'cannot pop stack',0, MA_ERR)

      end subroutine


c $Id$
