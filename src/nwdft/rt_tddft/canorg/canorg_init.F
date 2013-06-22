C     
C     canorg_init.F
C
C     Initializes canonical orthogonalization routines--calculates
C     number of linearly dependent eigenvectors and constructs the
C     truncated transformation matricies for F -> F', P->P', and P'->P
C     (note: prime notation means MO basis, ^+ means conjg transpose).
C
C     The transformations are:
C
C     F' = X^+  F   X
C     P  = X    P'  X^+
C     P' = Y^+  P   Y
C     F  = Y    F'  Y^+
C
C
C     X   = U s^{-1/2},        (n x m) matrix
C     X^+ = s^{-1/2} U^+,      (m x n) matrix
C     Y   = X s = U s^{+1/2}   (n x m) matrix
C     Y^+ = s^{+1/2} U^+       (m x n) matrix
C
C
C     n = number of AO basis functions
C     m = number of MO basis functions, m < n if we have linear deps.
C
C     U is the (n x m) matrix with eigenvectors of S as columns
C     and s is the (m x m) diagonal matrix of eigenvalues.
C
C     For spin orbit calculations, this routine will still make the
C     nbf_ao x nbf_mo matrices; they will later be "upconverted" to
C     super matrices by rt_tddft_so_upconvert().
C
C
      subroutine canorg_init (params, svals, g_svecs)
      implicit none
      
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "geom.fh"
C#include "cdft.fh"
#include "rtdb.fh"
#include "rt_tddft.fh"
      

C     == In/out ==
      type(rt_params_t), intent(inout) :: params

      
C     == Inputs ==
      double precision, intent(in) :: svals(params%nbf_ao)  !eigenvalues of S
      integer, intent(in)          :: g_svecs              !eigenvectors of S


C     == Parameters ==
      character(*), parameter :: pname = "canorg_init: "
      logical, parameter      :: oprint_conv_details = .false.


C     == External routines ==
      integer  dft_lindep
      external dft_lindep
      

C     == Variables ==
      character(20) str         !character version of numdep
      integer me
      integer i, j
      integer istart
      double complex zval
      double precision toll_s
      integer n_dep
      integer nbf_mo
      integer numdep
      integer g_zX, g_zY
      integer g_xpad, g_zXpad
      integer g_zs

C     (junk)
      integer g_tmp1, g_tmp2
      double precision junkval
      


      me = ga_nodeid ()

      
C
C     Get overlap tolerance.
C
      
C      if (.not. rtdb_get(params%rtdb, 'dft:toll_s', MT_DBL, 1, toll_s))
C     &     call errquit(pname//'failed to get dft:toll_s from rtdb',
C     $     0, RTDB_ERR)

      call util_tolls (params%rtdb, .false., toll_s,
     $     n_dep, svals, params%nbf_ao)


C
C     Determine number of molecular orbitals (< number of AO if we have
C     linear dependency).
C
      nbf_mo = dft_lindep (params%nbf_ao, toll_s, svals,
     $     oprint_conv_details, luout)
      
      numdep = params%nbf_ao - nbf_mo


C
C     Allocate trucated transformation matricies, size nbf_ao x nbf_mo, which
C     will be carried through entire simulation.  Note that I only make
C     complex-valued ones.
C
      if (.not. ga_create(mt_dcpl, params%nbf_ao, nbf_mo,
     $     "X trans. mat.", 0, 0, g_zX))
     $     call errquit (pname//"Failed to allocate g_zX", 0, GA_ERR)

      if (.not. ga_create(mt_dcpl, params%nbf_ao, nbf_mo,
     $     "Y trans. mat.", 0, 0, g_zY))
     $     call errquit (pname//"Failed to allocate g_zY", 0, GA_ERR)


C
C     Allocate local scratch space.
C
      if (.not. ga_create(mt_dbl, params%nbf_ao, params%nbf_ao,
     $     "X padded", 0,0, g_xpad))
     $     call errquit (pname//"alloc failed for g_xpad", 0, GA_ERR)

      if (.not. ga_create(mt_dcpl, params%nbf_ao, params%nbf_ao,
     $     "Complex X padded", 0,0, g_zXpad))
     $     call errquit (pname//"alloc failed for g_zXpad", 0, GA_ERR)

      if (.not. ga_create(mt_dcpl, nbf_mo, nbf_mo, "zs", 0, 0, g_zs))
     $     call errquit (pname//"alloc failed for zs", 0, GA_ERR)

      
C
C     Build truncated transformation matricies.  We first build the
C     transformation matrix using the dft canorg routine, which returns
C     the a nbf_ao x nbf_ao matrix padded with zeros.  We then take the nbf_ao x
C     nbf_mo patch and store it as our X matrix.  Finally, we multiply by s
C     to get Y.
C

C
C     Build X.  dft_canorg() builds a real X, so we have to convert to
C     complex before extracting the patch.
C
      call ga_zero (g_xpad)

      call dft_canorg (nbf_mo, svals, g_svecs, g_xpad)
      call convert_d2z (1d0, g_xpad, 0d0, g_xpad, g_zXpad)

      call ga_copy_patch ("N", g_zXpad, 1, params%nbf_ao, 1, nbf_mo,
     $     g_zX, 1, params%nbf_ao, 1, nbf_mo)


C
C     Build Y = X s and Z = s X.  s is diagonal m x m matrix of
C     eigenvalues.  Note that we start at numdep + 1, i.e., start at the
C     "non-zero" evals.
C
C     Note: This all assumes eigenvalues are ordered low to high...
C
      call ga_zero (g_zs)
      istart = numdep + 1


c$$$      if (me.eq.0) then
c$$$         do i = 1, params%nbf_ao
c$$$            write (luout, *) i , svals(i)
c$$$         enddo
c$$$      endif
c$$$      call halt ()

      
      if (me.eq.0) then
         j = 0 
         do i = istart, params%nbf_ao
            
            j = j + 1
            if (j > nbf_mo) call errquit (pname//"bad j > nbf_mo", 0, 0)
            
            zval = dcmplx (svals(i), 0d0)
C            write (luout, *) svals(i), 1d0/svals(i) !XXX TMP
            call ga_put (g_zs, j, j, j, j, zval, 1)
            
         enddo
      endif
      call ga_sync ()

      call ga_zgemm ("N", "N", params%nbf_ao, nbf_mo, nbf_mo,
     $     z1, g_zX, g_zs, z0, g_zY)  !Y transform matrix

      
C
C     Clean up.
C
      if (.not. ga_destroy (g_xpad))
     $     call errquit (pname//"failed to destroy g_xpad", 0, GA_ERR)

      if (.not. ga_destroy (g_zXpad))
     $     call errquit (pname//"failed to destroy g_zXpad", 0, GA_ERR)

      if (.not. ga_destroy (g_zs))
     $     call errquit (pname//"failed to destroy g_zs", 0, GA_ERR)

      
C
C     Store parameters and GA handles in params struct.
C
      params%nbf_mo = nbf_mo
      params%g_zX = g_zX
      params%g_zY = g_zY


C      if (numdep > 0) then
C         write (str, "(i0)") numdep

C         call rt_tddft_print_notice (pname//"detected "//trim(str)//
C     $        " linear dependencies.")

      end subroutine


      
C====================================================================
C
C     Diagonalize complex overlap matrix
C
            subroutine canorg_diag_zs (params)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params


C     == Parameters ==
      character(*), parameter :: pname = "canorg_diag_zs: "

      
C     == Variables ==
      integer g_zs
      integer lzevals, izevals
      integer g_zevecs
      integer i
      integer me


      me = ga_nodeid ()


C
C     Allocation
C      
      if (.not. ga_create(mt_dcpl, params%nbf_ao, params%nbf_ao,
     $     "complex S", 0, 0, g_zs))
     $     call errquit (pname//"Failed to allocate g_zs", 0, GA_ERR)

      if (.not. ga_create(mt_dcpl, params%nbf_ao, params%nbf_ao,
     $     "zS evecs", 0, 0, g_zevecs))
     $     call errquit (pname//"Failed to allocate g_zevecs", 0,GA_ERR)


      if (.not.ma_push_get(mt_dcpl,params%nbf_ao,'zevals',
     $     lzevals,izevals))
     &     call errquit(pname//'cannot allocate zevals',0, MA_ERR)
      


      call ga_zero (g_zs)
      call convert_d2z (1d0, params%g_s, 0d0, params%g_s, g_zs)
      call zmat_diag(g_zs, g_zevecs, dcpl_mb(izevals), .false., .false.)

      
      if (me.eq.0) then
         do i = 1, params%nbf_ao
            write (luout, *) i , dcpl_mb(izevals+i-1)
         enddo
      endif


C
C     Clean up
C      
      if (.not. ga_destroy (g_zs))
     $     call errquit (pname//"Failed to destroy g_zs", 0, GA_ERR)

      if (.not. ga_destroy (g_zevecs))
     $     call errquit (pname//"Failed to destroy g_zevecs", 0, GA_ERR)

      if (.not.ma_pop_stack(lzevals))
     $     call errquit (pname//"Failed to pop stack", 0, GA_ERR)

      call errquit ("XXX canorg_diag_zs", 0, 0)

      end subroutine
c $Id$
