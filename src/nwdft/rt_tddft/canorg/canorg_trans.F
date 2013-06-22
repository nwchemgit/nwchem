C
C     canorg_trans.F
C
C     Wrapper routine for transforming different matricies between
C     atomic orbital (AO) and molecular orbital (MO) basis using
C     canonical orthogonalization.
C     
C     Notation:  no prime means AO, prime means MO.
C
C     Overview of transforms:
C     
C     F' = X^+ F  X
C     S' = X^+ S  X    ( = 1_m)
C     P  = X   P' X^+
C     P' = Y^+ P  Y
C     F  = Y   F' Y^+
C
C     where X and Y are the complex-valued (put purely real) n x m
C     transformation matricies calculated in canorg_init().  All
C     conversions are for complex matricies only.  All matricies in AO
C     basis are ns_ao x ns_ao, while those in MO basis are ns_mo x ns_mo (where
C     ns_mo < ns_ao if we have linear dependencies).
C
C     g_a - input matrix
C     g_b - transformed matrix
C     mat - "F", "S", "P"
C     dir - "AO->MO", "MO->AO"
C
      subroutine canorg_trans (params, mat, dir, g_a, g_b)
      implicit none

#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_a               ! the input matrix
      character*1, intent(in)       :: mat
      character(*), intent(in)      :: dir


C     == Outputs ==
      integer, intent(in) :: g_b               ! transformed matrix returned in this GA


      
C     == Parameters ==
      character(*), parameter :: pname = "canorg_trans: "


C     == Variables ==
      integer dtypea, dtypeb
      integer n1a, n2a, n1b, n2b
      integer g_zs_nm         ! ns_ao x ns_mo scratch
      integer g_zs_mn         ! ns_mo x ns_ao scratch

      
C
C     Check GAs
C      
      call ga_check_handle (g_a,
     $     "fourth argument of "//pname//"not a valid GA")

      call ga_check_handle (g_b,
     $     "fifth argument of "//pname//"not a valid GA")

      call ga_inquire (g_a, dtypea, n1a, n2a)
      call ga_inquire (g_b, dtypeb, n1b, n2b)

      if (dtypea .ne. mt_dcpl)
     $     call errquit (pname//"input matrix must be complex", 0, 0)

      if (dtypeb .ne. mt_dcpl)
     $     call errquit (pname//"output matrix must be complex", 0, 0)


C
C     Allocate scratch space.
C
      if (.not. ga_create (mt_dcpl, params%ns_ao, params%ns_mo,
     $     "scr nxm", 0, 0, g_zs_nm))
     $     call errquit (pname//"alloc zscr_nm failed", 0, GA_ERR)

      if (.not. ga_create (mt_dcpl, params%ns_mo, params%ns_ao,
     $     "scr nxm", 0, 0, g_zs_mn))
     $     call errquit (pname//"alloc zscr_mn failed", 0, GA_ERR)

      

C
C     Check sizes of matricies and figure out transform routine to use.
C
      if (dir.eq."AO->MO") then ! g_a should be ns_ao x ns_ao, g_b ns_mo x ns_mo
         
         if ((n1a .ne. params%ns_ao).or.(n2a .ne. params%ns_ao))
     $        call errquit (pname//"input matrix in AO basis "//
     $        "should be ns_ao x ns_ao", 0, 0)
         
         if ((n1b .ne. params%ns_mo).or.(n2b .ne. params%ns_mo))
     $        call errquit (pname//"output matrix in MO basis "//
     $        "should be ns_mo x ns_mo", 0, 0)

         if ((mat .eq. "F").or.(mat .eq. "S")) then
            call canorg_trans_zfock_ao2mo (params, g_zs_nm, g_a, g_b)
         elseif (mat .eq. "P") then
            call canorg_trans_zdens_ao2mo (params, g_zs_nm, g_a, g_b)
         else
            call errquit(pname//'unknown matrix type: "'// mat//'"',0,0)
         endif
         

      elseif (dir.eq."MO->AO") then ! g_a: ns_mo x ns_mo, g_b: ns_ao x ns_ao
         
         if ((n1a .ne. params%ns_mo).or.(n2a .ne. params%ns_mo))
     $        call errquit (pname//"input matrix in MO basis "//
     $        "should be ns_mo x ns_mo", 0, 0)
         
         if ((n1b .ne. params%ns_ao).or.(n2b .ne. params%ns_ao))
     $        call errquit (pname//"output matrix in AO basis "//
     $        "should be ns_ao x ns_ao", 0, 0)

         if ((mat .eq. "F").or.(mat .eq. "S")) then
            call canorg_trans_zfock_mo2ao (params, g_zs_mn, g_a, g_b)
         elseif (mat .eq. "P") then
            call canorg_trans_zdens_mo2ao (params, g_zs_mn, g_a, g_b)
         else
            call errquit (pname//'unknown matrix type: "'//mat//'"',0,0)
         endif

      else
         call errquit (pname//'unknown transform: "'//dir//'"', 0, 0)
      endif


C
C     Clean up scratch.
C
      if (.not. ga_destroy (g_zs_nm))
     $     call errquit ("destroy scr_nm failed", 0, GA_ERR)

      if (.not. ga_destroy (g_zs_mn))
     $     call errquit ("destroy scr_mn failed", 0, GA_ERR)

      end subroutine

C====================================================================
C
C     Transform complex Fock matrix from AO to MO basis.  F' = X^+ F X.
C     Also applies for the overlap matrix.
C
      subroutine canorg_trans_zfock_ao2mo (params, g_zs_nm,
     $     g_zao, g_zmo)
      implicit none
      
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "rt_tddft.fh"


C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_zs_nm         ! ns_ao x ns_mo scratch
      integer, intent(in)           :: g_zao

      
C     == Outputs ==
      integer, intent(in) :: g_zmo


      call ga_zero (g_zs_nm)
      call ga_zero (g_zmo)

      call ga_zgemm ("N", "N", params%ns_ao, params%ns_mo, params%ns_ao,
     $     z1, g_zao, params%g_zX, z0, g_zs_nm)

      call ga_zgemm ("T", "N", params%ns_mo, params%ns_mo, params%ns_ao,
     $     z1, params%g_zX, g_zs_nm, z0, g_zmo)
      
      end subroutine


C====================================================================
C
C     Transform complex Fock matrix from MO to AO basis.  F = Y F' Y^+
C
      subroutine canorg_trans_zfock_mo2ao (params, g_zs_mn,
     $     g_zmo, g_zao)
      implicit none

#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "rt_tddft.fh"

C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_zs_mn       ! ns_mo x ns_ao scratch
      integer, intent(in)           :: g_zmo

C     == Outputs ==
      integer, intent(in) :: g_zao

      
C     == Parameters ==
      character(len=*), parameter ::
     $     pname = "canorg_trans_zfock_mo2ao: "


C     == Variables ==
      integer g_tmpao_re, g_tmpao_im


      call ga_zero (g_zs_mn)
      call ga_zero (g_zao)

      call ga_zgemm ("N", "T", params%ns_mo, params%ns_ao, params%ns_mo,
     $     z1, g_zmo, params%g_zY, z0, g_zs_mn)

      call ga_zgemm ("N", "N", params%ns_ao, params%ns_ao, params%ns_mo,
     $     z1, params%g_zY, g_zs_mn, z0, g_zao)

      end subroutine



C====================================================================
C
C     Transform complex Fock dens mat from AO to MO basis.
C     P' = Y^+ P Y
C     
      subroutine canorg_trans_zdens_ao2mo (params, g_zs_nm,
     $     g_zao, g_zmo)
      implicit none

#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "rt_tddft.fh"


C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_zs_nm         ! ns_ao x ns_mo scratch
      integer, intent(in)           :: g_zao

      
C     == Outputs ==
      integer, intent(in) :: g_zmo

C     == Parameters ==
      character(len=*), parameter ::
     $     pname = "canorg_trans_zdens_ao2mo: "


C     == Variables ==
      integer g_tmpmo_re, g_tmpmo_im


      call ga_zero (g_zs_nm)
      call ga_zero (g_zmo)

      call ga_zgemm ("N", "N", params%ns_ao, params%ns_mo, params%ns_ao,
     $     z1, g_zao, params%g_zY, z0, g_zs_nm)

      call ga_zgemm ("T", "N", params%ns_mo, params%ns_mo, params%ns_ao,
     $     z1, params%g_zY, g_zs_nm, z0, g_zmo)

      end subroutine



C====================================================================
C
C     Transform complex dens mat from MO to AO basis.  P = X P' X^+.
C     
      subroutine canorg_trans_zdens_mo2ao (params, g_zs_mn,
     $     g_zmo, g_zao)
      implicit none


#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "rt_tddft.fh"


C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_zs_mn       ! ns_mo x ns_ao scratch
      integer, intent(in)           :: g_zmo

      
C     == Outputs ==
      integer, intent(in) :: g_zao

      
C     == Parameters ==
      character(len=*), parameter ::
     $     pname = "canorg_trans_zdens_mo2ao: "


C     == Variables ==
      integer g_tmpao_re, g_tmpao_im


      call ga_zero (g_zs_mn)
      call ga_zero (g_zao)

      call ga_zgemm ("N", "T", params%ns_mo, params%ns_ao, params%ns_mo,
     $     z1, g_zmo, params%g_zX, z0, g_zs_mn)

      call ga_zgemm ("N", "N", params%ns_ao, params%ns_ao, params%ns_mo,
     $     z1, params%g_zX, g_zs_mn, z0, g_zao)


      end subroutine
c $Id$
