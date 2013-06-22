C
C     Projects a matrix in AO basis onto a particular geometry fragment
C     via mask (already created).  Can accept either double or double
C     complex GAs.
C     
      subroutine rt_tddft_geom_mask (params, g_full, ig, g_maskd)
      implicit none
      
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_full  ! handle for full mat to be projected (AO bas)
      integer, intent(in)           :: ig      ! index of geom to project onto
      integer, intent(in)           :: g_maskd ! handle for masked matrix (values set)

      
C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_geom_mask: "
      double precision, parameter :: thresh = 1d-8


C     == Variables ==
      integer dtype_full, n1_full, n2_full
      integer dtype_maskd, n1_maskd, n2_maskd
      integer g_fullre, g_fullim


C
C     Check input GA
C
      call ga_check_handle (g_full,
     $     "second argument of "//pname//"not a valid GA")
      
      call ga_inquire (g_full, dtype_full, n1_full, n2_full)

      if (n1_full .ne. n2_full)
     $     call errquit (pname//"n1 must equal n2 for g_full")

      if ( (dtype_full.ne.mt_dbl).and.(dtype_full.ne.mt_dcpl) )
     $     call errquit (pname//"g_full must be dbl or dcpl",0,0)


C
C     Check output GA
C      
      call ga_check_handle (g_maskd,
     $     "fourth argument of "//pname//"not a valid GA")
      
      call ga_inquire (g_maskd, dtype_maskd, n1_maskd, n2_maskd)

      if (n1_maskd .ne. n2_maskd)
     $     call errquit (pname//"n1 must equal n2 for g_maskd")

      if ( (dtype_maskd.ne.mt_dbl).and.(dtype_maskd.ne.mt_dcpl) )
     $     call errquit (pname//"g_maskd must be dbl or dcpl",0,0)
      

C
C     Check input and output compatible.
C
      if (n1_full.ne.n1_maskd)
     $     call errquit (pname//"g_full and g_maskd must be same size",
     $     0,0)
      
      if (dtype_full.ne.dtype_maskd)
     $     call errquit (pname//"g_full and g_maskd must be same type",
     $     0,0)



C
C     Check this is valid geom index
C
      if ( (ig.gt.params%ngeoms).or.(ig.lt.1) )
     $     call errquit (pname//"bad ig: ", ig, 0)



C
C     Now if the GA are double, we just elem mult by the mask.  If
C     complex, we make extract real and imag parts and mask
C     individually.
C     
      if (dtype_full .eq. mt_dbl) then !double
         call ga_zero (g_maskd)
         call ga_elem_multiply (g_full, params%geom(ig)%g_mask, g_maskd)
C         call ga_screen (g_maskd, thresh)
         
      elseif (dtype_full .eq. mt_dcpl) then !complex

         if (.not. ga_create (mt_dbl, n1_full, n1_full, "fullre",
     $        0, 0, g_fullre))
     $        call errquit (pname//"failed to create fullre",0, GA_ERR)
         
         if (.not. ga_create (mt_dbl, n1_full, n1_full, "fullim",
     $        0, 0, g_fullim))
     $        call errquit (pname//"failed to create fullim",0, GA_ERR)

         call convert_z2d (g_full, g_fullre, g_fullim)
         
         call ga_elem_multiply(g_fullre,params%geom(ig)%g_mask,g_fullre)
C         call ga_screen (g_fullre, thresh)

         call ga_elem_multiply(g_fullim,params%geom(ig)%g_mask,g_fullim)
C         call ga_screen (g_fullim, thresh)
 
         call convert_d2z (1d0, g_fullre, 1d0, g_fullim, g_maskd)
        
         if (.not. ga_destroy (g_fullre))
     $        call errquit (pname//"failed to destroy fullre", 0,GA_ERR)
         
         if (.not. ga_destroy (g_fullim))
     $        call errquit (pname//"failed to destroy fullim", 0,GA_ERR)

      else
         call errquit (pname//"bad data type for g_full",0,0)
      endif

      end subroutine
c $Id$
