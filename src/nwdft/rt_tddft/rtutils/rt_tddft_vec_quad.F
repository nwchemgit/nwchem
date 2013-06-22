C
C     rt_tddft_vec_quad.F
C     
C     Routines for dealing with vector and quadrupole data types.
C
      
      subroutine rt_tddft_vec_zero (vec)
      implicit none
#include "rt_tddft.fh"

C     == In/out ==
      type(rt_vector_t) vec

      vec%x = 0d0
      vec%y = 0d0
      vec%z = 0d0
      end subroutine

      
C====================================================================
      subroutine rt_tddft_quad_zero (quad)
      implicit none
#include "rt_tddft.fh"

C     == In/out ==
      type(rt_quad_t) quad

      quad%xx = 0d0
      quad%xy = 0d0
      quad%xz = 0d0
      quad%yy = 0d0
      quad%yz = 0d0
      quad%zz = 0d0
      
      end subroutine



C====================================================================
      subroutine rt_tddft_vec_add (vec1, vec2, sum)
      implicit none
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_vector_t), intent(in) :: vec1, vec2
      

C     == Outputs ==
      type(rt_vector_t), intent(out) :: sum


      sum%x = vec1%x + vec2%x
      sum%y = vec1%y + vec2%y
      sum%z = vec1%z + vec2%z
      
      end subroutine
c $Id$
