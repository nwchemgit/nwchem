C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_angrad (n_comp_a,n_cont_a,n_comp_b,n_cont_b,
     &    angint,radint,ecp_ints)
*                                                                      *
*   Combine angular and radial parts to produce final integrals        *
*                                                                      *
*   Argument (status) - description                                    *
*                                                                      *
*   n_comp_a (inp) - number of components on centre A (cart or sph)    *
*   n_cont_a (inp) - number of contracted functions on centre A        *
*   n_comp_b (inp) - number of components on centre B (cart or sph)    *
*   n_cont_b (inp) - number of contracted functions on centre B        *
*   angint - angular integrals                                         *
*   radint - radial integrals                                          *
*   ecp_ints - final combined integrals                                *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
      integer n_comp_a,n_cont_a,n_comp_b,n_cont_b
      integer i_a,i_b,j_a,j_b
      double precision angint(n_comp_a,n_comp_b),
     &    radint(n_cont_a,n_cont_b),
     &    ecp_ints(n_comp_a,n_cont_a,n_comp_b,n_cont_b)
*
      do i_b = 1,n_cont_b
        do j_b = 1,n_comp_b
          do i_a = 1,n_cont_a
            do j_a = 1,n_comp_a
              ecp_ints(j_a,i_a,j_b,i_b) = ecp_ints(j_a,i_a,j_b,i_b)
     &            +radint(i_a,i_b)*angint(j_a,j_b)
            end do
          end do
        end do
      end do
*
      return
      end
