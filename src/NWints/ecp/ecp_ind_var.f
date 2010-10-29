C $Id$
************************************************************************
*                                                                      *
      subroutine ecp_ind_var (l_c,n_blk,n_co1,n_co2,i_co1,i_co2,
     &    i_ze1,i_ze2,n_x,n_co_tot,n_co_max,i_off,n_pass,i_cont_c,
     &    n_cont_c,skip)
*                                                                      *
*   Define various index and length parameters for the selected set of *
*   integrals for a given ECP projector. These depend on the existence *
*   of the ECP parameters and the relationships between spin-free and  *
*   spin-orbit potentials.                                             *
*                                                                      *
*   Argument (status) - description                                    *
*                                                                      *
*   l_c (inp) - angular momentum of ECP projector                      *
*   n_blk (inp) - number of blocks, 1 for SF, 3 for SO, 4 for both     *
*   n_co1 (inp) - number of SF ECP coefficients                        *
*   n_co2 (inp) - number of SO ECP coefficients                        *
*   i_co1 (inp) - index of first SF ECP coefficient                    *
*   i_co2 (inp) - index of first SO ECP coefficient                    *
*   i_ze1 (inp) - index of first SF ECP exponent                       *
*   i_ze2 (inp) - index of first SO ECP exponent                       *
*   n_x (out) - actual number of blocks which can be calculated        *
*   n_co_tot (out) - total number of coefficients for one pass of the  *
*                    radial integrals                                  *
*   n_co_max (out) - maximum length of ECP contraction                 *
*   i_off (out) - offset of actual data in the final index of the ECP  *
*                 integral array                                       * 
*   n_pass (out) - number of passes of the radial integrals            *
*   i_cont_c (out) - address of first ECP contraction                  *
*   n_cont_c (out) - number of ECP contractions                        *
*   skip (out) - logical to skip if no integrals would be generated    *
*                                                                      *
*   Written by K. G. Dyall                                             *
*                                                                      *
************************************************************************
      implicit none
      integer l_c,n_blk,n_co1,n_co2,i_co1,i_co2,i_ze1,i_ze2,
     &    n_x,n_co_tot,n_co_max,i_off,n_pass,i_cont_c,n_cont_c
      logical skip
*
*   Define parameters for the different integral class cases
*
      n_co_tot = n_co1
      n_co_max = n_co1
      i_cont_c = 1
      n_cont_c = 1
      i_off = 0
      n_pass = 1
      if (n_blk .eq. 4) then
*
*     Both spin-free and spin-orbit integrals requested.
*
        if ((i_co1 .eq. 0) .or. (n_co1 .eq. 0)) then
*
*       There is no spin-free ECP: if there is also no spin-orbit ECP
*       then exit loop, otherwise set parameters for spin-orbit only.
*
          skip = (i_co2 .eq. 0) .or. (n_co2 .eq. 0)
          n_x = 3
          i_cont_c = 2
          i_off = 1
        else if ((i_co2 .eq. 0) .or. (n_co2 .eq. 0)) then
*
*       There is a spin-free ECP but there is no spin-orbit ECP:
*       set parameters for spin-free only.
*
          n_x = 1
        else
*
*       There are both spin-free and spin-orbit ECPs
*
          n_x = 4
          if ((i_ze2 .eq. i_ze1) .and. (n_co2 .eq. n_co1)) then
*
*         SF and SO ECPs have the same exponents: do together
*
            n_cont_c = 2
            n_co_tot = n_co1*2
          else
*
*         SF and SO ECPs have different exponents: do two passes
*
            n_co_tot = max(n_co1,n_co2)
            n_co_max = n_co_tot
            n_pass = 2
          end if
        end if
      else
*
*     Either spin-free or spin-orbit integrals requested.
*
        if (n_blk .eq. 3) then
          skip = (l_c .eq. 0) .or.
     &        (i_co2 .eq. 0) .or. (n_co2 .eq. 0)
          i_cont_c = 2
          n_co_tot = n_co2
          n_co_max = n_co2
        else
          skip = (i_co1 .eq. 0) .or. (n_co1 .eq. 0)
        end if
        n_x = n_blk
      end if
*
      return
      end
