C     
C     Computes -D.E (dipole coupling to a uniform external applied
C     electric field), which can be added to the Fock matrix in the AO
C     basis.
C
C     The input matrix g_ddote will be replaced with -D.E.
C
      subroutine dft_fock_efield (rtdb, efield, g_ddote, lscale)

      implicit none
      
#include "bas.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "cdft.fh"
#include "global.fh"
#include "util.fh"
#include "errquit.fh"


C     == Inputs ==
      integer, intent(in)            :: rtdb
      double precision, intent(in)   :: efield(3) !x,y,z components of electric field
      integer, intent(in)            :: g_ddote   !-D.E--value will be computed and stored in this GA
      logical, intent(in)            :: lscale    ! if .true., scale -D.E by norm of the dipole matrices

C     == Parameters ==
      character(*), parameter        :: pname = "dft_fock_efield: "

C     == Variables ==
      integer                        :: me
      integer                        :: g_dipx, g_dipy, g_dipz
      integer                        :: g_mask
      integer                        :: dtype, n1, n2
      integer                        :: ibf, cont, icont
      integer                        :: bf_start, bf_end
      integer                        :: iexp, lexp
      integer                        :: ival
      integer                        :: ncont
      integer                        :: l, nprim, ncontr, isphere
      double precision               :: norm, scale, exp_cutoff
      logical                        :: cut_this_cont


      me = ga_nodeid ()

C
C     Checks
C
      call ga_check_handle (g_ddote,
     $     pname//"third argument not a valid GA")
      
      call ga_inquire (g_ddote, dtype, n1, n2)

      if (dtype.ne.mt_dbl)
     $     call errquit (pname//"only valid for real GA",0,0)

      if (n1.ne.n2)
     $     call errquit (pname//"only valid for square GAs",0,0)

      if (n1.ne.nbf)
     $     call errquit (pname//"GA must be nbf x nbf",0,0)

      if (.not. bas_check_handle (ao_bas_han, pname))
     $     call errquit (pname//"invalid basis handle",0,BASIS_ERR)


C
C     Alloc dipole mats
C     
      if (.not. ga_duplicate (g_ddote, g_dipx, "dipx"))
     $     call errquit (pname//"alloc failed dipx", 0, GA_ERR)
      if (.not. ga_duplicate (g_ddote, g_dipy, "dipy"))
     $     call errquit (pname//"alloc failed dipy", 0, GA_ERR)
      if (.not. ga_duplicate (g_ddote, g_dipz, "dipz"))
     $     call errquit (pname//"alloc failed dipz", 0, GA_ERR)

      if (.not. ga_duplicate (g_ddote, g_mask, "mask"))
     $     call errquit (pname//"alloc failed mask", 0, GA_ERR)


C     (start mask off at all 1's)
      call ga_fill (g_mask, 1d0)


C
C     Mask out small exponents (if desired), this will add 0's to the mask
C
      if (rtdb_get (rtdb, "dft:efield_exp_cutoff",
     $     mt_dbl, 1, exp_cutoff)) then

         if (exp_cutoff .lt. 0d0)
     $        call errquit (pname//"invalid negative exp_cutoff")


         if (.not. bas_numcont (ao_bas_han, ncont))
     $        call errquit (pname//"ncont bas_numcont failed",
     $        0,BASIS_ERR)

         do icont = 1, ncont

            if (.not. bas_continfo(ao_bas_han, icont, l, nprim,
     &           ncontr, isphere))
     &           call errquit(pname//"bas_continfo failed",4, BASIS_ERR)

            if (.not.ma_push_get(mt_dbl,nprim,"exp",lexp,iexp))
     $           call errquit(pname//"make exp failed",0, MA_ERR)

            if (.not. bas_cn2bfr (ao_bas_han, icont, bf_start, bf_end))
     $           call errquit (pname//"bas_cn2bfr failed",0,BASIS_ERR)

            call dfill (nprim, 0d0, dbl_mb(iexp), 1)

            if (.not. bas_get_exponent(ao_bas_han, icont, dbl_mb(iexp)))
     $           call errquit (pname//"bas_get_exponent failed",
     $           0, BASIS_ERR)

c$$$            if (me.eq.0) then
c$$$               write (luout, *) ""            
c$$$               write (luout, *) ""
c$$$               write (luout, *) "icont   :", icont
c$$$               write (luout, *) "l       :", l
c$$$               write (luout, *) "nprim   :", nprim
c$$$               write (luout, *) "ncontr  :", ncontr
c$$$               write (luout, *) "isphere :", isphere
c$$$               write (luout, *) "bf_start:", bf_start
c$$$               write (luout, *) "bf_end  :", bf_end
c$$$               write (luout, *) "coefficients:"
c$$$               do ival = 1, nprim
c$$$                  write (luout, *) dbl_mb (iexp + ival - 1)
c$$$               enddo
c$$$               call util_flush (luout)
c$$$            endif


C     (if this contraction has at least one "small" coefficient we want
C     to cut, zero the parts of the mask that correspond to AOs in this
C     contraction)
            
            cut_this_cont = .false.
            do ival = 1, nprim
               if ( dbl_mb (iexp + ival - 1) .le. exp_cutoff )
     $              cut_this_cont = .true.
            enddo

            if (cut_this_cont) then
               do ibf = bf_start, bf_end
                  call ga_fill_patch (g_mask, 1, nbf_ao, ibf, ibf, 0d0)
                  call ga_fill_patch (g_mask, ibf, ibf, 1, nbf_ao, 0d0)
               enddo

               if (me.eq.0) 
     $           write (luout, "(1x,a,i0,a,es11.4,a,i0,a,i0)")
     $             "Contraction ", icont, " contains exponent below",
     $             exp_cutoff,  "; masked basis function(s): ",
     $             bf_start, " - ", bf_end
            endif

            if (.not. ma_pop_stack (lexp))
     $           call errquit (pname//"pop exp failed", 0, MA_ERR)
         enddo
      endif
      


C
C     Calculate dipole mats and scale by mask
C
      call ga_zero (g_dipx)
      call ga_zero (g_dipy)
      call ga_zero (g_dipz)

      call int_dip_ga (ao_bas_han, ao_bas_han,
     $     g_dipx, g_dipy, g_dipz)

      call ga_elem_multiply (g_dipx, g_mask, g_dipx)
      call ga_elem_multiply (g_dipy, g_mask, g_dipy)
      call ga_elem_multiply (g_dipz, g_mask, g_dipz)

c Daniel (10-10-12): Probably unnecessary...
c      call ga_print (g_dipx)

      
      if (lscale) then
         call ga_norm_infinity (g_dipx, norm)
         scale = 1d0 / norm
         call ga_scale (g_dipx, scale)

         call ga_norm_infinity (g_dipy, norm)
         scale = 1d0 / norm
         call ga_scale (g_dipy, scale)

         call ga_norm_infinity (g_dipz, norm)
         scale = 1d0 / norm
         call ga_scale (g_dipz, scale)
      endif
      
C
C     Compute -D.E
C
      call ga_zero (g_ddote)
      call ga_dadd (efield(1), g_dipx, 1d0, g_ddote, g_ddote)
      call ga_dadd (efield(2), g_dipy, 1d0, g_ddote, g_ddote)  
      call ga_dadd (efield(3), g_dipz, 1d0, g_ddote, g_ddote)

c Tested this with finite differentiation of pol. about the dipole.
c This line causes negative pol. w.r.t. the static response calc.
c      call ga_scale (g_ddote, -1d0)  ! so that we return -D.E. 


      
C
C     Cleanup
C
      if (.not. ga_destroy (g_dipx))
     $     call errquit (pname//"destroy failed dipx", 0, GA_ERR)
      if (.not. ga_destroy (g_dipy))
     $     call errquit (pname//"destroy failed dipy", 0, GA_ERR)
      if (.not. ga_destroy (g_dipz))
     $     call errquit (pname//"destroy failed dipz", 0, GA_ERR)
      if (.not. ga_destroy (g_mask))
     $     call errquit (pname//"destroy failed mask", 0, GA_ERR)

      end
