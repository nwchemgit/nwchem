c
c     $Id$                       
c

*  ************************************************************
*  *                                                          *
*  *             Band by Band Kohn-Sham Minimizer             *
*  *                    (fixed densities)                    *
*  *                                                          *
*  *                                                          *
*  ************************************************************
      subroutine bybminimize0()
      implicit none

#include "stdio.fh"
#include "bafdecls.fh"
#include "util.fh"

*     **** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)

      logical oprint
      integer iterations,it_out,it_ortho
      integer ito
      
*     **** external functions ****
      logical  control_precondition,control_print
      external control_precondition,control_print
      integer  control_ks_algorithm
      external control_ks_algorithm
      real*8   control_tole
      external control_tole
      integer  control_H1_it_in,control_H1_it_out, control_H1_it_ortho
      external control_H1_it_in,control_H1_it_out, control_H1_it_ortho

      iterations = control_H1_it_in()
      it_out     = control_H1_it_out()
      it_ortho   = control_H1_it_ortho()

            
      if (it_ortho.ge.0) then
         call Parallel_taskid(taskid)
         oprint = (taskid.eq.MASTER).and.control_print(print_medium)
         if (oprint) then
           write(luout,'(A,A,I0,A,I0,A,I0)') 
     >                 "     -  one_electron_guess,",
     >                 " it_in=",iterations,
     >                 " it_out=",it_out,
     >                 " it_ortho=",it_ortho
         end if

         call electron_gen_vall0()

         do ito=0,it_ortho

            if (ito.gt.0) then
*              **** diaganolize KS matrix - note that Hpsik set to last application of H(rho)|psi> ****
               call psi_1gen_hml()
               call psi_diagonalize_hml()
               call psi_1rotate2()
               call psi_2to1()
            end if

            call psi_KS_update00(1,it_out,iterations,
     >                        control_ks_algorithm(),
     >                        control_precondition(),
     >                        control_tole())
         end do
      end if
      return
      end
 

