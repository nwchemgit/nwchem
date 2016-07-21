c
c     print cosmo, smd output      
c
      subroutine cosmo_smd_output(egas, esol, do_cosmo_smd, gcds, 
     $  gstote, gspol) 
c
      implicit none
c
#include "stdio.fh"
#include "cosmo.fh"
#include "util_params.fh"
c
      double precision egas, esol, gcds
      double precision gstote, gspol
      logical do_cosmo_smd
c
c     Standard cosmo output
c
        if (.not.do_cosmo_smd) then
           write(luout,911)
           if (dabs(egas).gt.0.d0) then
                write(luout,912) egas
                write(luout,913) esol
                write(luout,914) (egas-esol),(egas-esol)*cau2kcal
           else
                write(luout,915)
                write(luout,913) esol
           endif
c         
        else  ! cosmo-smd output
c
           write(luout,819)
           if (dabs(egas).gt.0.d0) then
             write(luout,820) egas
             write(luout,821) (gstote-gspol)
             write(luout,822) (gstote-gspol-egas),
     $                   (gstote-gspol-egas)*cau2kcal
             write(luout,823) gstote
             write(luout,824) gspol,gspol*cau2kcal
             write(luout,825) esol
             write(luout,826) gcds,gcds*cau2kcal
             write(luout,827) (esol-egas),(esol-egas)*cau2kcal
           else
             write(luout,821) (gstote-gspol)
             write(luout,823) gstote
             write(luout,824) gspol,gspol*cau2kcal
             write(luout,825) esol
             write(luout,826) gcds,gcds*cau2kcal
           endif
        endif
c
c    Standard cosmo output
c
 911  format(
     $32x,'   COSMO solvation results',/,
     $32x,'   -----------------------',/)
 912  format('                 gas phase energy = ',f22.12)
 913  format('                 sol phase energy = ',f22.12)
 914  format(' (electrostatic) solvation energy = ',f22.12,
     $     ' (',f8.2,' kcal/mol)'                       )
 915  format('     skipped: no gas phase energy')
c
c     Custom cosmo-smd output
c
 819  format(
     $32x,'COSMO-SMD solvation results',/,
     $32x,'---------------------------',/
     $' Reference for the SMD model:',/,
     $' Marenich, A. V.; Cramer, C. J.; Truhlar, D. G.',
     $' J. Phys. Chem. B 2009, 113, 6378',/)
 820  format(1x,
     $'            internal energy in gas <Psi(g)|H|Psi(g)> = ',
     $ f22.12)
 821  format(1x,
     $'        internal energy in solvent <Psi(s)|H|Psi(s)> = ',
     $ f22.12)
 822  format(1x,
     $'                              delta internal energy  = ',
     $ f22.12,' (',f8.2,' kcal/mol)')
 823  format(1x,
     $'  total free energy in solvent <Psi(s)|H+V/2|Psi(s)> = ',
     $ f22.12)
 824  format(1x,
     $'polarization energy contribution <Psi(s)|V/2|Psi(s)> = ',
     $ f22.12,' (',f8.2,' kcal/mol)')
 825  format(1x,
     $'   total free energy in solvent including G(SMD-CDS) = ',
     $ f22.12)
 826  format(1x,
     $'                      G(SMD-CDS) energy contribution = ',
     $ f22.12,
     $     ' (',f8.2,' kcal/mol)')
 827  format(1x,
     $'    1 M fixed-concentration free energy of solvation = ',
     $ f22.12,' (',f8.2,' kcal/mol)'                       )
c
      end
