c
c     print vem output
c
      subroutine cosmo_vem_output(estote, espol, gstote, gspolneq, 
     $              gspoldyn)
c
      implicit none
c
#include "stdio.fh"
#include "cosmo.fh"
#include "util_params.fh"
c
      double precision estote, espol
      double precision gstote, gspolneq, gspoldyn
c
      write (luout,930)
      write (luout,931)
      write (luout,932) estote
      write (luout,933) espol,espol*cau2ev
      write (luout,934) gstote
      write (luout,935) gspolneq,gspolneq*cau2ev
      write (luout,936) gspoldyn,gspoldyn*cau2ev
c     write (luout,937) gspol,gspol*cau2ev
      write (luout,938) (estote-gstote),(estote-gstote)*cau2ev
c
 930  format(/,
     $'                          COSMO-VEM solvation results',/,
     $'                          ---------------------------',/
     $' Reference for the VEM model:',/,
     $' Marenich, A. V.; Cramer, C. J.; Truhlar, D. G.;',
     $' Guido, C. A.; Mennucci, B.;',/,' Scalmani, G.; Frisch, M. J.',
     $' Chem. Sci. 2011, 2, 2143',/)               
 931  format(1x,
     $'emission spectrum data: ES = initial state, GS = final state')
 932  format(1x,
     $'(1)   VEM ES equilibrium total free energy = ',
     $ f20.10)
 933  format(1x,
     $'(2)            ES polarization free energy = ',
     $ f20.10,' (',f10.4,' eV)')
 934  format(1x,
     $'(3)    GS nonequilibrium total free energy = ',
     $ f20.10)
 935  format(1x,
     $'(4)            GS polarization free energy = ',
     $ f20.10,' (',f10.4,' eV)')
 936   format(1x,
     $'(5)     fast polarization component of (4) = ',
     $ f20.10,' (',f10.4,' eV)')
 937  format(1x,
     $'(5a)        1/2 Vgs * ( Qes_in + Qgs_dyn ) = ',
     $ f20.10,' (',f10.4,' eV)')
 938  format(1x,
     $'(6) VEM vertical emission energy (1) - (3) = ',
     $ f20.10,' (',f10.4,' eV)')

      end 
