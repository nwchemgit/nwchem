      subroutine  vib_thermo(rtdb,nat,nat3,
     &    eigvals,masses,coords)
      implicit none
*
* $Id: vib_thermo.F,v 1.1 1999-07-02 23:32:01 d3e129 Exp $
*
* routine to compute the thermodynamic properties under the 
* rigid rotor harmonic approximation:
*
* See:
* 1. "Statistical Mechanics" by Donald A. McQuarrie, Chapters 5,6,8,
*     Publisher Harper and Row, (c) 1976
* 2. 
*     
*:passed::
      integer rtdb ! [input] rtdb handle
      integer nat  ! [input] number of atoms
      integer nat3 ! [input] 3*nat
      integer eigvals(nat3) ! [input] eigenvalues (in cm**-1)
      integer masses(nat)   ! [input] mass of atoms (in amu)
      integer coords(nat)   ! [input] atomic coordinates (in a.u.)
      return
      end
