      subroutine  vib_thermo(rtdb,nat,nat3,
     &    eigvals,masses,coords,scr)
      implicit none
#include "errquit.fh"
*
* $Id$
*
* routine to compute the thermodynamic properties under the 
* rigid rotor harmonic approximation:
*
* See:
* 1. "Statistical Mechanics" by Donald A. McQuarrie, Chapters 5,6,8,
*     Publisher Harper and Row, (c) 1976
* 2. 
*     
*::functions::
#include "rtdb.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "stdio.fh"
#include "inp.fh"
*::passed::
      integer rtdb ! [input] rtdb handle
      integer nat  ! [input] number of atoms
      integer nat3 ! [input] 3*nat
      double precision eigvals(nat3) ! [input] eigenvalues (in cm**-1)
      double precision masses(nat)   ! [input] mass of atoms (in amu)
      double precision coords(3,nat) ! [input] atomic coordinates (in a.u.)
      double precision scr(*)
*::local::
      integer i, j, count, error
      integer iat
      integer len_theory
      integer nopen, multiplicity
      logical linear
      double precision com(3) ! center of mass coordinates
      double precision Inertia(3,3), Inertia2(3,3), rx, ry, rz
      double precision VdI(3,3)
      double precision Ivals(3), Ivecs(3,3), Itri(3*4/2)
      double precision total_mass, xmass
      double precision H_trans , G_trans , S_trans , Cp_trans , Cv_trans
      double precision H_rot   , G_rot   , S_rot   , Cp_rot   , Cv_rot  
      double precision H_vib   , G_vib   , S_vib   , Cp_vib   , Cv_vib  
      double precision H_elec  , G_elec  , S_elec  , Cp_elec  , Cv_elec 
      double precision H_total , G_total , S_total , Cp_total , Cv_total
      double precision Temp
      double precision R_gas
      double precision Pi
      double precision Sigma
      character*32 theory
      character*50 my_string
*
      double precision thresh
      parameter (thresh=1.0d-3)
c::-statement function
      logical is_it_close_to  
      double precision value,test
      intrinsic abs
*---          is value close to test?
      is_it_close_to(value,test) = (abs(value-test).lt.thresh)
*
* -- until it is working
      if (nat.gt.0) return
*     
      Pi = 2.0d00*acos(0.0d00)
      H_trans  = 0.0d00
      H_rot    = 0.0d00
      H_vib    = 0.0d00
      H_elec   = 0.0d00
      H_total  = 0.0d00
      G_trans  = 0.0d00
      G_rot    = 0.0d00
      G_vib    = 0.0d00
      G_elec   = 0.0d00
      G_total  = 0.0d00
      S_trans  = 0.0d00
      S_rot    = 0.0d00
      S_vib    = 0.0d00
      S_elec   = 0.0d00
      S_total  = 0.0d00
      Cp_trans = 0.0d00
      Cp_rot   = 0.0d00
      Cp_vib   = 0.0d00
      Cp_elec  = 0.0d00
      Cp_total = 0.0d00
      Cv_trans = 0.0d00
      Cv_rot   = 0.0d00
      Cv_vib   = 0.0d00
      Cv_elec  = 0.0d00
      Cv_total = 0.0d00
*-- electronic component 
      my_string = ' '
      theory = ' '
      if (.not.rtdb_cget(rtdb, 'task:theory', 1, theory)) call errquit
     &      ('vib_thermo: could not read task:theory',911, RTDB_ERR)
      len_theory = inp_strlen(theory)
      my_string = theory(1:len_theory) // ':nopen'
      nopen = 0
      if (.not.rtdb_get(rtdb,my_string,mt_int,1,nopen)) then
        write(luout,*)' theory is :',theory(1:len_theory)
        call errquit('vib_thermo: could not read --theory--:nopen',911,
     &       RTDB_ERR)
      endif
      multiplicity = nopen + 1
      S_elec  = R_gas*log(dble(multiplicity))
      H_elec  = 0.0d00
      Cp_elec = 0.0d00
      Cv_elec = 0.0d00
      G_elec  = H_elec - Temp * S_elec
*-- rotational component
*- compute total mass and center of mass
      call dfill(3,0.0d00,com,1)
      total_mass = 0.0d00
      do iat = 1,nat
        xmass = masses(iat)
        total_mass = total_mass + xmass
        com(1) = com(1) + xmass*coords(1,iat)
        com(2) = com(2) + xmass*coords(2,iat)
        com(3) = com(3) + xmass*coords(3,iat)
      enddo
      call dscal(3,(1.0d00/total_mass),com,1)
*
      write(luout,*)' Coordinates '
      call output(coords,1,3,1,3,3,3,1)
      write(luout,*)' Center of Mass :',com
*
      call dfill((3*3),0.0d00,Inertia,1)
      do iat = 1,nat
        xmass = masses(iat)
        rx = coords(1,iat) - com(1)
        ry = coords(2,iat) - com(2)
        rz = coords(3,iat) - com(3)
        Inertia(1,1) = Inertia(1,1) +
     &        xmass*(ry*ry + rz*rz)
        Inertia(2,1) = Inertia(2,1) -
     &        xmass*ry*rx
        Inertia(3,1) = Inertia(3,1) -
     &        xmass*rz*rx
        Inertia(2,2) = Inertia(2,2) +
     &        xmass*(rx*rx + rz*rz)
        Inertia(3,2) = Inertia(3,2) -
     &        xmass*rz*ry
        Inertia(3,3) = Inertia(3,3) +
     &        xmass*(rx*rx + ry*ry)
      enddo
      Inertia(1,2) = Inertia(2,1)
      Inertia(1,3) = Inertia(3,1)
      Inertia(2,3) = Inertia(3,2)
      write(luout,*)' Raw Inertial Matrix'
      call output(Inertia,1,3,1,3,3,3,1)
* set up triangular matrix
      count = 0
      do i = 1,3
        do j = 1,i
          count = count + 1
          Itri(count) = Inertia(i,j)
        enddo
      enddo
* diagonilze Inertia Matrix (parinoia)
      call dfill(3,0.0d00,Ivals,1)
      call dfill((3*3),0.0d00,Ivecs,1)
      call vib_sjacobi(3,3,Itri,Ivals,Ivecs,error)
      if (error.ne.0) call errquit('vib_thermo: vib_sjacobi failed',911,
     &       UNKNOWN_ERR)
*
      write(luout,*) 'eigenvalues of Inertia Matrix'
      call output(Ivals,1,3,1,1,3,1,1)
      write(luout,*) 'eigenvectors of Inertia Matrix'
      call output(Ivecs,1,3,1,3,3,3,1)
* form "diagonal" matrix      
      call dfill((3*3),0.0d00,Inertia2,1)
      call dgemm('t','n',3,3,3,1.0d00,Ivecs,3,Inertia,3,0.0d00,VdI,3)
      call dgemm('n','n',3,3,3,1.0d00,VdI,3,Ivecs,3,0.0d00,Inertia2,3)
      write(luout,*)' Diagonalized Inertial Matrix'
      call output(Inertia2,1,3,1,3,3,3,1)
* 
      linear = is_it_close_to(Ivals(1),0.0d00) .and.
     &         is_it_close_to((abs(Ivals(2)-Ivals(3))),0.0d00)
*
      sigma = 1
      if (linear) then
        S_rot = R_gas*(1 + log(Ivals(2)/Sigma*1.0d00))
        if (rtdb_put(rtdb,'vib:linear',mt_log,1,.true.))
     &    call errquit('vib_thermo: failed to set linear',555, RTDB_ERR)
      else
      endif
*
      end
