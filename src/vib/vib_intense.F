* $Id$
      subroutine vib_intense(rtdb,eigvec,eigval,nat,ddpolx,ddpolq,
     &    intensity,first_pass)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "util_params.fh"
      integer nat
      integer rtdb
      double precision eigvec(3,nat,3,nat)  ! (cart,q)
      double precision eigval(3,nat)        ! (q)
      double precision ddpolx(3,3,nat)      ! (m,cart)
      double precision ddpolq(3,3,nat)      ! (m,q)
      double precision intensity(3,nat,4)   ! (q,3)  3->1 = au, 
* . . . . . . . . . . . . . . . . . . . . . !           2 = (debye/angs)**2,
* . . . . . . . . . . . . . . . . . . . . . !           3 = (KM/mole)
* . . . . . . . . . . . . . . . . . . . . . !           4 = arbitrary units 
* . . . . . . . . . . . . . . . . . . . . . !               normalized units
      logical first_pass
c
      integer iq, iqatom, iqxyz, iqm, atom, xyz
      double precision ddq, dint, ave_intensity, factor
      double precision dpdq1, dpdq2, dpdq3
      integer ivalx, ivala
      double precision val
      logical printit
      logical debug
c
      debug = util_print('everything',print_debug)
c
c fix phase of eigvec
c
      do iqatom = 1,nat
        do iqxyz = 1,3
          ivalx = -1
          ivala = -1
          val = -1.0d00
          do atom = 1,nat
            do xyz = 1,3
              if (abs(eigvec(xyz,atom,iqxyz,iqatom)).gt.val) then
                val = eigvec(xyz,atom,iqxyz,iqatom)
                ivala = atom
                ivalx = xyz
              endif
            enddo
          enddo
          if (eigvec(ivalx,ivala,iqxyz,iqatom).lt.0.0d00) then
            do atom = 1,nat
              do xyz = 1,3
                eigvec(xyz,atom,iqxyz,iqatom)=
     &              -eigvec(xyz,atom,iqxyz,iqatom)
              enddo
            enddo
          endif
        enddo
      enddo
c
*...  first form ddpolq from ddpolx
c
      do iqatom = 1,nat
        do iqxyz = 1,3
          do iqm = 1,3
            ddq = 0.0d00
            do atom = 1,nat
              do xyz = 1,3
                ddq = ddq +
*ddpolx(q,cart)eigvec(cart,q)
     &              ddpolx(iqm,xyz,atom)*eigvec(xyz,atom,iqxyz,iqatom)
              enddo
            enddo
            ddpolq(iqm,iqxyz,iqatom) = ddq
          enddo
        enddo
      enddo
c
c write dipole derivatives in normal modes 
c
      printit = util_print('dipole derivatives',print_default)
      printit = printit.or.debug
      if (printit) then
        write(luout,'(//)')
        write(luout,10000)
        if (first_pass) then
          write(luout,10001)
        else
          write(luout,10002)
        endif
        iq = 0
        do iqatom = 1,nat
          do iqxyz = 1,3
            iq = iq + 1
            dpdq1 = ddpolq(1,iqxyz,iqatom)*2.5417478d00/cau2ang
            dpdq2 = ddpolq(2,iqxyz,iqatom)*2.5417478d00/cau2ang
            dpdq3 = ddpolq(3,iqxyz,iqatom)*2.5417478d00/cau2ang
            write(luout,10003)
     &          iq,eigval(iqxyz,iqatom),dpdq1,dpdq2,dpdq3
          enddo
        enddo
        write(luout,10000)
        write(luout,'(//)')
      endif
c
*... form Intensity vector
c
      call dfill ((3*nat*4),0.0d00,intensity,1)
c
      ave_intensity = 0.0d00
      do iqatom = 1,nat
        do iqxyz = 1,3
          dint = 0.0d00
          do iqm = 1,3
            dint = dint +
     &          ddpolq(iqm,iqxyz,iqatom)*ddpolq(iqm,iqxyz,iqatom)
          enddo
          intensity(iqxyz,iqatom,1) = dint
          ave_intensity = ave_intensity + dint
*
* from A Physicists Desk Reference, The Second Edition of Physics Vade Mecum
*      Herbert L. Anderson, Editor in Chief
*      Copyright (C) 1989 American Institute of Physics
*      335 East 45th Street, New York, NY 10017
*
*1 debye = 10**(-18) esu cm * [1 e/4.8032068 x 10**(-10) esu]*[1 m /100cm]*[a0/5.29177249 x 10**(-11) m]
*1 debye = (1.0/4.8032068/5.29177249) * 10**(-18 + 10 - 2 + 11) e a0
*1 debye = (1.0/4.8032068/5.29177249) * 10**(1) e a0
*1 e a0  = (4.8032068*5.29177249) * 10**(-1) debye
*1 e a0  = 25.417477608 * 10**(-1) debye
*1 e a0  = 2.5417477608 debye
*
*use 1 e a0 = 2.541 7478 debye
*
          intensity(iqxyz,iqatom,2) =
     &        dint*2.5417478d00/0.529177249d00 ! bohr->angstrom matches current geom data
     &        *2.5417478d00/0.529177249d00
          intensity(iqxyz,iqatom,3) =
     &        intensity(iqxyz,iqatom,2)*42.2547  ! convert to KM/mol from debye**2/angs**2
        enddo
      enddo
c
      ave_intensity = ave_intensity/3.0d00/nat
c
      if (debug) then
        iq = 0
        do iqatom = 1,nat
          do iqxyz = 1,3
            iq = iq + 1
            write(luout,90000)
     &          iq,
     &          eigval(iqxyz,iqatom),
     &          intensity(iqxyz,iqatom,1),
     &          intensity(iqxyz,iqatom,2),
     &          intensity(iqxyz,iqatom,3)
          enddo
        enddo
      endif
c
*... normalize the average intensity to max = 10.0
c  10.0 = value*factor
c
      factor = 10.d00/ave_intensity
      call dcopy((3*nat),intensity(1,1,1),1,intensity(1,1,4),1)
      call dscal((3*nat),factor,intensity(1,1,4),1)
c
      if (debug) then
        iq = 0
        do iqatom = 1,nat
          do iqxyz = 1,3
            iq = iq + 1
            write(luout,90001)
     &          iq,
     &          eigval(iqxyz,iqatom),
     &          intensity(iqxyz,iqatom,4)
          enddo
        enddo
      endif
c
*... print information out.
c
      printit = util_print('IR intensities',print_low)
      printit = printit.or.debug
      if (printit) then
        write(luout,*)' '
        write(luout,*)' '
        write(luout,10000)
        if (first_pass) then
          write(luout,10004)
        else
          write(luout,10005)
        endif
        iq = 0
        do iqatom = 1,nat
          do iqxyz = 1,3
            iq = iq + 1
            write(luout,10006)
     &          iq,eigval(iqxyz,iqatom),
     &          (intensity(iqxyz,iqatom,iqm),iqm=1,4)
          enddo
        enddo
        write(luout,10000)
        write(luout,'(//)')
        call util_flush(luout)
      endif
c
*... print information to ecce file
c      
      if (first_pass) then
        call ecce_print1('intensities',
     &      mt_dbl,intensity(1,1,1),(3*nat))
        if(.not. rtdb_put(rtdb,'vib:intensities',MT_DBL,3*nat,
     &      intensity(1,1,1)))
     &      call errquit('vib_intense:rtdb_put of intensities failed',
     &      555, RTDB_ERR)
        call ecce_print1('intensities (debye/ang)2',
     &      mt_dbl,intensity(1,1,2),(3*nat))
        call ecce_print1('intensities (KM/mol)',
     &      mt_dbl,intensity(1,1,3),(3*nat))
      else
        call ecce_print1('projected intensities',
     &      mt_dbl,intensity(1,1,1),(3*nat))
        if(.not. rtdb_put(rtdb,'vib:projected intensities',
     &      MT_DBL,3*nat,intensity(1,1,1)))
     &      call errquit(
     &      'vib_intense:rtdb_put of projected intensities failed',
     &      555, RTDB_ERR)
        call ecce_print1('projected intensities (debye/ang)2',
     &      mt_dbl,intensity(1,1,2),(3*nat))
        call ecce_print1('projected intensities (KM/mol)',
     &      mt_dbl,intensity(1,1,3),(3*nat))
      endif
c
10000 format(1x,76('-'))
10001 format(
     &    1x,'Normal',1x,'Eigenvalue',1x,'||',
     &    9x,'Derivative Dipole Moments (debye/angs)',/,
     &    2x,'Mode',
     &    3x,'[cm**-1]',2x,'||',
     &    6x,'[d/dqX]',
     &   13x,'[d/dqY]',
     &   11x,'[d/dqZ]',/,
     &    1x,'------',
     &    1x,'----------',1x,'||',
     &    1x,'------------------',
     &    1x,'------------------',
     &    1x,'-----------------')
10002 format(
     &    1x,'Normal',1x,'Eigenvalue',1x,'||',
     &    4x,'Projected Derivative Dipole Moments (debye/angs)',/,
     &    2x,'Mode',
     &    3x,'[cm**-1]',2x,'||',
     &    6x,'[d/dqX]',
     &   13x,'[d/dqY]',
     &   11x,'[d/dqZ]',/,
     &    1x,'------',
     &    1x,'----------',1x,'||',
     &    1x,'------------------',
     &    1x,'------------------',
     &    1x,'-----------------')
10003 format(2x,i3,3x,f10.3,1x,'||',
     &    2x,f10.3,10x,f10.3,8x,f10.3)
10004 format(
     &    1x,'Normal',1x,'Eigenvalue',1x,'||',
     &   14x,'   Infra Red Intensities',/,
     &    2x,'Mode',
     &    3x,'[cm**-1]',2x,'||',
     &    1x,'[atomic units]',
     &    1x,'[(debye/angs)**2]',
     &    1x,'[(KM/mol)]',
     &    1x,'[arbitrary]',/,
     &    1x,'------',
     &    1x,'----------',1x,'||',
     &    1x,'--------------',
     &    1x,'-----------------',
     &    1x,'----------',
     &    1x,'-----------')
10005 format(
     &    1x,'Normal',1x,'Eigenvalue',1x,'||',
     &    8x,'   Projected Infra Red Intensities',/,
     &    2x,'Mode',
     &    3x,'[cm**-1]',2x,'||',
     &    1x,'[atomic units]',
     &    1x,'[(debye/angs)**2]',
     &    1x,'[(KM/mol)]',
     &    1x,'[arbitrary]',/,
     &    1x,'------',
     &    1x,'----------',1x,'||',
     &    1x,'--------------',
     &    1x,'-----------------',
     &    1x,'----------',
     &    1x,'-----------')
10006 format(2x,i3,3x,f10.3,1x,'||',
     &    1x,f11.6,6x,f10.3,4x,f10.3,2x,f10.3)
90000 format(1x,'Intensity for eigenvalue ',i3,' [',f10.3,
     &    ' cm**(-1)] is ',
     &    f10.3,' [atomic units] ',
     &    f10.3,'[(debey/angs)**2]',
     &    f10.3,'[(KM/mol)]')
90001 format(1x,'Normalized Intensity for eigenvalue ',i3,' [',f10.3,
     &    ' cm**(-1)] is ',f10.3,' arbitrary units ')
c
      end
