      subroutine vib_intense(eigvec,eigval,nat,ddpolx,ddpolq,
     &    intensity,first_pass)
      implicit none
#include "stdio.fh"
#include "util.fh"
      integer nat
      double precision eigvec(3,nat,3,nat)  ! (cart,q)
      double precision eigval(3,nat)        ! (q)
      double precision ddpolx(3,3,nat)      ! (m,cart)
      double precision ddpolq(3,3,nat)      ! (m,q)
      double precision intensity(3,nat,3)   ! (q,3)  3->1 = au, 
* . . . . . . . . . . . . . . . . . . . . . !           2 = (debye/angs)**2,
* . . . . . . . . . . . . . . . . . . . . . !           3 = arbitrary units 
* . . . . . . . . . . . . . . . . . . . . . !               normalized units
      logical first_pass
c
      integer iq, iqatom, iqxyz, iqm, atom, xyz
      double precision ddq, dint, ave_intensity, factor
      logical debug
c
      debug = .false.
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
*(cart,q)
     &              ddpolx(iqm,xyz,atom)*eigvec(xyz,atom,iqxyz,iqatom)
              enddo
            enddo
            ddpolq(iqm,iqxyz,iqatom) = ddq
          enddo
        enddo
      enddo
c
*... form Intensity vector
c
      call dfill ((3*nat*3),0.0d00,intensity,1)
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
            write(luout,10000)
     &          iq,
     &          eigval(iqxyz,iqatom),
     &          intensity(iqxyz,iqatom,1),
     &          intensity(iqxyz,iqatom,2)
          enddo
        enddo
      endif
c
*... normalize the average intensity to max = 10.0
c  10.0 = value*factor
c
      factor = 10.d00/ave_intensity
      call dcopy((3*nat),intensity(1,1,1),1,intensity(1,1,3),1)
      call dscal((3*nat),factor,intensity(1,1,3),1)
c
      if (debug) then
        iq = 0
        do iqatom = 1,nat
          do iqxyz = 1,3
            iq = iq + 1
            write(luout,10001)
     &          iq,
     &          eigval(iqxyz,iqatom),
     &          intensity(iqxyz,iqatom,3)
          enddo
        enddo
      endif
c
*... print information out.
c
      if (util_print('IR intensities',print_low)) then
        write(luout,*)' '
        write(luout,*)' '
        write(luout,10002)
        if (first_pass) then
          write(luout,10003)
        else
          write(luout,10004)
        endif
        iq = 0
        do iqatom = 1,nat
          do iqxyz = 1,3
            iq = iq + 1
            write(luout,10005)
     &          iq,eigval(iqxyz,iqatom),
     &          (intensity(iqxyz,iqatom,iqm),iqm=1,3)
          enddo
        enddo
        write(luout,10002)
        write(luout,*)' '
        write(luout,*)' '
        call util_flush(luout)
      endif
c
10000 format(1x,'Intensity for eigenvalue ',i3,' [',f10.3,
     &    ' cm**(-1)] is ',
     &    f10.3,' [atomic units] ',
     &    f10.3,'[(debey/angs)**2]' )
10001 format(1x,'Normalized Intensity for eigenvalue ',i3,' [',f10.3,
     &    ' cm**(-1)] is ',f10.3,' arbitrary units ')
10002 format(1x,65('-'))
10003 format(
     &    1x,'Normal',1x,'Eigenvalue',1x,'||',
     &    8x,'   Infra Red Intensities',/,
     &    2x,'Mode',
     &    3x,'[cm**-1]',2x,'||',
     &    1x,'[atomic units]',
     &    1x,'[(debye/angs)**2]',
     &    1x,'[arbitrary]',/,
     &    1x,'------',
     &    1x,'----------',1x,'||',
     &    1x,'--------------',
     &    1x,'-----------------',
     &    1x,'-----------')
10004 format(
     &    1x,'Normal',1x,'Eigenvalue',1x,'||',
     &    8x,'   Projected Intensities',/,
     &    2x,'Mode',
     &    3x,'[cm**-1]',2x,'||',
     &    1x,'[atomic units]',
     &    1x,'[(debye/angs)**2]',
     &    1x,'[arbitrary]',/,
     &    1x,'------',
     &    1x,'----------',1x,'||',
     &    1x,'--------------',
     &    1x,'-----------------',
     &    1x,'-----------')
10005 format(2x,i3,3x,f10.3,1x,'||',
     &    1x,f11.6,6x,f10.3,4x,f10.3)
c
      end
