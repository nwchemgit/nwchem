      subroutine bse_oscstr(pars,omega,totpoles)
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "bse.fh"
      type(bse_params_t) :: pars

      character(*), parameter :: pname = 'bse_oscstr: '
      character(len=7) :: key

      integer totpoles
      double precision omega(*)
      double precision oscstr,tmp,moms(3)

      integer llpole,ulpole,iroot,ixyz,ipole

      if (pars%singlet) then
        key = 'singlet'
      else
        key = 'triplet'
      endif

      ! Get the lower and upper limits
      llpole = -1
      if (pars%dowindow) then
        iroot = 0
        do ipole=1,totpoles
          if (omega(ipole).lt.pars%elower) cycle
          if (llpole.eq.-1) llpole = ipole
          iroot = iroot + 1
          if (omega(ipole).gt.pars%eupper) then
            ulpole = ipole - 1
            exit
          endif
          if (iroot.eq.pars%nroots) then
            ulpole = ipole
            exit
          endif
        enddo
      else
        llpole = 1
        ulpole = pars%nroots
      endif
          

      ! Compute the oscillator strengths around an energy window
      iroot = 0
      do ipole=llpole,ulpole
        iroot = iroot + 1
        oscstr = 0d0
        moms(:) = 0d0
        do ixyz=1,3
          moms(ixyz) = moms(ixyz) + nga_ddot_patch(pars%g_apb,'n',
     &          (/1,ipole/),(/pars%npoles(1),ipole/),
     &          pars%g_dipa(ixyz),'n',1,pars%npoles(1))
          oscstr = oscstr + moms(ixyz)**2
        enddo
        if (pars%ipol.gt.1) then
          do ixyz=1,3
            moms(ixyz) = moms(ixyz) + 
     &        nga_ddot_patch(pars%g_apb,'n',(/pars%npoles(1)+1,ipole/),
     &            (/totpoles,ipole/),
     &            pars%g_dipb(ixyz),'n',1,pars%npoles(2))
            oscstr = oscstr + moms(ixyz)**2
          enddo
        endif
        oscstr = 2d0/3d0*omega(ipole)*oscstr

        if (pars%me.eq.0) then
          write(luout,9110) iroot,key,'a1  ',omega(ipole),
     &                      omega(ipole)*ha2ev
          write(luout,9200) moms(1),moms(2),moms(3),oscstr
        endif
      enddo

 9110 format(2x,
     1  '---------------------------------------------------------------
     3-------------',
     2  /,2x,'Root',i4,1x,a7,1x,a4,f22.9,' a.u.',f22.4,' eV ',/,2x,
     1  '---------------------------------------------------------------
     3-------------')
 9200 format(5x,'Transition Moments    X',f9.5,'   Y',f9.5,'   Z',f9.5
     5    ,/,5x,'Dipole Oscillator Strength',18x,f14.10)

      end
