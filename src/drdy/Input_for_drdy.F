C $Id: Input_for_drdy.F,v 1.2 2001-01-08 15:44:17 windus Exp $
      Subroutine Input_for_drdy(rtdb)
*     Implicit Double Precision (A-H,O-Z)
      Implicit none
#include "stdio.fh"
#include "drdyP.fh"
*::passed:
      Integer rtdb ! [input] RTDB handle
*::locals:
      logical nwc_start
      logical nwc_continue
      logical nwc_restart
      Integer fu5,fu6
      Logical lexit
      Integer i, j, l, iop
C
C Read in data from unit fu5
      fu5 = luin
      fu6 = luout
C
C     Initialize values
C
      do i = 1,39
       lgs(i) = 0
      enddo
      do i = 1,5
       irepr(i) = 0  
       icode(i) = 4
      enddo
      lgs(1) = 1
      dels = 0.01d0
      delsv = 0.1d0
      slp = 1.0d0
      slm = -1.0d0
C     
****************************************************     call read5
      call read_drdy_rtdb(rtdb)
**************************************************** check for global NWChem restart
      nwc_start    = .false.
      nwc_continue = .false.
      nwc_restart  = .false.
      call util_get_rtdb_state(rtdb,nwc_start,nwc_continue,nwc_restart)
      if (nwc_continue.or.nwc_restart) then
        lgs(8) = 1
      endif
*
      if (irepr(1).eq.1.and.irepr(2).eq.1) then
         if (irepr(3).eq.1.and.irepr(4).eq.1) lgs(6) = 1
         if (irepr(3).eq.1.and.irepr(4).eq.0) lgs(6) = 2
      endif
      if (irepr(1).eq.1.and.irepr(2).eq.0) then
         if (irepr(3).eq.1.and.irepr(4).eq.1) lgs(6) = 3
         if (irepr(3).eq.1.and.irepr(4).eq.0) lgs(6) = 4
      endif
C    
      write (fu6,600) (lgs(i),i=1,39)
      if (lgs(1).eq.0) then
         write (fu6,601) lgs(1)
      else
         lgs(1) = 1
         write (fu6,602) lgs(1)
      endif
      if (lgs(6).le.1) then
         lgs(6) = 1
         write (fu6,603) lgs(6),2,2
      else if (lgs(6).eq.2) then
         write (fu6,603) lgs(6),2,1
      else if (lgs(6).eq.3) then
         write (fu6,603) lgs(6),1,2
      else
         lgs(6) = 4
         write (fu6,603) lgs(6),1,1
      endif
      if (lgs(8).le.0) then
         lgs(8) = 0
         write (fu6,604) lgs(8)
      else
         lgs(8) = 1
         write (fu6,605) lgs(8)
      endif
      if (lgs(27).eq.0) then
         write (fu6,606) lgs(27)
         dir = 1
      else
         lgs(27) = -1
         dir = -1
         write (fu6,607)
      endif
      if(lgs(31).le.1) then
         lgs(31) = 1
         write (fu6,608) lgs(31)
      else if (lgs(31).eq.2) then
         write (fu6,609) lgs(31)
      else if (lgs(31).eq.3) then
         write (fu6,610) lgs(31)
      else 
         write (fu6,611) lgs(31)
      endif
      write (fu6,612) ctitle
      n3 = 3*natom 
      write (fu6,613) natom,n3
      if (natom.gt.natoms) then
         lexit = .true.
         write (fu6,*) ' number of atoms is greater than max=',natoms
      endif
      write (fu6,614)
      write (fu6,615) (i,xmass(i),xmass(i)/cau,xmass(i)*ckg,i=1,natom)
      write (fu6,616) redm,redm/cau,redm*ckg
C
c Scale atomic masses to reduced mass of reactants
c  redefine xmass array to hold coordinate mass factors repeated for each coordinate
      do i = 1, natom
         xmass(i) = sqrt(xmass(i)/redm)
      enddo
      l = 0
      do i = 1, natom
         do j = 1, 3
            l = l+1
            amass(l) = xmass(i)
         enddo
      enddo
c input for reactant and product properties
C   icode = 1 for an atomic, 
C         = 2 for a diatomic treated as a Morse oscillator
C         = 3 for a diatomic or linear polyatomic,
C         = 4 for a nonlinear polyatomic GTS
      iop = 1
      write(fu6,620)
      write (fu6,623) 
      write (fu6,624) (i,(xr(3*i-3+j,iop),j=1,3),i=1,natom)
      write (fu6,625) iop,(iatsv(i,iop),i=1,nratom(iop))
      write (fu6,626) icode(iop)
C transform to mass-scaled coordinates
      call drdy_trans(1,n3,amass,xr(1,iop),dx)
C move center of mass to origin
      call drdy_center(natom,n3,amass,xr(1,iop),scr1)
      if (lgs(6).le.2) then
         iop = 2
         write(fu6,625) iop,(iatsv(i,iop),i=1,nratom(iop))
         write (fu6,626) icode(iop)
      endif
C
      iop = 3
      write(fu6,621)
      write(fu6,623)
      write(fu6,624) (i,(xr(3*i-3+j,iop),j=1,3),i=1,natom)
      write(fu6,625) iop,(iatsv(i,iop),i=1,nratom(iop))
      write (fu6,626) icode(iop)
C transform to mass-scaled coordinates
      call drdy_trans(1,n3,amass,xr(1,iop),dx)
C move center of mass to origin
      call drdy_center(natom,n3,amass,xr(1,iop),scr1)
      if (lgs(6).eq.1.or.lgs(6).eq.3) then
         iop = 4
         write(fu6,625) iop,(iatsv(i,iop),i=1,nratom(iop))
         write (fu6,626) icode(iop)
      endif
C Input for saddle point      
      iop = 5
      write(fu6,622)
      write(fu6,623)
      write(fu6,624) (i,(xr(3*i-3+j,iop),j=1,3),i=1,natom)
      write(fu6,626) icode(iop)
c transform to mass-scaled coordinates
      call drdy_trans(1,n3,amass,xr(1,iop),dx)
C move center of mass to origin
      call drdy_center(natom,n3,amass,xr(1,iop),scr1)
C rotate reactant coordinates to have prinicple axes same as saddle point
      call drdy_rotx(xr(1,5),xr,scr1,natom,n3)
C rotate product coordinates to have prinicple axes same as saddle point
      call drdy_rotx(xr(1,5),xr(1,3),scr1,natom,n3)
      write(fu6,627)
      write(fu6,628)
      write(fu6,624) (i,(xr(3*i-3+j,1),j=1,3),i=1,natom)
      write(fu6,629)
      write(fu6,624) (i,(xr(3*i-3+j,3),j=1,3),i=1,natom)
      write(fu6,631)
      write(fu6,624) (i,(xr(3*i-3+j,5),j=1,3),i=1,natom)
c Input for reaction path calculation
      write (fu6,630) dels,delsv,slp,slm
c
      return
500   format(80a1)
501   Format(2i5,4f10.6,i10,g10.3)
600   Format(
     &      //,' DRDY_NWChem - DiRect DYnamics with the NWChem elect',
     &      'ronic structure code',//,5x,'This code prepares the fil',
     &      'e30 input for POLYRATE from electronic structure',/,5x,
     &      'calculations of gradients hessians at the reactant, pro',
     &      'duct, and saddle ',/,5x,'point geometries, and along t',
     &      'he minimum energy path.  Cartesian geometries',/,5x,'f',
     &      'or the reactants, products, and saddle points need to b',
     &      'e input to this ',/,5x,'code; optimization for these g',
     &      'eometries is not performed in this code.  ',/,5x,'Poin',
     &      'ts along the minimum energy path are calculated here.',
     &      //,' Options, lgs=',/,2x,39i2,
     &      //,10x,' Note that only lgs 1,6,8,27,31 are used.')
601   Format(/,6x,'lgs(1)=',i2,', no saddle point; MEP started at',
     *  ' initial geometry along the gradient')
602   Format(/,6x,'lgs(1)=',i2,', saddle point; MEP started along',
     *  ' unbound eigenvector')
603   Format(6x,'lgs(6)=',i2,',',i3,', reactants,',i3,', products')
604   Format(6x,'lgs(8)=',i2,', no restart')
605   Format(6x,'lgs(8)=',i2,', restart')
606   Format(5x,'lgs(27)=',i2,', use original direction of unbound',
     * ' eigenvector at saddle point')
607   Format(5x,'lgs(27)=',i2,', flip direction of unbound',
     * ' eigenvector at saddle point')
608   Format(5x,'lgs(31)=',i2,', MEP calculation uses euler integrator')
609   Format(5x,'lgs(31)=',i2,', MEP calculation uses Page-McIver LQA', 
     *   ' algorithm')
610   Format(5x,'lgs(31)=',i2,', MEP calculation uses Page-McIver CLQA',
     *   ' algorithm')
611   Format(5x,'lgs(31)=',i2,', MEP calculation uses Page-McIver CUBE',
     *   ' algorithm')
612   Format(//,' Title:',/,1x,80a1)
613   Format(//,' Number of atoms =',i6,',',i6,' dimensions')
614   Format (//,' Masses',/,T35,'amu',13x,'a.u.',12x,'kg')
615   Format(5x,'Mass of atom',i3,T30,0pf10.6,f17.6,1pe15.6)
616   Format(5x,'Forward reduced mass',T30,0pf10.6,f17.6,1pe15.6)
620   Format(//,' Reactants')
621   Format(//,' Products')
622   Format(//,' Saddle point or starting geometry on MEP')
623   Format (/,5x,'Input geometry (cartesian coordinates in au)',/,6x,
     *   4hatom,12x,'x',13x,'y',13x,'z',/)
624   Format (6x,i3,5x,f12.5,2x,f12.5,2x,f12.5)
625   Format(/,5x,'Species',I2,' consists of atom(s):' ,T38,20i3,/,
     *   (T38,20i3))
626   Format(10x,'icode=',i5)
627   Format(/,' Mass scaled coordinates, rotated so principle axes',
     *   ' are aligned')
628   Format(/,5x,'Reactants')
629   Format(/,5x,'Products')
631   Format(/,5x,'Saddle point or starting geometry on MEP')
630   Format(//,' Parameters for MEP calculation',//,
     *   5x,'Step size, dels=',f10.6,/,
     *   5x,'Hessian save size, delsv =',f10.6,/,
     *   5x,'Path also stopped if s .gt.',f10.6,' or s .lt.',f10.6)
      end
