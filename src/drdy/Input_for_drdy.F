C $Id$
      Subroutine Input_for_drdy(rtdb)
C
C Options in lgs array
C   1 - is there a saddle point?
C       =0, no
C       =1, yes
C   2 - are frequencies and eigenvectors printed at save points along MEP?
C       =0, no
C       =1, yes
C   6 - what is the number of reactants and products?
C	=1,  2 reactants and 2 products
C	=2,  2 reactants and 1 product
C	=3,  1 reactant and 2 products
C	=4,  1 reactant and 1 product
C   8 - is this a job that is restarted from data on file 30?
C       =0, no
C       =1, yes (requires that file 30 exists)
C  27 - is direction of unbound eigenvector at saddle point switched?
C       =0, no 
C       =+1 or -1, yes, change sign of unbound eigenvector
C  31 - which integrator is used for MEP calculation?
C       =1, Euler integrator
C       =2, use Page-McIver LQA algorithm
C       =3, use Page-McIver CLQA algorithm
C       =4, use Page-McIver CUBE algorithm
C  38 - are electronic structure files save?
C       =0, no files saved
C       =1, yes, the movecs file is saved
C       =2, yes, the hessian file is saved
C       =3, yes, the movecs and hessian files are saved
C       =4, yes, the movecs and hessian files and the movecs from the single point calculation are saved 
C  39 - is a single point calculation performed at each save point?
C       =0, no
C	    =1, yes
C
C  All other lgs variables are not used in this version of the code
C
      Implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "drdyP.fh"
*::passed:
      Integer rtdb ! [input] RTDB handle
*::locals:
      logical nwc_start
      logical nwc_continue
      logical nwc_restart
      Integer i,ii,j,l,iop
      Logical  lexit
      character*32 theory
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
      redm = 1.0d0
      dels = 0.01d0
      delsv = 0.1d0
      delhss = 1.0d11
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
      if (lgs(2).eq.0) then
         write (fu6,661) lgs(2)
      else
         write (fu6,662) lgs(2)
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
      if(lgs(38).eq.0) then
         write (fu6,6110) lgs(38)
      else if (lgs(38).eq.1) then
         write (fu6,6111) lgs(38)
      else if (lgs(38).eq.2) then
         write (fu6,6112) lgs(38)
      else if (lgs(38).eq.3) then
         write (fu6,6113) lgs(38)
      else 
         write (fu6,6114) lgs(38)
      endif
      if(lgs(39).eq.0) then
         write (fu6,6115) lgs(39)
      else
         if (.not.rtdb_cget(rtdb,'drdy:sptheory',1,theory)) then
           write (fu6,6117)
           call errquit('Input_for_drdy: problem with singlepoint',
     &                   555, RTDB_ERR)
         endif
         write (fu6,6116) lgs(39)
      endif
      write (fu6,612) ctitle
      n3 = 3*natom 
      write (fu6,613) natom
      if (natom.gt.natoms) then
         lexit = .true.
         write (fu6,*) ' number of atoms is greater than max=',natoms
         call errquit('Input_for_drdy: too many atoms', 555, INPUT_ERR)
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
C convert redm to atomic units
      redm = redm/cau
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
      do iop = 1,4
         if (icode(iop).le.1) then
            icode(iop) = 1
            nf(iop) = 0
         elseif (icode(iop).eq.2) then
            write(fu6,6000) iop
            call errquit('This option is not supported',555, INPUT_ERR)
         elseif (icode(iop).eq.3) then
            nf(iop) = 3*nratom(iop) - 5
         elseif (icode(iop).ge.4) then
            nf(iop) = 3*nratom(iop) - 6
            icode(iop) = 4
         endif
      enddo
      iop = 5
      if (icode(iop).le.1) then
         icode(iop) = 1
         nf(iop) = 0
      elseif (icode(iop).eq.2) then
         write(fu6,6000) iop
         call errquit('This option is not supported',555, INPUT_ERR)
      elseif (icode(iop).eq.3) then
         nf(iop) = n3 - 6
      elseif (icode(iop).ge.4) then
         nf(iop) = n3 - 7
         icode(iop) = 4
      endif
c transform to mass-scaled coordinates
      call drdy_trans(1,n3,amass,xr(1,5),dx)
C move center of mass to origin
      call drdy_center(natom,n3,amass,xr(1,5),scr1)
      if (icode(5).eq.3) then
C if linear transition state then rotate saddle point to be along z axis
C  set up vector along z axis
         ii = 0
         do i = 1,natom
            x(ii+1) = 0.0d00
            x(ii+2) = 0.0d00
            x(ii+3) = i
            ii = ii + 3
         enddo
         call drdy_rotx(x,xr(1,5),scr1,natom,n3)
      endif
C rotate reactant coordinates to have prinicple axes same as saddle point
      call drdy_rotx(xr(1,5),xr,scr1,natom,n3)
C rotate product coordinates to have prinicple axes same as saddle point
      call drdy_rotx(xr(1,5),xr(1,3),scr1,natom,n3)
      write(fu6,627)
      write(fu6,628)
      write(fu6,624) (i,(xr(3*i-3+j,1),j=1,3),i=1,natom)
      write(fu6,629)
      write(fu6,624) (i,(xr(3*i-3+j,3),j=1,3),i=1,natom)
      write(fu6,630)
      write(fu6,624) (i,(xr(3*i-3+j,5),j=1,3),i=1,natom)
c Input for reaction path calculation
      if (delhss.ge.1.0d11) delhss = delsv
      write (fu6,631) dels,delsv,delhss,slp,slm
      if (lgs(31).eq.1 .and. delhss.lt.delsv) then
         delhss = delsv
         write (fu6,632) delhss
      endif
c
      return
501   Format(2i5,4f10.6,i10,g10.3)
600   Format(//,' DIRDYVTST - DiRect DYnamics for Variational',
     *  ' Transition State Theory with the NWChem ',
     *  ' electronic structure code',//,5x,
     *  'This code prepares the file30 input for POLYRATE from',
     *  ' electronic structure calculations of gradients',/,5x,
     *  'hessians at the reactant, product, and saddle point',
     *  ' geometries, and along the minimum energy path.',/,5x,
     *  'Cartesian geometries for the reactants, products, and saddle',
     *  ' points need to be input to this code;',/,5x,'optimization',
     *  ' for these geometries is not performed in the code.  Points',
     *  ' along the minimum energy path',/,5x,'are calculated here.',
     *  //,' Options, lgs=',/,2x,39i2,/,/10x,' Note that only lgs ',
     *  '1,2,6,8,27,31,38,39 are used.')
601   Format(/,5x,'lgs( 1)=',i2,', no saddle point; MEP started at',
     *  ' initial geometry along the gradient')
602   Format(/,5x,'lgs( 1)=',i2,', saddle point; MEP started along',
     *  ' unbound eigenvector')
661   Format(5x,'lgs( 2)=',i2,', frequencies and eigenvectors',
     *  ' at save points along MEP are not printed')
662   Format(5x,'lgs( 2)=',i2,', print frequencies and eigenvectors',
     *  ' at save points along MEP')
603   Format(5x,'lgs( 6)=',i2,',',i3,', reactants,',i3,', products')
604   Format(5x,'lgs( 8)=',i2,', no restart')
605   Format(5x,'lgs( 8)=',i2,', restart')
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
6110  Format(5x,'lgs(38)=',i2,', No electronic structure files saved')
6111  Format(5x,'lgs(38)=',i2,', Electronic structure output file',
     *   ' saved')
6112  Format(5x,'lgs(38)=',i2,', Electronic structure checkpoint file',
     *   ' saved')
6113  Format(5x,'lgs(38)=',i2,', Electronic structure output and',
     *   ' checkpoint files saved')
6114  Format(5x,'lgs(38)=',i2,', Electronic structure output, ',
     *   ' checkpoint, and single point calculation files saved')
6115  Format(5x,'lgs(39)=',i2,', No single point calculation ',
     *   ' performed')
6116  Format(5x,'lgs(39)=',i2,', Single point calculation performed')
6117  Format(5x,'ERROR: Single point calculation cannot be',
     *       ' performed since',/,5x,
     *       ' sptheory was not given in DIRDYVTST block')
612   Format(//,' Title:',/,1x,a80)
613   Format(//,' Number of atoms =',i6)
614   Format (//,' Masses',/,T35,'amu',13x,'a.u.',12x,'kg')
615   Format(5x,'Mass of atom',i3,T30,0pf10.6,f17.6,1pe15.6)
616   Format(5x,'Forward reduced mass',T30,0pf10.6,f17.6,1pe15.6)
620   Format(//,' Reactants')
621   Format(//,' Products')
622   Format(//,' Saddle point or starting geometry on MEP')
623   Format (/,5x,'Input geometry (cartesian coordinates in au)',/,6x,
     *   4hatom,12x,'x',13x,'y',13x,'z')
624   Format (6x,i3,5x,f12.5,2x,f12.5,2x,f12.5)
625   Format(/,5x,'Species',I2,' consists of atom(s):' ,T38,20i3,/,
     *   (T38,20i3))
626   Format(10x,'icode=',i5)
627   Format(/,' Mass scaled coordinates, rotated so principle axes',
     *   ' are aligned')
628   Format(/,5x,'Reactants',/,6x,'atom',12x,'x',13x,'y',13x,'z')
629   Format(/,5x,'Products',/,6x,'atom',12x,'x',13x,'y',13x,'z')
630   Format(/,5x,'Saddle point or starting geometry on MEP',/,6x,
     *   'atom',12x,'x',13x,'y',13x,'z')
631   Format(//,' Parameters for MEP calculation',//,
     *   5x,'Step size, dels=',f10.6,/,
     *   5x,'Interval for saving Hessian data, delsv =',f10.6,/,
     *   5x,'Interval for computing Hessian data, delhss =',f10.6,/,
     *   5x,'Path also stopped if s .gt.',f10.6,' or s .lt.',f10.6)
632   Format(/,5x,'For the Euler integrator delhss is reset to delsv',
     *   ' delhss=',f10.6)
6000  Format (' For iop=',i5,', icode=2, Morse diatomic, was chosen',
     *   /, 'This option is not supported, reset icode to 3 and try',
     *   ' again')
      end
