
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C         B E G I N   O F   P B C   P A R T 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c read coordinates in Angst and converts them to au 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine nwpwxc_pbcrdcoord(fname,lattice,n,xyz,iat,autoang)
      implicit none             
      interface
        subroutine nwpwxc_parse(str,delims,args,nargs)
        character(len=*),intent(inout) :: str
        character(len=*),intent(in)  :: delims
        character(len=*),dimension(:),intent(inout) :: args
        integer, intent(out) :: nargs
        end subroutine nwpwxc_parse
      end interface
       
      real*8                :: xyz(3,*)
      real*8, INTENT(OUT)   ::lattice(3,3)
      integer, INTENT(out)               :: iat(*) 
      integer, INTENT(in)               :: n 
      character*(*), INTENT(IN)          :: fname
      logical              :: selective=.FALSE. ! Selective dynamics
      logical              :: cartesian=.TRUE.  ! Cartesian or direct
      real*8, INTENT(IN)   ::autoang

      real*8 xx(10),scalar
      character*200 line
      character*80 args(90),args2(90)
      
      integer i,j,ich,nn,ntype,ntype2,atnum,i_dummy1,i_dummy2,ncheck


      lattice=0
      
      ich=142
      open(unit=ich,file=fname)
      rewind(ich)
      ncheck=0
      ntype=0
      read(ich,'(a)',end=200)line !first line must contain Element Info
      call nwpwxc_parse(line,' ',args,ntype)
      read(ich,'(a)',end=200)line !second line contains global scaling factor
      call nwpwxc_readl(line,xx,nn)
      scalar=xx(1)/autoang        !the Ang->au conversion is included in the scaling factor
c      write(*,'(F8.6)')scalar
      DO i=1,3            ! reading the lattice constants
        read(ich,'(a)',end=200)line
        call nwpwxc_readl(line,xx,nn)
        IF (nn < 3) 
     >     call nwpwxc_stoprun( 'Error reading unit cell vectors' )
        lattice(1,i)=xx(1)*scalar
        lattice(2,i)=xx(2)*scalar
        lattice(3,i)=xx(3)*scalar
      !  write(*,'(3F6.2)')lattice(1,i),lattice(2,i),lattice(3,i)
      ENDDO
      read(ich,'(a)',end=200)line !Ether here are the numbers of each element, or (>vasp.5.1) here are the element symbols
      line=adjustl(line)
      call nwpwxc_readl(line,xx,nn)
      IF (nn.eq.0) then      ! CONTCAR files have additional Element line here since vasp.5.1
        call nwpwxc_parse(line,' ',args,ntype)
        read(ich,'(a)',end=200)line
        line=adjustl(line)
        call nwpwxc_readl(line,xx,nn)
      ENDIF
!       call elem(args(1),i_dummy2)
!       IF (i_dummy2<1 .OR. i_dummy2>94) THEN
!          args=args2
!       ENDIF
      IF (nn.NE.ntype ) THEN
        call nwpwxc_stoprun( 'Error reading number of atomtypes')
      ENDIF
      ncheck=0
      DO i=1,nn
        i_dummy1=INT(xx(i))
        call nwpwxc_elem(args(i),i_dummy2)
        IF (i_dummy2<1 .OR. i_dummy2>94) 
     .   call nwpwxc_stoprun( 'Error: unknown element.')
        DO j=1,i_dummy1
          ncheck=ncheck+1
          iat(ncheck)=i_dummy2
        ENDDO
      ENDDO
      if (n.ne.ncheck) 
     >   call nwpwxc_stoprun('Error reading Number of Atoms')

      read(ich,'(a)',end=200)line
      line=adjustl(line)
      IF (line(:1).EQ.'s' .OR. line(:1).EQ.'S') THEN
        selective=.TRUE.
        read(ich,'(a)',end=200)line
        line=adjustl(line)
      ENDIF

c      write(*,*)line(:1)
      cartesian=(line(:1).EQ.'c' .OR. line(:1).EQ.'C' .OR. 
     .line(:1).EQ.'k' .OR. line(:1).EQ.'K')
      DO i=1,n
        read(ich,'(a)',end=200)line
        call nwpwxc_readl(line,xx,nn)
        IF (nn.NE.3) call nwpwxc_stoprun( 'Error reading coordinates.')

        IF (cartesian) THEN
          xyz(1,i)=xx(1)*scalar
          xyz(2,i)=xx(2)*scalar
          xyz(3,i)=xx(3)*scalar
        ELSE
          xyz(1,i)=lattice(1,1)*xx(1)+lattice(1,2)*
     .    xx(2)+lattice(1,3)*xx(3)
          xyz(2,i)=lattice(2,1)*xx(1)+lattice(2,2)*xx(2)+lattice(2,3)*
     .    xx(3)
          xyz(3,i)=lattice(3,1)*xx(1)+lattice(3,2)*xx(2)+lattice(3,3)*
     .    xx(3)
        ENDIF
        
c      write(*,'(3F20.10,1X,I3)')xyz(:,i),iat(i)   !debug printout
      
      ENDDO
      
      
 200  continue

      close(ich)
      end subroutine nwpwxc_pbcrdcoord


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C compute coordination numbers by adding an inverse damping function
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine nwpwxc_pbcncoord(natoms,rcov,iz,xyz,cn,
     >                            lat,rep_cn,crit_cn)
      implicit none  
      !include 'param'
      real*8 k1,k2,k3

c global ad hoc parameters
      parameter (k1=16.0)
      parameter (k2=4./3.) 

c reasonable choices are between 3 and 5
c this gives smoth curves with maxima around the integer values
c k3=3 give for CN=0 a slightly smaller value than computed
c for the free atom. This also yields to larger CN for atoms
c in larger molecules but with the same chem. environment
c which is physically not right
c values >5 might lead to bumps in the potential
      parameter (k3=-4.) 


      integer,intent(in) :: natoms,iz(*)
      real*8,intent(in)  :: rcov(94)

      integer i,max_elem,rep_cn(3)
      real*8 xyz(3,*),cn(*),lat(3,3)

      integer iat,taux,tauy,tauz    
      real*8 dx,dy,dz,r,damp,xn,rr,rco,tau(3)
      real*8, INTENT(IN) :: crit_cn

      do i=1,natoms
      xn=0.0d0
      do iat=1,natoms
        do taux=-rep_cn(1),rep_cn(1)
        do tauy=-rep_cn(2),rep_cn(2)
        do tauz=-rep_cn(3),rep_cn(3)
          if(iat.eq.i .and. taux.eq.0 .and. tauy.eq.0 .and. 
     .      tauz.eq.0)        cycle
          tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)
          dx=xyz(1,iat)-xyz(1,i)+tau(1)
          dy=xyz(2,iat)-xyz(2,i)+tau(2)
          dz=xyz(3,iat)-xyz(3,i)+tau(3)
          r=(dx*dx+dy*dy+dz*dz)
          if (r.gt.crit_cn) cycle
          r=sqrt(r)
c covalent distance in Bohr
          rco=rcov(iz(i))+rcov(iz(iat))
          rr=rco/r
c counting function exponential has a better long-range behavior than MHGs inverse damping
          damp=1.d0/(1.d0+exp(-k1*(rr-1.0d0)))
          xn=xn+damp
c            print '("cn(",I2,I2,"): ",E14.8)',i,iat,damp

        enddo !tauz
        enddo !tauy
        enddo !taux
      enddo !iat
      cn(i)=xn  
      enddo !i

      end subroutine nwpwxc_pbcncoord

      subroutine nwpwxc_pbcrdatomnumber(fname,n)
      implicit none             
      interface
        subroutine nwpwxc_parse(str,delims,args,nargs)
        character(len=*),intent(inout) :: str
        character(len=*),intent(in)  :: delims
        character(len=*),dimension(:),intent(inout) :: args
        integer, intent(out) :: nargs
        end subroutine nwpwxc_parse
      end interface
       
      integer, INTENT(out)               :: n 
      character*(*), INTENT(IN)          :: fname
      logical              :: selective=.FALSE. ! Selective dynamics
      logical              :: cartesian=.TRUE.  ! Cartesian or direct

      real*8 xx(10),scalar,fdum
      character*80 line,args(90),args2(90)
      
      integer i,j,ich,nn,ntype,ntype2,atnum,i_dummy1,i_dummy2
      
      ich=142
      open(unit=ich,file=fname)
      n=0
      ntype=0
      read(ich,'(a)',end=200)line !first line must contain Element Info
      call nwpwxc_parse(line,' ',args,ntype)
      read(ich,'(a)',end=200)line !second line contains global scaling factor
      call nwpwxc_readl(line,xx,nn)
c      write(*,'(F8.6)')scalar
      DO i=1,3            ! reading the lattice constants
        read(ich,'(a)',end=200)line
        call nwpwxc_readl(line,xx,nn)
        IF (nn < 3) 
     >     call nwpwxc_stoprun( 'Error reading unit cell vectors' )
      !  write(*,'(3F6.2)')lattice(1,i),lattice(2,i),lattice(3,i)
      ENDDO
      read(ich,'(a)',end=200)line !Ether here are the numbers of each element, or (>vasp.5.1) here are the element symbols
      line=adjustl(line)
      call nwpwxc_readl(line,xx,nn)
      IF (nn.eq.0) then      ! CONTCAR files have additional Element line here since vasp.5.1
        call nwpwxc_parse(line,' ',args,ntype)
        read(ich,'(a)',end=200)line
        line=adjustl(line)
        call nwpwxc_readl(line,xx,nn)
      ENDIF
!       call elem(args(1),i_dummy2)
!       IF (i_dummy2<1 .OR. i_dummy2>94) THEN
!          args=args2
!       ENDIF
      IF (nn.NE.ntype ) THEN
!         IF(nn.NE.ntype2) THEN
        call nwpwxc_stoprun( 'Error reading number of atomtypes')
!         ELSE
!           ntype=ntype2
!         ENDIF
      ENDIF
      n=0
      DO i=1,nn
        i_dummy1=INT(xx(i))
          n=n+i_dummy1
      ENDDO

 200  continue

      close(ich)
      end subroutine nwpwxc_pbcrdatomnumber



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C compute energy
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
      subroutine nwpwxc_pbcedisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,
     >           r2r4,r0ab,
     .           rcov,rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,
     .           e6,e8,e10,e12,e63,lat,rthr,rep_vdw,cn_thr,rep_cn)
      implicit none  
      integer max_elem,maxc
      real*8 r2r4(max_elem),rcov(max_elem)
      real*8 rs6,rs8,rs10,alp6,alp8,alp10
      real*8 rthr,cn_thr,crit_cn
      integer rep_vdw(3),rep_cn(3)
      integer n,iz(*),version,mxc(max_elem)
!      integer rep_v(3)=rep_vdw!,rep_cn(3)
      real*8 xyz(3,*),r0ab(max_elem,max_elem),lat(3,3)!,r2r4(*)
!      real*8 rs6,rs8,rs10,alp6,alp8,alp10,rcov(max_elem)
      real*8 c6ab(max_elem,max_elem,maxc,maxc,3)
      real*8 e6, e8, e10, e12, e63!,crit_vdw,crit_cn
      logical noabc
 
      integer iat,jat,kat
      real*8 r,r2,r6,r8,tmp,dx,dy,dz,c6,c8,c10,ang,rav,R0
      real*8 damp6,damp8,damp10,rr,thr,c9,r42,c12,r10,c14
      real*8 cn(n),rxyz(3),dxyz(3)
      real*8 r2ab(n*n),cc6ab(n*n),dmp(n*n),d2(3),t1,t2,t3,tau(3)
      integer nwpwxc_lin,ij,ik,jk
      integer taux,tauy,tauz,counter
      real*8 a1,a2  !BJ-parameter
      real*8 bj_dmp6,bj_dmp8
      real*8 tmp1,tmp2

#include "global.fh"
#include "msgtypesf.h"

      integer taskid,np,pcount

      np     = GA_Nnodes()
      taskid = GA_Nodeid()
      pcount = 0


      e6 =0
      e8 =0
      e10=0
      e12=0
      e63=0
      tau=(/0.0,0.0,0.0/)
      counter=0
      crit_cn=cn_thr
c Becke-Johnson parameters
      a1=rs6
      a2=rs8      
      cc6ab = 0.0d0
      


C DFT-D2
      if(version.eq.2)then


      do iat=1,n-1
        do jat=iat+1,n
          if (mod(pcount,np).eq.taskid) then
          c6=c6ab(iz(jat),iz(iat),1,1,1)
          do taux=-rep_vdw(1),rep_vdw(1)
          do tauy=-rep_vdw(2),rep_vdw(2)
          do tauz=-rep_vdw(3),rep_vdw(3)
            tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)
            dx=xyz(1,iat)-xyz(1,jat)+tau(1)
            dy=xyz(2,iat)-xyz(2,jat)+tau(2)
            dz=xyz(3,iat)-xyz(3,jat)+tau(3)
            r2=dx*dx+dy*dy+dz*dz
            if(r2.gt.rthr) cycle
            r=sqrt(r2)
            damp6=1./(1.+exp(-alp6*(r/(rs6*r0ab(iz(jat),iz(iat)))-1.)))
            r6=r2**3      
            e6 =e6+c6*damp6/r6
          enddo !taux
          enddo !tauy
          enddo !tauz
        end if 
        pcount = pcount + 1
        enddo
      enddo
      
      do iat=1,n
        if (mod(pcount,np).eq.taskid) then
        jat=iat
        c6=c6ab(iz(jat),iz(iat),1,1,1)
        do taux=-rep_vdw(1),rep_vdw(1)
        do tauy=-rep_vdw(2),rep_vdw(2)
        do tauz=-rep_vdw(3),rep_vdw(3)
          if (taux.eq.0 .and. tauy.eq.0 .and. tauz.eq.0) cycle
          tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)
          dx=tau(1)
          dy=tau(2)
          dz=tau(3)
          r2=dx*dx+dy*dy+dz*dz
           if(r2.gt.rthr) cycle
          r=sqrt(r2)
          damp6=1./(1.+exp(-alp6*(r/(rs6*r0ab(iz(jat),iz(iat)))-1.)))
          r6=r2**3      
          e6 =e6+c6*damp6/r6*0.50d0
        enddo
        enddo
        enddo
      end if 
      pcount = pcount + 1
      enddo !iat
      if (np.gt.1) call GA_DGOP(9+MSGDBL,e6,1,'+')
      
      

      else if ((version.eq.3).or.(version.eq.5)) then
C DFT-D3(zero-damping)

        call nwpwxc_pbcncoord(n,rcov,iz,xyz,cn,lat,rep_cn,crit_cn)

        do iat=1,n-1
          do jat=iat+1,n
          if (mod(pcount,np).eq.taskid) then
c get C6
            call nwpwxc_getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),
     .                                  cn(iat),cn(jat),c6)

            if(.not.noabc)then
              ij=nwpwxc_lin(jat,iat)
c store C6 for C9, calc as sqrt
              cc6ab(ij)=sqrt(c6)
            endif
            do taux=-rep_vdw(1),rep_vdw(1)
            do tauy=-rep_vdw(2),rep_vdw(2)
            do tauz=-rep_vdw(3),rep_vdw(3)
              tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)

              dx=xyz(1,iat)-xyz(1,jat)+tau(1)
              dy=xyz(2,iat)-xyz(2,jat)+tau(2)
              dz=xyz(3,iat)-xyz(3,jat)+tau(3)
              r2=dx*dx+dy*dy+dz*dz
c cutoff

              if(r2.gt.rthr) cycle
              r =sqrt(r2)
              R0=r0ab(iz(jat),iz(iat))
              rr=R0/r
c damping
              if(version.eq.3)then
C DFT-D3 zero-damp
                tmp=rs6*rr   
                damp6 =1.d0/( 1.d0+6.d0*tmp**alp6 )
                tmp=rs8*rr     
                damp8 =1.d0/( 1.d0+6.d0*tmp**alp8 )
              else
C DFT-D3M zero-damp
                tmp=(r/(rs6*R0))+rs8*R0
                damp6 =1.d0/( 1.d0+6.d0*tmp**(-alp6) )
                tmp=(r/R0)+rs8*R0
                damp8 =1.d0/( 1.d0+6.d0*tmp**(-alp8) )
              endif


              r6=r2**3      
              e6 =e6+damp6/r6* c6
c             write(*,*)'e6: ',c6*damp6/r6*autokcal

c stored in main as sqrt
              c8 =3.0d0*r2r4(iz(iat))*r2r4(iz(jat))*c6
              r8 =r6*r2

              e8 =e8+c8*damp8/r8

            enddo !tauz
            enddo !tauy
            enddo !taux
          end if 
          pcount = pcount + 1
          enddo !jat
        enddo !iat
      
        do iat=1,n
        if (mod(pcount,np).eq.taskid) then
          jat=iat
c get C6
          call nwpwxc_getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),
     .                                  cn(iat),cn(jat),c6)
         
          if(.not.noabc)then
            ij=nwpwxc_lin(jat,iat)
c store C6 for C9, calc as sqrt
          cc6ab(ij)=sqrt(c6)
          endif
          do taux=-rep_vdw(1),rep_vdw(1)
          do tauy=-rep_vdw(2),rep_vdw(2)
          do tauz=-rep_vdw(3),rep_vdw(3)
            if (taux.eq.0 .and. tauy.eq.0 .and. tauz.eq.0) cycle
            tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)

            dx=tau(1)
            dy=tau(2)
            dz=tau(3)
            r2=dx*dx+dy*dy+dz*dz
c cutoff
            if(r2.gt.rthr) cycle
            r =sqrt(r2)
            R0=r0ab(iz(jat),iz(iat))
            rr=R0/r

            if(version.eq.3)then
C DFT-D3 zero-damp
              tmp=rs6*rr   
              damp6 =1.d0/( 1.d0+6.d0*tmp**alp6 )
              tmp=rs8*rr     
              damp8 =1.d0/( 1.d0+6.d0*tmp**alp8 )
            else
C DFT-D3M zero-damp
              tmp=(r/(rs6*R0))+rs8*R0
              damp6 =1.d0/( 1.d0+6.d0*tmp**(-alp6) )
              tmp=(r/R0)+rs8*R0
              damp8 =1.d0/( 1.d0+6.d0*tmp**(-alp8) )
            endif


            r6=r2**3      

            e6 =e6+damp6/r6*0.50d0 *C6

c stored in main as sqrt
            c8 =3.0d0*r2r4(iz(iat))*r2r4(iz(jat)) *C6
            r8 =r6*r2

            e8 =e8+c8*damp8/r8*0.50d0
            counter=counter+1

          enddo !tauz
          enddo !tauy
          enddo !taux
        end if 
        pcount = pcount + 1
        enddo !iat
        if (np.gt.1) call GA_DGOP(9+MSGDBL,e6,1,'+')
        if (np.gt.1) call GA_DGOP(9+MSGDBL,e8,1,'+')
        if (np.gt.1) call GA_DGOP(9+MSGDBL,cc6ab,n*n,'+')
!      write(*,*)'counter(edisp): ',counter
      else if((version.eq.4).or.(version.eq.6)) then


C DFT-D3(BJ-damping)
        call nwpwxc_pbcncoord(n,rcov,iz,xyz,cn,lat,rep_cn,crit_cn)

        do iat=1,n
        if (mod(pcount,np).eq.taskid) then
          do jat=iat+1,n
c get C6
            call nwpwxc_getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),
     .                                  cn(iat),cn(jat),c6)

            rxyz=xyz(:,iat)-xyz(:,jat)
            r42=r2r4(iz(iat))*r2r4(iz(jat))
            bj_dmp6=(a1*dsqrt(3.0d0*r42)+a2)**6
            bj_dmp8=(a1*dsqrt(3.0d0*r42)+a2)**8

            if(.not.noabc)then
              ij=nwpwxc_lin(jat,iat)
c store C6 for C9, calc as sqrt
              cc6ab(ij)=sqrt(c6)
            endif
            do taux=-rep_vdw(1),rep_vdw(1)
            do tauy=-rep_vdw(2),rep_vdw(2)
            do tauz=-rep_vdw(3),rep_vdw(3)
              tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)
            
              dxyz=rxyz+tau

              r2=sum(dxyz*dxyz)
c cutoff
              if(r2.gt.rthr) cycle
              r =sqrt(r2)
              rr=r0ab(iz(jat),iz(iat))/r


              r6=r2**3      

              e6 =e6+c6/(r6+bj_dmp6)

c stored in main as sqrt
              c8 =3.0d0*c6*r42
              r8 =r6*r2

              e8 =e8+c8/(r8+bj_dmp8)

              counter=counter+1

            enddo !tauz
            enddo !tauy
            enddo !taux
          enddo !jat

! Now the self interaction
          jat=iat
c get C6
          call nwpwxc_getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),
     .                                  cn(iat),cn(jat),c6)
          r42=r2r4(iz(iat))*r2r4(iz(iat))
          bj_dmp6=(a1*dsqrt(3.0d0*r42)+a2)**6
          bj_dmp8=(a1*dsqrt(3.0d0*r42)+a2)**8
           
          if(.not.noabc)then
            ij=nwpwxc_lin(jat,iat)
c store C6 for C9, calc as sqrt
            cc6ab(ij)=dsqrt(c6)
          endif

          do taux=-rep_vdw(1),rep_vdw(1)
          do tauy=-rep_vdw(2),rep_vdw(2)
          do tauz=-rep_vdw(3),rep_vdw(3)
            if (taux.eq.0 .and. tauy.eq.0 .and. tauz.eq.0) cycle
            tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)

            r2=sum(tau*tau)
c cutoff
            if(r2.gt.rthr) cycle
            r =sqrt(r2)
            rr=r0ab(iz(jat),iz(iat))/r


            r6=r2**3      

            e6 =e6+c6/(r6+bj_dmp6)*0.50d0

c stored in main as sqrt
            c8 =3.0d0*c6*r42
            r8 =r6*r2

            e8 =e8+c8/(r8+bj_dmp8)*0.50d0
            counter=counter+1

          enddo !tauz
          enddo !tauy
          enddo !taux
        end if 
        pcount = pcount + 1
        enddo !iat
        if (np.gt.1) call GA_DGOP(9+MSGDBL,e6,1,'+')
        if (np.gt.1) call GA_DGOP(9+MSGDBL,e8,1,'+')
        if (np.gt.1) call GA_DGOP(9+MSGDBL,cc6ab,n*n,'+')


      endif !version


      if(noabc)return

C compute non-additive third-order energy using averaged C6
       call nwpwxc_pbcthreebody(max_elem,xyz,lat,n,iz,rep_cn,crit_cn,
     .                   cc6ab,r0ab,e63)

      end subroutine nwpwxc_pbcedisp


      SUBROUTINE nwpwxc_pbcthreebody(max_elem,xyz,lat,n,iz,repv,
     >            cnthr,cc6ab,
     .            r0ab,eabc)
      IMPLICIT NONE
      integer max_elem
      INTEGER         :: n,i,j,k,jtaux,jtauy,jtauz,iat,jat,kat
      INTEGER         :: ktaux,ktauy,ktauz,counter,ij,ik,jk,idum
      REAL*8          :: dx,dy,dz,rij2,rik2,rjk2,c9,rr0ij,rr0ik
      REAL*8          :: rr0jk,geomean,fdamp,rik,rjk,rij
      REAL*8          :: r0ij,r0ik,r0jk
      REAL*8,INTENT(OUT)::eabc
      REAL*8          :: tmp,tmp1,tmp2,tmp3,tmp4,ang
                                           
      REAL*8 ,DIMENSION(3,3),INTENT(IN)::lat
      REAL*8 ,DIMENSION(3,*),INTENT(IN) :: xyz
      INTEGER,DIMENSION(*),INTENT(IN)::iz
      REAL*8,DIMENSION(3):: jtau,ktau,jxyz,kxyz,ijvec,ikvec,jkvec,dumvec
      INTEGER,DIMENSION(3):: repv
      REAL*8,INTENT(IN) ::cnthr
      REAL*8,DIMENSION(n*n),INTENT(IN)::cc6ab
      REAL*8,DIMENSION(max_elem,max_elem),INTENT(IN):: r0ab
      REAL*8,PARAMETER::sr9=0.75d0    !reciprocal radii scaling parameter for damping function (s_r=4/3)
      REAL*8,PARAMETER::alp9=-16.0d0  !alpha saved with "-" sign
      INTEGER,EXTERNAL :: nwpwxc_lin
      REAL*8 :: abcthr
      INTEGER, DIMENSION(3) :: repmin,repmax
!      REAL*8   :: time1,time2

#include "global.fh"
#include "msgtypesf.h"

      integer taskid,np,pcount

      np     = GA_Nnodes()
      taskid = GA_Nodeid()
      pcount = 0

      counter=0
      eabc=0.0d0
      abcthr=cnthr
!      abcthr=1.0d99
!      write(*,*)'thr:',(abcthr)

!      call cpu_time(time1)

      do iat=3,n
        do jat=2,iat-1
          ijvec=xyz(:,jat)-xyz(:,iat)
          ij=nwpwxc_lin(iat,jat)
          r0ij=r0ab(iz(iat),iz(jat))
          do kat=1,jat-1
          if (mod(pcount,np).eq.taskid) then
            ik=nwpwxc_lin(iat,kat)
            jk=nwpwxc_lin(jat,kat)
            ikvec=xyz(:,kat)-xyz(:,iat)
            jkvec=xyz(:,kat)-xyz(:,jat)
            c9=-1.0d0*(cc6ab(ij)*cc6ab(ik)*cc6ab(jk))

            r0ik=r0ab(iz(iat),iz(kat))
            r0jk=r0ab(iz(jat),iz(kat))


            do jtaux=-repv(1),repv(1)
              repmin(1)=max(-repv(1),jtaux-repv(1))
              repmax(1)=min(repv(1),jtaux+repv(1))
            do jtauy=-repv(2),repv(2)
              repmin(2)=max(-repv(2),jtauy-repv(2))
              repmax(2)=min(repv(2),jtauy+repv(2))
            do jtauz=-repv(3),repv(3)
              repmin(3)=max(-repv(3),jtauz-repv(3))
              repmax(3)=min(repv(3),jtauz+repv(3))
              jtau=jtaux*lat(:,1)+jtauy*lat(:,2)+jtauz*lat(:,3)
              dumvec=ijvec+jtau
              dumvec=dumvec*dumvec
              rij2=SUM(dumvec)
              if(rij2.gt.abcthr)cycle

              rr0ij=DSQRT(rij2)/r0ij


              do ktaux=repmin(1),repmax(1)
              do ktauy=repmin(2),repmax(2)
              do ktauz=repmin(3),repmax(3)
                ktau=ktaux*lat(:,1)+ktauy*lat(:,2)+ktauz*lat(:,3)
                dumvec=ikvec+ktau
                dumvec=dumvec*dumvec
                rik2=SUM(dumvec)
                if(rik2.gt.abcthr)cycle
                rr0ik=DSQRT(rik2)/r0ik

                dumvec=jkvec+ktau-jtau
                rjk2=SUM(dumvec*dumvec)
                if(rjk2.gt.abcthr)cycle
                rr0jk=DSQRT(rjk2)/r0jk


                geomean=(rr0ij*rr0ik*rr0jk)**(1.0d0/3.0d0)
!               write(*,*)'geomean:',geomean
                fdamp=1./(1.+6.*(sr9*geomean)**alp9)  !alp9 is already saved with "-"
                tmp1 = (rij2+rjk2-rik2)
                tmp2 = (rij2+rik2-rjk2)
                tmp3 = (rik2+rjk2-rij2)
                tmp4=rij2*rjk2*rik2
                ang=(0.375d0*tmp1*tmp2*tmp3/tmp4+1.0d0)/tmp4**1.5d0
                
                eabc=eabc+ang*c9*fdamp

              ENDDO !ktauz
              ENDDO !ktauy
              ENDDO !ktaux

            ENDDO !jtauz
            ENDDO !jtauy
            ENDDO !jtaux

          end if 
          pcount = pcount + 1
          ENDDO !kat
        ENDDO !jat
      ENDDO !iat

      DO iat=2,n
        jat=iat
          ij=nwpwxc_lin(iat,jat)
          ijvec=0.0d0
          r0ij=r0ab(iz(iat),iz(jat))
        DO kat=1,iat-1
        if (mod(pcount,np).eq.taskid) then
          jk=nwpwxc_lin(jat,kat)
          ik=jk
          ikvec=xyz(:,kat)-xyz(:,iat)
          jkvec=ikvec
          c9=-(cc6ab(ij)*cc6ab(ik)*cc6ab(jk))

          r0ik=r0ab(iz(iat),iz(kat))
          r0jk=r0ab(iz(jat),iz(kat))
          do jtaux=-repv(1),repv(1)
            repmin(1)=max(-repv(1),jtaux-repv(1))
            repmax(1)=min(repv(1),jtaux+repv(1))
          do jtauy=-repv(2),repv(2)
            repmin(2)=max(-repv(2),jtauy-repv(2))
            repmax(2)=min(repv(2),jtauy+repv(2))
          do jtauz=-repv(3),repv(3)
            repmin(3)=max(-repv(3),jtauz-repv(3))
            repmax(3)=min(repv(3),jtauz+repv(3))
            IF (jtaux.eq.0 .and. jtauy.eq.0 .and. jtauz.eq.0) cycle
            jtau=jtaux*lat(:,1)+jtauy*lat(:,2)+jtauz*lat(:,3)
            dumvec=ijvec+jtau
            dumvec=dumvec*dumvec
            rij2=SUM(dumvec)
            if(rij2.gt.abcthr)cycle

            rr0ij=DSQRT(rij2)/r0ij
       
            do ktaux=repmin(1),repmax(1)
            do ktauy=repmin(2),repmax(2)
            do ktauz=repmin(3),repmax(3)
! every result * 0.5
              ktau=ktaux*lat(:,1)+ktauy*lat(:,2)+ktauz*lat(:,3)
              dumvec=ikvec+ktau
              dumvec=dumvec*dumvec
              rik2=SUM(dumvec)
              if(rik2.gt.abcthr)cycle
              rr0ik=DSQRT(rik2)/r0ik

              dumvec=jkvec+ktau-jtau
              dumvec=dumvec*dumvec
              rjk2=SUM(dumvec)
              if(rjk2.gt.abcthr)cycle
              rr0jk=DSQRT(rjk2)/r0jk


              geomean=(rr0ij*rr0ik*rr0jk)**(1./3.)
              fdamp=1./(1.+6.*(sr9*geomean)**alp9)
              tmp1 = (rij2+rjk2-rik2)
              tmp2 = (rij2+rik2-rjk2)
              tmp3 = (rik2+rjk2-rij2)
              tmp4=rij2*rjk2*rik2
              ang=(0.375d0*tmp1*tmp2*tmp3/tmp4+1.0d0)/tmp4**1.5d0

              eabc=eabc+c9*fdamp*ang/2.0
            ENDDO !ktauz
            ENDDO !ktauy
            ENDDO !ktaux
 
          ENDDO !jtauz
          ENDDO !jtauy
          ENDDO !jtaux
        end if 
        pcount = pcount + 1
        ENDDO !kat
      ENDDO !iat
      ! And now kat=jat, but cycling throug all imagecells without jtau=ktau. and jat>iat going though all cells
      ! But this counts only 1/2

      DO iat=2,n
        DO jat=1,iat-1
        if (mod(pcount,np).eq.taskid) then
          kat=jat
          ij=nwpwxc_lin(iat,jat)
          jk=nwpwxc_lin(jat,kat)
          ik=ij
          ikvec=xyz(:,kat)-xyz(:,iat)
          ijvec=ikvec
          jkvec=0.0d0
          c9=-(cc6ab(ij)*cc6ab(ik)*cc6ab(jk))

          r0ij=r0ab(iz(iat),iz(jat))
          r0ik=r0ij
          r0jk=r0ab(iz(jat),iz(kat))

            do jtaux=-repv(1),repv(1)
              repmin(1)=max(-repv(1),jtaux-repv(1))
              repmax(1)=min(repv(1),jtaux+repv(1))
            do jtauy=-repv(2),repv(2)
              repmin(2)=max(-repv(2),jtauy-repv(2))
              repmax(2)=min(repv(2),jtauy+repv(2))
            do jtauz=-repv(3),repv(3)
              repmin(3)=max(-repv(3),jtauz-repv(3))
              repmax(3)=min(repv(3),jtauz+repv(3))
            jtau=jtaux*lat(:,1)+jtauy*lat(:,2)+jtauz*lat(:,3)
            dumvec=ijvec+jtau
            dumvec=dumvec*dumvec
            rij2=SUM(dumvec)
            if(rij2.gt.abcthr)cycle

            rr0ij=DSQRT(rij2)/r0ij
       
            do ktaux=repmin(1),repmax(1)
            do ktauy=repmin(2),repmax(2)
            do ktauz=repmin(3),repmax(3)
! every result * 0.5
              IF (jtaux.eq.ktaux .and. jtauy.eq.ktauy 
     .            .and. jtauz.eq.ktauz) cycle
              ktau=ktaux*lat(:,1)+ktauy*lat(:,2)+ktauz*lat(:,3)
              dumvec=ikvec+ktau
              dumvec=dumvec*dumvec
              rik2=SUM(dumvec)
              if(rik2.gt.abcthr)cycle
              rr0ik=DSQRT(rik2)/r0ik

              dumvec=jkvec+ktau-jtau
              dumvec=dumvec*dumvec
              rjk2=SUM(dumvec)
              if(rjk2.gt.abcthr)cycle
              rr0jk=DSQRT(rjk2)/r0jk


              geomean=(rr0ij*rr0ik*rr0jk)**(1./3.)
              fdamp=1./(1.+6.*(sr9*geomean)**alp9)
              tmp1 = (rij2+rjk2-rik2)
              tmp2 = (rij2+rik2-rjk2)
              tmp3 = (rik2+rjk2-rij2)
              tmp4=rij2*rjk2*rik2
              ang=(0.375d0*tmp1*tmp2*tmp3/tmp4+1.0d0)/tmp4**1.5d0

              eabc=eabc+c9*fdamp*ang/2.0
            ENDDO !ktauz
            ENDDO !ktauy
            ENDDO !ktaux
 
          ENDDO !jtauz
          ENDDO !jtauy
          ENDDO !jtaux
        end if
        pcount = pcount + 1
        ENDDO !kat
      ENDDO !iat


! And finally the self interaction iat=jat=kat all 

      idum=0
      DO iat=1,n
      if (mod(pcount,np).eq.taskid) then
      jat=iat
      kat=iat
      ijvec=0.0d0
      ij=nwpwxc_lin(iat,iat)
      ik=ij
      jk=ij
      ikvec=ijvec
      jkvec=ikvec
          c9=-(cc6ab(ij)*cc6ab(ik)*cc6ab(jk))

          r0ij=r0ab(iz(iat),iz(iat))
          r0ik=r0ij
          r0jk=r0ij
          do jtaux=-repv(1),repv(1)
            repmin(1)=max(-repv(1),jtaux-repv(1))
            repmax(1)=min(repv(1),jtaux+repv(1))
          do jtauy=-repv(2),repv(2)
            repmin(2)=max(-repv(2),jtauy-repv(2))
            repmax(2)=min(repv(2),jtauy+repv(2))
          do jtauz=-repv(3),repv(3)
            repmin(3)=max(-repv(3),jtauz-repv(3))
            repmax(3)=min(repv(3),jtauz+repv(3))
          IF (jtaux.eq.0 .and. jtauy.eq.0 .and. jtauz.eq.0) cycle
          jtau=jtaux*lat(:,1)+jtauy*lat(:,2)+jtauz*lat(:,3)
          dumvec=jtau
          dumvec=dumvec*dumvec
          rij2=SUM(dumvec)
          if(rij2.gt.abcthr)cycle
          rr0ij=DSQRT(rij2)/r0ij

            do ktaux=repmin(1),repmax(1)
            do ktauy=repmin(2),repmax(2)
            do ktauz=repmin(3),repmax(3)
           IF ((ktaux.eq.0) .and.( ktauy.eq.0) .and.( ktauz.eq.0))cycle !IF iat and kat are the same then cycle
            IF ((ktaux.eq.jtaux) .and. (ktauy.eq.jtauy) 
     .         .and. (ktauz.eq.jtauz)) cycle      !If kat and jat are the same then cycle

! every result * 1/6 becaues every triple is counted twice due to the two loops jtau and ktau going from -repv to repv -> *1/2
! 
!plus 1/3 becaues every triple is three times in each unitcell
              ktau=ktaux*lat(:,1)+ktauy*lat(:,2)+ktauz*lat(:,3)
              dumvec=ktau
              dumvec=dumvec*dumvec
              rik2=SUM(dumvec)
              if(rik2.gt.abcthr)cycle
              rr0ik=DSQRT(rik2)/r0ik

              dumvec=jkvec+ktau-jtau
              dumvec=dumvec*dumvec
              rjk2=SUM(dumvec)
              if(rjk2.gt.abcthr)cycle
              rr0jk=DSQRT(rjk2)/r0jk


              geomean=(rr0ij*rr0ik*rr0jk)**(1./3.)
              fdamp=1./(1.+6.*(sr9*geomean)**alp9)
              tmp1 = (rij2+rjk2-rik2)
              tmp2 = (rij2+rik2-rjk2)
              tmp3 = (rik2+rjk2-rij2)
              tmp4=rij2*rjk2*rik2
              ang=(0.375d0*tmp1*tmp2*tmp3/tmp4+1.0d0)/tmp4**1.5d0

              eabc=eabc+c9*fdamp*ang/6.0d0
 
          ENDDO !ktauz
          ENDDO !ktauy
          ENDDO !ktaux
        ENDDO !jtauz
        ENDDO !jtauy
        ENDDO !jtaux


      end if
      pcount = pcount + 1
      ENDDO !iat
      if (np.gt.1) call GA_DGOP(9+MSGDBL,eabc,1,'+')

      END SUBROUTINE nwpwxc_pbcthreebody


c     Input Geometry sanity check for pbc (to avoid au/Angtstrom mixups) J.M. 
      subroutine nwpwxc_pbccheckrcov(n,iz,rcov,xyz,lat)
      implicit none
      logical check
      integer iz(*),n,i,j,taux,tauy,tauz
      real*8 rcov(94),dist,dx,dy,dz,thr,xyz(3,*),r,lat(3,3),tau(3)
      check=.false.
      do i=1,n-1
       do j=i+1,n
         thr=0.5*(rcov(iz(i))+rcov(iz(j)))
         do taux=-1,1
         do tauy=-1,1
         do tauz=-1,1
            tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)

            dx=xyz(1,i)-xyz(1,j)+tau(1)
            dy=xyz(2,i)-xyz(2,j)+tau(2)
            dz=xyz(3,i)-xyz(3,j)+tau(3)
 
           r=sqrt(dx*dx+dy*dy+dz*dz)
           if (r.lt.thr) then
             check=.true.
!      write(*,*)'short distance',i,'(',iz(i),') and ',j,'(',iz(j),'):'
!           write(*,*)r,' < ',thr
!             write(*,*)
           endif
         enddo !tauz
         enddo !tauy
         enddo !taux
       enddo !j
      enddo !i
      if (check) then
          write(*,*)'--------------------------------------------------'
          write(*,*)'!! SOME DISTANCES VERY SHORT. CHECK COORDINATES !!'
          write(*,*)'--------------------------------------------------'
      endif
      end subroutine nwpwxc_pbccheckrcov

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C compute gradient
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine nwpwxc_pbcgdisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,
     >                           r2r4,r0ab,
     .            rcov,s6,s18,rs6,rs8,rs10,alp6,alp8,alp10,noabc,num,
     .                 version,g,disp,gnorm,stress,lat,rep_v,rep_cn,
     .                 crit_vdw,echo,crit_cn)
      
      implicit none  
      !include  'param'
      real*8 k1,k2,k3

c global ad hoc parameters
      parameter (k1=16.0)
      parameter (k2=4./3.) 

c reasonable choices are between 3 and 5
c this gives smoth curves with maxima around the integer values
c k3=3 give for CN=0 a slightly smaller value than computed
c for the free atom. This also yields to larger CN for atoms
c in larger molecules but with the same chem. environment
c which is physically not right
c values >5 might lead to bumps in the potential
      parameter (k3=-4.) 

      integer n,iz(*),max_elem,maxc,version,mxc(max_elem)
      real*8 xyz(3,*),r0ab(max_elem,max_elem),r2r4(*)
      real*8 c6ab(max_elem,max_elem,maxc,maxc,3)
      real*8 g(3,*),s6,s18,rcov(max_elem)
      real*8 rs6,rs8,rs10,alp10,alp8,alp6        
      real*8 a1,a2 !BJ-parameters
      real*8 bj_dmp6,bj_dmp8 ! precalculated dampingterms
      logical noabc,num,echo
c coversion factors
      REAL*8, parameter ::autoang =0.52917726d0
      REAL*8, parameter ::autokcal=627.509541d0
      REAL*8, parameter ::autoev=27.211652d0   

      integer iat,jat,i,j,kat,my,ny,a,b,idum,tau2,drijsize
      real*8 R0,C6,alp,R42,disp,x1,y1,z1,x2,y2,z2,rr,e6abc,fdum  
      real*8 dx,dy,dz,r2,r,r4,r6,r8,r10,r12,t6,t8,t10,damp1
      real*8 damp6,damp8,damp9,e6,e8,e10,e12,gnorm,tmp1
      real*8 s10,s8,gC6(3),term,step,dispr,displ,r235,tmp2
      real*8 cn(n),gx1,gy1,gz1,gx2,gy2,gz2,rthr,testsum
      real*8,  DIMENSION(3,3) :: lat,stress,sigma,virialstress,lat_1
      real*8,  DIMENSION(3,3) :: gC6_stress
      integer, DIMENSION(3)   :: rep_v,rep_cn
      real*8 crit_vdw,crit_cn
      integer taux,tauy,tauz
      real*8, DIMENSION(3) :: tau,vec12,dxyz,dxyz0
      real*8,external  ::nwpwxc_volume
      real*8           ::outpr(3,3)
      real*8, DIMENSION(3,3):: outerprod

      real*8 rij(3),rik(3),rjk(3),r7,r9
      real*8 rik_dist,rjk_dist
      real*8 drik,drjk
      real*8 rcovij
      real*8 dc6,c6chk !d(C6ij)/d(r_ij)
      real*8 expterm,dcni
      real*8, allocatable,dimension(:,:,:,:) ::  drij  !d(E)/d(r_ij) derivative wrt. dist. iat-jat
      !real*8, allocatable,dimension(:,:,:,:) :: dcn    !dCN(iat)/d(r_ij) is equal to
      real*8 dcnn
                                                       !dCN(jat)/d(r_ij)     
      real*8 :: dc6_rest  
      integer,external :: nwpwxc_lin
      real*8,external :: nwpwxc_vectorsize
      real*8 vec(3),vec2(3),dummy
      real*8 dc6i(n)       !dC6i(iat) saves dE_dsp/dCN(iat) 
      real*8 dc6ij(n,n)    !dC6(iat,jat)/cCN(iat) in dc6ij(i,j) for ABC-grad
c      real*8 dc6_rest_sum(n*(n+1)/2)
      integer linij,linik,linjk
      real*8 abc(3,n)

      real*8 eabc
      real*8 gabc(3,n),glatabc(3,3) !threebody gradient
      real*8 sigma_abc(3,3)
      real*8 labc,rabc
      real*8 ,dimension(3) ::ijvec,ikvec,jkvec,jtau,ktau,dumvec
      integer jtaux,jtauy,jtauz,ktaux,ktauy,ktauz,mtaux,mtauy,mtauz
      integer,dimension(3) :: taumin,taumax
      integer mat,linim,linjm,linkm
      real*8 rij2,rik2,rjk2,c9,c6ij,c6ik,c6jk,geomean,geomean3
      real*8 rr0ij,rr0jk,rr0ik,dc6iji,dc6ijj
      real*8 :: sr9=0.75d0 !inverse of 4/3
      real*8, parameter :: alp9=-16.0d0
      real*8,DIMENSION(n*(n+1)) ::c6save
      real*8 abcthr,time1,time2,geomean2,r0av,dc9,dfdmp,dang,ang
      integer,dimension(3) ::repv,repmin,repmax,repmin2,repmax2

#include "global.fh"
#include "msgtypesf.h"

      integer taskid,np,pcount

      np     = GA_Nnodes()
      taskid = GA_Nodeid()
      pcount = 0


c R^2 cut-off 
      rthr=crit_vdw
      abcthr=crit_cn
!      write(*,*)'abcthr:', abcthr**(1./1.)
      sigma=0.0d0
      virialstress=0.0d0
      stress=0.0d0
      gabc=0.0d0
      glatabc=0.0d0
      
c      testsum=0.0d0

      if(echo)write(*,*) 

      if(num) then
      if (echo)
     .  write(*,*) 'doing numerical gradient O(N^3) ...'

      call nwpwxc_pbcedisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,
     .           rcov,rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,
     .           e6,e8,e10,e12,e6abc,lat,rthr,rep_v,crit_cn,rep_cn)
 

          disp=-s6*e6-s18*e8-e6abc

      step=2.d-5

      do i=1,n
        do j=1,3
          xyz(j,i)=xyz(j,i)+step        
          call nwpwxc_pbcedisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,
     >           r2r4,r0ab,
     .           rcov,rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,
     .           e6,e8,e10,e12,e6abc,lat,rthr,rep_v,crit_cn,rep_cn)
 
          dispr=-s6*e6-s18*e8-e6abc
          rabc=e6abc
          xyz(j,i)=xyz(j,i)-2*step      
          call nwpwxc_pbcedisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,
     >           r2r4,r0ab,
     .           rcov,rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,
     .           e6,e8,e10,e12,e6abc,lat,rthr,rep_v,crit_cn,rep_cn)
 
          displ=-s6*e6-s18*e8-e6abc
          labc=e6abc
          gabc(j,i)=0.5*(rabc-labc)/step
          g(j,i)=0.5*(dispr-displ)/step  
          xyz(j,i)=xyz(j,i)+step        
        enddo !jat
      enddo   !iat
      IF (echo) write(*,*)'Doing numerical stresstensor...'

      call nwpwxc_xyz_to_abc(xyz,abc,lat,n)
      step=2.d-5
      if (echo) write(*,*)'step: ',step
      do i=1,3
        do j=1,3
          lat(j,i)=lat(j,i)+step
          call nwpwxc_abc_to_xyz(abc,xyz,lat,n)
          !call edisp...dum1
          call nwpwxc_pbcedisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,
     >           r2r4,r0ab,
     .           rcov,rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,
     .           e6,e8,e10,e12,e6abc,lat,rthr,rep_v,crit_cn,rep_cn)

         dispr=-s6*e6-s18*e8-e6abc
         labc=e6abc


          lat(j,i)=lat(j,i)-2*step
          call nwpwxc_abc_to_xyz(abc,xyz,lat,n)
          !call edisp...dum2
          call nwpwxc_pbcedisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,
     >           r2r4,r0ab,
     .           rcov,rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,
     .           e6,e8,e10,e12,e6abc,lat,rthr,rep_v,crit_cn,rep_cn)

         displ=-s6*e6-s18*e8-e6abc
         rabc=e6abc
          stress(j,i)=(dispr-displ)/(step*2.0)
          glatabc(j,i)=(rabc-labc)/(step*2.0d0)

          lat(j,i)=lat(j,i)+step
          call nwpwxc_abc_to_xyz(abc,xyz,lat,n)
          
        enddo !j
      enddo !i

      sigma=0.0d0
      call nwpwxc_inv_cell(lat,lat_1)
      do a=1,3
        do b=1,3
           do my=1,3
              sigma(a,b)=sigma(a,b)-stress(a,my)*lat(b,my)
           enddo
        enddo !b
      enddo !a

      goto 999

      endif !num


      if(version.eq.2)then
      if(echo)write(*,*) 'doing analytical gradient D-old O(N^2) ...'
      disp=0
      stress=0.0d0
      do iat=1,n-1
         do jat=iat+1,n
         if (mod(pcount,np).eq.taskid) then
           R0=r0ab(iz(jat),iz(iat))*rs6
           c6=c6ab(iz(jat),iz(iat),1,1,1)*s6
           do taux=-rep_v(1),rep_v(1)
           do tauy=-rep_v(2),rep_v(2)
           do tauz=-rep_v(3),rep_v(3)
            tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)
              dxyz=xyz(:,iat)-xyz(:,jat)+tau            
            r2  =sum(dxyz*dxyz)
           if (.not.(r2.gt.rthr)) then
           !if(r2.gt.rthr) cycle
            r235=r2**3.5                       
            r   =dsqrt(r2)
            damp6=exp(-alp6*(r/R0-1.0d0))
            damp1=1.+damp6           
            tmp1=damp6/(damp1*damp1*r235*R0)
            tmp2=6./(damp1*r*r235)

            term=alp6*tmp1-tmp2
              g(:,iat)=g(:,iat)-term*dxyz*c6
              g(:,jat)=g(:,jat)+term*dxyz*c6
            disp=disp+c6*(1./damp1)/r2**3

            do ny=1,3
            do my=1,3
              sigma(my,ny)=sigma(my,ny)+term*dxyz(ny)*dxyz(my)*c6
            enddo !my
            enddo !ny

            end if

           enddo !tauz
           enddo !tauy
           enddo !taux
         end if 
         pcount = pcount + 1
         enddo !jat
      enddo !iat
c and now the self interaction, only for convenient energy in dispersion
      do iat=1,n
      if (mod(pcount,np).eq.taskid) then
         jat=iat
           R0=r0ab(iz(jat),iz(iat))*rs6
           c6=c6ab(iz(jat),iz(iat),1,1,1)*s6
           do taux=-rep_v(1),rep_v(1)
           do tauy=-rep_v(2),rep_v(2)
           do tauz=-rep_v(3),rep_v(3)
            if (.not.(taux.eq.0 .and. tauy.eq.0 .and. tauz.eq.0)) then
            !if (taux.eq.0 .and. tauy.eq.0 .and. tauz.eq.0) cycle
            tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)

            dxyz=tau
!             vec12=(/ dx,dy,dz /)
            r2  =sum(dxyz*dxyz)
            if (.not.(r2.gt.rthr)) then
            !if(r2.gt.rthr) cycle
            r235=r2**3.5                       
            r   =dsqrt(r2)
            damp6=exp(-alp6*(r/R0-1.0d0))
            damp1=1.+damp6           
            tmp1=damp6/(damp1*damp1*r235*R0)
            tmp2=6./(damp1*r*r235)
            disp=disp+(c6*(1./damp1)/r2**3)*0.50d0
            term=alp6*tmp1-tmp2
            do ny=1,3
            do my=1,3
             sigma(my,ny)=sigma(my,ny)+term*dxyz(ny)*dxyz(my)*c6*0.5d0
            enddo !my
            enddo !ny
            
           end if
           end if

           enddo !tauz
           enddo !tauy
           enddo !taux
      end if 
      pcount = pcount + 1
      enddo !iat
      if (np.gt.1) call GA_DGOP(9+MSGDBL,disp,1,'+')
      if (np.gt.1) call GA_DGOP(9+MSGDBL,sigma,9,'+')
      if (np.gt.1) call GA_DGOP(9+MSGDBL,g,3*n,'+')
      
      call nwpwxc_inv_cell(lat,lat_1)
      do a=1,3
        do b=1,3
           do my=1,3
              stress(a,b)=stress(a,b)-sigma(a,my)*lat_1(b,my)
           enddo
        enddo !b
      enddo !a

      disp=-disp
!       sigma=virialstress
      goto 999
      endif !version==2

      if ((version.eq.3).or.(version.eq.5)) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          
!    begin ZERO DAMPING GRADIENT         
!          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (echo)   
     . write(*,*) 'doing analytical gradient O(N^2) ...'
c precompute for analytical part
      call nwpwxc_pbcncoord(n,rcov,iz,xyz,cn,lat,rep_cn,crit_cn)


      s8 =s18
      s10=s18
      allocate(drij(-rep_v(3):rep_v(3),-rep_v(2):rep_v(2),
     .              -rep_v(1):rep_v(1),n*(n+1)/2))

      drijsize =(2*rep_v(3)+1)*(2*rep_v(2)+1)*(2*rep_v(1)+1)*(n*(n+1)/2)
      !write(*,*) "n=",n
      !write(*,*) "rep_v=",rep_v
      !write(*,*) "drijsize1=",drijsize
      disp=0

      drij=0.0d0
      dc6_rest=0.0d0
c      dc6_rest_sum=0.0d0
      c6save=0.0d0
      kat=0
      dc6i=0.0d0


      do iat=1,n
      if (mod(pcount,np).eq.taskid) then  !* drij,disp,dc6i,dc6ij,c6save,dc6_rest_sum
        call nwpwxc_get_dC6_dCNij(maxc,max_elem,c6ab,mxc(iz(iat)),
     .          mxc(iz(iat)),cn(iat),cn(iat),iz(iat),iz(iat),iat,iat,
     .          c6,dc6iji,dc6ijj)

        c6save(nwpwxc_lin(iat,iat))=c6
        dc6ij(iat,iat)=dc6iji
        r0=r0ab(iz(iat),iz(iat))
        r42=r2r4(iz(iat))*r2r4(iz(iat))
        rcovij=rcov(iz(iat))+rcov(iz(iat))


        do taux=-rep_v(1),rep_v(1)
        do tauy=-rep_v(2),rep_v(2)
        do tauz=-rep_v(3),rep_v(3)
          tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)


!first dE/d(tau) saved in drij(i,i,counter)
          rij=tau
          r2=sum(rij*rij)
!          if (r2.gt.rthr) cycle

          if (r2.gt.0.1.and.r2.lt.rthr) then


          r=dsqrt(r2)
          r6=r2*r2*r2
          r7=r6*r
          r8=r6*r2
          r9=r8*r

!
!  Calculates damping functions:
C DGAS Edit
        if (version.eq.3) then 
          t6 = (r/(rs6*R0))**(-alp6)
          damp6 =1.d0/( 1.d0+6.d0*t6 )
          t8 = (r/(rs8*R0))**(-alp8)
          damp8 =1.d0/( 1.d0+6.d0*t8 )

          drij(tauz,tauy,taux,nwpwxc_lin(iat,iat))
     >     =drij(tauz,tauy,taux,nwpwxc_lin(iat,iat))
     .        +(-s6*(6.0/(r7)*C6*damp6)  ! d(r^(-6))/d(tau)
     .        -s8*(24.0/(r9)*C6*r42*damp8))*0.5d0


          drij(tauz,tauy,taux,nwpwxc_lin(iat,iat))
     >    =drij(tauz,tauy,taux,nwpwxc_lin(iat,iat))
     .        +(s6*C6/r7*6.d0*alp6*t6*damp6*damp6     !d(f_dmp)/d(tau)
     .        +s8*C6*r42/r9*18.d0*alp8*t8*damp8*damp8)*0.5d0
        else !version.eq.5
          t6 = (r/(rs6*R0)+R0*rs8)**(-alp6)
          damp6 =1.d0/( 1.d0+6.d0*t6 )
          t8 = (r/(R0)+R0*rs8)**(-alp8)
          damp8 =1.d0/( 1.d0+6.d0*t8 )

          tmp1=s6*6.d0*damp6*C6/r7
          tmp2=s8*6.d0*C6*r42*damp8/r9
          drij(tauz,tauy,taux,nwpwxc_lin(iat,iat))
     >     =drij(tauz,tauy,taux,nwpwxc_lin(iat,
     .     iat)) - (tmp1  ! d(r^(-6))/d(r_ij)
     .                +4.d0*tmp2)*0.5d0
      
      
          drij(tauz,tauy,taux,nwpwxc_lin(iat,iat))
     >     =drij(tauz,tauy,taux,nwpwxc_lin(iat,
     .     iat))
     .      +(tmp1*alp6*t6*damp6*r/(r+rs6*R0*R0*rs8)  !d(f_dmp)/d(r_ij)
     .      +3.d0*tmp2*alp8*t8*damp8*r/(r+R0*R0*rs8))*0.5d0
        endif
!
!      in dC6_rest all terms BUT C6-term is saved for the kat-loop
!          
          dc6_rest=
     .        (s6/r6*damp6+3.d0*s8*r42/r8*damp8)*0.50d0


          disp=disp-dc6_rest*c6  ! calculate E_disp for sanity check

          dc6i(iat)=dc6i(iat)+dc6_rest*(dc6iji+dc6ijj)
!          if (r2.lt.crit_cn)
c          dc6_rest_sum(nwpwxc_lin(iat,iat))
c     >    =dc6_rest_sum(nwpwxc_lin(iat,iat))+dc6_rest


          else !r2 < 0.1>rthr
             drij(tauz,tauy,taux,nwpwxc_lin(iat,iat))=0.0d0
          endif


        ENDDO !tauz
        ENDDO !tauy
        ENDDO !taux

!!!!!!!!!!!!!!!!!!!!!!!!!!
! B E G I N   jat  L O O P 
!!!!!!!!!!!!!!!!!!!!!!!!!!         
        do jat=1,iat-1
!
!      get_dC6_dCNij calculates the derivative dC6(iat,jat)/dCN(iat) and
!      dC6(iat,jat)/dCN(jat). these are saved in dC6ij for the kat loop
!
          call nwpwxc_get_dC6_dCNij(maxc,max_elem,c6ab,mxc(iz(iat)),
     .          mxc(iz(jat)),cn(iat),cn(jat),iz(iat),iz(jat),iat,jat,
     .          c6,dc6iji,dc6ijj)

          r0=r0ab(iz(jat),iz(iat))
          r42=r2r4(iz(iat))*r2r4(iz(jat))
          rcovij=rcov(iz(iat))+rcov(iz(jat))
          linij=nwpwxc_lin(iat,jat)
 
          dc6ij(iat,jat)=dc6iji
          dc6ij(jat,iat)=dc6ijj
          c6save(linij)=c6
            do taux=-rep_v(1),rep_v(1)
            do tauy=-rep_v(2),rep_v(2)
            do tauz=-rep_v(3),rep_v(3)
              tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)
  
  
            rij=xyz(:,jat)-xyz(:,iat)+tau
            r2=sum(rij*rij)
            if (.not.(r2.gt.rthr)) then
            !if (r2.gt.rthr) cycle
  
 
            r=dsqrt(r2)
            r6=r2*r2*r2
            r7=r6*r
            r8=r6*r2
            r9=r8*r
  
!
!  Calculates damping functions:
            if (version.eq.3) then
              t6 = (r/(rs6*R0))**(-alp6)
              damp6 =1.d0/( 1.d0+6.d0*t6 )
              t8 = (r/(rs8*R0))**(-alp8)
              damp8 =1.d0/( 1.d0+6.d0*t8 )
    
              drij(tauz,tauy,taux,linij)=drij(tauz,tauy,taux,
     .           linij)
     .          -s6*(6.0/(r7)*C6*damp6)  ! d(r^(-6))/d(r_ij)
     .          -s8*(24.0/(r9)*C6*r42*damp8)
  
              drij(tauz,tauy,taux,linij)=drij(tauz,tauy,taux,
     .           linij)
     .          +s6*C6/r7*6.d0*alp6*t6*damp6*damp6     !d(f_dmp)/d(r_ij)
     .          +s8*C6*r42/r9*18.d0*alp8*t8*damp8*damp8
            else !version.eq.5
              t6 = (r/(rs6*R0)+R0*rs8)**(-alp6)
              damp6 =1.d0/( 1.d0+6.d0*t6 )
              t8 = (r/(R0)+R0*rs8)**(-alp8)
              damp8 =1.d0/( 1.d0+6.d0*t8 )

              tmp1=s6*6.d0*damp6*C6/r7
              tmp2=s8*6.d0*C6*r42*damp8/r9
              drij(tauz,tauy,taux,linij)=drij(tauz,tauy,taux,
     .           linij) - (tmp1  ! d(r^(-6))/d(r_ij)
     .                    +4.d0*tmp2)
      
      
              drij(tauz,tauy,taux,linij)=drij(tauz,tauy,taux,
     .           linij)
     .          +(tmp1*alp6*t6*damp6*r/(r+rs6*R0*R0*rs8)  !d(f_dmp)/d(r_ij)
     .          +3.d0*tmp2*alp8*t8*damp8*r/(r+R0*R0*rs8))
            endif
!
!      in dC6_rest all terms BUT C6-term is saved for the kat-loop
!          
              dc6_rest=
     .        (s6/r6*damp6+3.d0*s8*r42/r8*damp8)

 
              disp=disp-dc6_rest*c6  ! calculate E_disp for sanity check

              dc6i(iat)=dc6i(iat)+dc6_rest*dc6iji
              dc6i(jat)=dc6i(jat)+dc6_rest*dc6ijj
!            if (r2.lt.crit_cn) 
c              dc6_rest_sum(linij)=dc6_rest_sum(linij)
c     .          +dc6_rest 

            end if

            enddo !tauz
            enddo !tauy
            enddo !taux
    
          enddo !jat

        end if
        pcount = pcount + 1
        enddo !iat



      elseif ((version.eq.4).or.(version.eq.6)) then

!!!!!!!!!!!!!!!!!!!!!!!
! NOW THE BJ Gradient !
!!!!!!!!!!!!!!!!!!!!!!!

      if (echo) write(*,*) 'doing analytical gradient O(N^2) ...'
      call nwpwxc_pbcncoord(n,rcov,iz,xyz,cn,lat,rep_cn,crit_cn)

      a1 =rs6
      a2 =rs8
      s8 =s18

      allocate(drij(-rep_v(3):rep_v(3),-rep_v(2):rep_v(2),
     .              -rep_v(1):rep_v(1),n*(n+1)/2))
      drijsize =(2*rep_v(3)+1)*(2*rep_v(2)+1)*(2*rep_v(1)+1)*(n*(n+1)/2)
      !write(*,*) "drijsize2=",drijsize
      disp=0
      drij=0.0d0
      dc6_rest=0.0d0
c      dc6_rest_sum=0.0d0
      kat=0

      do iat=1,n
      if (mod(pcount,np).eq.taskid) then !* drij,disp,dc6i,dc6ij,c6save,dc6_rest_sum
        call nwpwxc_get_dC6_dCNij(maxc,max_elem,c6ab,mxc(iz(iat)),
     .          mxc(iz(iat)),cn(iat),cn(iat),iz(iat),iz(iat),iat,iat,
     .          c6,dc6iji,dc6ijj)

        dc6ij(iat,iat)=dc6iji
        c6save(nwpwxc_lin(iat,iat))=c6
        r42=r2r4(iz(iat))*r2r4(iz(iat))
        rcovij=rcov(iz(iat))+rcov(iz(iat))

        R0=a1*sqrt(3.0d0*r42)+a2

        do taux=-rep_v(1),rep_v(1)
        do tauy=-rep_v(2),rep_v(2)
        do tauz=-rep_v(3),rep_v(3)
          tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)

!first dE/d(tau) saved in drij(i,i,counter)
          rij=tau
          r2=sum(rij*rij)
!          if (r2.gt.rthr) cycle

!          if (r2.gt.0.1) then
          if (r2.gt.0.1.and.r2.lt.rthr) then
!
!      get_dC6_dCNij calculates the derivative dC6(iat,jat)/dCN(iat) and
!      dC6(iat,jat)/dCN(jat). these are saved in dC6ij for the kat loop
!
          r=dsqrt(r2)
          r4=r2*r2
          r6=r4*r2
          r7=r6*r
          r8=r6*r2
          r9=r8*r

!
!  Calculates damping functions:
 
          t6=(r6+R0**6)
          t8=(r8+R0**8)

          drij(tauz,tauy,taux,nwpwxc_lin(iat,iat))
     >    =drij(tauz,tauy,taux,nwpwxc_lin(iat,
     .     iat))
     .        -s6*C6*6.0d0*r4*r/(t6*t6)*0.5d0 ! d(1/(r^(6)+R0^6)/d(r)
     .        -s8*C6*24.0d0*r42*r7/(t8*t8)*0.5d0


!
!      in dC6_rest all terms BUT C6-term is saved for the kat-loop
!          
          dc6_rest=
     .        (s6/t6+3.d0*s8*r42/t8)*0.50d0


          disp=disp-dc6_rest*c6  ! calculate E_disp for sanity check

          dc6i(iat)=dc6i(iat)+dc6_rest*(dc6iji+dc6ijj)
!          if (r2.lt.crit_cn)
c          dc6_rest_sum(nwpwxc_lin(iat,iat))
c     >    =dc6_rest_sum(nwpwxc_lin(iat,iat))+dc6_rest


          else !r2 < 0.1>rthr
             drij(tauz,tauy,taux,nwpwxc_lin(iat,iat))=0.0d0
          endif


        ENDDO !tauz
        ENDDO !tauy
        ENDDO !taux

!!!!!!!!!!!!!!!!!!!!!!!!!!
! B E G I N   jat  L O O P 
!!!!!!!!!!!!!!!!!!!!!!!!!!         
        do jat=1,iat-1
!
!      get_dC6_dCNij calculates the derivative dC6(iat,jat)/dCN(iat) and
!      dC6(iat,jat)/dCN(jat). these are saved in dC6ij for the kat loop
!
          call nwpwxc_get_dC6_dCNij(maxc,max_elem,c6ab,mxc(iz(iat)),
     .          mxc(iz(jat)),cn(iat),cn(jat),iz(iat),iz(jat),iat,jat,
     .          c6,dc6iji,dc6ijj)

          r42=r2r4(iz(iat))*r2r4(iz(jat))
          rcovij=rcov(iz(iat))+rcov(iz(jat))
 
          R0=a1*dsqrt(3.0d0*r42)+a2

          linij=nwpwxc_lin(iat,jat)
          dc6ij(iat,jat)=dc6iji
          dc6ij(jat,iat)=dc6ijj
          c6save(linij)=c6
            do taux=-rep_v(1),rep_v(1)
            do tauy=-rep_v(2),rep_v(2)
            do tauz=-rep_v(3),rep_v(3)
              tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)
  
  
            rij=xyz(:,jat)-xyz(:,iat)+tau
            r2=sum(rij*rij)

            if (.not.(r2.gt.rthr)) then
            !if (r2.gt.rthr) cycle
  
 
            r=dsqrt(r2)
            r4=r2*r2
            r6=r4*r2
            r7=r6*r
            r8=r6*r2
            r9=r8*r
  
!
!  Calculates damping functions:
            t6=(r6+R0**6)
            t8=(r8+R0**8)

 
            drij(tauz,tauy,taux,linij)=drij(tauz,tauy,taux,
     .           linij)
     .        -s6*C6*6.0d0*r4*r/(t6*t6)
     .        -s8*C6*24.0d0*r42*r7/(t8*t8)

!
!      in dC6_rest all terms BUT C6-term is saved for the kat-loop
!          
            dc6_rest=
     .        (s6/t6+3.d0*s8*r42/t8)

 
            disp=disp-dc6_rest*c6  ! calculate E_disp for sanity check

          dc6i(iat)=dc6i(iat)+dc6_rest*dc6iji
          dc6i(jat)=dc6i(jat)+dc6_rest*dc6ijj
!            if (r2.lt.crit_cn) 
c          dc6_rest_sum(nwpwxc_lin(iat,jat))=dc6_rest_sum(linij)
c     >             +dc6_rest 


            end if

          enddo !tauz
          enddo !tauy
          enddo !taux
 
         enddo !jat

      end if
      pcount = pcount + 1
      enddo !iat

      endif ! version=3 or 4
      if (np.gt.1) call GA_DGOP(9+MSGDBL,disp,1,'+')
      if (np.gt.1) call GA_DGOP(9+MSGDBL,c6save,n*(n+1)/2,'+')
c      if (np.gt.1) call GA_DGOP(9+MSGDBL,dc6_rest_sum,n*(n+1)/2,'+')
      if (np.gt.1) call GA_DGOP(9+MSGDBL,dc6ij,n*n,'+')


!
!!!!!!!!!!!!!!!!!!!!!!!
!!    BEGIN Threebody gradient
!!!!!!!!!!!!!!!!!!!!!!!      
      if (.not.noabc) then

!        write(*,*)'!!!!!!!!!!    THREEBODY  GRADIENT  !!!!!!!!!!'
        sr9=0.75d0
        eabc=0.0d0
        abcthr=crit_cn
        repv=rep_cn
!        write(*,*)'thr:',sqrt(abcthr)

        !call nwpwxc_cpu_time(time1)
        do iat=3,n
          do jat=2,iat-1
            linij=nwpwxc_lin(iat,jat)
            ijvec=xyz(:,jat)-xyz(:,iat)
          
            c6ij=c6save(linij)
            do kat=1,jat-1
            if (mod(pcount,np).eq.taskid) then !* using - c6save,  adding - drij,dc6i,
              linik=nwpwxc_lin(iat,kat)
              linjk=nwpwxc_lin(jat,kat)
              ikvec=xyz(:,kat)-xyz(:,iat)
              jkvec=xyz(:,kat)-xyz(:,jat)

              c6ik=c6save(linik)
              c6jk=c6save(linjk)
              c9=-1.0d0*dsqrt(c6ij*c6ik*c6jk)

              do jtaux=-rep_cn(1),rep_cn(1)
                repmin(1)=max(-rep_cn(1),jtaux-rep_cn(1))
                repmax(1)=min(rep_cn(1),jtaux+rep_cn(1))
              do jtauy=-rep_cn(2),rep_cn(2)
                repmin(2)=max(-rep_cn(2),jtauy-rep_cn(2))
                repmax(2)=min(rep_cn(2),jtauy+rep_cn(2))
              do jtauz=-rep_cn(3),rep_cn(3)
                repmin(3)=max(-rep_cn(3),jtauz-rep_cn(3))
                repmax(3)=min(rep_cn(3),jtauz+rep_cn(3))
                jtau=jtaux*lat(:,1)+jtauy*lat(:,2)+jtauz*lat(:,3)
                rij2=SUM((ijvec+jtau)*(ijvec+jtau))

                if (.not.(rij2.gt.abcthr)) then
                !if(rij2.gt.abcthr)cycle

                rr0ij=DSQRT(rij2)/r0ab(iz(iat),iz(jat))


                do ktaux=repmin(1),repmax(1)
                do ktauy=repmin(2),repmax(2)
                do ktauz=repmin(3),repmax(3)
                  ktau=ktaux*lat(:,1)+ktauy*lat(:,2)+ktauz*lat(:,3)
                  rik2=SUM((ikvec+ktau)*(ikvec+ktau))
                  if (.not.(rik2.gt.abcthr)) then
                  !if(rik2.gt.abcthr)cycle

                  dumvec=jkvec+ktau-jtau
                  rjk2=SUM(dumvec*dumvec)

                  if (.not.(rjk2.gt.abcthr)) then
                  !if(rjk2.gt.abcthr)cycle

                  rr0ik=dsqrt(rik2)/r0ab(iz(iat),iz(kat))
                  rr0jk=dsqrt(rjk2)/r0ab(iz(jat),iz(kat))
                  geomean2=(rij2*rjk2*rik2)
! first calculate the three components for the energy calculation fdmp
! and ang 
                  r0av=(rr0ij*rr0ik*rr0jk)**(1.0d0/3.0d0)
                  damp9=1./(1.+6.*(sr9*r0av)**alp9)  !alp9 is already saved with "-"

                  geomean=dsqrt(geomean2)
                  geomean3=geomean*geomean2
                  ang=0.375d0*(rij2+rjk2-rik2)*(rij2-rjk2+rik2)
     .                *(-rij2+rjk2+rik2)/(geomean3*geomean2)
     .                +1.0d0/(geomean3)

                  dc6_rest=ang*damp9
                  eabc=eabc+dc6_rest*c9
!
!start calculating the gradient components dfdmp, dang and dc9                  
                  
!dfdmp is the same for all three distances 
                 dfdmp=2.d0*alp9*(0.75d0*r0av)**(alp9)*damp9*damp9

!start calculating the derivatives of each part w.r.t. r_ij              
                  r=dsqrt(rij2)


                 dang=-0.375d0*(rij2**3+rij2**2*(rjk2+rik2)
     .               +rij2*(3.0d0*rjk2**2+2.0*rjk2*rik2+3.0*rik2**2)
     .               -5.0*(rjk2-rik2)**2*(rjk2+rik2))
     .               /(r*geomean3*geomean2)

                 tmp1=-dang*c9*damp9+dfdmp/r*c9*ang
                 drij(jtauz,jtauy,jtaux,linij)=
     .             drij(jtauz,jtauy,jtaux,linij)-tmp1
             
!start calculating the derivatives of each part w.r.t. r_ik              

                 r=dsqrt(rik2)


                 dang=-0.375d0*(rik2**3+rik2**2*(rjk2+rij2)
     .               +rik2*(3.0d0*rjk2**2+2.0*rjk2*rij2+3.0*rij2**2)
     .               -5.0*(rjk2-rij2)**2*(rjk2+rij2))
     .               /(r*geomean3*geomean2)

                 tmp1=-dang*c9*damp9+dfdmp/r*c9*ang
!                 tmp1=-dc9
                 drij(ktauz,ktauy,ktaux,linik)=
     .             drij(ktauz,ktauy,ktaux,linik)-tmp1

!
!start calculating the derivatives of each part w.r.t. r_jk              

                r=dsqrt(rjk2)

              dang=-0.375d0*(rjk2**3+rjk2**2*(rik2+rij2)
     .             +rjk2*(3.0d0*rik2**2+2.0*rik2*rij2+3.0*rij2**2)
     .             -5.0*(rik2-rij2)**2*(rik2+rij2))
     .               /(r*geomean3*geomean2)

              tmp1=-dang*c9*damp9+dfdmp/r*c9*ang
              drij(ktauz-jtauz,ktauy-jtauy,ktaux-jtaux,linjk)=
     .  drij(ktauz-jtauz,ktauy-jtauy,ktaux-jtaux,linjk)-tmp1

!calculating the CN derivative dE_disp(ijk)/dCN(i)

           dc9=dc6ij(iat,jat)/c6ij+dc6ij(iat,kat)/c6ik
              dc9=0.5d0*c9*dc9
              dc6i(iat)=dc6i(iat)+dc6_rest*dc9
             
           dc9=dc6ij(jat,iat)/c6ij+dc6ij(jat,kat)/c6jk
              dc9=0.5d0*c9*dc9
              dc6i(jat)=dc6i(jat)+dc6_rest*dc9

           dc9=dc6ij(kat,iat)/c6ik+dc6ij(kat,jat)/c6jk
              dc9=0.5d0*c9*dc9
              dc6i(kat)=dc6i(kat)+dc6_rest*dc9


                  end if
                  end if
                enddo !ktauz
                enddo !ktauy
                enddo !ktauz

                end if

              enddo !jtauz
              enddo !jtauy
              enddo !jtaux
            end if
            pcount = pcount + 1
            enddo !kat
          enddo !jat
        enddo !iat

! Now the interaction with jat=iat of the triples iat,iat,kat
      DO iat=2,n
        jat=iat
        linij=nwpwxc_lin(iat,jat)
        ijvec=0.0d0

        c6ij=c6save(linij)
        DO kat=1,iat-1
        if (mod(pcount,np).eq.taskid) then  !*** Adding-drij    Using- c6save
          linjk=nwpwxc_lin(jat,kat)
          linik=linjk

          c6ik=c6save(linik)
          c6jk=c6ik  
          ikvec=xyz(:,kat)-xyz(:,iat)
          jkvec=ikvec
          c9=-dsqrt(c6ij*c6ik*c6jk)
        do jtaux=-repv(1),repv(1)
          repmin(1)=max(-repv(1),jtaux-repv(1))
          repmax(1)=min(repv(1),jtaux+repv(1))
        do jtauy=-repv(2),repv(2)
          repmin(2)=max(-repv(2),jtauy-repv(2))
          repmax(2)=min(repv(2),jtauy+repv(2))
        do jtauz=-repv(3),repv(3)
          repmin(3)=max(-repv(3),jtauz-repv(3))
          repmax(3)=min(repv(3),jtauz+repv(3))

            IF(.not.(jtaux.eq.0 .and. jtauy.eq.0 .and. jtauz.eq.0)) then
            !IF (jtaux.eq.0 .and. jtauy.eq.0 .and. jtauz.eq.0) cycle

            jtau=jtaux*lat(:,1)+jtauy*lat(:,2)+jtauz*lat(:,3)
            dumvec=jtau
            rij2=SUM(dumvec*dumvec)

            if (.not.(rij2.gt.abcthr)) then
            !if(rij2.gt.abcthr)cycle

            rr0ij=DSQRT(rij2)/r0ab(iz(iat),iz(jat))
       
            do ktaux=repmin(1),repmax(1)
            do ktauy=repmin(2),repmax(2)
            do ktauz=repmin(3),repmax(3)
! every result * 0.5

              ktau=ktaux*lat(:,1)+ktauy*lat(:,2)+ktauz*lat(:,3)
              dumvec=ikvec+ktau
              dumvec=dumvec*dumvec
              rik2=SUM(dumvec)

              if (.not.(rik2.gt.abcthr)) then
              !if(rik2.gt.abcthr)cycle

              dumvec=jkvec+ktau-jtau
              dumvec=dumvec*dumvec
              rjk2=SUM(dumvec)

              if (.not.(rjk2.gt.abcthr)) then
              !if(rjk2.gt.abcthr)cycle

              rr0ik=DSQRT(rik2)/r0ab(iz(iat),iz(kat))
              rr0jk=DSQRT(rjk2)/r0ab(iz(jat),iz(kat))


              geomean2=(rij2*rjk2*rik2)
              r0av=(rr0ij*rr0ik*rr0jk)**(1.0d0/3.0d0)
              damp9=1./(1.+6.*(sr9*r0av)**alp9)  !alp9 is already saved with "-"

              geomean=dsqrt(geomean2)
              geomean3=geomean*geomean2
              ang=0.375d0*(rij2+rjk2-rik2)*(rij2-rjk2+rik2)
     .           *(-rij2+rjk2+rik2)/(geomean3*geomean2)
     .           +1.0d0/(geomean3)


              dc6_rest=ang*damp9/2.0d0   !factor 1/2 for doublecounting
              eabc=eabc+dc6_rest*c9

!              iat=jat
              dfdmp=2.d0*alp9*(0.75d0*r0av)**(alp9)*damp9*damp9

!start calculating the derivatives of each part w.r.t. r_ij              
              r=dsqrt(rij2)

              dang=-0.375d0*(rij2**3+rij2**2*(rjk2+rik2)
     .            +rij2*(3.0d0*rjk2**2+2.0*rjk2*rik2+3.0*rik2**2)
     .            -5.0*(rjk2-rik2)**2*(rjk2+rik2))
     .               /(r*geomean3*geomean2)

                 tmp1=-dang*c9*damp9+dfdmp/r*c9*ang
                 drij(jtauz,jtauy,jtaux,linij)=
     .             drij(jtauz,jtauy,jtaux,linij)-tmp1/2.0
             
!start calculating the derivatives of each part w.r.t. r_ik              
                 r=dsqrt(rik2)


                 dang=-0.375d0*(rik2**3+rik2**2*(rjk2+rij2)
     .               +rik2*(3.0d0*rjk2**2+2.0*rjk2*rij2+3.0*rij2**2)
     .               -5.0*(rjk2-rij2)**2*(rjk2+rij2))
     .               /(r*geomean3*geomean2)

                 tmp1=-dang*c9*damp9+dfdmp/r*c9*ang
                 drij(ktauz,ktauy,ktaux,linik)=
     .             drij(ktauz,ktauy,ktaux,linik)-tmp1/2.0
!
!start calculating the derivatives of each part w.r.t. r_ik              
              r=dsqrt(rjk2)

              dang=-0.375d0*(rjk2**3+rjk2**2*(rik2+rij2)
     .             +rjk2*(3.0d0*rik2**2+2.0*rik2*rij2+3.0*rij2**2)
     .             -5.0*(rik2-rij2)**2*(rik2+rij2))
     .               /(r*geomean3*geomean2)

              tmp1=-dang*c9*damp9+dfdmp/r*c9*ang

              drij(ktauz-jtauz,ktauy-jtauy,ktaux-jtaux,linjk)=
     .  drij(ktauz-jtauz,ktauy-jtauy,ktaux-jtaux,linjk)-tmp1/2.0

           dc9=dc6ij(iat,jat)/c6ij+dc6ij(iat,kat)/c6ik
              dc9=0.5d0*c9*dc9
              dc6i(iat)=dc6i(iat)+dc6_rest*dc9
             
           dc9=dc6ij(jat,iat)/c6ij+dc6ij(jat,kat)/c6jk
              dc9=0.5d0*c9*dc9
              dc6i(jat)=dc6i(jat)+dc6_rest*dc9

           dc9=dc6ij(kat,iat)/c6ik+dc6ij(kat,jat)/c6jk
              dc9=0.5d0*c9*dc9
              dc6i(kat)=dc6i(kat)+dc6_rest*dc9



              end if
              end if

            ENDDO !ktauz
            ENDDO !ktauy
            ENDDO !ktaux
 
            end if 
            end if 

          ENDDO !jtauz
          ENDDO !jtauy
          ENDDO !jtaux
        end if
        pcount = pcount + 1
        ENDDO !kat
      ENDDO !iat
      ! And now kat=jat, but cycling throug all imagecells without jtau=ktau. and jat>iat going though all cells    (iat,jat,jat)
      ! But this counts only 1/2

      DO iat=2,n
        DO jat=1,iat-1
        if (mod(pcount,np).eq.taskid) then   !*Adding - eabc,dc6i,drij
          kat=jat
          linij=nwpwxc_lin(iat,jat)
          linjk=nwpwxc_lin(jat,kat)
          linik=linij

          c6ij=c6save(linij)
          c6ik=c6ij

          c6jk=c6save(linjk)
          ikvec=xyz(:,kat)-xyz(:,iat)
          ijvec=ikvec
          jkvec=0.0d0

          c9=-1.0d0*dsqrt(c6ij*c6ik*c6jk)
          do jtaux=-repv(1),repv(1)
            repmin(1)=max(-repv(1),jtaux-repv(1))
            repmax(1)=min(repv(1),jtaux+repv(1))
          do jtauy=-repv(2),repv(2)
            repmin(2)=max(-repv(2),jtauy-repv(2))
            repmax(2)=min(repv(2),jtauy+repv(2))
          do jtauz=-repv(3),repv(3)
            repmin(3)=max(-repv(3),jtauz-repv(3))
            repmax(3)=min(repv(3),jtauz+repv(3))

            jtau=jtaux*lat(:,1)+jtauy*lat(:,2)+jtauz*lat(:,3)
            dumvec=ijvec+jtau
            dumvec=dumvec*dumvec
            rij2=SUM(dumvec)

            if (.not.(rij2.gt.abcthr)) then
            !if(rij2.gt.abcthr)cycle

            rr0ij=SQRT(rij2)/r0ab(iz(iat),iz(jat))
       
          do ktaux=repmin(1),repmax(1)
          do ktauy=repmin(2),repmax(2)
          do ktauz=repmin(3),repmax(3)
! every result * 0.5
              IF (.not.(jtaux.eq.ktaux .and. jtauy.eq.ktauy 
     .            .and. jtauz.eq.ktauz)) then
c              IF (jtaux.eq.ktaux .and. jtauy.eq.ktauy 
c     .            .and. jtauz.eq.ktauz) cycle

              ktau=ktaux*lat(:,1)+ktauy*lat(:,2)+ktauz*lat(:,3)
              dumvec=ikvec+ktau
              dumvec=dumvec*dumvec
              rik2=SUM(dumvec)

              if (.not.(rik2.gt.abcthr)) then
              !if(rik2.gt.abcthr)cycle

              rr0ik=SQRT(rik2)/r0ab(iz(iat),iz(kat))

              dumvec=jkvec+ktau-jtau
              dumvec=dumvec*dumvec
              rjk2=SUM(dumvec)

              if (.not.(rjk2.gt.abcthr)) then
              !if(rjk2.gt.abcthr)cycle

              rr0jk=SQRT(rjk2)/r0ab(iz(jat),iz(kat))

!              if (rij*rjk*rik.gt.abcthr)cycle

              geomean2=(rij2*rjk2*rik2)
              r0av=(rr0ij*rr0ik*rr0jk)**(1.0d0/3.0d0)
              damp9=1./(1.+6.d0*(sr9*r0av)**alp9)  !alp9 is already saved with "-"

              geomean=dsqrt(geomean2)
              geomean3=geomean*geomean2
              ang=0.375d0*(rij2+rjk2-rik2)*(rij2-rjk2+rik2)
     .           *(-rij2+rjk2+rik2)/(geomean2*geomean3)
     .           +1.0d0/(geomean3)
              dc6_rest=ang*damp9/2.0d0   !factor 1/2 for doublecounting
              eabc=eabc+dc6_rest*c9


!              jat=kat
              dfdmp=2.d0*alp9*(0.75d0*r0av)**(alp9)*damp9*damp9
!start calculating the derivatives of each part w.r.t. r_ij              
                 r=dsqrt(rij2)

                 dang=-0.375d0*(rij2**3+rij2**2*(rjk2+rik2)
     .               +rij2*(3.0d0*rjk2**2+2.0d0*rjk2*rik2+3.0d0*rik2**2)
     .               -5.0d0*(rjk2-rik2)**2*(rjk2+rik2))
     .               /(r*geomean3*geomean2)

                 tmp1=-dang*c9*damp9+dfdmp/r*c9*ang
                 drij(jtauz,jtauy,jtaux,linij)=
     .             drij(jtauz,jtauy,jtaux,linij)-tmp1/2.0d0
             
!start calculating the derivatives of each part w.r.t. r_ik              
                 r=dsqrt(rik2)


                 dang=-0.375d0*(rik2**3+rik2**2*(rjk2+rij2)
     .               +rik2*(3.0d0*rjk2**2+2.0*rjk2*rij2+3.0*rij2**2)
     .               -5.0*(rjk2-rij2)**2*(rjk2+rij2))
     .               /(r*geomean3*geomean2)

                 tmp1=-dang*c9*damp9+dfdmp/r*c9*ang
!                 tmp1=-dc9
                 drij(ktauz,ktauy,ktaux,linik)=
     .             drij(ktauz,ktauy,ktaux,linik)-tmp1/2.0d0
!
!start calculating the derivatives of each part w.r.t. r_jk              
              r=dsqrt(rjk2)

              dang=-0.375d0*(rjk2**3+rjk2**2*(rik2+rij2)
     .             +rjk2*(3.0d0*rik2**2+2.0*rik2*rij2+3.0*rij2**2)
     .             -5.0d0*(rik2-rij2)**2*(rik2+rij2))
     .               /(r*geomean3*geomean2)

              tmp1=-dang*c9*damp9+dfdmp/r*c9*ang
              drij(ktauz-jtauz,ktauy-jtauy,ktaux-jtaux,linjk)=
     .  drij(ktauz-jtauz,ktauy-jtauy,ktaux-jtaux,linjk)-tmp1/2.0d0

!calculating the CN derivative dE_disp(ijk)/dCN(i)

           dc9=dc6ij(iat,jat)/c6ij+dc6ij(iat,kat)/c6ik
              dc9=0.5d0*c9*dc9
              dc6i(iat)=dc6i(iat)+dc6_rest*dc9
             
           dc9=dc6ij(jat,iat)/c6ij+dc6ij(jat,kat)/c6jk
              dc9=0.5d0*c9*dc9
              dc6i(jat)=dc6i(jat)+dc6_rest*dc9

           dc9=dc6ij(kat,iat)/c6ik+dc6ij(kat,jat)/c6jk
              dc9=0.5d0*c9*dc9
              dc6i(kat)=dc6i(kat)+dc6_rest*dc9



              end if
              end if
              end if

            ENDDO !ktauz
            ENDDO !ktauy
            ENDDO !ktaux
 
            end if
          ENDDO !jtauz
          ENDDO !jtauy
          ENDDO !jtaux
        end if
        pcount = pcount + 1
        ENDDO !kat
      ENDDO !iat


! And finally the self interaction iat=jat=kat all 

      idum=0
      DO iat=1,n
      if (mod(pcount,np).eq.taskid) then !** adding drij,eabc,dc6i
      jat=iat
      kat=iat
      ijvec=0.0d0
      linij=nwpwxc_lin(iat,jat)
      linik=nwpwxc_lin(iat,kat)
      linjk=nwpwxc_lin(jat,kat)
      ikvec=ijvec
      jkvec=ikvec
          c6ij=c6save(linij)
          c6ik=c6ij
          c6jk=c6ij
          c9=-(DSQRT(c6ij*c6ij*c6ij))

        do jtaux=-repv(1),repv(1)
          repmin(1)=max(-repv(1),jtaux-repv(1))
          repmax(1)=min(repv(1),jtaux+repv(1))
        do jtauy=-repv(2),repv(2)
          repmin(2)=max(-repv(2),jtauy-repv(2))
          repmax(2)=min(repv(2),jtauy+repv(2))
        do jtauz=-repv(3),repv(3)
          repmin(3)=max(-repv(3),jtauz-repv(3))
          repmax(3)=min(repv(3),jtauz+repv(3))
         IF (.not.((jtaux.eq.0).and.(jtauy.eq.0).and.(jtauz.eq.0))) then
         !IF ((jtaux.eq.0) .and.(jtauy.eq.0) .and.(jtauz.eq.0))cycle !IF iat and jat are the same then cycle
          jtau=jtaux*lat(:,1)+jtauy*lat(:,2)+jtauz*lat(:,3)
          dumvec=jtau
          dumvec=dumvec*dumvec
          rij2=SUM(dumvec)

          if (.not.(rij2.gt.abcthr)) then
          !if(rij2.gt.abcthr)cycle

          rr0ij=SQRT(rij2)/r0ab(iz(iat),iz(jat))

          do ktaux=repmin(1),repmax(1)
          do ktauy=repmin(2),repmax(2)
          do ktauz=repmin(3),repmax(3)
           IF (.not.((ktaux.eq.0) .and.( ktauy.eq.0) 
     >         .and.( ktauz.eq.0))) then
           !IF ((ktaux.eq.0) .and.( ktauy.eq.0) .and.( ktauz.eq.0))cycle !IF iat and kat are the same then cycle

            IF (.not.((ktaux.eq.jtaux) .and. (ktauy.eq.jtauy) 
     .         .and. (ktauz.eq.jtauz))) then
c            IF ((ktaux.eq.jtaux) .and. (ktauy.eq.jtauy) 
c     .         .and. (ktauz.eq.jtauz)) cycle      !If kat and jat are the same then cycle

! every result * 1/6 becaues every triple is counted twice due to the two loops jtau and ktau going from -repv to repv -> *1/2
! 
!plus 1/3 becaues every triple is three times in each unitcell
              ktau=ktaux*lat(:,1)+ktauy*lat(:,2)+ktauz*lat(:,3)
              dumvec=ktau
              dumvec=dumvec*dumvec
              rik2=SUM(dumvec)

              if (.not.(rik2.gt.abcthr)) then
              !if(rik2.gt.abcthr)cycle

              rr0ik=SQRT(rik2)/r0ab(iz(iat),iz(kat))

              dumvec=jkvec+ktau-jtau
              dumvec=dumvec*dumvec
              rjk2=SUM(dumvec)

              if (.not.(rjk2.gt.abcthr)) then
              !if(rjk2.gt.abcthr)cycle

              rr0jk=SQRT(rjk2)/r0ab(iz(jat),iz(kat))

              geomean2=(rij2*rjk2*rik2)
              r0av=(rr0ij*rr0ik*rr0jk)**(1.0d0/3.0d0)
              damp9=1./(1.+6.*(sr9*r0av)**alp9)  !alp9 is already saved with "-"

              geomean=dsqrt(geomean2)
              geomean3=geomean*geomean2
              ang=0.375d0*(rij2+rjk2-rik2)*(rij2-rjk2+rik2)
     .            *(-rij2+rjk2+rik2)/(geomean2*geomean3)
     .            +1.0d0/(geomean3)
              dc6_rest=ang*damp9/6.0d0
              eabc=eabc+c9*dc6_rest
 
!                          iat=jat=kat
              dfdmp=2.d0*alp9*(0.75d0*r0av)**(alp9)*damp9*damp9
!start calculating the derivatives of each part w.r.t. r_ij              

                 r=dsqrt(rij2)
                 dang=-0.375d0*(rij2**3+rij2**2*(rjk2+rik2)
     .               +rij2*(3.0d0*rjk2**2+2.0*rjk2*rik2+3.0*rik2**2)
     .               -5.0*(rjk2-rik2)**2*(rjk2+rik2))
     .               /(r*geomean3*geomean2)


                 tmp1=-dang*c9*damp9+dfdmp/r*c9*ang
                 drij(jtauz,jtauy,jtaux,linij)=
     .             drij(jtauz,jtauy,jtaux,linij)-tmp1/6.0d0
             
!start calculating the derivatives of each part w.r.t. r_ik              

                 r=dsqrt(rik2)

                 dang=-0.375d0*(rik2**3+rik2**2*(rjk2+rij2)
     .               +rik2*(3.0d0*rjk2**2+2.0*rjk2*rij2+3.0*rij2**2)
     .               -5.0*(rjk2-rij2)**2*(rjk2+rij2))
     .               /(r*geomean3*geomean2)

                 tmp1=-dang*c9*damp9+dfdmp/r*c9*ang
                 drij(ktauz,ktauy,ktaux,linik)=
     .             drij(ktauz,ktauy,ktaux,linik)-tmp1/6.0d0
!
!start calculating the derivatives of each part w.r.t. r_jk

              r=dsqrt(rjk2)
              dang=-0.375d0*(rjk2**3+rjk2**2*(rik2+rij2)
     .             +rjk2*(3.0d0*rik2**2+2.0*rik2*rij2+3.0*rij2**2)
     .             -5.0*(rik2-rij2)**2*(rik2+rij2))
     .               /(r*geomean3*geomean2)

              tmp1=-dang*c9*damp9+dfdmp/r*c9*ang
              drij(ktauz-jtauz,ktauy-jtauy,ktaux-jtaux,linjk)=
     .  drij(ktauz-jtauz,ktauy-jtauy,ktaux-jtaux,linjk)-tmp1/6.0d0


!calculating the CN derivative dE_disp(ijk)/dCN(i)

           dc9=dc6ij(iat,jat)/c6ij+dc6ij(iat,kat)/c6ik
              dc9=0.5d0*c9*dc9
              dc6i(iat)=dc6i(iat)+dc6_rest*dc9
             
           dc9=dc6ij(jat,iat)/c6ij+dc6ij(jat,kat)/c6jk
              dc9=0.5d0*c9*dc9
              dc6i(jat)=dc6i(jat)+dc6_rest*dc9

           dc9=dc6ij(kat,iat)/c6ik+dc6ij(kat,jat)/c6jk
              dc9=0.5d0*c9*dc9
              dc6i(kat)=dc6i(kat)+dc6_rest*dc9


              end if
              end if
              end if
              end if

          ENDDO !ktauz
          ENDDO !ktauy
          ENDDO !ktaux

          end if 
          end if 

        ENDDO !jtauz
        ENDDO !jtauy
        ENDDO !jtaux


      end if
      pcount = pcount + 1
      ENDDO !iat
      if (np.gt.1) call GA_DGOP(9+MSGDBL,eabc,1,'+')


        !call nwpwxc_cpu_time(time2)

!        write(*,*)' eabc(gdisp):  ',eabc
!      write(*,'('' time(abc)  '',f6.1)')time2-time1
        disp=disp-eabc  ! cbannwarth: the minus is correct here, because the angular term is calculated with 0.375d0 instead of -0.375d0 (c.f., non-periodic code)
!        write(*,*)'gdisp:',disp
      endif !.not.noabc


451   continue

      if (np.gt.1) call GA_DGOP(9+MSGDBL,dc6i,n,'+')
      if (np.gt.1) call GA_DGOP(9+MSGDBL,drij,drijsize,'+')
      !write(*,*) "after noabc drijsize=",drijsize
      !write(*,*) "disp=",disp,n
      !write(*,*) "DC6I=",dc6i
      !write(*,*) "dc6ij=",dc6ij
      !write(*,*) "drij=",drij

      sigma_abc=0.0d0
      sigma=0.0d0

! After calculating all derivatives dE/dr_ij w.r.t. distances,
! the grad w.r.t. the coordinates is calculated dE/dr_ij * dr_ij/dxyz_i       
      do iat=2,n
        do jat=1,iat-1
        if (mod(pcount,np).eq.taskid) then
          linij=nwpwxc_lin(iat,jat)
!       write(*,'(3E17.6,XX,2I2)'),drij(0,0,-1:1,lin(iat,jat)),iat,jat
          rcovij=rcov(iz(iat))+rcov(iz(jat))
          do taux=-rep_v(1),rep_v(1)
          do tauy=-rep_v(2),rep_v(2)
          do tauz=-rep_v(3),rep_v(3)
            tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)

            rij=xyz(:,jat)-xyz(:,iat)+tau
            r2=sum(rij*rij)  
            if (.not.(r2.gt.rthr.or.r2.lt.0.5)) then
            !if (r2.gt.rthr.or.r2.lt.0.5) cycle
            r=dsqrt(r2)

            if (r2.lt.crit_cn) then
              expterm=exp(-k1*(rcovij/r-1.d0))
              dcnn=-k1*rcovij*expterm/
     .                 (r2*(expterm+1.d0)*(expterm+1.d0))
            else
              dcnn=0.0d0
            endif

            x1=drij(tauz,tauy,taux,linij)+dcnn*(dc6i(iat)+dc6i(jat))

            vec=x1*rij/r
            g(:,iat)=g(:,iat)+vec
            g(:,jat)=g(:,jat)-vec
            do i=1,3
            do j=1,3
              sigma(j,i)=sigma(j,i)+vec(j)*rij(i)
            enddo !j
            enddo !i



            end if
          enddo !tauz
          enddo !tauy
          enddo !taux
        end if
        pcount = pcount + 1
        enddo !jat
      enddo !iat
      if (np.gt.1) call GA_DGOP(9+MSGDBL,g,3*n,'+')

      do iat=1,n
      if (mod(pcount,np).eq.taskid) then
        rcovij=rcov(iz(iat))+rcov(iz(iat))
          do taux=-rep_v(1),rep_v(1)
          do tauy=-rep_v(2),rep_v(2)
          do tauz=-rep_v(3),rep_v(3)
            if (.not.(taux.eq.0.and.tauy.eq.0.and.tauz.eq.0)) then
            !if (taux.eq.0.and.tauy.eq.0.and.tauz.eq.0) cycle

            tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)
            r2=(sum(tau*tau))
            r=dsqrt(r2)
            if (r2.lt.crit_cn) then
              expterm=exp(-k1*(rcovij/r-1.d0))
              dcnn=-k1*rcovij*expterm/
     .                   (r2*(expterm+1.d0)*(expterm+1.d0))
            else
              dcnn=0.0d0
            endif
            x1=drij(tauz,tauy,taux,nwpwxc_lin(iat,iat))+dcnn*dc6i(iat)
            vec=x1*tau/r
            vec2(1)=taux
            vec2(2)=tauy
            vec2(3)=tauz
            do i=1,3
            do j=1,3
              sigma(j,i)=sigma(j,i)+vec(j)*tau(i)
            enddo !j
            enddo !i

            end if

          enddo !tauz
          enddo !tauy
          enddo !taux



      end if
      pcount = pcount + 1
      enddo !iat
      if (np.gt.1) call GA_DGOP(9+MSGDBL,sigma,9,'+')
      !write(*,*) "SIGMA=",sigma



      stress=0.0d0
      glatabc=0.0d0
      call nwpwxc_inv_cell(lat,lat_1)
      do a=1,3
        do b=1,3
           do my=1,3
              stress(a,b)=stress(a,b)-sigma(a,my)*lat_1(b,my)
           enddo
        enddo !b
      enddo !a



!          write(*,*)'drij:',drij(lin(iat,jat),:)
!          write(*,*)'g:',g(1,1:3)
!          write(*,*)'dcn:',sum(dcn(lin(2,1),:))


        !write(*,*) "deallocate drij"
        deallocate(drij)
        !deallocate(drij,dc6_rest,dcn)


      

 999  continue
!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!      
!This is where the D2 gradient and the numerical gradient jump.
!      
!!!!!!!!!!!!!!!!!!!!!!!!!!      
!      do i=1,n
!        write(*,'(83F17.12)') g(1:3,i)
!      enddo
      gnorm=sum(abs(g(1:3,1:n)))
      if(echo)then
c      write(*,*)'testsum:',testsum*autoev/autoang
      write(*,*)'|G(force)| =',gnorm
      gnorm=sum(abs(stress(1:3,1:3)))
      write(*,*)'|G(stress)|=',gnorm
      endif

      end subroutine nwpwxc_pbcgdisp


      subroutine nwpwxc_pbcwregrad(nat,g,g_lat)
      implicit none
      integer nat,i
      real*8 g(3,nat)
      real*8 g_lat(3,3)

      open(unit=142,file='dftd3_gradient')

!      write(*,*)'Gradient:' !Jonas
!      write(*,*)            !Jonas
      do i=1,nat
         write(142,'(3E22.14)')g(1:3,i)
!         write(*,'(3D22.14)')g(1:3,i) !Jonas
      enddo

      close(142)

      open(unit=143,file='dftd3_cellgradient')

!      write(*,*)'Gradient:' !Jonas
!      write(*,*)            !Jonas
      do i=1,3
         write(143,'(3E22.14)')g_lat(1:3,i)
!         write(*,'(3D22.14)')g(1:3,i) !Jonas
      enddo

      close(143)
      end subroutine nwpwxc_pbcwregrad

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C analyse all pairs
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine nwpwxc_pbcadisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,
     >           r2r4,r0ab,
     .           rcov,rs6,rs8,rs10,alp6,alp8,alp10,version,autokcal,
     .           autoang,rthr,rep_v,cn_thr,rep_cn,s6,s18,etot,lat)
      implicit none  
      integer n,iz(*),max_elem,maxc,version,mxc(max_elem) 
      real*8 xyz(3,*),r0ab(max_elem,max_elem),r2r4(*),s6
      real*8 rs6,rs8,rs10,alp6,alp8,alp10,autokcal,etot,s18,autoang
      real*8 c6ab(max_elem,max_elem,maxc,maxc,3),rcov(max_elem)
      real*8 lat(3,3)
      integer rep_v(3),rep_cn(3)
 
      integer iat,jat,i,j,k,nbin,taux,tauy,tauz
      real*8 R0,r,r2,r6,r8,tmp,alp,dx,dy,dz,c6,c8,c10
      real*8 damp6,damp8,damp10,r42,rr,check,rthr,cn_thr,rvdw
      real*8 cn(n),i6,e6,e8,e10,edisp                   
      real*8,allocatable ::  dist(:),li(:,:)
      real*8 xx(500),eg(10000)
      integer grplist(500,20)
      integer grpn(20),at(n)
      integer ngrp,dash
      integer nwpwxc_lin, iiii, jjjj, iii, jjj, ii, jj, ni, nj 
      integer iout(500)
      logical ex
      character*80 atmp
      real*8 tau(3)
 
      real*8,dimension(:,:), allocatable :: ed
      allocate(ed(n,n))


c distance bins
      nbin=17
      allocate(dist(0:nbin))
      allocate(li(0:nbin,2))

      li(0,1)=0   
      li(0,2)=1.5 
      li(1,1)=1.5
      li(1,2)=2
      li(2,1)=2
      li(2,2)=2.3333333333
      li(3,1)=2.3333333333
      li(3,2)=2.6666666666
      li(4,1)=2.6666666666
      li(4,2)=3.0
      li(5,1)=3.0          
      li(5,2)=3.3333333333
      li(6,1)=3.3333333333
      li(6,2)=3.6666666666
      li(7,1)=3.6666666666
      li(7,2)=4.0
      li(8,1)=4.0
      li(8,2)=4.5
      li(9,1)=4.5
      li(9,2)=5.0
      li(10,1)=5.0
      li(10,2)=5.5
      li(11,1)=5.5
      li(11,2)=6.0
      li(12,1)=6.0
      li(12,2)=7.0           
      li(13,1)=7.0           
      li(13,2)=8.0           
      li(14,1)=8.0           
      li(14,2)=9.0           
      li(15,1)=9.0           
      li(15,2)=10.0          
      li(16,1)=10.0           
      li(16,2)=20.0          
      li(17,1)=20.0           
      li(17,2)=dsqrt(rthr)*autoang


      call nwpwxc_pbcncoord(n,rcov,iz,xyz,cn,lat,rep_cn,cn_thr)

      write(*,*)
      write(*,*)'analysis of pair-wise terms (in kcal/mol)'
      write(*,'(''pair'',2x,''atoms'',9x,''C6'',14x,''C8'',12x,
     .''E6'',7x,''E8'',7x,''Edisp'')')
      e8=0
      ed=0
      dist=0
      check=0
      do iat=1,n
         do jat=iat,n

             do taux=-rep_v(1),rep_v(1)
             do tauy=-rep_v(2),rep_v(2)
             do tauz=-rep_v(3),rep_v(3)
              tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)
            dx=xyz(1,iat)-xyz(1,jat)+tau(1)
            dy=xyz(2,iat)-xyz(2,jat)+tau(2)
            dz=xyz(3,iat)-xyz(3,jat)+tau(3)
            r2=(dx*dx+dy*dy+dz*dz)
CTHR
            if(r2.gt.rthr.or.r2.lt.0.5) cycle
            r =sqrt(r2)
            R0=r0ab(iz(jat),iz(iat))
            rr=R0/r
            r6=r2**3       

            if(version.eq.3)then
              tmp=rs6*rr   
              damp6 =1.d0/( 1.d0+6.d0*tmp**alp6 )
              tmp=rs8*rr     
              damp8 =1.d0/( 1.d0+6.d0*tmp**alp8 )
            else
              tmp=(r/(R0*rs6)+R0*rs8)**(-alp6) 
              damp6 =1.d0/( 1.d0+6.d0*tmp )
              tmp=(r/(R0)+R0*rs8)**(-alp8) 
              damp8 =1.d0/( 1.d0+6.d0*tmp )
            endif

            if (version.eq.2)then
              c6=c6ab(iz(jat),iz(iat),1,1,1)
              damp6=1.d0/(1.d0+exp(-alp6*(r/(rs6*R0)-1.0d0)))
              if (iat.eq.jat) then
                e6 =s6*autokcal*c6*damp6/r6
              else
                e6 =s6*autokcal*c6*damp6/r6
              endif
              e8=0.0d0
            else
              call nwpwxc_getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),
     .                                      cn(iat),cn(jat),c6)
            endif

            if((version.eq.3).or.(version.eq.5))then
              r8 =r6*r2
              r42=r2r4(iz(iat))*r2r4(iz(jat))
              c8 =3.0d0*c6*r42
              if (iat.eq.jat) then
                e6 =s6*autokcal*c6*damp6/r6*0.5
                e8 =s18*autokcal*c8*damp8/r8*0.5
              else
                e6 =s6*autokcal*c6*damp6/r6
                e8 =s18*autokcal*c8*damp8/r8
              endif
            endif

            if((version.eq.4).or.(version.eq.6))then
              r42=r2r4(iz(iat))*r2r4(iz(jat))
              c8 =3.0d0*c6*r42
c use BJ radius
              R0=dsqrt(c8/c6)              
              rvdw=rs6*R0+rs8
              r8 =r6*r2
              if (iat.eq.jat) then
                e6 =s6*autokcal*c6/(r6+rvdw**6)*0.5
                e8 =s18*autokcal*c8/(r8+rvdw**8)*0.5
              else
                e6 =s6*autokcal*c6/(r6+rvdw**6)
                e8 =s18*autokcal*c8/(r8+rvdw**8)
              endif
            endif

            edisp=-(e6+e8)
            ed(iat,jat)=edisp
            ed(jat,iat)=edisp

!           write(*,'(2i4,2x,2i3,2D16.6,2F9.4,F10.5)')
!     .     iat,jat,iz(iat),iz(jat),c6,c8,
!     .    -e6,-e8,edisp

            check=check+edisp
            rr=r*autoang
            do i=0,nbin
               if(rr.gt.li(i,1).and.rr.le.li(i,2)) dist(i)=dist(i)+edisp
            enddo
          enddo !tauz
          enddo !tauy
          enddo !taux
         enddo !jat
      enddo !iat

      write(*,'(/''distance range (Angstroem) analysis'')')
      write(*,'( ''writing histogram data to <histo.dat>'')')
      open(unit=11,file='histo.dat')
      do i=0,nbin
         write(*,'(''R(low,high), Edisp, %tot :'',2f5.1,F12.5,F8.2)')
     .   li(i,1),li(i,2),dist(i),100.*dist(i)/etot
         write(11,*)(li(i,1)+li(i,2))*0.5,dist(i)
      enddo
      close(11)

      write(*,*) 'checksum (Edisp) ',check
      if(abs(check-etot).gt.1.d-3)stop'something is weired in adisp'

      deallocate(dist,li)
      return








      inquire(file='fragment',exist=ex)
      if(ex) return
      write(*,'(/''fragment based analysis'')')
      write(*,'( ''reading file <fragment> ...'')')
      open(unit=55,file='fragment')
      i=0
      at=0
 111  read(55,'(a)',end=222) atmp
      call nwpwxc_readfrag(atmp,iout,j)
      if(j.gt.0)then
         i=i+1
         grpn(i)=j
         do k=1,j
            grplist(k,i)=iout(k)      
            at(grplist(k,i))=at(grplist(k,i))+1
         enddo
      endif
      goto 111
 222  continue
      ngrp=i  
      k=0
      do i=1,n
         if(at(i).gt.1) stop 'something is weird in file <fragment>'
         if(at(i).eq.0)then
            k=k+1
            grplist(k,ngrp+1)=i
         endif
      enddo
      if(k.gt.0) then
         ngrp=ngrp+1
         grpn(ngrp)=k
      endif
c Implemented display of atom ranges instead of whole list of atoms
      write(*,*)'group #        atoms '
      dash=0
      do i=1,ngrp
       write(*,'(i4,3x,i4)',advance='no')i,grplist(1,i)
       do j=2,grpn(i)
        if(grplist(j,i).eq.(grplist(j-1,i)+1)) then
         if(dash.eq.0)then
          write(*,'(A1)',advance='no')'-'
          dash=1
         endif
        else
         if(dash.eq.1)then
          write(*,'(i4)',advance='no') grplist(j-1,i)
          dash=0
         endif
         write(*,'(i4)',advance='no') grplist(j,i)
        endif
       enddo 
       if(dash.eq.1)then
        write(*,'(i4)',advance='no') grplist(j-1,i)
        dash=0
       endif
      write(*,*)''
      enddo

c old display list code
c      write(*,*)'group #        atoms '
c      do i=1,ngrp      
c         write(*,'(i4,3x,100i3)')i,(grplist(j,i),j=1,grpn(i))
c      enddo

      eg=0
      iii=0
      do i=1,ngrp
         ni=grpn(i)
         iii=iii+1
         jjj=0
         do j=1,ngrp
            nj=grpn(j)
            jjj=jjj+1
            do ii=1,ni
               iiii=grplist(ii,i)
               do jj=1,nj
                  jjjj=grplist(jj,j)
                  if(jjjj.lt.iiii)cycle
                  eg(nwpwxc_lin(iii,jjj)) 
     >             = eg(nwpwxc_lin(iii,jjj))+ed(iiii,jjjj)
               enddo
            enddo
         enddo
      enddo

c     call prmat(6,eg,ngrp,0,'intra- + inter-group dispersion energies')
      write(*,*)' group i      j     Edisp'
      k=0
      check=0
      do i=1,ngrp
      do j=1,i    
      k=k+1
      check=check+eg(k) 
      write(*,'(5x,i4,'' --'',i4,F8.2)')i,j,eg(k)
      enddo
      enddo
      write(*,*) 'checksum (Edisp) ',check

      deallocate(dist,li)
      end subroutine nwpwxc_pbcadisp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE NWPWXC_SET_CRITERIA(rthr,lat,tau_max)

        REAL*8 :: r_cutoff,rthr
        REAL*8 :: lat(3,3)
        REAL*8 :: tau_max(3)
        REAL*8 :: norm1(3),norm2(3),norm3(3)
        REAL*8 :: cos10,cos21,cos32
        real*8,external :: nwpwxc_vectorsize

        r_cutoff=sqrt(rthr)
c          write(*,*) 'lat',lat
          !c find normal to the plane...
        call nwpwxc_kreuzprodukt(lat(:,2),lat(:,3),norm1)
        call nwpwxc_kreuzprodukt(lat(:,3),lat(:,1),norm2)
        call nwpwxc_kreuzprodukt(lat(:,1),lat(:,2),norm3)
c        write(*,*) 'norm2',norm2
        !c ...normalize it...
        norm1=norm1/NWPWXC_VECTORSIZE(norm1)
        norm2=norm2/NWPWXC_VECTORSIZE(norm2)
        norm3=norm3/NWPWXC_VECTORSIZE(norm3)
c        write(*,*) 'norm2_',norm2
          !c cos angles between normals and lattice vectors
        cos10=SUM(norm1*lat(:,1))
        cos21=SUM(norm2*lat(:,2))
        cos32=SUM(norm3*lat(:,3))
          !write(*,*) 'cos32',cos32
          !tau_max(1)=abs(2*r_cutoff/cos10)
          !tau_max(2)=abs(2*r_cutoff/cos21)
          !tau_max(3)=abs(2*r_cutoff/cos32)
          !write(*,*) 'r_cutoff',r_cutoff
        tau_max(1)=abs(r_cutoff/cos10)
        tau_max(2)=abs(r_cutoff/cos21)
        tau_max(3)=abs(r_cutoff/cos32)
c        write(*,'(3f8.4)')tau_max(1),tau_max(2),tau_max(3)
      END SUBROUTINE NWPWXC_SET_CRITERIA


      SUBROUTINE nwpwxc_kreuzprodukt(A,B,C)
        IMPLICIT NONE
  
        REAL*8 :: A(3),B(3)
        REAL*8 :: X,Y,Z
        REAL*8 :: C(3)
        
        X=A(2)*B(3)-B(2)*A(3)
        Y=A(3)*B(1)-B(3)*A(1)
        Z=A(1)*B(2)-B(1)*A(2)
        C=(/X,Y,Z/)
      END SUBROUTINE nwpwxc_kreuzprodukt

       FUNCTION NWPWXC_VECTORSIZE(VECT)

         REAL*8 :: VECT(3)
         REAL*8 :: SVECT(3)
         REAL*8 :: NWPWXC_VECTORSIZE

         SVECT=VECT*VECT
         NWPWXC_VECTORSIZE=SUM(SVECT)
         NWPWXC_VECTORSIZE=NWPWXC_VECTORSIZE**(0.5)
       END FUNCTION NWPWXC_VECTORSIZE

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine nwpwxc_inv_cell(x,a) !x is normal lat, a is lat^(-1)
      IMPLICIT NONE
      real*8, intent(in)   :: x(3,3) !unitcell vectors in direct space
      real*8, intent(out)  :: a(3,3) !unitcell vectors in reciprocal space
      integer i
      real*8 det
      
      a=0.0
      det=x(1,1)*x(2,2)*x(3,3)+x(1,2)*x(2,3)*x(3,1)+x(1,3)*x(2,1)*
     .    x(3,2)-x(1,3)*x(2,2)*x(3,1)-x(1,2)*x(2,1)*x(3,3)-x(1,1)*
     .    x(2,3)*x(3,2)
!      write(*,*)'Det:',det
      a(1,1)=x(2,2)*x(3,3)-x(2,3)*x(3,2)
      a(2,1)=x(2,3)*x(3,1)-x(2,1)*x(3,3)
      a(3,1)=x(2,1)*x(3,2)-x(2,2)*x(3,1)
      a(1,2)=x(1,3)*x(3,2)-x(1,2)*x(3,3)
      a(2,2)=x(1,1)*x(3,3)-x(1,3)*x(3,1)
      a(3,2)=x(1,2)*x(3,1)-x(1,1)*x(3,2)
      a(1,3)=x(1,2)*x(2,3)-x(1,3)*x(2,2)
      a(2,3)=x(1,3)*x(2,1)-x(1,1)*x(2,3)
      a(3,3)=x(1,1)*x(2,2)-x(1,2)*x(2,1)
      a=a/det
      end subroutine nwpwxc_inv_cell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE nwpwxc_xyz_to_abc(xyz,abc,lat,n)
      IMPLICIT NONE
      integer,intent(in) :: n
      real*8, INTENT(in) :: xyz(3,n)
      real*8, intent(in) :: lat(3,3)
      real*8, intent(out) :: abc(3,n)
      
      real*8 lat_1(3,3)
      integer i,j,k

      call nwpwxc_inv_cell(lat,lat_1)
      
      abc(:,:n)=0.0d0
      do i=1,n
        do j=1,3
          do k=1,3
            abc(j,i)=abc(j,i)+lat_1(j,k)*xyz(k,i)            
          enddo !k
          abc(j,i)=dmod(abc(j,i),1.0d0)
        enddo !j
      enddo !i
  
      END SUBROUTINE nwpwxc_xyz_to_abc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE nwpwxc_abc_to_xyz(abc,xyz,lat,n)
      IMPLICIT NONE
      real*8, INTENT(in) :: abc(3,*)
      real*8, intent(in) :: lat(3,3)
      real*8, intent(out) :: xyz(3,*)
      integer,intent(in) :: n
      
      integer i,j,k

      xyz(:,:n)=0.0d0
      do i=1,n
        do j=1,3
          do k=1,3
            xyz(j,i)=xyz(j,i)+lat(j,k)*abc(k,i)            
          enddo !k
        enddo !j
      enddo !i
  
      END SUBROUTINE nwpwxc_abc_to_xyz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      REAL*8 FUNCTION nwpwxc_volume(lat)
        IMPLICIT NONE
        REAL*8, INTENT(in) ::lat(3,3)
        REAL*8 zwerg
        
        zwerg=lat(1,1)*lat(2,2)*lat(3,3)+lat(1,2)*lat(2,3)*lat(3,1)+
     .        lat(1,3)*lat(2,1)*lat(3,2)-lat(1,3)*lat(2,2)*lat(3,1)-
     .        lat(1,2)*lat(2,1)*lat(3,3)-lat(1,1)*lat(2,3)*lat(3,2)
        nwpwxc_volume=abs(zwerg)
      END FUNCTION nwpwxc_volume







ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       
c            string pars procedures
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine nwpwxc_parse(str,delims,args,nargs)

! Parses the string 'str' into arguments args(1), ..., args(nargs) based on
! the delimiters contained in the string 'delims'. Preceding a delimiter in
! 'str' by a backslash (\) makes this particular instance not a delimiter.
! The integer output variable nargs contains the number of arguments found.
      interface 
        subroutine nwpwxc_split(str,delims,before,sep)
         character(len=*),intent(inout) :: str,before
         character(len=*),intent(in) :: delims
         character,optional,intent(inout) :: sep
       end subroutine nwpwxc_split
      end interface

      character(len=*),intent(inout) :: str
      character(len=*),intent(in) :: delims
      character(len=len_trim(str)) :: strsav
      character(len=*),dimension(:),intent(inout) :: args
      integer, intent(out) :: nargs
      
      strsav=str
      call nwpwxc_compact(str)
      na=size(args)
      do i=1,na
        args(i)=' '
      end do  
      nargs=0
      lenstr=len_trim(str)
      if(lenstr==0) return
      k=0

      do
         if(len_trim(str) == 0) exit
         nargs=nargs+1
         call nwpwxc_split(str,delims,args(nargs))
         call nwpwxc_removebksl(args(nargs))
      end do   
      str=strsav

      end subroutine nwpwxc_parse

!**********************************************************************

      subroutine nwpwxc_compact(str)

! Converts multiple spaces and tabs to single spaces; deletes control characters;
! removes initial spaces.

      character(len=*):: str
      character(len=1):: ch
      character(len=len_trim(str)):: outstr
      
      str=adjustl(str)
      lenstr=len_trim(str)
      outstr=' '
      isp=0
      k=0

      do i=1,lenstr
        ch=str(i:i)
        ich=iachar(ch)
  
        select case(ich)
  
          case(9,32)     ! space or tab character
            if(isp==0) then
              k=k+1
              outstr(k:k)=' '
            end if
            isp=1
            
          case(33:)      ! not a space, quote, or control character
            k=k+1
            outstr(k:k)=ch
            isp=0
      
        end select
        
      end do

      str=adjustl(outstr)

      end subroutine nwpwxc_compact

!**********************************************************************

      subroutine nwpwxc_removesp(str)

      ! Removes spaces, tabs, and control characters in string str

      character(len=*):: str
      character(len=1):: ch
      character(len=len_trim(str))::outstr

      str=adjustl(str)
      lenstr=len_trim(str)
      outstr=' '
      k=0

      do i=1,lenstr
        ch=str(i:i)
        ich=iachar(ch)
        select case(ich)    
          case(0:32)  ! space, tab, or control character
               cycle       
          case(33:)  
            k=k+1
            outstr(k:k)=ch
        end select
      end do
      
      str=adjustl(outstr)
      
      end subroutine nwpwxc_removesp



#if defined(XLFLINUX)
@PROCESS NOESCAPE
#endif
      subroutine nwpwxc_split(str,delims,before,sep)


! Routine finds the first instance of a character from 'delims' in the
! the string 'str'. The characters before the found delimiter are
! output in 'before'. The characters after the found delimiter are
! output in 'str'. The optional output character 'sep' contains the 
! found delimiter. A delimiter in 'str' is treated like an ordinary 
! character if it is preceded by a backslash (\). If the backslash 
! character is desired in 'str', then precede it with another backslash.

      character(len=*),intent(inout) :: str,before
      character(len=*),intent(in) :: delims
      character,optional :: sep
      logical :: pres
      character :: ch,cha

      pres=present(sep)
      str=adjustl(str)
      call nwpwxc_compact(str)
      lenstr=len_trim(str)
      if(lenstr == 0) return        ! string str is empty
      k=0
      ibsl=0                        ! backslash initially inactive
      before=' '
      do i=1,lenstr
         ch=str(i:i)
         if(ibsl == 1) then          ! backslash active
            k=k+1
            before(k:k)=ch
            ibsl=0
            cycle
         end if
         if(ch == '\') then          ! backslash with backslash inactive
            k=k+1
            before(k:k)=ch
            ibsl=1
            cycle
         end if
         ipos=index(delims,ch)         
         if(ipos == 0) then          ! character is not a delimiter
            k=k+1
            before(k:k)=ch
            cycle
         end if
         if(ch /= ' ') then          ! character is a delimiter that is not a space
            str=str(i+1:)
            if(pres) sep=ch
            exit
         end if
         cha=str(i+1:i+1)            ! character is a space delimiter
         iposa=index(delims,cha)
         if(iposa > 0) then          ! next character is a delimiter
            str=str(i+2:)
            if(pres) sep=cha
            exit
         else
            str=str(i+1:)
            if(pres) sep=ch
            exit
         end if
      end do
      if(i >= lenstr) str=''
      str=adjustl(str)              ! remove initial spaces
      return

      end subroutine nwpwxc_split

!**********************************************************************

#if defined(XLFLINUX)
@PROCESS NOESCAPE
#endif
      subroutine nwpwxc_removebksl(str)

! Removes backslash (\) characters. Double backslashes (\\) are replaced
! by a single backslash.

      character(len=*):: str
      character(len=1):: ch
      character(len=len_trim(str))::outstr

      str=adjustl(str)
      lenstr=len_trim(str)
      outstr=' '
      k=0
      ibsl=0                        ! backslash initially inactive
      
      do i=1,lenstr
        ch=str(i:i)
        if(ibsl == 1) then          ! backslash active
         k=k+1
         outstr(k:k)=ch
         ibsl=0
         cycle
        end if
        if(ch == '\') then          ! backslash with backslash inactive
         ibsl=1
         cycle
        end if
        k=k+1
        outstr(k:k)=ch              ! non-backslash with backslash inactive
      end do
      
      str=adjustl(outstr)
      
      end subroutine nwpwxc_removebksl



