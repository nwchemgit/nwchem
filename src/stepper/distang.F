C> \ingroup stpr_priv
C> @{
      subroutine distan(nat,coord)
C     $Id$
      implicit none
      integer nat,i,j
      double precision dist
      external dist
      double precision coord(3,*),rij
      double precision au2ang
      write(6,1000)
      au2ang = 0.529177249d00
      do 20 i=1,nat-1 
         write(6,1100)
         do 10 j=i+1,nat 
            rij=dsqrt(dist(i,j,coord))       
            if (rij.le.5.0d00)
     $           write(6,1200) ' ', i, ' ', j, rij, rij*au2ang
 10      continue
 20   continue
      write(6,1100)
      return        
 1000 format(//,6x,'internuclear distances      (atomic units)',    
     .     '     angstroms')        
 1100 format(6x, 58(1h-))
 1200 format('        ',a6,i2,' <---> ',a6,i2,'   (', f12.6,')',f14.6,
     .     '  ') 
      end
      subroutine angle(nat,coord)
      integer nat,i,j,k
      double precision coord(3,nat)
      double precision dij,rij,dik,djk,rjk
      double precision acosarg,thet
      double precision dist
      external dist
c 
      write(6,1000)
      do i=1,nat-1 
         write(6,1100)
         do j=i+1,nat 
            dij=dist(i,j,coord) 
            rij=sqrt(dij) 
            if (rij .lt. 5.0d0) then
               do k=j+1,nat
                  dik=dist(i,k,coord)      
                  djk=dist(j,k,coord)
                  if (max(dij,dik,djk).le.25.0d0) then
                     rjk =sqrt(djk) 
                     acosarg =  (dij+djk-dik)/(2*rij*rjk)
                     acosarg = int(acosarg*1.0d+06)
                     acosarg = acosarg*1.0d-6
                     thet=(dacos(acosarg))*57.29578d0
                     write(6,1200) ' ', i, ' ', j, ' ', k, thet 
                  endif
               enddo
            endif
         enddo
      enddo
      return        
 1000 format(//,'      internuclear angles',30x,'degrees') 
 1100 format(6x, 58(1h-))
 1200 format('        ',a6,i2,2(' <---> ',a6,i2),f16.6,'  ')        
      end
      double precision function  dist(at1,at2,coord)
      implicit none
      integer at1,at2
      double precision coord(3,*)
c.....................................................................
c     w a r n i n g. this routine computes the square of the distance.
c.....................................................................
      dist = (coord(1,at1)-coord(1,at2))**2 + 
     +     (coord(2,at1)-coord(2,at2))**2 + 
     +     (coord(3,at1)-coord(3,at2))**2 
      return        
      end
C> @}
