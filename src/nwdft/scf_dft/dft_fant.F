c
      double precision function dft_fant(geom,natoms,a,d)
c
C     Bonacic-Kouteck, V.; Cespiva, L.; Fantucci, P.; Pittner, J.; Kouteck, J. J. 
C     Chem. Phys. 1994, 98, 490.
C     Mitric, R.; Hartmann, M.; Stanca, B.; Bonacic-Koutecky, V.; Fantucci, P.;
C     J. Phys. Chem. A.; (Article); 2001; 105(39); 8892-8905
C     The core-core repulsion has been corrected according to
C     (CC(rij) ) l/rij + D exp(-arij)). 
C     Constants D and a obtained for 1e-RECP's from the fitting procedure
C     described in Appendix A necessery for corection of core-core potential: 
C     DAg-Ag ) 1619.887392, aAg-Ag ) 2.49630105, 
C     DAu-Au ) 1911.131090, aAu-Au ) 2.46590129, 
C     DAu-Ag ) 1765.509532, aAu-Ag ) 2.48110117. 
c
      implicit none
c
#include "errquit.fh"
#include "geom.fh"
#include "global.fh"
#include "stdio.fh"
c
      integer geom,natoms
      double precision a,d
c
      integer i,j
      double precision qi,qj
      double precision xi(3),xj(3),rij
      character*16 tag
c     
      dft_fant=0d0
      do i=1,natoms
         if (.not. geom_cent_get(geom, i, tag, 
     &        xi, qi))call errquit
     &        ('dft_fant: geom_cent_get failed', 0, GEOM_ERR)
         do j=1,i-1
            if (.not. geom_cent_get(geom, j, tag, 
     &           xj, qj))call errquit
     &           ('dft_fant: geom_cent_get failed', 0, GEOM_ERR)
            rij=sqrt((xi(1)-xj(1))**2+(xi(2)-xj(2))**2+
     +           (xi(3)-xj(3))**2)
            dft_fant=dft_fant + d*exp(-a*rij)

         enddo
      enddo
c
      if(ga_nodeid().eq.0)
     W     write(luout,*) 
     W     ' Bonacic-Koutecky Fantucci-Repulsive Term ',dft_fant
c
      return
      end
c $Id$
