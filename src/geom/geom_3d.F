c***********************************************************************
c
c         Direct space geometry stuff for the 3-d Bravias lattices
c
c 1) computes the a-matrix used to transform from crystallographic to 
c    cartesian coordinates
c
c 2) computes the metric matrix for the lattice
c
c 3) computes the direct space cell volume
c
c 4) computes reciprocal space lattice constants & volume 
c
c --> Important Conventions
c
c     The lattice constants are input and stored in Angstroms and Degrees
c***********************************************************************
      subroutine geom_3d(geom)
C$Id: geom_3d.F,v 1.3 1995-02-02 23:12:44 d3g681 Exp $
      implicit none
#include "inp.fh"
#include "geom.fh"
#include "geomP.fh"

      integer geom,i,j,iang,indx(3)
      real*8 cdist(3),cang(3),gmat(3,3),amat(3,3),dperm,vol,
     &       c1,c2,c3,s1,s2,s3,rad,ainv(3,3),rdist(3),
     &       rang(3)
      parameter(rad=57.295779513082343d0)
c
c--> get direct space lattice vectors
c
      do 100 i=1,3
        cdist(i)=lattice_vectors(i,geom)
        cang(i) =lattice_angles(i,geom)
  100 continue
c
c--> convert angles to radians
c
      do 110 i=1,3
        cang(i)=cang(i)/rad
  110 continue
c     
c--------> build the metrical matrix
c     
      
      do 200 i=1,3
        gmat(i,i)=cdist(i)**2
  200 continue
      iang=3
      do 210 i=1,3
        do 220 j=i+1,3
          gmat(i,j)=cdist(i)*cdist(j)*dcos(cang(iang))
          gmat(j,i)=gmat(i,j)
          iang=iang-1
  220   continue
  210 continue
c      
      do 230 i=1,3
        do 240 j=1,3
          metric_matrix(i,j,geom)=gmat(i,j)
  240   continue
  230 continue
c     
c----------> compute the volume
c     
      call ludcmp(gmat,3,3,indx,dperm)
      do 245 i=1,3
        dperm=dperm*gmat(i,i)
  245 continue
      vol=dsqrt(dperm)             
      volume_direct(geom)=vol
c     
c-----------> build a-matrix
c     
c Chosen convention of two axis systems.
c    Orient oblique crystal system to external rectilinear laboratory
c    coordinate system such that
c    b-axis is along j (unit rectilinear vector), ab-plane contains i (rect),
c    i.e., j along b, i in  ab-plane
c     
c  Note this matrix is implicitly built to convert the input coordinates
c  from angstroms to bohr.
c
      c1=dcos(cang(1))
      c2=dcos(cang(2))
      c3=dcos(cang(3))
      s3=dsin(cang(3))
      amat(1,1) = cdist(1)*s3*angstrom_to_au
      amat(1,2) = 0.0d+00
      amat(1,3) = (cdist(3)*(c2-c1*c3)/s3)*angstrom_to_au
      amat(2,1) = cdist(1)*c3*angstrom_to_au
      amat(2,2) = cdist(2)*angstrom_to_au
      amat(2,3) = cdist(3)*c1*angstrom_to_au
      amat(3,1) = 0.0d+00
      amat(3,2) = 0.0d+00
      amat(3,3) = (vol/(cdist(1)*cdist(2)*s3))*angstrom_to_au
c     
c     load it into common
c
      do 250 i=1,3
        do 260 j=1,3
          amatrix(i,j,geom)=amat(i,j)
  260   continue
  250 continue
c
c--> build amatrix_inv for 3-d systems
c
      c1=dcos(cang(1))
      c2=dcos(cang(2))
      c3=dcos(cang(3))
      s3=dsin(cang(3))

      ainv(1,1) = (1.0d+00/(cdist(1)*s3))/angstrom_to_au
      ainv(1,2) = 0.0d+00
      ainv(1,3) = (cdist(2)*cdist(3)*(c1*c3-c2)/(vol*s3))/angstrom_to_au
      ainv(2,1) = (-c3/(cdist(2)*s3))/angstrom_to_au
      ainv(2,2) = (1.0d+00/cdist(2))/angstrom_to_au
      ainv(2,3) = (cdist(1)*cdist(3)*(c2*c3-c1)/(vol*s3))/angstrom_to_au
      ainv(3,1) = 0.0d+00
      ainv(3,2) = 0.0d+00
      ainv(3,3) = (cdist(1)*cdist(2)*s3/vol)/angstrom_to_au
c     
c     load it into common
c
      do 270 i=1,3
        do 280 j=1,3
          amatrix_inv(i,j,geom)=ainv(i,j)
  280   continue
  270 continue
c
c--> construct reciprocal space vectors and angles
c
      c1=dcos(cang(1))
      c2=dcos(cang(2))
      c3=dcos(cang(3))
c
      s1=dsin(cang(1))
      s2=dsin(cang(2))
      s3=dsin(cang(3))
c
      rdist(1)=((cdist(2)*cdist(3))/vol)*s1
      rdist(2)=((cdist(1)*cdist(3))/vol)*s2
      rdist(3)=((cdist(1)*cdist(2))/vol)*s3
c
      rang(1)=dacos((c2*c3-c1)/(s2*s3))
      rang(2)=dacos((c1*c3-c2)/(s1*s3))
      rang(3)=dacos((c1*c2-c3)/(s1*s2))
c
c--> load
c
      do 2000 i=1,3
        recip_lat_vectors(i,geom)=rdist(i)
        recip_lat_angles(i,geom)=rang(i)*rad
 2000 continue

c
c***angles in radians
c
c      write(*,10)
c      write(*,20) (rdist(i), i=1,3)
c      write(*,40)
c      write(*,30) (rang(i), i=1,3)
c
c***convert to degrees (angles in dang/cang in radians)
c
c      rang(1)=rang(1)*rad
c      rang(2)=rang(2)*rad
c      rang(3)=rang(3)*rad
c      write(*,50)
c      write(*,30) (rang(i), i=1,3)
c10    format(/,"reciprocal Lattice Paramaters")
c20    format("distances:a=",f8.4,3x,"b= ",f8.4,3x,"c= ",f8.4)
c30    format("angles:  alpha=",f8.4,3x,"beta= ",f8.4,3x,"gamma= "
c     &,f8.4)
c40    format(/,"angles in radians")
c50    format(/,"angles in degrees")
c
      return
      end

      
