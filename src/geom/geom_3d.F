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
c     The lattice angles are input in Degrees
c***********************************************************************
      subroutine geom_3d(geom,c2au)
C$Id$
      implicit none
#include "global.fh"
#include "inp.fh"
#include "geom.fh"
#include "stdio.fh"
#include "nwc_const.fh"
#include "geomP.fh"
c
      double precision deter3
      external deter3
      integer geom,i,j,iang
      integer me
*      integer indx(3)
      double precision c2au
      double precision cdist(3),cang(3),gmat(3,3),amat(3,3),dperm,vol
      double precision c1,c2,c3,s1,s2,s3,rad,ainv(3,3),rdist(3),rang(3)
      double precision bmat(3,3),scale,pi
c
      parameter(rad=57.295779513082343d0)
c
      pi=acos(0.0d0)*2.0d0
c
      me  = ga_nodeid()
c
c--> get direct space lattice vectors
c
      do 100 i=1,3
        cdist(i)=lattice_vectors(i,geom)*c2au
        cang(i) =lattice_angles(i,geom)/rad
  100 continue
c
c--> convert angles to radians
c
c      do 110 i=1,3
c        cang(i)=cang(i)/rad
c  110 continue
c     
c--------> build the metrical matrix (atomic units)
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
c----------> compute the volume (atomic units)
c     
*      call ludcmp(gmat,3,3,indx,dperm)
*      do 245 i=1,3
*        dperm=dperm*gmat(i,i)
*  245 continue
      if (me.eq.0) then 
        write(LuOut,*)
        write(LuOut,*) ' !!!!!!!!!  geom_3d NEEDS TESTING !!!!!!!!!! '
      endif
      dperm = deter3(gmat)
*
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
c  from angstroms to bohr; i.e., amat is in atomic units.
c
c  the Jth column of this matrix contains the jth lattice vector
c aij is the ith Cartesian component of the jth lattice vector.
c
c
      c1=dcos(cang(1))
      c2=dcos(cang(2))
      c3=dcos(cang(3))
      s3=dsin(cang(3))
      amat(1,1) = cdist(1)*s3
      amat(1,2) = 0.0d+00
      amat(1,3) = (cdist(3)*(c2-c1*c3)/s3)
      amat(2,1) = cdist(1)*c3
      amat(2,2) = cdist(2)
      amat(2,3) = cdist(3)*c1
      amat(3,1) = 0.0d+00
      amat(3,2) = 0.0d+00
      amat(3,3) = (vol/(cdist(1)*cdist(2)*s3))
c     
c     load it into common
c
      do 250 i=1,3
        do 260 j=1,3
          amatrix(i,j,geom)=amat(i,j)
  260   continue
  250 continue
c
c compute direct space volume in atomic bohr**3
c
c      vol = abs(amat(1,1)*(amat(2,2)*amat(3,3)-amat(3,2)*amat(2,3)) +
c     $     amat(1,2)*(amat(2,3)*amat(3,1)-amat(3,3)*amat(2,1)) +
c     $     amat(1,3)*(amat(2,1)*amat(3,2)-amat(3,1)*amat(2,2)))
c      volume_direct(geom)=vol
c
c--> build amatrix_inv for 3-d systems (units=inverse a.u.)
c
      c1=dcos(cang(1))
      c2=dcos(cang(2))
      c3=dcos(cang(3))
      s3=dsin(cang(3))
c
      ainv(1,1) = (1.0d+00/(cdist(1)*s3))
      ainv(1,2) = 0.0d+00
      ainv(1,3) = (cdist(2)*cdist(3)*(c1*c3-c2)/(vol*s3))
      ainv(2,1) = (-c3/(cdist(2)*s3))
      ainv(2,2) = (1.0d+00/cdist(2))
      ainv(2,3) = (cdist(1)*cdist(3)*(c2*c3-c1)/(vol*s3))
      ainv(3,1) = 0.0d+00
      ainv(3,2) = 0.0d+00
      ainv(3,3) = (cdist(1)*cdist(2)*s3/vol)
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
c  Cartesian components of reciprocal lattice vectors in a.u.
c
      scale=(2.0d0*pi/vol)
c
      bmat(1,1)=scale*(amat(2,2)*amat(3,3)-amat(3,2)*amat(2,3))
      bmat(1,2)=scale*(amat(2,3)*amat(3,1)-amat(3,3)*amat(2,1))
      bmat(1,3)=scale*(amat(2,1)*amat(3,2)-amat(3,1)*amat(2,2))
      bmat(2,1)=scale*(amat(3,2)*amat(1,3)-amat(1,2)*amat(3,3))
      bmat(2,2)=scale*(amat(3,3)*amat(1,1)-amat(1,3)*amat(3,1))
      bmat(2,3)=scale*(amat(3,1)*amat(1,2)-amat(1,1)*amat(3,2))
      bmat(3,1)=scale*(amat(1,2)*amat(2,3)-amat(2,2)*amat(1,3))
      bmat(3,2)=scale*(amat(1,3)*amat(2,1)-amat(2,3)*amat(1,1))
      bmat(3,3)=scale*(amat(1,1)*amat(2,2)-amat(2,1)*amat(1,2))
c
c  load bmatrix into common
c
      do 290 i=1,3
         do 300 j=1,3
            bmatrix(i,j,geom)=bmat(i,j)
 300     continue
 290  continue
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

      
      subroutine geom_3d_amatrix(geom,scale)
      implicit none
      integer geom
      double precision scale

#include "errquit.fh"
#include "nwc_const.fh"
#include "geomP.fh"

*     !**** local variables ****
      integer i,j    
      double precision amat(3,3)   
      double precision c(3,3), vol
c
c
      do i=1,3
        do j=1,3
          amat(i,j) = amatrix(i,j,geom)
        end do
      end do
c
c     Mmmm ... the original code only set this stuff from the input
c     using the a,b,c,alpha,beta,gamma, but now we have changed
c     the amatrix ... need to update ainv and also recompute the
c     other crap ... for now just set the other crap to crap so that
c     we'll know if it is used
c
      do i = 1,3
         do j = 1,3
            metric_matrix(i,j,geom) = 1d300
            bmatrix(i,j,geom) = 1d300
         end do
         recip_lat_vectors(i,geom) = 1d300
         recip_lat_angles(i,geom) = 1d300
      end do
c
c     HERE SHOULD RECOMPUTE AMATRIX WITH STANDARD ORIENTATION
c     SINCE IF THE GEOMETRY IS STORED AND RELOADED THE
c     STANDARD ORIENTATION IS IMPOSED.
c
c     Update the amatrix inverse
c      - Since amat=[a1,a2,a3]
c              ainv=[b1,b2,b3]^t
c
      call dfill(9,0.0d0,c,1)
      c(1,1) = amat(2,2)*amat(3,3) - amat(3,2)*amat(2,3)  ! = b(1,1)
      c(1,2) = amat(3,2)*amat(1,3) - amat(1,2)*amat(3,3)  ! = b(2,1)
      c(1,3) = amat(1,2)*amat(2,3) - amat(2,2)*amat(1,3)  ! = b(3,1)
      c(2,1) = amat(2,3)*amat(3,1) - amat(3,3)*amat(2,1)  ! = b(1,2)
      c(2,2) = amat(3,3)*amat(1,1) - amat(1,3)*amat(3,1)  ! = b(2,2)
      c(2,3) = amat(1,3)*amat(2,1) - amat(2,3)*amat(1,1)  ! = b(3,2)
      c(3,1) = amat(2,1)*amat(3,2) - amat(3,1)*amat(2,2)  ! = b(1,3)
      c(3,2) = amat(3,1)*amat(1,2) - amat(1,1)*amat(3,2)  ! = b(2,3)
      c(3,3) = amat(1,1)*amat(2,2) - amat(2,1)*amat(1,2)  ! = b(3,3)
      vol = amat(1,1)*c(1,1)
     >    + amat(2,1)*c(1,2)
     >    + amat(3,1)*c(1,3)
      volume_direct(geom) = vol
c
      call dscal(9,1.0d0/vol,c,1)
c
      call dcopy(9,c,1,amatrix_inv(1,1,geom),1)
c
c     Ooops ... must also update the pesky lattice parameters
c
      call xlattice_abc_abg(
     $     lattice_vectors(1,geom),
     $     lattice_vectors(2,geom),
     $     lattice_vectors(3,geom),
     $     lattice_angles(1,geom),
     $     lattice_angles(2,geom),
     $     lattice_angles(3,geom),amat)

      lattice_vectors(1,geom) = lattice_vectors(1,geom)/scale
      lattice_vectors(2,geom) = lattice_vectors(2,geom)/scale
      lattice_vectors(3,geom) = lattice_vectors(3,geom)/scale
c
      return
      end

