c***********************************************************************
c
c
c   VERY IMPORTANT: Note that the same data structures for amatrix,
c                   amatrixinv, gmatrix, etc. are used for the 2
c                   and 3-dimensional cases. The last row and column
c                   of these matrices have been set equal to zero EXCEPT
c                   element (3,3) which is equal to 1. This allows a large
c                   amount of code reuse for the various cases. However,
c                   Beware of doing mindless matrix algebra with 
c                   these quantities since only the (1,1),(1,2),(2,1) and
c                   (2,2) elements have any physical meaning.
c
c
c         Direct space geometry stuff for the 2-d Bravias lattices
c
c 1) Computes the a-matrix used to transform from crystallographic to 
c    cartesian coordinates
c    note the following defn. of the amatrix to reduce confusion
c
c      A = a11 i + a21 j 
c      B = 0   i + a22 j 
c
c      thus: Amatrix =[ a11 0   ]
c                     [ a21 a22 ]
c
c     where: {A,B} are crystallographic axis, {i,j} unit Cartesian 
c            vectors. There is no concept of Crystallographic C-axis
c
c     So Amatrix premultiplies a fractional vector (a vector defined in
c     the crystal basis). Each successive row of a (when mult. times 
c     the vector) will give the particular cartesian component. 
c     Atranspose*A=metrix matrix. An alternative way's to think about 
c     its definition are:
c     
c**    Premultiplying a fractional vector by amatrix will result in a
c**    Cartesian vector in atomic units.
c
c    It represents a chosen convention of two axis systems.Where the 
c    oblique crystal system is oriented relative to an external 
c    rectilinear laboratory coordinate system such that: the b-axis is 
c    along j (unit rectilinear vector) ab-plane contains i (rect)
c    ,i.e., j along b, i in  ab-plane. This convention will result in
c    the z-direction of Cartesian space being perpendicular to the 
c    periodic repeat direction of the 2-d slab. 
c
c
c 2) computes the metric matrix for the lattice.
c
c 3) computes the direct space cell volume
c
c 4) computes reciprocal space lattice constants & volume 
c
c 5) computes the transformation matrix from reciprocal lattice vectors
c    to Cartesians (i.e, finds cartesian components of a recip lat vec)
c
c 6) computes some constants needed for Ewald routines that depend on
c    the geometry
c
c --> Important Conventions
c
c     The lattice angles are input in Degrees
c
c                                                   AC Hess
c***********************************************************************
      subroutine geom_2d(geom,c2au)
*
* $Id$
*
      implicit none
#include "inp.fh"
#include "geom.fh"
#include "nwc_const.fh"
#include "geomP.fh"
c
      integer geom,i,j
      double precision c2au
      double precision cdist(3),cang(3),gmat(3,3),amat(3,3),vol
      double precision c3,s3,rad,ainv(3,3)
      double precision pi,bmat(3,3)
c
      parameter(rad=57.295779513082343d0)
c
      pi=acos(0.0d0)*2.0d0
c
c--> get direct space lattice vectors
c
      do 100 i=1,3
        cdist(i)=lattice_vectors(i,geom)*c2au
  100 continue
        cang(3) =lattice_angles(3,geom)/rad
c
c        write(*,*) 'lattice vectors and angles'
c        write(*,*) 'vectors',(lattice_vectors(i,geom),i=1,3)
c        write(*,*) 'angles',(lattice_angles(i,geom),i=1,3)
c
c--------> build the metrical matrix (atomic units)
c     
      do 200 i=1,2
        gmat(i,i)=cdist(i)**2
  200 continue
      gmat(1,2)=cdist(1)*cdist(2)*cos(cang(3))
      gmat(2,1)=gmat(1,2)
      gmat(1,3)=0.0d+00
      gmat(2,3)=0.0d+00
      gmat(3,1)=0.0d+00
      gmat(3,2)=0.0d+00
      gmat(3,3)=0.0d+00
c      
      do 230 i=1,3
        do 240 j=1,3
          metric_matrix(i,j,geom)=gmat(i,j)
  240   continue
  230 continue
c
c      write(*,*) 'metrix matrix'
c      do i=1,3
c         write(*,*) (gmat(i,j),j=1,3)
c      enddo
c     
c-----------> build a-matrix
c     
c***j along b, i in  ab-plane
c
      s3=dsin(cang(3))
      c3=dcos(cang(3))
      amat(1,1) = cdist(1)*s3
      amat(1,2) = 0.0d+00
      amat(1,3) = 0.0d+00
      amat(2,1) = cdist(1)*c3
      amat(2,2) = cdist(2)
      amat(2,3) = 0.0d+00
      amat(3,1) = 0.0d+00
      amat(3,2) = 0.0d+00
      amat(3,3) = cdist(3)
c
c      write(*,*) 'a- matrix'
c      do i=1,3
c         write(*,*) (amat(i,j),j=1,3)
c      enddo
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
      vol = abs(amat(1,1)*amat(2,2)-amat(1,2)*amat(2,1))
c
      volume_direct(geom)=vol
c
c--> build amatrix_inv for 2-d systems
c
      ainv(1,1)= amat(2,2)/vol
      ainv(1,2)=-amat(1,2)/vol
      ainv(1,3)= 0.0d0
      ainv(2,1)=-amat(2,1)/vol
      ainv(2,2)= amat(1,1)/vol
      ainv(2,3)= 0.0d0
      ainv(3,1)= 0.0d0
      ainv(3,2)= 0.0d0
      ainv(3,3)= 1.0d0/cdist(3)
c
c      write(*,*) 'ainv- matrix'
c      do i=1,3
c         write(*,*) (ainv(i,j),j=1,3)
c      enddo
c
c     load it into common
c
      do 270 i=1,3
        do 280 j=1,3
          amatrix_inv(i,j,geom)=ainv(i,j)
  280   continue
  270 continue
c
c--> construct bmatrix (used to transform reciprocal vectors
c    to Cartesian form). 2pi*(At)-1
c
      bmat(1,1)=(2.0d0*pi/vol)*amat(2,2)
      bmat(1,2)=(2.0d0*pi/vol)*(-amat(2,1))
      bmat(1,3)=0.0d0
      bmat(2,1)=(2.0d0*pi/vol)*(-amat(1,2))
      bmat(2,2)=(2.0d0*pi/vol)*amat(1,1)
      bmat(2,3)=0.0d0
      bmat(3,1)=0.0d0
      bmat(3,2)=0.0d0
c      bmat(3,3)=1.0d0
      bmat(3,3)=2.0d0*pi/cdist(3)
c
c      write(*,*) 'b-matrix'
c      do i=1,3
c         write(*,*) (bmat(i,j),j=1,3)
c      enddo
c
c load in into common
c
      do 290 i=1,3
         do 300 j=1,3
            bmatrix(i,j,geom)=bmat(i,j)
 300     continue
 290  continue
c
      return
      end

      
      subroutine geom_2d_amatrix(geom,scale)
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








