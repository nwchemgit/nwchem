c***********************************************************************
*
* $Id$
*
c
c                            dctr.f
c
c   This program converts the space group symmetry operators of the centered
c    lattices from the non-primitive basis to the primitive basis.
c
c              Written by Daryl G. Clerc   May 22, 1995
c
c***********************************************************************
c
c   The conventions used to relate the centered cells (a'b'c') to the 
c  primitive cells (abc) are as follows:
c
c    A-centered: The author's (DGC)
c                a'=a, b'=(b+c)/2, c'=(-b+c)/2
c    C-centered: Orthohexagonal cell convention (from Int. Tables, Fig. 5.8, p.70, cell #C2)
c                a'=(a+b)/2, b'=(-a+b)/2, c'=c
c    R-centered: Obverse Int. Tables, Table 5.1, p.78 - Cell #1
c                See also Figs. 5.7a & c, p.79)
c                   a' = 2a/3 + b/3 + c/3
c                   b' = -a/3 + b/3 + c/3
c                   c' = -a/3 -2b/3 + c/3
c    I-centered: Int. Tables, Table 5.1, p.76. See also Fig. 5.4, p.77.
c                   a' = -a/2 + b/2 + c/2
c                   b' =  a/2 - b/2 + c/2
c                   c' =  a/2 + b/2 - c/2
c    F-centered: Int. Tables, Table 5.1, p.77. See also Fig. 5.5, p.77.
c                a' = (b+c)/2, b' = (a+c)/2, c' = (a+b)/2
c***********************************************************************
c
c  Each space group operator is defined as {R|v - Rs + s}, where v = the translation
c  associated with the point group operator R (e.g. the translation part of a glide
c  or screw operation) and s = shift of origin.
c  
c  The operators in the de-centered coordinate systems (primed) are related
c  to the original (centered, non-primed) operators as follows:
c  
c        { R' | v' - R's' + s' } = { P^-1RP | P^-1(v - Rs + s) },
c
c  where P is a 3x3 matrix which transforms coordinates as Pr'=r,
c        r  = vector r components in the non-primitive (centered) basis,
c        r' = vector r components in primitive (de-centered) basis.
c
c   Note that the de-centered crystal possesses the same symmetry as the
c    centered crystal except for the centering vectors. The centering vectors
c    form at least two of the lattice translations in the de-centered system.
c   Because the orientation of the de-centered lattice vectors relative to the
c    symmetry elements is different from that of the centered lattice vectors,
c    the symmetry elements must be transformed as described above.
c
c***********************************************************************
c
c---> The following is a list of the centered 3D groups:
c
c      A-centered (4 total): 38, 39, 40, 41
c
c      C-centered (16 total): 5, 8, 9, 12, 15, 20, 21, 35, 36, 37, 63,
c                              64, 65, 66, 67, 68
c
c      R-centered (7 total):146,148,155,160,161,166,167
c
c      I-centered (38 total):
c          23, 24, 44, 45, 46, 71, 72, 73, 74, 79, 80, 82, 87, 88, 97, 98,
c           107,108,109,110,119,120,121,122,139,140,141,142,197,199,204,206,
c           211,214,217,220,229,230
c
c      F-centered (16 total)
c            22, 42, 43, 69, 70,196,202,203,209,210,216,219,225,226,227,228
c  
c*******************************************************************************************
c
c---> The following is a list of the centered 2D groups:
c      C-centered (9 total):10, 13, 16, 22, 26, 35, 36, 47, 48
c
c*******************************************************************************************
c
      subroutine dctr( symops, nops_dctr, numgrp, sol_name, numvec,
     & cntvec, numset,geom)
c
      integer i,j,maxops,ipointer,numgrp,nops_dctr,numvec,numset,geom
      double precision tf1(90),tf(6,3),temp(3,4),temp2(3,4),trans(3)
      character *10 sol_name
      character *1 ctr_type
c
      parameter(maxops=192)
      double precision symops(maxops*3,4),cntvec(3,3),xcnttot,one,tol
      double precision lattice(6)

      logical  value
      logical  geom_convert_to_primitive
      external geom_convert_to_primitive
c
c.......................................................................
c--->Assign transformation matrices to 1-dim array in successive rows as follows:
c    P1^-1,P1; P2^-1,P2;...,P5^-1,P5.  Pr'=r and the order is A-C-R-I-F.
c.......................................................................
      data (tf1(i),i=1,90)/
     & 1.0000000000000000d0, 0.0000000000000000d0, 0.0000000000000000d0,
     & 0.0000000000000000d0, 1.0000000000000000d0, 1.0000000000000000d0,
     & 0.0000000000000000d0,-1.0000000000000000d0, 1.0000000000000000d0,
     & 1.0000000000000000d0, 0.0000000000000000d0, 0.0000000000000000d0,
     & 0.0000000000000000d0, 0.5000000000000000d0,-0.5000000000000000d0,
     & 0.0000000000000000d0, 0.5000000000000000d0, 0.5000000000000000d0,
     & 1.0000000000000000d0, 1.0000000000000000d0, 0.0000000000000000d0,
     &-1.0000000000000000d0, 1.0000000000000000d0, 0.0000000000000000d0,
     & 0.0000000000000000d0, 0.0000000000000000d0, 1.0000000000000000d0,
     & 0.5000000000000000d0,-0.5000000000000000d0, 0.0000000000000000d0,
     & 0.5000000000000000d0, 0.5000000000000000d0, 0.0000000000000000d0,
     & 0.0000000000000000d0, 0.0000000000000000d0, 1.0000000000000000d0,
     & 1.0000000000000000d0, 0.0000000000000000d0, 1.0000000000000000d0,
     &-1.0000000000000000d0, 1.0000000000000000d0, 1.0000000000000000d0,
     & 0.0000000000000000d0,-1.0000000000000000d0, 1.0000000000000000d0,
     & 0.6666666666666667d0,-0.3333333333333333d0,-0.3333333333333333d0,
     & 0.3333333333333333d0, 0.3333333333333333d0,-0.6666666666666667d0,
     & 0.3333333333333333d0, 0.3333333333333333d0, 0.3333333333333333d0,
     & 0.0000000000000000d0, 1.0000000000000000d0, 1.0000000000000000d0,
     & 1.0000000000000000d0, 0.0000000000000000d0, 1.0000000000000000d0,
     & 1.0000000000000000d0, 1.0000000000000000d0, 0.0000000000000000d0,
     &-0.5000000000000000d0, 0.5000000000000000d0, 0.5000000000000000d0,
     & 0.5000000000000000d0,-0.5000000000000000d0, 0.5000000000000000d0,
     & 0.5000000000000000d0, 0.5000000000000000d0,-0.5000000000000000d0,
     &-1.0000000000000000d0, 1.0000000000000000d0, 1.0000000000000000d0,
     & 1.0000000000000000d0,-1.0000000000000000d0, 1.0000000000000000d0,
     & 1.0000000000000000d0, 1.0000000000000000d0,-1.0000000000000000d0,
     & 0.0000000000000000d0, 0.5000000000000000d0, 0.5000000000000000d0,
     & 0.5000000000000000d0, 0.0000000000000000d0, 0.5000000000000000d0,
     & 0.5000000000000000d0, 0.5000000000000000d0, 0.0000000000000000d0/
c
c.......................................................................
c       Assign tf(i,j) row pointer based on the type of centering
c.......................................................................
c
      one = 1.00d+00
      tol = 1.00d-12
c
      ctr_type='Z'
      xcnttot=cntvec(1,1)+cntvec(1,2)+cntvec(1,3)
c
      if(xcnttot.eq.1.0.and.cntvec(1,1).eq.0.0.and.numvec.eq.1) then
         ctr_type='A'
         ipointer=0
      endif
      if(xcnttot.eq.1.0.and.cntvec(1,3).eq.0.0.and.numvec.eq.1) then
         ctr_type='C'
         ipointer=18
      endif
      if(xcnttot.gt.1.33.and.xcnttot.lt.1.34) then
         ctr_type='R'
         ipointer=36
      endif
      if(xcnttot.eq.1.5) then
         ctr_type='I'
         ipointer=54
      endif
      if(numvec.eq.3) then
         ctr_type='F'
         ipointer=72
      endif
      if(ctr_type.eq.'Z') then
         write(*,985)
         stop
      endif
985   format(5x,'ERROR IN SYMMETRY/DCTR.F ... INAPPROPRIATE CENTERING')
c
c.......................................................................
c--->Put 1-dim array tf1() data into two 3-dim tf() arrays
c.......................................................................
c
      nrows=6
      do 100 i=1,nrows
         do 110 j=1,3
            tf(i,j)=tf1(ipointer+3*(i-1)+j)
110      continue
c         write(*,1000) (tf(i,j),j=1,3)
100   continue
      write(*,*) ' '
c
1000  format(5x,3(f7.4,3x))
c
c.......................................................................
c          Print transformation matrices & number of symops
c.......................................................................
      write(*,989)
989   format('**********************************************************
     &**************************************')
      write(*,*) ' '
      write(*,970)
970   format('DECENTERING   ...')
      write(*,*) ' '
      write(*,987) ctr_type,sol_name,numgrp,numset
987   format(3x,'The lattice is ',a1,'-centered:  Group name = ',a10,', 
     & Group No. = ',i3,',  Setting = ',i1)
      write(*,*) ' '
      write(*,990) numgrp,nops_dctr+1
990   format('The number of symmetry operations (incl. E) in group no. '
     &,i3,' after decentering is ',i3)
c
      write(*,994)
994   format(/,22x,'The basis transformation Prp = r is:',/)
      write(*,997) (tf(4,j),j=1,3)
      write(*,996) (tf(5,j),j=1,3)
      write(*,995) (tf(6,j),j=1,3)
997   format(20x,'x  = ',f7.4,' xp  +  ',f7.4,' yp  +  ',f7.4,' zp')
996   format(20x,'y  = ',f7.4,' xp  +  ',f7.4,' yp  +  ',f7.4,' zp')
995   format(20x,'z  = ',f7.4,' xp  +  ',f7.4,' yp  +  ',f7.4,' zp')
c
      write(*,998)
998   format(/,22x,'The inverse transformation P^-1r = rp is:',/)
      write(*,993) (tf(1,j),j=1,3)
      write(*,992) (tf(2,j),j=1,3)
      write(*,991) (tf(3,j),j=1,3)
993   format(20x,'xp = ',f7.4,' x   +  ',f7.4,' y   +  ',f7.4,' z')
992   format(20x,'yp = ',f7.4,' x   +  ',f7.4,' y   +  ',f7.4,' z')
991   format(20x,'zp = ',f7.4,' x   +  ',f7.4,' y   +  ',f7.4,' z')
      write(*,*) ' '
      write(*,989)
      write(*,*) ' '
      write(*,*) ' '
c
c.......................................................................
c                          Begin Transformations
c.......................................................................
      do 190 iop=1,nops_dctr
c         write(*,999)
999   format('**********************************************************
     &******************')
c
c--->Zero out temp arrays
      do 10 i=1,3
         do 20 j=1,4
            temp(i,j)=0d0
            temp2(i,j)=0d0
            trans(j)=0d0
20       continue
10    continue
c
c--->Begin similarity transformation: M=P^-1RP ... compute RP
      do 200 i=1,3
         do 220 j=1,3
            do 240 k=1,3
               temp(i,j)=temp(i,j)+symops(3*(iop-1)+i,k)*tf(k+3,j)
240         continue
220      continue
c         write(*,1004) 3*(iop-1)+i,(symops(3*(iop-1)+i,k),k=1,4)
1004     format('symops(',i3,',1-4)=',3x,4(f7.4,3x))
200   continue
c      write(*,*) ' '
c
c--->Compute P^-1*RP
      do 300 i=1,3
         do 320 j=1,3
            do 340 k=1,3
               temp2(i,j)=temp2(i,j)+tf(i,k)*temp(k,j)
340         continue
320      continue
300   continue
c      write(*,*) ' '
c
c--->Compute P^-1(v - Rs + s) (i.e. transform translational part)
      do 400 i=1,3
         do 420 j=1,3
            trans(i)=trans(i)+tf(i,j)*symops(3*(iop-1)+j,4)
420      continue
         temp2(i,4)=trans(i)
c         write(*,1001) i,(temp2(i,j),j=1,4)
1001     format('temp2(',i3,',1-4)=',3x,4(f7.4,3x))
400   continue
c      write(*,999)
c
c--->Remove integral translations
      do 424 i=1,3
         if(temp2(i,4).ge.one-tol.and.temp2(i,4).le.one+tol) then
            temp2(i,4) = temp2(i,4) - one
         endif
424   continue
c.......................................................................
c                  Set symops=temp2
c.......................................................................
      do 500 i=3*(iop-1)+1,3*iop
         do 510 j=1,4
            symops(i,j)=temp2(i-3*(iop-1),j)
510      continue
500   continue
c
c--->Transform next operator in symops list
190   continue
c.......................................................................
c                      End Transformations
c.......................................................................
c
      value =  geom_convert_to_primitive(geom,ctr_type,tf)

      return
      end
