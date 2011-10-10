      subroutine smd_pot_compute_generic(np,
     >                           nl,
     >                           na,
     >                           rcutsq,
     >                           it,
     >                           p1,
     >                           p2,
     >                           p3,
     >                           pt,
     >                           ccc,
     >                           fff,
     >                           point,
     >                           list,
     >                           e)

      implicit none

#include "smd_const_data.fh"


c       np         [in] number of potential parameters
c       na         [in] number of atoms
c       nl         [in] size of verlet list
c       rcutsq     [in] cutoff radius squared
c       it(na)     [in] type index
c       p1(np)     [in] potential parameter 1
c       p2(np)     [in] potential parameter 2
c       p3(np)     [in] potential parameter 3
c       pt(np)     [in] potential type
c       point(na)  [in] verlet pointer
c       list(nl)   [in] verlet list
c       ccc(nl,3)  [in] verlet pair distance vector
c       fff(na,3)  [inout] force
c       e          [out]   energy

      integer np
      integer na
      integer nl
      double precision rcutsq
      integer it(na)
      double precision p1(np)
      double precision p2(np)
      double precision p3(np)
      integer pt(np)
      integer point(na)
      integer list(nl)
      double precision ccc(3,nl)
      double precision fff(3,na)
      double precision e,etmp
c
      integer i,j,k,jnab
      integer jbeg,jend
      integer nlist
      integer it1,it2,nt
      integer potindex

      double precision  dr,rsq
      double precision  erfxc,force

      double precision x,y,z

      e=0

c     figure out number of types
c     from nt(nt+1)/2=np
c     ------------------
      nt=int(sqrt(1+8.0*np))
      nt=(-1+nt)/2
      nlist = 0
      do i=1,na-1

       it1 = it(i)
       jbeg=point(i)
       jend=point(i+1)-1


      do jnab=jbeg,jend

       etmp = 0.0d0
       j=list(jnab)

       nlist = nlist + 1
       x=ccc(1,nlist)
       y=ccc(2,nlist)
       z=ccc(3,nlist)

       rsq=x*x+y*y+z*z

       if(rsq.lt.rcutsq)then

        it2 = it(j)
        dr=sqrt(rsq)


        potindex = ((it1-1)*(2*nt+2-it1))/2+it2-it1+1

        if(pt(potindex).eq.1)then

         etmp=(p1(potindex)/dr**12
     $             -p2(potindex)/dr**6)

         e = e + etmp
         force=(12*p1(potindex)/dr**12
     $          -6*p2(potindex)/dr**6)/rsq

         fff(1,i)=fff(1,i)+convfct2*force*x
         fff(2,i)=fff(2,i)+convfct2*force*y
         fff(3,i)=fff(3,i)+convfct2*force*z

         fff(1,j)=fff(1,j)-convfct2*force*x
         fff(2,j)=fff(2,j)-convfct2*force*y
         fff(3,j)=fff(3,j)-convfct2*force*z


        elseif(pt(potindex).eq.2)then

         etmp=(p1(potindex)
     $              *exp(-dr/p2(potindex))
     $              -p3(potindex)/dr**6)
         e = e + etmp

         force=(dr*p1(potindex)
     $         *exp(-dr/p2(potindex))
     $         /p2(potindex)-6*p3(potindex)/dr**6)
     $         /rsq

         fff(1,i)=fff(1,i)+convfct2*force*x
         fff(2,i)=fff(2,i)+convfct2*force*y
         fff(3,i)=fff(3,i)+convfct2*force*z

         fff(1,j)=fff(1,j)-convfct2*force*x
         fff(2,j)=fff(2,j)-convfct2*force*y
         fff(3,j)=fff(3,j)-convfct2*force*z

        endif

       endif


      end do
      end do

      return

      end

c $Id$
