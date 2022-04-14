      subroutine cosmo_screen(
     &               Axyz,Aprims,Acoefs,NPA,NCA,La,
     &               Bxyz,Bprims,Bcoefs,NPB,NCB,Lb,
     &               Cxyz, zan, nefc, skip )

      Implicit none

      double precision pi,explim
      Parameter (PI=3.1415926535898d0)
      Parameter (EXPLIM=100.d0)
#include "apiP.fh"

      Logical skip

      integer npa,nca,la,npb,ncb,lb,nefc
      double precision axyz(3),bxyz(3),cxyz(3,nefc)
      double precision aprims(npa),acoefs(npa,nca)
      double precision bprims(npb),bcoefs(npb,ncb)
      double precision zan(nefc)

      double precision n3over2,logpix15
      double precision eps_small
      parameter (eps_small=1d-32)
      double precision const

      integer mpa,mpb
      double precision a,b,abi,beta,r2,over

      if(val_int_acc.lt.eps_small)  then
         const = log(pi)-log(eps_small)
      else
         const = log(pi)-log(val_int_acc)
      endif

      R2 = (Axyz(1) - Bxyz(1))**2 +
     $     (Axyz(2) - Bxyz(2))**2 +
     $     (Axyz(3) - Bxyz(3))**2

      skip = .true.

      do mpa = 1,NPA
        a = Aprims(mpa)
        do mpb = 1,NPB
          b = Bprims(mpb)
         
          abi = 1.0d0/(a+b)
          beta = a*b*abi

          if (beta*r2 .lt. const+log(abi)) then
            skip = .false.
            return
          endif

        enddo
      enddo

      end subroutine
