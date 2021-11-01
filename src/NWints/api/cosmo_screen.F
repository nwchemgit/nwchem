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
      parameter (LOGPIx15=1.7170948287741d0)
      parameter (N3OVER2=-1.5d0)
      double precision eps_small
      parameter (eps_small=1d-32)
      double precision const
      double precision qmax

      integer mpa,mpb,mca,mcb
      double precision a,b,abi,beta,r2,over,ac,bc
      double precision est_int(nca,ncb)

      qmax = maxval(abs(zan))

      if(val_int_acc.lt.eps_small)  then
         const = LOGPIx15-log(eps_small)+log(nefc*qmax/3d0)
      else
         const = LOGPIx15-log(val_int_acc)+log(nefc*qmax/3d0)
      endif

      R2 = (Axyz(1) - Bxyz(1))**2 +
     $     (Axyz(2) - Bxyz(2))**2 +
     $     (Axyz(3) - Bxyz(3))**2

      skip = .true.

      est_int(:,:) = 0d0
      do mpa = 1,NPA
        a = Aprims(mpa)
        do mpb = 1,NPB
          b = Bprims(mpb)
         
          abi = 1.0d0/(a+b)
          beta = a*b*abi

          over = exp(-beta*r2)*abi**n3over2

          do mcb=1,ncb
            bc = bcoefs(mpb,mcb)
            if (bc.eq.0d0) cycle
            do mca=1,nca
              ac = acoefs(mpa,nca)
              if (ac.eq.0d0) cycle
              est_int(mca,mcb) = est_int(mca,mcb) + ac*bc*over
            enddo
          enddo
         
        enddo
      enddo

      do mcb=1,ncb
        do mca=1,nca
          if (abs(est_int(mca,mcb)*qmax*0.5d0*acos(-1d0)).gt.
     $        val_int_acc) skip = .false.     
        enddo
      enddo

      end subroutine
