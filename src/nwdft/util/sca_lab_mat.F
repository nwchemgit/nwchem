      Subroutine sca_lab_mat(basis_hand,isha,ishb,ig_a,vec)

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      logical status
c**** 
c**** nwchem basis handle
c**** 
      integer basis_hand
c**** 

#include <bas.fh>
#include <dimensions.fh>
 
      Dimension vec(*)
c     
c     Given a handle to a Gaussian basis set and two shell indices, return the 
c     composite index for the pair of basis functions.
c     
c******************************************************************************

      if (.not. bas_cn2bfr( basis_hand, isha, ifirsta, ilasta)) then
        call errquit('Exiting in sca_lab_mat',1)
      endif 
      if (.not. bas_cn2bfr( basis_hand, ishb, ifirstb, ilastb)) then
        call errquit('Exiting in sca_lab_mat',2)
      endif 

      nshbfa=ilasta-ifirsta+1

      nshbfb=ilastb-ifirstb+1
      
      call ga_acc(ig_a,
     *     ifirstb,ilastb,
     *     ifirsta,ilasta,
     *     vec,nshbfb,1.d0)
      if(isha.ne.ishb) then
        n=1
        do ii=ifirsta,ilasta
             call ga_acc(ig_a,
     *       ii,ii,
     *       ifirstb,ilastb,
     *       vec(n),1,1.d0)
          n=n+nshbfb
        enddo
      endif


      return
      end

