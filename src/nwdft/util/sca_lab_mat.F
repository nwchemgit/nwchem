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
      
c            call ga_put(ig_a,
c     *           ifirsta,ilasta,
c     *           ifirstb,ilastb,
c     *           vec,nshbfa)
c            if(isha.ne.ishb)
c     *      call ga_put(ig_a,
c     *           ifirstb,ilastb,
c     *           ifirsta,ilasta,
c     *           vec,nshbfb)
      ioffseta = ifirsta - 1
      ioffsetb = ifirstb - 1
        ii=ifirsta
        call ga_put(ig_a,
     *       ii,ii,
     *       ifirstb,ioffsetb+nshbfb,
     *       vec,1)
        if(isha.ne.ishb)
     *  call ga_put(ig_a,
     *       ifirstb,ioffsetb+nshbfb,
     *       ii,ii,
     *       vec,1)
        if(nshbfa.gt.1) then
          n = nshbfb
          do 20 ma = 2,nshbfa
            ii=ii+1
            call ga_put(ig_a,
     *           ii,ii,
     *           ifirstb,ioffsetb+nshbfb,
     *           vec(n+1),1)
        if(isha.ne.ishb)
     *       call ga_put(ig_a,
     *           ifirstb,ioffsetb+nshbfb,
     *           ii,ii,
     *           vec(n+1),1)

            n=n+nshbfb
            
  20     continue
       endif


      return
      end

