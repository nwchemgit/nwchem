      Subroutine sca_lab_mat(basis_hand,isha,ishb,
     &     nshbfa,nshbfb,ig_a,vec,vec2)

C$Id: sca_lab_mat.F,v 1.6 1995-02-02 23:20:08 d3g681 Exp $
      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      logical status
c**** 
c**** nwchem basis handle
c**** 
      integer basis_hand
c**** 

#include <bas.fh>
 
      Dimension vec(nshbfb,*),vec2(nshbfa,*)
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

       call ga_acc(ig_a,
     *     ifirstb,ilastb,
     *     ifirsta,ilasta,
     *     vec,nshbfb,1.d0)
      if(isha.ne.ishb) then
         do jj=1,nshbfb
            do ii=1,nshbfa
               vec2(ii,jj)=vec(jj,ii)
            enddo
         enddo
         call ga_acc(ig_a,
     *        ifirsta,ilasta,
     *        ifirstb,ilastb,
     *        vec2,nshbfa,1.d0)

      endif


      return
      end

