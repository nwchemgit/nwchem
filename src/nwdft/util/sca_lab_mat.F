      Subroutine sca_lab_mat(basis_hand,isha,ishb,
     &     nshbfa,nshbfb,g_a,vec,vec2)

C$Id: sca_lab_mat.F,v 1.9 1996-10-18 21:50:15 d3h449 Exp $
      implicit none
      
      integer basis_hand
      integer isha,ishb
      integer nshbfa,nshbfb
      integer g_a
      double precision vec(nshbfb,*),vec2(nshbfa,*)
c**** 

#include "bas.fh"
c 
      integer ifirsta, ilasta, ifirstb, ilastb,jj,ii
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

       call ga_acc(g_a,
     *     ifirstb,ilastb,
     *     ifirsta,ilasta,
     *     vec,nshbfb,1.d0)
      if(isha.ne.ishb) then
         do jj=1,nshbfb
            do ii=1,nshbfa
               vec2(ii,jj)=vec(jj,ii)
            enddo
         enddo
         call ga_acc(g_a,
     *        ifirsta,ilasta,
     *        ifirstb,ilastb,
     *        vec2,nshbfa,1.d0)

      endif


      return
      end

