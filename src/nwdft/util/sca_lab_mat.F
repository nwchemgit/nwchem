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
c     removed <AO_basis_set.fh>
c     removed <CD_basis_set.fh>
c     removed <XC_basis_set.fh>

      Dimension vec(*)
c     
c     Given a handle to a Gaussian basis set and two shell indices, return the 
c     composite index for the pair of basis functions.
c     
c******************************************************************************
c**** 
c**** if( name.eq.'AO' )then
c**** ioffseta = kloc0(isha) - 1
c**** ioffsetb = kloc0(ishb) - 1
c**** else if( name.eq.'CD' )then
c**** ioffseta = kloc0_cd(isha) - 1
c**** ioffsetb = kloc0_cd(ishb) - 1
c**** else if( name.eq.'XC' )then
c**** ioffseta = kloc0_xc(isha) - 1
c**** ioffsetb = kloc0_xc(ishb) - 1
c**** end if

      if (.not. bas_cn2bfr( basis_hand, isha, ifirsta, ilast)) then
        call errquit('Exiting in sca_lab_mat',1)
      endif 
      ioffseta = ifirsta - 1
      if (.not. bas_cn2bfr( basis_hand, ishb, ifirstb, ilast)) then
        call errquit('Exiting in sca_lab_mat',1)
      endif 
      ioffsetb = ifirstb - 1

c**** nshbfa = igb_shell_bf(name,isha)
      status = bas_continfo(basis_hand, isha,
     &     idum1, idum2, nshbfa, idum3)

      if( isha.eq.ishb )then

        ic=ioffseta+1
        call ga_put(ig_a,
     *       ic,ic,
     *       ic,ic,
     *       vec,1)
        if(nshbfa.gt.1) then
          ii=ic
          n = 1
          do 10 ma = 2,nshbfa
            ii=ii+1
            call ga_put(ig_a,
     *           ii,ii,
     *           ic,ii,
     *           vec(n+1),1)
            call ga_put(ig_a,
     *           ic,ii-1,
     *           ii,ii,
     *           vec(n+1),1)
            n = n + ma

   10     continue
        endif


      else 

c****   nshbfb = igb_shell_bf(name,ishb)
        status = bas_continfo(basis_hand, ishb,
     &       idum1, idum2, nshbfb, idum3)
        ii=ioffseta+1
        call ga_put(ig_a,
     *       ii,ii,
     *       ioffsetb+1,ioffsetb+nshbfb,
     *       vec,1)
        call ga_put(ig_a,
     *       ioffsetb+1,ioffsetb+nshbfb,
     *       ii,ii,
     *       vec,1)
        if(nshbfa.gt.1) then
          n = nshbfb
          do 20 ma = 2,nshbfa
            ii=ii+1
            call ga_put(ig_a,
     *           ii,ii,
     *           ioffsetb+1,ioffsetb+nshbfb,
     *           vec(n+1),1)
            call ga_put(ig_a,
     *           ioffsetb+1,ioffsetb+nshbfb,
     *           ii,ii,
     *           vec(n+1),1)

            n=n+nshbfb
            
   20     continue
        endif


      end if

      end

