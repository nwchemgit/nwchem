      Subroutine nbf_to_mbf(  geom_hand, basis_hand,
     &     rq0,mbf,mcenters,
     &     zprim, nzprim)

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      integer basis_hand, geom_hand

      Logical LResult

#include <dimensions.fh>
#include <bas.fh>

c--   > Distance Squared between Sampling Points and Centers

      Dimension rq0(mcenters)

      Dimension zprim(nzprim)

      if( .not. bas_numcont(basis_hand,ncontrset) )
     &          call errquit('Exiting in int_eval_gbsets.',5)

        n1 = 0
        n2 = 0
        do 400 icset = 1,ncontrset

          LResult = bas_cn2ce( basis_hand, icset, ictr )
c         
c         get info about current contraction set
c         
          LResult = bas_continfo(basis_hand, icset,
     &         itype, nprimo, ncontr, isphere)
            if (.not. bas_cn2bfr( basis_hand,icset,
     &                      ifirst,ilast))
     &        call errquit('Exiting in int_eval_gbsets.',33)
                nshbf=ilast-ifirst+1
c         
c         angular momentum
c         
          if (itype .lt. 0) then
            call errquit(
     &           'int_eval_gbsets: sp-type orbital not coded', 1)
          else 
            l = itype
          endif 
c         
c         cartesian/spherical harmonic
c         
          if ( isphere .eq. 0 ) then !  cartesian set
            nang = (l+1)*(l+2)/2
          elseif ( isphere .eq. 1 ) then !  spherical harmonic
            nang = 2*l+1
          else 
            call errquit( 'int_cgtf:  illegal isphere value', 1 )
          endif 
c         
c         get exponents and contraction coefficients for this contraction set
c         
          LResult = bas_get_exponent( basis_hand, icset, zprim )

c         Determine the minimum Gaussian exponent.

          zmin = 1.D+06
          do 300 iprimo = 1,nprimo
            zmin = min(zprim(iprimo),zmin)
  300     continue

c         Only include those basis functions that are "non-zero" for at least one 
c         point in the sampling set.

          if( zmin*rq0(ictr).gt.2.D+01 )then

            n2 = n2 + nshbf
            go to 400

          else

            n2 = n2 + nshbf

          end if

          if( l.eq.0 )then

c           =============>  S Contractions  <============= 

            n1 = n1 + ncontr

          else if( l.eq.1 )then

c           =============>  P Contractions  <============= 

            n1 = n1+ncontr*3

          else if( l.eq.2 )then

c           =============>  D Contractions  <============= 

            n1 = n1 + ncontr*6

          else

c           =============>  General Case  <============= 

            n1 = n1 + ncontr*nang

          end if

  400   continue

        mbf = n1
        return
      end
