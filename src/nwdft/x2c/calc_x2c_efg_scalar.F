      subroutine calc_x2c_efg_scalar (rtdb, nexc,
     &   g_dscf, g_dens_at, g_u)

c     -----------------------------------------------------------------
c     purpose: calculate one-electron property EFG using the 
c     two-component exact-decoupling Hamiltonian

c     For the time being we MUST use an uncontracted basis !

c     Integrals are supposed to be initialized by the calling
c     routine.

c     Driver (c) 2012, 2022 J. Autschbach, SUNY, jochena@buffalo.edu

c     Approach: The transformation matrices g_u come from calc_x2c_1e 
c     if certain options are set.
c     
c     -----------------------------------------------------------------

      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"
#include "apiP.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "bas.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "dftpara.fh"
#include "cdft.fh"
#include "geom.fh"
c#include "rel_consts.fh"
#include "cosmo.fh"

c     subroutine arguments:

      integer rtdb          ! [in] runtime database GA handle
      integer nexc          ! [in] no. of XC terms
      integer g_dscf(2)     ! [in] density matrix from SCF
      integer g_dens_at(2)  ! [in] atom densities for model pot.
      integer g_u(2)        ! [in] decoupling matrices

c     local variables:

      integer basis  
      integer nbas, nbas2
      integer me, type, i, j
      double precision sum
      character*2  symbol
      character*16 element, at_tag

c     GA handles used for temp storage in this routine:
      integer g_P(2), g_efgint, g_temp, g_efg, g_w(6)

c     dimensioning arrays
      integer alo(3), ahi(3), blo(3), bhi(3), clo(3), chi(3)

c     other MA handles
      integer l_xyzpt, k_xyzpt, l_zanpt, k_zanpt

c     bq variables (MV)
      logical dobq
      integer bq_ncent
      integer i_cbq
      integer i_qbq
      double precision elpotbq

      double precision ga_trace_diag
      external ga_trace_diag
      double precision tol2e, accval
      logical master, analytic_pvp, status
      integer dbg
      
      integer ndmem, ntemp
      integer iat, nat
      
      double precision xp, yp, zp, xn, yn, zn, zan
      double precision rr, rr5
      double precision efgel(6), efgnuc(6), etmp(6)
      double precision origin(3)

      double precision pifac, pi, deg, third
      parameter (third=1.0d0/3.0d0)
c      parameter(froth=4.0d0/3.0d0) 
c      parameter(pifac=froth*3.14159265358979323846264338327950288419d0)

c     =================================================================

c     ---------------
c     initializations
c     ---------------

c     options / settings

      dbg=0

      me=ga_nodeid()
      master = me.eq.0

      if (dbg>0.and.master) write(luout,*) 'entering x2c_efg'

      call ga_sync()

c     constants:

      tol2e=10.d0**(-itol2e)

      call int_acc_get(accval)

      if (dbg>0.and.master) write(luout,*) 'accval = ',accval
c      accval = 1d-25
c      call int_acc_set(accval)

      pi  = acos(-1.0d0)
      deg = 180.0d0/pi
      pifac = (4.0d0/3.0d0)*pi

      basis = ao_bas_han        ! AO basis handle

c     make sure the input GAs have the correct dimension

      call ga_inquire(g_u(1),type,nbas,nbas)
      if (nbas.ne.nbf_ao) then
        call errquit('x2c_efg: g_u has wrong dimension',nbas,GA_ERR)
      end if

      if (master.and.dbg>0) write (luout,*) 'x2c: nbas =',nbas
      nbas2 = nbas * nbas

c     in light of the recent updates of the x2c ham. code, it is 
c     important to check that we are working with an uncontracted
c     basis here. 
c     eventually, this will be fixed

      if (.not. bas_numbf(x2c_bas,ntemp))
     &  call errquit('x2c_efg: ntemp error',ntemp, BASIS_ERR)

      if (ntemp.ne.nbas)
     &  call errquit(
     &  'x2c_efg: cannot currently use a contracted basis',
     &  ntemp, BASIS_ERR)

c     -------------
c     Calculate EFG
c     -------------

c     initialize integrals that were not already initialized

c      call hnd_giao_init(basis,1)

c     ----- define points for calculation -----
c           1. nuclei: for now, we take ALL nuclei
c
      status=geom_ncent(geom,nat)
c
      if (.not. ma_push_get(mt_dbl,3*nat,'xyz pnt',l_xyzpt,k_xyzpt))
     &    call errquit('x2c_efg: ma failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,nat,'zan pnt',l_zanpt,k_zanpt))
     &    call errquit('x2c_efg: ma failed',911,MA_ERR)
c
      do  iat=1,nat
        status=geom_cent_get(geom,iat,at_tag,dbl_mb(k_xyzpt+3*(iat-1)),
     &     dbl_mb(k_zanpt+iat-1))
       end do

      alo(1) = nbas
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbas
      ahi(2) = nbas
      ahi(3) = 6*nat
      if (.not.nga_create(MT_DBL,3,ahi,'g_efgint',alo,g_efgint)) call
     &   errquit('x2c_efg: nga_create failed g_efgint',11,GA_ERR)
      call ga_zero(g_efgint)

c     calculate field gradient integrals (same as FC+SD in NMR ...)
      call int_giao_1ega(basis,basis,g_efgint,'sd+fc',
     &   dbl_mb(k_xyzpt),nat, oskel)

c     comment: The electronic integrals contain the isotropic
c     `contact' piece. However, when constructing the traceless
c     EFG tensor, those contributions disappear. We might as well
c     take them out explicitly

c     calculate FC integrals separately
      ahi(3) = nat
      if (.not.nga_create(MT_DBL,3,ahi,'temp',alo,g_temp)) call
     &   errquit('add_x2c: nga_create failed g_temp',31,GA_ERR)
      call ga_zero(g_temp)

      call int_giao_1ega(basis,basis,g_temp,'fc',
     &   dbl_mb(k_xyzpt),nat,oskel)

c
c     Take out FC part from SD+FC (i.e. from xx, yy and zz components)
c
      alo(1) = 1
      ahi(1) = nbas
      alo(2) = 1
      ahi(2) = nbas
      blo(1) = 1
      bhi(1) = nbas
      blo(2) = 1
      bhi(2) = nbas
      do iat = 1, nat
        alo(3) = 1+(iat-1)*6
        ahi(3) = 1+(iat-1)*6
        blo(3) = iat
        bhi(3) = iat
        do i = 1, 3 ! runs over xx, yy, zz
          call nga_add_patch(1.0d0,g_efgint,alo,ahi,
     &       pifac,g_temp,blo,bhi,
     &       g_efgint,alo,ahi)
          alo(3) = alo(3)+1
          ahi(3) = ahi(3)+1
         enddo
      enddo

      if (.not. ga_destroy(g_temp)) call errquit
     &   ('x2c_efg: temp ga corrupt?',31, GA_ERR)

c     the field gradient integrals carry a factor of 3
c     that we don't want:

      call ga_scale(g_efgint, third)

      if (dbg>1 .and.master) write (luout,*) 'x2c_efg before iat loop'

c     ---------------------------------------
c     calculate EFGs in a loop over the atoms
c     ---------------------------------------

      do iat = 1,nat

        xp = dbl_mb(k_xyzpt  +3*(iat-1))
        yp = dbl_mb(k_xyzpt+1+3*(iat-1))
        zp = dbl_mb(k_xyzpt+2+3*(iat-1))

        if (dbg>1 .and.master) write (luout,*) 'x2c_efg: iat=',iat

c       --------------------
c       nuclear contribution
c       --------------------

        do i = 1,6
          efgnuc(i) = 0d0
        end do

        do i = 1,nat
          xn  = dbl_mb(k_xyzpt  +3*(i-1)) - xp
          yn  = dbl_mb(k_xyzpt+1+3*(i-1)) - yp
          zn  = dbl_mb(k_xyzpt+2+3*(i-1)) - zp
          zan = dbl_mb(k_zanpt+i-1)
          rr = sqrt(xn*xn + yn*yn + zn*zn)
          if (rr.ge.1.0d-3) then 
c           (skip if it is the same nucleus or nuclei too close)
            rr5=rr*rr*rr*rr*rr
            efgnuc(1) = efgnuc(1) - zan*xn*xn/rr5
            efgnuc(2) = efgnuc(2) - zan*yn*yn/rr5
            efgnuc(3) = efgnuc(3) - zan*zn*zn/rr5
            efgnuc(4) = efgnuc(4) - zan*xn*yn/rr5
            efgnuc(5) = efgnuc(5) - zan*xn*zn/rr5
            efgnuc(6) = efgnuc(6) - zan*yn*zn/rr5
          end if 
        end do                  ! i = loop over all other nuclei

c       create traceless version of nuclear contribution to EFG:
        etmp(1) = 2.0d0*efgnuc(1) - efgnuc(2) - efgnuc(3)
        etmp(2) = 2.0d0*efgnuc(2) - efgnuc(1) - efgnuc(3)
        etmp(3) = 2.0d0*efgnuc(3) - efgnuc(1) - efgnuc(2)
        etmp(4) = 3.0d0*efgnuc(4)
        etmp(5) = 3.0d0*efgnuc(5)
        etmp(6) = 3.0d0*efgnuc(6)        
        do i = 1,6
          efgnuc(i) = etmp(i)
        end do
        
c       -----------------------
c       electronic contribution
c       -----------------------
        
c       First, we do the case without picture change corrections
        
c       copy electronic EFG integrals into temp array:
c       we also rearrange them such that we get the traceless
c       versions of the integrals right away
        
        alo(1) = nbas
        alo(2) = -1
        alo(3) = -1
        ahi(1) = nbas
        ahi(2) = nbas
        ahi(3) = 6
       
       if (.not.nga_create(MT_DBL,3,ahi,'efg',alo,g_efg)) call
     &     errquit('add_x2c: nga_create failed g_efg',32,GA_ERR)
        call ga_zero(g_efg)
 
        ahi(3) = 1
        if (.not.nga_create(MT_DBL,2,ahi,'temp',alo,g_temp)) call
     &     errquit('add_x2c: nga_create failed g_temp',31,GA_ERR)
        call ga_zero(g_temp)
               
        alo(1) = 1          
        ahi(1) = nbas
        alo(2) = 1           
        ahi(2) = nbas

        blo(1) = 1
        bhi(1) = nbas
        blo(2) = 1           
        bhi(2) = nbas

        clo(1) = 1
        chi(1) = nbas
        clo(2) = 1
        chi(2) = nbas
        clo(3) = 1
        chi(3) = 1

c       create traceless version of EFG matrices in g_efg
c       for the selected reference point

c       step 1: accumulate 3xx, 3yy, 3zz, 
c       and accumulate xx+yy+zz in g_temp
        do i = 1,3
          alo(3) = (iat-1)*6 + i 
          ahi(3) = alo(3)
          blo(3) = i
          bhi(3) = blo(3)
          call nga_add_patch(3.0d0,g_efgint,alo,ahi,
     &       1d0,g_efg,blo,bhi,
     &       g_efg,blo,bhi)
          call nga_add_patch(1.0d0,g_efgint,alo,ahi,
     &       1d0,g_temp,clo,chi,
     &       g_temp,clo,chi)          
        end do
c       step 2: subtract xx+yy+zz from diagonal 
        do i = 1,3
          blo(3) = i
          bhi(3) = blo(3)
          call nga_add_patch(1.0d0,g_efg,blo,bhi,
     &       -1d0,g_temp,clo,chi,
     &       g_efg,blo,bhi)
        end do

c       indices 1, 2, 3 of g_efg are now the components
c       3xx-rr, 3yy-rr, 3zz-rr

c       step 3: assign 3xy, 3xz, 3yz to the remaining elements 4, 5, 6
        do i = 4,6
          alo(3) = (iat-1)*6 + i 
          ahi(3) = alo(3)
          blo(3) = i
          bhi(3) = blo(3)
          call nga_add_patch(3.0d0,g_efgint,alo,ahi,
     &       1d0,g_efg,blo,bhi,
     &       g_efg,blo,bhi)
        end do

c       calculate electronic EFG components:

        alo(3) = 1 ! for density matrix
        ahi(3) = 1
               
        do i = 1,6              ! loop over EFG components
          
          if (dbg>1 .and.master) write (luout,*) 'x2c_efg: i=',i

          efgel(i) = 0d0

c         if we use array g_efg, use this:
          blo(3) = i            ! select this section from array g_efg
          bhi(3) = blo(3)

c         if we use g_efgint, use the boundaries below:
c$$$          blo(3) = i + (iat-1)*6
c$$$          bhi(3) = blo(3)

                   
c         tr[P*EFG]
          alo(1) = 1           
          ahi(1) = nbas
          alo(2) = 1
          ahi(2) = nbas
          call ga_zero(g_temp)
          call nga_matmul_patch('n','n',1d0,0d0,
     &       g_dscf(1),alo,ahi,
     &       g_efg, blo,bhi,
     &       g_temp,clo,chi)
          sum = ga_trace_diag(g_temp)
          
          efgel(i) = efgel(i) + sum
          
        end do                  ! i = terms in EFG

c       form the traceless EFG tensor elements
c       3x**2 - r**2, etc, for diagonal elements

c$$$        etmp(1) = 2.0d0*efgel(1) - efgel(2) - efgel(3)
c$$$        etmp(2) = 2.0d0*efgel(2) - efgel(1) - efgel(3)
c$$$        etmp(3) = 2.0d0*efgel(3) - efgel(1) - efgel(2)
c$$$        etmp(4) = 3.0d0*efgel(4)
c$$$        etmp(5) = 3.0d0*efgel(5)
c$$$        etmp(6) = 3.0d0*efgel(6)
c$$$
c$$$        do i = 1,6
c$$$          efgel(i) = etmp(i)
c$$$        end do

c        print results:
        
        if (master) then
          write (luout,*)
          write (luout,*) 'X2C EFG: Before picture-change corr.'
          write (luout,'(1x,a,i4,6e15.7)') 'x2cefg: inuc, efgnuc =',
     &       iat, (efgnuc(i), i=1,6)
          write (luout,'(1x,a,i4,6e15.7)') 'x2cefg: inuc, efgel  =',
     &       iat, (efgel(i), i = 1,6)
          write (luout,'(1x,a,i4,6e15.7)') 'x2cefg: inuc, efgtot =',
     &       iat, (efgel(i)+efgnuc(i), i = 1,6)
          write (luout,*)
        end if

c       ---------------------------------------------------
c       calculate picture-change corrected EFG operator
c       for the selected reference point (nuclear position)
c       ---------------------------------------------------

c       allocate arrays to hold the AO matrices of the EFG
c       operator in the small-component basis
c       (we'll calculate all 6 matrices in one go)

        alo(1) = nbas
        alo(2) = -1
        alo(3) = -1
        ahi(1) = nbas
        ahi(2) = nbas
        ahi(3) = 1
        do j = 1,6
          if (.not.nga_create(MT_DBL,2,ahi,'g_w',alo,g_w(j)))
     &       call errquit(
     &       'x2cefg: nga_create failed g_w',1+10*j,GA_ERR)  
          call ga_zero(g_w(j))
        end do
        
        origin(1) = xp
        origin(2) = yp
        origin(3) = zp
        if (dbg>0.and.master) write (luout,*) 'calling getw'
        call x2c_getw_efg_scalar(
     &     rtdb, g_dens_at, origin,
     &     g_w,
     &     nexc)
        if (dbg>0.and.master) write (luout,*) 'return from getw'
        

c       --------------------
c       allocate operator GA
c       --------------------

        if(.not.ga_create(mt_dbl, nbas, nbas,'g_P',0,0,
     &     g_P))
     &     call errquit(
     &     'x2c_efg: error creating g_P',666, GA_ERR)

c

        
c       ------------------------------------------------------------
c       in a loop over the EFG components, transform the 
c       operator to two-component form, and calculate the electronic
c       EFG component again
c       ------------------------------------------------------------
        
        do i = 1,6              ! loop over EFG components
          
          call ga_zero(g_P)

c         array bounds needed to extract EFG component from GA:

          alo(1) = 1
          ahi(1) = nbas
          alo(2) = 1
          ahi(2) = nbas
          alo(3) = i
          ahi(3) = i

c         g_temp always uses the same dimensions clo, chi
c         we use the same ones for g_u, g_P
          clo(1) = 1
          chi(1) = nbas
          clo(2) = 1
          chi(2) = nbas

c         transform EFG operator

c         UL' EFG UL, store in g_P

          call nga_matmul_patch('t','n',1d0,0d0,
     &      g_u(1),clo,chi,
     &      g_efg,alo,ahi,
     &      g_temp,clo,chi)
          
          call nga_matmul_patch('n','n',1d0,0d0,
     &      g_temp,clo,chi,
     &      g_u(1),clo,chi,
     &      g_P,clo,chi)
          
c         US' pEFGp US, add to g_P
          
          call nga_matmul_patch('t','n',1d0,0d0,
     &      g_u(2),clo,chi,
     &      g_w(i),clo,chi,
     &      g_temp,clo,chi)
          
          call nga_matmul_patch('n','n',1d0,1d0,
     &      g_temp,clo,chi,
     &      g_u(2),clo,chi,
     &      g_P,clo,chi)

c          call ga_print(g_P)

c         -----------------------------------------
c         calculate electronic EFG again, this time 
c         using the transformed operator
c         -----------------------------------------

          efgel(i) = 0.0d0
          
          alo(1) = 1                
          ahi(1) = nbas
          alo(2) = 1                
          ahi(2) = nbas
          blo(1) = 1               
          bhi(1) = nbas
          blo(2) = 1
          bhi(2) = nbas
          
c         tr[P*EFG]
          call ga_zero(g_temp)
          call nga_matmul_patch('n','n',1d0,0d0,
     &       g_dscf(1),alo,ahi,
     &       g_P,blo,bhi,
     &       g_temp,clo,chi)
          sum = ga_trace_diag(g_temp)
          efgel(i) = efgel(i) + sum          
          
        end do                  ! i = loop over EFG components

c       print results:
        
        if (master) then
          write (luout,*)
          write (luout,*) 'X2C EFG: After picture-change corr.'
          write (luout,'(1x,a,i4,6e15.7)') 'x2cefg: inuc, efgnuc =',
     &       iat, (efgnuc(i), i=1,6)
          write (luout,'(1x,a,i4,6e15.7)') 'x2cefg: inuc, efgel  =',
     &       iat, (efgel(i), i = 1,6)
          write (luout,'(1x,a,i4,6e15.7)') 'x2cefg: inuc, efgtot =',
     &       iat, (efgel(i)+efgnuc(i), i = 1,6)
          write (luout,*)
        end if

c       --------------------
c       clean up some GA memory
c       --------------------

        if (.not. ga_destroy(g_efg)) call errquit
     &     ('x2cefg: ga corrupt?',32, GA_ERR)
        if (.not. ga_destroy(g_temp)) call errquit
     &     ('x2cefg: ga corrupt?',31, GA_ERR)

        if (.not. ga_destroy(g_P)) call errquit
     &     ('x2cefg: ga P corrupt?',666, GA_ERR)

        do j = 1,6
          if (.not. ga_destroy(g_w(j))) call errquit
     &       ('x2cefg: failed destroying g_w',1+10*j, GA_ERR)
        end do

c       --------------------------------------------------------
c       end loop over refernce points at which EFG is calculated 
c       (typically the list of nuclei)
c       --------------------------------------------------------

      end do ! iat, loop over nuclei

      if (.not. ga_destroy(g_efgint)) call errquit
     &   ('x2cefg: ga corrupt?',11, GA_ERR)

c     --------
c     all done
c     --------

      call ga_sync()     

      if (dbg>0.and.master) write(luout,*) 'leaving x2c_efg'

c     =================================================================

      return
      end

