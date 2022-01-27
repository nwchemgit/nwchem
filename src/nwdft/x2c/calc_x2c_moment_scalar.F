      subroutine calc_x2c_moment_scalar (rtdb, nexc,
     &   g_dscf, g_dens_at, g_u, type, ncomp)

c     -----------------------------------------------------------------
c     purpose: calculate some multipole moments using the 
c     two-component exact-decoupling Hamiltonian

c     For the time being we MUST use an uncontracted basis !

c     Integrals are supposed to be initialized by the calling
c     routine.

c     Driver (c) 2012,2022 J. Autschbach, SUNY, jochena@buffalo.edu

c     Approach: The transformation matrices g_u come from calc_x2c_1e 
c     if certain options are set.
c     
c     -----------------------------------------------------------------

      implicit none
c
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
#include "msgids.fh"
#include "dra.fh"
#include "inp.fh"
c
c     subroutine arguments:

      integer rtdb          ! [in] runtime database GA handle
      integer nexc          ! [in] no. of XC terms
      integer g_dscf(2)     ! [in] density matrix from SCF
      integer g_dens_at(2)  ! [in] atom densities for model pot.
      integer g_u(2)        ! [in] decoupling matrices
      integer ncomp         ! [in] number of components for moment
      character*(*) type    ! type of moment

c     local variables:

      integer basis  
      integer nbas, nbas2
      integer me, idum, i, j
      double precision sum
      character*2  symbol
      character*16 element, at_tag

c     GA handles used for temp storage in this routine:
      integer g_P, g_int, g_temp, g_moment, g_w(ncomp)

c     dimensioning arrays
      integer alo(3), ahi(3), blo(3), bhi(3), clo(3), chi(3)



      double precision ga_trace_diag
      external ga_trace_diag
      double precision tol2e, accval
      logical master, analytic_pvp, status
      integer dbg
      
      integer ndmem, ntemp

      double precision moment(ncomp)
      double precision origin(3)
c
      integer i_comp
      character*3 ch_comp
      logical dosavepcdip
      character*(nw_max_path_len) fn_pcdip
c
      logical dmat_to_file
      external dmat_to_file
c
c     =================================================================

c     ---------------
c     initializations
c     ---------------

c     options / settings

      dbg=1

      me=ga_nodeid()
      master = me.eq.0

      if (dbg>0.and.master) write(luout,*) 'entering x2c_moment'

      call ga_sync()
c
      dosavepcdip = .false.
      if(.not.(rtdb_get(rtdb,'x2c:savepcdip',mt_log,1,dosavepcdip)))
     &      dosavepcdip = .false.

c     sanity check: type of integral, and dimensioning

      if (type.eq.'dipole') then
        if (ncomp.ne.3)  call errquit
     &     ('x2c_moment: ncomp has wrong value',666,INPUT_ERR)  
        if(master) write(luout,*) 'x2c: dipole moment calculation'
      else 
        if (master) write(luout,*)
     &     'x2c_moment: integral type not supported. Aborting'
        return
      end if

      tol2e=10.d0**(-itol2e)

      call int_acc_get(accval)

      if (dbg>0.and.master) write(luout,*) 'accval = ',accval
c      accval = 1d-25
c      call int_acc_set(accval)


      basis = ao_bas_han        ! AO basis handle

c     make sure the input GAs have the correct dimension

      call ga_inquire(g_u(1),idum,nbas,nbas)
      if (nbas.ne.nbf_ao) then
        call errquit('x2c_moment: g_u has wrong dimension',nbas,GA_ERR)
      end if

      if (master.and.dbg>0) write (luout,*) 'nbas =',nbas
      nbas2 = nbas * nbas

c     in light of the recent updates of the x2c ham. code, it is 
c     important to check that we are working with an uncontracted
c     basis here. 
c     eventually, this will be fixed

      if (.not. bas_numbf(x2c_bas,ntemp))
     &  call errquit('x2c_moment: ntemp error',ntemp, BASIS_ERR)

      if (ntemp.ne.nbas)
     &  call errquit(
     &  'x2c_moment: cannot currently use a contracted basis',
     &  ntemp, BASIS_ERR)


c     coordinate origin

      origin(1) = 0d0
      origin(2) = 0d0
      origin(3) = 0d0

c     --------------------------
c     Calculate electric moments
c     --------------------------

c     initialize integrals that were not already initialized

c      call hnd_giao_init(basis,1)

c     allocate integral array and calculate AO matrix.
c     the integral type in variable 'type' must match
c     the available types and no. of components in 
c     routine int_mpole_1ega

      alo(1) = nbas
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbas
      ahi(2) = nbas
      ahi(3) = ncomp
      if (.not.nga_create(MT_DBL,3,ahi,'g_int',alo,g_int)) call
     &   errquit('x2c_moment: nga_create failed g_int',11,GA_ERR)
      call ga_zero(g_int)
      call int_mpole_1ega(basis,basis,g_int,type,origin,
     &     oskel)

      if(.not.ga_create(mt_dbl, nbas, nbas,'g_temp',0,0,g_temp))
     &   call errquit(
     &   'x2c_moment: error creating g_temp',12, GA_ERR)
      
c     First, we do the case without picture change corrections
                              
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
      
c     calculate electronic moment components:
      
      alo(3) = 1                ! for density matrix
      ahi(3) = 1
      
      do i = 1,ncomp            ! loop over moment components
        
        moment(i) = 0.0d0
        
c       
        blo(3) = i              ! select this section from array g_int
        bhi(3) = blo(3)        

c       tr[P*Moment]
        moment(i) = 0.0d0
        
        call ga_zero(g_temp)
        call nga_matmul_patch('n','n',1d0,0d0,
     &     g_dscf(1),alo,ahi,
     &     g_int, blo,bhi,
     &     g_temp,clo,chi)
        sum = ga_trace_diag(g_temp)
        
        moment(i) = moment(i) + sum
        
      end do                    ! i = terms in moment

c     print results:
      
      if (master) then
        write (luout,*)
        write (luout,*) 'X2C Moment: Before picture-change corr.'
        write (luout,'(1x,a,50e15.7)') 'moment = ',
     &     (moment(i), i=1,ncomp)
        write (luout,*)
      end if
      
c     allocate arrays to hold the AO matrices of the EFG
c     operator in the small-component basis
c     (we'll calculate all matrices in one go)
      
      alo(1) = nbas
      alo(2) = -1
      alo(3) = -1
      ahi(1) = nbas
      ahi(2) = nbas
      ahi(3) = 1
      do j = 1,ncomp
        if (.not.nga_create(MT_DBL,2,ahi,'g_w',alo,g_w(j)))
     &     call errquit(
     &     'x2cefg: nga_create failed g_w',1+10*j,GA_ERR)  
        call ga_zero(g_w(j))
      end do
      
      if (dbg>0.and.master) write (luout,*) 'calling getw'
      call x2c_getw_moment_scalar(
     &   rtdb, g_dens_at, origin,
     &   g_w,
     &   nexc, type, ncomp)
      if (dbg>0.and.master) write (luout,*) 'return from getw'
      
      
c     --------------------
c     allocate operator GA
c     --------------------
      
      if(.not.ga_create(mt_dbl, nbas, nbas,'g_P',0,0,
     &   g_P))
     &   call errquit(
     &   'x2c_efg: error creating g_P',666, GA_ERR)
      
      
c     ------------------------------------------------------------
c     in a loop over the moment components, transform the 
c     operator to two-component form, and recalculate moment
c     ------------------------------------------------------------
      
      do i = 1,ncomp            ! loop over components
        
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

c         transform moment operator

c         UL' MomenT UL, store in g_P

          call nga_matmul_patch('t','n',1d0,0d0,
     &      g_u(1),clo,chi,
     &      g_int,alo,ahi,
     &      g_temp,clo,chi)
          
          call nga_matmul_patch('n','n',1d0,0d0,
     &      g_temp,clo,chi,
     &      g_u(1),clo,chi,
     &      g_P,clo,chi)
          
c         US' pMomenTp US, add to g_P
          
          call nga_matmul_patch('t','n',1d0,0d0,
     &      g_u(2),clo,chi,
     &      g_w(i),clo,chi,
     &      g_temp,clo,chi)
          
          call nga_matmul_patch('n','n',1d0,1d0,
     &      g_temp,clo,chi,
     &      g_u(2),clo,chi,
     &      g_P,clo,chi)

c          call ga_print(g_P)        

c       ---------------------------------
c       calculate moment again, this time 
c       using the transformed operator
c       ---------------------------------
        
        moment(i) = 0.0d0
        
        alo(1) = 1                
        ahi(1) = nbas
        alo(2) = 1                
        ahi(2) = nbas
        blo(1) = 1               
        bhi(1) = nbas
        blo(2) = 1
        bhi(2) = nbas
        
c       tr[P*Moment]
        call ga_zero(g_temp)
        call nga_matmul_patch('n','n',1d0,0d0,
     &     g_dscf(1),alo,ahi,
     &     g_P,blo,bhi,
     &     g_temp,clo,chi)
        sum = ga_trace_diag(g_temp)
        moment(i) = moment(i) + sum     

c       The GA g_P now has the picture-change transformed
c       electric moment AO matrix for component no. i
c       (as returned from int_mpole_1ega).
cng     check if we need to save the g_P matrices for use in TDDFT
c
        if (dosavepcdip) call write_x2c_pc_dip(i,g_P) ! save components
c
      end do                    ! i = loop over moment components
      
c     print results:
      
      if (master) then
        write (luout,*)
        write (luout,*) 'X2C Moment: After picture-change corr.'
        write (luout,'(1x,a,50e15.7)') 'moment = ',
     &     (moment(i), i=1,ncomp)
        write (luout,*)
      end if

c     --------------------
c     clean up some GA memory
c     --------------------
     
      
      if (.not. ga_destroy(g_P)) call errquit
     &   ('x2c_moment: ga P corrupt?',666, GA_ERR)
      
      do j = 1,ncomp
        if (.not. ga_destroy(g_w(j))) call errquit
     &     ('x2c_moment: failed destroying g_w',1+10*j, GA_ERR)
      end do

      if (.not. ga_destroy(g_temp)) call errquit
     &   ('x2c_moment: ga corrupt?',12, GA_ERR)
      
      
      if (.not. ga_destroy(g_int)) call errquit ('x2cefg: ga corrupt?'
     &   ,11, GA_ERR)

      
c     --------
c     all done
c     --------

      call ga_sync()     

      if (dbg>0.and.master) write(luout,*) 'leaving x2c_moment'

c     =================================================================

      return
      end

