      subroutine dimqm_Dens2Dipoles(rtdb, nbf, geom, basis,
     $                              g_dens_r, g_dens_im, lstatic,
     $                              lifetime, dir)
c
c     Determines the induced dipoles for the given density matrix 
c    
c     Called from aoresponse.F
c
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "dimqm_constants.fh"
#include "dimqm.fh"
c
c     Input Variables
      integer rtdb                ! RTDB handle
      integer g_dens_r            ! density matrix handle (real part)
      integer g_dens_im           ! density matrix handle (imaginary part)
      integer nbf                 ! number of basis functions
      integer geom                ! geometry handle
      integer basis               ! basis set handle
      logical lstatic   ! static or FD?
      logical lifetime   ! Lifetime?
      integer dir
c
c     Local variables
      integer i3
      integer icmplx
      integer l_fld, k_fld
      integer l_dimxyz, k_dimxyz
      integer l_muind, k_muind
      character*1 direction(3)
      data direction /'x', 'y', 'z'/
c
c     Common variables used
c
c     logical lcmplx
      if(ldebug) then
        write(luout,*) "Start Dens2Dipoles"
      end if
      i3 = nDIM * 3
      icmplx = COMPLEXMATH
      if(lstatic) icmplx = REALMATH
c      if(lfreq .and. lcmplx) icmplx = COMPLEXMATH
c      call ga_print(g_dens_r)
c
c     Allocate memory for arrays
      if(.not.ma_push_get(mt_dbl,i3*icmplx,'dimrsp:fld',l_fld,k_fld))
     $  call errquit('malloc dimrsp:fld failed',1,MA_ERR)
c
      if(.not.ma_push_get(mt_dbl,i3*icmplx,'dimrsp:muind',
     $                                            l_muind,k_muind))
     $  call errquit('malloc dimrsp:muind failed',1,MA_ERR)
c
      if(.not.ma_push_get(mt_dbl,i3*icmplx,'dimrsp:muold',
     $                                            l_muold,k_muold))
     $  call errquit('malloc dimrsp:muold failed',1,MA_ERR)
c
c      if(.not.ma_push_get(mt_dbl,i3*icmplx,'dimrsp:muold',
c     $                                            l_muold,k_muold))
c     $  call errquit('malloc dimrsp:muold failed',1,MA_ERR)
c
      if(.not.ma_push_get(mt_dbl,i3,'dimrsp:xyz',l_dimxyz,k_dimxyz))
     $  call errquit('malloc dimrsp:xyz failed',1,MA_ERR)
c      write(luout,*) "Arrays Allocated"
c
c     Initizalize arrays
      call dfill(i3*icmplx, ZERO, dbl_mb(k_muind), 1)
      call dfill(i3*icmplx, ZERO, dbl_mb(k_muold), 1)
      if(.not.rtdb_get(rtdb,'dimpar:coords',mt_dbl,i3,
     $                                dbl_mb(k_dimxyz)))
     $  call errquit('get dimpar:coords failed', 1, RTDB_ERR)
c
c       call dimqm_screening(rtdb, geom, basis, dbl_mb(k_dimxyz),
c     $                      screen)      
c      write(luout,*) "Arrays initialized"
c
c     Initialize the integral routines since they get close out earlier
c      call int_init(rtdb,1,basis)
c      write(luout,*) "Integrals initialized"
c
c     Calculate the real part of the field
      call dimqm_EqmE(rtdb, g_dens_r, geom, basis, dbl_mb(k_fld),
     $                                   dbl_mb(k_dimxyz))
c      write(luout,*) "Real Field Done"
c
c     Calculate the imaginary part of the field (maybe?)
      if(lifetime) then
c        write(luout,*) "Complex Field Needed"
        call dimqm_EqmE(rtdb, g_dens_im, geom, basis,
     $             dbl_mb(k_fld+i3), dbl_mb(k_dimxyz))
      endif
c      write(luout,*) "Field Created"
c
c     Close the integral routines again
c      call int_terminate()
c      write(luout,*) "Integrals closed"
c
c     Calculate the dipoles induced by the QM density
c      lfirst = .true.
      call dimqm_f2d(rtdb, dbl_mb(k_fld), 
     $               dbl_mb(k_muind), dbl_mb(k_muold),
     $           dbl_mb(k_dimxyz), icmplx, direction(dir), ' ',.false.)
c      write(luout,*) "Dipoles Created"
c
c     Write induced dipoles to RTDB
      if(.not.rtdb_put(rtdb,'dimqm:muind',mt_dbl,i3*icmplx,
     $                                        dbl_mb(k_muind)))
     $  call errquit('put dimqm:muind failed',1,RTDB_ERR)
c
c      write(luout,*) "Dipoles Written"
c
c     Print total dipole moment
      call dimqm_printDipole(rtdb, dbl_mb(k_muind), 
     $                       icmplx, .false.)
c     Clean up arrays
      if(.not.ma_chop_stack(l_fld))
     $  call errquit('ma chop at l_fld failed',1,MA_ERR)
      if(ldebug) then
        write(luout,*) "End Dens2Dipoles"
      end if
c
      end subroutine dimqm_Dens2Dipoles
