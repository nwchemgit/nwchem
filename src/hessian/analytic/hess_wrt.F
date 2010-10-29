      subroutine hess_wrt
*     $Id$
      implicit none
#include "mafdecls.fh"
#include "hess_info.fh"
c
      character*255 filehess,fileddipole
      logical dipole_okay
c
      dipole_okay=.true.
      call util_file_name('fd_ddipole', .false., .false.,fileddipole)
      call util_file_name('hess', .false., .false.,filehess)
      call stpr_wrt_fd_from_sq(dbl_mb(k_exy),n3xyz,filehess)
      if (dipole_okay) then
         call stpr_wrt_fd_dipole(dbl_mb(k_dder),ncent,fileddipole)
      else
         call util_file_unlink(fileddipole)
      endif
      return
      end
c
      subroutine hess_read(nxyz, e_xyz)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
c
      integer nxyz, nvar2, n3xyz
      integer l_exyt, k_exyt, nhesst, iii
      double precision dbl_tmp
      double precision e_xyz(nxyz*3,nxyz*3)
      character*255 filehess
c
c This routine reads the hessian from disk.  This was essentially 
c taken from stepper/vib_vib.F.
c
      n3xyz = 3 * nxyz
      nvar2 = n3xyz * n3xyz
      if (.not.ma_push_get(MT_DBL,nvar2,'hessian: triangle exy',
     &     l_exyt,k_exyt))
     &     call errquit('hess_read: cannot allocate',555, MA_ERR)
      nhesst = n3xyz*(n3xyz+1)/2
      call util_file_name('hess', .false., .false.,filehess)
      open(unit=69,file=filehess,form='formatted',status='old',
     &    err=99900,access='sequential')
      do iii = 0,(nhesst-1)
        read(69,*,err=99901,end=99902) dbl_tmp
        dbl_mb(k_exyt+iii) = dbl_tmp
      enddo
      close(unit=69,status='keep')
      call vib_dsquar(dbl_mb(k_exyt),e_xyz,n3xyz,n3xyz)
      if (.not.ma_pop_stack(l_exyt))
     *   call errquit('hess_read: cannot deallocate Exyt',555, MA_ERR)
c
      return
c
99900 continue
      write(luout,*)'hess_file => ',filehess
      call errquit('hess_read: error opening file: "filehess"',
     & 911,0)
99901 continue
      write(luout,*)'hess_file => ',filehess
      call errquit('hess_read: error reading file: "filehess"',
     & 911,0)
99902 continue
      write(luout,*)'hess_file => ',filehess
      call errquit
     & ('hess_read: unexpected EOF when reading file: "filehess"',
     & 911,0)
c
      end
c
      subroutine ddip_read(nxyz, ddipder)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "hess_info.fh"
#include "stdio.fh"
c
      integer nxyz, iii
      double precision dbl_tmp
      double precision ddipder(nxyz*9)
      character*255 fileddipole
      logical dipole_file_exists
c
c This routine reads the dipole derivative terms
c from disk.  This was essentially taken from stepper/vib_vib.F.
c
      call util_file_name('fd_ddipole',.false., .false.,fileddipole)
      dipole_file_exists = .false.
      inquire(file=fileddipole,exist=dipole_file_exists)
      if (dipole_file_exists) then
        open(unit=70,file=fileddipole,form='formatted',status='old',
     &      err=89900,access='sequential')
        do iii = 1, 9*nxyz
          read(70,*,err=89901,end=89902) dbl_tmp
          ddipder(iii) = dbl_tmp
        enddo
        close(unit=70,status='keep')
      endif
c
      return
c
89900 continue
      write(luout,*)'dipole_file => ',fileddipole
      call errquit('hess_read: error opening file: "fileddipole"',811,
     &       DISK_ERR)
89901 continue
      write(luout,*)'dipole_file => ',fileddipole
      call errquit('hess_read: error reading file: "fileddipole"',811,
     &       DISK_ERR)
89902 continue
      write(luout,*)'dipole_file => ',fileddipole
      call errquit
     & ('hess_read: unexpected EOF when reading file: "fileddipole"',
     &   811, DISK_ERR)
c
      end
