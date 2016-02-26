c
c     Saves the solution of the linear systems of equations in its raw form.
c     HvD 10/2007
c
      subroutine tddft_grad_solve_save_x(g_x,ipol,nao,nfc,naoc,nocc,
     +                                   nav,nfv,nroot,g_mo,fn_x,x_when,
     +                                   geom,basis,oskel)
c
c     This routine saves the solution of the linear systems of equations
c     in its raw form. The purpose of this is that it may be beneficial
c     to use this solution as a starting point for subsequent linear
c     system solves in a geometry optimisation or a numerical hessian
c     evaluation.
c
c     The data is stored in AO basis to make it less susceptible to
c     changes in the MOs as the geometry changes. The AO basis data
c     needs to be multiplied by the overlap matrix to allow the 
c     transformation back to MO basis to succeed. This multiplication
c     is done here because the number of saves is always less than or
c     equal to the number of loads.
c
c     The solution g_x is saved only if x_when equals "once" or 
c     "always". The other legal value of x_when is "never", all other
c     values are considered to be an error. The management of x_when
c     to control the appropriate time and frequency of saving g_x
c     is not taken care of in this routine but is expected to happen
c     elsewhere.
c
c     Also any failure to produce the file fn_x is considered to be
c     an error.
c
c     Finally, this routine is a companion of tddft_grad_solve_load_x
c     which reads the data produced by this routine. So any changes
c     need to be consistent across both routines.
c
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
c
c     Input:
c
      integer g_x          ! 2D global array for linear system solutions
      integer ipol         ! the number of spin components
      integer nao          ! the number of AO's
      integer nfc(ipol)    ! the number of frozen core orbitals
      integer naoc(ipol)   ! the number of active occupied orbitals
      integer nocc(ipol)   ! the number of occupied orbitals
      integer nav(ipol)    ! the number of active virtual orbitals
      integer nfv(ipol)    ! the number of active frozen virtual
                           ! orbitals
      integer nroot        ! the number of columns in g_x
      integer g_mo(ipol)   ! the molecular orbitals
c
      character*(*) fn_x   ! filename for solutions file
      character*(*) x_when ! when is g_x to be saved?
c
      integer geom         ! the geometry handle
      integer basis        ! the basis set handle
      logical oskel        ! 
c
c     Local:
c
      integer maxi                     ! the leading dimension of the
                                       ! buffer
      parameter(maxi=100)
      double precision buf(maxi) ! a local buffer for the data
      integer i, j, k
      integer ir                 ! counter over roots
      integer ip, ipm1           ! counter over spin components
      integer vlo(3), vhi(3)     ! patch limits
      integer mlo(3), mhi(3)     ! patch limits
      integer g_mat_mo(ipol) ! temporary global array in MO basis
      integer g_mat_ao       ! temporary global array in AO basis
      integer g_ovlp         ! overlap matrix
      integer g_work         ! scratch matrix
c
      integer idim(3) ! dimensions of temporary global arrays
      integer ichk(3) ! chunking for matrix distribution
c
      character*132 msg          ! error message
      character*32 pname
c
      pname = "tddft_grad_solve_save_x: "
c
c     Code:
c
      if (x_when.eq."never") then
c
c       Nothing to do here.
c
      elseif (x_when.eq."once".or.x_when.eq."always") then
c
c       Save g_x.
c
        call ga_sync
        do ip = 1, ipol
          idim(1) = 1
          idim(2) = naoc(ip)
          idim(3) = nav(ip)
          ichk(1) = 1
          ichk(2) = -1
          ichk(3) = -1
          if (.not.nga_create(mt_dbl,3,idim,'g_mat_mo',ichk,
     +        g_mat_mo(ip))) 
     +      call errquit(pname//'failed to create g_mat_mo',0,GA_ERR)
        enddo
        idim(1) = ipol
        idim(2) = nao
        idim(3) = nao
        ichk(1) = ipol
        if (.not.nga_create(mt_dbl,3,idim,'g_mat_ao',ichk,
     +      g_mat_ao)) 
     +    call errquit(pname//'failed to create g_mat_ao',0,GA_ERR)
c
        if (ga_nodeid().eq.0) then
          open(unit=69,file=fn_x,form='unformatted',err=90000)
          write(69,err=90010)ipol,nao,
     +         (nfc(i),i=1,ipol),
     +         (naoc(i),i=1,ipol),
     +         (nav(i),i=1,ipol),
     +         (nfv(i),i=1,ipol),
     +         nroot
        endif
        do ir = 1, nroot
          do ip = 1, ipol
            vlo(1) = 1
            vhi(1) = naoc(ip)*nav(ip)
            do ipm1 = 1, ip - 1
              vlo(1) = vlo(1) + naoc(ipm1)*nav(ipm1)
              vhi(1) = vhi(1) + naoc(ipm1)*nav(ipm1)
            enddo
            vlo(2) = ir
            vhi(2) = ir
            mlo(1) = 1
            mhi(1) = 1
            mlo(2) = 1
            mhi(2) = naoc(ip)
            mlo(3) = 1
            mhi(3) = nav(ip)
            call nga_copy_patch('n',g_x,vlo,vhi,g_mat_mo(ip),mlo,mhi)
          enddo
cDEBUG
c         write(6,*)'*** tddft_grad_solve_save_x: A'
cDEBUG
          call tddft_grad_trans_mo2ao(ipol,nao,nfc,naoc,nocc,nav,
     +                                nfv,1,1.0d0,0.0d0,"ib",g_mo,
     +                                g_mat_mo,"ib",g_mat_ao)
c
c         Multiply from left and right by the overlap matrix
c
          if (.not.ga_create(mt_dbl,nao,nao,'AO overlap',-1,-1,g_ovlp))
     +       call errquit(
     +       'tddft_grad_solve_save_x: failed to create g_ovlp',0,
     +       GA_ERR)
          call ga_zero(g_ovlp)
          call int_1e_ga(basis,basis,g_ovlp,'overlap',oskel)
          if (oskel) call sym_symmetrize(geom,basis,.false.,g_ovlp)
          if (.not.ga_create(mt_dbl,nao,nao,'work',-1,-1,g_work))
     +       call errquit(pname//'failed to create g_work',0,GA_ERR)
          do ip = 1, ipol
            vlo(1) = 1
            vhi(1) = nao
            vlo(2) = 1
            vhi(2) = nao
            mlo(1) = ip
            mhi(1) = ip
            mlo(2) = 1
            mhi(2) = nao
            mlo(3) = 1
            mhi(3) = nao
cDEBUG
c         write(6,*)'*** tddft_grad_solve_save_x: B'
cDEBUG
            call nga_matmul_patch('N','N',1.0d0,0.0d0,g_ovlp,vlo,vhi,
     +           g_mat_ao,mlo,mhi,g_work,vlo,vhi)
cDEBUG
c         write(6,*)'*** tddft_grad_solve_save_x: C'
cDEBUG
            call nga_matmul_patch('N','N',1.0d0,0.0d0,g_work,vlo,vhi,
     +           g_ovlp,vlo,vhi,g_mat_ao,mlo,mhi)
          enddo
          if (.not.ga_destroy(g_work)) 
     +      call errquit(pname//'failed to destroy g_work',0,GA_ERR)
          if (.not.ga_destroy(g_ovlp)) 
     +      call errquit(pname//'failed to destroy g_ovlp',0,GA_ERR)
c
c         Overlap matrix done
c
          if (ga_nodeid().eq.0) then
            do ip = 1, ipol
              do j = 1, nao
                do i = 1, nao, maxi
                  mlo(1) = ip
                  mhi(1) = ip
                  mlo(2) = i
                  mhi(2) = min(nao,i-1+maxi)
                  mlo(3) = j
                  mhi(3) = j
                  idim(1) = 1
                  idim(2) = maxi
                  call nga_get(g_mat_ao,mlo,mhi,buf,idim)
                  do k = mlo(2), mhi(2)
                    write(69,err=90010)k,j,buf(k-mlo(2)+1)
                  enddo
                enddo
              enddo
            enddo
          endif
        enddo
        if (ga_nodeid().eq.0) then
          close(69)
        endif
        do ip = 1, ipol
          if (.not.ga_destroy(g_mat_mo(ip))) 
     +      call errquit(pname//'failed to destroy g_mat_mo',0,GA_ERR)
        enddo
        if (.not.ga_destroy(g_mat_ao)) 
     +      call errquit(pname//'failed to destroy g_mat_ao',0,GA_ERR)
      else
        call errquit(pname//"illegal value of x_when",0,CALC_ERR)
      endif
cDEBUG
c         write(6,*)'*** tddft_grad_solve_save_x: Z'
cDEBUG
      return
c
c     Error handling
c
90000 continue
      msg = 'tddft_grad_solve_save_x: error opening file: '//fn_x
      call errquit(pname//msg,0,DISK_ERR)
      return
90010 continue
      msg = 'tddft_grad_solve_save_x: write error on: '//fn_x
      call errquit(pname//msg,0,DISK_ERR)
      return
c
      end

c $Id$
