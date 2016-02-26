c
c     Loads the solution of the linear systems of equations
c     HvD 11/2007
c
      logical function tddft_grad_solve_load_x(g_x,ipol,nao,nfc,naoc,
     +                 nocc,nav,nfv,nroot,g_mo,fn_x)
c
c     The purpose of this is that it may be beneficial
c     to use this solution as a starting point for linear system solves
c     in a geometry optimisation or a numerical hessian evaluation.
c
c     Note that on disk the solutions are stored in AO basis. The 
c     reason is that in geometry optimisations the changes to the 
c     geometry cause the orbitals to change so much that the MO basis
c     data becomes useless (i.e. it would be better to start from
c     scratch then to try and reuse it). Storing the data in AO basis
c     makes the data much less susceptible to changes in the MOs. Hence
c     the recovered data leads to genuine performance benefits.
c
c     If the data is loaded successfully the function returns .true.
c     otherwise .false. is returned.
c
c     Finally, this routine is a companion of tddft_grad_solve_save_x
c     which writes the data produced by this routine. So any changes
c     need to be consistent across both routines.
c
      implicit none
c
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
c
c     Input:
c
      integer g_x          ! 2D global array for linear system solutions
      integer ipol         ! the number of spin components
      integer nao          ! the number of AO's
      integer nfc(ipol)    ! the number of frozen core orbitals
      integer naoc(ipol)   ! the number of active occupied orbitals
      integer nocc(ipol)   ! the number of occupied orbitals
      integer nav(ipol)    ! the number of virtual orbitals
      integer nfv(ipol)    ! the number of active frozen virtual
                           ! orbitals
      integer nroot        ! the number of columns in g_x
      integer g_mo(ipol)   ! the molecular orbitals
c
      character*(*) fn_x   ! filename for solutions file
c
c     Local:
c
      integer ir                       ! counter over roots
      integer ip, ipm1                 ! counter over spins
      integer maxi                     ! the leading dimension of the 
                                       ! buffer
      parameter(maxi=100)
      double precision buf(maxi) ! a local buffer for the data
      integer i, j, k, ik, ij
      integer vlo(3), vhi(3)     ! patch limits
      integer mlo(3), mhi(3)     ! patch limits
c
      integer idim(3) ! dimensions of temporary global arrays
      integer ichk(3) ! chunking for matrix distribution
      integer g_mat_mo(ipol) ! temporary global array in MO basis
      integer g_mat_ao       ! temporary global array in AO basis
c
      integer iostatus            ! IOSTAT variable
      logical success             ! the overall success of the operation
      logical ioerror             ! has an ioerror occurred?
      logical x_exists            ! .true. if the file exists
      integer itpol, itnao        ! temporaries for checking
      integer itnfc(2), itnaoc(2) ! temporaries for checking
      integer itnav(2), itnfv(2)  ! temporaries for checking
      integer itnroot             ! temporaries for checking
c
      character*32 pname
c
      pname="tddft_grad_solve_load_x: "
c
c     Code:
c
      ioerror = .false.
      success = .false.
      call ga_sync
      do ip = 1, ipol
        idim(1) = 1
        idim(2) = naoc(ip)
        idim(3) = nav(ip)
        ichk(1) = 1
        ichk(2) = -1
        ichk(3) = -1
        if (.not.nga_create(mt_dbl,3,idim,'g_mat_mo',ichk,
     +      g_mat_mo(ip))) 
     +   call errquit(pname//'failed to create g_mat_mo',0,GA_ERR)
      enddo
      idim(1) = ipol
      idim(2) = nao
      idim(3) = nao
      ichk(1) = ipol
      if (.not.nga_create(mt_dbl,3,idim,'g_mat_ao',ichk,
     +    g_mat_ao)) 
     +   call errquit(pname//'failed to create g_mat_ao',0,GA_ERR)
c
      if (ga_nodeid().eq.0) then
        inquire(file=fn_x,exist=x_exists)
      endif
      call ga_brdcst(1001,x_exists,ma_sizeof(MT_LOG,1,MT_BYTE),0)
      if (x_exists) then
        if (ga_nodeid().eq.0) then
          open(unit=69,file=fn_x,form='unformatted',iostat=iostatus,
     +         status='unknown')
          ioerror = ioerror .or. iostatus.gt.0
          rewind(69)
          read(69,iostat=iostatus)itpol,itnao,
     +        (itnfc(i),i=1,ipol),
     +        (itnaoc(i),i=1,ipol),
     +        (itnav(i),i=1,ipol),
     +        (itnfv(i),i=1,ipol),
     +        itnroot
          ioerror = ioerror .or. iostatus.gt.0
        endif
        call ga_brdcst(1002,itpol,ma_sizeof(MT_INT,1,MT_BYTE),0)
        call ga_brdcst(1003,itnao,ma_sizeof(MT_INT,1,MT_BYTE),0)
        call ga_brdcst(1004,itnfc,ma_sizeof(MT_INT,ipol,MT_BYTE),0)
        call ga_brdcst(1005,itnaoc,ma_sizeof(MT_INT,ipol,MT_BYTE),0)
        call ga_brdcst(1006,itnav,ma_sizeof(MT_INT,ipol,MT_BYTE),0)
        call ga_brdcst(1007,itnfv,ma_sizeof(MT_INT,ipol,MT_BYTE),0)
        call ga_brdcst(1008,itnroot,ma_sizeof(MT_INT,1,MT_BYTE),0)
        if (itpol.eq.ipol.and.itnao.eq.nao.and.itnroot.eq.nroot.and.
     +      itnfc(1).eq.nfc(1).and.itnaoc(1).eq.naoc(1).and.
     +      itnav(1).eq.nav(1).and.itnfv(1).eq.nfv(1).and.
     +      (ipol.eq.1.or.(
     +       itnfc(2).eq.nfc(2).and.itnaoc(2).eq.naoc(2).and.
     +       itnav(2).eq.nav(2).and.itnfv(2).eq.nfv(2))))
     +  then
          do ir = 1, nroot
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
                  if (ga_nodeid().eq.0) then
                    do k = mlo(2), mhi(2)
                      read(69,iostat=iostatus)ik,ij,buf(k-mlo(2)+1)
                      ioerror = ioerror .or. iostatus.gt.0
                      ioerror = ioerror .or. (ik.ne.k.or.ij.ne.j)
                    enddo
                    call nga_put(g_mat_ao,mlo,mhi,buf,idim)
                  endif
                enddo
              enddo
            enddo
            call tddft_grad_trans_ao2mo(ipol,nao,nfc,naoc,nocc,nav,
     +                                  nfv,1,1.0d0,0.0d0,"ib",g_mo,
     +                                  g_mat_ao,g_mat_mo,"ib")
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
              call nga_copy_patch('n',g_mat_mo(ip),mlo,mhi,
     +                            g_x,vlo,vhi)
            enddo
          enddo
        endif
        if (ga_nodeid().eq.0) then
          close(69,iostat=iostatus)
          ioerror = ioerror .or. iostatus.gt.1
        endif
        success = .not.ioerror
      endif
      do ip = 1, ipol
        if (.not.ga_destroy(g_mat_mo(ip))) 
     +    call errquit(pname//'failed to destroy g_mat_mo',0,GA_ERR)
      enddo
      if (.not.ga_destroy(g_mat_ao)) 
     +    call errquit(pname//'failed to destroy g_mat_ao',0,GA_ERR)
c
      call ga_brdcst(maxi,success,ma_sizeof(MT_LOG,1,MT_BYTE),0)
c
      tddft_grad_solve_load_x = success
c
      return
      end
c $Id$
