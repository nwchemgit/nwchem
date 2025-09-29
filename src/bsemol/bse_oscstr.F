      subroutine bse_oscstr(pars,omega,totpoles)
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "bse.fh"
#include "cdft.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "util_parms.fh"
#include "inp.fh"
      type(bse_params_t) :: pars

      character(*), parameter :: pname = 'bse_oscstr: '
      character(len=7) :: key

      integer totpoles,npoles
      double precision omega(*)
      double precision oscstr(3),oscstr_tot,tmp,moms(20),maxr,r(1)

      integer llpole,ulpole,iroot,ixyz,ipole,mag_transdip(3)

      integer isym,jsym,i,j,k,l,pol,occ,vir,ioff
      integer g_td(2),g_xy(2),g_dipmag, g_ovlp, g_work, g_corr
      integer g_tdtot(2)
      integer k_corr, l_corr, ulloop
      character*4 oname,vname

      integer g_dipx, g_dipy, g_dipz
      integer g_qdrxx, g_qdrxy, g_qdrxz, g_qdryy, g_qdryz, g_qdrzz
      integer g_work2, g_work3, g_x(2), g_y(2)

      logical ltransden,stat
      double precision origin(3)
      double precision na,nb,s2,s2_full
      double precision espread
      data origin/0d0,0d0,0d0/

      character*(nw_max_path_len) fn_civecs, fn_transden
      logical,external :: dmat_to_file
      integer,external :: ga_create_atom_blocked

      character*5 spin(2)
      double precision x(totpoles),y(totpoles)
      double precision thresh  ! Threshold for printing
      parameter (thresh=5.0d-2)
      double precision xe(3), xh(3), xe2_exp, xh2_exp, xe_exp_2
      double precision xh_exp_2, xhxe_exp, dhtoe, dexc, sigma_h, sigma_e
      double precision omega_norm, cov, pcc, dcd
      spin(1)='alpha'
      spin(2)='beta '

      if(.not.rtdb_get(pars%rtdb,'bse:ltransden',mt_log,1,ltransden))
     $    ltransden = .false.
      if(.not.nga_create(mt_dbl,3,(/nbf_ao,nbf_ao,3/),'mag_dipole',
     $                   (/nbf_ao,-1,-1/),g_dipmag))
     $   call errquit('bse_oscstr: nga_create failed g_dipmag',0,GA_ERR)
      call ga_zero(g_dipmag)
      call int_giao_1ega(ao_bas_han,ao_bas_han,g_dipmag,'angmom',origin,
     $                   1,.false.) 
      call ga_scale(g_dipmag,0.5d0)

      g_ovlp = ga_create_atom_blocked(geom, ao_bas_han, 'Overlap')
      call ga_zero(g_ovlp)
      call int_1e_ga(ao_bas_han, ao_bas_han, g_ovlp, 'overlap', oskel)
      if(oskel) call sym_symmetrize(geom, ao_bas_han, .false., g_ovlp)

      if (pars%ipol.eq.2) then
        if(.not.ma_push_get(mt_dbl,pars%nmo*pars%nmo,'corr',
     $                      l_corr,k_corr))
     $    call errquit('bse_oscstr: failed to create MA',0,MA_ERR)
        if(.not.ga_create(mt_dbl,nbf_ao,pars%nmo,'work',0,0,g_work))
     $    call errquit('bse_oscstr: failed to create GA',0,GA_ERR)
        if(.not.ga_create(mt_dbl,pars%nmo,pars%nmo,'work',
     $                    pars%nmo,0,g_corr))
     $    call errquit('bse_oscstr: failed to create GA',0,GA_ERR)
        call ga_dgemm('n','n',nbf_ao,pars%nmo,nbf_ao,1d0,g_ovlp,
     $                 pars%g_movecs(2),0d0,g_work)
        call ga_dgemm('t','n',pars%nmo,pars%nmo,nbf_ao,1d0,
     $                 pars%g_movecs(1),g_work,0d0,g_corr)
        call ga_get(g_corr,1,pars%nmo,1,pars%nmo,dbl_mb(k_corr),
     $              pars%nmo)
        stat = ga_destroy(g_work).and.ga_destroy(g_corr)
        if(.not.stat)
     $    call errquit('bse_osctr: failed to destroy GA',58,GA_ERR)   
      endif

      g_dipx = ga_create_atom_blocked(geom, ao_bas_han, 'X')
      g_dipy = ga_create_atom_blocked(geom, ao_bas_han, 'Y')
      g_dipz = ga_create_atom_blocked(geom, ao_bas_han, 'Z')
      g_qdrxx = ga_create_atom_blocked(geom, ao_bas_han, 'XX')
      g_qdrxy = ga_create_atom_blocked(geom, ao_bas_han, 'XY')
      g_qdrxz = ga_create_atom_blocked(geom, ao_bas_han, 'XZ')
      g_qdryy = ga_create_atom_blocked(geom, ao_bas_han, 'YY')
      g_qdryz = ga_create_atom_blocked(geom, ao_bas_han, 'YZ')
      g_qdrzz = ga_create_atom_blocked(geom, ao_bas_han, 'ZZ')
      g_work = ga_create_atom_blocked(geom, ao_bas_han, 'work1')
      g_work2 = ga_create_atom_blocked(geom, ao_bas_han, 'work2')
      g_work3 = ga_create_atom_blocked(geom, ao_bas_han, 'work3')

      call int_dip_ga(ao_bas_han, ao_bas_han, g_dipx, g_dipy, g_dipz)
      call int_qdr_ga(ao_bas_han, ao_bas_han, g_qdrxx, g_qdrxy, g_qdrxz,
     $                g_qdryy, g_qdryz, g_qdrzz)



      npoles = pars%npoles(1)
      if (pars%ipol.gt.1) npoles = npoles + pars%npoles(2)


      if (pars%singlet) then
        key = 'singlet'
      else
        key = 'triplet'
      endif

      ! Get the lower and upper limits
      llpole = -1
      if (pars%dowindow .and. (.not.pars%davidson)) then
        iroot = 0
        do ipole=1,npoles
          if (omega(ipole).lt.pars%elower) cycle
          if (llpole.eq.-1) llpole = ipole
          iroot = iroot + 1
          if (omega(ipole).gt.pars%eupper) then
            ulpole = ipole - 1
            exit
          endif
          if (iroot.eq.pars%nroots) then
            ulpole = ipole
            exit
          endif
        enddo
      else
        llpole = 1
        ulpole = pars%nroots
      endif
        

      ! Get symmetry of ground state
      isym = 1
      if (pars%ipol.eq.2) then
        do i=1,pars%ipol
          do j=1,pars%nocc(i)
            isym = ieor(isym-1,int_mb(pars%k_irs(i)+j-1)-1)+1
          enddo
        enddo
      endif
      call sym_irrepname(geom, isym, oname)


      ! Allocate transition density matrices
      do i=1,pars%ipol
        if(.not.ga_create(mt_dbl,nbf_ao,nbf_ao,'Transition density',
     $                    -1,-1,g_td(i))) call errquit
     $    ('bse_oscstr: failed to create g_td',0,GA_ERR)
        if(.not.ga_create(mt_dbl,nbf_ao,nbf_ao,'Transition density',
     $                    -1,-1,g_tdtot(i))) call errquit
     $    ('bse_oscstr: failed to create g_tdtot',0,GA_ERR)   
        if(.not.ga_create(mt_dbl,pars%npoles(i),1,'X+Y vector',-1,1,
     $                    g_xy(i))) call errquit
     $    ('bse_oscstr: failed to create g_xy',0,GA_ERR)
        if(.not.ga_create(mt_dbl,pars%npoles(i),1,'X vector',-1,1,
     $                    g_x(i))) call errquit
     $    ('bse_oscstr: failed to create g_x',0,GA_ERR)   
        if(.not.ga_create(mt_dbl,pars%npoles(i),1,'Y vector',-1,1,
     $                    g_y(i))) call errquit
     $    ('bse_oscstr: failed to create g_y',0,GA_ERR)   
      enddo

      ! Compute the oscillator strengths around an energy window
      iroot = 0
      do ipole=llpole,ulpole
        iroot      = iroot + 1
        oscstr_tot = 0d0
        moms(:)    = 0d0
        oscstr(:)  = 0d0
        mag_transdip = 0d0

        call ga_copy_patch('n',pars%g_apb,1,pars%npoles(1),ipole,ipole,
     &                         g_x(1),1,pars%npoles(1),1,1)
        if(pars%ipol.gt.1) then
          call ga_copy_patch('n',pars%g_apb,pars%npoles(1)+1,npoles,
     &                       ipole,ipole,g_x(2),1,pars%npoles(2),1,1)
        endif
        call tddft_transfm(1, g_x, pars%g_movecs, nbf_ao, pars%nocc,
     &                     (/pars%nmo, pars%nmo/), (/0,0/), (/0,0/),
     &                     pars%ipol, g_td, 1)
        call ga_transpose(g_td(1), g_tdtot(1))
        if(pars%ipol.gt.1) call ga_transpose(g_td(2), g_tdtot(2))

        if(pars%davdison) then
          call ga_copy_patch('n',pars%g_amb,1,pars%npoles(1),ipole,ipole,
     &                         g_y(1),1,pars%npoles(1),1,1)
          if(pars%ipol.gt.1) then
            call ga_copy_patch('n',pars%g_amb, pars%npoles(1)+1, npoles,
     &                          ipole, ipole, g_y(2), 1, pars%npoles(2),
     &                          1, 1)
          endif
          call tddft_transfm(1, g_y, pars%g_movecs, nbf_ao, pars%nocc,
     &                       (/pars%nmo, pars%nmo/), (/0,0/), (/0,0/),
     &                       pars%ipol, g_td, 1)
          call ga_add(1.0d0, g_td(1), 1.0d0, g_tdtot(1), g_td(1))
          if(pars%ipol.gt.1)
     &      call ga_add(1.0d0, g_td(2), 1.0d0, g_tdtot(2), g_td(2))
        else
          call ga_copy(g_tdtot(1), g_td(1))
          if (pars%ipol.gt.1) call ga_copy(g_tdtot(2), g_td(2))
        endif

        do i=1,pars%ipol
          call multipole_density(ao_bas_han,origin,3,g_td(i),moms,20)
          do ixyz=1,3
            mag_transdip(ixyz) = mag_transdip(ixyz) + 
     $        nga_ddot_patch(g_td(i),'n',(/1,1/),(/nbf_ao,nbf_ao/),
     $                       g_dipmag,'n',(/1,1,ixyz/),
     $                       (/nbf_ao,nbf_ao,ixyz/))
          enddo
        enddo

        if(pars%ipol.gt.1) call ga_add(1d0,g_td(1),1d0,g_td(2),g_td(1))

        ! Indicators
        call ga_dgemm('t','n',nbf_ao,nbf_ao,nbf_ao,1d0,g_td(1),g_ovlp,
     $                 0d0,g_work)
        call ga_dgemm('n','n',nbf_ao,nbf_ao,nbf_ao,1d0,g_work,g_td(1),
     $                 0d0,g_work2) 
        omega_norm = ga_ddot(g_work2,g_ovlp)
        xe(1) = ga_ddot(g_work2,g_dipx)/omega_norm
        xe(2) = ga_ddot(g_work2,g_dipy)/omega_norm
        xe(3) = ga_ddot(g_work2,g_dipz)/omega_norm
        xe_exp_2 = xe(1)**2 + xe(2)**2 + xe(3)**2
        xe2_exp  = (ga_ddot(g_work2,g_qdrxx) + 
     $             ga_ddot(g_work2,g_qdryy) +
     $             ga_ddot(g_work2,g_qdrzz))/omega_norm

          
        call ga_dgemm('n','n',nbf_ao,nbf_ao,nbf_ao,1d0,g_td(1),g_ovlp,
     $                 0d0,g_work)
        call ga_dgemm('n','t',nbf_ao,nbf_ao,nbf_ao,1d0,g_work,g_td(1),
     $                 0d0,g_work2) 
        xh(1) = ga_ddot(g_work2,g_dipx)/omega_norm
        xh(2) = ga_ddot(g_work2,g_dipy)/omega_norm
        xh(3) = ga_ddot(g_work2,g_dipz)/omega_norm
        xh_exp_2 = xh(1)**2 + xh(2)**2 + xh(3)**2
        xh2_exp  = (ga_ddot(g_work2,g_qdrxx) + 
     $             ga_ddot(g_work2,g_qdryy) +
     $             ga_ddot(g_work2,g_qdrzz))/omega_norm


        xhxe_exp = 0d0 
        call ga_dgemm('n','n',nbf_ao,nbf_ao,nbf_ao,1d0,g_td(1),g_dipx,
     $                 0d0,g_work)
        call ga_dgemm('n','t',nbf_ao,nbf_ao,nbf_ao,1d0,g_work,g_td(1),
     $                 0d0,g_work2)
        xhxe_exp = xhxe_exp + ga_ddot(g_work2,g_dipx)/omega_norm
        call ga_dgemm('n','n',nbf_ao,nbf_ao,nbf_ao,1d0,g_td(1),g_dipy,
     $                 0d0,g_work)
        call ga_dgemm('n','t',nbf_ao,nbf_ao,nbf_ao,1d0,g_work,g_td(1),
     $                 0d0,g_work2)
        xhxe_exp = xhxe_exp + ga_ddot(g_work2,g_dipy)/omega_norm
        call ga_dgemm('n','n',nbf_ao,nbf_ao,nbf_ao,1d0,g_td(1),g_dipz,
     $                 0d0,g_work)
        call ga_dgemm('n','t',nbf_ao,nbf_ao,nbf_ao,1d0,g_work,g_td(1),
     $                 0d0,g_work2)
        xhxe_exp = xhxe_exp + ga_ddot(g_work2,g_dipz)/omega_norm


        sigma_h = dsqrt(xh2_exp - xh_exp_2)
        sigma_e = dsqrt(xe2_exp - xe_exp_2)

        dhtoe   = dsqrt((xe(1)-xh(1))**2 + 
     $                  (xe(2)-xh(2))**2 +
     $                  (xe(3)-xh(3))**2)
        dexc    = dsqrt(xe2_exp + xh2_exp - 2d0*xhxe_exp)
        cov     = xhxe_exp - xe(1)*xh(1) - xe(2)*xh(2) - xe(3)*xh(3)
        pcc     = cov/sigma_e/sigma_h
        dcd     = dhtoe - 0.5d0*(sigma_e + sigma_h)

        if (ltransden) then
          call util_file_name('tdens',.false.,.false.,fn_transden)
          call createuniquefilename(fn_transden,ipole)
          call util_file_name_resolve(fn_transden,.false.)
          if (pars%ipol.gt.1)
     &      call ga_add(1d0,g_td(1),1d0,g_td(2),g_td(1))

          call ga_symmetrize(g_td(1))

          if(.not. dmat_to_file(g_td(1),fn_transden))
     &      call errquit('bse_oscstr: dmat_to_file error',0,CALC_ERR)
        endif

        if (pars%ipol.eq.1) then
          do i=1,20
            moms(i) = moms(i)*dsqrt(2.0d0)
          enddo
        endif

        if (pars%davidson) then
          call ga_get(pars%g_apb,1,npoles,ipole,ipole,x,npoles)
          call ga_get(pars%g_amb,1,npoles,ipole,ipole,y,npoles)
        endif

        if (pars%davidson.and.pars%ipol.eq.2) then
          call tddft_s2(pars%tda, (/0,0/), pars%nocc, 
     &                    (/pars%nmo,pars%nmo/), (/0,0/),
     &                    pars%npoles,x,x(pars%npoles(1)+1),y,
     &                    y(pars%npoles(1)+1),dbl_mb(k_corr),s2)
          call bse_s2(pars%nocc,pars%nmo,pars%npoles,x,
     &     x(pars%npoles(1)+1), y, y(pars%npoles(1)+1), dbl_mb(k_corr),
     &     s2_full) 
          if (abs(s2_full).lt.0.25d0) then
            key = 'singlet'
          elseif (abs(abs(s2_full)-2.0d0).lt.0.25) then
            key = 'triplet'
          else
            key = '       '
          endif
        endif

        if (pars%me.eq.0) then
          call tddft_oscstr(moms, omega(ipole), mag_transdip,
     &                      oscstr, oscstr_tot)

          l = 0
          maxr = 0.0d0
          pol = 1
          occ = 1
          vir = pars%nocc(1) + 1
          do i=1,pars%ipol
            if (pars%npoles(i).eq.0) cycle
            do j=1,pars%nocc(i)
              do k=pars%nocc(i)+1,pars%nmo
                l = l + 1
                call ga_get(pars%g_apb,l,l,ipole,ipole,r,1)
                if (dabs(r(1)).gt.maxr) then
                  maxr = dabs(r(1))
                  pol = i
                  occ = j
                  vir = k
                endif
              enddo
            enddo
          enddo
          jsym = ieor(isym-1,int_mb(pars%k_irs(pol)+occ-1)-1)+1
          jsym = ieor(jsym-1,int_mb(pars%k_irs(pol)+vir-1)-1)+1
          call sym_irrepname(geom, jsym, vname)

          write(luout,9110) iroot,key,vname,omega(ipole),
     &                        omega(ipole)*ha2ev

          if (pars%davidson.and.pars%ipol.eq.2) then
            write(luout,*)
            write(luout,9180) s2
            write(luout,9190) s2_full
          endif

          if (pars%singlet) then
            write(luout,9200) moms(2),moms(3),moms(4),
     &                        moms(5),moms(6),moms(7),
     &                        moms(8),moms(9),moms(10),
     &                        oscstr(1),oscstr(2),oscstr(3),
     &                        oscstr_tot
          else
            write(luout,9210)
          endif

          write(luout,*)
          write(luout,9340) omega_norm
          write(luout,9300) sigma_e*cau2ang
          write(luout,9310) sigma_h*cau2ang
          write(luout,9320) dhtoe*cau2ang
          write(luout,9330) dexc*cau2ang
          write(luout,9360) dcd*cau2ang
          write(luout,9350) pcc


          if (pars%davidson) then

          write(luout,*)
          if (ipol.eq.1) then
            l=0
            do j=1,pars%nocc(1)
              do k=pars%nocc(1)+1,pars%nmo
                l=l+1
                if (pars%tda) then
                  r=x(l)
                  call sym_irrepname
     1              (geom,int_mb(pars%k_irs(1)+j-1),oname)
                  call sym_irrepname
     1              (geom,int_mb(pars%k_irs(1)+k-1),vname)
                  if ((dabs(r(1)).gt.thresh))
     2              write(LuOut,9120) j,oname,k,vname,r
                else
                  r=x(l)
                  call sym_irrepname
     1              (geom,int_mb(pars%k_irs(1)+j-1),oname)
                  call sym_irrepname
     1              (geom,int_mb(pars%k_irs(1)+k-1),vname)
                  if ((dabs(r(1)).gt.thresh))
     2              write(LuOut,9130) j,oname,k,vname,r
                  r=y(l)
                  if ((dabs(r(1)).gt.thresh))
     2              write(LuOut,9140) j,oname,k,vname,r
                endif
              enddo
            enddo
          else
            do i=1,ipol
              ioff = (i-1)*pars%npoles(1)
              l=0
              do j=1,pars%nocc(i)
                do k=pars%nocc(i)+1,pars%nmo
                  l=l+1
                  if (pars%tda) then
                    r = x(ioff+l)
                    call sym_irrepname
     1                (geom,int_mb(pars%k_irs(i)+j-1),oname)
                    call sym_irrepname
     1                (geom,int_mb(pars%k_irs(i)+k-1),vname)
                    if ((dabs(r(1)).gt.thresh))
     2                write(LuOut,9150) j,spin(i),oname,
     3                k,spin(i),vname,r
                  else
                    r = x(ioff+l)
                    call sym_irrepname
     1                (geom,int_mb(pars%k_irs(i)+j-1),oname)
                    call sym_irrepname
     1                (geom,int_mb(pars%k_irs(i)+k-1),vname)
                    if ((dabs(r(1)).gt.thresh))
     2                write(LuOut,9160) j,spin(i),oname,
     3                k,spin(i),vname,r
                    r=y(ioff+l)
                    if ((dabs(r(1)).gt.thresh))
     2                write(LuOut,9170) j,spin(i),oname,
     3                k,spin(i),vname,r
                  endif
                enddo
              enddo
            enddo      
          endif
          call util_flush(LuOut)

          endif

        endif

        call ga_sync()

      enddo

      if (.not.ga_destroy(g_dipmag))
     &  call errquit('bse_oscstr: ga_destroy failed',0,GA_ERR)
      do i=1,pars%ipol
        if (.not.ga_destroy(g_td(i)))
     &  call errquit('bse_oscstr: ga_destroy failed',0,GA_ERR)
        if (.not.ga_destroy(g_xy(i)))
     &  call errquit('bse_oscstr: ga_destroy failed',0,GA_ERR)
        if (.not.ga_destroy(g_tdtot(i)))
     &  call errquit('bse_oscstr: ga_destroy failed',0,GA_ERR)
        if (.not.ga_destroy(g_x(i)))
     &  call errquit('bse_oscstr: ga_destroy failed',0,GA_ERR)
        if (.not.ga_destroy(g_y(i)))
     &  call errquit('bse_oscstr: ga_destroy failed',0,GA_ERR)
      enddo
      if(.not.(ga_destroy(g_dipx).and.ga_destroy(g_dipy).and.
     $         ga_destroy(g_dipz).and.ga_destroy(g_qdrxx).and. 
     $         ga_destroy(g_qdrxy).and.ga_destroy(g_qdrxz).and. 
     $         ga_destroy(g_qdryy).and.ga_destroy(g_qdryz).and. 
     $         ga_destroy(g_qdrzz).and.ga_destroy(g_ovlp).and.
     $         ga_destroy(g_work).and.ga_destroy(g_work2).and.
     $         ga_destroy(g_work3)))
     &  call errquit('bse_oscstr: ga_destroy failed',0,GA_ERR)


 9110 format(2x,
     1  '---------------------------------------------------------------
     3-------------',
     2  /,2x,'Root',i4,1x,a7,1x,a4,f22.9,' a.u.',f22.4,' eV ',/,2x,
     1  '---------------------------------------------------------------
     3-------------')
 9200 format(5x,'Transition Moments    X',f9.5,'   Y',f9.5,'   Z',f9.5
     1    ,/,5x,'Transition Moments   XX',f9.5,'  XY',f9.5,'  XZ',f9.5
     2    ,/,5x,'Transition Moments   YY',f9.5,'  YZ',f9.5,'  ZZ',f9.5
     5    ,/,5x,'Dipole Oscillator Strength',18x,f14.10
     5    ,/,5x,'Electric Quadrupole       ',18x,f14.10
     5    ,/,5x,'Magnetic Dipole           ',18x,f14.10
     5    ,/,5x,'Total Oscillator Strength ',18x,f14.10)
 9210 format(5x,'Transition Moments                    Spin forbidden'
     1    ,/,5x,'Oscillator Strength                   Spin forbidden')
 9120 format(5x,'Occ.',i5,2x,a4,'---  Virt.',i5,2x,a4,f12.5)
 9130 format(5x,'Occ.',i5,2x,a4,'---  Virt.',i5,2x,a4,f8.5,' X')
 9140 format(5x,'Occ.',i5,2x,a4,'---  Virt.',i5,2x,a4,f8.5,' Y')
 9150 format(5x,'Occ.',i5,1x,a5,1x,a4,'---  Virt.',
     1  i5,1x,a5,1x,a4,f12.5)
 9160 format(5x,'Occ.',i5,1x,a5,1x,a4,'---  Virt.',
     1  i5,1x,a5,1x,a4,f8.5,' X')
 9170 format(5x,'Occ.',i5,1x,a5,1x,a4,'---  Virt.',
     1  i5,1x,a5,1x,a4,f8.5,' Y')
 9180 format(2x,'    <S2> = ',f8.4)
 9190 format(2x,'   D<S2> = ',f8.4)
 9300 format(5x,'     sigma_e = ',f8.4, ' Angstrom')
 9310 format(5x,'     sigma_h = ',f8.4, ' Angstrom')
 9320 format(5x,'    |d_h->e| = ',f8.4, ' Angstrom')
 9330 format(5x,'       d_exc = ',f8.4, ' Angstrom')
 9360 format(5x,'        d_CD = ',f8.4, ' Angstrom')
 9350 format(5x,'Pearson C.C. = ',f8.4)
 9340 format(5x,'Exciton norm = ',f8.4)



      end
