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
#include "inp.fh"
      type(bse_params_t) :: pars

      character(*), parameter :: pname = 'bse_oscstr: '
      character(len=7) :: key

      integer totpoles
      double precision omega(*)
      double precision oscstr(3),oscstr_tot,tmp,moms(20),maxr,r(1)

      integer llpole,ulpole,iroot,ixyz,ipole,mag_transdip(3)

      integer isym,jsym,i,j,k,l,pol,occ,vir
      integer g_td(2),g_xy(2),g_dipmag
      character*4 oname,vname

      logical ltransden
      double precision origin(3)
      data origin/0d0,0d0,0d0/

      character*(nw_max_path_len) fn_civecs, fn_transden
      logical,external :: dmat_to_file

      if(.not.rtdb_get(pars%rtdb,'bse:ltransden',mt_log,1,ltransden))
     $    ltransden = .false.
      if(.not.nga_create(mt_dbl,3,(/nbf_ao,nbf_ao,3/),'mag_dipole',
     $                   (/nbf_ao,-1,-1/),g_dipmag))
     $   call errquit('bse_oscstr: nga_create failed g_dipmag',0,GA_ERR)
      call ga_zero(g_dipmag)
      call int_giao_1ega(ao_bas_han,ao_bas_han,g_dipmag,'angmom',origin,
     $                   1,.false.) 
      call ga_scale(g_dipmag,0.5d0)


      if (pars%singlet) then
        key = 'singlet'
      else
        key = 'triplet'
      endif

      ! Get the lower and upper limits
      llpole = -1
      if (pars%dowindow) then
        iroot = 0
        do ipole=1,totpoles
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
        if(.not.ga_create(mt_dbl,pars%npoles(i),1,'X+Y vector',-1,1,
     $                    g_xy(i))) call errquit
     $    ('bse_oscstr: failed to create g_xy',0,GA_ERR)
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
     &                         g_xy(1),1,pars%npoles(1),1,1)
        if (pars%ipol.gt.1) then
          call ga_copy_patch('n',pars%g_apb,pars%npoles(1)+1,totpoles,
     &                        ipole,ipole,g_xy(2),1,pars%npoles(2),1,1)
        endif
        call tddft_transfm(1,g_xy,pars%g_movecs,nbf_ao,pars%nocc,
     &                    (/pars%nmo,pars%nmo/),0,0,pars%ipol,g_td,1)

        do i=1,pars%ipol
          call multipole_density(ao_bas_han,origin,3,g_td(i),moms,20)
          do ixyz=1,3
            mag_transdip(ixyz) = mag_transdip(ixyz) + 
     $        nga_ddot_patch(g_td(i),'n',(/1,1/),(/nbf_ao,nbf_ao/),
     $                       g_dipmag,'n',(/1,1,ixyz/),
     $                       (/nbf_ao,nbf_ao,ixyz/))
          enddo
        enddo

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

        if (pars%me.eq.0) then
          call tddft_oscstr(moms, omega(ipole), mag_transdip,
     &                      oscstr, oscstr_tot)

c        do ixyz=1,3
c          moms(ixyz) = moms(ixyz) + nga_ddot_patch(pars%g_apb,'n',
c     &          (/1,ipole/),(/pars%npoles(1),ipole/),
c     &          pars%g_dipa(ixyz),'n',1,pars%npoles(1))
c          oscstr = oscstr + moms(ixyz)**2
c        enddo
c        if (pars%ipol.gt.1) then
c          do ixyz=1,3
c            moms(ixyz) = moms(ixyz) + 
c     &        nga_ddot_patch(pars%g_apb,'n',(/pars%npoles(1)+1,ipole/),
c     &            (/totpoles,ipole/),
c     &            pars%g_dipb(ixyz),'n',1,pars%npoles(2))
c            oscstr = oscstr + moms(ixyz)**2
c          enddo
c        endif
c        oscstr = 2d0/3d0*omega(ipole)*oscstr

          l = 0
          maxr = 0.0d0
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

          if (pars%singlet) then
            write(luout,9200) moms(2),moms(3),moms(4),
     &                        moms(5),moms(6),moms(7),
     &                        moms(8),moms(9),moms(10),
     &                        oscstr(1),oscstr(2),oscstr(3),
     &                        oscstr_tot
          else
            write(luout,9210)
          endif


        endif

      enddo

      if (.not.ga_destroy(g_dipmag))
     &  call errquit('bse_oscstr: ga_destroy failed',0,GA_ERR)
      do i=1,pars%ipol
        if (.not.ga_destroy(g_td(i)))
     &  call errquit('bse_oscstr: ga_destroy failed',0,GA_ERR)
        if (.not.ga_destroy(g_xy(i)))
     &  call errquit('bse_oscstr: ga_destroy failed',0,GA_ERR)
      enddo

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


      end
