      subroutine bse_dipole(pars)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "bse.fh"
#include "cdft.fh"
      type(bse_params_t) :: pars

      character(*),parameter :: pname = 'bse_dipole: '

      integer ixyz,isp,ilo,ihi,jlo,jhi,amo,imo
      integer g_dip(3),handle(3),g_tmp1,g_tmp2
      integer k,klocal,kproc,kglobal,tovern,modtn,ktmp,ld
      integer ga_create_atom_blocked
      external ga_create_atom_blocked

      ! Create dipole matrices
      do ixyz=1,3
        g_dip(ixyz) = ga_create_atom_blocked(geom,ao_bas_han,'Dipoles')
        if(.not.nga_create(mt_dbl,1,pars%npoles(1),'dip a',0,
     &                     pars%g_dipa(ixyz)))
     &   call errquit(pname//'could not create g_dipa',0,GA_ERR)
        if (pars%ipol.gt.1) then
          if(.not.nga_create(mt_dbl,1,pars%npoles(2),'dip b',0,
     &                       pars%g_dipb(ixyz)))
     &     call errquit(pname//'could not create g_dipb',0,GA_ERR)
        endif
      enddo

      ! Compute dipole AO integrals for the 3 components
      call int_dip_ga(ao_bas_han, ao_bas_han, g_dip(1), g_dip(2),
     &                g_dip(3))

      ! Transform dipole matrices
      do isp=1,pars%ipol
        if (.not.ga_create(mt_dbl,nbf_ao,pars%nocc(isp),'tmp1',
     &                     0,pars%nocc(isp),g_tmp1))
     &    call errquit(pname//'could not create g_tmp1',0,GA_ERR)
        if (.not.ga_create(mt_dbl,pars%nvir(isp),pars%nocc(isp),'tmp2',
     &                     pars%nvir(isp),0,g_tmp2))
     &    call errquit(pname//'could not create g_tmp2',0,GA_ERR)

        if (isp.eq.1) then
          handle(:) = pars%g_dipa(:)
        else
          handle(:) = pars%g_dipb(:)
        endif

        tovern = pars%npoles(isp)/pars%nprocs
        modtn  = mod(pars%npoles(isp),pars%nprocs)
        call ga_distribution(g_tmp2,pars%me,ilo,ihi,jlo,jhi)
        do ixyz=1,3
          call ga_matmul_patch('n','n',1d0,0d0,
     &                    g_dip(ixyz),1,nbf_ao,1,nbf_ao,
     &                    pars%g_movecs(isp),1,nbf_ao,1,pars%nocc(isp),
     &                    g_tmp1,1,nbf_ao,1,pars%nocc(isp))
          call ga_matmul_patch('t','n',1d0,0d0,
     &                    pars%g_movecs(isp),1,pars%nvir(isp),1,nbf_ao,
     &                    g_tmp1,1,nbf_ao,1,pars%nocc(isp),
     &                    g_tmp2,1,pars%nvir(isp),1,pars%nocc(isp))
          call nga_copy_patch('n',g_tmp2,(/1,1/),
     &         (/pars%nvir(isp),pars%nocc(isp)/),handle(ixyz),1,
     &         pars%npoles(isp))
        enddo
        if ((.not.ga_destroy(g_tmp1)).and.(.not.ga_destroy(g_tmp2)))
     &    call errquit(pname//'could not destroy GAs',0,GA_ERR)   
      enddo

      do ixyz=1,3
        if (.not.ga_destroy(g_dip(ixyz)))
     &      call errquit(pname//'could not destroy GAs',0,GA_ERR)   
      enddo

      end
