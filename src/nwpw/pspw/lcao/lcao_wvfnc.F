*
* $Id: lcao_wvfnc.F,v 1.1 2001-11-30 01:55:13 bylaska Exp $
*
      logical function lcao_wvfnc(oprint,it_max,
     >                      wavefunction_filename,
     >             			version,
     >						ngrid,
     >						unita,
     >						ispin,
     >						ne)

      implicit none
      logical oprint
      integer it_max
      character*30 wavefunction_filename
      integer      version
      integer      ngrid(3)
      real*8       unita(3,3)
      integer      ispin,ne(2)

#include "mafdecls.fh"
 
*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      logical value
      integer ierr

*     **** initialize atomic orbitals ****
      call aorbs_init()
      call aorbs_readall()

*     ***** allocate memory from heap memory ****
      value = MA_alloc_get(mt_dbl,(2*n2ft3d),'dn',dn(2),dn(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(nbasis2),'smatrix',
     >                     smatrix(2),smatrix(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(nbasis2),'kematrix',
     >                     kematrix(2),kematrix(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(nbasis2),'vlmatrix',
     >                     vlmatrix(2),vlmatrix(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(nbasis2),'vnlmatrix',
     >                     vnlmatrix(2),vnlmatrix(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(2*nbasis2),'vhxcmatrix',
     >                     vhxcmatrix(2),vhxcmatrix(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(2*nbasis2),'hmatrix',
     >                     hmatrix(2),hmatrix(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(2*nbasis2),'psimatrix',
     >                     psimatrix(2),psimatrix(1))
      if (.not. value) then
        go to 100 
      end if



*     **** initialize lcao matrices ****
      call lcao_generate_smatrix(dbl_mb(smatrix(1)))
      call lcao_generate_kematrix(dbl_mb(kematrix(1)))
      call lcao_generate_vlmatrix(dbl_mb(vlmatrix(1)))
      call lcao_generate_vnlmatrix(dbl_mb(vnlmatrix(1)))
      call lcao_init_dn(ispin,dbl_mb(rho(1)))

      do it=1,it_max
         call lcao_generate_vhxcmatrix(ispin,
     >           dbl_mb(rho(1)),
     >           dbl_mb(vhxcmatrix(1)))
         call lcao_generate_hmatrix(ispin,
     >                              dbl_mb(vhxcmatrix(1)),
     >                              dbl_mb(kematrix(1)),
     >                              dbl_mb(vlmatrix(1)),
     >                              dbl_mb(vnlmatrix(1)),
     >                              dbl_mb(hmatrix(1)))

          call lcao_solve(ispin,
     >                    dbl_mb(hmatrix(1)),
     >                    dbl_mb(smatrix(1)),
     >                    dbl_mb(psimatrix(1)),
     >                    ierr)
          if (ierr.ne.0) then
             value = .false.
             go to 100
          end if
          call D3dB_r_copy(1,dbl_mb(rho(1)),dbl_mb(rho_old(1)))
          call lcao_new_dn(ispin,
     >                    dbl_mb(psimatrix(1)),
     >                    dbl_mb(rho(1)))
          call lcao_mix_dn(ispin,
     >                    dbl_mb(rho_old(1)),
     >                    dbl_mb(rho(1)))
        end do
      end do

      call lcao_write_psi(dbl_mb(psimatrix(1)),
     >                    wavefunction_filename,
                          version,
                          ngrid(3),
                          unita(3,3),
                          ispin,ne(2))

*     ***** free heap memory ****
  100 continue
      value2 = MA_free_heap(rho(2))
      value2 = MA_free_heap(smatrix(2))
      value2 = MA_free_heap(kematrix(2))
      value2 = MA_free_heap(vlmatrix(2))
      value2 = MA_free_heap(vnlmatrix(2))
      value2 = MA_free_heap(vhxcmatrix(2))
      value2 = MA_free_heap(hmatrix(2))
      value2 = MA_free_heap(psimatrix(2))
      call aorbs_end()

      lcao_wvfnc = value
      return
      end



