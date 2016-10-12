c
c     Calculate or read the density matrix
c
      subroutine esp_denmat(occ,ene)
c
c $Id$
c
      implicit none
c
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
#include "esp_params.fh"
#include "esp_common.fh"
#include "util.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "inp.fh"
c  
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
      logical movecs_read,movecs_read_header
      external movecs_read,movecs_read_header
c
      real*8 ene(2*nbf),occ(2*nbf), scale
c
      integer ilo,ihi,jlo,jhi,lthi,lthj,i
      integer lg_s,k_scr,l_scr,icol,iset, lg_p(2)
      integer ioff
c
      character*20 scftyp
      character*255 dmat_file
c
      logical ao_1prdm_read
      external ao_1prdm_read
c
c     allocate global arrays
c
      lg_d    =ga_create_atom_blocked(igeom,ibasis,'density')
      lg_v    =ga_create_atom_blocked(igeom,ibasis,'vector')
      lg_s    =ga_create_atom_blocked(igeom,ibasis,'scr')
      lg_p(1)   =ga_create_atom_blocked(igeom,ibasis,'density')
c
      call ga_zero(lg_d)
      call ga_zero(lg_v)
      call ga_zero(lg_s)
      call ga_zero(lg_p(1))
c
      if (denmat.ne.0) then
c
c       read density matrix
c
        call util_file_name('dmat',.false.,.false.,dmat_file) !
        if (.not. ao_1prdm_read(nbf,lg_d,dmat_file))
     &     call errquit('esp_denmat: ao_1prdm_read failed',0,0)
c
      else ! denmat
c
c       calculate density matrix from movecs
c
        if (movfil.eq.' ') then
           call util_file_name('movecs',.false.,.false.,movfil)
        else
           call util_file_name_resolve(movfil,.false.)
        end if
        
        if(ga_nodeid().eq.0) write(luout,*) ' using M.O. file = ',
     C       movfil(1:inp_strlen(movfil))

        if(.not.movecs_read_header(movfil(1:inp_strlen(movfil)),
     R       title,basnam,scftyp,nbf,nset,
     +       nmo,2)) call md_abort('esp: movecs_read_header failed',me)
c
c       if nset=2, put the first set's occupation numbers and eigenvalues
c       in 1:nbf, put the second set's in nbf+1:2*nbf, and allocate
c       beta density matrix.
c
        if (nset.gt.1) then
          ioff = nbf
          lg_p(2)   =ga_create_atom_blocked(igeom,ibasis,'beta density')
          call ga_zero(lg_p(2))
        else
          ioff = 0
        endif
c
        nocc=0
        charge=chgnuc
c
c     read orbital occupations, energies and vectors from movecs
c
      do 1 iset=1,nset
c
      if (.not. movecs_read(movfil(1:inp_strlen(movfil)), 
     I        iset, dbl_mb(ioptr+ioff),
     &     dbl_mb(ieptr+ioff), lg_v))
     & call md_abort('esp_denmat:movecs_read failed ',0)
c
      do 2 i=1,nbf
      if(dbl_mb(ioptr+i-1+ioff).gt.small) then
         charge=charge-dbl_mb(ioptr+i-1+ioff)
      endif
      occ(i+ioff) = dbl_mb(ioptr+i-1+ioff)
      ene(i+ioff) = dbl_mb(ieptr+i-1+ioff)
    2 continue
c
      call ga_copy(lg_v,lg_s)
      call ga_sync()
c
      call ga_distribution(lg_s,me,ilo,ihi,jlo,jhi)
      lthi=ihi-ilo+1
      lthj=jhi-jlo+1
c
      if(lthi.ne.0.and.lthj.ne.0) then
      if(lthi.lt.0.or.lthj.lt.0) 
     + call md_abort('esp_denmat: distribution error',me)
c
      if(.not.ma_push_get(mt_dbl,lthi,'vector',l_scr,k_scr))
     + call md_abort('esp_denmat: failed to allocate vector',me)
c
      do 3 icol=jlo,jhi
      call ga_get(lg_s,ilo,ihi,icol,icol,dbl_mb(k_scr),1)
cx
cx      write(*,6666) icol,(dbl_mb(k_scr+i),i=0,nbf-1)
cx 6666 format(//,'Vector',i5,/,(10f10.6))
cx
      scale=dbl_mb(ioptr+icol-1+ioff)
cx
cx      write(*,'(a,i5,f12.6)') 'Scale ',icol,scale
cx
      call dscal(lthi,scale,dbl_mb(k_scr),1)
      call ga_put(lg_s,ilo,ihi,icol,icol,dbl_mb(k_scr),1)
cx      call ga_get(lg_s,ilo,ihi,icol,icol,dbl_mb(k_scr),1)
cx
cx      write(*,6667) (dbl_mb(k_scr+i),i=0,nbf-1)
cx 6667 format(10f10.6)
cx
    3 continue
c
      if(.not.ma_pop_stack(l_scr))
     + call md_abort('esp_denmat: failed to deallocate vector',me)
c
      endif
c
      call ga_sync()
      call ga_dgemm('n','t',nbf,nbf,nbf,1.d0,lg_v,lg_s,0.d0,lg_p(iset))
c
    1 continue 
c
      if (nset.eq.1) then
        call ga_copy(lg_p(1),lg_d)
      else 
        call ga_dadd(1.0d0,lg_p(1),1.0d0,lg_p(2),lg_d) ! consolidate the density matrices
      endif
c
      endif !denmat
c
      call ga_sync()
c
c     deallocate global arrays lg_s, lg_v and lg_d
c
      if(.not.ga_destroy(lg_p(1)))
     + call md_abort('esp_denmat: failed to destroy lg_p',me)
      if (nset.eq.2) then
        if(.not.ga_destroy(lg_p(2)))
     +  call md_abort('esp_denmat: failed to destroy lg_p',me)
      endif
      if(.not.ga_destroy(lg_s))
     + call md_abort('esp_denmat: failed to destroy lg_s',me)
      if(.not.ga_destroy(lg_v))
     + call md_abort('esp_denmat: ga_destroy lg_v failed',me)
c
      return
      end
