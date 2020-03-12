      subroutine ducc(rtdb,d_t1,k_t1_offset,d_t2,k_t2_offset,
     &    d_v2,k_v2_offset,ene_orb,nos,nvs,noas,nobs,nvas,nvbs,nactv)
c
c all varaible below have to be defined in the tce_energy.f
c the nact should be defined in the input set tce:nactv
c
c nactv - number of active virtual spinorbitals
c nos   - number of occupied spinorbitals
c nvs   - number of virtual spinorbitals
c noas  - number of occupied alpha spinorbitals
c nobs  - number of occupied beta spinorbitals
c nvas  - number of virtual  alhoa spinorbitals 
c nvbs  - number of virtual beta spinorbitals
c
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer d_t1    ! GA handle for t1
      integer d_t2    ! GA handle for t2
      integer d_v2    ! GA handle for v2 
      integer nos,nvs ! # of occupied/virtual spinorbitals
      integer noas,nobs,nvas,nvbs ! # of occupied/virtual alphas betas
      integer nactv
      double precision ene_orb(nos+nvs)  ! a copy of dbl_mb(k_sorted) 
      integer p1b,p2b,h3b,h4b
      integer p1,p2,h3,h4
      integer k_t1_offset
      integer k_t2_offset
      integer k_v2_offset
      integer size
      integer rtdb
      integer l_aux1,k_aux1 ! for local memory allocator: loc. mem. buffer 1
      integer i,j,k,l         ! auxiliary indices
      integer ia,ib,ic,id,ie,if !auxiliary indices
      integer nact              !number of active virtual orbitals
c matrices with original integrals
c h created here from "orbital" matrix horb
      double precision h(nos+nvs,nos+nvs)
c horb valid only for RHF case
      double precision horb((nos+nvs)/2,(nos+nvs)/2)
      double precision v(nos+nvs,nos+nvs,nos+nvs,nos+nvs)
c t2 in a nice representation : be careful here (you may want to reindex
c virtuals
      double precision t1(nos,nos+1:nos+nvs)
      double precision t2(nos,nos,nos+1:nos+nvs,nos+1:nos+nvs)
c transformed matrices
      double precision ht(nos+nvs,nos+nvs)
      double precision vt(nos+nvs,nos+nvs,nos+nvs,nos+nvs)  
c correlation energies
      double precision eccsd,etfnt
      double precision xxx
c
      logical nodezero
c
c
      nodezero = (ga_nodeid().eq.0)
c
c one electron integrals (original core Hamiltonian)    
c      call kinetic_hcore_1(rtdb,horb)
c       call kinetic_hcore(rtdb,horb,noas+nvas,spust)
c mapping from horb-->h --------------
c
c
      do i=1,(nos+nvs)/2
      do j=1,(nos+nvs)/2
        horb(i,j)=0.0d0
      enddo
      enddo
c
      do i=1,nos+nvs
      do j=1,nos+nvs
        h(i,j)=0.0d0
        ht(i,j)=0.0d0
      enddo
      enddo
c
      call kinetic_hcore_1(rtdb,horb,h,nos,nvs,noas,nobs,nvas,nvbs)
c
      if(nodezero) then
c      xxx=0.0d0
      do i=1,nos+nvs
      do j=1,nos+nvs
       if(dabs(h(i,j)).gt.1.0d-11)  then
           write(6,10) h(i,j),i,j
           call util_flush(6)
       endif
      enddo
      enddo
c      write(6,*)' max h ',xxx
      call util_flush(6)
      endif

c     
c forming v matrix
c
      do i=1,nos+nvs
      do j=1,nos+nvs
      do k=1,nos+nvs
      do l=1,nos+nvs
        v(i,j,k,l)=0.0d0 
        vt(i,j,k,l)=0.0d0
      enddo
      enddo
      enddo
      enddo
      call mapping_v2_m(rtdb,d_v2,k_v2_offset,v,nos,nvs)
c *** debug ****
      if(nodezero) then 
      xxx=0.0d0
      do i=1,nos+nvs
      do j=1,nos+nvs
      do k=1,nos+nvs
      do l=1,nos+nvs
       if(dabs(v(i,j,k,l)).gt.xxx) xxx=dabs(v(i,j,k,l))
      enddo
      enddo
      enddo
      enddo
      write(6,*)'max v ',xxx
      call util_flush(6)
      end if !nodezero
c **************
c
      do i=1,nos
      do ia=nos+1,nos+nvs
        t1(i,ia)=0.0d0
      enddo
      enddo
c
      do i=1,nos
        do j=1,nos
          do ia=nos+1,nos+nvs
            do ib=nos+1,nos+nvs
               t2(i,j,ia,ib)=0.0d0
            enddo
          enddo
        enddo
      enddo
c
      call mapping_t1(rtdb,d_t1,k_t1_offset,t1,nos,nvs)
c
c *** debug ****
      if(nodezero) then 
      write(6,*) 'T1 from ducc'
      do i=1,nos
      do ia=nos+1,nos+nvs
        if(dabs(t1(i,ia)).gt.1.0d-10) then
         write(6,10) t1(i,ia),i,ia
        endif
      enddo
      enddo
      call util_flush(6)
      endif !nodezero
 10   format(f17.10,3x,2i4)
 11   format(f17.10,3x,4i4)
c ***************
c
      call mapping_t2(rtdb,d_t2,k_t2_offset,t2,nos,nvs)
c
c ***  debug ***
      if(nodezero) then
      write(6,*) 'T2 from ducc'
      do i=1,nos
        do j=1,nos
          do ia=nos+1,nos+nvs
            do ib=nos+1,nos+nvs
      if(dabs(t2(i,j,ia,ib)).gt.1.0d-10) then
        write(6,11) t2(i,j,ia,ib),i,j,ia,ib
        call util_flush(6)
      endif
            enddo
          enddo
        enddo
      enddo
      endif !nodezero
c **************
c
c check point: calculating CCSD correlation energy 
c  *** degub **
      write(6,*)'test eccsd'
      call util_flush(6)
c *************
 
      eccsd=0.0d0
      do i=1,nos
      do j=1,nos
      do ia=nos+1,nos+nvs
      do ib=nos+1,nos+nvs
        eccsd=eccsd+0.250d0*v(i,j,ia,ib)*t2(i,j,ia,ib)
        eccsd=eccsd+0.50d0*v(i,j,ia,ib)*t1(i,ia)*t1(j,ib)
c *** debug ***
       write(6,11) v(i,j,ia,ib),i,j,ia,ib
       write(6,11) t2(i,j,ia,ib)
       write(6,10) t1(i,ia)
       write(6,10) t1(j,ib)
       write(6,*)'   '
       call util_flush(6)
c *************
c   
      enddo
      enddo
      enddo
      enddo
c 
      if(nodezero) then 
        write(6,*)'From DUCC CCSD corr. ene.',eccsd
        call util_flush(6)
      endif
c transformation t2-->t2e
      do i=1,nos
      do j=1,nos
c virt a a
        do ia=nos+1,nos+nactv
        do ib=nos+1,nos+nactv
         t2(i,j,ia,ib)=0.0d0
        enddo
        enddo
c virt a b  & b a
        do ia=nos+1,nos+nactv
        do ib=nos+nvas+1,nos+nvas+nactv
         t2(i,j,ia,ib)=0.0d0
         t2(i,j,ib,ia)=0.0d0
        enddo
        enddo
c virt b b 
        do ia=nos+nvas+1,nos+nvas+nactv
        do ib=nos+nvas+1,nos+nvas+nactv
         t2(i,j,ia,ib)=0.0d0
        enddo
        enddo
c
       enddo  !j t2->t2e
       enddo  !i t2->t2e
c calculating <T2e+  F_N  T2e>
c *** debug ***
      write(6,*)'orbital energies'
      do i=1,nos+nvs
       write(6,12) ene_orb(i)
      enddo
      call util_flush(6)
 12   format(f17.10,3x,i3)
c *************
      etfnt=0.0d0
      do i=1,nos
      do j=1,nos
      do ia=nos+1,nos+nvs
      do ib=nos+1,nos+nvs
        etfnt=etfnt+0.50d0*ene_orb(ia)*t2(i,j,ia,ib)*t2(i,j,ia,ib)
        etfnt=etfnt-0.50d0*ene_orb(i)*t2(i,j,ia,ib)*t2(i,j,ia,ib)
      enddo
      enddo
      enddo
      enddo
      if(nodezero) then
        write(6,*)'From DUCC CCSD corr. ene.-e-part',etfnt
        call util_flush(6)
      endif
c
c work starts here 
c
c  calculate ht  (print only ht(alpha,alpha)

c  calculate vt  (print only vt(alpha,beta,alpha,beta)

c 
      return
      end
c
c
c
c
c
c
c

