      SUBROUTINE mapping_v2_m(rtdb,
     1             d_v2,k_v2_offset,v,nos,nvs)
      IMPLICIT NONE
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer rtdb             ! Run-time database
      integer g1b,g2b,g3b,g4b
      integer g1b_1,g2b_1,g3b_1,g4b_1
      integer size_loc,k_a,l_a
      integer d_v2,k_v2_offset
ccx      INTEGER next
ccx      INTEGER nprocs
ccx      INTEGER count
ccx      integer nxtask
ccx      external nxtask
      logical nodezero
      integer igi1,igi2,igi3,igi4,ipos1,ipos2,ipos3,ipos4
      integer i,j,k,l
      integer nos,nvs
      double precision v(nos+nvs,nos+nvs,nos+nvs,nos+nvs)
      double precision xxx
c
c
      nodezero=(ga_nodeid().eq.0)
c
c
ccx      nprocs = GA_NNODES()
ccx      count = 0
ccx      next = NXTASK(nprocs, 1)
c 
      DO g3b = 1,noab+nvab
      DO g4b = g3b,noab+nvab
      DO g1b = 1,noab+nvab
      DO g2b = g1b,noab+nvab
      IF (int_mb(k_spin+g3b-1)+int_mb(k_spin+g4b-1) .eq. int_mb(k_spin+g
     &1b-1)+int_mb(k_spin+g2b-1)) THEN
      IF (ieor(int_mb(k_sym+g3b-1),ieor(int_mb(k_sym+g4b-1),ieor(int_mb(
     &k_sym+g1b-1),int_mb(k_sym+g2b-1)))) .eq. irrep_v) THEN
ckk      IF ((.not.restricted).or.(int_mb(k_spin+g3b-1)+int_mb(k_spin+g4b-1
ckk     &)+int_mb(k_spin+g1b-1)+int_mb(k_spin+g2b-1).ne.8)) THEN ! spin case 8
ccx      IF (next.eq.count) THEN
c
        size_loc = int_mb(k_range+g3b-1)*int_mb(k_range+g4b-1)*
     &             int_mb(k_range+g1b-1)*int_mb(k_range+g2b-1)
       IF (.not.MA_PUSH_GET(mt_dbl,size_loc,'a',l_a,k_a))
     &              CALL ERRQUIT('map1',1,MA_ERR)
       CALL TCE_RESTRICTED_4(g3b,g4b,g1b,g2b,g3b_1,g4b_1,g1b_1,g2b_1)
       call  GET_HASH_BLOCK_I(d_v2,dbl_mb(k_a),size_loc,
     1               int_mb(k_v2_offset),(g2b_1 - 1 + (noab+nvab) *
     2               (g1b_1 - 1 + (noab+nvab) * (g4b_1 - 1 +
     3               (noab+nvab) * (g3b_1 - 1)))),
     4               g2b_1,g1b_1,g4b_1,g3b_1)
c mapping starts here 
       i=0
       do igi3=1,int_mb(k_range+g3b-1)
       do igi4=1,int_mb(k_range+g4b-1)
       do igi1=1,int_mb(k_range+g1b-1)
       do igi2=1,int_mb(k_range+g2b-1)
        i=i+1
        ipos3=int_mb(k_offset+g3b-1)+igi3
        ipos4=int_mb(k_offset+g4b-1)+igi4
        ipos1=int_mb(k_offset+g1b-1)+igi1
        ipos2=int_mb(k_offset+g2b-1)+igi2
        xxx=dbl_mb(k_a+i-1)
        if((ipos1.le.ipos2).and.(ipos3.le.ipos4)) then
         v(ipos3,ipos4,ipos1,ipos2)=xxx
         v(ipos3,ipos4,ipos2,ipos1)=-xxx
         v(ipos4,ipos3,ipos1,ipos2)=-xxx
         v(ipos4,ipos3,ipos2,ipos1)=xxx
         v(ipos1,ipos2,ipos3,ipos4)=xxx
         v(ipos1,ipos2,ipos4,ipos3)=-xxx
         v(ipos2,ipos1,ipos3,ipos4)=-xxx
         v(ipos2,ipos1,ipos4,ipos3)=xxx
        endif
       enddo
       enddo
       enddo
       enddo

c mapping done
ccx       call put_hash_block(d_v2aux,dbl_mb(k_a),size_loc,
ccx     1             int_mb(k_v2aux_offset),(g2b - 1 + (noab+nvab) *
ccx     2               (g1b - 1 + (noab+nvab) * (g4b - 1 +
ccx     3               (noab+nvab) * (g3b - 1)))))
       IF (.not.MA_POP_STACK(l_a))
     1          CALL ERRQUIT('ccsd_t2_1',4,MA_ERR)
c
ccx      next = NXTASK(nprocs, 1)
ccx      END IF
ccx      count = count + 1
ckk      END IF !spin case 8 
      END IF
      END IF
      END DO
      END DO
      END DO
      END DO
ccx      next = NXTASK(-nprocs, 1)
      call ga_sync()
      return 
      end 
c
c
c
c
c
      SUBROUTINE mapping_t1(rtdb,d_t1,k_t1_offset,t1,nos,nvs,
     &noas,nobs,nvas,nvbs)
      IMPLICIT NONE
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer rtdb             ! Run-time database
      integer d_t1
      integer k_t1_offset
      integer nos,nvs,noas,nobs,nvas,nvbs
      integer p1b,h2b
      integer p1b_1,h2b_1
      integer i,iph2,ipp1
      integer p1,h2
      integer size,k_aux,l_aux
      integer ispin
      double precision t1(nos,nos+1:nos+nvs)
      double precision xxx
      do p1b = noab+1,noab+nvab
        do h2b = 1,noab
          if (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+h2b-1)) then
ckk            if ((.not.restricted).or.(int_mb(k_spin+p1b-1)
ckk     1        +int_mb(k_spin+h2b-1).ne.4)) then  ! spin case 4
              if (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+h2b-1))
     1          .eq.irrep_t) then
                ispin=int_mb(k_spin+p1b-1)+int_mb(k_spin+h2b-1)
                size = int_mb(k_range+p1b-1) * int_mb(k_range+h2b-1)
                if (.not.ma_push_get(mt_dbl,size,'taux',l_aux,k_aux))
     1            call errquit('t1map: MA problem',0,MA_ERR)
               CALL TCE_RESTRICTED_2(p1b,h2b,p1b_1,h2b_1) 
                call get_hash_block(d_t1,dbl_mb(k_aux),size,
     1            int_mb(k_t1_offset),((p1b_1-noab-1)*noab+h2b_1-1))
                i = 0
                do p1 = 1,int_mb(k_range+p1b-1)
                  do h2 = 1,int_mb(k_range+h2b-1)
                    i = i + 1
                    xxx=dbl_mb(k_aux+i-1)
                    if(dabs(xxx).gt.1.0d-11) then !filling t1
                       iph2=int_mb(k_offset+h2b-1)+h2
                       ipp1=int_mb(k_offset+p1b-1)+p1
                       t1(iph2,ipp1)=xxx
c only RHF case is covered here beta->beta
c                        if(ispin.eq.2)
c                          iph2=iph2+noas
c                          ipp1=ipp1+nvas 
c                          t1(iph2,ipp1)=xxx
c                        endif
                     endif ! tol 1.0d-11
                  enddo
                enddo
c
                if (.not.ma_pop_stack(l_aux))
     1            call errquit('t1map2: MA problem',1,MA_ERR)
              endif
ckk            endif !spin case 4
          endif
        enddo
      enddo
c
      call ga_sync()
c
      RETURN
      END
c
c
c
c
      SUBROUTINE mapping_t2(rtdb,d_t2,k_t2_offset,t2,nos,nvs)
      IMPLICIT NONE
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer rtdb             ! Run-time database
      integer d_t2
      integer k_t2_offset
      integer nos,nvs
      integer p1b,p2b,h3b,h4b
      integer p1b_1,p2b_1,h3b_1,h4b_1
      integer i,iph4,iph3,ipp2,ipp1
      integer p1,p2,h3,h4
      integer size,k_aux,l_aux
      double precision t2(nos,nos,nos+1:nos+nvs,nos+1:nos+nvs)
      double precision xxx
c
      do p1b = noab+1,noab+nvab
        do p2b = p1b,noab+nvab
          do h3b = 1,noab
            do h4b = h3b,noab
              if (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)
     1          .eq. int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1)) then
ckk              if ((.not.restricted).or.
ckk     1          (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+
ckk     2          int_mb(k_spin+h3b-1)+int_mb(k_spin+h4b-1).ne.8)) then !spin case 8
              if (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),
     1          ieor(int_mb(k_sym+h3b-1),int_mb(k_sym+h4b-1))))
     2          .eq. irrep_t) then
                size = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1)
     1               * int_mb(k_range+h3b-1) * int_mb(k_range+h4b-1)
                if (.not.ma_push_get(mt_dbl,size,'t2aux',l_aux,k_aux))
     1            call errquit('t2map: MA problem',0,MA_ERR)
                call tce_restricted_4(p1b,p2b,h3b,h4b,p1b_1,p2b_1,
     1               h3b_1,h4b_1)
                call get_hash_block(d_t2,dbl_mb(k_aux),size,
     1            int_mb(k_t2_offset),
     1            ((((p1b_1-noab-1)*nvab+p2b_1-noab-1)
     2            *noab+h3b_1-1)*noab+h4b_1-1))
                i = 0
                do p1 = 1,int_mb(k_range+p1b-1)
                  do p2 = 1,int_mb(k_range+p2b-1)
                    do h3 = 1,int_mb(k_range+h3b-1)
                      do h4 = 1,int_mb(k_range+h4b-1)
                        i = i + 1
                        xxx=dbl_mb(k_aux+i-1)
                        if(dabs(xxx).gt.1.0d-11) then !filling t2
                          iph4=int_mb(k_offset+h4b-1)+h4
                          iph3=int_mb(k_offset+h3b-1)+h3
                          ipp2=int_mb(k_offset+p2b-1)+p2
                          ipp1=int_mb(k_offset+p1b-1)+p1
                        if((iph3.lt.iph4).and.(ipp1.lt.ipp2)) then
                          t2(iph3,iph4,ipp1,ipp2)=xxx
                          t2(iph3,iph4,ipp2,ipp1)=-xxx
                          t2(iph4,iph3,ipp1,ipp2)=-xxx
                          t2(iph4,iph3,ipp2,ipp1)=xxx
                        endif !ordering iph3<iph4 ipp1<ipp2
                        endif ! filling t2
                      enddo
                    enddo
                  enddo
                enddo
                if (.not.ma_pop_stack(l_aux))
     1            call errquit('t2map2: MA problem',1,MA_ERR)
              endif
ckk              endif !spin case 8
              endif
c
            enddo  !h4b
          enddo  !h3b
        enddo  !p2b
      enddo  !p1b
c
      call ga_sync()
      RETURN 
      END
c
c
c
      SUBROUTINE kinetic_hcore_1(rtdb,f2m,h2m,
     &  nos,nvs,noas,nobs,nvas,nvbs)
C     $Id: tce_mo2e_zones_4a_disk_ga_N5.F 25705 2014-06-02 18:58:56Z
C     kowalski $
      IMPLICIT NONE
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "util.fh"
#include "bas.fh"
#include "schwarz.fh"
#include "sym.fh"
#include "sf.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer rtdb
      integer g_hcore
      integer nbf2
      integer nos,nvs,noas,nobs,nvas,nvbs
      integer k_core,l_core
      integer iha,ihb,ipa,ipb
      integer isigma1
      integer l_movecs_orb,k_movecs_orb
      integer l_aux,k_aux
      integer dim_mol
      integer l_final,k_final
      integer i,j,k,ipos1,ipos2
      double precision xxx
      double precision f2m((nos+nvs)/2,(nos+nvs)/2)
      double precision h2m(nos+nvs,nos+nvs)
cc      logical oskel
      logical nodezero
      oskel=.false.
      nodezero=(ga_nodeid().eq.0)
c
c f2m & h2m asummed to be zero 
c
      do i=1,(nos+nvs)/2
      do j=1,(nos+nvs)/2
       f2m(i,j)=0.0d0
      enddo
      enddo
c
      do i=1,nos+nvs
      do j=1,nos+nvs
        h2m(i,j)=0.0d0
      enddo
      enddo
c
c
      if (.not. ga_create(mt_dbl,nbf,nbf, 'AO ovl',
     &        0, 0, g_hcore))
     &        call errquit('Error creating ga_hcore',0,GA_ERR)
c
      call ga_zero(g_hcore)
c 
      call int_1e_ga(ao_bas_han,ao_bas_han,g_hcore,'kinetic',oskel)
      call int_1e_ga(ao_bas_han,ao_bas_han,g_hcore,'potential',oskel)
c create MA l_core
      nbf2=nbf*nbf
      if (.not.ma_push_get(mt_dbl,nbf2,'h_core',l_core,k_core))
     1    call errquit('tce_core: MA problem',0,MA_ERR)
      call dfill(nbf2, 0.0d0, dbl_mb(k_core), 1)
c h_hcore ==> l_core ! attention - the last parameter (leading
c             dimension of local memory buffer)
        call ga_get(g_hcore,1,nbf,1,nbf,dbl_mb(k_core),nbf)
c 
c alpha orbitals only
c
        iha = nocc(1)-nfc(1)
        ihb = nocc(ipol)-nfc(ipol)
        ipa = nmo(1)-nocc(1)-nfv(1)
        ipb = nmo(ipol)-nocc(ipol)-nfv(ipol)
c *** debug ****
        if(nodezero) then 
        write(6,*)'from h1'
        write(6,*)'iha ',iha
        write(6,*)'ihb ',ihb
        write(6,*)'ipa ',ipa
        write(6,*)'ipb ',ipb
        call util_flush(6)
        endif
c **************
       iha=noas
       ihb=nobs
       ipa=nvas
       ipb=nvbs
c
        dim_mol=iha+ipa
c
      if (.not.ma_push_get(mt_dbl,nbf*(iha+ipa)
     1  ,"sorted MO coeffs",
     2  l_movecs_orb,k_movecs_orb))
     3  call errquit('tce_mo2e_zone: MA problem 1',0,
     2    BASIS_ERR)
      call dfill(nbf*(iha+ipa),0.0d0, dbl_mb(k_movecs_orb), 1)
      do i=1,iha
      do isigma1=1,nbf
       dbl_mb(k_movecs_orb+(i-1)*nbf+isigma1-1)=
     & dbl_mb(k_movecs_sorted+(i-1)*nbf+isigma1-1)
      enddo
      enddo
      do i=iha+1,iha+ipa
      do isigma1=1,nbf
       dbl_mb(k_movecs_orb+(i-1)*nbf+isigma1-1)=
     & dbl_mb(k_movecs_sorted+(i+ihb-1)*nbf+isigma1-1)
      enddo
      enddo
c
      if (.not.ma_push_get(mt_dbl,nbf*dim_mol,'l_aux',l_aux,k_aux))
     1    call errquit('l_aux_core: MA problem',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,dim_mol*dim_mol,'l_final',
     1    l_final,k_final))
     1    call errquit('tce_core: MA problem',0,MA_ERR)
c   l_aux (nbf,dim_mol)= k_core(nbf,nbf) * k_movecs_orb(nbf,dim_mol)
       call dgemm('N','N',nbf,dim_mol,nbf,1.0d0,dbl_mb(k_core),nbf,
     1            dbl_mb(k_movecs_orb),nbf,0.0d0,dbl_mb(k_aux),nbf)
c   l_final= k_movecs_orb(nbf,dim_mol)^T * l_aux(nbf,dim_mol)
c   l_final= k_movecs_orb_T(dim_mol,nbf) * l_aux(nbf,dim_mol)
       call dgemm('T','N',dim_mol,dim_mol,nbf,1.0d0,
     1            dbl_mb(k_movecs_orb),nbf,dbl_mb(k_aux),
     1            nbf,0.0d0,dbl_mb(k_final),dim_mol)
c
       if(nodezero) then
        write(6,*)'h_1 matrix elements'
        k=0
        do j=1,dim_mol
        do i=1,dim_mol
         k=k+1
         xxx=dabs(dbl_mb(k_final+k-1))
         if(xxx.gt.1.0d-10) then
          write(6,100) i,j,dbl_mb(k_final+k-1)
         end if 
        enddo
        enddo
        call util_flush(6)
       end if !nodezero
c spust f2m --------
        k=0
        do j=1,dim_mol
        do i=1,dim_mol
         k=k+1
         f2m(i,j)=dbl_mb(k_final+k-1)
        enddo
        enddo
c spust f2m --------
c delete k_final
      if (.not.ma_pop_stack(l_final))
     1  call errquit('tcc_hcore_final: MA problem',15,MA_ERR)
c delete k_aux
      if (.not.ma_pop_stack(l_aux))
     1  call errquit('tcc_hcore_aux: MA problem',15,MA_ERR)
c delete k_movecs_orb
      if (.not.ma_pop_stack(l_movecs_orb))
     1  call errquit('tcc_hcore: MA problem',15,MA_ERR)
c delete MA l_core
      if (.not.ma_pop_stack(l_core))
     1     call errquit('l_core',15,MA_ERR)
c
        if (.not.ga_destroy(g_hcore))
     1    call errquit('deletefile: hcore GA problem',0,GA_ERR)
c
c
c
c mapping f2m --> h2m (valid only for RHF).
c
c block h-h
       do i=1,iha
       do j=1,iha
        h2m(i,j)=f2m(i,j)
        h2m(i+iha,j+iha)=f2m(i,j)       
       enddo
       enddo
c block h-p
       do i=1,iha
       do j=1,ipa
         h2m(i,j+iha+ihb)=f2m(i,j+iha)
         h2m(i+iha,j+iha+ihb+ipa)=f2m(i,j+iha)
       enddo
       enddo
c block p-h
       do i=1,ipa
       do j=1,iha
         h2m(i+iha+ihb,j)=f2m(i+iha,j)
         h2m(i+iha+ihb+ipa,j+iha)=f2m(i+iha,j)
       enddo
       enddo
c block p-p
       do i=1,ipa
       do j=1,ipa
        h2m(i+iha+ihb,j+iha+ihb)=f2m(i+iha,j+iha)
        h2m(i+iha+ihb+ipa,j+iha+ihb+ipa)=f2m(i+iha,j+iha)
       enddo
       enddo
c
c
c
 100  format(2i5,3x,f18.10)
c
      call ga_sync()
      return 
      end

