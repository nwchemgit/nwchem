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
      INTEGER next
      INTEGER nprocs
      INTEGER count
      integer nxtask
      external nxtask
      logical nodezero
      integer igi1,igi2,igi3,igi4,ipos1,ipos2,ipos3,ipos4
      integer i,j,k,l,ld
      integer nos,nvs
      double precision v(nos+nvs,nos+nvs,nos+nvs,nos+nvs)
      double precision xxx
      integer howmany,l_elem,k_elem
      external dabssum
      double precision dabssum
      integer mmap_elem
      mmap_elem(i,j,k,l,ld)=(((l-1)*ld+k-1)*ld+j-1)*ld+i-1
c
c
      nodezero=(ga_nodeid().eq.0)
c
c
      nprocs = GA_NNODES()
      count = 0
      howmany=0
c     array to avoid double (or multiple-) counting after dgop
c     sinc elements update is not unique
      IF (.not.MA_PUSH_GET(mt_int,(nos+nvs)**4,'elem',l_elem,k_elem))
     &     CALL ERRQUIT('map1',0,MA_ERR)
      call ifill((nos+nvs)**4,0,int_mb(k_elem),1)
      next = NXTASK(nprocs, 1)
c 
      DO g3b = 1,noab+nvab
      DO g4b = g3b,noab+nvab
      DO g1b = 1,noab+nvab
         IF (next.eq.count) THEN
            howmany=howmany+1
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
           int_mb(k_elem+mmap_elem(ipos3,ipos4,ipos1,ipos2,nos+nvs))=1
           v(ipos3,ipos4,ipos2,ipos1)=-xxx
           int_mb(k_elem+mmap_elem(ipos3,ipos4,ipos2,ipos1,nos+nvs))=1
           v(ipos4,ipos3,ipos1,ipos2)=-xxx
           int_mb(k_elem+mmap_elem(ipos4,ipos3,ipos1,ipos2,nos+nvs))=1
           v(ipos4,ipos3,ipos2,ipos1)=xxx
           int_mb(k_elem+mmap_elem(ipos4,ipos3,ipos2,ipos1,nos+nvs))=1
           v(ipos1,ipos2,ipos3,ipos4)=xxx
           int_mb(k_elem+mmap_elem(ipos1,ipos2,ipos3,ipos4,nos+nvs))=1
           v(ipos1,ipos2,ipos4,ipos3)=-xxx
           int_mb(k_elem+mmap_elem(ipos1,ipos2,ipos4,ipos3,nos+nvs))=1
           v(ipos2,ipos1,ipos3,ipos4)=-xxx
           int_mb(k_elem+mmap_elem(ipos2,ipos1,ipos3,ipos4,nos+nvs))=1
           v(ipos2,ipos1,ipos4,ipos3)=xxx
           int_mb(k_elem+mmap_elem(ipos2,ipos1,ipos4,ipos3,nos+nvs))=1
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
ccx      count = count + 1
ckk      END IF !spin case 8 
      END IF
      END IF
      END DO
      next = NXTASK(nprocs, 1)
      END IF
      count = count + 1
      END DO
      END DO
      END DO
      next = NXTASK(-nprocs, 1)
c     global sum
c      call ga_dgop(2019,v,(nos+nvs)**4,'+')
      call util_mygadgop(2019,v,(nos+nvs)**4,'+')
      call util_mygaigop(2020,int_mb(k_elem),(nos+nvs)**4,'+')
      call map_elemd(v,int_mb(k_elem),(nos+nvs)**4)
       IF (.not.MA_POP_STACK(l_elem))
     1          CALL ERRQUIT('map_hcore',14,MA_ERR)
c      call ga_sync()
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
      SUBROUTINE mapping_t3(rtdb,d_t3,k_t3_offset,t3,nos,nvs)
      IMPLICIT NONE
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
       integer rtdb             ! Run-time database
       integer d_t3
       integer k_t3_offset
       integer nos,nvs
       integer p1b,p2b,p3b,h4b,h5b,h6b
       integer p1b_1,p2b_1,p3b_1,h4b_1,h5b_1,h6b_1
       integer i,key
       integer iph6,iph5,iph4,ipp3,ipp2,ipp1
       integer p1,p2,p3,h4,h5,h6
       integer size,k_aux,l_aux
       double precision t3(nos,nos,nos,
     &  nos+1:nos+nvs,nos+1:nos+nvs,nos+1:nos+nvs)
      double precision xxx
c
      DO p1b = noab+1,noab+nvab
      DO p2b = p1b,noab+nvab
      DO p3b = p2b,noab+nvab
      DO h4b = 1,noab
      DO h5b = h4b,noab
      DO h6b = h5b,noab
      IF (int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1)+int_mb(k_spin+p3b-1)
     & .eq. int_mb(k_spin+h4b-1)+int_mb(k_spin+h5b-1)+int_mb(k_spin+h6b-
     &1)) THEN
      IF (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p2b-1),ieor(int_mb(
     &k_sym+p3b-1),ieor(int_mb(k_sym+h4b-1),ieor(int_mb(k_sym+h5b-1),int
     &_mb(k_sym+h6b-1)))))) .eq. irrep_t) THEN
cuu      IF ((.not.restricted).or.(int_mb(k_spin+p2b-1)+int_mb(k_spin+p3b-1
cuu     &)+int_mb(k_spin+p4b-1)+int_mb(k_spin+h1b-1)+int_mb(k_spin+h5b-1)+i
cuu     &nt_mb(k_spin+h6b-1).ne.12)) THEN
c
      size =   int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1)
     1       * int_mb(k_range+p3b-1) * int_mb(k_range+h4b-1)
     1       * int_mb(k_range+h5b-1) * int_mb(k_range+h6b-1)
c
       if (.not.ma_push_get(mt_dbl,size,'t3aux',l_aux,k_aux))
     1      call errquit('t2map: MA problem',0,MA_ERR)
c
      key= h6b - 1 + noab * (h5b - 1 + noab * (h4b
     &-1+noab * (p3b - noab - 1 + nvab * (p2b - noab - 1 + nvab * (p1b
     & - noab - 1)))))
c
      call tce_restricted_6(p1b,p2b,p3b,h4b,h5b,h6b,p1b_1,p2b_1,p3b_1,
     &h4b_1,h5b_1,h6b_1)
c
      call get_hash_block(d_t3,dbl_mb(k_aux),size,
     1     int_mb(k_t3_offset),key)
c
       i=0   
       do p1 = 1,int_mb(k_range+p1b-1)
       do p2 = 1,int_mb(k_range+p2b-1)
       do p3 = 1,int_mb(k_range+p3b-1)
       do h4 = 1,int_mb(k_range+h4b-1)
       do h5 = 1,int_mb(k_range+h5b-1)
       do h6 = 1,int_mb(k_range+h6b-1)
          i = i + 1
          xxx=dbl_mb(k_aux+i-1)
            if(dabs(xxx).gt.1.0d-11) then !filling t2
               ipp1=int_mb(k_offset+p1b-1)+p1
               ipp2=int_mb(k_offset+p2b-1)+p2
               ipp3=int_mb(k_offset+p3b-1)+p3
               iph4=int_mb(k_offset+h4b-1)+h4
               iph5=int_mb(k_offset+h5b-1)+h5
               iph6=int_mb(k_offset+h6b-1)+h6
            if((ipp1.lt.ipp2).and.(ipp2.lt.ipp3).and.
     &         (iph4.lt.iph5).and.(iph5.lt.iph6)) then
c
c ipp1,ipp2,ipp3 ->  1
c ipp1,ipp3,ipp2 -> -1
c ipp2,ipp1,ipp3 -> -1
c ipp2,ipp3,ipp1 ->  1
c ipp3,ipp1,ipp2 ->  1
c ipp3,ipp2,ipp1 -> -1
c
c iph4,iph5,iph6 ->  1
c iph4,iph6,iph5 -> -1
c iph5,iph4,iph6 -> -1
c iph5,iph6,iph4 ->  1
c iph6,iph4,iph5 ->  1
c iph6,iph5,iph4 -> -1
c 
c iph1,iph2,iph3 base +
           t3(iph4,iph5,iph6,ipp1,ipp2,ipp3) =  xxx
           t3(iph4,iph5,iph6,ipp1,ipp3,ipp2) = -xxx
           t3(iph4,iph5,iph6,ipp2,ipp1,ipp3) = -xxx
           t3(iph4,iph5,iph6,ipp2,ipp3,ipp1) =  xxx
           t3(iph4,iph5,iph6,ipp3,ipp1,ipp2) =  xxx
           t3(iph4,iph6,iph6,ipp3,ipp2,ipp1) = -xxx
c iph1,iph3,iph2 base -
           t3(iph4,iph6,iph5,ipp1,ipp2,ipp3) = -xxx
           t3(iph4,iph6,iph5,ipp1,ipp3,ipp2) =  xxx
           t3(iph4,iph6,iph5,ipp2,ipp1,ipp3) =  xxx
           t3(iph4,iph6,iph5,ipp2,ipp3,ipp1) = -xxx
           t3(iph4,iph6,iph5,ipp3,ipp1,ipp2) = -xxx
           t3(iph4,iph6,iph5,ipp3,ipp2,ipp1) =  xxx
c iph2,iph1,iph3 base -
           t3(iph5,iph4,iph6,ipp1,ipp2,ipp3) = -xxx
           t3(iph5,iph4,iph6,ipp1,ipp3,ipp2) =  xxx
           t3(iph5,iph4,iph6,ipp2,ipp1,ipp3) =  xxx
           t3(iph5,iph4,iph6,ipp2,ipp3,ipp1) = -xxx
           t3(iph5,iph4,iph6,ipp3,ipp1,ipp2) = -xxx
           t3(iph5,iph4,iph6,ipp3,ipp2,ipp1) =  xxx
c iph2,iph3,iph1 base +
           t3(iph5,iph6,iph4,ipp1,ipp2,ipp3) =  xxx
           t3(iph5,iph6,iph4,ipp1,ipp3,ipp2) = -xxx
           t3(iph5,iph6,iph4,ipp2,ipp1,ipp3) = -xxx
           t3(iph5,iph6,iph4,ipp2,ipp3,ipp1) =  xxx
           t3(iph5,iph6,iph4,ipp3,ipp1,ipp2) =  xxx
           t3(iph5,iph6,iph4,ipp3,ipp2,ipp1) = -xxx
c iph3,iph1,iph2 base +
           t3(iph6,iph4,iph5,ipp1,ipp2,ipp3) =  xxx
           t3(iph6,iph4,iph5,ipp1,ipp3,ipp2) = -xxx
           t3(iph6,iph4,iph5,ipp2,ipp1,ipp3) = -xxx
           t3(iph6,iph4,iph5,ipp2,ipp3,ipp1) =  xxx
           t3(iph6,iph4,iph5,ipp3,ipp1,ipp2) =  xxx
           t3(iph6,iph4,iph5,ipp3,ipp2,ipp1) = -xxx
c iph3,iph2,iph1 base -
           t3(iph6,iph5,iph4,ipp1,ipp2,ipp3) = -xxx
           t3(iph6,iph5,iph4,ipp1,ipp3,ipp2) =  xxx
           t3(iph6,iph5,iph4,ipp2,ipp1,ipp3) =  xxx
           t3(iph6,iph5,iph4,ipp2,ipp3,ipp1) = -xxx
           t3(iph6,iph5,iph4,ipp3,ipp1,ipp2) = -xxx
           t3(iph6,iph5,iph4,ipp3,ipp2,ipp1) =  xxx
c
            endif !ordering iph4<iph5<iph6 etc.
            endif ! filling t3
           enddo
           enddo
           enddo
           enddo
           enddo
           enddo
c
      if (.not.ma_pop_stack(l_aux))
     1   call errquit('t2map2: MA problem',1,MA_ERR)
c
ccx      END IF ! one less
      END IF
      END IF
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO


c
      return 
      end
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
       if((noas.ne.iha).or.(nobs.ne.ihb).or.(nvas.ne.ipa).or.
     &    (nvbs.ne.ipb)) 
     & call errquit('tce_hcore:inconsistency1',0,MA_ERR)
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

      subroutine map_elemd(v,k_elem,n)
      implicit none
      double precision v(*)
      integer k_elem(*)
      integer n
c
      integer i
c
      do i=1,n
         if(dabs(v(i)).gt.1d-30) v(i)=v(i)/dble(k_elem(i))
      enddo
      return
      end
c
c
      subroutine hn_form(rtdb,
     &    d_v2,k_v2_offset,ene_orb,nos,nvs,noas,nobs,nvas,nvbs,
     &    ehfa)
c
c  !!!!!!!!  Only for RHF !!!!!!!!
c
c all varaible below have to be defined in the tce_energy.f
c
c nos   - number of occupied spinorbitals
c nvs   - number of virtual spinorbitals
c noas  - number of occupied alpha spinorbitals
c nobs  - number of occupied beta spinorbitals
c nvas  - number of virtual  alhoa spinorbitals
c nvbs  - number of virtual beta spinorbitals
c
c spinorbital convention
c
c | noas | nobs | nvas | nvbs |
c
c ene_orb contains orbital energies
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
      integer d_v2    ! GA handle for v2
      integer k_v2_offset
      integer nos,nvs ! # of occupied/virtual spinorbitals
      integer noas,nobs,nvas,nvbs ! # of occupied/virtual alphas betas
      double precision ene_orb(nos+nvs)  ! a copy of dbl_mb(k_sorted)
      integer size
      integer rtdb
      integer i,j,k,l         ! auxiliary indices
      integer ia,ib,ic,id,ie,if !auxiliary indices
      integer p,q,r,s,m,n
c
c matrices with original integrals
c h created here from "orbital" matrix horb
      double precision h(nos+nvs,nos+nvs)
c horb valid only for RHF case
      double precision horb((nos+nvs)/2,(nos+nvs)/2)
c
      double precision v(nos+nvs,nos+nvs,nos+nvs,nos+nvs)
c
      double precision xxx
c shift energies
      double precision ehfa
c
      logical nodezero
c
c
      nodezero = (ga_nodeid().eq.0)
c
      do i=1,nos+nvs
      do j=1,nos+nvs
        h(i,j)=0.0d0
      enddo
      enddo
c
      do i=1,(nos+nvs)/2
      do j=1,(nos+nvs)/2
        horb(i,j)=0.0d0
      enddo
      enddo
c
c forming v matrix
c
      do i=1,nos+nvs
      do j=1,nos+nvs
      do k=1,nos+nvs
      do l=1,nos+nvs
        v(i,j,k,l)=0.0d0 
      enddo
      enddo
      enddo
      enddo
c
      call mapping_v2_m(rtdb,d_v2,k_v2_offset,v,nos,nvs)
c
      do p=1,nos+nvs
       h(p,p)=ene_orb(p)
      enddo
c
      do p=1,nos+nvs
      do q=1,nos+nvs
      do m=1,nos
       h(p,q)=h(p,q)-v(p,m,q,m)
      enddo
      enddo
      enddo
c
      ehfa=0.0d0
c
      do m=1,nos
       ehfa=ehfa+ene_orb(m)
      enddo
c
      do m=1,nos
      do n=1,nos
       ehfa=ehfa-(0.50d0)*v(m,n,m,n) 
      enddo
      enddo
c
c mapping h --> horb (valid only for RHF).
c
c block h-h
       do i=1,noas 
       do j=1,noas
        horb(i,j)=h(i,j)
       enddo
       enddo
c block h-p
       do i=1,noas
       do j=1,nvas
         horb(i,j+noas)=h(i,j+noas+nobs)
       enddo
       enddo
c block p-h
       do i=1,nvas
       do j=1,noas
         horb(i+noas,j)=h(i+noas+nobs,j)
       enddo
       enddo
c block p-p
       do i=1,nvas
       do j=1,nvas
        horb(i+noas,j+noas)=h(i+noas+nobs,j+noas+nobs)
       enddo
       enddo
c
      if(nodezero) then 
        do j=1,noas+nvas
        do i=1,noas+nvas
         xxx=dabs(horb(i,j))
         if(xxx.gt.1.0d-10) then
            write(6,100) i,j,horb(i,j)
         endif
        enddo
        enddo
      endif ! nodezero
c
 100  format(2i5,3x,f18.10)
c
      return
      end
c
c
c
c
      SUBROUTINE mapping_v2(rtdb,d_v2aux,k_v2aux_offset,
     1             d_v2,k_v2_offset)
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
      integer d_v2aux,k_v2aux_offset
      integer d_v2,k_v2_offset
      INTEGER next
      INTEGER nprocs
      INTEGER count
      integer nxtask
      external nxtask
      logical nodezero
c *** debug ***
      integer ig1,ig2,ig3,ig4,i,ipos1,ipos2,ipos3,ipos4
c *************
c
c
      nodezero=(ga_nodeid().eq.0)
c
c
      nprocs = GA_NNODES()
      count = 0
      next = NXTASK(nprocs, 1)
c 
c *** debug ****
c       write(6,*)"2-electron integrals mapping xxx"
c **************
      DO g3b = 1,noab+nvab
      DO g4b = g3b,noab+nvab
      DO g1b = 1,noab+nvab
      DO g2b = g1b,noab+nvab
      IF (int_mb(k_spin+g3b-1)+int_mb(k_spin+g4b-1) .eq. int_mb(k_spin+g
     &1b-1)+int_mb(k_spin+g2b-1)) THEN
      IF (ieor(int_mb(k_sym+g3b-1),ieor(int_mb(k_sym+g4b-1),ieor(int_mb(
     &k_sym+g1b-1),int_mb(k_sym+g2b-1)))) .eq. irrep_v) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+g3b-1)+int_mb(k_spin+g4b-1
     &)+int_mb(k_spin+g1b-1)+int_mb(k_spin+g2b-1).ne.8)) THEN
      IF (next.eq.count) THEN
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
c *** debug ***
c         write(6,*)"2-electron integrals mapping"
c         if(int_mb(k_spin+g3b-1)+int_mb(k_spin+g4b-1).eq.3) then 
c         i=1
c         do ig3=1,int_mb(k_range+g3b-1)
c         do ig4=1,int_mb(k_range+g4b-1)
c         do ig1=1,int_mb(k_range+g1b-1)
c         do ig2=1,int_mb(k_range+g2b-1)
c          ipos3=int_mb(k_offset+g3b-1)+ig3
c          ipos4=int_mb(k_offset+g4b-1)+ig4
c          ipos1=int_mb(k_offset+g1b-1)+ig1
c          ipos2=int_mb(k_offset+g2b-1)+ig2
c          write(6,1000) dbl_mb(k_a+i-1),ipos3,ipos4,ipos1,ipos2
c          i=i+1 
c         enddo
c         enddo
c         enddo
c         enddo
c         call util_flush
c         end if ! spin alpha beta alpha beta
c 1000    format(f15.6,2x,4i6)
c *************
       call put_hash_block(d_v2aux,dbl_mb(k_a),size_loc,
     1             int_mb(k_v2aux_offset),(g2b - 1 + (noab+nvab) *
     2               (g1b - 1 + (noab+nvab) * (g4b - 1 +
     3               (noab+nvab) * (g3b - 1)))))
       IF (.not.MA_POP_STACK(l_a))
     1          CALL ERRQUIT('ccsd_t2_1',4,MA_ERR)
c
      next = NXTASK(nprocs, 1)
      END IF
      count = count + 1
      END IF
      END IF
      END IF
      END DO
      END DO
      END DO
      END DO
      next = NXTASK(-nprocs, 1)
      call ga_sync()
      return
      end

c
c
c
      subroutine pw_tce(rtdb,d_v2,k_v2_offset,d_f1,k_f1_offset,
     &ene_orb,nos,nvs,noas,nobs,nvas,nvbs,mov,ehfa)
c
c Currently for RHF reference
c
c tilesize=1 (assumed)!!!!!!
c
c no-frozen occupied orbitals
c
c nos  - total # of occpuied spinorbitals
c nvs  - total # of virtual  spinorbitals
c noas - total # of alpha occupied spinorbitals
c nobs - total # of beta  occupied spinorbitals
c nvas - total # of alpha virtual  spinorbitals
c nvbs - total # of beta  virtual  spinorbitals
c mov  - total # of correlated ORBITALS = noas+nvas
c
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "inp.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
c
c
      integer rtdb             ! Run-time database
      integer d_v2
      integer k_v2_offset
      integer d_f1
      integer k_f1_offset
      integer nos,nvs,mov
      integer noas,nobs,nvas,nvbs
c vm - 2-electron integrals in the Mulliken convention 
c vm(i,j,k,l)=(ij|kl)
      double precision vm(mov,mov,mov,mov)
      double precision hm(mov,mov)
      double precision fm(mov,mov)
      double precision fms(nos+nvs,nos+nvs)
      double precision ene_orb(nos+nvs)
      logical nodezero         ! True if node 0
!QuantumComputing
      integer l_aux,k_aux
      integer g3b,g4b,g1b,g2b,g5b,g6b
      integer irow,icol
      integer index_pair
      double precision xxx
      double precision enrep_tce
      double precision ehfa
      integer i,j,k,l,m,n,p,q,r,s
      integer iha,ihb,ipa,ipb
      integer hlimit,vlimit
      integer offset_test,size_test
!QuantumComputing
c
      nodezero=ga_nodeid().eq.0
      ehfa = 0.0d0
c
      do g3b=1,noab+nvab
       if(int_mb(k_range+g3b-1).ne.1)
     &      call errquit('TCE_PW:tile eq 1 ',geom,
     &      GEOM_ERR)
      enddo
c
      iha=noas
      ihb=nobs
      ipa=nvas
      ipb=nvbs
c
c ehfa
           if (.not. geom_nuc_rep_energy(geom, enrep_tce))
     &      call errquit('TCE:unable to get nuclear rep ',geom, 
     &      GEOM_ERR)
c ene_orb initialization
           do i=1,nos+nvs
            ene_orb(i)=0.0d0
           enddo
c hm initialization 
           do i=1,mov
           do j=1,mov
            hm(i,j)=0.0d0
            fm(i,j)=0.0d0
           enddo
           enddo
c
           do i=1,nos+nvs
           do j=1,nos+nvs
            fms(i,j)=0.0d0
           enddo
           enddo
c vm
           do i=1,mov
           do j=1,mov
           do k=1,mov
           do l=1,mov
            vm(i,j,k,l)=0.0d0
           enddo
           enddo
           enddo
           enddo
c
c
c opening of the files : reading from nodezero
c   
      if(nodezero) then ! ------------- 987
       open(120,file='hm.dat',status='unknown') 
         read(120,*) hlimit
         do m=1,hlimit
          read(120,*) i,j,xxx
          hm(i,j)=xxx
         enddo
       close(120)
       open(121,file='vm.dat',status='unknown')
         read(121,*) vlimit
         do m=1,vlimit
          read(121,*) i,j,k,l,xxx
          vm(i,j,k,l)=xxx
          vm(i,j,l,k)=xxx
          vm(j,i,k,l)=xxx
          vm(j,i,l,k)=xxx
          vm(k,l,i,j)=xxx
          vm(k,l,j,i)=xxx
          vm(l,k,i,j)=xxx
          vm(l,k,j,i)=xxx
         enddo
       close(121)
c orbital fock matrix
      do i=1,mov
      do j=1,mov
       fm(i,j)=hm(i,j)
      enddo
      enddo
c
      do p=1,mov
      do q=1,mov
       do m=1,noas
        fm(p,q)=fm(p,q)+(2.0d0)*vm(p,q,m,m)-vm(p,m,q,m)
       enddo
      enddo
      enddo
c
c       do p=1,mov
c       do q=1,mov
c        xxx=fm(p,q)
c        if(dabs(xxx).gt.1.0d-10) then
c          write(6,100) p,q,fm(p,q)
c        endif
c       enddo
c       enddo
c
c mapping fm --> fms
c block h-h
       do i=1,iha
       do j=1,iha
        fms(i,j)=fm(i,j)
        fms(i+iha,j+iha)=fm(i,j)
       enddo
       enddo
c block h-p
       do i=1,iha
       do j=1,ipa
         fms(i,j+iha+ihb)=fm(i,j+iha)
         fms(i+iha,j+iha+ihb+ipa)=fm(i,j+iha)
       enddo
       enddo
c block p-h
       do i=1,ipa
       do j=1,iha
         fms(i+iha+ihb,j)=fm(i+iha,j)
         fms(i+iha+ihb+ipa,j+iha)=fm(i+iha,j)
       enddo
       enddo
c block p-p
       do i=1,ipa
       do j=1,ipa
        fms(i+iha+ihb,j+iha+ihb)=fm(i+iha,j+iha)
        fms(i+iha+ihb+ipa,j+iha+ihb+ipa)=fm(i+iha,j+iha)
       enddo
       enddo
c
c
c *** debug ***
       write(6,*)'Fock matrix'
       do i=1,nos+nvs
       do j=1,nos+nvs
        write(6,*) fms(i,j),i,j
       enddo
       enddo
       write(6,*)'--- end of F matrix ---'
c *************
c
c pw - tce works only for rhf and tilesize=1
c      further improvements possible
c
      offset_test=0
      DO g3b = 1,noa+nva   !j
      DO g4b = g3b,noa+nva !i
      DO g1b = 1,noa+nva   !l
      DO g2b = g1b,noa+nva !k
      IF (int_mb(k_spin_alpha+g3b-1)+int_mb(k_spin_alpha+g4b-1).eq.
     &int_mb(k_spin_alpha+g1b-1)+int_mb(k_spin_alpha+g2b-1)) THEN
      IF (ieor(int_mb(k_sym_alpha+g3b-1),ieor(int_mb(k_sym_alpha+g4b-1),
     &    ieor(int_mb(k_sym_alpha+g1b-1),int_mb(k_sym_alpha+g2b-1))))
     &    .eq. irrep_v) THEN
      IROW=INDEX_PAIR(g4b,g3b)
      ICOL=INDEX_PAIR(g2b,g1b)
      IF(IROW.GE.ICOL) THEN
      size_test=int_mb(k_range_alpha+g3b-1) *
     &int_mb(k_range_alpha+g4b-1) * int_mb(k_range_alpha+g1b-1) *
     &int_mb(k_range_alpha+g2b-1)
c
       if(.not.ma_push_get(mt_dbl,size_test,'test',l_aux,k_aux))
     1     call errquit('tce_aux2: MA problem',0,MA_ERR)
c
       dbl_mb(k_aux)=vm(g4b,g3b,g2b,g1b)
c
        call ga_put(d_v2,offset_test+1,offset_test+size_test,1,1,
     1    dbl_mb(k_aux),size_test)
c
      if(.not.ma_pop_stack(l_aux))
     1     call errquit('l_aux2',15,MA_ERR)
c
      offset_test=offset_test+size_test
c
      END IF
      END IF
      END IF
      END DO
      END DO
      END DO
      END DO
c
c
c mapping f1
c
c
      offset_test=0
      DO g6b = 1,noab+nvab
      DO g5b = 1,noab+nvab
      IF (int_mb(k_spin+g6b-1) .eq. int_mb(k_spin+g5b-1)) THEN
      IF (ieor(int_mb(k_sym+g6b-1),int_mb(k_sym+g5b-1)) .eq. irrep_f) TH
     &EN
      IF ((.not.restricted).or.(int_mb(k_spin+g6b-1)+int_mb(k_spin+g5b-1
     &).ne.4)) THEN
c
      size_test = int_mb(k_range+g6b-1) * int_mb(k_range+g5b-1)
c
       if(.not.ma_push_get(mt_dbl,size_test,'test',l_aux,k_aux))
     1     call errquit('tce_aux2: MA problem',0,MA_ERR)
c
       dbl_mb(k_aux)=fms(g6b,g5b)
       call ga_put(d_f1,offset_test+1,offset_test+size_test,1,1,
     1    dbl_mb(k_aux),size_test)
c
      if (.not.ma_pop_stack(l_aux))
     1     call errquit('l_aux2',15,MA_ERR)
c
      offset_test=offset_test+size_test
c
      END IF
      END IF
      END IF
      END DO
      END DO
c
c
c mapping ene_orb -->k_evl_sorted
c
c
      do i=1,nos+nvs
       ene_orb(i)=fms(i,i)
      enddo
c *** debug ***
      write(6,*)'diag fock matrix'
      do i=1,nos+nvs
        write(6,*) ene_orb(i),i
      enddo
      write(6,*)'repulsion energy: ',enrep_tce
c *************
c
c
c     ehfa for RHF reference
c
c
      ehfa=0.0d0
      ehfa=ehfa+enrep_tce 
      do i=1,noas
       ehfa=ehfa+(2.0d0)*hm(i,i)
      enddo
      do i=1,noas
      do j=1,noas
       ehfa=ehfa+2.0d0*vm(i,i,j,j)
       ehfa=ehfa-vm(i,j,i,j)
      enddo
      enddo

c *** debug ***
         write(6,*)'from pw_tce routine ehfa: ',ehfa
         call util_flush(6)
c *************
c
c only nodezer has a value of ehfa !!!!
c
        if (.not.rtdb_put(rtdb,'tce:pwref',mt_dbl,1,ehfa))
     1  call errquit('pwtce: failed to write ehfa',0,RTDB_ERR)
c
       call util_flush(6)
      endif !nodezero ---------------- 987 
c
      call ga_sync()
c
      return
      end



