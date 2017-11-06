      program wfn1_2rdm
      implicit none
c
c     Given an 2RDM Gamma_{abcd} consider every possible ordering 
c     of the labels a, b, c, and d and work the sign out for that
c     ordering relative to the ordering abcd. Note that an 
c     interchange of labels 1 and 2, or 3 and 4, or 1 and 4, or 
c     2 and 3 change the sign (assuming real valued 2RDMs). Other
c     interchanges leave the sign unchanged.
c
      integer nmax
      parameter(nmax = 24)
      character*1 labels(4,nmax)
      character*1 word(4)
      character*1 temp
      character*1 q
      parameter(q="'")
      integer perm(4,4)
      integer sign(nmax)
      integer ncol ! the column number
      integer icol ! the column counter
      integer ia   ! counter for label "a"
      integer ib   ! counter for label "b"
      integer ic   ! counter for label "c"
      integer id   ! counter for label "d"
      integer ii,jj
c
      sign = 1
      perm = 1
      perm(1,3) = -1
      perm(1,4) = -1
      perm(2,3) = -1
      perm(2,4) = -1
      perm(3,1) = -1
      perm(3,2) = -1
      perm(4,1) = -1
      perm(4,2) = -1
      ncol = 0
      do ia = 1, 4
        do ib = 1, 4
          if (ib.ne.ia) then
            do ic = 1, 4
              if (ic.ne.ia.and.ic.ne.ib) then
                do id = 1, 4
                  if (id.ne.ia.and.id.ne.ib.and.id.ne.ic) then
                    ncol = ncol + 1
                    labels(ia,ncol) = "a"
                    labels(ib,ncol) = "b"
                    labels(ic,ncol) = "c"
                    labels(id,ncol) = "d"
                  endif
                enddo ! id
              endif
            enddo ! ic
          endif
        enddo ! ib
      enddo ! ia
c
      do icol = 1, ncol
        word(1:4) = labels(1:4,icol)
        do ii = 1, 3
          do jj = ii+1, 4
            if (word(jj).lt.word(ii)) then
              temp     = word(jj)
              word(jj) = word(ii)
              word(ii) = temp
              sign(icol) = sign(icol)*perm(ii,jj)
            endif
          enddo
        enddo
      enddo
c
      do icol = 1, ncol
        write(*,*)"\Gamma^{\sigma\sigma}_{abcd} &=&",sign(icol),
     &            "\Gamma^{\sigma\sigma}_{",labels(1:4,icol),"} \\\\"
        if (icol.eq.1) cycle
        write(80,'("              if(abs(rdm2(ia,ib,ic,id)+(",i2,
     &             ")*rdm2(i",a,",i",a,",i",a,",i",a,")).gt.tol)")')
     &             -sign(icol),labels(1,icol),labels(2,icol),
     &              labels(3,icol),labels(4,icol)
        write(80,'("     &        then")')
        write(80,'("                write(*,*)",a,"failed ",
     &             i2," ",5a)')q,sign(icol),labels(1:4,icol),q
        write(80,'("                write(*,*)ia,ib,ic,id")')
        write(80,'("                write(*,*)rdm2(ia,ib,ic,id),",
     &             "rdm2(i",a,",i",a,",i",a,",i",a,")")')
     &             labels(1:4,icol)
        write(80,'("              endif")')
      enddo
c
      end
