      subroutine tce_mrcc_create_cas1(rtdb)
        implicit none
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_mrcc.fh"
#include "tce_main.fh"
#include "geom.fh"

        integer iactel,iactorb
        integer rtdb
ckbn gf -2      integer nodezero,nprocs
      logical nodezero
        character*4 irrepname
        integer i,j,irrep
c        integer nref,iref
        integer iref
        integer k,n
        integer id,ia,ib
        integer mrcc_factorial
        external mrcc_factorial
        integer itotal
        integer iacounts(maxorb)
        integer ibcounts(maxorb)
        integer iapos(maxorb+1)
        integer ibpos(maxorb+1)
c        integer itotreal,idbocc(10240,maxorb)
        integer itotreal
        character*10240 s
c        character*maxref s
        integer iactshift
        character*8 t1
        character*3 t2
        integer ic
        integer l_idbocc,k_idbocc
c        integer l_sidbocc,k_sidbocc

        nodezero = (ga_nodeid().eq.0)

        do k=1,10240
          s(k:k)=''
        enddo

        t1 = 'bwcc:ref'

      if (.not.rtdb_get(rtdb,'mrcc:iactel',mt_int,1,iactel))
     2   call  errquit('tce_mrcc_cas',2,UNKNOWN_ERR)
      if (.not.rtdb_get(rtdb,'mrcc:iactorb',mt_int,1,iactorb))
     2    call  errquit('tce_mrcc_cas',3,UNKNOWN_ERR)

      
      if (.not.ma_push_get(mt_int,maxorb*maxref,"idbocc",
     +  l_idbocc,k_idbocc))
     +  call errquit("tce_mrcc_energy: MA problem",0,MA_ERR)
c       do i=0,maxorb*10240
c        int_mb(k_idbocc + i-1 ) = 0
c       enddo
c      call ifill(maxorb*maxref,0,int_mb(k_idbocc),1)

c      if (.not.ma_push_get(mt_int,maxorb*maxref,"sidbocc",
c      if (.not.ma_push_get(mt_byte,maxorb*maxref,"sidbocc",
c     +  l_sidbocc,k_sidbocc))
c     +  call errquit("tce_mrcc_energy: MA problem",0,MA_ERR)
c      call ifill(maxorb*maxref,0,int_mb(k_sidbocc),1)

        if(nodezero) then

          write(LuOut,"(/,'Active elec :',I4)")iactel
          write(LuOut,"('Active orbs :',I4,/)")iactorb

          if(mod(iactel,2).eq.1)
     1 write(LuOut,"('Odd number of elec.',I4,/)")iactel/2

        endif

        itotal = 0

        do i=1,min(iactel/2,iactorb)
          id=iactel/2-i+1
          ia=i-1+mod(iactel,2)
          ib=i-1
          itotal = itotal +
     1 mrcc_factorial(iactorb)/(mrcc_factorial(id)*
     4 mrcc_factorial(ia)*mrcc_factorial(ib)*
     5 mrcc_factorial(iactorb-id-ia-ib))

c        if(nodezero)
c     1  write(LuOut,"('Category ',I5,I5,I5,' has ',I8,
c     2 ' possibilities')")id,ia,ib,
c     3 mrfactorial(iactorb)/(mrfactorial(id)*
c     4 mrfactorial(ia)*mrfactorial(ib)*
c     5 mrfactorial(iactorb-id-ia-ib))

        enddo


        do k=1,iactorb
            iacounts(k) = 0
            ibcounts(k) = 0
          if(k.le.id) then
            iacounts(k) = 1
            ibcounts(k) = 1
          else
           if(k.le.(id+ia)) iacounts(k) = 1
           if(k.le.(id+ib)) ibcounts(k) = 1
          endif
        enddo

        if(nodezero) then
          write(LuOut,"(/,'Alpha ',100I1)")(iacounts(k),k=1,iactorb)
          write(LuOut,"('Beta  ',100I1,/)")(ibcounts(k),k=1,iactorb)
        endif

        do k=1,id+ib
          ibpos(k) = k   
        enddo

          ibpos(id+ib+1) = iactorb+1

        itotreal = 0
 8746   continue

        do i=1,id+ib
             if(i.gt.1) then
               if(ibpos(id+ib-i+1).lt.(iactorb-i+1)) then
             ibpos(id+ib-i+1)=ibpos(id+ib-i+1)+1
                 do n=1,i-1
                   ibpos(id+ib-n+1)=ibpos(id+ib-i+1)+i-n
                 enddo
               goto 8746
               else
               goto 8745
               endif
             endif
         do k=ibpos(id+ib-i+1),ibpos(id+ib-i+2)
           if(k.eq.(ibpos(id+ib-i+2))) then
            if((ibpos(id+ib-i+1)+1).lt.(ibpos(id+ib-i+2))) then
             ibpos(id+ib-i+1)=ibpos(id+ib-i+1)+1
             if(i.gt.1) goto 8746
            endif
            goto 8745
           endif
           do n=1,iactorb
            ibcounts(n) = 0
           enddo
           do n=1,id+ib
            if(n.eq.(id+ib+1-i)) then
            ibcounts(k) = 1
            else
            ibcounts(ibpos(n)) = 1
            endif
           enddo
c           if(nodezero)write(LuOut,"(100I1)")(ibcounts(n),n=1,iactorb)
           itotreal = itotreal + 1
           do n=1,iactorb
c             idbocc(itotreal,n)=ibcounts(n)
            int_mb(k_idbocc+(n-1)*maxorb+itotreal-1)=ibcounts(n)
           enddo
         enddo
 8745    continue
        enddo

corg        if(ibpos(1).lt.(iactorb-id-ib)) goto 8746
         do i=1,(nocc(1)+nocc(2)-iactel)/2
           s(i:i)='2'
c           byte_mb(k_sidbocc+(i-1)*maxorb+i-1)='2'
c       write(*,'(A,A,A)') "s",s(i:i),byte_mb(k_sidbocc+(i-1)*maxorb+i-1)
c       write(*,'(A,A)') "s",byte_mb(k_sidbocc+(i-1)*maxorb+i-1)
c           write(*,'(A,A,A)') "s",s
         enddo

         iactshift = (nocc(1)+nocc(2)-iactel)/2

c         write(6,*)iactshift,nocc(1),nocc(2)

         ic = 0

         do i=1,itotreal
          do k=1,itotreal
          irrep = 0
           do n=1,iactorb
c          if((idbocc(i,n)+idbocc(k,n)).eq.2) then
           if(((int_mb(k_idbocc+(n-1)*maxorb+i-1))+
     +            int_mb(k_idbocc+(n-1)*maxorb+k-1)).eq.2) then
           s(n+iactshift:n+iactshift)='2'
c         byte_mb(k_sidbocc+((n+iactshift)-1)*maxorb+(n+iactshift)-1)='2'
c             elseif(idbocc(i,n).eq.1) then
             elseif(int_mb(k_idbocc+(n-1)*maxorb+i-1).eq.1) then
               s(n+iactshift:n+iactshift)='a'
c         byte_mb(k_sidbocc+((n+iactshift)-1)*maxorb+(n+iactshift)-1)='a'
               irrep = ieor(irrep,int_mb(k_irs(1)+n+iactshift-1))
c             elseif(idbocc(k,n).eq.1) then
             elseif((int_mb(k_idbocc+(n-1)*maxorb+k-1)).eq.1) then
               s(n+iactshift:n+iactshift)='b'
c         byte_mb(k_sidbocc+((n+iactshift)-1)*maxorb+(n+iactshift)-1)='b'
               irrep = ieor(irrep,int_mb(k_irs(2)+n+iactshift-1))
             else
               s(n+iactshift:n+iactshift)='0'
c         byte_mb(k_sidbocc+((n+iactshift)-1)*maxorb+(n+iactshift)-1)='0'
             endif
           enddo
           call sym_irrepname(geom,irrep+1,irrepname)
        if(targetsym.eq.irrepname) then

          ic = ic + 1

             if(nodezero)write(LuOut,*)s(1:iactorb+iactshift),' ',
     1 ic
c      if(nodezero) write(LuOut,'(A,A,I5)')
c     +byte_mb(k_sidbocc+((iactorb+iactshift)-1)*maxorb+(1)-1),'ii ',ic

      write(t2,"(I3.3)")ic
      if (.not.rtdb_cput(rtdb,t1//t2,1,s(1:iactorb+iactshift)))
c      if (.not.rtdb_cput(rtdb,t1//t2,1,byte_mb(k_sidbocc)))
     1   call errquit('tce_mrcc_cas: failed writing to rtdb',0,
     2   RTDB_ERR)

           endif

           enddo
         enddo

        if(nodezero)
     1  write(LuOut,"(/,'Total no. of references:',I5,/)")ic

      if (.not.rtdb_put(rtdb,'bwcc:nref',mt_int,1,ic))
     1   call errquit('tce_mrcc_input: failed writing to rtdb',0,
     2   RTDB_ERR)

c       this routine has to be rewritten at some point

ckbn    check whether number of references go beyond maxref
        if(ic .gt. maxref)  then
        if(nodezero)  
     +   write(LuOut,'(A,I5,A,I5)') "Number of references: ",ic,
     +   " greater than maximum number of references: " , maxref
        if(nodezero) call util_flush(LuOut)
        call errquit('tce_mrcc_create_cas: nref > maxref',3,RTDB_ERR)
        endif

c      if (.not.ma_pop_stack(l_sidbocc))
c     +    call errquit("tce_mrcc_energy: MA problem",615,MA_ERR)

      if (.not.ma_pop_stack(l_idbocc))
     +    call errquit("tce_mrcc_energy: MA problem",615,MA_ERR)

        return
        end

        integer function mrcc_factorial(input)
         implicit none
#include "errquit.fh"
         integer input,res
         integer i

c         if(input .lt.  0) 
c     + call errquit("mrcc_factorial: cas input",30,INPUT_ERR)
c
c         if(input .gt. 20) 
c     + call errquit("mrcc_factorial: cas input",30,INPUT_ERR)


         res = 1
         do i=1,input
           res = res*i
         enddo

         if(res .le. 0) 
     + call errquit("mrcc_factorial: cas input",30,INPUT_ERR)
         mrcc_factorial = res

        end


c $Id$
