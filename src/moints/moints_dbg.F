#ifdef DEBUG_PRINT
       subroutine moints_print_opermatrix(noper,nbf,g_a)
       implicit none
#include "global.fh"
#include "mafdecls.fh"       
       integer noper,nbf
       integer g_a
       integer j,k,jk
       integer my_id, ilo, ihi, jlo, jhi
       integer k_local, ld_local

       my_id = ga_nodeid()
       call ga_distribution(g_a,my_id,ilo,ihi,jlo,jhi)
       do j=1,noper
         do k=1,j
           jk = (j*(j-1))/2 + k
           if ((jk.ge.jlo).and.(jk.le.jhi)) then
             write(6,901) j,k
 901         format(//,'Operator: [',i2,',',i2,']',/)
             print*,my_id,'  ',ilo,ihi,'  ',jk
             call ga_access(g_a,ilo,ihi,jk,jk,k_local,ld_local)
             call moints_print_square(nbf,dbl_mb(k_local))
           endif
           call ga_sync()
         enddo
       enddo

       return
       end
       



       subroutine moints_print_square(n,x)
       implicit none
       integer n,i,j
       double precision x(n,n)
       
       do i=1,n
         write(6,901) (x(i,j),j=1,n)
 901     format(15f9.4)
       enddo
       call util_flush(6)
       end
       





       subroutine moints_print_1(ilo,ihi,jlo,jhi,
     $                           klo,khi,mo_lo,mo_hi,x)
       implicit none
       integer ilo,ihi,jlo,jhi,klo,khi,mo_lo,mo_hi
       double precision x(klo:khi,jlo:jhi,ilo:ihi,mo_lo:mo_hi)
       integer m,i,j,k

       do m=mo_lo,mo_hi
         do i=ilo,ihi
           write(6,933) ((x(k,j,i,m),k=klo,khi),j=jlo,jhi)
 933       format(15f9.4)
         enddo
         write(6,*)
       enddo

       return
       end



#else
       subroutine moints_print_opermatrix(noper,nbf,g_a)
       implicit none
       integer noper,nbf
       integer g_a

       return
       end

#endif

