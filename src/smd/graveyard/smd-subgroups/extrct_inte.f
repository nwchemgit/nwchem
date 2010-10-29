c
c $Id$
c

       SUBROUTINE ex_inter(n,ist,irecord,jbuf,ierr)

       implicit none

       integer i,j,n,ifact,ist,nstart,ierr
       integer jbuf,irec_len

       character*1 cchar,irecord

       parameter(irec_len=100)

       dimension cchar(15),irecord(irec_len)

       data cchar/'0','1','2','3','4','5','6','7','8','9','+','&','/',
     $            '-','.'/

       jbuf=0
       ierr=0
       nstart=ist+n-1
       ifact=1
       do i = 1,n
        do j=1,14
          if(cchar(j).eq.irecord(nstart))goto 120
        enddo
110     ierr=1
        return
120     if(j.lt.11)goto 130
        if(nstart.ne.ist)goto 110
        if(j.ge.14)jbuf=-jbuf
        goto 150
130     jbuf=jbuf+(j-1)*ifact
        ifact=ifact*10
        nstart=nstart-1
       enddo

150    return

       END
