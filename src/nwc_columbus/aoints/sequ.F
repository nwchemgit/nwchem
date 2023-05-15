      integer function indx(a,b,c,d,lena,lenb,lenc,lend,s,t,u,v,dbg,
     &  iimn, jjmn, kkmn, llmn)
c     safe but inefficient brute force shell index computer
c     replace with direct computation of index
      integer a, b, c, d
      integer lena, lenb, lenc, lend
      integer s, t, u, v
      integer iimn, jjmn, kkmn, llmn
      integer ii, jj, kk, ll
      integer i, j, k, l, imn, imx, jmn, jmx, lmx, kmx
      logical dbg
      integer ij, kl
      integer st, uv

*     write(26,*) 'in indx'

      st=(s*(s-1))/2+t
      uv=(u*(u-1))/2+v

      if(dbg) then
        write(26,100)'calculating local shell quartet index for ',
     &    a,b,c,d
        write(26,100)'   shell sizes: ', lena, lenb, lenc, lend
      endif
  100 format(a,4i2)

      if(st.lt.uv) then
        write(26,*) 'u == v and d > c'
        indx=-1
        goto 1000
      endif

      if(u.eq.v.and.d.gt.c) then
        write(26,*) 'u == v and d > c'
        indx=-1
        goto 1000
      endif
*     if(t.eq.v.and.d.gt.b) then
*       write(26,*) 't == v and d > b'
*       indx=-1
*       goto 1000
*     endif
*     if(s.eq.u.and.c.gt.a) then
*       write(26,*) 's == u and c > a'
*       indx=-1
*       goto 1000
*     endif
      if(s.eq.t.and.b.gt.a) then
        write(26,*) 's == t and b > a'
        indx=-1
        goto 1000
      endif

      if(t.gt.s) then
        write(26,*) 't > s'
        indx=-1
        goto 1000
      endif
   
      if(v.gt.u) then
        write(26,*) 'v > u'
        indx=-1
        goto 1000
      endif


      indx=0
*     i=1
*     j=1
*     k=1
*     l=1
      do k=1,c
        if (v.eq.u) then
*         if(k.lt.c.or.j.lt.b.or.i.lt.a) then
          if(k.lt.c) then
            lmx=k
          else
            lmx=d
          endif
        else
*         if(k.lt.c.or.j.lt.b.or.i.lt.a) then
          if(k.lt.c) then
            lmx=lend
          else
            lmx=d
          endif
        endif
        kk=k+kkmn-1
        do l=1,lmx
*         if (u.eq.s) then
*           imn=k
*         else
            imn=1
*         endif
          if(l.lt.d.or.k.lt.c) then 
            imx=lena
          else 
            imx=a
          endif
          ll=l+llmn-1
          kl=((kk-1)*kk)/2+ll
          do i=imn,imx
*           write(*,*)'i = ', i, ' b = ', b, ' lenb = ', lenb
            if(s.eq.t) then 
             if(l.lt.d.or.k.lt.c.or.i.lt.a) then
              jmx=i
             else
              jmx=b
             endif
            else
             if(l.lt.d.or.k.lt.c.or.i.lt.a) then
              jmx=lenb
             else
              jmx=b
             endif
            endif
            jmn=1
*           write(*,*)'jmn=',jmn,' jmx=',jmx
            ii=iimn+i-1
            do j=jmn,jmx
              jj=jjmn+j-1
              ij=((ii-1)*ii)/2+jj
*       write(26,2323) '[',i,',',j,';',k,',',l,']'
 2323   format(4(a,i2),a1)
*       write(26,'(2(a,i3))') '(ij)=',ij,' (kl)=',kl

              if(st.eq.uv) then
               if(ij.ge.kl) then
                indx=indx+1
                if(dbg) then
                 write(26,'(4i2,a,i3)') i,j,k,l, ' - ', indx
                endif
               endif
              else
               indx=indx+1
               if(dbg) then
                write(26,'(4i2,a,i3)') i,j,k,l, ' - ', indx
               endif
              endif

            enddo
          enddo
        enddo
      enddo
      if(dbg) then
        write(26,'(a,i4)') 'local index=', indx
      endif
*     write(26,*) 'leaving indx'
 1000 continue

      return
      end
     

