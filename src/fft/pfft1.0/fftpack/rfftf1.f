      subroutine rfftf1 (n,c,ch,wa,ifac)
      dimension       ch(*)      ,c(*)       ,wa(*)      ,ifac(*)
      nf = ifac(2)
      na = 1
      l2 = n
      iw = n
      do 111 k1=1,nf
         kh = nf-k1
         ip = ifac(kh+3)
         l1 = l2/ip
         ido = n/l2
         idl1 = ido*l1
         iw = iw-(ip-1)*ido
         na = 1-na
         if (ip .ne. 4) go to 102
         ix2 = iw+ido
         ix3 = ix2+ido
         if (na .ne. 0) go to 101
         call radf4 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
         go to 110
  101    call radf4 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
         go to 110
  102    if (ip .ne. 2) go to 104
         if (na .ne. 0) go to 103
         call radf2 (ido,l1,c,ch,wa(iw))
         go to 110
  103    call radf2 (ido,l1,ch,c,wa(iw))
         go to 110
  104    if (ip .ne. 3) go to 106
         ix2 = iw+ido
         if (na .ne. 0) go to 105
         call radf3 (ido,l1,c,ch,wa(iw),wa(ix2))
         go to 110
  105    call radf3 (ido,l1,ch,c,wa(iw),wa(ix2))
         go to 110
  106    if (ip .ne. 5) go to 108
         ix2 = iw+ido
         ix3 = ix2+ido
         ix4 = ix3+ido
         if (na .ne. 0) go to 107
         call radf5 (ido,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 110
  107    call radf5 (ido,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
         go to 110
  108    if (ido .eq. 1) na = 1-na
         if (na .ne. 0) go to 109
         call radfg (ido,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
         na = 1
         go to 110
  109    call radfg (ido,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
         na = 0
  110    l2 = l1
  111 continue
      if (na .eq. 1) return
      do 112 i=1,n
         c(i) = ch(i)
  112 continue
      return
      end
