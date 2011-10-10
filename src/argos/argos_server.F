      subroutine argos_server
c
      implicit none
c
#include "argos.fh"
c
#if !defined(WIN32)
      integer*4 create_client_socket
      integer client_socket_write
      external create_client_socket,client_socket_write
c
      integer*4 iiport,lens
c
      character*255 string
      integer numbyt
c
      if(iport.le.0.or.me.gt.0) return
c
c     open socket to server
c
      iiport=iport
      if(.not.lserver) then
      write(*,1) server(1:index(server,' ')-1),iiport
    1 format('Attempt to open socket to ',a,' port ',i5)
      isocket=create_client_socket(server,iiport)
      lserver=isocket.gt.0
      endif
c
      if(lserver) then
      print*,'server socket open'
      else
      print*,'server socket error'
      endif
c
      if(.not.lserver) return
c
      write(string,1000) stime,etot,temp,pres*1.0125d-5,volume
 1000 format("tETpV",4f12.3,f12.6)
c
      lens=66
      string(66:66)=char(13)
c
      numbyt=client_socket_write(isocket,string,lens)
c
      print*,'Bytes written to socket is ',numbyt
c
#endif
      return
      end
c $Id$
