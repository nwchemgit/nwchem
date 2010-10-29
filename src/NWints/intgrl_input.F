      subroutine intgrl_input(rtdb)
      implicit none
* $Id$
#include "stdio.fh"
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
*     if (rtdb_get(rtdb,'int:cando_txs',MT_LOG,1,user_status)) then
*     if (rtdb_get(rtdb,'int:cando_sp',MT_LOG,1,user_status)) then
*     if (rtdb_get(rtdb,'int:cando_hnd',MT_LOG,1,user_status)) then
*     if (.not.(rtdb_get(rtdb,'int:full_mem',Mt_Log,1,o_do_full)))
*     if (.not. rtdb_get(rtdb,'int:mem2emax',Mt_Int,1,mem2emax))
*     if (.not.(rtdb_get(rtdb,'int:full_mem',Mt_Log,1,o_do_full)))
*     if (rtdb_get(rtdb,'int:acc_std',mt_dbl,1,val_int_acc_std)) then
*     if (rtdb_get(rtdb,'int:acc_high',mt_dbl,1,val_int_acc_high)) then 
*     if (.not. rtdb_get(rtdb,'dyall_mod_dir',mt_log,1,dyall_mod_dir))
*     if (.not. rtdb_get(rtdb,'atomic_nesc',mt_log,1,atomic_nesc))
*     if (.not. rtdb_get(rtdb,'nesc_1e_approx',mt_log,1,nesc_1e_approx))
*     if (.not. rtdb_get(rtdb,'ss_one_cent',mt_log,1,ss_one_cent))
*     if (.not. rtdb_get(rtdb,'do_ssss',mt_log,1,do_ssss)) 
*     if (.not. rtdb_get(rtdb,'doug_kroll',MT_LOG,1,doug_kroll)) 
*     if (.not. rtdb_get(rtdb,'doug_kroll:type',MT_INT,1,
*     if (.not. rtdb_get(rtdb,'doug_kroll:eventemp',MT_LOG,1,
*     if (.not. rtdb_get(rtdb,'doug_kroll:etratio',MT_DBL,1,
*     if (rtdb_get(rtdb,'int:cando_sp',MT_LOG,1,status)) then
*     if (rtdb_get(rtdb,'int:cando_nw',MT_LOG,1,status)) then
*     if (rtdb_get(rtdb,'int:cando_txs',MT_LOG,1,status)) then
*     if (rtdb_get(rtdb,'int:cando_hnd',MT_LOG,1,status)) then
*     if (.not. rtdb_get(rtdb,'dyall_mod_dir',mt_log,1,dyall_mod_dir))
*     if (.not. rtdb_get(rtdb,'atomic_nesc',mt_log,1,atomic_nesc))
*     if (.not. rtdb_get(rtdb,'nesc_1e_approx',mt_log,1,nesc_1e_approx))
*     if (.not. rtdb_get(rtdb,'ss_one_cent',mt_log,1,ss_one_cent))
*     if (.not. rtdb_get(rtdb,'do_ssss',mt_log,1,do_ssss)) 
*     if (.not. rtdb_get(rtdb,'doug_kroll',MT_LOG,1,doug_kroll))
*     if (.not. rtdb_get(rtdb,'doug_kroll:type',MT_INT,1,
*     if (.not. rtdb_get(rtdb,'doug_kroll:eventemp',MT_LOG,1,
*     if (.not. rtdb_get(rtdb,'doug_kroll:etratio',MT_DBL,1,
*     if (rtdb_get(rtdb,'intd:cando_txs',MT_LOG,1,status))
*     if (rtdb_get(rtdb,'int:cando_sp',MT_LOG,1,status)) then
*     if (rtdb_get(rtdb,'int:cando_nw',MT_LOG,1,status)) then
*     if (rtdb_get(rtdb,'int:cando_txs',MT_LOG,1,status)) then
*     if (rtdb_get(rtdb,'int:cando_hnd',MT_LOG,1,status)) then
*     if (rtdb_get(rtdb,'intd:approxmem',MT_INT,1,intd_memtmp)) then
*     if (rtdb_get(rtdb,'intdd:cando_txs',MT_LOG,1,status))
*     if (rtdb_get(rtdb,'int:cando_sp',MT_LOG,1,status)) then
*     if (rtdb_get(rtdb,'int:cando_nw',MT_LOG,1,status)) then
*     if (rtdb_get(rtdb,'int:cando_txs',MT_LOG,1,status)) then
*     if (rtdb_get(rtdb,'int:cando_hnd',MT_LOG,1,status)) then
*     if (rtdb_get(rtdb,'intdd:approxmem',MT_INT,1,intd_memtmp)) then
*     status = rtdb_get(rtdb, 'doug_kroll', mt_log, 1, dk_run)
*     status = rtdb_get(rtdb, 'dyall_mod_dir', mt_log, 1, dmd_run)
*     if (rtdb_get(rtdb,'int:cando_txs',MT_LOG,1,user_status)) then
*     if (rtdb_get(rtdb,'intd:cando_txs',MT_LOG,1,user_status)) then
*     if (rtdb_get(rtdb,'intdd:cando_txs',MT_LOG,1,user_status)) then
*     if (rtdb_get(rtdb,'int:cando_sp',MT_LOG,1,user_status)) then
*     if (rtdb_get(rtdb,'int:cando_hnd',MT_LOG,1,user_status)) then
*     if (rtdb_get(rtdb,'int:cando_nw',MT_LOG,1,user_status)) then
*emem if (.not.(rtdb_get(rtdb,'int:full_mem',Mt_Log,1,o_do_full)))
*     if (rtdb_get(rtdb,'intd:approxmem',MT_INT,1,intd_memtmp)) then
*     if (rtdb_get(rtdb,'intdd:approxmem',MT_INT,1,intd_memtmp)) then
*     if (rtdb_get(rtdb,'intcsum:all',MT_LOG,1,dummyL)) then
*     if (rtdb_get(rtdb,'intcsum:1e',MT_LOG,1,dummyL)) then
*     if (rtdb_get(rtdb,'intcsum:2e',MT_LOG,1,dummyL)) then
*     if (rtdb_get(rtdb,'intcsum:ov',MT_LOG,1,dummyL)) then
*     if (rtdb_get(rtdb,'intcsum:ke',MT_LOG,1,dummyL)) then
*     if (rtdb_get(rtdb,'intcsum:pe',MT_LOG,1,dummyL)) then
*     if (rtdb_get(rtdb,'intcsum:h1',MT_LOG,1,dummyL)) then
*     if (rtdb_get(rtdb,'intcsum:3ov',MT_LOG,1,dummyL)) then
*     if (rtdb_get(rtdb,'intcsum:b2e4c',MT_LOG,1,dummyL)) then
*     if (rtdb_get(rtdb,'intcsum:2e4c',MT_LOG,1,dummyL)) then
*     if (rtdb_get(rtdb,'intcsum:2e3c',MT_LOG,1,dummyL)) then
*     if (rtdb_get(rtdb,'intcsum:2e2c',MT_LOG,1,dummyL)) then
*     if (rtdb_get(rtdb,'intcsum:ovprint',MT_LOG,1,dumm69)) then
*     if (rtdb_get(rtdb,'intcsum:keprint',MT_LOG,1,dumm69)) then
*     if (rtdb_get(rtdb,'intcsum:peprint',MT_LOG,1,dumm69)) then
*     if (rtdb_get(rtdb,'intcsum:h1print',MT_LOG,1,dumm69)) then
*     if (rtdb_get(rtdb,'intcsum:3ovprint',MT_LOG,1,dumm69)) then
*     if (rtdb_get(rtdb,'intcsum:2e4cprint',MT_LOG,1,dumm69)) then
*     if (rtdb_get(rtdb,'intcsum:2e3cprint',MT_LOG,1,dumm69)) then
*     if (rtdb_get(rtdb,'intcsum:2e2cprint',MT_LOG,1,dumm69)) then
*     if (rtdb_get(rtdb,'intcsum:b2e4cprint',MT_LOG,1,dumm69)) then
*     if (rtdb_get(rtdb,'intcsum:b2e4c:all_at_once',MT_LOG,1,status)) then
*     if (rtdb_get(rtdb,'intcsum:b2e4c:atom_size',MT_INT,1,int_tmp)) then
*     if (rtdb_get(rtdb,'intcsum:b2e4c:quartet_size',MT_INT,1,int_tmp)) then
*     if (rtdb_get(rtdb,'intcsum:b2e4c:buffer_size',MT_INT,1,int_tmp)) then
*     if(rtdb_get(rtdb,'intcsum:b2e4c:scratch_size',MT_INT,1,int_tmp)) then
*     if(.not.rtdb_get(rtdb,'int:txs:iprint',MT_INT,1,iprint)) then
*     if (rtdb_get(rtdb,'int:txs:icheck',MT_INT,1,int_tmp)) then
*     if (rtdb_get(rtdb,'int:txs:ncache',MT_INT,1,int_tmp)) then
*     else if (rtdb_get(rtdb,'int:txs:ncac',MT_INT,1,int_tmp)) then
*     if (rtdb_get(rtdb,'int:txs:maxprice',MT_INT,1,int_tmp)) then
*     else if(rtdb_get(rtdb,'int:txs:ipay',MT_INT,1,int_tmp)) then
*     if (rtdb_get(rtdb,'int:txs:iprint',MT_INT,1,int_tmp)) then
*     else if(rtdb_get(rtdb,'int:txs:prin',MT_INT,1,int_tmp)) then
*     if (rtdb_get(rtdb,'int:txs:limi',MT_INT,3,rtdblim)) then
*     if (rtdb_get(rtdb,'int:txs:limxmem',MT_INT,1,int_tmp)) then
*     if (rtdb_get(rtdb,'int:txs:limblks',MT_INT,1,int_tmp)) then
*     if (rtdb_get(rtdb,'int:txs:limpair',MT_INT,1,int_tmp)) then
*     if (rtdb_get(rtdb,'int:txs:iroute',MT_INT,1,int_tmp)) then
*     else if(rtdb_get(rtdb,'int:txs:rout',MT_INT,1,int_tmp)) then
*     if(rtdb_get(rtdb,'int:txs:lost_limit',MT_INT,1,int_tmp)) then
*     if (rtdb_get(rtdb,'intcsum:all',MT_LOG,1,dummyL)) then
*     if (rtdb_get(rtdb,'intcsum:1e',MT_LOG,1,dummyL)) then
*     if (rtdb_get(rtdb,'intcsum:2e',MT_LOG,1,dummyL)) then
*     if (rtdb_get(rtdb,'intcsum:ov',MT_LOG,1,dummyL)) then
*     if (rtdb_get(rtdb,'intcsum:ke',MT_LOG,1,dummyL)) then
*     if (rtdb_get(rtdb,'intcsum:pe',MT_LOG,1,dummyL)) then
*     if (rtdb_get(rtdb,'intcsum:h1',MT_LOG,1,dummyL)) then
*     if (rtdb_get(rtdb,'intcsum:3ov',MT_LOG,1,dummyL)) then
*     if (rtdb_get(rtdb,'intcsum:b2e4c',MT_LOG,1,dummyL)) then
*     if (rtdb_get(rtdb,'intcsum:2e4c',MT_LOG,1,dummyL)) then
*     if (rtdb_get(rtdb,'intcsum:2e3c',MT_LOG,1,dummyL)) then
*     if (rtdb_get(rtdb,'intcsum:2e2c',MT_LOG,1,dummyL)) then
*rak: if (rtdb_get(rtdb,'int:ov69',MT_LOG,1,dumm69)) then
*rak: if (rtdb_get(rtdb,'int:ke69',MT_LOG,1,dumm69)) then
*rak: if (rtdb_get(rtdb,'int:pe69',MT_LOG,1,dumm69)) then
*rak: if (rtdb_get(rtdb,'int:h169',MT_LOG,1,dumm69)) then
*rak: if (rtdb_get(rtdb,'int:4c69',MT_LOG,1,dumm69)) then
*rak: if (rtdb_get(rtdb,'int:3c69',MT_LOG,1,dumm69)) then
*rak: if (rtdb_get(rtdb,'int:2c69',MT_LOG,1,dumm69)) then
*     if (rtdb_get(rtdb,'intcsum:ovprint',MT_LOG,1,dumm69)) then
*     if (rtdb_get(rtdb,'intcsum:keprint',MT_LOG,1,dumm69)) then
*     if (rtdb_get(rtdb,'intcsum:peprint',MT_LOG,1,dumm69)) then
*     if (rtdb_get(rtdb,'intcsum:h1print',MT_LOG,1,dumm69)) then
*     if (rtdb_get(rtdb,'intcsum:3ovprint',MT_LOG,1,dumm69)) then
*     if (rtdb_get(rtdb,'intcsum:2e4cprint',MT_LOG,1,dumm69)) then
*     if (rtdb_get(rtdb,'intcsum:2e3cprint',MT_LOG,1,dumm69)) then
*     if (rtdb_get(rtdb,'intcsum:2e2cprint',MT_LOG,1,dumm69)) then
*     if (.not. rtdb_get(rtdb,'dyall_mod_dir',mt_log,1,dyall_mod_dir))
*     if (.not. rtdb_get(rtdb,'nesc_1e_approx',mt_log,1,nesc_1e_approx))
*     if (.not. rtdb_get(rtdb,'int:mem2emax',Mt_Int,1,mem2c))
*     if (rtdb_get(rtdb,'intcsum:b2e4cprint',MT_LOG,1,dumm69)) then
*     if (rtdb_get(rtdb,'intcsum:b2e4c:all_at_once',
*     if(rtdb_get(rtdb,'intcsum:b2e4c:atom_size',
*     if(rtdb_get(rtdb,'intcsum:b2e4c:quartet_size',
*     if(rtdb_get(rtdb,'intcsum:b2e4c:buffer_size',
*     if(rtdb_get(rtdb,'intcsum:b2e4c:scratch_size',
*     if (.not. rtdb_get(rtdb,'atomic_nesc',mt_log,1,atomic_nesc))
*     if (.not. rtdb_get(rtdb,'nesc_1e_approx',mt_log,1,one_e_approx))
*     if (.not. rtdb_get(rtdb,'ss_one_cent',mt_log,1,one_cent)) 
*     if (.not. rtdb_get(rtdb,'do_ssss',mt_log,1,do_ssss)) 
*     if (.not. rtdb_get(rtdb,'atomic_nesc',mt_log,1,atomic_nesc))
*     if (.not. rtdb_get(rtdb,'nesc_1e_approx',mt_log,1,nesc_1e_approx))
*     if (rtdb_get(rtdb,'relativistic:clight',mt_dbl,1,clight))
*     if(.not.rtdb_get(rtdb,'int:txs:iprint',MT_INT,1,iprint)) then
*     if (rtdb_get(rtdb,'int:txs:icheck',MT_INT,1,int_tmp)) then
*     if (rtdb_get(rtdb,'int:txs:ncache',MT_INT,1,int_tmp)) then
*     else if (rtdb_get(rtdb,'int:txs:ncac',MT_INT,1,int_tmp)) then
*     if (rtdb_get(rtdb,'int:txs:maxprice',MT_INT,1,int_tmp)) then
*     else if(rtdb_get(rtdb,'int:txs:ipay',MT_INT,1,int_tmp)) then
*     if (rtdb_get(rtdb,'int:txs:iprint',MT_INT,1,int_tmp)) then
*     else if(rtdb_get(rtdb,'int:txs:prin',MT_INT,1,int_tmp)) then
*     if (rtdb_get(rtdb,'int:txs:limi',MT_INT,3,rtdblim)) then
*     if (rtdb_get(rtdb,'int:txs:limxmem',MT_INT,1,int_tmp)) then
*     if (rtdb_get(rtdb,'int:txs:limblks',MT_INT,1,int_tmp)) then
*     if (rtdb_get(rtdb,'int:txs:limpair',MT_INT,1,int_tmp)) then
*     if (rtdb_get(rtdb,'int:txs:iroute',MT_INT,1,int_tmp)) then
*     else if(rtdb_get(rtdb,'int:txs:rout',MT_INT,1,int_tmp)) then
*     if(rtdb_get(rtdb,'int:txs:lost_limit',MT_INT,1,int_tmp)) then
* int
*   [sp|rotated_axis] [on|off]
*   hondo [on|off]
*   [nwchem|mmdav] [on|off] 
*   texas [on|off] [1stoff] [2ndoff] [iprint <integer>] [icheck <integer>] \
*         [ncache <integer>] [route <integer>]\
*         [limit [xmem <integer>] [blocks <integer>] [pairs <integer>] lost <integer>]
*   accuracy [<real> | tol <integer>] [sleazy|default|tight|maximum] 
*   memory [approximate|full] [1stmem <integer>] [2ndmem <integer>]
*   checksum [all|1e|2e] || [ov|ke|pe|h1|3ov|2e4c|2e3c|2e2c|b2e4c] \
*            [[buffer|scratch] <integer>] [[quartet|atom] <integer>] 
*   print [ov|ke|pe|h1|3ov|2e4c|2e3c|2e2c|b2e4c] 
*   profile|timing_stats [on|off]
* end
      integer rtdb
*
      integer flen
      character*255 field
*
      if (ga_nodeid().ne.0) return
      call inp_set_field(0) ! goto start of the line
      if (.not. inp_a(field)) call errquit
     &    ('intgrl_input: no input present',911, UNKNOWN_ERR)
      if (.not.(inp_compare(.false.,'intgrl',field).or.
     &          inp_compare(.false.,'integral',field))) call errquit
     &    ('intgrl_input: not INTeGRaL input block',911, UNKNOWN_ERR)
00100 continue
      if (inp_read()) then
        if (.not.inp_a(field)) call errquit
     &      ('intgrl_input:any-keyword could not read field',911,
     &       UNKNOWN_ERR)
        if (inp_compare(.false.,'end',field)) then
          goto 90001
        else if (inp_compare(.false.,'rotated_axis',field).or.
     &           inp_compare(.false.,'sp',field)) then
          if (.not.inp_a(field)) call errquit
     &      ('intgrl_input:sp|rotated_axis: could not read on|off',911, 
     &      UNKNOWN_ERR)
          if (inp_compare(.false.,'on',field)) then
            if (.not.rtdb_put(rtdb,'int:cando_sp',MT_LOG,1,.true.))
     &          call errquit('intgrl_input: rtdb_put failed',911,
     &       RTDB_ERR)
          else if (inp_compare(.false.,'off',field)) then
            if (.not.rtdb_put(rtdb,'int:cando_sp',MT_LOG,1,.false.))
     &          call errquit('intgrl_input: rtdb_put failed',911,
     &       RTDB_ERR)
          else
            flen = inp_strlen(field)
            write(luout,*)' expected on|off got: ',field(1:flen)
            call errquit('intgrl_input: fatal error ',911, UNKNOWN_ERR)
          endif
          goto 00100
        else 
          goto 00100
        endif
      endif
90001 continue
      end
