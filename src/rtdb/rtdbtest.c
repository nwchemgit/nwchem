/*$Id$*/
#include <stdio.h>
#include "rtdb.h"
#include "misc.h"
#include "macdecls.h"
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[])
{
  int rtdb;
  char buf[4];
  int ibuf;
  double dbuf;

  PBEGIN_(argc, argv);
  
  (void) MA_init(MT_CHAR, -1, -1);
  (void) MA_set_auto_verify(1);
  
  if (!rtdb_open("test.db", "unknown", &rtdb))
    error("rtdbtest: open failed on %s\n", "test.db");

  /* Check character */

  if (!rtdb_put(rtdb, "hello", MT_CHAR, 4, "bob"))
    error("rtdbtest: c put failed\n", 0);

  if (!rtdb_get(rtdb, "hello", MT_CHAR, 4, buf))
    error("rtdbtest: c get failed\n", 0);

  printf("buf=%s\n", buf);

  /* Check integer */

  ibuf = 1;
  if (!rtdb_put(rtdb, "ibuf", MT_INT, 1, &ibuf))
    error("rtdbtest: i put failed\n", 0);

  ibuf = 0;
  if (!rtdb_get(rtdb, "ibuf", MT_INT, 1, &ibuf))
    error("rtdbtest: i get failed\n", 0);

  printf("ibuf=%d\n", ibuf);


  /* Check double */

  dbuf = 1.0;
  if (!rtdb_put(rtdb, "dbuf", MT_DBL, 1, &dbuf))
    error("rtdbtest: d put failed\n", 0);

  dbuf = 0.0;
  if (!rtdb_get(rtdb, "dbuf", MT_DBL, 1, &dbuf))
    error("rtdbtest: d get failed\n", 0);

  printf("dbuf=%f\n", dbuf);

  /* Check the info on all the above */

  {
    static const char *names[] = {"hello", "ibuf", "dbuf"};
    int n = 3;
    
    while (n--) {
      int ma_type, nelem;
      char date[26];
      
      if (!rtdb_get_info(rtdb, names[n], &ma_type, &nelem, date))
	error("get info failed n=%d\n", n);
      
      printf("\"%s\" is %d elements of type %d, dated %s",
	     names[n], nelem, ma_type, date);
    }
  }
  
  /* Now loop all processes thru the entire db doing parallel ma gets */
  
  {
    char name[256];
    int available;
    
    for (available=rtdb_first(rtdb, sizeof name, name);
	 available;
	 available=rtdb_next(rtdb, sizeof name, name)) {
      int ma_type, nelem, ma_handle;
      void *data;
      
      if (!rtdb_ma_get(rtdb, name, &ma_type, &nelem, &ma_handle))
	error("rtdbtest: ma_get failed\n", 0);

      if (!MA_get_pointer(ma_handle, &data))
	error("rtdbtest: MA_get_pointer failed\n", 0);
      
      (void) printf("%2ld: name=%s, t=%d, n=%d, ",
		    NODEID_(), name, ma_type, nelem);
      ma_print(stdout, ma_type, nelem, data);
      (void) fflush(stdout);
      (void) MA_free_heap(ma_handle);
    }
  }
  
  /* Now print the entire data base out */
  
  rtdb_print(rtdb, 1);
  
  if (!rtdb_close(rtdb, "delete"))
    error("rtdbtest: close failed\n", 0);
  
  pend_();
  return 0;
}
