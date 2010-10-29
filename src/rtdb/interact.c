/*$Id$*/
#include <stdio.h>
#include "macdecls.h"
#include "rtdb.h"

extern const char *ma_typename(int);

#define MAX_NELEM 1024

int ReadNelem()
{
  int nelem = -1;
  
  printf("Enter nelem -> "); fflush(stdout);
  if (scanf("%d", &nelem) != 1)
    exit(1);
  if (nelem > MAX_NELEM) {
    printf("too many elements\n");
    nelem = -1;
  }
  return nelem;
}

int read_string(char *a, int len)
{
  int n;

  while (1) {
    int c = getc(stdin);
    if (c == EOF)
      return 0;
    else if (!(c==' ' || c=='\t' || c=='\n')) {
      ungetc(c,stdin);
      break;
    }
  }
  for (n=0; n<(len-1); n++) {
    int c = getc(stdin);
    if (c == EOF)
      return 0; 
    else if (c == '\n') {
      a[n] = 0;
      return 1;
    }
    else
      a[n] = c;
  }
  return 0;
}

int main(int argc, char *argv[])
{
  int rtdb = -1, ma_type, nelem;
  char filename[256], mode[256], name[1024], date[26];
  void *data;

  char carray[MAX_NELEM+1];
  int iarray[MAX_NELEM];
  double darray[MAX_NELEM];
  
#define N_OPTS 9
  const char *opts[N_OPTS] = {"quit", "open", "close", "info", "put", "get", 
			  "first", "next", "print"};

  if (!MA_init(MT_DBL, -1, -1))
    exit(1);

  while (1) {
     int opt;

     printf("\n\n\n      Interactive RTDB\n      ----------------\n\n");
     for (opt=0; opt<N_OPTS; opt++)
       printf("        %-8s %3d\n", opts[opt], opt);

     printf("\nEnter option number -> "); fflush(stdout);
     if (scanf("%d", &opt) != 1) 
       break;

     switch (opt) {
     case 0:
       exit(0); break;

     case 1:
       printf("\nOpen database\n-----------\n\n");
       printf("Enter filename -> "); fflush(stdout);
       if (scanf("%s",filename) != 1) 
	 exit(1);
       printf("Enter mode (new, old, unknown, empty, scratch) -> ");
       fflush(stdout);
       if (scanf("%s",mode) != 1) 
	 exit(1);
       if (!rtdb_seq_open(filename, mode, &rtdb))
	 printf("\nOpen of %s with mode %s failed\n", filename, mode);
       break;

     case 2:
       printf("\nClose database\n-------------\n\n"); fflush(stdout);
       if (rtdb < 0) 
	 printf("database is not open\n");
       else {
	 printf("Enter mode (keep, delete) -> "); fflush(stdout);
	 if (scanf("%s",mode) != 1) 
	   exit(1);
	 if (!rtdb_seq_close(rtdb, mode))
	   printf("Close of %s with mode %s failed\n", filename, mode);
	 else
	   rtdb = -1;
       }
       break;

     case 3:
       printf("\nInformation on entry\n--------------------\n\n");
       if (rtdb < 0) 
	 printf("database is not open\n");
       else {
	 printf("Enter name -> "); fflush(stdout);
  	 if (!read_string(name, sizeof name))
	   exit(1);

	 if (!rtdb_seq_get_info(rtdb, name, &ma_type, &nelem, date))
	   printf("Get info on \"%s\" failed\n", name);
	 else
	   printf("%s -> type=%s, nelem=%d, date=%s\n", 
		  name, ma_typename(ma_type), nelem, date);
       }
       break;

     case 4:
       printf("\nPut entry\n---------\n\n");
       if (rtdb < 0) 
	 printf("database is not open\n");
       else {
	 printf("Enter name -> "); fflush(stdout);
	 if (!read_string(name, sizeof name))
	   exit(1);
	 printf("Enter type (int, char, double) -> "); fflush(stdout);
	 if (scanf("%s", mode) != 1)
	   exit(1);
	 if (!strcmp(mode,"int")) {
	   int i;
	   if ((nelem = ReadNelem()) <= 0) break;
	   ma_type = MT_INT;
	   for (i=0; i<nelem; i++)
	     if (scanf("%d", iarray+i) != 1)
	       exit(1);
	   data = (void *) iarray;
	 }
	 else if (!strcmp(mode, "char")) {
	   int i;
	   ma_type = MT_CHAR;
	   if (!read_string(carray, MAX_NELEM)) break;
	   nelem = strlen(carray) + 1;
	   data = (void *) carray;
	 }
	 else if (!strcmp(mode, "double")) {
	   int i;
	   if ((nelem = ReadNelem()) <= 0) break;
	   ma_type = MT_DBL;
	   for (i=0; i<nelem; i++)
	     if (scanf("%lf", darray+i) != 1)
	       exit(1);
	   data = (void *) darray;
	 }
	 else {
	   printf("invalid type\n");
	   break;
	 }
	 if(!rtdb_seq_put(rtdb, name, ma_type, nelem, data))
	   printf("put %s, nelem=%d, type=%s failed", name, 
		  nelem, ma_typename(ma_type)); fflush(stdout);
       }
       break;


     case 5:
       printf("\nGet entry\n---------\n\n");
       if (rtdb < 0) 
	 printf("database is not open\n");
       else {
	 int ma_handle;

	 if(rtdb_seq_ma_get(rtdb, name, &ma_type, &nelem, &ma_handle)) {
	   ma_print(stdout, ma_type, nelem, MA_get_pointer(ma_handle, &data));
	   MA_free_heap(ma_handle);
	 }
	 else 
	   printf("Get of %s failed\n", name);
       }
       break;

     case 8:
       printf("\nPrint database\n--------------\n\n");
       if (rtdb < 0) 
	 printf("database is not open\n");
       else {
	 int values;
	 printf("Print values (0=no, 1=yes) -> "); fflush(stdout);
	 if (scanf("%d", &values) != 1)
	   exit(1);
	 rtdb_seq_print(rtdb, values);
       }
       break;
       
     default:
       break;
     }
   }
  return 0;
}


