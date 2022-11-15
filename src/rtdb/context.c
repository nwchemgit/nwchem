/*$Id$*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rtdb.h"
#include "macdecls.h"
#include "misc.h"

#define MAX_CLEN 4096
static char context[MAX_CLEN];

int context_set(const char *string)
{
  if (strlen(string) < sizeof(context)) {
    (void) strcpy(context, string);
    return 1;
  }
  else {
    fprintf(stderr, "context_set: string too long? %s\n", string);
    fflush(stderr);
    return 0;
  }
}

char *context_get(void)
{
  return strdup(context);
}

int context_rtdb_store(int rtdb)
{
  return rtdb_put(rtdb, "Context", MT_CHAR, strlen(context)+1, context);
}

int context_rtdb_load(int rtdb)
{
  return rtdb_get(rtdb, "Context", MT_CHAR, sizeof(context), context);
}

int context_push(const char *string)
{
  int clen = strlen(context);
  int slen = strlen(string);
  
  if (slen+clen+2 >= sizeof(context)) {
    fprintf(stderr, "context_push: static dimension of context too small\n");
    fprintf(stderr, "context_push: current = %s\n", context);
    fprintf(stderr, "context_push: pushing = %s\n", string);
    return 0;
  }
  else {
    (void) strcpy(context+clen, string);
    (void) strcpy(context+clen+slen, ":");
    return 1;
  }
}

int context_pop(const char *string)
{
  int clen = strlen(context);
  int slen = strlen(string);
  
  if (clen)
    clen--;			/* Trailing colon */
  
  if (slen <= clen && strncmp(context+clen-slen, string, slen) == 0) {
    context[clen-slen] = 0;
    return 1;
  }
  else {
    fprintf(stderr, "context_pop: current = %s\n", context);
    fprintf(stderr, "context_pop: popping = %s\n", string);
    return 0;
  }
}

int context_rtdb_match(int rtdb, const char *name, int reslen,
		       char *result)
{
  char buf[MAX_CLEN];
  int blen = strlen(context);
  
  if (blen+strlen(name)+1 > sizeof(buf)) {
    fprintf(stderr, "context_rtdb_match: buffer size exceeded\n");
    fprintf(stderr, "context_rtdb_match: current = %s\n", context);
    fprintf(stderr, "context_rtdb_match: pushing = %s\n", name);
    return 0;
  }
  
  strcpy(buf, context);
  
  while (1) {
    int ma_type, nelem;
    char date[26];
    
    /* Append name to current context */
    
    (void) strcpy(buf+blen, name);
    
    if (rtdb_get_info(rtdb, buf, &ma_type, &nelem, date)) {
      if (ma_type == MT_CHAR) {
	if (!rtdb_get(rtdb, buf, ma_type, reslen, result)) {
	  fprintf(stderr, "context_rtdb_match: rtdb_get failed?\n");
	  return 0;
	}
	reslen = strlen(result);
	if (result[reslen-1] == '\n') /* Fortran cput appends an unwanted CR */
	  result[reslen-1] = 0;
	return 1;
      }
      else {
	fprintf(stderr, "context_rtdb_match: found %s but is wrong type\n",
		name);
	return 0;
      }
    }
    else {
      
      /* Did not find entry ... pop the context stack */
      
      if (!blen)
	return 0;		/* Stack is alredy empty */
      
      blen--;
      while (--blen > 0)
	if (buf[blen] == ':')
	  break;
    }
  }
  
  return 1;			/* Never executed */
}



int context_prefix(const char *name, char *result, int result_len)
{
  if ((strlen(name)+strlen(context)+1) > result_len) {
    fprintf(stderr, "constant_prefix: result too short\n");
    return 0;
  }
  strcpy(result,context);
  strcpy(result+strlen(context),name);
  
  return 1;
}

/*
static void context_print()
{
  printf("context = -%s-\n", context);
}
int main()
{
  int rtdb;
  char *cntx;

  (void) MA_initialize(MT_CHAR, -1, -1);

  if (!rtdb_open("test.db", "unknown", &rtdb))
    error("testcontext: open failed on %s\n", "test.db");

  context_print();
  if (!context_push("optimize"))
    error("context push failed %d\n", 0);
    context_print();
  if (!context_push("scf"))
    error("context push failed %d\n", 0);
    context_print();
  if (!context_push("rhf"))
    error("context push failed %d\n", 0);
    context_print();
  if (!context_push("pcg"))
    error("context push failed %d\n", 0);
    context_print();

  (void) context_store(rtdb);

  (void) context_set("");

  (void) context_print();

  if (!context_load(rtdb))
    error("context_load: failed %d\n", 0);

  (void) context_print();

  cntx = context_get();
  printf("context from get = %s\n", cntx);

  if (context_pop("scf"))
    error("context pop succeeded %d\n", 0);
  if (!context_pop("pcg"))
    error("context pop failed %d\n", 0);
  context_print();
  if (!context_pop("rhf"))
    error("context pop failed %d\n", 0);
    context_print();
  if (!context_pop("scf"))
    error("context pop failed %d\n", 0);
    context_print();
  if (!context_pop("optimize"))
    error("context pop failed %d\n", 0);
    context_print();

  (void) rtdb_close(rtdb, "delete");

  return 0;
}
*/
