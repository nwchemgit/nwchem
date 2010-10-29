/*
 $Id$
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "hdbm.h"

int random(void);
int srandom(int);

#define INNER 1000

int datum_to_int(d)
     datum d;
{
  int tmp;
  (void) memcpy((void *) &tmp, d.dptr, sizeof(tmp));
  return (int) tmp;
 }

void insert_random_pair(hdbm h)
{
  int j, keylen, datalen;
  char *key, *data;
  datum k, v;

  keylen = random()%127 + 1;
  if (!(key = malloc((size_t) keylen))) exit(43);
  for (j=0; j<keylen; j++)
    key[j] = j & 127;
  datalen = random()%32768 + 1;
  if (!(data = malloc((size_t) datalen))) exit(43);
  for (j=0; j<datalen; j++)
    data[j] = j & 127;
  datum_wrap(key, keylen, &k);
  datum_wrap(data, datalen, &v);
  if (!hdbm_replace(h, k, v)) exit(13);
  free(key); free(data);

  (void) printf("inserted key=%d data=%d\n", keylen, datalen);

  if (!hdbm_file_flush(h)) exit(51);
}

void MixedLengthTest(const char *filename)
{
  hdbm h;
  int n = 10, i, ntime = 10, available;
  datum k, v;

  if (!hdbm_open(filename, 1, &h)) exit(1);

  printf("before initial insert\n");
  hdbm_print_usage();

  /* Insert with random size for keys and values */

  for (i=0; i<n; i++) {
    insert_random_pair(h);
  }

  printf("before random delete\n");
  hdbm_print_usage();

  /* Loop thru database multiple times randomly deleting and 
     writing keys */

  while (ntime--) {
    for (available=hdbm_first_key(h, &k);
	 available;
	 available=hdbm_next_key(h, &k)) {
      if ( (random() & 1) ) {
	if (!hdbm_delete(h, k)) exit (31);
	insert_random_pair(h);
      }
      datum_free(k);
    }
  }

  printf("Before compress \n\n");
  if (!hdbm_print_table(h, 0, 0)) exit(77);
    
  printf("before close\n");
  hdbm_print_usage();

  if (!hdbm_close(h)) exit(2);

  printf("after close\n");
  hdbm_print_usage();

  if (!hdbm_file_compress(filename)) exit(3);

  if (!hdbm_open(filename, 1, &h)) exit(1);
  
  printf("After compress \n\n");
  if (!hdbm_print_table(h, 0, 0)) exit(77);

  if (!hdbm_close(h)) exit(2);

  hdbm_print_usage();
}

void Test1(int argc, char **argv)
{
  int i, outer;
  int use_random = 1;
  hdbm h, h2;

  if (!hdbm_open("test.db", 1, &h)) exit(1);

  if (argc != 1) use_random = 0;

#define INDEX(a) (use_random ? (random()&0x0f) : (a))

/*   PrintTableStats(); */

  for (outer=0; outer<200; outer++) {
    int ilo = outer*INNER;
    int ihi = ilo + INNER;

    srandom((int) (outer + 5));

    for (i=ilo; i<ihi; i++) {
      int index = INDEX(i);
      int i3 = index*3;
      datum key;
      datum value;
      if (!datum_make((char *) &i3, sizeof i3, &value))
	exit(1);
      if (!datum_make((char *) &index, sizeof index, &key))
	exit(1);
      if (!hdbm_insert(h, key, value))
	exit(1);
      datum_free(value);
      datum_free(key);
    }
/*     PrintTableStats(); */
  }
  
  hdbm_print_stats(h);

  hdbm_close(h);
  if (!hdbm_file_copy("test.db", "test2.db"))
    exit(1);
  if (!(hdbm_open("test.db", 1, &h) && hdbm_open("test2.db", 1, &h2)))
    exit(1);

  for (outer=190; outer>=0; outer--) {
    int ilo = outer*INNER;
    int ihi = ilo + INNER;

    srandom((int) (outer + 5));

    for (i=ilo; i<ihi; i++) {
      int index = INDEX(i);
      datum key;
      datum value;
      int i3;
      if (!datum_make((char *) &index, sizeof index, &key))
	exit(1);
      if (!hdbm_extract(h, key, &value))
	exit(1);
      if (!hdbm_delete(h2, key))
	exit(1);
      i3 = datum_to_int(value);
      if (i3 != (index*3)) {
	(void) printf("i=%d, index=%d, i3=%d\n", i, index, i3);
	exit(1);
      }
      else if ((i%1000) == 0)
	printf("%d, ", i);
      datum_free(value);
      datum_free(key);
    }
    printf("\n");
/*     PrintTableStats(); */
  }

  hdbm_close(h);
  hdbm_close(h2);
  if (!hdbm_file_copy("test.db", "test2.db")) /* test2 should be empty */
    exit(1);
}

/*ARGSUSED*/
int main(argc, argv)
     int argc;
     char **argv;
{
  /* Test1(argc, argv); */

  MixedLengthTest("test.db");

  return 0;
}
