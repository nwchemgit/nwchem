/*
 $Id$
 */

#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <string.h>
#include "hdbm.h"

void uppercase(char *s)
{
  int c;

  while ((c = *s++))
    if (c >= 'a' && c <= 'z')
      *(s-1) = c + 'A' - 'a';
}

void lowercase(char *s)
{
  int c;

  while ((c = *s++))
    if (c >= 'A' && c <= 'Z')
      *(s-1) = c + 'a' - 'A';
}

int read_and_insert_word(FILE *file, hdbm db)
/*
  Read a word from the file and insert it into the database
  with the key being lowercase and the value uppercase.
*/
{
  int status;
  char word[80];
  datum key, value;

  if (fscanf(file, "%s", word) != 1)
    return 0;
  
  lowercase(word);
  if (!datum_make(word, strlen(word)+1, &key))
    return 0;

  uppercase(word);
  if (!datum_make(word, strlen(word)+1, &value))
    return 0;

  status = hdbm_insert(db, key, value);

  datum_free(key); datum_free(value);

  return status;
}

int main()
{
  FILE *file = fopen("/usr/dict/words", "r");
  hdbm db;
  int nentry = 0;
  char word[80];

  if (!file)
    exit(1);
  if (!hdbm_open("/tmp/words.db", 1, &db))
    exit(1);

  while (read_and_insert_word(file, db) && nentry < 20000)
    nentry++;

  printf("Inserted %d words\n", nentry);

  hdbm_print_stats(db);
  hdbm_print_usage();

  (void) hdbm_close(db);

  if (!hdbm_open("/tmp/words.db", 1, &db))
    exit(1);

  printf(" enter :");
  while (scanf("%s",word) == 1) {
      datum key;

      lowercase(word);

      datum_wrap(word, strlen(word)+1, &key);

      if (hdbm_probe(db, key)) 
	  printf("%s is found\n", word);
      else
	  printf("%s is not found\n", word);
      printf(" enter :");
  }

  (void) hdbm_close(db);

  (void) unlink("/tmp/words.db");

  return 0;
}





  
