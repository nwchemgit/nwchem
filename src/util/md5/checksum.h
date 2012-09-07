/*
 $Id$
 */
#ifndef CHECKSUM_H
#define CHECkSUM_H

/**
\brief Initialize the checksum
*/

extern void checksum_init(void);

/**
\brief Include new data in the checksum
*/

extern void checksum_update(int len, const void *buf);

/**
\brief Retrieve the final checksum
*/

extern void checksum_final(char sum[33]);

/**
\brief Calculate the checksum of a single buffer
*/

extern void checksum_simple(int len, const void *buf, char sum[33]);

#endif
