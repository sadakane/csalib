/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#ifndef _WT_H_
#define _WT_H_

#include "typedef.h"
//#include "comparray.h"
//#include "huffman.h"

typedef struct wavelet3 {
  i64 n; // length
  i64 sigma; // alphabet size
  int m; // number of distinct characters actually appeared in the string
  int k; // number of bytes in an integer
//  Huffman2 *huffman;
  void *huffman;
//  comparray **da;
  void *da;
  int L;
  i64 size;
  int opt;

//  MMAP *map;
  void *map;

  int (*access)(struct wavelet3 *wt, i64 i);
  i64 (*rank)(struct wavelet3 *wt, i64 i, int c);
  i64 (*select)(struct wavelet3 *wt, i64 i, int c);
  i64 (*rank_access)(struct wavelet3 *wt, i64 i, int *c);
  int (*enumerate)(struct wavelet3 *wt, i64 l, i64 r, int *C, i64 *freq);
} wavelet3;

i64 wt3_makeindex(i64 k2, int sigma, char *fname, int opt);
void wt3_read(wavelet3 *wt, char *fname);

#endif // _WT_H_
