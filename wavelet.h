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

typedef struct wavelet {
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

  int (*access)(struct wavelet *wt, i64 i);
  i64 (*rank)(struct wavelet *wt, i64 i, int c);
  i64 (*select)(struct wavelet *wt, i64 i, int c);
  i64 (*rank_access)(struct wavelet *wt, i64 i, int *c);
  int (*enumerate)(struct wavelet *wt, i64 l, i64 r, int *C, i64 *freq);
} wavelet;

i64 wt_makeindex(i64 k2, int sigma, char *fname, int opt);
void wt_read(wavelet *wt, char *fname);

#endif // _WT_H_
