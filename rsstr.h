/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#ifndef _RSSTR_H_
#define _RSSTR_H_

#ifdef __SSE4_2__
#include <smmintrin.h>
#define POPCOUNT4(x) _mm_popcnt_u64(x)
#else
#define POPCOUNT4(x) popcount4(x)
#endif


//#ifndef _BITVEC_T_
//#define _BITVEC_T_
//typedef u64 bitvec_t;
//#endif
#define DD (sizeof(bitvec_t)*8)

#define logSIGMADNA 4
//#define SIGMADNA (1<<logSIGMADNA)

#define SIGMADNA 9 // acgtACGT$

//#define SB 32
//#define MB 256
#define logLB 16
#define LB (1<<logLB)


typedef struct rsstr {
  i64 n;
  int k; // number of bytes in an integer
  int l;
  int lmb;
  bitvec_t *BW;
  uchar *RL;
  u16 *RM;
  i64 *C;
  i64 size;
  MMAP *mapbwt;

  int (*access)(struct rsstr *wt, i64 i);
  i64 (*rank)(struct rsstr *wt, i64 i, int c);
  i64 (*select)(struct rsstr *wt, i64 i, int c);
  i64 (*succ)(struct rsstr *wt, i64 i, int c);
//  i64 (*rank_access)(struct rsstr *wt, i64 i, int *c);

} rsstr;

i64 rsstr_makeindex(char *fname, int L);
void rsstr_read(rsstr *sa, char *fname);

#endif // _RSSTR_H_
