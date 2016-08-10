/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#ifndef _LF_DNA2_H_
#define _LF_DNA2_H_

#ifdef __SSE4_2__
#include <smmintrin.h>
#define POPCOUNT3(x) _mm_popcnt_u64(x)
#else
#define POPCOUNT3(x) popcount(x)
#endif



typedef struct {
  i64 n;
  i64 last;
  int k; // number of bytes in an integer
  int lb, mb;
  bitvec_t *BW; /* pointer to the bit stream encoding psi */
  uchar *RL;
  u16 *RM;
  i64 *C;
  i64 psize;
  i64 id; // type of encoding
  MMAP *mapbwt,*mapidx;
} lf_dna2;

i64 lf_dna2_makeindex(CSA *csa, char *fname);
void lf_dna2_read(CSA *sa, char *fname);
//static int lf_dna_BW_sub(lf_dna *lf,i64 i);
//static int lf_dna_BW(CSA *csa,i64 i);
//static i64 lf_dna_rankc_sub(lf_dna *lf, i64 i, int c);
//static i64 lf_dna_rankc(CSA *csa, i64 i, int c);
void lf_dna2_options(CSA *csa, char *p);



#endif // _LF_DNA2_H_
