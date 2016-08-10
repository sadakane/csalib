/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#ifndef _DENSEARRAY4_H_
#define _DENSEARRAY4_H_

#include "typedef.h"
#include "mman.h"

#ifdef __SSE4_2__
#include <smmintrin.h>
#define POPCOUNT(x) _mm_popcnt_u64(x)
#else
#define POPCOUNT(x) popcount(x)
#endif

#define SDARRAY_RANK1 (1<<0)
//#define SDARRAY_RANK0 (1<<1)
//#define SDARRAY_RR (1<<1)
#define SDARRAY_SPARSE (1<<1)
#define SDARRAY_SELECT1 (1<<2)
#define SDARRAY_SELECT0 (1<<3)
//#define SDARRAY_LENHUF (1<<4)
#define SDARRAY_NOBUF (1<<5)
//#define SDARRAY_COMPRANK (1<<6)
//#define SDARRAY_COMPPTR (1<<7)
//#define SDARRAY_SUC (1<<8)


#ifndef _BITVEC_T_
#define _BITVEC_T_
typedef u64 bitvec_t;
#endif

typedef struct densearray4 {
  i64 n; // ベクトルの長さ
  i64 m; // 1の数
  i64 k; // 整数のバイト数
  i64 size; // 索引サイズ (ベクトルは含まない)
  bitvec_t *buf; // ベクトル
  int opt; // サポートする操作
  i64 ml[2],ms[2]; // selectの表の要素数
  i64 rrr; // small blockのサイズ

// for select 1,0
  uchar *lp[2];
  uchar *sl[2];
  word *ss[2];
  i64 *p[2];

//for rank
  uchar *rl;
  word *rs;

// for sparsearray
  int low_width;
  bitvec_t *low;

  int (*access)(struct densearray4 *da, i64 i);
  void (*setbit)(struct densearray4 *da, i64 i, int x);
  u64 (*getbits)(struct densearray4 *da, i64 i, int d);
  void (*setbits)(struct densearray4 *da, i64 i, int d, u64 x);

  i64 (*rank)(struct densearray4 *da, i64 i, int c);
  i64 (*select)(struct densearray4 *da, i64 i, int c);
  i64 (*succ)(struct densearray4 *da, i64 i, int c);
  i64 (*pred)(struct densearray4 *da, i64 i, int c);


} densearray4;

densearray4 *densearray4_new(i64 n);
i64 densearray4_makeindex(densearray4 *da, int opt);
i64 densearray4_write(densearray4 *da, FILE *f);
densearray4 *densearray4_read(uchar **map);
void densearray4_free(densearray4 *da);

#endif // _DENSEARRAY_H_
