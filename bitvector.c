/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/

// supports select1 and select0

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/timeb.h>
#include "densearray4.h"

#define ID_DENSEARRAY 0x10

#define CHECK
//#define RANDOM

#define logD 6

#define PBS (sizeof(bitvec_t)*8)
#define D (1<<logD)
#define logM 5
#define M (1<<logM)
#define logP 8
#define P (1<<logP)
#define logLL 16    // size of word
#define LL (1<<logLL)
//#define logLLL 7
//#define LLL 128
//#define LLL 32
#define logLLL 5
//#define logLLL 2
#define LLL (1<<logLLL)
//#define logL 10
//#define logL (logLL-3)
#define logL (logLL-1-5)
#define L (1<<logL)

#define logRR 16
#define RR (1<<logRR)
//#define logRRR 9
//#define RRR (1<<logRRR)

#ifndef min
 #define min(x,y) ((x)<(y)?(x):(y))
#endif


static int msize=0;
#define mymalloc(p,n,f) {p = malloc((n)*sizeof(*p)); msize += (f)*(n)*sizeof(*p); /* if (f) printf("malloc %d bytes at line %d total %d\n",(n)*sizeof(*p),__LINE__,msize);  */ if ((p)==NULL) {printf("not enough memory (%d bytes) in line %d\n",msize,__LINE__); exit(1);};}

static int blog(i64 x)
{
i64 l;
  l = -1;
  while (x>0) {
    x>>=1;
    l++;
  }
  return l;
}


static void writeuint(int k,u64 x,FILE *f)
{
  int i;
  for (i=k-1; i>=0; i--) {
    fputc(x & 0xff,f);
    x >>= 8;
  }
}

static u64 getuint(uchar *s, i64 i, i64 w)
{
  u64 x;
  i64 j;
  s += i*w;
  x = 0;
  for (j=0; j<w; j++) {
    x += ((u64)(*s++)) << (j*8);
  }
  return x;
}

static void putuint(uchar *s, i64 i, i64 x, i64 w)
{
  i64 j;
  s += i*w;
  for (j=0; j<w; j++) {
    *s++ = x & 0xff;
    x >>= 8;
  }
}


static int setbit(bitvec_t *B, i64 i,int x)
{
  i64 j,l;

  j = i / D;
  l = i & (D-1);
  if (x==0) B[j] &= (~(1L<<(D-1-l)));
  else if (x==1) B[j] |= (1L<<(D-1-l));
  else {
    printf("error setbit x=%d\n",x);
    exit(1);
  }
  return x;
}

static void densearray4_setbit(densearray4 *da, i64 i,int x)
{
  i64 j,l;
  bitvec_t *B;

  B = da->buf;
  j = i / D;
  l = i & (D-1);
  if (x==0) B[j] &= (~(1L<<(D-1-l)));
  else if (x==1) B[j] |= (1L<<(D-1-l));
  else {
    printf("error setbit x=%d\n",x);
    exit(1);
  }
}

static void setbits(bitvec_t *B, i64 i, int d, u64 x)
{
  int j;

  for (j=0; j<d; j++) {
    setbit(B,i+j,(x>>(d-j-1))&1);
  }
}

static void densearray4_setbits(densearray4 *da, i64 i, int d, u64 x)
{
  u64 y,m;
  int d2;
  bitvec_t *B;
  
  B = da->buf;
  B += (i>>logD);
  i &= (D-1);

  while (i+d > D) {
    d2 = D-i; // x の上位 d2 ビットを格納
    y = x >> (d-d2);
    m = (1<<d2)-1;
    *B = (*B & (~m)) | y;
    B++;  i=0;
    d -= d2;
    x &= (1<<d)-1; // x の上位ビットを消去
  }
  m = (1<<d)-1;
  y = x << (D-i-d);
  m <<= (D-i-d);
  *B = (*B & (~m)) | y;

}

static int getbit(bitvec_t *B, i64 i)
{
  i64 j,l;

  j = i >> logD;
  l = i & (D-1);
  return (B[j] >> (D-1-l)) & 1;
}

static int densearray4_getbit(densearray4 *da, i64 i)
{
  i64 j,l;

  j = i >> logD;
  l = i & (D-1);
  return (da->buf[j] >> (D-1-l)) & 1;
}

static bitvec_t densearray4_getbits(densearray4 *da, i64 i, int d)
{
  bitvec_t x,z;
  bitvec_t *B;

  B = da->buf;
  B += (i >>logD);
  i &= (D-1);
  if (i+d <= D) {
    x = B[0];
    x <<= i;
    x >>= (D-d);  // D==64, d==0 だと動かない
  } else {
    x = B[0] << i;
    x >>= D-d;
    z = B[1] >> (D-(i+d-D));
    x += z;
  }
  return x;
}

#if 0
static const unsigned int popCount[] = {
0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8
};
#endif

static unsigned int selecttbl[8*256];

static unsigned int popcount(bitvec_t x)
{
  bitvec_t r;
//  __uint128_t rr;
#if 1
  r = x;
  r = ((r & 0xaaaaaaaaaaaaaaaa)>>1) + (r & 0x5555555555555555);
  r = ((r & 0xcccccccccccccccc)>>2) + (r & 0x3333333333333333);
  r = ((r>>4) + r) & 0x0f0f0f0f0f0f0f0f;
//  r = (r>>8) + r;
//  r = (r>>16) + r;
//  r = ((r>>32) + r) & 127;

  r *= 0x0101010101010101;
  r >>= 64-8;
//  printf("r1 %016lx\n",r);
//  rr = r;
//  rr *= 0x0101010101010101;
//  r = rr >> 56;
//  printf("r2 %016lx\n",r);
//  r &= 0xff;
#else
  r = popCount[x & 0xff];
  x >>= 8;
  r += popCount[x & 0xff];
  x >>= 8;
  r += popCount[x & 0xff];
  x >>= 8;
  r += popCount[x & 0xff];
#if 1
  x >>= 8;
  r += popCount[x & 0xff];
  x >>= 8;
  r += popCount[x & 0xff];
  x >>= 8;
  r += popCount[x & 0xff];
  x >>= 8;
  r += popCount[x & 0xff];
#endif
#endif
  return r;
}

#if 0
static int select_sub(bitvec_t x, int s)
{
  bitvec_t r, w;
  __uint128_t rr;
  int p;

  p = 0;

  r = x;
  r = ((r & 0xaaaaaaaaaaaaaaaa)>>1) + (r & 0x5555555555555555);
  r = ((r & 0xcccccccccccccccc)>>2) + (r & 0x3333333333333333);
  r = ((r>>4) + r) & 0x0f0f0f0f0f0f0f0f;

  rr = r;
  rr *= 0x0101010101010101;
  w = rr >> 56;
  w += 0x8080808080808080;
  w -= s * 0x0101010101010101;
  w &= 0x8080808080808080;
  if ((w & 0xffffffff00000000) == 0) {
    p += 32;
    x <<= 32;
  }
  if ((w & 0xffff000000000000) == 0) {
    p += 16;
    x <<= 16;
  }
  if ((w & 0xff00000000000000) == 0) {
    p += 8;
    x <<= 8;
  }
  return p;
}
#endif

///////////////////////////////////////////////////////////
// rank(i, c)
//   returns the number of c's in [0,i]
///////////////////////////////////////////////////////////
static i64 densearray4_rank(densearray4 *da, i64 i, int c)
{
  i64 r,j;
  bitvec_t *p;
  i64 rrr;

#ifdef DEBUG
   if (c < 0 || c > 1) {
    printf("densearray4_rank: error c = %d\n", c);
    exit(1);
  }
#endif
  if (i < 0) return 0;
  if (i >= da->n) i = da->n-1;
  
  rrr = da->rrr;
  r = getuint(da->rl,i>>logRR,da->k) + da->rs[i/rrr];
  p = da->buf + ((i & (~(rrr-1))) >> logD);
  j = i & (rrr-1);
  while (j >= D) {
    r += POPCOUNT(*p++);
    j -= D;
  }
  r += POPCOUNT(*p >> (D-1-j));

  if (c == 0) {
    r = (i+1) - r;
  }

  return r;
}

static i64 densearray4_rank_and_bit(densearray4 *da, i64 i, int *c)
{
  i64 r,j;
  bitvec_t *p,x;
  i64 rrr;

  if (i < 0) return 0;
  if (i >= da->n) i = da->n-1;

  rrr = da->rrr;
  r = getuint(da->rl,i>>logRR,da->k) + da->rs[i/rrr];
  p = da->buf + ((i & (~(rrr-1))) >> logD);
  j = i & (rrr-1);
  while (j >= D) {
    r += POPCOUNT(*p++);
    j -= D;
  }
  x = *p >> (D-1-j);
  r += POPCOUNT(x);
  *c = (int)(x & 1);

  return r;
}

static i64 densearray4_rank0(densearray4 *da, i64 i)
{
  if (i < 0) return 0;
  if (i >= da->n) i = da->n-1;
  return i+1 - densearray4_rank(da,i, 1);
}

static i64 densearray4_select0_naive(densearray4 *da, i64 i)
{
  i64 l,r,m;

  l = 0; r = da->n-1;
  while (l <= r) {
    m = (l+r)/2;
    if (i <= densearray4_rank0(da,m)) r = m-1;  else l = m+1;
  }
  return l;
}

static i64 densearray4_select_by_rank(densearray4 *da, i64 ith, int c)
{
  bitvec_t x, *buf;
  
  i64 j,k;
  i64 ofs, ofs2;
  i64 r, r2;
  int rr;
  i64 d,runlen;

  i64 ll, rl, ml, pl; // 大ブロックの2分探索用
  i64 p; // 答え
  i64 ii;
  i64 rrr;
  bitvec_t *q;
  i64 p0;
  i64 l0;

  rrr = da->rrr;

  ii = ith;

  ll = 0;  rl = (da->n-1) >> logRR;
  pl = ll;  r2 = 0;
  while (ll <= rl) {
    ml = (ll+rl) >> 1;
    r = getuint(da->rl,ml,da->k);
    if (c == 0) {r = (ml<<logRR) - r;} // select0 の場合
    if (r < ith) {
      pl = ml;
      ll = ml+1;
      r2 = r;
    } else {
      rl = ml-1;
    }
  }
  ith -= r2; // 大ブロック内のランク
  ll = (pl<<logRR)/rrr;
  rl = min(((pl+1)<<logRR)/rrr-1, (da->n-1)/rrr);

//  r = da->rs[rl];
//  if (c == 0) {r = ((rl-ll)*rrr) - r;} // select0 の場合
//  if (r < ith) {
//    printf("ith = %ld r = %ld\n", ith, r);
//  }

  pl = ll;  r2 = 0;
  l0 = ll;
  while (ll <= rl) {
    ml = (ll+rl) >> 1;
    r = da->rs[ml];
    if (c == 0) {r = ((ml-l0)*rrr) - r;} // select0 の場合
    if (r < ith) {
      pl = ml;
      ll = ml+1;
      r2 = r;
    } else {
      rl = ml-1;
    }
  }
  ith -= r2; // 小ブロック内のランク

  p = pl * rrr; // ith を含む小ブロックの最初のビット位置
  rr = 0;

  p0 = p;

  q = &(da->buf[p>>logD]);

  while (1) {
    x = *q;
    if (c == 0) x = ~x;
    rr = POPCOUNT(x);
    if (rr >= ith) break;
    ith -= rr;
    p += D;
    q++;
  }
      
  x = *q;
  if (c == 0) x = ~x;
  while (1) {
    rr = POPCOUNT(x >> (D-8));
    if (rr >= ith) break;
    ith -= rr;
    p += 8;
    x <<= 8;
  }
  p += selecttbl[((ith-1)<<8)+(x>>(D-8))];

  if (p >= p0 + rrr) {
    printf("i = %ld p = %ld p0 = %ld\n", ii, p, p0);
  }

#if 0
  if (p != densearray4_select0_naive(da, ii)) {
    printf("select0: i = %ld p = %ld naive = %ld\n", 
      ii, p, densearray4_select0_naive(da, ii));
    exit(1);
  }
#endif
  return p;
  
}

///////////////////////////////////////////////////////////
// select(i, f)
//   returns the position of i-th f. The position is in [0,n-1].
//   if i == 0, returns -1
//   if f is larger than the number of f's in the vector, returns n
///////////////////////////////////////////////////////////
static i64 densearray4_select(densearray4 *da, i64 i,int f)
{
  
  //  dword *s;
  i64 p,r;
  i64 il;
  i64 rr;
  bitvec_t x;
  bitvec_t *q;
  i64 m;
  i64 stmp;

  if (i == 0) return -1;
  m = (f == 1) ? da->m : da->n - da->m;
  if (i > m) return da->n;

  if (f == 0 && (da->opt & SDARRAY_SELECT0) == 0) {
    return densearray4_select_by_rank(da,i,0);
  }
  if (f == 1 && (da->opt & SDARRAY_SELECT1) == 0) {
    return densearray4_select_by_rank(da,i,1);
  }

  i--;

  il = da->p[f][i>>logL];
  if (il < 0) {
    il = -il-1;
    p = getuint(da->sl[f],(il<<logL)+(i & (L-1)),da->k);
  } else {
    p = getuint(da->lp[f],i>>logL,da->k);
    p += da->ss[f][(il<<(logL-logLLL))+(i & (L-1))/LLL];
    r = i - (i & (LLL-1));

    q = &(da->buf[p>>logD]);

    if (f == 1) {
      rr = p & (D-1);
      r -= POPCOUNT(*q >> (D-1-rr));
      p = p - rr;
      
      while (1) {
        rr = POPCOUNT(*q);
        if (r + rr >= i) break;
        r += rr;
        p += D;
        q++;
      }
      
      x = *q;
      while (1) {
        rr = POPCOUNT(x >> (D-8));
        if (r + rr >= i) break;
        r += rr;
        p += 8;
        x <<= 8;
      }
      p += selecttbl[((i-r-1)<<8)+(x>>(D-8))];
    } else {
      rr = p & (D-1);
      r -= POPCOUNT((~(*q))  >> (D-1-rr));
      p = p - rr;
      
      while (1) {
        rr = POPCOUNT(~(*q));
        if (r + rr >= i) break;
        r += rr;
        p += D;
        q++;
      }
      
      x = ~(*q);

      while (1) {
        rr = POPCOUNT(x >> (D-8));
        if (r + rr >= i) break;
        r += rr;
        p += 8;
        x <<= 8;
      }
      p += selecttbl[((i-r-1)<<8)+(x>>(D-8))];
    }
  }
  
  return p;
}

///////////////////////////////////////////////////////////
// pred(x, c)
//   returns max{i | B[i]=c and i >= x}
//   if there is no such i, returns -1
//   (if x >= n, error, or the position of last c?)
///////////////////////////////////////////////////////////
static i64 densearray4_pred(densearray4 *da, i64 x, int c)
{
  i64 s0;
  i64 r;
  i64 m;
  i64 s;

  r = densearray4_rank(da, x, c);
  s = densearray4_select(da, r, c);
  
  return s;
}


// rank(min{i | B[i]=c and i >= x})
static i64 densearray4_succ_rank(densearray4 *da, i64 x)
{
  i64 r;

  r = (x > 0) ? densearray4_rank(da, x-1, 1) + 1 : 1;

  return r;
}

// min{i | B[i]=c and i >= x}
static i64 densearray4_succ(densearray4 *da, i64 x, int c)
{
  i64 s0;
  i64 r;
  i64 m;
  i64 s;

  bitvec_t *q, z;
  int i,d, rr;

  if (x < 0) x = 0; // debug

  s = x;

  d = x & (D-1);
  q = &(da->buf[x>>logD]);

  z = *q++;
  if (c == 0) z = ~z;
  z <<= d;
  z >>= d; // 上位 d ビットを 0 にする
  s &= ~(D-1);

  for (i=0; i < 4; i++) {
    rr = POPCOUNT(z);
    if (rr > 0) break;
    s += D;
    z = *q++;
    if (c == 0) z = ~z;
  }

  if (i < 4) {
    while (1) {
      rr = POPCOUNT(z >> (D-8));
      if (rr > 0) break;
      s += 8;
      z <<= 8;
    }
    s += selecttbl[((1-1)<<8)+(z>>(D-8))];

    if (s > da->n) s = da->n;
/*
    s0 = s;
    r = (x > 0) ? densearray4_rank(da, x-1, c) + 1 : 1;
    m = (c == 1) ? da->m : da->n - da->m;
    s = (r <= m) ? densearray4_select(da, r, c) : da->n;

    if (s0 != s) {
      printf("succ: x = %ld s = %ld s0 = %ld\n", x, s, s0);
    }
*/
  } else {
    r = (x > 0) ? densearray4_rank(da, x-1, c) + 1 : 1;
    m = (c == 1) ? da->m : da->n - da->m;
    s = (r <= m) ? densearray4_select(da, r, c) : da->n;
  }
  return s;
}

static void densearray4_make_selecttbl(void)
{
  i64 i,x,r;
  bitvec_t buf[1];

  for (x = 0; x < 256; x++) {
    setbits(buf,0,8,x);
    for (r=0; r<8; r++) selecttbl[(r<<8)+x] = -1;
    r = 0;
    for (i=0; i<8; i++) {
      if (getbit(buf,i)) {
        selecttbl[(r<<8)+x] = i;
        r++;
      }
    }
  }

}

static int densearray4_construct_init(densearray4 *da, i64 n)
{
  i64 i;

  da->n = n;
  mymalloc(da->buf, (n+D-1)/D, 0);
  for (i=0; i<(n+D-1)/D; i++) da->buf[i] = 0;

  return 0;
}

static void densearray4_notsupported(void)
{
  printf("densearray4: not supported.\n");
  exit(1);
}

densearray4 *densearray4_new(i64 n)
{
  densearray4 *da;
  i64 i;

  mymalloc(da, 1, 0);
  da->n = n;
  mymalloc(da->buf, (n+D-1)/D, 0);
  for (i=0; i<(n+D-1)/D; i++) da->buf[i] = 0;

  da->rl = NULL;
  da->rs = NULL;
  da->lp[0] = da->lp[1] = NULL;
  da->sl[0] = da->sl[1] = NULL;
  da->ss[0] = da->ss[1] = NULL;
  da->p[0] = da->p[1] = NULL;

  da->access = densearray4_getbit;
  da->setbit = densearray4_setbit;
  da->getbits = densearray4_getbits;
  da->setbits = densearray4_setbits;

  da->rank = (i64 (*)(struct densearray4 *, i64, int))densearray4_notsupported;
  da->select = (i64 (*)(struct densearray4 *, i64, int))densearray4_notsupported;
  da->succ = (i64 (*)(struct densearray4 *, i64, int))densearray4_notsupported;
  da->pred = (i64 (*)(struct densearray4 *, i64, int))densearray4_notsupported;

  return da;
}

void densearray4_free(densearray4 *da)
{
  int i;
  free(da->buf);
  if (da->rl) free(da->rl);
  if (da->rs) free(da->rs);
  for (i=0; i<2; i++) {
    if (da->lp[i]) free(da->lp[i]);
    if (da->sl[i]) free(da->sl[i]);
    if (da->ss[i]) free(da->ss[i]);
    if (da->p[i]) free(da->p[i]);
  }
  free(da);
}

#if 0
int densearray4_construct_set(densearray4 *da, i64 i, int x)
{
  if (x > 0) setbit(da->buf,i,x);
  return 0;
}

i64 densearray4_construct_end(densearray4 *da, ushort l, int opt)
{
  densearray4_construct(da, da->n, da->buf, opt);
//  return da->size;
  return da->size + sizeof(*da->buf) * ((da->n+D-1)/D);
}
#endif

static void densearray4_construct(densearray4 *da, i64 n, bitvec_t *buf, int opt)
{
  i64 i,j,m,k;
  i64 nl;
  i64 p,pp;
  i64 il,is,ml[2],ms[2];
  i64 r;
  i64 size;
//  dword *s;
  i64 p1,m1,p2,m2;
  i64 rrr;

i64 s1,s2,s3,s4;
  
#if 0
  i64 b;
  i64 freq[RRR+1];
#endif

  size = sizeof(densearray4);
//  printf("densearray-size:0 %ld\n",size);

  densearray4_make_selecttbl(); // 最初に1回だけでいいので消す?

  da->rrr = rrr = opt >> 16;
  opt &= 0xffff;

  if (L/LLL == 0) {
    printf("ERROR: L=%d LLL=%d\n",L,LLL);
    exit(1);
  }

  m = 0;
  for (i=0; i<n; i++) m += getbit(buf,i);
  da->n = n;
  da->m = m;

//  da->k = k = (blog(m+1)+1+8-1)/8;
  da->k = k = (blog(n+1)+1+8-1)/8;


//  printf("n=%d m=%d (%f)\n",n,m,(double)m/n);

  da->buf = buf;

//  mymalloc(s,m,0);
  //da->s = s;
  m = 0;
  for (i=0; i<n; i++) {
    if (getbit(buf,i)) {
//      s[m] = i;
      m++;
    }
  }

  if (opt & SDARRAY_SELECT1) {
    nl = (m-1) / L + 1;
    mymalloc(da->lp[1],(nl+1)*k,1);  size += sizeof(*(da->lp[1]))*(nl+1)*k;
//    for (il = 0; il < nl+1; il++) da->lp[1][il] = 0;
    for (il = 0; il < nl+1; il++) putuint(da->lp[1],il,0,k);

    mymalloc(da->p[1],nl+1,1);  size += sizeof(*(da->p[1]))*(nl+1);
    for (il = 0; il < nl+1; il++) da->p[1][il] = 0;
//    printf("densearray-size:1 %ld\n",size);

    for (r = 0; r < 2; r++) {
      ml[1] = ms[1] = 0;
      p1 = p2 = -1;
      m1 = m2 = -1;
      for (il = 0; il < nl; il++) {
  //      pp = s[il*L];
        while (m1 < il*L) m1 += getbit(buf,++p1);
        if (p1 >= n) printf("???1 p1 = %ld\n",p1);
  //      if (pp != p1) printf("pp = %ld  p1 = %ld\n",pp,p1);
        pp = p1;
//        da->lp[1][il] = pp;
        putuint(da->lp[1],il,pp,k);
        i = min((il+1)*L-1,m-1);
  //      p = s[i];
        while (m2 < i) m2 += getbit(buf,++p2);
        if (p2 >= n) printf("???2 p2 = %ld\n",p2);
  //      if (p != p2) printf("p = %ld  p2 = %ld\n",p,p2);
        p = p2;
        //printf("%d ",p-pp);
        if (p - pp >= LL) {
          if (r == 1) {
            for (is = 0; is < L; is++) {
              if (il*L+is >= m) break;
              while (m1 < il*L+is) m1 += getbit(buf,++p1);
              if (p1 >= n) printf("???3 p1 = %ld\n",p1);
  //          if (s[il*L+is] != p1) printf("1:s = %ld  p1 = %ld\n",s[il*L+is],p1);
//            da->sl[ml*L+is] = s[il*L+is];
//              da->sl[1][ml[1]*L+is] = p1;
              putuint(da->sl[1],ml[1]*L+is,p1,k);
            }
          }
          da->p[1][il] = -(ml[1]+1);
          ml[1]++;
        } else {
          if (r == 1) {
            for (is = 0; is < L/LLL; is++) {
              if (il*L+is*LLL >= m) break;
              while (m1 < il*L+is*LLL) m1 += getbit(buf,++p1);
              if (p1 >= n) printf("???4 p1 = %ld\n",p1);
//            if (s[il*L+is*LLL] != p1) 
//              printf("2:s = %ld  p1 = %ld\n",s[il*L+is*LLL],p1);
//            da->ss[ms*(L/LLL)+is] = s[il*L+is*LLL] - pp;
              da->ss[1][ms[1]*(L/LLL)+is] = p1 - pp;
            }
          }
          da->p[1][il] = ms[1];
          ms[1]++;
        }
      }
      if (r == 0) {
        da->ml[1] = ml[1];  da->ms[1] = ms[1];
//        mymalloc(da->sl[1],ml[1]*L+1,1);  size += sizeof(*(da->sl[1]))*(ml[1]*L+1);
//        for (il = 0; il < ml[1]*L+1; il++) da->sl[1][il] = 0;
        mymalloc(da->sl[1],(ml[1]*L+1)*k,1);  size += sizeof(*(da->sl[1]))*(ml[1]*L+1)*k;
        for (il = 0; il < ml[1]*L+1; il++) putuint(da->sl[1],il,0,k);
        mymalloc(da->ss[1],ms[1]*(L/LLL)+1,1);
        for (il = 0; il < ms[1]*(L/LLL)+1; il++) da->ss[1][il] = 0;
        size += sizeof(*(da->ss[1]))*(ms[1]*(L/LLL)+1);
//        printf("densearray-size:2 %ld\n",size);
      }
    }
  } else {
    da->lp[1] = NULL;  da->p[1] = NULL;
    da->sl[1] = NULL;  da->ss[1] = NULL;
  }

  if (opt & SDARRAY_SELECT0) {
    i64 m0;
    m0 = n - m;
    nl = (m0-1) / L + 1;
    mymalloc(da->lp[0],(nl+1)*k,1);  size += sizeof(*(da->lp[0]))*(nl+1)*k;
//    for (il = 0; il < nl+1; il++) da->lp[0][il] = 0;
    for (il = 0; il < nl+1; il++) putuint(da->lp[0],il,0,k);
    mymalloc(da->p[0],nl+1,1);  size += sizeof(*(da->p[0]))*(nl+1);
    for (il = 0; il < nl+1; il++) da->p[0][il] = 0;
//    printf("densearray-size:1 %ld\n",size);

    for (r = 0; r < 2; r++) {
      ml[0] = ms[0] = 0;
      p1 = p2 = -1;
      m1 = m2 = -1;
      for (il = 0; il < nl; il++) {
  //      pp = s[il*L];
        while (m1 < il*L) m1 += 1-getbit(buf,++p1);
        if (p1 >= n) printf("???1 p1 = %ld\n",p1);
  //      if (pp != p1) printf("pp = %ld  p1 = %ld\n",pp,p1);
        pp = p1;
//        da->lp[0][il] = pp;
        putuint(da->lp[0],il,pp,k);
        i = min((il+1)*L-1,m0-1);
  //      p = s[i];
        while (m2 < i) m2 += 1-getbit(buf,++p2);
        if (p2 >= n) printf("???2 p2 = %ld\n",p2);
  //      if (p != p2) printf("p = %ld  p2 = %ld\n",p,p2);
        p = p2;
        //printf("%d ",p-pp);
        if (p - pp >= LL) {
          if (r == 1) {
            for (is = 0; is < L; is++) {
              if (il*L+is >= m0) break;
              while (m1 < il*L+is) m1 += 1-getbit(buf,++p1);
              if (p1 >= n) printf("???3 p1 = %ld\n",p1);
  //          if (s[il*L+is] != p1) printf("1:s = %ld  p1 = %ld\n",s[il*L+is],p1);
//            da->sl[ml*L+is] = s[il*L+is];
//              da->sl[0][ml[0]*L+is] = p1;
              putuint(da->sl[0],ml[0]*L+is,p1,k);
            }
          }
          da->p[0][il] = -(ml[0]+1);
          ml[0]++;
        } else {
          if (r == 1) {
            for (is = 0; is < L/LLL; is++) {
              if (il*L+is*LLL >= m0) break;
              while (m1 < il*L+is*LLL) m1 += 1-getbit(buf,++p1);
              if (p1 >= n) printf("???4 p1 = %ld\n",p1);
//            if (s[il*L+is*LLL] != p1) 
//              printf("2:s = %ld  p1 = %ld\n",s[il*L+is*LLL],p1);
//            da->ss[ms*(L/LLL)+is] = s[il*L+is*LLL] - pp;
              da->ss[0][ms[0]*(L/LLL)+is] = p1 - pp;
            }
          }
          da->p[0][il] = ms[0];
          ms[0]++;
        }
      }
      if (r == 0) {
        da->ml[0] = ml[0];  da->ms[0] = ms[0];
//        mymalloc(da->sl[0],ml[0]*L+1,1);  size += sizeof(*(da->sl[0]))*(ml[0]*L+1);
//        for (il = 0; il < ml[0]*L+1; il++) da->sl[0][il] = 0;
        mymalloc(da->sl[0],(ml[0]*L+1)*k,1);  size += sizeof(*(da->sl[0]))*(ml[0]*L+1)*k;
        for (il = 0; il < ml[0]*L+1; il++) putuint(da->sl[0],il,0,k);
        mymalloc(da->ss[0],ms[0]*(L/LLL)+1,1);
        for (il = 0; il < ms[0]*(L/LLL)+1; il++) da->ss[0][il] = 0;
        size += sizeof(*(da->ss[0]))*(ms[0]*(L/LLL)+1);
//        printf("densearray-size:2 %ld\n",size);
      }
    }
  } else {
    da->lp[0] = NULL;  da->p[0] = NULL;
    da->sl[0] = NULL;  da->ss[0] = NULL;
  }


// rank index
  if (opt & SDARRAY_RANK1) {
//    mymalloc(da->rl,(n+RR-1)/RR,1);  
//      size += sizeof(*(da->rl))*((n+RR-1)/RR);
    mymalloc(da->rl,(n+RR-1)/RR*k,1);  
      size += sizeof(*(da->rl))*((n+RR-1)/RR)*k;
    mymalloc(da->rs,(n+rrr-1)/rrr,1);
      size += sizeof(*(da->rs))*((n+rrr-1)/rrr);
//  printf("densearray-size:3 %ld\n",size);
    r = 0;
    for (i=0; i<n; i+=RR) {
//      da->rl[i/RR] = r;
      putuint(da->rl,i/RR,r,k);
      m = 0;
      for (j=0; j<RR; j++) {
        if (j % rrr == 0 && i+j < n) {
          da->rs[(i+j)/rrr] = m;
        }
        if (i+j < n && getbit(buf,i+j)==1) m++;
      }
      r += m;
    }
  } else {
    da->rl = NULL;  da->rs = NULL;
  }
//
  da->opt = opt;
  da->size = size;


//  da->access = densearray4_getbit;
  da->rank = densearray4_rank;
  da->select = densearray4_select;
  da->succ = densearray4_succ;
  da->pred = densearray4_pred;

  return;
}

i64 densearray4_makeindex(densearray4 *da, int opt)
{
  densearray4_construct(da, da->n, da->buf, opt);
  return da->size + sizeof(*da->buf) * ((da->n+D-1)/D);
}

i64 densearray4_write(densearray4 *da, FILE *f)
{
  i64 i,k;
  i64 nl;
  i64 rrr;

  i64 size;

  k = da->k;
  rrr = da->rrr;
  size = 0;
  writeuint(1,ID_DENSEARRAY,f);
  writeuint(1, da->k, f);
  writeuint(sizeof(da->n), da->n, f);
  writeuint(sizeof(da->m), da->m, f);
  writeuint(4, da->opt, f);
  size += 1 + sizeof(da->n) + sizeof(da->m) + 1;
  if (da->opt & SDARRAY_RANK1) {
    writeuint(sizeof(int), RR, f);
    writeuint(sizeof(int), rrr, f);
    size += 2*sizeof(int);
  }
  if (da->opt & SDARRAY_SELECT1) {
    writeuint(sizeof(da->ml[1]), da->ml[1], f);  size += sizeof(da->ml[1]);
    writeuint(sizeof(da->ms[1]), da->ms[1], f);  size += sizeof(da->ms[1]);
  }
  if (da->opt & SDARRAY_SELECT0) {
    writeuint(sizeof(da->ml[0]), da->ml[0], f);  size += sizeof(da->ml[0]);
    writeuint(sizeof(da->ms[0]), da->ms[0], f);  size += sizeof(da->ms[0]);
  }

//  writeuint(sizeof(i64), size, f);

  if (da->opt & SDARRAY_RANK1) {
    for (i=0; i<(da->n+RR-1)/RR; i++) {
      writeuint(k, getuint(da->rl,i,k), f); size += k;
    }
    for (i=0; i<(da->n+rrr-1)/rrr; i++) {
      writeuint(sizeof(*da->rs), da->rs[i], f);
      size += sizeof(*da->rs);
    }
  }

  if (da->opt & SDARRAY_SELECT1) {
    nl = (da->m-1) / L + 1;
    for (i=0; i<nl+1; i++) {
//      writeuint(sizeof(*da->lp[1]), da->lp[1][i], f);  size += sizeof(*da->lp[1]);
      writeuint(k, getuint(da->lp[1],i,k), f);  size += k;
    }
    for (i=0; i<nl+1; i++) {
      writeuint(sizeof(*da->p[1]), da->p[1][i], f);  size += sizeof(*da->p[1]);
    }
    for (i=0; i<da->ml[1]*L+1; i++) {
//      writeuint(sizeof(*da->sl[1]), da->sl[1][i], f);  size += sizeof(*da->sl[1]);
      writeuint(k, getuint(da->sl[1],i,k), f);  size += k;
    }
    for (i=0; i<da->ms[1]*(L/LLL)+1; i++) {
      writeuint(sizeof(*da->ss[1]), da->ss[1][i], f);  size += sizeof(*da->ss[1]);
    }
  }
  if (da->opt & SDARRAY_SELECT0) {
    nl = (da->n-da->m-1) / L + 1;
    for (i=0; i<nl+1; i++) {
//      writeuint(sizeof(*da->lp[0]), da->lp[0][i], f);  size += sizeof(*da->lp[0]);
      writeuint(k, getuint(da->lp[0],i,k), f);  size += k;
    }
    for (i=0; i<nl+1; i++) {
      writeuint(sizeof(*da->p[0]), da->p[0][i], f);  size += sizeof(*da->p[0]);
    }
    for (i=0; i<da->ml[0]*L+1; i++) {
//      writeuint(sizeof(*da->sl[0]), da->sl[0][i], f);  size += sizeof(*da->sl[0]);
      writeuint(k, getuint(da->sl[0],i,k), f);  size += k;
    }
    for (i=0; i<da->ms[0]*(L/LLL)+1; i++) {
      writeuint(sizeof(*da->ss[0]), da->ss[0][i], f);  size += sizeof(*da->ss[0]);
    }
  }

  if (!(da->opt & SDARRAY_NOBUF)) {
    for (i=0; i<(da->n+D-1)/D; i++) {
      writeuint(sizeof(*da->buf), da->buf[i], f);
      size += sizeof(*da->buf);
    }
  }
  return size;
}

densearray4 *densearray4_read(uchar **map)
{
  densearray4 *da;
  i64 nl;
  uchar *p;

  mymalloc(da, 1, 0);

  
  i64 rr, rrr;

//  make_selecttbl();

  p = *map;

  rr = getuint(p,0,1);  p += 1;
  if (rr != ID_DENSEARRAY) {
    printf("densearray4_read: id = %ld is not supported.\n",rr);
    exit(1);
  }

  da->k = getuint(p,0,1);  p += 1;
  da->n = getuint(p,0,sizeof(da->n));  p += sizeof(da->n);
  da->m = getuint(p,0,sizeof(da->m));  p += sizeof(da->m);
  da->opt = getuint(p,0,4);  p += 4;
//  printf("densearray_read: n=%ld m=%ld opt=%d\n",da->n,da->m,da->opt);
  if (da->opt & SDARRAY_RANK1) {
    rr = getuint(p,0,sizeof(int));  p += sizeof(int);
    if (rr != RR) {
      printf("error2 RR=%ld must be %d\n",rr,RR);
    }
    rrr = getuint(p,0,sizeof(int));  p += sizeof(int);
    da->rrr = rrr;
//    printf("RRR = %ld\n",rrr);
#if 0
    if (rrr != RRR) {
      printf("error RRR=%ld must be %d\n",rrr,RRR);
    }
#endif
  }
  if (da->opt & SDARRAY_SELECT1) {
    da->ml[1] = getuint(p,0,sizeof(da->ml[1]));  p += sizeof(da->ml[1]);
    da->ms[1] = getuint(p,0,sizeof(da->ms[1]));  p += sizeof(da->ms[1]);
  }
  if (da->opt & SDARRAY_SELECT0) {
    da->ml[0] = getuint(p,0,sizeof(da->ml[0]));  p += sizeof(da->ml[0]);
    da->ms[0] = getuint(p,0,sizeof(da->ms[0]));  p += sizeof(da->ms[0]);
  }
//  size = getuint(p,0,sizeof(size));  p += sizeof(size);
//  printf("size %ld\n",size);

  if (da->opt & SDARRAY_RANK1) {
    da->rl = p;
    p += sizeof(*da->rl) * ((da->n+rr-1)/rr) * da->k;
    da->rs = (word *)p;
    p += sizeof(*da->rs) * ((da->n+rrr-1)/rrr);
  }

  if (da->opt & SDARRAY_SELECT1) {
    nl = (da->m-1) / L + 1;
//    da->lp[1] = (dword *)p;
//    p += sizeof(*da->lp[1]) * (nl+1);
    da->lp[1] = p;
    p += sizeof(*da->lp[1]) * (nl+1) * da->k;
    da->p[1] = (i64 *)p;
    p += sizeof(*da->p[1]) * (nl+1);
//    da->sl[1] = (dword *)p;
//    p += sizeof(*da->sl[1]) * (da->ml[1]*L+1);
    da->sl[1] = p;
    p += sizeof(*da->sl[1]) * (da->ml[1]*L+1) * da->k;
    da->ss[1] = (word *)p;
    p += sizeof(*da->ss[1]) * (da->ms[1]*(L/LLL)+1);
  }
  if (da->opt & SDARRAY_SELECT0) {
    nl = (da->n - da->m-1) / L + 1;
//    da->lp[0] = (dword *)p;
//    p += sizeof(*da->lp[0]) * (nl+1);
    da->lp[0] = p;
    p += sizeof(*da->lp[0]) * (nl+1) * da->k;
    da->p[0] = (i64 *)p;
    p += sizeof(*da->p[0]) * (nl+1);
//    da->sl[0] = (dword *)p;
//    p += sizeof(*da->sl[0]) * (da->ml[0]*L+1);
    da->sl[0] = p;
    p += sizeof(*da->sl[0]) * (da->ml[0]*L+1) * da->k;
    da->ss[0] = (word *)p;
    p += sizeof(*da->ss[0]) * (da->ms[0]*(L/LLL)+1);
  }
  
  if (!(da->opt & SDARRAY_NOBUF)) {
    da->buf = (bitvec_t *)p;
    p += sizeof(*da->buf) * ((da->n+D-1)/D);
  }
  *map = p;

  da->access = densearray4_getbit;
  da->setbit = densearray4_setbit;
  da->getbits = densearray4_getbits;
  da->setbits = densearray4_setbits;

  da->rank = densearray4_rank;
  da->select = densearray4_select;
  da->succ = densearray4_succ;
  da->pred = densearray4_pred;

  densearray4_make_selecttbl();

  return da;
}



#ifdef DENSEMAIN
typedef struct timeb mytimestruct;

void mygettime(mytimestruct *t)
{
  ftime(t);
}

double mylaptime(mytimestruct *before,mytimestruct *after)
{
  double t;
  t = after->time - before->time;
  t += (double)(after->millitm - before->millitm)/1000;
  return t;
}


#define N 10485760

int main(int argc, char *argv[])
{
  densearray4 s1,s2;
//  densearray2 s2;

  i64 i,r,n,m,rr;
  u64 hoge,sum;
  FILE *infp = NULL;
  FILE *out;

  bitvec_t *B;
  byte *B2;
  dword *S,*R;
  MMAP *map;
  uchar *mapp;

  double t;
  mytimestruct before,after;

#ifdef __SSE4_2__
  printf("SSE4.2 is available.\n");
#else
  printf("no SSE4.2\n");
#endif

  srand(2);

  n = N; // length of bit vector
  r = 2; // ratio of ones
  m = 0; // number of ones
  if (argc >= 2){
    r = atoi(argv[1]);
    if (r == 0) {
      infp = fopen(argv[2],"rb");
      if (infp == NULL){
        printf("cannot open %s\n",argv[2]);
        return -1;
      }
      fseek(infp,0,SEEK_END);
      n = ftell(infp);
      rewind(infp);
      //printf("n: %d\n",n);
    } else if (argc >= 3) {
      n = atoi(argv[2]);
    }
  }
  rr = r;

  mymalloc(B,(n+PBS-1+P)/PBS,0);

  if (!infp) {
    m = 0;
    for (i = 0; i < n; i++) {
      if (rand() % 100 < r) {
        setbit(B,i,1);
        m++;
      } else {
        setbit(B,i,0);
      }
    }
  } else {
    m = 0;
    for (i = 0; i < n; i++){
      int c = fgetc(infp);
      if (c == EOF){
        printf("unexpected error at reading from %s\n",argv[2]);
        return -1;
      }
      if (c == '1' || c == '(') {
        setbit(B,i,1);
        m++;
      } else if (c == '0' || c == ')') {
        setbit(B,i,0);
      } else {
        printf("unexpected error (2) at reading from %s\n",argv[2]);
        return -1;
      }
    }
  }


  densearray4_construct(&s2,n,B, (128<<16)+(SDARRAY_RANK1 | SDARRAY_SELECT1));
  printf("da: used memory: %d bytes (%lf bpc)\n",s2.size,(double)s2.size*8/n);

#if 1
  out = fopen("sparsearraytmp.dat","w");
  densearray4_write(&s2, out);
  fclose(out);

  map = mymmap("sparsearraytmp.dat");
  if (map->addr==NULL) {
    perror("mmap2\n");
    exit(1);
  }
  mapp = (uchar *)map->addr;

  densearray4_read(&s1, &mapp);
#endif



#ifdef CHECK
  mymalloc(S,n+1,0);
  mymalloc(R,n+1,0);
  r = 0;
  S[r] = -1;
  for (i=0; i<n; i++) {
    if (getbit(B,i)) {
      r++;
      S[r] = i;
    }
    R[i] = r;
  }
#endif

#if 1
  srand(3);

  mygettime(&before);
  hoge = rand();
  sum = 0;
  for (i = 0; i < 100000000; i++) {
    i64 j;
    //j = (rand() % r)+1;
#ifdef RANDOM
    j = hoge % m + 1;
#else
    j = i % m + 1;
#endif
#ifdef CHECK
    if (densearray4_select(&s1,j,1) != S[j]) {
      printf("ERROR: S[%d] = %d, s = %d\n",j,S[j],densearray4_select(&s1,j,1));
    }
    sum += densearray4_select(&s1,j,1);
#else
    sum += densearray4_select(&s1,j,1);
#endif
    hoge = hoge * 1566083941U + 1;
  }
  mygettime(&after);
  t = mylaptime(&before,&after);
  printf("%ld %f densearray4 select time sum=%ld\n",rr,t,sum);
#endif


  return 0;
}

#endif //  DENSEMAIN
