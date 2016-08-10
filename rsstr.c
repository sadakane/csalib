/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "typedef.h"
#include "mman.h"
#include "rsstr.h"



#ifndef min
 #define min(x,y) ((x)<(y)?(x):(y))
#endif

#define mymalloc(p,n) {p = malloc((n)*sizeof(*p)); if ((p)==NULL) {printf("not enough memory (%d bytes) in line %d\n",(n)*sizeof(*p),__LINE__); exit(1);};}

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

static void writeint(int k,i64 x,FILE *f)
{
  int i;
  for (i=k-1; i>=0; i--) {
    fputc(x & 0xff,f);
    x >>= 8;
  }
}

static void writeuint(int k,u64 x,FILE *f)
{
  int i;
  for (i=k-1; i>=0; i--) {
    fputc(x & 0xff,f);
    x >>= 8;
  }
}



static int setbit(bitvec_t *B, i64 i,bitvec_t x)
{
  i64 j,l;
  j = i / DD;
  l = i % DD;
  if (x==0) B[j] &= (~(1L<<(DD-1-l)));
  else if (x==1) B[j] |= (1L<<(DD-1-l));
  else {
    printf("error setbit x=%ld\n",x);
    exit(1);
  }
  return x;
}

static int setbits(bitvec_t *B, i64 i, bitvec_t x, int w)
{
  i64 j;
  
  for (j=0; j<w; j++) {
    setbit(B,i+j,(x>>(w-1-j)) & 1);
  }
  return 0;
}


static int convertchar(uchar t)
{
  int c = 0;
  switch (t) {
  case 'a':  c = 0;  break;
  case 'c':  c = 1;  break;
  case 'g':  c = 2;  break;
  case 't':  c = 3;  break;
  case 'A':  c = 4;  break;
  case 'C':  c = 5;  break;
  case 'G':  c = 6;  break;
  case 'T':  c = 7;  break;
  case '$':  c = 8;  break;
  default:  printf("error char = %c [%02x]\n",t,t);
  }
  return c;
}

static int CC2(int c)
{
  int CC2[9] = {'a', 'c', 'g', 't', 'A', 'C', 'G', 'T', '$'};
  return CC2[c];
}

static int access(rsstr *lf,i64 i)
{
  bitvec_t x;
  x = lf->BW[i/(DD/logSIGMADNA)];
  x = (x >> (DD-logSIGMADNA-(i % (DD/logSIGMADNA))*logSIGMADNA)) & ((1<<logSIGMADNA)-1);
  return CC2(x);
}

#if 1
static i64 rankc_slow(rsstr *lf, i64 i, int c)
{
  int j;
  i64 r;
  i64 i2;
  int lmb,mb;
  int c2;

  c2 = convertchar(c);

//  printf("rankc(%ld)\n",i);

  lmb = lf->lmb;
  mb = 1 << lmb;
  r = getuint(lf->RL,(i/LB)*(SIGMADNA-1) + c2,lf->k);
  r += lf->RM[(i>>lmb)*(SIGMADNA-1) + c2];

  i2 = (i/mb)*mb;
  for (j=0; j <= (i % mb); j++) {
    if (access(lf,i2+j) == c) r++;
  }
  return r;
}
#endif

static bitvec_t popcount4(bitvec_t x)
{
  bitvec_t r;
  r = x;
  r = ((r>>4) + r) & 0x0f0f0f0f0f0f0f0fL;
  r = (r>>8) + r;
  r = (r>>16) + r;
  r = (r>>32) + r;
  r = r & 127;
  return r;
}


static i64 rankc(rsstr *lf, i64 i, int c)
{
  i64 j,r;
  i64 i2;
  i64 i3;
  int lmb;
  bitvec_t x,m,*p;
  static bitvec_t masktbl[9] = {0,0x1111111111111111L,0x2222222222222222L,0x3333333333333333L,
              0x4444444444444444L,0x5555555555555555L,0x6666666666666666L,0x7777777777777777L,
              0x8888888888888888L};
  int c2;

  if (i < 0) return 0;

  c2 = convertchar(c);

  
  i3 = i;

  lmb = lf->lmb;
  r = getuint(lf->RL,(i>>logLB)*(SIGMADNA-1) + c2,lf->k);
  r += lf->RM[(i>>lmb)*(SIGMADNA-1) + c2];
  i2 = (i>>lmb) << lmb;
  p = &lf->BW[i2/(DD/logSIGMADNA)];

  for (j=0; j+(DD/logSIGMADNA)-1 <= (i-i2); j+=(DD/logSIGMADNA)) {
    x = (*p++) ^ masktbl[c2];
    x = (x | (x>>1) | (x>>2) | (x>>3)) & 0x1111111111111111L;
    r += (DD/logSIGMADNA) - POPCOUNT4(x);
  }
  x = (*p) ^ masktbl[c2];
  x = (x | (x>>1) | (x>>2) | (x>>3)) & 0x1111111111111111L;
  m = 0x1111111111111111L >> (((i-i2) - j + 1) * logSIGMADNA);
  x |= m;
  r += (DD/logSIGMADNA) - POPCOUNT4(x);

#if 0
  if (r != rankc_slow(lf,i3,c)) {
	  printf("rankc: error i=%d c=%d rank=%d (%d)\n",i3,c,r,rankc_slow(lf,i3,c));
  }
#endif
  return r;
}

static void make_tbl(rsstr *lf)
{
  i64 n;
  i64 i;
  i64 C[SIGMADNA];
  int c,mb;
  bitvec_t *BW;
  int k;

  n = lf->n;
  k = lf->k;
  BW = lf->BW;
  mb = 1 << lf->lmb;

  mymalloc(lf->RL,k*((n+LB-1)/LB+1)*(SIGMADNA-1));
  mymalloc(lf->RM,sizeof(lf->RM[0])*((n+mb-1)/mb+1)*(SIGMADNA-1));

  for (i=0;i<((n+LB-1)/LB)*(SIGMADNA-1);i++) {
    putuint(lf->RL,i,0,k);
  }
  for (i=0;i<((n+mb-1)/mb)*(SIGMADNA-1);i++) lf->RM[i] = 0;

  for (c = 0; c < SIGMADNA; c++) C[c] = 0;

  fprintf(stderr,"making tables...\n");
  for (i = 0; i < n; i++) {
    if (i % 1000000 == 0) {
      fprintf(stderr,"%ld\r",i/1000);
      fflush(stderr);
    }
    if (i % LB == 0) {
      for (c = 0; c < SIGMADNA-1; c++) {
        putuint(lf->RL,(i / LB)*(SIGMADNA-1) + c,C[c],k);
      }
    }
    if (i % mb == 0) {
      for (c = 0; c < SIGMADNA-1; c++) {
        lf->RM[(i / mb)*(SIGMADNA-1) + c]
            = (u16)(C[c] - getuint(lf->RL,(i / LB)*(SIGMADNA-1) + c,k));
      }
    }
    c = convertchar(access(lf,i));
    C[c]++;
  }
  
}

static i64 select_naive(rsstr *wt, i64 i, int c)
{
  i64 l,r,m;

  l = 0; r = wt->n-1;
  while (l <= r) {
    m = (l+r)/2;
    if (i <= wt->rank(wt,m,c)) r = m-1;  else l = m+1;
  }
  return l;
}

static i64 select_by_rank(rsstr *da, i64 ith, int c)
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

  int lmb;
  bitvec_t m;
  static bitvec_t masktbl[9] = {0,0x1111111111111111L,0x2222222222222222L,0x3333333333333333L,
              0x4444444444444444L,0x5555555555555555L,0x6666666666666666L,0x7777777777777777L,
              0x8888888888888888L};
  int c2;
  i64 i3;
  if (ith < 0) return 0;

  c2 = convertchar(c);

  
  i3 = ith;

  lmb = da->lmb;
  ll = 0;  rl = (da->n-1) >> logLB;
  pl = ll;  r2 = 0;
  while (ll <= rl) {
    ml = (ll+rl) >> 1;
    r = getuint(da->RL,ml*(SIGMADNA-1)+c2, da->k);
    if (r < ith) {
      pl = ml;
      ll = ml+1;
      r2 = r;
    } else {
      rl = ml-1;
    }
  }
  ith -= r2; // 大ブロック内のランク
  ll = (pl<<logLB)>>lmb;
  rl = min((((pl+1)<<logLB)>>lmb)-1, (da->n-1)>>lmb);

  pl = ll;  r2 = 0;
  l0 = ll;
  while (ll <= rl) {
    ml = (ll+rl) >> 1;
    r = da->RM[ml*(SIGMADNA-1)+c2];
    if (r < ith) {
      pl = ml;
      ll = ml+1;
      r2 = r;
    } else {
      rl = ml-1;
    }
  }
  ith -= r2; // 小ブロック内のランク

  p = pl << lmb; // ith を含む小ブロックの最初のビット位置

  if (p >= da->n) {
    p = da->n;
    goto end;
  }


#if 0
  while (1) {
    if (p >= da->n) break;
    if (access(da,p) == c) {
      ith--;
      if (ith == 0) break;
    }
    p++;
  }
#else
  q = &da->BW[p/(DD/logSIGMADNA)];

  while (1) {
    x = (*q++) ^ masktbl[c2];
    x = (x | (x>>1) | (x>>2) | (x>>3)) & 0x1111111111111111L;
    r = (DD/logSIGMADNA) - POPCOUNT4(x);
    if (r >= ith) break;
    ith -= r;
    p += DD/logSIGMADNA;
    if (p >= da->n) {
      p = da->n;  
      goto end;
    }
  }
  while (1) {
    if (p >= da->n) break;
    if (((x >> (DD-logSIGMADNA)) & ((1<<logSIGMADNA)-1)) == 0) {
      ith--;
      if (ith == 0) break;
    }
    x <<= logSIGMADNA;
    p++;
  }
  
#endif

end:
#if 0
  if (p != select_naive(da,i3,c)) {
	  printf("select: error i=%ld c=%d select=%ld (%ld)\n",i3,c,p,select_naive(da,i3,c));
  }
#endif

  return p;
  
}

////////////////////////////////////////
// returns select(rank(p-1,c)+1,c)
////////////////////////////////////////
static i64 succ(rsstr *da, i64 p, int c)
{
  bitvec_t x;
  bitvec_t *q;
  
  int r,s,i;
  static bitvec_t masktbl[9] = {0,0x1111111111111111L,0x2222222222222222L,0x3333333333333333L,
              0x4444444444444444L,0x5555555555555555L,0x6666666666666666L,0x7777777777777777L,
              0x8888888888888888L};
  int c2;
  i64 i3;

  c2 = convertchar(c);

  i3 = p;

  q = &da->BW[p/(DD/logSIGMADNA)];

  for (s=0; s<2; s++) {
    r = p % (DD/logSIGMADNA);
    x = (*q++) ^ masktbl[c2];
    x = (x | (x>>1) | (x>>2) | (x>>3)) & 0x1111111111111111L;
    x <<= (r*logSIGMADNA);
    while ((r < DD/logSIGMADNA) && (p < da->n)) {
      if (((x >> (DD-logSIGMADNA)) & ((1<<logSIGMADNA)-1)) == 0) break;
      x <<= logSIGMADNA;
      p++;
      r++;
    }
    if ((r < DD/logSIGMADNA) || (p >= da->n)) break;
  }
  if (s == 2) { // not found
    i64 d;
    if (i3 > 0) {
      d = da->rank(da, i3-1,c) + 1;
    } else {
      d = 1;
    }
    p = da->select(da, d, c);
  }

  return p;
  
}

i64 rsstr_succ2(rsstr *da, i64 s, i64 t, int c)
{
  bitvec_t x,x3;
  bitvec_t *q;
  
  int r,l,i;
  static bitvec_t masktbl[9] = {0,0x1111111111111111L,0x2222222222222222L,0x3333333333333333L,
              0x4444444444444444L,0x5555555555555555L,0x6666666666666666L,0x7777777777777777L,
              0x8888888888888888L};
  int c2,c3;
  i64 i3;

  c2 = convertchar(tolower(c));
  c3 = convertchar(toupper(c));

  i3 = s;

  q = &da->BW[s/(DD/logSIGMADNA)];

  for (l=0; l<2; l++) {
    r = s % (DD/logSIGMADNA);
    x = (*q) ^ masktbl[c2];
    x = (x | (x>>1) | (x>>2) | (x>>3)) & 0x1111111111111111L;
    x <<= (r*logSIGMADNA);
    x3 = (*q) ^ masktbl[c3];
    x3 = (x3 | (x3>>1) | (x3>>2) | (x3>>3)) & 0x1111111111111111L;
    x3 <<= (r*logSIGMADNA);
    q++;
    while ((r < DD/logSIGMADNA) && (s <= t)) {
      if (((x >> (DD-logSIGMADNA)) & ((1<<logSIGMADNA)-1)) == 0) break;
      if (((x3 >> (DD-logSIGMADNA)) & ((1<<logSIGMADNA)-1)) == 0) break;
      x <<= logSIGMADNA;
      x3 <<= logSIGMADNA;
      s++;
      r++;
    }
    if ((r < DD/logSIGMADNA) || (s > t)) break;
  }
  if (l == 2) { // not found
    s = -1;
  }

  return s;
  
}


i64 rsstr_makeindex(char *fname, int L)
{
  FILE *in,*out;
  bitvec_t *BW;
  char *fbw;
  int k,c,c2;
  i64 size;
  i64 n,i;
  rsstr *lf;

  lf = malloc(sizeof(*lf));
  if (lf == NULL) {
    printf("rsstr_makeindex: malloc failed.\n");
    exit(1);
  }

  size = 0;

  lf->lmb = blog(L);
  L = 1 << lf->lmb;
  lf->l = L;

  k = strlen(fname);
  mymalloc(fbw,k+5);
  sprintf(fbw,"%s.rs",fname);

  in = fopen(fname,"rb");
  if (in == NULL) {
    printf("rsstr_makeindex: cannot open %s\n",fbw);
    exit(1);
  }
  fseek(in,0,SEEK_END);
  n = ftell(in);
  fseek(in,0,SEEK_SET);
  printf("n=%ld\n",n);
  lf->n = n;

  mymalloc(BW,sizeof(*BW) * (n/(DD/logSIGMADNA)+1));
  lf->BW = BW;
  BW[n/(DD/logSIGMADNA)] = 0;

  fprintf(stderr,"packing...\n");
  for (i = 0; i < n; i++) {
    if (i % 1000000 == 0) {
      fprintf(stderr,"%ld\r",i/1000);
      fflush(stderr);
    }
    c = fgetc(in);
    c2 = convertchar(c);
    setbits(BW,i*logSIGMADNA,c2,logSIGMADNA);
  }
  fclose(in);

  out = fopen(fbw,"w");
  if (out == NULL) {
    printf("rsstr_makeindex: cannot open %s\n",fbw);
  }

  lf->k = k = (blog(n+1)+1+8-1)/8;

  make_tbl(lf);

//  writeint(1,ID_LF,out);
  writeint(1,k,out); /* #bytes of integer */
  writeint(k,n,out);
  writeint(sizeof(lf->l),lf->l,out);
  size += 1+k+sizeof(lf->l);

  for (i=0; i<n/(DD/logSIGMADNA)+1; i++) {
    writeuint(sizeof(*BW),BW[i],out);
    size += sizeof(*BW);
  }

  for (i=0;i<((n+LB-1)/LB)*(SIGMADNA-1);i++) {
    writeint(k,getuint(lf->RL,i,k),out);
    size += k;
  }
  for (i=0;i<((n+L-1)/L)*(SIGMADNA-1);i++) {
    writeint(sizeof(lf->RM[0]),lf->RM[i],out);
    size += sizeof(lf->RM[0]);
  }
  fclose(out);

  free(fbw);

  lf->size = size;

  lf->access = access;
  lf->rank = rankc;
//  lf->select = select_naive;
  lf->select = select_by_rank;
  lf->succ = succ;

  return size;
}


void rsstr_read(rsstr *lf, char *fname)
{
  int k,l,id;

  i64 psize;
  i64 n;
  uchar *p, *q;

  lf->mapbwt = mymmap(fname);
  if (lf->mapbwt->addr==NULL) {
    perror("rsstr_read: mmap2\n");
    exit(1);
  }
  p = q = (uchar *)lf->mapbwt->addr;
  psize = lf->mapbwt->len;

//  id = getuint(p,0,1);  p += 1;
//  if (id != ID_LF) {
//    printf("lf_dna_read: id = %d is not supported.\n",id);
//    exit(1);
//  }
  lf->k = k = getuint(p,0,1);  p += 1;
  lf->n = n = getuint(p,0,k);  p += k;
  lf->l = l = getuint(p,0,sizeof(lf->l));  p += sizeof(lf->l);
  lf->lmb = blog(l);
  if ((1 << lf->lmb) != l) {
    printf("L=%d must be a power of 2.\n",l);
    exit(1);
  }

  lf->BW = (bitvec_t *)p;  p += sizeof(*lf->BW) * (n/(DD/logSIGMADNA)+1);

  lf->RL = (uchar *)p;  p += ((n+LB-1)/LB)*(SIGMADNA-1) * k;
  lf->RM = (u16 *)p;  p += ((n+l-1)/l)*(SIGMADNA-1) * sizeof(lf->RM[0]);

  lf->size = psize;

  lf->access = access;
  lf->rank = rankc;
//  lf->rank = rankc_slow;
//  lf->select = select_naive;
  lf->select = select_by_rank;
  lf->succ = succ;


}

