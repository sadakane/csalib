/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "densearray.h"
#include "diskbuf.h"
#include "wavelet2.h"
#include "huffman.h"

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

static void writeint(int k,i64 x,FILE *f)
{
  int i;
  for (i=k-1; i>=0; i--) {
    fputc(x & 0xff,f);
    x >>= 8;
  }
}

static u64 readuint(int k, FILE *f)
{
  int j, c;
  u64 x;
  x = 0;
  for (j=0; j<k; j++) {
    c = fgetc(f);
    if (c == -1) {
      printf("readuint: c = %d\n", c);
      exit(1);
    }
    x += ((u64)c) << (j*8);
  }
  return x;
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




static i64 wt2_rank_access(wavelet2 *wt, i64 i, int *cc)
{
  int node;
  int c,n,m;
  i64 r;
  Huffman2 *huf;
  densearray **da;

  if (i < 0 || i >= wt->n) {
    printf("wt_access: i = %ld n = %ld\n", i, wt->n);
    exit(1);
  }

  huf = (Huffman2 *)wt->huffman;
  da = (densearray **)wt->da;
  n = huf->n;
  m = huf->m;
  node = n+m - 2;

  while (node >= n && i >= 0) {
    r = densearray_rank_and_bit(da[node-n], i, &c);
    if (c == 0) {
      i = (i+1) - r - 1;
      node = huf->left[node];
    } else {
      i = r - 1;
      node = huf->right[node];
    }
  }
  *cc = node;
  r = i+1;
  return r;
}

static int wt2_access(wavelet2 *wt, i64 i)
{
  int c;

  wt2_rank_access(wt, i, &c);
  return c;
}


static i64 wt2_rank(wavelet2 *wt, i64 i, int c)
{
  int node;
  u64 x;
  int d,n,m;
  Huffman2 *huf;
  densearray **da;

  if (i < 0 || i >= wt->n) {
    printf("wt_rank: i = %ld n = %ld\n", i, wt->n);
    exit(1);
  }

  huf = (Huffman2 *)wt->huffman;
  da = (densearray **)wt->da;
  n = huf->n;
  m = huf->m;
  node = n+m - 2;
  if (huf->clen[c] == 0) {
//    printf("wt_rank: no c = %d, i = %ld\n", c, i);
    return 0;
  }
  x = huf->code[c];

  while (node >= n && i >= 0) {
    d = (x >> (sizeof(x)*8-1)) & 1;
    if (d == 0) {
      i = densearray_rank0(da[node-n], i) - 1;
      node = huf->left[node];
    } else {
      i = densearray_rank(da[node-n], i) - 1;
      node = huf->right[node];
    }
    x <<= 1;
  }
  return i+1;
}

static i64 wt2_select_naive(wavelet2 *wt, i64 i, int c)
{
  i64 l,r,m;

  l = 0; r = wt->n-1;
  while (l <= r) {
    m = (l+r)/2;
    if (i <= wt->rank(wt,m,c)) r = m-1;  else l = m+1;
  }
//  if (l != wt2_select_fast(wt,i,c)) {
//    printf("select(i=%ld, c=%d) naive = %ld fast = %ld\n", i, c, l, wt2_select_fast(wt,i,c));
//  }
  return l;
}


static i64 wt2_select_fast(wavelet2 *wt, i64 i, int c)
{
  int node;
  u64 x;
  int d,n,m;
  Huffman2 *huf;
  densearray **da;
  int node2[64];
  char bit[64];
  int depth;
  i64 r,ii;

  if (i < 0 || i > wt->n) {
    printf("wt_select: i = %ld n = %ld\n", i, wt->n);
    exit(1);
  }
  if (i == 0) return -1;
  ii = i;

  huf = (Huffman2 *)wt->huffman;
  da = (densearray **)wt->da;
  n = huf->n;
  m = huf->m;
  node = n+m - 2;
  if (huf->clen[c] == 0) {
//    printf("wt_select: no c = %d, i = %ld\n", c, i);
    return -1;
  }
  x = huf->code[c];

  depth = 0;
  while (node >= n) {
    d = (x >> (sizeof(x)*8-1)) & 1;
    node2[depth] = node;
    bit[depth] = d;
    if (d == 0) {
      node = huf->left[node];
    } else {
      node = huf->right[node];
    }
    x <<= 1;
    depth++;
  }

  depth--;
  if (bit[depth] == 0) {
    r = da[node2[depth]-n]->n - da[node2[depth]-n]->m; // number of zeros
  } else {
    r = da[node2[depth]-n]->m; // number of ones
  }
  if (i > r) return wt->n; // no i-th character

  i--;
  while (depth >= 0) {
//    i = densearray_select(da[node-n], i+1, bit[depth]);
    i = densearray_select(da[node2[depth]-n], i+1, bit[depth]);
    depth--;
  }
#if 0
  if (i != wt2_select_naive(wt, ii, c)) {
    printf("wt2_select_fast: i = %ld ans = %ld naive = %ld\n",
      ii, i, wt2_select_naive(wt, ii, c));
    exit(1);
  }
#endif
  return i;
}




static void wt2_enum_sub(wavelet2 *wt, i64 l, i64 r, int node, int *num, int *C)
{
  int n;
  i64 rank, rank1, rank2;
  Huffman2 *huf;
  densearray **da;

  huf = (Huffman2 *)wt->huffman;
  da = (densearray **)wt->da;
  n = huf->n;

  if (node < n) {
//    putuint(C, *num, node, wt->k);
    C[*num] = node;
    *num = *num+1;
    return;
  }

  rank1 = densearray_rank(da[node-n], r);
  if (l > 0) rank2 = densearray_rank(da[node-n], l-1); else rank2 = 0;
  rank = rank1 - rank2;

  if (rank < r-l+1) { // there exists at least one 0
    wt2_enum_sub(wt, l-rank2, r-rank1, huf->left[node], num, C);
  }
  if (rank > 0) { // there exists at least one 1
    wt2_enum_sub(wt, rank2, rank1-1, huf->right[node], num, C);
  }
}

static int wt2_enum(wavelet2 *wt, i64 l, i64 r, int *C, i64 *freq)
{
  int node;
  int i,c,n,m, num;
  int *Ctmp;
  Huffman2 *huf;

  huf = (Huffman2 *)wt->huffman;

  n = huf->n;
  m = huf->m;
  node = n+m - 2;

  Ctmp = mymalloc(m * sizeof(*Ctmp));

  num = 0;
  wt2_enum_sub(wt, l, r, node, &num, Ctmp);

  for (i=0; i<num; i++) {
    C[i] = Ctmp[i];
    freq[i] = wt->rank(wt, r, C[i]);
    if (l>0) freq[i] -= wt->rank(wt, l-1, C[i]);
  }

  free(Ctmp);

  return num;
}




static void make_wavelet2_sub(wavelet2 *wt, i64 n, int k2, int depth, int node)
{
  i64 i,nl,nr;
  int c,d,m,w;
  FILE *in,*fl,*fr;
  u64 x;
  i64 stmp;
  int sigma;
  Huffman2 *h;
  densearray **da;

  sigma = wt->sigma;

  printf("sub node=%d depth=%d n=%ld\n",node,depth,n);
  fflush(stdout);
  h = (Huffman2 *)wt->huffman;
  da = (densearray **)wt->da;
  m = h->n;

  if (node < m) {
    remove_tmp(node);
    return;
  }
  w = sizeof(x)*8;

  da[node-m] = mymalloc(sizeof(densearray));
  densearray_construct_init(da[node-m], n);

  in = open_tmp(node);

  fl = create_tmp(h->left[node]);
  fr = create_tmp(h->right[node]);

  nl = nr = 0;
  for (i=0; i<n; i++) {
    c = readuint(k2, in);
    x = h->code[c];
    d = (x >> (w-1-depth)) & 1;
    densearray_construct_set(da[node-m], i, d);
    if (d == 0) {
      writeint(k2, c, fl);
      nl++;
    } else {
      writeint(k2, c, fr);
      nr++;
    }
  }
//  printf("nl = %ld nr = %ld\n", nl, nr);
  fclose(in);
  fclose(fl);  fclose(fr);
  remove_tmp(node);

  stmp = densearray_construct_end(da[node-m], wt->L, wt->opt);
  printf("densearray size=%ld (%1.3f bpc)\n",stmp,(double)stmp*8/n);
  fflush(stdout);
  wt->size += stmp;
  make_wavelet2_sub(wt, nl, k2, depth+1, h->left[node]);
  make_wavelet2_sub(wt, nr, k2, depth+1, h->right[node]);
}

static wavelet2 *make_wavelet2(int k2, int sigma, char *fwt, int opt)
{
  FILE *in,*out;
  i64 i,n;
  int root;
  u64 c;
  i64 *C;
//  double *freq;
  wavelet2 *wt;
  Huffman2 *huf;
  densearray **da;
  
  wt = mymalloc(sizeof(*wt));
  wt->sigma = sigma;
  wt->L = 512; //

  opt |= SDARRAY_RANK1;
  opt |= SDARRAY_SELECT1;
//    opt |= SDARRAY_LENHUF;
//    opt |= SDARRAY_RR;
//    opt |= SDARRAY_COMPRANK | SDARRAY_COMPPTR;
//    opt |= SDARRAY_SUC;

  wt->opt = opt;

  in = fopen(fwt,"rb");
  if (in == NULL) {
    printf("wt:make_wavelet: cannot open %s\n",fwt);
    exit(1);
  }
  fseek(in,0,SEEK_END);
  n = ftell(in);
  fseek(in,0,SEEK_SET);
  if (n % k2 != 0) {
    printf("??? file length = %ld k2 = %d\n", n, k2);
    exit(1);
  }
  n /= k2;
  printf("n=%ld\n",n);
  wt->n = n;
  wt->k = (blog(n+1)+1+8-1)/8;

 
  fprintf(stderr,"counting...\n");
  C = mymalloc(sizeof(*C) * sigma);
  for (i=0; i<sigma; i++) C[i] = 0;
  for (i = 0; i < n; i++) {
    if (i % 1000000 == 0) {
      fprintf(stderr,"%ld\r",i/1000);
      fflush(stderr);
    }
//    c = fgetc(in);
    c = readuint(k2, in);
    if (c >= sigma) {
      printf("wt:make_wavelet: c = %d >= sigma = %d\n",c, sigma);
      exit(1);
    }
    C[c]++;
  }
  fseek(in,0,SEEK_SET);

#if 0
  freq = mymalloc(sigma*sizeof(*freq));
  for (i=0; i<sigma; i++) {
    freq[i] = (double)C[i] / n;
//    printf("freq[%d] = %lf\n",i,freq[i]);
  }
#endif
redo:
//  wt->huffman = huf = MakeHuffman2Tree(sigma, freq);
//  wt->huffman = huf = MakeHuffman2Tree2(sigma, freq, 8);
  wt->huffman = huf = MakeHuffman2Tree2(sigma, C, 8);
  wt->m = huf->m;
  c = 0;
  for (i=0; i<sigma; i++) {
    if (huf->clen[i] > c) c = huf->clen[i];
  }
  if (c > Huffman_MAXLEN) {
    printf("maxlen = %d redo\n", c);
    for (i=0; i<sigma; i++) {
//      if (freq[i] > 0) freq[i] = (freq[i]+1)/2;
      if (C[i] > 0) C[i] = (C[i]+1)/2;
    }
    freeHuffman2(huf);
    goto redo;
  }

  free(C);
//  free(freq);

  root = huf->n + huf->m - 2;
  out = create_tmp(root);
  fprintf(stderr,"packing...\n");
  for (i = 0; i < n; i++) {
    if (i % 1000000 == 0) {
      fprintf(stderr,"%ld\r",i/1000);
      fflush(stderr);
    }
//    c = fgetc(in);
//    fputc(csa->CtoA[c],out);
    c = readuint(k2, in);
    writeint(k2, c, out);
  }
  fclose(out);
  
  wt->da = da = mymalloc(sizeof(densearray)*sigma);
  for (i=0; i<sigma; i++) da[i] = NULL;

  wt->size = 0;
  make_wavelet2_sub(wt, n, k2, 0, root);
  printf("wavelet: size = %ld (%1.3f bpc) \n",wt->size,
          (double)wt->size*8/n);
  fclose(in);

  return wt;
}

i64 wt2_makeindex(i64 k2, int sigma, char *fname, int opt)
{
  FILE *in,*out;
  char *fwt;
  i64 size;
  i64 i, d;
  wavelet2 *wt;
  int k, nn;
  Huffman2 *huf;
  densearray **da;

  d = strlen(fname);
  fwt = mymalloc(d+5);
  sprintf(fwt,"%s.wt",fname);

  densearray_make_selecttbl();

  wt = make_wavelet2(k2, sigma, fname, opt);
  huf = (Huffman2 *)wt->huffman;
  da = (densearray **)wt->da;
  k = wt->k;

  out = fopen(fwt,"w");
  if (out == NULL) {
    printf("wt_makeindex: cannot open %s\n",fwt);
  }

  size = 0;

  writeint(1,ID_WAVELET2,out);
  writeint(1,k,out); /* #bytes of integer */
  writeint(k,wt->n,out);
  writeint(k,wt->sigma-1,out);
  size += 1+1+k*2;

//  writeint(1,wt->id,out);
//  size += 1;

  Huffman2_write(wt->huffman, out);

#if 0
  nn = huf->n;
  for (i=nn; i<=2*nn-2; i++) {
    if (huf->left[i] < 2*nn) {
      size += densearray_write(da[i-nn], out);
    }
  }
#else
  nn = huf->m;
  for (i=0; i<nn-1; i++) {
    size += densearray_write(da[i], out);
  }
#endif

  fclose(out);

  free(fwt);

  return size;
}

void wt2_read(wavelet2 *wt, char *fname)
{
  int k;
  i64 psize1,bsize;
  i64 n,i,id;
  int nn;
  uchar *p, *q;
  Huffman2 *huf;
  densearray **da;
  MMAP *map;

  wt->map = map = mymmap(fname);
  if (map->addr==NULL) {
    perror("wt_read: mmap2\n");
    exit(1);
  }
  p = (uchar *)map->addr;
  q = p;

  id = getuint(p,0,1);  p += 1;
  if (id != ID_WAVELET2) {
    printf("wt2_read: id = %ld is not supported.\n",id);
    exit(1);
  }
  wt->k = k = getuint(p,0,1);  p += 1;
  wt->n = n = getuint(p,0,k);  p += k;
  wt->sigma = getuint(p,0,k)+1;  p += k;

//  id = getuint(p,0,1);  p += 1;
//  wt->id = id;

  wt->huffman = huf = Huffman2_read(&p);
  wt->m = huf->m;

  densearray_make_selecttbl();

  nn = huf->m;
  wt->da = da = mymalloc(sizeof(*da)*nn);
  for (i=0; i<nn; i++) da[i] = NULL;

  for (i=0; i<nn-1; i++) {
    da[i] = mymalloc(sizeof(densearray));
    densearray_read(da[i], &p);
  }
  wt->size = p - q;

  wt->access = (int (*)(wavelet2 *, i64))wt2_access;
  wt->rank = (i64 (*)(wavelet2 *, i64, int))wt2_rank;
  wt->rank_access = (i64 (*)(wavelet2 *, i64, int *))wt2_rank_access;
//  wt->select = (i64 (*)(wavelet2 *, i64, int))wt2_select_naive;
  wt->select = (i64 (*)(wavelet2 *, i64, int))wt2_select_fast;
  wt->enumerate = (int (*)(wavelet2 *, i64, i64, int *, i64 *))wt2_enum;

}

#ifdef MAIN
int main(int argc, char *argv[])
{
  i64 i,n;
  int c, cc, k2, sigma, j;
  wavelet2 wt;
  char fwt[100];
  FILE *in;
  i64 *C;
  i64 r;
  int opt;

  k2 = atoi(argv[1]);
  sigma = 1 << (8*k2);
  opt = 0;
  if (argc >= 4) {
    opt = atoi(argv[3]);
    printf("option = %d\n", opt);
  }
  wt2_makeindex(k2, sigma, argv[2], opt);

#if 1
  sprintf(fwt,"%s.wt",argv[2]);
  wt2_read(&wt, fwt);
  n = wt.n;

  C = mymalloc(sizeof(*C) * sigma);
  for (i=0; i<sigma; i++) C[i] = 0;

  in = fopen(argv[2], "r");
  for (i=0; i<n; i++) {
    c = readuint(k2, in);
    C[c]++;
    r = wt2.rank_access(&wt, i, &cc);
    if (c != cc) {
      printf("i=%ld c = %d cc = %d\n", i, c, cc);
    }
    if (r != C[c] || 0) {
      printf("i=%ld r = %d c = %d C = %ld\n", i, r, c, C[c]);
    }
#if 0
    for (cc=0; cc<sigma; cc++) {
      r = wt2.rank(&wt, i, cc);
      if (r != C[cc]) {
        printf("i=%ld c = %d r = %d C = %ld\n", i, cc, r, C[cc]);
      }
    }
#endif
  }
#endif
}

#endif
