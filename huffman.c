/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdio.h>
#include <stdlib.h>
#include "huffman.h"
#include "heap.h"

#define mymalloc(p,n) {p = malloc((n)*sizeof(*p)); if ((p)==NULL) {printf("not enough memory in line %d\n",__LINE__); exit(1);};}

struct freqs {
  double freq;
  int idx;
};

struct freqs2 {
//  int freq;
  i64 freq;
  int idx;
};

static int tmpf_cmp(double *s, double *t, i64 w)
{
  if (*s < *t) return -1;
  if (*s > *t) return  1;
  return 0;
}

static int tmpf_cmp2(int *s, int *t, i64 w)
{
  if (*s < *t) return -1;
  if (*s > *t) return  1;
  return 0;
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

static void writeuint(int k,u64 x,FILE *f)
{
  int i;
  for (i=k-1; i>=0; i--) {
    fputc((int)(x & 0xff),f);
    x >>= 8;
  }
}

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

static Huffman *newHuffman(int n)
{
  int i;
  Huffman *p;
  mymalloc(p,1);
  p->n = n;
//  mymalloc(p->v,n);
  mymalloc(p->left,2*n+3);
  mymalloc(p->right,2*n+3);
//  mymalloc(p->clen,n*2);
//  mymalloc(p->code,n*2);
  mymalloc(p->clen,2*n+3);
  for (i=0; i<2*n+3; i++) p->clen[i] = 0;
  mymalloc(p->code,2*n+3);
  return p;
}

static Huffman2 *newHuffman2(int n, int m)
{
  int i;
  Huffman2 *p;
  mymalloc(p,1);
  p->n = n;
  p->m = m;
//  mymalloc(p->v,n);
  mymalloc(p->left,2*n+3);
  mymalloc(p->right,2*n+3);
//  mymalloc(p->clen,n*2);
//  mymalloc(p->code,n*2);
  mymalloc(p->clen,2*n+3);
  for (i=0; i<2*n+3; i++) p->clen[i] = 0;
  mymalloc(p->code,2*n+3);
  return p;
}


void freeHuffman(Huffman *p)
{
//  free(p->v);
  free(p->left);
  free(p->right);
  free(p->clen);
  free(p->code);
  free(p);
}

void freeHuffman2(Huffman2 *p)
{
//  free(p->v);
  free(p->left);
  free(p->right);
  free(p->clen);
  free(p->code);
  if (p->tbl_width > 0) free(p->tbl);
  free(p);
}

static void maketree(int n, int r, int d, Huffman *h, u64 c)
{
  h->clen[r] = d;
//  if (h->left[r] >= 0) {
  if (h->left[r] < 2*n) {
    maketree(n,h->left[r],d+1,h, c);
  } else {
    //h->v[r] = r;
    h->code[r] = c;
  }
//  if (h->right[r] >= 0) {
  if (h->right[r] < 2*n) {
    maketree(n,h->right[r],d+1,h, c + (1L<<(sizeof(c)*8-1-d)));
  } else {
    //h->v[r] = r;
    h->code[r] = c;
  }
}

static void maketree2(int n, int r, int d, Huffman2 *h, u64 c)
{
  h->clen[r] = d;
//  if (h->left[r] >= 0) {
  if (h->left[r] < 2*n) {
    maketree2(n,h->left[r],d+1,h, c);
  } else {
    //h->v[r] = r;
    h->code[r] = c;
  }
//  if (h->right[r] >= 0) {
  if (h->right[r] < 2*n) {
    maketree2(n,h->right[r],d+1,h, c + (1L<<(sizeof(c)*8-1-d)));
  } else {
    //h->v[r] = r;
    h->code[r] = c;
  }
}

void mkhufdecodetable(Huffman *h)
{
  u64 i;
  u64 v,w;

  for (i = 0; i < HUFTBLSIZ; i++) {
    w = i << (sizeof(u64)*8-HUFTBLWIDTH); // bit stream
    v = DecodeHuffman(h, w);
    if (h->clen[v] <= HUFTBLWIDTH) {
      h->tbl[i] = v;
    } else {
      h->tbl[i] = -1;
    }
  }
}

void mkhufdecodetable2(Huffman2 *h)
{
  u64 i;
  u64 v,w;
  int width;

//  if (h->n <= 256) width = 8;
//  else if (h->n <= 65536) width = 16;
//  else width = 0;

  width = h->tbl_width;
  mymalloc(h->tbl, 1<<width);

  for (i = 0; i < (1<<width); i++) {
    w = i << (sizeof(u64)*8-width); // bit stream
    v = DecodeHuffman2(h, w);
    if (h->clen[v] <= width) {
      h->tbl[i] = v;
    } else {
      h->tbl[i] = -1;
    }
  }
}

int DecodeHuffman(Huffman *h, u64 x)
{
  unsigned int i = sizeof(x)*8 - 1;
  int p = h->n * 2 - 2;
  while (1) {
    if ((x >> i) & 1) p = h->right[p]; else p = h->left[p];
    i--;
    //if (p < h->n) return h->v[p];
    if (p < h->n) return p;
  }
}

int DecodeHuffman2(Huffman2 *h, u64 x)
{
  unsigned int i = sizeof(x)*8 - 1;
  int p = h->n + h->m - 2;
  while (1) {
    if ((x >> i) & 1) p = h->right[p]; else p = h->left[p];
    i--;
    //if (p < h->n) return h->v[p];
    if (p < h->n) return p;
  }
}

int DecodeHuffman_tbl(Huffman *h, u64 x)
{
  u64 w;

  w = x >> (sizeof(u64)*8-HUFTBLWIDTH);
  if (h->tbl[w] >= 0) return h->tbl[w];
  return DecodeHuffman(h,x);
}

int DecodeHuffman2_tbl(Huffman2 *h, u64 x)
{
  u64 w;
  
  if (h->tbl_width > 0) {
    w = x >> (sizeof(u64)*8-h->tbl_width);
    if (h->tbl[w] >= 0) return h->tbl[w];
  }
  return DecodeHuffman2(h,x);
}

Huffman *MakeHuffmanTree(int n, double *freq)
{
  int i,j;
  Huffman *h;
  double *tmpf;
  int l,r;
  int m1,m2;
  int m; // ïpìxÇ™ 0 Ç≈ÇÕÇ»Ç¢ï∂éöÇÃêî

  m = 0;
  for (i=0; i<n; i++) {
    if (freq[i] > 0) m++;
  }
//  m = n; //

  h = newHuffman(n);
  mymalloc(tmpf,2*n+3);

//  for (i=0; i<n*2+3; i++) h->left[i] = h->right[i] = -1;
  for (i=0; i<n*2+3; i++) h->left[i] = h->right[i] = n*2;
  for (i=0; i<n; i++) tmpf[i] = freq[i];
//  for (i=n; i<n*2-1+2; i++) tmpf[i] = 0.0;
  for (i=n; i<=n+m-2; i++) tmpf[i] = 0.0;
//  tmpf[n*2-1] = tmpf[n*2] = 1000000.0;
  tmpf[n+m-1] = tmpf[n+m] = 1000000.0;

  l = 0; r = n-1;
//  for (j=0; j<n-1; j++) {
  for (j=0; j<m-1; j++) {
    printf("%d\r", j);  fflush(stdout);
//    m1 = n*2-1;  m2 = n*2;
    m1 = n+m-1;  m2 = n+m;
    for (i=l; i<=r; i++) {
      if ((tmpf[i] > 0) && (tmpf[i] < tmpf[m2])) m2 = i;
      if ((tmpf[i] > 0) && (tmpf[i] < tmpf[m1])) {m2 = m1; m1 = i;}
    }
    if (j == m-2) {
      r = n*2-3;
    }
    h->left[r+1] = m1;  h->right[r+1] = m2;
    tmpf[r+1] = tmpf[m1] + tmpf[m2];
    tmpf[m1] = tmpf[m2] = 0.0;
    r++;
  }
  maketree(n, r, 0, h, 0);

  free(tmpf);

  mkhufdecodetable(h);
  return h;
}

Huffman *MakeHuffmanTree2(int n, double *freq)
{
  int i,j;
  Huffman *h;
  struct freqs *tmpf, tmin1, tmin2, tnew;
  HEAP H;
  int l,r;
  int m1,m2;
  int m; // ïpìxÇ™ 0 Ç≈ÇÕÇ»Ç¢ï∂éöÇÃêî

  m = 0;
  for (i=0; i<n; i++) {
    if (freq[i] > 0) m++;
  }
//  m = n; //

  h = newHuffman(n);
  mymalloc(tmpf,2*n+3+1);

  for (i=0; i<n*2+3; i++) h->left[i] = h->right[i] = n*2;

  m = 0;
  for (i=0; i<n; i++) {
    if (freq[i] > 0) {
      tmpf[m+1].idx = i;
      tmpf[m+1].freq = freq[i];
      m++;
    }
  }

  heap_build(&H, m, (uchar *)tmpf, n+m, sizeof(struct freqs), sizeof(double),
            (void *)&(tmpf[0].freq) - (void *)&tmpf[0], (int (*)(uchar *, uchar *, i64))tmpf_cmp);

  l = 0; r = n-1;
  for (j=0; j<m-1; j++) {
//  	printf("j = %d  size = %d\n", j, H.size);
//  	heap_print(&H);
    heap_extract(&H, (uchar *)&tmin1);
    heap_extract(&H, (uchar *)&tmin2);
//  	printf("tmin1 (%d, %f) tmin2 (%d, %f)\n", tmin1.idx, tmin1.freq, tmin2.idx, tmin2.freq);
    if (j == m-2) {
      r = n*2-3;
    }
    r++;
    h->left[r] = tmin1.idx;  h->right[r] = tmin2.idx;
    tnew.idx = r;
    tnew.freq = tmin1.freq + tmin2.freq;
    heap_insert(&H, (uchar *)&tnew);
  }

  maketree(n, r, 0, h, 0);

  free(tmpf);

  mkhufdecodetable(h);
  return h;
}



Huffman2 *MakeHuffman2Tree(int n, double *freq)
{
  int i,j;
  Huffman2 *h;
  double *tmpf;
  int l,r;
  int m1,m2;
  int m; // ïpìxÇ™ 0 Ç≈ÇÕÇ»Ç¢ï∂éöÇÃêî

  m = 0;
  for (i=0; i<n; i++) {
    if (freq[i] > 0) m++;
  }

  h = newHuffman2(n,m);
  mymalloc(tmpf,2*n+3);

//  for (i=0; i<n*2+3; i++) h->left[i] = h->right[i] = -1;
  for (i=0; i<n*2+3; i++) h->left[i] = h->right[i] = n*2;
  for (i=0; i<n; i++) tmpf[i] = freq[i];
//  for (i=n; i<n*2-1+2; i++) tmpf[i] = 0.0;
  for (i=n; i<=n+m-2; i++) tmpf[i] = 0.0;
//  tmpf[n*2-1] = tmpf[n*2] = 1000000.0;
  tmpf[n+m-1] = tmpf[n+m] = 1000000.0;

  l = 0; r = n-1;
//  for (j=0; j<n-1; j++) {
  for (j=0; j<m-1; j++) {
//    m1 = n*2-1;  m2 = n*2;
    m1 = n+m-1;  m2 = n+m;
    for (i=l; i<=r; i++) {
      if ((tmpf[i] > 0) && (tmpf[i] < tmpf[m2])) m2 = i;
      if ((tmpf[i] > 0) && (tmpf[i] < tmpf[m1])) {m2 = m1; m1 = i;}
    }
//    printf("%d (%d,%d)\n", r+1, m1, m2);
    h->left[r+1] = m1;  h->right[r+1] = m2;
    tmpf[r+1] = tmpf[m1] + tmpf[m2];
    tmpf[m1] = tmpf[m2] = 0.0;
    r++;
  }
  maketree2(n, r, 0, h, 0);

  free(tmpf);

  h->tbl_width = HUFTBLWIDTH;
  mkhufdecodetable2(h);
  return h;
}

//Huffman2 *MakeHuffman2Tree2(int n, double *freq, int tbl_width)
Huffman2 *MakeHuffman2Tree2(int n, i64 *freq, int tbl_width)
{
  int i,j;
  Huffman2 *h;
  struct freqs2 *tmpf, tmin1, tmin2, tnew;
  HEAP H;
  int l,r;
  int m1,m2;
  int m; // ïpìxÇ™ 0 Ç≈ÇÕÇ»Ç¢ï∂éöÇÃêî

  m = 0;
  for (i=0; i<n; i++) {
    if (freq[i] > 0) m++;
  }

  h = newHuffman2(n,m);
  mymalloc(tmpf,2*n+3);

  for (i=0; i<n*2+3; i++) h->left[i] = h->right[i] = n*2;

  m = 0;
  for (i=0; i<n; i++) {
    if (freq[i] > 0) {
      tmpf[m+1].idx = i;
      tmpf[m+1].freq = (freq[i]+3)/4;
      m++;
    }
  }

//  heap_build(&H, m, (uchar *)tmpf, n+m, sizeof(struct freqs), sizeof(double),
//            (void *)&(tmpf[0].freq) - (void *)&tmpf[0], (int (*)(uchar *, uchar *, i64))tmpf_cmp);
  heap_build(&H, m, (uchar *)tmpf, n+m, sizeof(tmpf[0]), sizeof(tmpf[0].freq),
            (void *)&(tmpf[0].freq) - (void *)&tmpf[0], (int (*)(uchar *, uchar *, i64))tmpf_cmp2);

  l = 0; r = n-1;
  for (j=0; j<m-1; j++) {
//  	printf("j = %d  size = %d\n", j, H.size);
//  	heap_print(&H);
    heap_extract(&H, (uchar *)&tmin1);
    heap_extract(&H, (uchar *)&tmin2);

    r++;
    h->left[r] = tmin1.idx;  h->right[r] = tmin2.idx;
//    printf("%d (%d,%d)\n", r, tmin1.idx, tmin2.idx);
    tnew.idx = r;
    tnew.freq = tmin1.freq + tmin2.freq;
    heap_insert(&H, (uchar *)&tnew);
  }
  maketree2(n, r, 0, h, 0);

  free(tmpf);

  h->tbl_width = tbl_width;
  mkhufdecodetable2(h);
  return h;
}



void Huffman_write(Huffman *h, FILE *out)
{
  int k,n,i,len,d;
  u64 x;
//  printf("Huffman_write\n");
  writeuint(1,ID_HUFFMAN,out);
  n = h->n;
  k = (blog(2*n+4)+1+8-1)/8;
  writeuint(1,k,out);
  writeuint(k,n,out);

  for (i=n; i<=2*n-2; i++) {
//    printf("i=%d left=%d right=%d\n",i,h->left[i],h->right[i]);
    writeuint(k,h->left[i],out);
    writeuint(k,h->right[i],out);
  }
  for (i=0; i<n; i++) {
    len = h->clen[i];
    writeuint(1,len,out);
    if (len > 0) {
      d = (len+7)/8;
      x = h->code[i];
      x >>= (sizeof(x)-d)*8;
//    printf("i=%d len=%d d=%d code=%lx\n",i,len,d,x);
      writeuint(d,x,out);
    }
  }
}

void Huffman2_write(Huffman2 *h, FILE *out)
{
  int k,n,i,len,d;
  u64 x;
  int m;
//  printf("Huffman_write\n");
  writeuint(1,ID_HUFFMAN2,out);
  n = h->n;
  k = (blog(2*n+4)+1+8-1)/8;
  writeuint(1,k,out);
  writeuint(k,n,out);
  m = h->m;
  writeuint(k,m,out);

  for (i=n; i<=n+m-2; i++) {
//    printf("i=%d left=%d right=%d\n",i,h->left[i],h->right[i]);
    writeuint(k,h->left[i],out);
    writeuint(k,h->right[i],out);
  }
  for (i=0; i<n; i++) {
    len = h->clen[i];
    if (len > 0) {
      writeuint(1,len,out);
      d = (len+7)/8;
      x = h->code[i];
      x >>= (sizeof(x)-d)*8;
//    printf("i=%d len=%d d=%d code=%lx\n",i,len,d,x);
      writeuint(d,x,out);
    }
  }
}

Huffman *Huffman_read(uchar **map)
{
  int id;
  int k,n,i,len,d;
  u64 x;
  Huffman *h;
  uchar *p;
  
  p = *map;

  if ((id = getuint(p,0,1)) != ID_HUFFMAN) {
//    printf("Huffman_read: id = %d\n",id);
//    exit(1);
  }
  p += 1;
  k = getuint(p,0,1);  p += 1;
//  printf("Huffman_read k=%d\n",k);
  n = getuint(p,0,k);  p += k;
//  printf("Huffman_read n=%d\n",n);

  h = newHuffman(n);

  for (i=n; i<=2*n-2; i++) {
    h->left[i] = getuint(p,0,k);  p += k;
    h->right[i] = getuint(p,0,k);  p += k;
//    printf("i=%d left=%d right=%d\n",i,h->left[i],h->right[i]);
  }
  for (i=0; i<n; i++) h->code[i] = 0;
  for (i=0; i<n; i++) {
    len = h->clen[i] = getuint(p,0,1); p += 1;
    if (len > 0) {
      d = (len+7)/8;
      x = getuint(p,0,d);  p += d;
//    printf("i=%d len=%d d=%d code=%lx\n",i,len,d,x);
      x <<= (sizeof(x)-d)*8;
      h->code[i] = x;
    }
  }
  mkhufdecodetable(h);
  *map = p;
  return h;
}


Huffman2 *Huffman2_read(uchar **map)
{
  int id;
  int k,n,i,len,d;
  u64 x;
  Huffman2 *h;
  uchar *p;
  int m;
  
  p = *map;

  if ((id = getuint(p,0,1)) != ID_HUFFMAN2) {
//    printf("Huffman_read: id = %d\n",id);
//    exit(1);
  }
  p += 1;
  k = getuint(p,0,1);  p += 1;
//  printf("Huffman_read k=%d\n",k);
  n = getuint(p,0,k);  p += k;
//  printf("Huffman_read n=%d\n",n);
  m = getuint(p,0,k);  p += k;

  h = newHuffman2(n,m);

  for (i=0; i<n; i++) h->clen[i] = 0;
  for (i=n; i<=2*n-2; i++) {
    h->left[i] = h->right[i] = 2*n;
  }

  for (i=n; i<=n+m-2; i++) {
    h->left[i] = getuint(p,0,k);  p += k;
    h->right[i] = getuint(p,0,k);  p += k;
//    printf("i=%d left=%d right=%d\n",i,h->left[i],h->right[i]);
    if (h->left[i] < n) {
      h->clen[h->left[i]] = 1;
    }
    if (h->right[i] < n) {
      h->clen[h->right[i]] = 1;
    }
  }
  for (i=0; i<n; i++) h->code[i] = 0;
  for (i=0; i<n; i++) {
    if (h->clen[i] > 0) {
      len = h->clen[i] = getuint(p,0,1); p += 1;
      d = (len+7)/8;
      x = getuint(p,0,d);  p += d;
//      printf("i=%d len=%d d=%d code=%lx\n",i,len,d,x);
      x <<= (sizeof(x)-d)*8;
      h->code[i] = x;
    }
  }
  h->tbl_width = 0;
  mkhufdecodetable2(h);
  *map = p;
  return h;
}

i64 Huffman2_usedmemory(Huffman2 *h)
{
  i64 size;
  size = sizeof(*h);
//  size += h->n * sizeof(*h->v);
  size += (2*h->n+3) * sizeof(*h->left);
  size += (2*h->n+3) * sizeof(*h->right);
  size += (2*h->n+3) * sizeof(*h->clen);
  size += (2*h->n+3) * sizeof(*h->code);
  if (h->tbl_width > 0) size += (1 << h->tbl_width) * sizeof(*h->tbl);
  return size;
}
