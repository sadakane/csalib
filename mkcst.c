/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/

#include <stdio.h>
#include <stdlib.h>

#define display_progressbar(str,i,n) if (i % 10000000 == 0) {fprintf(stderr,"%s %ld/%ld                       \r",str,i/10000000,n/10000000);  fflush(stderr); }

#define SIGMA 256

#ifndef uchar
typedef unsigned char uchar;
#endif
typedef long i64;

int blog(unsigned long x) // [0,n] の数を格納するには blog(n)+1 ビット必要
{
int l;
  l = -1;
  while (x>0) {
    x>>=1;
    l++;
  }
  return l;
}

void *mymalloc(size_t n)
{
  void *p;

//  printf("allocating %ld bytes (%ld bytes available)\n",n,available_memory());

  p = malloc(n);
  if (p == NULL) {
    printf("malloc failed.\n");
    exit(1);
  }
  if (n == 0) {
    printf("warning: 0 bytes allocated p=%p\n",p);
  }

  return p;
}


typedef int suf_t;

int fast_hgt(suf_t *ISA, suf_t *SA, uchar *T, suf_t *H, i64 n)
{
  i64 p,q,h,i;
  FILE *f;

  for (i=0; i<=n; i++) ISA[SA[i]] = i;

  f = fopen("output.hgt","wb");

  h = 0;
  for (p=0; p<=n; p++) {
    i = ISA[p];
//    printf("isa[%ld] = %ld\n",p,i);
    if (i == 0) {
      H[i] = 0;
    } else {
      q = SA[i-1];
//      printf("h = %ld cmp T[%ld] T[%ld]\n",h,p+h,q+h);
      if (h == 0) {
        while ((p+h < n) && (q+h < n) && (T[p+h] == T[q+h])) h++;
      } else {
        h--;
        while ((p+h < n) && (q+h < n) && (T[p+h] == T[q+h])) h++;
      }
      H[i] = h;
//      printf("H[%ld] = %ld\n",i,h);
      i = h;
      while (i>0) {
        fputc('1',f);
        i--;
      }
      fputc('0',f);
    }
  }

  fclose(f);
  return 0;
}

int TB(uchar c, int r, int w)
{
  return (c >> (w-1 - r)) & 1;
}

int fast_hgt2(suf_t *ISA, suf_t *SA, uchar *T, suf_t *H, i64 n)
{
  i64 p,q,h,i,r;
  FILE *f;
  int w,j;
  int CtoB[SIGMA+1], BtoC[SIGMA+1], C[SIGMA+1];
  
  for (i=0; i<=SIGMA; i++) C[i] = 0;
  for (i=0; i<n; i++) C[T[i]+1] = 1;
  C[0] = 1;
  CtoB[-1+1] = 0;  BtoC[0] = -1;
  r = 1; // 出現する文字の数
  for (i=0; i<SIGMA; i++) {
    if (C[i+1]>0) { // 文字 i が T に出現
      CtoB[i+1] = r; // 文字 i が r にマップされる
      BtoC[r] = i; // r の元の文字は i ($ は -1)
      r++;
    }
  }
  w = blog(r)+1; // r 個の文字を格納するのに必要なビット数
  printf("w = %d\n",w);
  for (i=0; i<r; i++) {
    printf("%3d (%c) -> ",BtoC[i],isprint(BtoC[i])?BtoC[i]:' ');
    for (j=0; j<w; j++) putchar('0' + ((i>>(w-1-j))&1));
    printf("\n");
  }

  for (i=0; i<=n; i++) ISA[SA[i]] = i;

  f = fopen("output.hb","wb");

  h = 0;
  for (p=0; p<=n; p++) {
    int c1,c2;
    r = 0;
    i = ISA[p];
//    printf("isa[%ld] = %ld\n",p,i);
    if (i == 0) {
      H[i] = 0;
    } else {
      q = SA[i-1];
//      printf("h = %ld cmp T[%ld] T[%ld]\n",h,p+h,q+h);
      if (h == 0) {
        while ((p+h < n) && (q+h < n) && (T[p+h] == T[q+h])) h++;
        if (p+h < n) c1 = CtoB[T[p+h]+1]; else c1 = 0;
        if (q+h < n) c2 = CtoB[T[q+h]+1]; else c2 = 0;
        while (r < w && TB(c1,r,w) == TB(c2,r,w)) r++;
      } else {
        h--;
        while ((p+h < n) && (q+h < n) && (T[p+h] == T[q+h])) h++;
        if (p+h < n) c1 = CtoB[T[p+h]+1]; else c1 = 0;
        if (q+h < n) c2 = CtoB[T[q+h]+1]; else c2 = 0;
        while (r < w && TB(c1,r,w) == TB(c2,r,w)) r++;
      }
      H[i] = h * w + r;
//      printf("H[%ld] = %ld\n",i,h);
      r = H[i];
      while (r>0) {
        fputc('1',f);
        r--;
      }
      fputc('0',f);
    }
  }

  fclose(f);
  return 0;
}


int outputtree(uchar *fname, suf_t *H, i64 n)
{
  i64 i,j,h,h2;
  suf_t *s;
  suf_t *L,*R;
  FILE *f;

  s = malloc((n+2) * sizeof(*s));
  L = malloc((n+2) * sizeof(*L));
  R = malloc((n+2) * sizeof(*R));

  if (s == NULL || L == NULL || R == NULL) {
    printf("not enough memory. s=%p L=%p R=%p\n",s,L,R);
    exit(1);
  }

  for (i=0; i<=n+1; i++) L[i] = R[i] = 0;

  h = 0;
  s[h++] = -1;
  for (i=0; i<=n; i++) {
    h2 = h;
    while ((s[h-1] > H[i])) {
      h--;
    }
    R[i] = h2-h;
    //printf("R[%d] = %d\n",i,R[i]);
    if (s[h-1] == H[i]) {
      //R[i] = 0;
    }
    if (s[h-1] < H[i]) {
      s[h++] = H[i];
    }
  }
  h2 = h;
  while (s[h-1] >= 0) {
    h--;
  }
  R[n+1] = h2-h;
  //printf("R[%d] = %d\n",i,R[i]);

  h = 0;
  s[h++] = -1;
  for (i=n; i>=0; i--) {
    h2 = h;
    while ((s[h-1] > H[i])) {
      h--;
    }
    L[i+1] = h2-h;
    //printf("L[%d] = %d\n",i,L[i]);
    if (s[h-1] == H[i]) {
      //R[i] = 0;
    }
    if (s[h-1] < H[i]) {
      s[h++] = H[i];
    }
  }
  h2 = h;
  while (s[h-1] >= 0) {
    h--;
  }
  L[-1+1] = h2-h + 1;
  //printf("L[%d] = %d\n",i,L[i+1]);

  f = fopen(fname,"wb");

  fprintf(f,"(");
  for (i=0; i<=n; i++) {
    for (j=0; j<R[i]; j++) fprintf(f,")");
    for (j=0; j<L[i+1]; j++) fprintf(f,"(");
    fprintf(f,"()");
  }
  for (j=0; j<R[n+1]; j++) fprintf(f,")");

  fclose(f);

  free(R);
  free(L);
  free(s);
}

int outputtree2(uchar *fname, suf_t *H, i64 n)
{
  i64 i,j,h,h2;
  suf_t *s;
  suf_t *L,*R;
  FILE *f;

  s = malloc((n+2) * sizeof(*s));
  L = malloc((n+2) * sizeof(*L));
  R = malloc((n+2) * sizeof(*R));

  if (s == NULL || L == NULL || R == NULL) {
    printf("not enough memory. s=%p L=%p R=%p\n",s,L,R);
    exit(1);
  }

  for (i=0; i<=n+1; i++) L[i] = R[i] = 0;

  h = 0;
  s[h++] = -1;
  for (i=0; i<=n; i++) {
    h2 = h;
    while ((s[h-1] > H[i])) {
      h--;
    }
    R[i] = h2-h;
    //printf("R[%d] = %d\n",i,R[i]);
    if (s[h-1] == H[i]) {
      //R[i] = 0;
    }
    if (s[h-1] < H[i]) {
      s[h++] = H[i];
    }
  }
  h2 = h;
  while (s[h-1] >= 0) {
    h--;
  }
  R[n+1] = h2-h;
  //printf("R[%d] = %d\n",i,R[i]);

  h = 0;
  s[h++] = -1;
  for (i=n; i>=0; i--) {
    h2 = h;
    while ((s[h-1] > H[i])) {
      h--;
    }
    L[i+1] = h2-h;
    //printf("L[%d] = %d\n",i,L[i]);
    if (s[h-1] == H[i]) {
      //R[i] = 0;
    }
    if (s[h-1] < H[i]) {
      s[h++] = H[i];
    }
  }
  h2 = h;
  while (s[h-1] >= 0) {
    h--;
  }
  L[-1+1] = h2-h + 1;
  //printf("L[%d] = %d\n",i,L[i+1]);

  f = fopen(fname,"wb");

  fprintf(f,"(");
  for (i=0; i<=n; i++) {
//    for (j=0; j<R[i]; j++) fprintf(f,")");
    for (j=0; j<L[i+1]; j++) fprintf(f,"(");
    fprintf(f,")");
  }
//  for (j=0; j<R[n+1]; j++) fprintf(f,")");

  fclose(f);

  free(R);
  free(L);
  free(s);
}

void bw_to_psi(i64 *n0, i64 *last0, suf_t **psi0, char *fbw, char *flst)
{
  FILE *in;
  i64 i,j,n,last;
  suf_t *psi;
  i64 C[SIGMA];
  i64 w,c,k;

  in = fopen(flst,"r");
  if (in == NULL) {
    perror("bw_to_psi:");  exit(1);
  }
  fscanf(in,"%ld",&last);
  fclose(in);

  for (c=0; c<SIGMA; c++) {
    C[c] = 0;
  }

  in = fopen(fbw,"r");
  if (in == NULL) {
    perror("bw_to_psi:");  exit(1);
  }
  n = 0;
  while (1) {
    display_progressbar("reading ",n,(long int)0);
    c = fgetc(in);
    if (c == EOF) break;
    C[c]++;
    n++;
  }
  rewind(in);
  printf("n = %ld last = %ld\n",n,last);

  psi = (suf_t *)mymalloc((n+1)*sizeof(suf_t));

  for (j=1,c=0; c<SIGMA; c++) {
    k = C[c];
    C[c] = j;
    j += k;
  }

  for (i = 0; i<=n; i++) {
    display_progressbar("computing psi ",i,n);
    if (i == last) {
      psi[0] = i;
    } else {
      c = fgetc(in);
//      printf("psi[%ld] = %ld\n",C[c],i);
      psi[C[c]++] = i;
    }
  }
  fclose(in);
  *n0 = n;
  *last0 = last;
  *psi0 = psi;
}

void psi_to_sa(i64 n, i64 last, suf_t *psi, suf_t *sa)
{
  i64 i,j,v;

  j = last;
  for (i=0; i<=n; i++) {
    display_progressbar("making sa ",i,n);
    v = psi[j];
    sa[j] = i;
    j = v;
  }
  sa[0] = n;
}

int main(int argc, char *argv[])
{
  suf_t *sa,*isa, *H;
  i64 n, last;
  i64 i;
  uchar *T;
  FILE *in;

  if (argc < 3) {
    printf("mkcst file file.bw file.lst\n");
    return 0;
  }

  bw_to_psi(&n, &last, &sa, argv[2], argv[3]);
  psi_to_sa(n, last, sa, sa);
#if 1
  for (i=0; i<=n; i++) {
    printf("sa[%ld] = %ld\n",i,sa[i]);
  }
#endif
  isa = (suf_t *)mymalloc((n+2)*sizeof(suf_t));
  T = (uchar *)mymalloc((n+1)*sizeof(uchar));
  H = (suf_t *)mymalloc((n+1)*sizeof(suf_t));
  in = fopen(argv[1],"r");
  if (in == NULL) {printf("cannot open %s\n",argv[1]);  exit(1);}
  for (i=0; i<n; i++) T[i] = fgetc(in);
  fclose(in);

  fast_hgt(isa, sa, T, H, n);
#if 1
  for (i=0; i<=n; i++) {
    printf("H[%ld] = %d\n",i,H[i]);
  }
#endif
  outputtree("output.bp", H, n);

  fast_hgt2(isa, sa, T, H, n);
#if 1
  for (i=0; i<=n; i++) {
    printf("H[%ld] = %d\n",i,H[i]);
  }
#endif
  outputtree2("output.bp2", H, n);

}
