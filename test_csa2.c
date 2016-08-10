/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#if 1
 #include <sys/timeb.h>
#else
 #include <sys/time.h>
 #include <sys/resource.h>
#endif
#include "csa.h"

#ifndef min
#define min(x,y) ((x)<(y)?(x):(y))
#endif

#ifndef _MYTIMESTRUCT_
#define _MYTIMESTRUCT_
#if 1
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
#else
typedef mytimestruct struct rusage;
void mygettime(mytimestruct *t)
{
  getrusage(RUSAGE_SELF,t);
}
double mylaptime(mytimestruct *before,mytimestruct *after)
{
  double t;
  t = after->ru_utime.tv_sec - before->ru_utime.tv_sec;
  t += (double)(after->ru_utime.tv_usec - before->ru_utime.tv_usec)
               /1000000;
  return t;
}
#endif
#endif

int strlen2(unsigned char *p)
{
  int l;
  l = 0;
  while (*p++ != 0) l++;
  return l;
}

void traverse2(CSA *csa, int max_depth, i64 l, i64 r, int depth, uchar *key, FILE *out)
{
  int f,c,cc;
  i64 ll, rr;
  i64 num;
  int C[4] = {'a', 'c', 'g', 't'};

#if 1
  if (l == r) {
    ll = csa->BW_LF(csa, l, &c);
    if (c == '$' || c == -1) return;
    fputc('(', out);
    fputc(c, out);
    key[depth] = c;
    if (depth+1 < max_depth) {
      traverse2(csa, max_depth, ll, ll, depth+1, key, out);
    } else {
      num = 1;
      fprintf(out, "%ld", num);
    }
    fputc(')', out);
    return;
  }
#endif

  for (c=0; c<4; c++) {
    ll = l;  rr = r;
    csa->searchsub(C[c], csa, &ll, &rr);
    if (ll <= rr) {
      fputc('(', out);
      fputc(C[c], out);
      key[depth] = c;
      if (depth+1 < max_depth) {
        traverse2(csa, max_depth, ll, rr, depth+1, key, out);
      } else {
        num = rr - ll + 1;
        fprintf(out, "%ld", num);
      }
      fputc(')', out);
    }
  }
}


unsigned char key[25600];
int main(int argc, char *argv[])
{
  i64 i,n;
  CSA SA;
  CSA SA2;
  mytimestruct before,after;
  double t;
  int k;
  FILE *in;
//  i64 count[10002];
  int C[4] = {'a', 'c', 'g', 't'};

   if (argc<2) {
      fprintf(stderr, "syntax: suftest file\n");
      return 1;
   }
   
   k = atoi(argv[1]);
//  in = fopen(argv[4], "rb");

   csa_read(&SA, 2, argv+2);
   n = SA.n;


#if 0
  for (i=0; i<n; i++) {
    int c1, c2;
    if (i != 20095) c1 = fgetc(in);
    c2 = SA.BW(&SA, i);
    if (c1 != c2) {
      printf("%d %c %c\n", i, c1, c2);
    }
  }
#endif

#if 0
  csa_read(&SA2, 2, argv+4);
  for (i=0; i<n; i++) {
    int r1, r2, c;
    if (i == 132) {
      printf("hoge\n");
    }
    c = SA2.BW(&SA2, i);
    printf("%d %c \n", i, c);
    for (c=0; c<4; c++) {
      r1 = SA.rankc(&SA, i, SA2.CtoA[C[c]]);
      r2 = SA2.rankc(&SA2, i, SA2.CtoA[C[c]]);
      if (r1 != r2) {
        printf("%d %c %d %d\n", i, C[c], r1, r2);
      }
    }
  }
#endif

{
  i64 th;
  uchar key[100];
  n = SA.n;
  traverse2(&SA, k, 0, n, 0, key, stdout);
  exit(0);
}

   return 0;
}
