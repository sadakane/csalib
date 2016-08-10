/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdlib.h>
#include <stdio.h>
#include "csa.h"

int gap;

int utf_head(int c)
{
  return (0 <= c && c <= 0x7f);
}

int utf8_len(int c)
{
  int k;
  if (0 <= c && c <= 0x7f) k = 1;
  else if (0xc2 <= c && c <= 0xdf) k = 2;
  else if (0xe0 <= c && c <= 0xef) k = 3;
  else if (0xf0 <= c && c <= 0xf7) k = 4;
  else if (0xf8 <= c && c <= 0xfb) k = 5;
  else if (0xfc <= c && c <= 0xfd) k = 6;
  else k = -1;
  return k;
}

int valid_utf8(uchar *p, int len)
{
  int c, k, i;

  while (len > 0) {
    c = *p++;
    k = utf8_len(c);
    if (k <= 0) return 0;
    if (k > len) return 0;

    for (i=1; i<k; i++) {
      c = *p++;
      if (!(0x80 <= c && c <= 0xbf)) return 0;
    }
    len -= k;
  }
  return 1;
}

int i64cmp(const void *p, const void *q)
{
  i64 x, y;
  x = *(i64 *)p;  y = *(i64 *)q;
  if (x < y) return -1;
  if (x > y) return 1;
  return 0;
}


typedef struct {
  i64 l, r; // パタンのランクの範囲
//  int len;  // パタンの長さ
} pattern;

typedef struct {
  int n; // 要素数
  int m; // 1要素のバイト数
  int size; // 最大要素数
  int step; // サイズ更新時に増やす要素数
  void *list; // 要素を格納する配列
} resizable_array;

resizable_array *resizable_array_init(int initial_size, int step, int m)
{
  resizable_array *array;
  if ((array = malloc(sizeof(resizable_array))) == NULL) {
    printf("init_resizable_array: not enough memory.\n");
    exit(1);
  }
  array->n = 0;
  array->m = m;
  array->size = initial_size;
  array->step = step;
  if ((array->list = malloc(initial_size*m)) == NULL) {
    printf("init_resizable_array: not enough memory.\n");
    exit(1);
  }
  return array;
}

void resizable_array_enlarge(resizable_array *array)
{
  array->size += array->step;
  array->list = realloc(array->list, array->size * array->m);
  if (array->list == NULL) {
    printf("resizable_array_enlarge: not enough memory.\n");
    exit(1);
  }
}

char *resizable_array_get(resizable_array *array, int i)
{
  char *q;
  q = (char *)array->list;
  q += i * array->m;
  return q;
}

void resizable_array_append(resizable_array *array, char *p)
{
  int i;
  char *q;

  if (array->n == array->size) resizable_array_enlarge(array);

  q = resizable_array_get(array, array->n);
  for (i=0; i < array->m; i++) *q++ = *p++;
  array->n++;
}

void resizable_array_truncate(resizable_array *array)
{
  array->size = array->n;
  array->list = realloc(array->list, array->size * array->m);
  if (array->list == NULL) {
    printf("resizable_array_truncate: not enough memory.\n");
    exit(1);
  }
}

void resizable_array_free(resizable_array *array)
{
  free(array->list);
  free(array);
}


void find_proximity_patterns_sub(CSA *csa, i64 l, i64 r, int gap, i64 minfreq, resizable_array *array)
{
  int i,c,k;
  i64 ll, rr, num;
  uchar head[SIGMA];
  i64 L[SIGMA], R[SIGMA];
  pattern pat;
  
  k = csa->child_l(csa,l,r,head, L, R);
  for (i=0; i<k; i++) {
    c = head[i];
    ll = L[i];  rr = R[i];
    if (rr - ll + 1 >= minfreq || 1) {
      pat.l = ll;  pat.r = rr;
      resizable_array_append(array, (char *)&pat);
    }
    if (gap > 0) find_proximity_patterns_sub(csa, ll, rr, gap-1, minfreq, array);
  }
}

void find_proximity_patterns2(CSA *csa, i64 l, i64 r, int gap, i64 minfreq)
{
  resizable_array *array;

  array = resizable_array_init(1000, 500, sizeof(pattern));
  find_proximity_patterns_sub(csa, l, r, gap, minfreq, array);
  printf("freq = %ld #patterns = %ld\n", r-l+1, array->n);
  resizable_array_free(array);
}

void find_proximity_patterns(CSA *csa, i64 l, i64 r, int gap, i64 minfreq)
{
  i64 *A;
  i64 i,k;
  i64 rank;
  i64 size;
  
  size = (r-l+1)*gap;
  A = malloc(size*sizeof(*A));
  if (A == NULL) {
    printf("find_proximity_patterns: not enough memory num = %ld gap = %ld\n", r-l+1, gap);
    exit(1);
  }
  for (i=0; i<r-l+1; i++) {
    rank = l+i;
    for (k=0; k<gap; k++) {
      rank = csa->LF(csa, rank);
      A[i*gap+k] = rank;
    }
  }
  qsort(A, size, sizeof(*A), i64cmp);
  
  free(A);
}

////////////////////////////////////////////////////////////////
// 左極大な頻出パタンを列挙する
// パタンの長さは maxlen 以下
// テキストがUnicode (utf8) と仮定し，文字の切れ目が正しいものだけ求める
////////////////////////////////////////////////////////////////
void frequent_leftmaximal(CSA *csa, i64 l, i64 r, i64 th, uchar *key, int len, int maxlen)
{
  int i,c,k;
  int maximal;
  i64 ll, rr, num;
  uchar head[SIGMA];
  i64 L[SIGMA], R[SIGMA];

  len++;  if (len > maxlen) return; // パタン長の最大値を超えたら終了

  k = csa->child_l(csa,l,r,head, L, R); // パタン P の左側に出現する文字集合を求める

  maximal = 1;
  for (i=0; i<k; i++) {
    c = head[i];
    ll = L[i];  rr = R[i];
    if (rr - ll + 1 >= th) { // パタン cP の頻度が閾値以上ならば探索を続ける
      *(key-1) = c;
      frequent_leftmaximal(csa, ll, rr, th, key-1, len, maxlen);
      maximal = 0; // cP が頻出なので，P は極大ではない
      // 注: cP はvalidではないとき，cPもPも表示されなくなってしまう
    }
  }

  if (maximal) {
    printf("[%ld,%ld] key = %s\n",l,r,key);
    find_proximity_patterns2(csa, l, r, gap, th);
  }
}

#define MAXLEN 100
int main(int argc, char *argv[])
{
  i64 n;
  CSA SA;
  double t;
  i64 th;
  uchar key[MAXLEN+1];

  if (argc<4) {
    fprintf(stderr, "syntax: %s threshold gap {indexfiles}\n", argv[0]);
    return 1;
  }

  t = atof(argv[1]);
  gap = atoi(argv[2]);

  csa_read(&SA,argc-3, argv+3);
  n = SA.n;
  if (t < 1.0) th = n * t; else th = t;
  printf("n = %ld threshold = %ld gap = %d\n", n, th, gap);

  key[MAXLEN] = 0;
  frequent_leftmaximal(&SA, 0, n, th, key+MAXLEN, 0, MAXLEN);

  return 0;
}
