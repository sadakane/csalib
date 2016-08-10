/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdlib.h>
#include <stdio.h>
#include "csa.h"

typedef struct {
  i64 l, r; // パタンのランクの範囲
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

void resizable_array_delete(resizable_array *array, int n)
{
  if (n > array->n) {
    printf("resizable_array_delete: n = %d > %d\n", n, array->n);
    exit(1);
  }
  array->n = n;
}

void resizable_array_free(resizable_array *array)
{
  free(array->list);
  free(array);
}


////////////////////////////////////////////////////////////////
// 左極大な頻出パタンを列挙する
// パタンの長さは maxlen 以下
////////////////////////////////////////////////////////////////
void frequent_leftmaximal(CSA *csa, i64 th, uchar *key, int len, int maxlen, int minlen, int error, int maxerror, resizable_array *array)
{
  int i,j,c,k;
  int maximal;
  i64 ll, rr, num;
  uchar head[SIGMA];
  i64 L[SIGMA], R[SIGMA];
  pattern newpat, *pat;
  resizable_array *newarray[SIGMA], *new2;
  i64 size[SIGMA], size2;
  int nl;
  int on;

  len++;  if (len > maxlen) return; // パタン長の最大値を超えたら終了
#if 0
  printf("sub |array| = %d\n", array->n);
  for (i=0; i<array->n; i++) {
    pat = (pattern *)resizable_array_get(array,i);
    printf("[%ld,%ld]\n", pat->l, pat->r);
  }
#endif
  on = array->n;


  for (c=0; c<csa->m; c++) {
    newarray[c] = resizable_array_init(array->n, 1, sizeof(pattern));
    size[c] = 0;
  }

  size2 = 0;  nl = 0;
  for (i=0; i<array->n; i++) {
    pat = (pattern *)resizable_array_get(array, i);
    k = csa->child_l(csa,pat->l,pat->r,head, L, R); // 各パタンに対し cP のランクを求める
    for (j=0; j<k; j++) {
      newpat.l = L[j];  newpat.r = R[j];
//      c = csa->head(csa,L[j]);
      c = csa->CtoA[head[j]];
      resizable_array_append(newarray[c], (char *)&newpat);
      size[c] += R[j] - L[j] + 1;
      size2 += R[j] - L[j] + 1;
      nl++;
    }
  }

  maximal = 1;
  for (c=0; c<csa->m; c++) {
    if (size[c] >= th) {
      *(key-1) = csa->AtoC[c];
      frequent_leftmaximal(csa, th, key-1, len, maxlen, minlen, error, maxerror, newarray[c]);
      maximal = 0; // cP が頻出なので，P は極大ではない
    }
  }

//  if (maximal == 1 && len > 1 && size2 >= th && error < maxerror) {
  if (len > 1 && size2 >= th && error < maxerror) { // ?P の頻度が閾値以上
    new2 = resizable_array_init(nl, 1, sizeof(pattern));
    *(key-1) = '*'; // wild card
    for (c=0; c<csa->m; c++) {
      if (size[c] > 0) {
        for (i=0; i<newarray[c]->n; i++) {
          resizable_array_append(new2, resizable_array_get(newarray[c],i));
        }
      }
    }
    frequent_leftmaximal(csa, th, key-1, len, maxlen, minlen, error+1, maxerror, new2);
//    maximal = 0;
    resizable_array_free(new2);
  }



  if (maximal && len >= minlen && size2 >= th) {
    printf("pattern = [%s] freq = %ld\n", key, size2);
  }

  for (c=0; c<csa->m; c++) resizable_array_free(newarray[c]);

}


#define MAXLEN 100
int main(int argc, char *argv[])
{
  i64 n;
  CSA SA;
  double t;
  i64 th;
  uchar key[MAXLEN+1];
  pattern pat, *p;
  resizable_array *array, *ans;
  int minlen, maxerror;

  if (argc<5) {
    fprintf(stderr, "syntax: %s threshold minlen maxerror {indexfiles}\n", argv[0]);
    return 1;
  }

  t = atof(argv[1]);
  minlen = atof(argv[2]);
  maxerror = atof(argv[3]);

  csa_read(&SA,argc-4, argv+4);
  n = SA.n;
  if (t < 1.0) th = n * t; else th = t;
  printf("n = %ld threshold = %ld\n", n, th);

  array = resizable_array_init(1000, 500, sizeof(pattern));

  key[MAXLEN] = 0;
  pat.l = 0;  pat.r = n;
  resizable_array_append(array, (char *)&pat);
  frequent_leftmaximal(&SA, th, key+MAXLEN, 0, MAXLEN, minlen, 0, maxerror, array);

  return 0;
}
