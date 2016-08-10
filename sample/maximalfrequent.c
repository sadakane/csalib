/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdlib.h>
#include <stdio.h>
#include "csa.h"

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

int valid_utf8(uchar *p, int len, uchar **b, uchar **e)
{
  int c, k, i;

  *b = *e = NULL;
  while (len > 0) {
    c = *p;
    k = utf8_len(c);
    if (k > 0) break;
    p++;
    len--;
  }
  if (len == 0) return 0;
  *b = p;
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
    *e = p-1;
  }
  return 1;
}

typedef struct {
  i64 l, r; // パタンのランクの範囲
  int len;  // パタンの長さ
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

int sign(i64 x)
{
  if (x < 0) return -1;
  if (x > 0) return 1;
  return 0;
}

int cmppattern(const void *p, const void *q) // 範囲の比較に使う関数
{
  pattern *pp, *qq;
  pp = (pattern *)p;  qq = (pattern *)q;
  if (pp->l != qq->l) return sign(pp->l - qq->l); // まず範囲の左端の小さい順にソート
  if (pp->r != qq->r) return sign(qq->r - pp->r); // 左端が同じときは左端の大きい順にソート
  return pp->len - qq->len;                 // 範囲が同じときは短い順にソート
}


////////////////////////////////////////////////////////////////
// 左極大なパタンの中から右極大なものを見つける
// パタンの長さは maxlen 以下
////////////////////////////////////////////////////////////////
resizable_array *find_rightmaximal(CSA *csa, resizable_array *array, int maxlen)
{
  int i, j, d;
  pattern *pat;
  pattern *stack;
  resizable_array *ans;
  
  if ((stack = malloc((maxlen+2) * sizeof(pattern))) == NULL) {
    printf("find_rightmaximal: not enough memory.\n");
    exit(1);
  }

  qsort(array->list, array->n, array->m, cmppattern); // 左極大パタンをソート

  ans = resizable_array_init(1000, 500, sizeof(pattern)); // 答えを入れる配列

  d = 0;
  stack[d].l = 0;  stack[d].r = csa->n;  stack[d].len = 0; // スタックの底に範囲全体を入れておく
  for (i=0; i<array->n; i++) {
    pat = (pattern *)resizable_array_get(array, i);
    if (stack[d].l <= pat->l && pat->r <= stack[d].r) { // 入れ子になっていればスタックに積む
      d++;
      if (d > maxlen+1) printf("??? d = %d\n",d); // 起こらないはず
      stack[d].l = pat->l;  stack[d].r = pat->r;  stack[d].len = pat->len;
    } else if (stack[d].r < pat->l) { // スタックの一番上と重なりがないとき
      resizable_array_append(ans, (char *)&stack[d]); // スタックの一番上の要素は極大
      while (stack[d].r < pat->l) { // pat と入れ子になるまでスタックの中身を取り出す
        d--;
        if (d < 0) printf("????? d = %d\n",d); // 起こらないはず
      }
    } else { // 範囲に重なりがあるとき（起こらないはず）
      printf("??? stack top [%ld,%ld] pat [%ld,%ld]\n",stack[d].l,stack[d].r,pat->l,pat->r);
    }
  }
  
  free(stack);
  return ans;
}

////////////////////////////////////////////////////////////////
// 左極大な頻出パタンを列挙する
// パタンの長さは maxlen 以下
////////////////////////////////////////////////////////////////
void frequent_leftmaximal(CSA *csa, pattern pat, i64 th, uchar *key, int maxlen, resizable_array *array)
{
  int i,c,k,len;
  int maximal;
  i64 ll, rr, num;
  uchar head[SIGMA];
  i64 L[SIGMA], R[SIGMA];
  pattern newpat;

  len = pat.len;
  len++;  if (len > maxlen) return; // パタン長の最大値を超えたら終了

  k = csa->child_l(csa,pat.l,pat.r,head, L, R); // パタン P の左側に出現する文字集合を求める

  maximal = 1;
  for (i=0; i<k; i++) {
    c = head[i];
    ll = L[i];  rr = R[i];
    if (rr - ll + 1 >= th) { // パタン cP の頻度が閾値以上ならば探索を続ける
      *(key-1) = c;
      newpat.l = ll;  newpat.r = rr;  newpat.len = len;
      frequent_leftmaximal(csa, newpat, th, key-1, maxlen, array);
      maximal = 0; // cP が頻出なので，P は極大ではない
    }
  }

  if (maximal) {
    resizable_array_append(array, (char *)&pat);
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
  pattern pat, *p;
  resizable_array *array, *ans;
  int i;
  uchar *str;

  if (argc<3) {
    fprintf(stderr, "syntax: %s threshold {indexfiles}\n", argv[0]);
    return 1;
  }

  t = atof(argv[1]);

  csa_read(&SA,argc-2, argv+2);
  n = SA.n;
  if (t < 1.0) th = n * t; else th = t;
  printf("n = %ld threshold = %ld\n", n, th);

  array = resizable_array_init(1000, 500, sizeof(pattern));

  key[MAXLEN] = 0;
  pat.l = 0;  pat.r = n;  pat.len = 0;
  frequent_leftmaximal(&SA, pat, th, key+MAXLEN, MAXLEN, array);

  ans = find_rightmaximal(&SA, array, MAXLEN);

  if ((str = malloc((MAXLEN+1) * sizeof(*str))) == NULL) {
    printf("not enough memory.\n");
    exit(1);
  }

  for (i=0; i<ans->n; i++) {
    uchar *b, *e;
    p = (pattern *)resizable_array_get(ans,i);
    SA.substring(str, &SA, p->l, p->len);
    valid_utf8(str, p->len, &b, &e);
    if (b != NULL && e != NULL) {
      e[1] = 0;
      printf("[%ld,%ld] len = %d [%s]\n",p->l,p->r,e-b+1,b);
    }
  }
  printf("number of left-maximal patterns = %d\n",array->n);
  printf("number of maximal patterns = %d\n",ans->n);

  return 0;
}
