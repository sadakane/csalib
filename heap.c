/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdio.h>
#include <stdlib.h>
#include "csa.h"
#include "heap.h"

#define LEFT(i) (2*i)
#define RIGHT(i) (2*(i)+1)
#define PARENT(i) (i/2)

static int default_cmp(uchar *s, uchar *t, i64 key_w)
{
  i64 i;
  for (i=0; i<key_w; i++) {
    if (s[i] < t[i]) return -1;
    if (s[i] > t[i]) return  1;
  }
  return 0;
}

static void copy(uchar *t, uchar *s, i64 item_w)
{
  i64 i;
  for (i=0; i<item_w; i++) *t++ = *s++;
}



void heapify(HEAP *H, i64 i)
{
  i64 l, r, largest, heap_size;
  uchar *A;
  i64 w,o,k;
  uchar tmp[64];

  A = H->A;  heap_size = H->size;
  w = H->item_w;  o = H->key_ofs;  k = H->key_w;
  if (w > 64) {
    printf("heapify: w = %ld\n",w);  exit(1);
  }

  while (1) {
    l = LEFT(i);   r = RIGHT(i);
    if (l <= heap_size && H->cmp(&A[l*w+o],&A[i*w+o],k) < 0) largest = l; // A[i] と左の子で小さい
    else  largest = i;                                                 // 方をlargestに
    if (r <= heap_size && H->cmp(&A[r*w+o],&A[largest*w+o],k) < 0)        // 右の子の方が小さい
      largest = r;
    if (largest == i) break;
    copy(tmp,&A[i*w],w);
    copy(&A[i*w],&A[largest*w],w);
    copy(&A[largest*w],tmp,w);           // A[i]を子供と入れ替える
    i = largest;
  }
}

void heap_build(HEAP *H, i64 n, uchar *A, i64 max_n, i64 item_w, i64 key_w, i64 key_ofs, int (*cmp)(uchar *, uchar *, i64))
{
  i64 i;
  H->max_n = max_n;
  H->size = n;
  H->item_w = item_w;  H->key_w = key_w;  H->key_ofs = key_ofs;
  H->A = A;
  if (cmp != NULL) H->cmp = cmp;  else H->cmp = default_cmp;

  for (i = H->max_n/2; i >= 1; i--) {
    heapify(H,i);
  }
}

void heap_print(HEAP *H)
{
  i64 i,j;
  for (i=1; i<=H->size; i++) {
    printf("%ld: ",i);
    for (j=0; j<H->item_w; j++) printf("%02x ",H->A[i*H->item_w+j]);
    printf("\n");
  }
}

#if 0
void heapsort(i64 n, data *A)
{
  i64 i;
  data tmp;
  HEAP H;
  uchar tmp[64];
  i64 w;
  
  w = A->item_w;

  BUILD_HEAP(&H,n,A,n);
  for (i = n; i >= 2; i--) {
    copy(tmp,&A[1*w],w);
    copy(&A[1*w],&A[i*w],w);
    copy(&A[i*w],tmp,w);         // 根と最後の要素を交換
    H.size = H.size - 1;
    HEAPIFY(&H,1);
  }
}
#endif

void heap_extract(HEAP *H, uchar *MIN)  // O(lg n) 時間
{
  uchar *A;
  i64 w;

  A = H->A;
  w = H->item_w;

  if (H->size < 1) {
    printf("ERROR ヒープのアンダーフロー\n");
    exit(1);
  }
  copy(MIN,&A[1*w],w);
  copy(&A[1*w],&A[H->size*w],w);
  H->size = H->size - 1;
  heapify(H,1);
  return;
}

void heap_insert(HEAP *H, uchar *key)  // O(lg n) 時間
{
  i64 i,w,o,k;
  uchar *A;

  A = H->A;
  w = H->item_w;  o = H->key_ofs;  k = H->key_w;

  H->size = H->size + 1;
  if (H->size > H->max_n) {
    printf("ERROR ヒープのオーバーフロー\n");
    exit(1);
  }

  i = H->size;
  while (i > 1 && H->cmp(&A[PARENT(i)*w+o],key+o,k) > 0) {
    copy(&A[i*w],&A[PARENT(i)*w],w);
    i = PARENT(i);
  }
  copy(&A[i*w],key,w);
}

void heap_delete(HEAP *H, i64 i)  // O(lg n) 時間
{
  i64 w;
  uchar *A;

  A = H->A;

  if (i < 1 || i > H->size) {
    printf("ERROR 範囲エラー (%d,%d)\n",i,H->size);
    exit(1);
  }

  while (i > 1) {
    copy(&A[i*w],&A[PARENT(i)*w],w);
    i = PARENT(i);
  }
  copy(&A[1*w],&A[H->size*w],w);
  H->size = H->size - 1;
  heapify(H,1);
}

#if 0
int main(int argc, char *argv[])
{
  data A[14+1] = {-1,27,17,3,16,13,10,1,5,7,12,4,8,9,0};
  HEAP H;
  int i,n;
  

  //  HEAPSORT(14,A);
  //  for (i=1;i<=14;i++) printf("%d ",A[i]);
  //  printf("\n");

  heap_build(&H,14,A,14);
  n = H.size+1;
  for (i=1; i<=n; i++) {
    printf("min = %d\n",heap_extract(&H));
    for (i=1;i<=H.size;i++) printf("%d ",A[i]);
    printf("\n");
  }

}
#endif
