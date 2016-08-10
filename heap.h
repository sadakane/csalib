/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
typedef struct {
  i64 max_n;
  i64 size;
  i64 item_w;
  i64 key_w;
  i64 key_ofs;
  uchar *A;
  int (*cmp)(uchar *s, uchar *t, i64 key_w);
} HEAP;


void heapify(HEAP *H, i64 i);
void heap_build(HEAP *H, i64 n, uchar *A, i64 max_n,i64 item_w,i64 key_w,i64 key_ofs, int (*cmp)(uchar *, uchar *, i64));
//void heapsort(int n, data *A);
void heap_extract(HEAP *H, uchar *MAX);
void heap_insert(HEAP *H, uchar *key);
void heap_delete(HEAP *H, i64 i);
void heap_print(HEAP *H);
