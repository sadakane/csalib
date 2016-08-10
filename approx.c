/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdio.h>
#include <stdlib.h>
#include "csa.h"
#include "approx.h"

#define INDEX(x,y) ((x)*(1+keylen)+(y))
#define GAP_PENALTY (1)
#define MISS_PENALTY (1)

#define LEFT 0
#define RIGHT 1

static int nnodes;

#define BLOCKSIZE (1<<20)
approx_list *approx_list_init(void)
{
  approx_list *list;
  list = malloc(sizeof(approx_list));
  if (list == NULL) {
    printf("approx_list_init: malloc failed.\n");
    exit(1);
  }
  list->p = malloc(sizeof(approx_set)*BLOCKSIZE);
  if (list->p == NULL) {
    printf("approx_list_init: malloc failed.\n");
    exit(1);
  }
  list->len = BLOCKSIZE;
  list->num = 0;
  return list;
}

void approx_list_append(approx_list *list, int error, int len, i64 l, i64 r)
{
  if (list->num == list->len) {
    list->len += BLOCKSIZE;
//    printf("list->p = %p ",list->p);
    list->p = realloc(list->p, sizeof(approx_set)*list->len);
//    printf("realloc list->p = %p\n",list->p);
  }
  list->p[list->num].error = error;
  list->p[list->num].len = len;
  list->p[list->num].l = l;
  list->p[list->num].r = r;
  list->num++;
}

#define MATCH 0
#define GAP 1
#define MISS 2
int dp(uchar *T, int *tlenp, int *shift, uchar *P, int keylen, int prn)
{
  int i,j,k;
  int *score;
  int *e1, *e2;
  int d;
  int tlen;

  tlen = *tlenp;

  score = malloc(sizeof(int)*(1+keylen)*(1+tlen));
  for (j=0; j<=keylen; j++) score[INDEX(0,j)] = j * GAP_PENALTY;
  for (i=0; i<=tlen; i++) score[INDEX(i,0)] = i * GAP_PENALTY;

  for (i=1; i<=tlen; i++) {
    for (j=1; j<=keylen; j++) {
      score[INDEX(i,j)]
       = min(score[INDEX(i-1,j-1)] + (T[i-1] != P[j-1])*MISS_PENALTY,  // MATCH or MISS
         min(score[INDEX(i-1,j)] +GAP_PENALTY, // GAP in T
             score[INDEX(i,j-1)] +GAP_PENALTY  // GAP in P
            ));
    }
  }
  k = score[INDEX(tlen,keylen)];

  e1 = malloc(sizeof(int)*(tlen+keylen+2));
  e2 = malloc(sizeof(int)*(tlen+keylen+2));

  i = tlen;  j = keylen;
  d = 0;
  while (i>0 || j>0) {
    if (i>0 && j>0 && score[INDEX(i,j)] == score[INDEX(i-1,j-1)] + (T[i-1] != P[j-1])*MISS_PENALTY) {
      if (T[i-1] == P[j-1]) {
        e1[d] = MATCH;
        e2[d] = MATCH;
      } else {
        e1[d] = MISS;
        e2[d] = MISS;
      }
      i--;
      j--;
    }
    if (j>0 && (score[INDEX(i,j)] == score[INDEX(i,j-1)]+GAP_PENALTY)) {
      e1[d] = GAP;
      e2[d] = MATCH;
      j--;
    }
    if (i>0 && (score[INDEX(i,j)] == score[INDEX(i-1,j)]+GAP_PENALTY)) {
      e1[d] = MATCH;
      e2[d] = GAP;
      i--;
    } 
    d++;
  }
  if (prn) {
    printf("T ");
    for (j=0,i=d-1; i>=0; i--) {
      switch(e1[i]) {
      case MISS:
        printf("[%c]", T[j]);
        j++;
        break;
      case GAP:
        printf("-");
        break;
      case MATCH:
        printf("%c", T[j]);
        j++;
        break;
      }
    }
    printf("\n");
    printf("P ");
    for (i=0,j=d-1; j>=0; j--) {
      switch(e2[j]) {
      case MISS:
        printf("[%c]", P[i]);
        i++;
        break;
      case GAP:
        printf("-");
        break;
      case MATCH:
        printf("%c", P[i]);
        i++;
        break;
      }
    }
    printf("\n");
  }
  *shift = 0;
#if 1
  for (i=0,j=d-1; j>=0; j--) {
    if (e2[j] != GAP) break;
    i++;
  }
  *shift = i;
  tlen -= i;

  for (i=0,j=0; j<d; j++) {
    if (e2[j] != GAP) break;
    i++;
  }
  tlen -= i;
  *tlenp = tlen;
  if (prn) printf("shift %d tlen %d\n", *shift, *tlenp);
#endif
  free(e2);  free(e1);

  free(score);
  return k;
}

void approx_set_print(CSA *csa, approx_set *T, uchar *key, int keylen)
{
  int i, j, len, shift;
  uchar *text;

  len = T->len;
  text = mymalloc(len+1);
  csa->substring(text, csa, T->l, len);
  text[len] = 0;
  dp(text, &len, &shift, key, keylen, 1);
  free(text);
}

void approx_list_fix(CSA *csa, approx_list *list, uchar *key, int keylen)
{
  int i, j, len, shift;
  uchar *text;
  for (i=0; i<list->num; i++) {
    len = list->p[i].len;
    text = mymalloc(len+1);
    csa->substring(text, csa, list->p[i].l, len);
    text[len] = 0;
    dp(text, &len, &shift, key, keylen, 0);
    for (j=0; j<shift; j++) {
      list->p[i].l = csa->psi(csa, list->p[i].l);
      list->p[i].r = csa->psi(csa, list->p[i].r);
    }
    list->p[i].len = len;
    csa->substring(text, csa, list->p[i].l, len);
    text[len] = 0;
    list->p[i].error = dp(text, &len, &shift, key, keylen, 0);
    free(text);
  }
}

void approx_list_print(CSA *csa, approx_list *list, uchar *key, int keylen)
{
  int i, j, len, shift;
  uchar *text;
  printf("list num=%d\n", list->num);
  for (i=0; i<list->num; i++) {
    len = list->p[i].len;
    printf("[%ld,%ld] error=%d \n",list->p[i].l,list->p[i].r,list->p[i].error);
    text = mymalloc(len+1);
    csa->substring(text, csa, list->p[i].l, len);
    text[len] = 0;
    dp(text, &len, &shift, key, keylen, 1);
    free(text);
  }
}

void approx_list_clear(approx_list *list)
{
  list->num = 0;
}

void approx_list_free(approx_list *list)
{
  free(list->p);
  free(list);
}

static int sign(i64 x)
{
  if (x < 0) return -1;
  if (x > 0) return 1;
  return 0;
}


int approx_set_cmp(const void *p, const void *q)
{
  approx_set *pp, *qq;
  pp = (approx_set *)p;
  qq = (approx_set *)q;
  if (pp->l != qq->l) return sign(pp->l - qq->l);
  if (pp->error != qq->error) return pp->error - qq->error;
  if (pp->r != qq->r) return sign(pp->r - qq->r);
  if (pp->len != qq->len) return pp->len - qq->len;
  return 0;
}

void approx_list_sort(approx_list *list)
{
  qsort(list->p, list->num, sizeof(*(list->p)), approx_set_cmp);
}

void approx_list_uniq(approx_list *list)
{
  int i,j;
  approx_set *p;
  p = list->p;
  j = 0;
  for (i=1; i<list->num; i++) {
    if (!(p[i].error == p[j].error && p[i].l == p[j].l && p[i].r == p[j].r && p[i].len == p[j].len)){
      p[j+1].error = p[i].error;
      p[j+1].l = p[i].l;
      p[j+1].r = p[i].r;
      p[j+1].len = p[i].len;
      j++;
    }
  }
  list->num = j+1;
}

static int csa_approxsearch_sub(
  uchar *key,int keylen, /* ŒŸõ‚·‚éƒpƒ^ƒ“ P ‚Æ‚»‚Ì’·‚³ */
  int k, /* Å‘å‹——£ */
  int *score, /* DP•\ */
  i64 tl, i64 tr, /* T ‚ÌŽ«‘‡ */
  int tlen, /* T ‚Ì’·‚³ */
  uchar *T, /* T[0..tlen-1] */
  int lr, /* P ‚Ì‰E‚©‚ç’T‚· (RIGHT==1) ¶‚©‚ç’T‚· (LEFT==0) */
  CSA *csa, /* T ‚ÌCSA */
  approx_list *list /* “š‚¦‚ðŠi”[‚·‚éƒŠƒXƒg */
  )
{
  i64 i,j;
  int c,cc;
  i64 lll, rrr;
  uchar set[csa->m];
  rank_t ll[csa->m], rr[csa->m];
  int mm;
  int num;

  tlen++;
  score[INDEX(tlen,0)] = score[INDEX(tlen-1,0)]+1;
#if 0
  if (tl == tr) { // Œó•â‚ª1‚Â‚µ‚©‚È‚¢
    if (lr == LEFT) {
      lll = csa->BW_LF(csa, l, &c);
    } else {
      lll = l;
      for (i=0; i<tlen-1; i++) lll = csa->psi(csa, lll);
      c = csa->T(csa, lll);
      lll = l;
    }
    if (c == -1) {
      num = 0;
    } else {
      num = 1;
      set[0] = c;
      ll[0] = rr[0] = lll;
    }
  } else {
    if (lr == LEFT) {
      num = csa->child_l(csa, tl, tr, set, ll, rr);
    } else {
      num = csa_child_rs(csa, T, tlen-1, set, ll, rr);
    }
  }
#else
  if (lr == LEFT) {
    num = csa->child_l(csa, tl, tr, set, ll, rr);
  } else {
    num = csa_child_rs(csa, T, tlen-1, set, ll, rr);
  }
#endif
  for (i=0; i<num; i++) {
    nnodes++;
    c = set[i];
    lll = ll[i];  rrr = rr[i];
    T[tlen-1] = c;
    mm = score[INDEX(tlen,0)];
    for (j=1; j<=keylen; j++) {
      if (lr == LEFT) cc = key[keylen-j]; else cc = key[j-1];
      score[INDEX(tlen,j)]
        = min(score[INDEX(tlen-1,j-1)] + (c != cc)*MISS_PENALTY, // MATCH or MISS
          min(score[INDEX(tlen-1,j)] +GAP_PENALTY, // GAP in T
              score[INDEX(tlen,j-1)] +GAP_PENALTY));  // GAP in P
      mm = min(mm, score[INDEX(tlen,j)]);
    }
#if 1
    if (mm == k) { // ‚±‚ÌŒã‚ÍŠ®‘Sˆê’v‚¾‚¯
      for (j=0; j<keylen; j++) {
        if (score[INDEX(tlen,j)] == k) {
          int jj;
          i64 l3, r3;
          if (lr == LEFT) {
            l3 = lll;  r3 = rrr;
            for (jj=j+1; jj<=keylen; jj++) {
              nnodes++;
              cc = key[keylen-jj];
              if (csa->searchsub(cc, csa, &l3, &r3) != 0) break;
            }
            if (jj > keylen) { // ‘S•”ˆê’v‚µ‚½
              approx_list_append(list, k, tlen+keylen-j, l3, r3);
            }
          } else {
            jj = csa->search(&key[j], keylen-j, csa, &l3, &r3);
            nnodes += jj;
            if (jj == keylen-j) {
              for (jj = tlen; jj > 0; jj--) {
                nnodes++;
                cc = T[jj-1];
                if (csa->searchsub(cc, csa, &l3, &r3) != 0) break;
              }
              if (jj == 0) { // ‘S•”ˆê’v‚µ‚½
                approx_list_append(list, k, tlen+keylen-j, l3, r3);
              }
            }
          }
        }
      }
    }
#endif
    if (mm < k) { // P ‘S‘Ì‚ª‹——£ k ˆÈ“à‚Åˆê’v‚·‚é‰Â”\«‚ª‚ ‚é
      csa_approxsearch_sub(key,keylen,k,score,lll,rrr,tlen,T,lr,csa,list);
    }
    if (score[INDEX(tlen,keylen)] <= k) { // P ‘S‘Ì‚ª‹——£ k ˆÈ“à‚Åˆê’v‚µ‚½
      approx_list_append(list, score[INDEX(tlen,keylen)], tlen, lll, rrr);
    }
  }
  return 0;
}

approx_list *csa_approxsearch(unsigned char *key,int keylen,int k,CSA *csa)
{
  i64 s1,s2;
  int c,s,t;
  i64 i,j,h,l,r,n;
  int *score;
  uchar *buf;
  int buflen;
  approx_list *list, *list0;

  score = malloc(sizeof(int)*(1+keylen)*(1+keylen+k+1));
  buflen = keylen+k+1;
  buf = malloc(buflen);
//  printf("buf %d\n",buflen);
  for (i=0; i<=keylen; i++) score[INDEX(0,i)] = i * GAP_PENALTY;
  for (i=0; i<=keylen+k; i++) score[INDEX(i,0)] = i * GAP_PENALTY;
  for (i=0; i<=keylen+k; i++) buf[i] = 0;

  l = 0;  r = csa->n;

  nnodes = 0;
  list = approx_list_init();
  csa_approxsearch_sub(key, keylen, k, score, l, r, 0, buf, LEFT, csa, list);
  printf("#visited nodes = %d\n",nnodes);
  approx_list_fix(csa, list, key, keylen);
  approx_list_sort(list);
  approx_list_uniq(list);
//  approx_list_print(csa, list, key, keylen);

  free(score);
  free(buf);

  return list;
}

approx_list *csa_approxsearch2(unsigned char *key,int keylen,int k,CSA *csa)
{
  i64 s1,s2;
  int c,s,t;
  i64 i,j,h,l,r,n;
  int *score;
  uchar *buf;
  int buflen;
  approx_list *list, *list1, *list2;
  
  int l1, l2, e1, e2;

  e1 = k/2;
  e2 = k-e1;
  l1 = keylen/2;
  l2 = keylen-l1;

  score = malloc(sizeof(int)*(1+keylen)*(1+keylen+k+1));
  buflen = keylen+k+1;
  buf = malloc(buflen);
//  printf("buf %d\n",buflen);
  for (i=0; i<=keylen; i++) score[INDEX(0,i)] = i;
  for (i=0; i<=keylen+k; i++) score[INDEX(i,0)] = i;
  for (i=0; i<=keylen+k; i++) buf[i] = 0;

  l = 0;  r = csa->n;

  nnodes = 0;
  list1 = approx_list_init();
  for (j=0; j<=l1; j++) score[INDEX(0,j)] = j;
  csa_approxsearch_sub(key, l1, e1, score, l, r, 0, buf, LEFT, csa, list1);

  list2 = approx_list_init();
  for (i=0; i<list1->num; i++) {
    int e, tlen, keylen;
    i64 ll, rr;
    e = list1->p[i].error;
    tlen = list1->p[i].len;
    ll = list1->p[i].l;
    rr = list1->p[i].r;
    keylen = l2;
    csa->substring(buf, csa, ll, tlen);
    for (j=0; j<=keylen; j++) score[INDEX(tlen,j)] = e+j;
    csa_approxsearch_sub(key+l1, keylen, k, score, ll, rr, tlen, buf, RIGHT, csa, list2);
  }
  approx_list_clear(list1);

  for (j=0; j<=l2; j++) score[INDEX(0,j)] = j;
  csa_approxsearch_sub(key+l1, l2, e2, score, l, r, 0, buf, LEFT, csa, list1);

  for (i=0; i<list1->num; i++) {
    int e, tlen, keylen;
    i64 ll,rr;
    e = list1->p[i].error;
    tlen = list1->p[i].len;
    ll = list1->p[i].l;
    rr = list1->p[i].r;
    keylen = l1;
    for (j=0; j<=keylen; j++) score[INDEX(tlen,j)] = e+j;
    csa_approxsearch_sub(key, keylen, k, score, ll, rr, tlen, buf, LEFT, csa, list2);
  }

  approx_list_fix(csa, list2, key, keylen);
  approx_list_sort(list2);
  approx_list_free(list1);

  approx_list_uniq(list2);

  printf("#visited nodes = %d\n",nnodes);

  free(score);
  free(buf);
  
  return list2;
}
