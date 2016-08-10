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
#include "approx.h"

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
  t += (double)(after->ru_utime.tv_usec
                - before->ru_utime.tv_usec)/1000000;
  return t;
}
#endif

void test_approxsearch(CSA *SA)
{
   long i,n_pos;
   int keylen;
   double t, t1, t2;
   unsigned char key[4096];
   char *buf,*c_ptr;
   mytimestruct before,after;
   approx_list *list, *list0;
   int max_error;

   while (!feof(stdin)) {
#if 1
      printf("\ninput error len pos ");  fflush(stdout);
      scanf(" %d %d %ld",&max_error, &keylen, &i);
//      i = rand() % (SA->n-keylen-1) + 1;
      SA->substring(key, SA, i, keylen);
//      i = 200000;
      printf("error = %d len = %d pos = %ld\n", max_error, keylen, i);
      SA->text(key,SA,i,i+keylen-1);
      key[keylen] = 0;
      printf("key %s\n",key);
#endif
#if 0
      printf("\ninput key error ");  fflush(stdout);
      scanf(" %s %d",key, &max_error);
      keylen = strlen(key);
      for (i=0; i<keylen; i++) key[i] = tolower(key[i]);
      printf("key %s\n",key);
#endif
#if 1
      mygettime(&before);
      list = csa_approxsearch(key, keylen, max_error, SA);
      mygettime(&after);
      t1 = mylaptime(&before,&after);
      printf("search1 %f sec\n",t1);
      approx_list_print(SA, list, key, keylen);
#endif
#if 1
      mygettime(&before);
      list0 = csa_approxsearch2(key, keylen, max_error, SA);
      mygettime(&after);
      t2 = mylaptime(&before,&after);
//      approx_list_print(SA, list0, key, keylen);
      printf("search2 %f sec\n",t2);
#endif
#if 0
  if (list0->num != list->num) {
    printf("list0 %d list %d\n", list0->num, list->num);
  }
  for (i=0; i<list0->num && i<list->num; i++) {
    approx_set *p0, *p;
    p0 = &list0->p[i];  p = & list->p[i];
    if (p0->error != p->error || p0->len != p->len || p0->l != p->l || p0->r != p->r) {
      printf("list0 error = %d len = %d [%ld,%ld]\n", p0->error, p0->len, p0->l, p0->r);
      approx_set_print(SA, p0, key, keylen);
      printf("list  error = %d len = %d [%ld,%ld]\n", p->error, p->len, p->l, p->r);
      approx_set_print(SA, p, key, keylen);
    }
  }
  if (list0->num > list->num) {
    approx_set *p0;
    for (i=list->num; i<list0->num; i++) {
      p0 = &list0->p[i];
      printf("list0[%d] error = %d len = %d [%ld,%ld]\n", i, p0->error, p0->len, p0->l, p0->r);
      approx_set_print(SA, &list0->p[i], key, keylen);
    }
  }
  if (list0->num < list->num) {
    approx_set *p;
    for (i=list0->num; i<list->num; i++) {
      p = &list->p[i];
      printf("list[%d] error = %d len = %d [%ld,%ld]\n", i, p->error, p->len, p->l, p->r);
      approx_set_print(SA, &list->p[i], key, keylen);
    }
  }
#endif


//      printf("search results\n");
//      approx_list_print(SA, list, key, keylen);

//      approx_list_free(list);
      //break;
   }
}


int main(int argc, char *argv[])
{
  CSA SA;

  csa_read(&SA,argc-1,argv+1);

  test_approxsearch(&SA);
}
