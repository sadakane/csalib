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
#include "cst.h"

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

#if 0
void traverse(CSA *csa, i64 l, i64 r, i64 th, uchar *key)
{
  int f,c;
  i64 ll, rr;
  i64 num;

  f = 0;
  for (c=0; c<=255; c++) {
    ll = l;  rr = r;
    csa->searchsub(c, csa, &ll, &rr);
    if (ll <= rr) {
      num = rr - ll + 1;
      if (num >= th) {
        key[-1] = c;
        traverse(csa, ll, rr, th, key-1);
        f = 1;
      }
    }
  }
        
  if (f == 0) {
    printf("[%ld,%ld] key = %s\n",l,r,key);
  }
}
#else
void traverse(CSA *csa, i64 l, i64 r, i64 th, uchar *key)
{
  int i,f,c,k, k2;
  i64 ll, rr;
  i64 num;
  uchar head[csa->m];
  i64 lll[csa->m], rrr[csa->m];
  uchar head2[csa->m];
  i64 lll2[csa->m], rrr2[csa->m];

  f = 0;
  k = csa_child_l(csa,l,r,head, lll, rrr);
  k2 = lf_wt_child_l(csa,l,r,head2, lll2, rrr2);

  printf("k = %d ",k);
  for (i=0; i<k; i++) {printf("(%d [%ld,%ld]) ",head[i],lll[i],rrr[i]);}
  printf("\n");
  printf("k2 = %d ",k2);
  for (i=0; i<k2; i++) {printf("(%d [%ld,%ld]) ",head2[i],lll2[i],rrr2[i]);}
  printf("\n");

  for (i=0; i<k; i++) {
    c = head[i];
    ll = lll[i];  rr = rrr[i];
    if (ll <= rr) {
      num = rr - ll + 1;
      if (num >= th) {
        key[-1] = c;
        traverse(csa, ll, rr, th, key-1);
        f = 1;
      }
    }
  }

  if (f == 0) {
    printf("[%ld,%ld] key = %s\n",l,r,key);
  }
}
#endif

void cst_traverse(cst_node node)
{
  cst_node child, parent;
  uchar *pathlabel;

  pathlabel = cst_pathlabel(node);
  printf("cst_traverse node = (%ld, [%ld, %ld]) \"%s\" \n",node.depth, node.l, node.r, pathlabel);
  free(pathlabel);

  if (cst_isleaf(node)) return;

  child = cst_firstchild(node);
  parent = cst_parent(child);
  if (!cst_eq(parent,node)) {
    printf("fistchild = (%ld, [%ld, %ld]) \n",child.depth, child.l, child.r);
    printf("parent = (%ld, [%ld, %ld]) \n",parent.depth, parent.l, parent.r);
  }
  while (!cst_isnull(child)) {
    parent = cst_parent(child);
    if (!cst_eq(parent,node)) {
      printf("nextchild = (%ld, [%ld, %ld]) \n",child.depth, child.l, child.r);
      printf("parent = (%ld, [%ld, %ld]) \n",parent.depth, parent.l, parent.r);
    }
    cst_traverse(child);
    child = cst_nextchild(node, child);
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
//  i64 count[10002];

   if (argc<2) {
      fprintf(stderr, "syntax: suftest file\n");
      return 1;
   }

   //if (argc>=3) sprintf(fname2,"%s",argv[2]);
//   csa_read(&SA,argc, argv);
   csa_read(&SA,2, argv+1);
   n = SA.n;

#if 0
   csa_read(&SA2,2, argv+3);

  for (i=0; i<n; i++) {
    rank_t p1, p2;
    int c1, c2;
    if (i % 100000 == 0) {printf("%ld\r", i);  fflush(stdout);}
    p1 = SA.BW_LF(&SA, i, &c1);
    p2 = SA2.BW_LF(&SA2, i, &c2);
    if (p1 != p2 || c1 != c2) {
      printf("i=%ld LF = %ld LF2 = %ld c1 = %d c2 = %d\n", i, p1, p2, c1, c2);
    }
  }

#endif
  

#if 1
  {
    CSAFILE *csafile;
    rank_t rank, rank2;
    pos_t pos, pos2;
    int c, c2;
    csafile = csa_fdopen(&SA, NULL);
    pos = 0;
    rank = SA.inverse(&SA, pos);
    csa_fseek(csafile, pos, SEEK_SET);
    mygettime(&before);
    for (i=0; i<n/10; i++) {
      if (i % 10000 == 0) {printf("%ld\r", i);  fflush(stdout);}
      c2 = csa_fgetc(csafile);
      pos2 = csafile->pos;
      rank2 = csafile->rank;
#if 1
      c = SA.T(&SA, rank);
      pos++;
      rank = SA.psi(&SA, rank);

      if (c != c2 || pos != pos2 || rank != rank2) {
        printf("pos  = %ld rank  = %ld c  = %d\n", pos, rank, c);
        printf("pos2 = %ld rank2 = %ld c2 = %d\n", pos2, rank2, c2);
      }

#endif

    }
     mygettime(&after);
     t = mylaptime(&before,&after);
     printf("fgetc %f sec\n",t);

    pos = n;
    rank = SA.inverse(&SA, pos);
    csa_fseek(csafile, pos, SEEK_SET);
    mygettime(&before);
    for (i=0; i<n/10; i++) {
      if (i % 10000 == 0) {printf("%ld\r", i);  fflush(stdout);}
      c2 = csa_fgetbw(csafile);
      pos2 = csafile->pos;
      rank2 = csafile->rank;
#if 1
      c = SA.BW(&SA, rank);
      pos--;
      rank = SA.LF(&SA, rank);

      if (c != c2 || pos != pos2 || rank != rank2) {
        printf("pos  = %ld rank  = %ld c  = %d\n", pos, rank, c);
        printf("pos2 = %ld rank2 = %ld c2 = %d\n", pos2, rank2, c2);
      }

#endif

    }
     mygettime(&after);
     t = mylaptime(&before,&after);
     printf("fgetbw %f sec\n",t);


    pos = 0;
    rank = SA.inverse(&SA, pos);
    mygettime(&before);
    for (i=0; i<n/10; i++) {
      if (i % 10000 == 0) {printf("%ld\r", i);  fflush(stdout);}
      c = SA.T(&SA, rank);
      pos++;
      rank = SA.psi(&SA, rank);
    }
     mygettime(&after);
     t = mylaptime(&before,&after);
     printf("psi %f sec\n",t);

    pos = n;
    rank = SA.inverse(&SA, pos);
    mygettime(&before);
    for (i=0; i<n/1000; i++) {
      if (i % 10000 == 0) {printf("%ld\r", i);  fflush(stdout);}
//      c = SA.BW(&SA, rank);
      pos--;
      rank = SA.LF(&SA, rank);
    }
     mygettime(&after);
     t = mylaptime(&before,&after);
     printf("LF %f sec\n",t);

  }
#endif

#if 1
  {
    uchar key[200];
    uchar buf[200];
    i64 keylen;
    i64 l,r;
    i64 s,t,j,k;
    while (1) {
      printf("key ");
      fgets(key,100,stdin);
      keylen = (strlen2(key)-1)/SA.k2;
      s = SA.search(key,keylen,&SA,&l,&r);
      if (s == keylen) {
        printf("number of occurrences = %ld\n",r-l+1);
        if (r-l+1 > 20) r = l+20;
#if 1
        for (i=l; i<=r; i++) {
          j = SA.lookup(&SA,i);
          printf("%12ld: ",j);
//          s = max(j-20,0);
//          t = min(j+20+keylen-1,SA.n-1);
          s = j;
          t = min(j+20,SA.n-1);
          SA.text(buf,&SA,s,t);
          for (k=0; k<(t-s+1)*SA.k2; k++) putchar(isspace(buf[k])?' ':buf[k]);
          printf("\n");
        }
#endif
#if 0
        for (i=l; i<=r; i++) {
          i64 l1,l2;
          l1 = SA.substring_lf(buf+20,&SA,i,20);
          l2 = SA.substring(buf+20,&SA,i,20+keylen);
          for (k=20-l1; k<20+l2; k++) putchar(isspace(buf[k])?' ':buf[k]);
          printf("\n");
        }
#endif
      } else {
        printf("no match max match = %ld\n",s);
      }
    }
  }
#endif



#if 0
  {
    cst_node node;
    node = cst_root(&SA);
    cst_traverse(node);
  }
#endif

#if 0
{
  i64 th;
  uchar key[100];
  n = SA.n;
  th = n / 10000;
  traverse(&SA, 0, n, th, key+99);
  exit(0);
}
#endif

#if 0
  {
    i64 x,j, sum;
     if (SA.psi != csa_psi_by_rankc_naive) {
       mygettime(&before);
       x = 0;
       sum = 0;
       for (i = 1; i <= n; i++) {
         if (i % 1000000 == 0) {fprintf(stderr,"%ld \r",i);  fflush(stderr);}
         SA.head(&SA,x);
         x = SA.psi(&SA,x);
//         printf("i=%ld x=%ld\n",i,x);
         sum += x*i;
       }
       mygettime(&after);
       t = mylaptime(&before,&after);
//       printf("psi %f sec\n",t);
       printf("%f %f\n",(double)SA.psize*8/n, t);
       printf("# sum = %ld\n",sum);
     } else {
       mygettime(&before);
       //x = SA.inverse(&SA,n);
       x = 0;
       sum = 0;
       for (i = n; i > 0; i--) {
         x = SA.LF(&SA,x);
//       printf("LF[%ld] = %ld\n",j,i);
         if (i % 1000000 == 0) {
           fprintf(stderr,"%10u\r",i);
            fflush(stderr);
         }
//         printf("i=%ld x=%ld\n",i,x);
         sum += x*i;
       }
       mygettime(&after);
       t = mylaptime(&before,&after);
//       printf("LF %f sec\n",t);
       printf("%f %f\n",(double)SA.psize*8/n, t);
       printf("# sum = %ld\n",sum);
     }
     return 0;
  }
#endif

#if 0
   {
     i64 j;
     mygettime(&before);
     i = 0;
     for (j=0; j<=min(n,n+100); j++) {
       i = SA.LF(&SA,i);
//       printf("LF[%ld] = %ld\n",j,i);
       if (j % 1000000 == 0) {
         fprintf(stderr,"%10u\r",j);
          fflush(stderr);
       }
     }
     mygettime(&after);
     t = mylaptime(&before,&after);
     printf("LF %f sec\n",t);
   }
#endif


#if 0
   {
     i64 j,m,c;
     i64 freq[SIGMA];
     m = SA.m;
     for (i=0; i<SIGMA; i++) freq[i] = 0;
     for (j=0; j<=n; j++) {
       c = SA.BW(&SA,j);
       printf("BW[%ld] = %ld (%d)\n",j,c,SA.AtoC[c]);
       if (c != -1) {
         freq[c]++;
       }
#if 1
       for (c=0; c<m; c++) {
         if (SA.rankc(&SA,j,c) != freq[c]) {
           printf("rankc(%ld,%d) = %ld  freq = %ld\n",j,c,
                  SA.rankc(&SA,j,c),freq[c]);
         }
       }
#endif
     }
     return 0;
   }
#endif

#if 0
   {
     i64 j,jj;
     mygettime(&before);
     i = 0;
     for (j=0; j<=min(n,n+20); j++) {
       i = SA.psi(&SA,j);
       if (j % 1000000 == 0) {
         fprintf(stderr,"%10u\r",j);
          fflush(stderr);
       }
       jj = SA.LF(&SA,i);
       if (jj != j) {
         printf("error\n");
         printf("psi[%ld] = %ld\n",j,i);
         printf("LF[%ld] = %ld\n",i,jj);
       }
     }
     mygettime(&after);
     t = mylaptime(&before,&after);
     fprintf(stderr,"psi %f sec\n",t);
   }
   return 0;
#endif

#if 0
   {
     i64 j;
     i64 hoge;
     mygettime(&before);
     hoge = 0;
     for (i=0; i<=min(n+100,n); i++) {
       j = SA.lookup(&SA,i);
       printf("sa[%ld] = %ld\n",i,j);
       if (i % 1000 == 0) {
         printf("%10u\r",i);
          fflush(stdout);
       }
       hoge += j;
     }
     mygettime(&after);
     t = mylaptime(&before,&after);
     printf("lookup %f sec sum = %ld\n",t,hoge);
   }
#endif
#if 0
   mygettime(&before);
   for (i = 1; i <= 10000; i++) {
     int c;
     c = csa_T(&SA,i);
   }
   mygettime(&after);
   t = mylaptime(&before,&after);
   printf("T %f sec\n",t);
#endif
#if 0
   mygettime(&before);
   for (i = 0; i <= n; i++) {
     i64 j;
     if (i % 1000 == 0) {
       printf("%u\r",i);
       fflush(stdout);
     }
     j = SA.inverse(&SA,i);
#if 1
     if (SA.lookup(&SA,j) != i) {
       printf("error SA[%ld] = %ld  ISA[%ld] = %ld\n",j,SA.lookup(&SA,j),i,j);
     }
#endif
   }
   mygettime(&after);
   t = mylaptime(&before,&after);
   printf("inverse %f sec\n",t);
#endif

   return 0;
}
